##### Libs #####
rm(list=ls())
library(corrplot)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(viridis)
library(TMB)
library(haven)
library(car)
library(gbm)
library(Matrix)
library(lme4)
library(caret)
library(mgcv)
library(nlme)
library(e1071)
library(minpack.lm)
library(ggpubr)
library(MASS)
library(patchwork)
library(tidyr)
library(readr)
library(grid)


###=====================================================================================================
###------------------------- Natural origin Jack pine -----------------------------------------------

# data

n <- read.csv("Natural origin jack pine data file")
############################## Parametric modeling ########################

set.seed(1001)

predictor_sets <- list(
  DBH                 = c("DBH"),
  DBH_SHT             = c("DBH", "StandHtDomOrCodom"),
  DBH_SHT_BAH         = c("DBH", "StandHtDomOrCodom", "TreeBALiveTotalPHa"),
  DBH_SHT_BAH_TPH     = c("DBH", "StandHtDomOrCodom", "TreeBALiveTotalPHa", "TreeDensityLiveTotalPHa"),
  DBH_BAH_TPH         = c("DBH", "TreeBALiveTotalPHa", "TreeDensityLiveTotalPHa")
)

cols_needed_all <- unique(c("PlotName","HtLiveBranch","HtTot", "CrownClassCode" , unlist(predictor_sets)))
n0 <- na.omit(n[, cols_needed_all])

#### Data summaries
names(n0)
unique(n0$PlotName)
summary_n0 <- n0 %>%
  summarise(
    n_Plots = n_distinct(PlotName),
    n_Trees = n(),
    HTLCB_mean = mean(HtLiveBranch, na.rm = TRUE),
    HTLCB_min  = min(HtLiveBranch, na.rm = TRUE),
    HTLCB_max  = max(HtLiveBranch, na.rm = TRUE),
    HTLCB_sd = sd (HtLiveBranch, na.rm = TRUE),
    DBH_mean = mean(DBH, na.rm = TRUE),
    DBH_min  = min(DBH, na.rm = TRUE),
    DBH_max  = max(DBH, na.rm = TRUE),
    DBH_sd  = sd(DBH, na.rm = TRUE),
    Ht_mean = mean(HtTot, na.rm = TRUE),
    Ht_min  = min(HtTot, na.rm = TRUE),
    Ht_max  = max(HtTot, na.rm = TRUE),
    Ht_sd  = sd(HtTot, na.rm = TRUE),
    SHT_mean = mean(StandHtDomOrCodom, na.rm = TRUE),
    SHT_min  = min(StandHtDomOrCodom, na.rm = TRUE),
    SHT_max  = max(StandHtDomOrCodom, na.rm = TRUE),
    SHT_sd  = sd(StandHtDomOrCodom, na.rm = TRUE),
    BAH_mean = mean(TreeBALiveTotalPHa, na.rm = TRUE),
    BAH_min  = min(TreeBALiveTotalPHa, na.rm = TRUE),
    BAH_max  = max(TreeBALiveTotalPHa, na.rm = TRUE),
    BAH_sd  = sd(TreeBALiveTotalPHa, na.rm = TRUE),
    TPH_mean = mean(TreeDensityLiveTotalPHa, na.rm = TRUE),
    TPH_min  = min(TreeDensityLiveTotalPHa, na.rm = TRUE),
    TPH_max  = max(TreeDensityLiveTotalPHa, na.rm = TRUE),
    TPH_sd  = sd(TreeDensityLiveTotalPHa, na.rm = TRUE)
  )

summary_n0

clamp_ratio <- function(df, eps = 1e-8) {
  r <- pmin(1 - eps, pmax(eps, df$HtLiveBranch / pmax(eps, df$HtTot)))
  df$._ratio_ <- r
  df
}

metrics_fullIC <- function(obs, pred, k_params) {
  n   <- length(obs)
  rss <- sum((obs - pred)^2)
  r2   <- cor(obs, pred)^2
  rmse <- sqrt(rss / n)
  bias <- mean(pred - obs)
  neg2loglik <- n * (log(2*pi) + 1 + log(rss / n))
  aic <- neg2loglik + 2 * k_params
  bic <- neg2loglik + log(n) * k_params
  c(R2 = r2, RMSE = rmse, Bias = bias, AIC = aic, BIC = bic)
}

bind_rows_safe <- function(dflist) {
  if (!length(dflist)) return(data.frame())
  all_cols <- unique(unlist(lapply(dflist, names)))
  padded <- lapply(dflist, function(df) {
    miss <- setdiff(all_cols, names(df))
    if (length(miss)) df[miss] <- NA
    df[all_cols]
  })
  do.call(rbind, padded)
}

.link_eta_start <- function(df, mt, vars) {
  eps <- 1e-8
  r <- pmin(1 - eps, pmax(eps, df$HtLiveBranch / pmax(eps, df$HtTot)))
  X <- as.data.frame(df[, vars, drop = FALSE])
  X <- cbind(Intercept = 1, X)
  
  y <- switch(mt,
              M1 = -log(1 - r),             
              M2 = log(r / (1 - r)),        
              M3 = -log(1/r - 1)            
  )
  ols <- try(lm(y ~ . - 1, data = as.data.frame(X)), silent = TRUE)
  if (inherits(ols, "try-error") || any(!is.finite(coef(ols)))) return(NULL)
  
  co <- as.numeric(coef(ols)); names(co) <- colnames(X)
  b0 <- co["Intercept"]; if (is.na(b0)) b0 <- 0
  bj <- co[setdiff(names(co), "Intercept")]
  names(bj) <- paste0("b", seq_along(vars))
  list(b0 = b0, bj = bj)
}

build_formula_M <- function(model_type, vars) {
  et <- paste(c("b0", paste0("b", seq_along(vars), " * ", vars)), collapse = " + ")
  switch(model_type,
         M1 = as.formula(paste0("HtLiveBranch ~ HtTot * (1 - exp(-(", et, ")))")),
         M2 = as.formula(paste0("HtLiveBranch ~ HtTot / (1 + exp(-(", et, ")))")),
         M3 = as.formula(paste0("HtLiveBranch ~ HtTot / ((1 + exp(-(", et, ")))^(1/m))"))
  )
}

fit_M_once <- function(train, test) {
  out <- list()
  for (mt in c("M1","M2","M3")) {
    for (ps_name in names(predictor_sets)) {
      vars <- predictor_sets[[ps_name]]
      model_id <- paste(mt, ps_name, sep = "_")
      fml  <- build_formula_M(mt, vars)
      
      st_lin  <- .link_eta_start(train, mt, vars)
      k       <- length(vars)
      st_base <- c(b0 = -0.05, setNames(rep(0.05, k), paste0("b", seq_len(k))))
      if (!is.null(st_lin)) {
        st_base["b0"] <- st_lin$b0
        for (j in seq_along(vars)) {
          nm <- paste0("b", j)
          if (!is.na(st_lin$bj[nm])) st_base[nm] <- st_lin$bj[nm]
        }
      }
      if (mt == "M3") st_base <- c(st_base, m = 1.0)
      

      start_sets <- list(
        st_base,
        { tmp <- st_base; tmp["b0"] <- tmp["b0"] * 0.5; tmp },
        { tmp <- st_base; tmp[grep("^b\\d+$", names(tmp))] <- 0; tmp },
        { tmp <- st_base; tmp[grep("^b\\d+$", names(tmp))] <- tmp[grep("^b\\d+$", names(tmp))] * 0.5; tmp }
      )
      
      fit <- NULL
      for (st in start_sets) {
        if (mt == "M3" && !("m" %in% names(st))) st <- c(st, m = 1.0)
        fit <- try(
          nlsLM(fml, data = train, start = st,
                control = nls.lm.control(maxiter = 600, ftol = 1e-10, ptol = 1e-10)),
          silent = TRUE
        )
        if (!inherits(fit, "try-error")) break
      }
      if (inherits(fit, "try-error")) {
        row <- data.frame(
          Model = model_id, Type = "Fixed", Predictors = paste(vars, collapse = "+"),
          R2 = NA_real_, RMSE = NA_real_, Bias = NA_real_, AIC = NA_real_, BIC = NA_real_,
          R2_test = NA_real_, RMSE_test = NA_real_, Bias_test = NA_real_,
          b0 = NA_real_, check.names = FALSE
        )
        for (j in seq_along(vars)) row[[paste0("b", j)]] <- NA_real_
        if (mt == "M3") row[["m"]] <- NA_real_ 
        out[[length(out)+1]] <- row
        next
      }
      
      pred_tr <- as.numeric(predict(fit, newdata = train))
      met_tr  <- c(R2 = cor(train$HtLiveBranch, pred_tr)^2,
                   RMSE = sqrt(mean((train$HtLiveBranch - pred_tr)^2)),
                   Bias = mean(pred_tr - train$HtLiveBranch),
                   AIC = AIC(fit), BIC = BIC(fit))
      
      cf <- coef(fit)
      eta_fun <- function(d) {
        eta <- cf["b0"]
        for (j in seq_along(vars)) eta <- eta + cf[paste0("b", j)] * d[[vars[j]]]
        eta
      }
      pred_te <- switch(mt,
                        M1 = test$HtTot * (1 - exp(-eta_fun(test))),
                        M2 = test$HtTot / (1 + exp(-eta_fun(test))),
                        M3 = { mpar <- cf["m"]; test$HtTot / ((1 + exp(-eta_fun(test)))^(1/mpar)) }
      )
      met_te <- c(R2_test = cor(test$HtLiveBranch, pred_te)^2,
                  RMSE_test = sqrt(mean((test$HtLiveBranch - pred_te)^2)),
                  Bias_test = mean(pred_te - test$HtLiveBranch))
      
      row <- data.frame(
        Model = model_id, Type = "Fixed", Predictors = paste(vars, collapse = "+"),
        t(met_tr), t(met_te), check.names = FALSE
      )
      param_names <- c("b0", paste0("b", seq_along(vars)), if (mt=="M3") "m")
      for (nm in param_names) row[[nm]] <- if (nm %in% names(cf)) cf[[nm]] else NA_real_
      out[[length(out)+1]] <- row
    }
  }
  bind_rows_safe(out)
}

grouped_kfold_indices <- function(df, K = 10) {
  pl <- sample(unique(df$PlotName))
  split(pl, cut(seq_along(pl), K, labels = FALSE))
}

run_grouped_kfold <- function(df, K = 10) {
  folds <- grouped_kfold_indices(df, K)
  per_fold <- list()
  for (k in seq_along(folds)) {
    test_plots  <- folds[[k]]
    train_plots <- setdiff(unique(df$PlotName), test_plots)
    train <- clamp_ratio(df[df$PlotName %in% train_plots, ])
    test  <- clamp_ratio(df[df$PlotName %in% test_plots, ])
    res   <- fit_M_once(train, test)         # << ONLY M1/M2/M3
    res$Fold <- k
    per_fold[[k]] <- res
  }
  bind_rows_safe(per_fold)
}

# Run CV
K <- 10
cv_results <- run_grouped_kfold(n0, K)

agg <- do.call(rbind, lapply(split(cv_results, list(cv_results$Model, cv_results$Predictors), drop = TRUE), function(d) {
  data.frame(
    Model = d$Model[1],
    Predictors = d$Predictors[1],
    Folds = nrow(d),
    # Train 
    R2_mean = mean(d$R2, na.rm = TRUE),
    RMSE_mean = mean(d$RMSE, na.rm = TRUE),
    Bias_mean = mean(d$Bias, na.rm = TRUE),
    AIC_mean = mean(d$AIC, na.rm = TRUE),
    BIC_mean = mean(d$BIC, na.rm = TRUE),
    # Test
    R2_test_mean   = mean(d$R2_test, na.rm = TRUE),
    R2_test_sd     = sd(d$R2_test, na.rm = TRUE),
    RMSE_test_mean = mean(d$RMSE_test, na.rm = TRUE),
    RMSE_test_sd   = sd(d$RMSE_test, na.rm = TRUE),
    Bias_test_mean = mean(d$Bias_test, na.rm = TRUE),
    Bias_test_sd   = sd(d$Bias_test, na.rm = TRUE),
    Family = sub("_.*$", "", d$Model[1]),
    stringsAsFactors = FALSE
  )
}))
agg <- agg[order(agg$RMSE_test_mean), ]
print(agg)


predict_row_on_data <- function(row, data) {
  fam  <- sub("_.*$", "", row$Model)   # M1/M2/M3
  vars <- strsplit(row$Predictors, "\\+")[[1]]
  
  b0 <- suppressWarnings(as.numeric(row$b0)); if (is.na(b0)) b0 <- 0
  bj <- sapply(seq_along(vars), function(j) suppressWarnings(as.numeric(row[[paste0("b", j)]])))
  bj[is.na(bj)] <- 0
  
  eta <- rep(b0, nrow(data))
  for (j in seq_along(vars)) eta <- eta + bj[j] * data[[vars[j]]]
  
  switch(fam,
         M1 = data$HtTot * (1 - exp(-eta)),
         M2 = data$HtTot / (1 + exp(-eta)),
         M3 = {
           mpar <- suppressWarnings(as.numeric(row$m)); if (is.na(mpar)) mpar <- 1.0
           data$HtTot / ((1 + exp(-eta))^(1/mpar))
         },
         NA_real_
  )
}

set.seed(1001)
folds <- grouped_kfold_indices(n0, K)

pred_plot_df <- list()
for (k in seq_along(folds)) {
  test_plots <- folds[[k]]
  test_df    <- n0[n0$PlotName %in% test_plots, , drop = FALSE]
  rows_k <- cv_results[cv_results$Fold == k, , drop = FALSE]
  if (!nrow(rows_k)) next
  
  for (i in seq_len(nrow(rows_k))) {
    rowi <- rows_k[i, ]
    pred <- try(predict_row_on_data(rowi, test_df), silent = TRUE)
    if (inherits(pred, "try-error") || all(!is.finite(pred))) next
    pred_plot_df[[length(pred_plot_df)+1]] <- data.frame(
      Model     = rowi$Model,
      Observed  = test_df$HtLiveBranch,
      Predicted = as.numeric(pred),
      Height = test_df$HtTot,
      CrownClass = test_df$CrownClassCode,
      plotname = test_df$PlotName,
      dbh = test_df$DBH,
      SHT = test_df$StandHtDomOrCodom,
      BAH = test_df$TreeBALiveTotalPHa,
      TPH = test_df$TreeDensityLiveTotalPHa
    )
  }
}
pred_plot_df <- if (length(pred_plot_df)) do.call(rbind, pred_plot_df) else data.frame()

######  Split Model into Family and Set _ PLOTS ####
pred_plot_df2 <- pred_plot_df %>%
  mutate(
    Family = sub("_.*$", "", Model),             # M1/M2/M3
    Set    = sub("^[^_]+_", "", Model)
  )

family_levels <- c("M1","M2","M3")
set_levels    <- c("DBH","DBH+SHT","DBH+SHT+BAH","DBH+SHT+BAH+TPH","DBH+BAH+TPH")
pred_plot_df2 <- pred_plot_df2 %>%
  filter(Family %in% family_levels) %>%
  mutate(Family = factor(Family, levels = family_levels),
         Set    = factor(Set,    levels = set_levels))

rng <- range(c(pred_plot_df2$Observed, pred_plot_df2$Predicted), na.rm = TRUE, finite = TRUE)
pad <- 0.05 * diff(rng); if (!is.finite(pad)) pad <- 0
lims <- c(rng[1] - pad, rng[2] + pad)


names(pred_plot_df2)

# changed _ with +

stats_df0 <- pred_plot_df2 %>%  group_by(Family, Set) %>%
  summarise(    R2 = round(cor(Observed, Predicted, use = "complete.obs")^2, 3),
                RMSE = round(sqrt(mean((Observed - Predicted)^2, na.rm = TRUE)), 3),
                .groups = "drop"  )


tiff(filename = "ObsPred_N_2.jpg", width = 240, height = 180,
     units = "mm", res = 400, 
     compression = "lzw") 
ggplot(pred_plot_df2, aes(x = Observed, y = Predicted)) +
  geom_bin2d(bins = 35) +  scale_fill_gradientn(
    colours = c("#c6dbef", "#6baed9", "#2171b9"),
    trans = "sqrt",
    limits = range(1, 55),  
    breaks = c(1, 55),      
    labels = c("1", "50"), 
    name = NULL,             
    guide = guide_colorbar(      direction = "horizontal",
                                 barheight = unit(0.3, "cm"),      barwidth = unit(4, "cm"),
                                 ticks = FALSE,      label.position = "bottom",
                                 title.position = "top"    )  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_smooth(
    method = "lm",    formula = y ~0+ x,    se = TRUE,    color = "red",
    fill = "red",    alpha = 0.10,    linewidth = 0.5  ) +
  facet_grid(Family ~ Set, drop = FALSE) +
  coord_equal(xlim = lims, ylim = lims, expand = FALSE) +
  labs(title = "Natural origin stand",
       x = "Observed HCB (m)",       y = "Predicted HCB (m)") +
  theme_bw(base_size = 11) +
  theme(    panel.grid.minor = element_blank(),
            strip.text = element_text(face = "bold"),
            panel.spacing = unit(0.6, "lines"),
            legend.position = "bottom",  # move legend below
            legend.justification = "left"  )+
  geom_text(    data = stats_df0,
                aes(x = Inf, y = -Inf,
                    label = paste0("R² = ", R2, "\nRMSE = ", RMSE)),
                hjust = 1.1, vjust = -0.5,
                size = 3.2,
                color = "black",
                inherit.aes = FALSE  )

dev.off()


qn<- pred_plot_df2
names(qn)
qn$HTLCBr <- (qn$Predicted *100) / qn$Height
qn$HTLCBrOBS <- (qn$Observed *100) / qn$Height
qn$res <- qn$Observed - qn$Predicted
qn$res_std <- (qn$Observed - qn$Predicted) / sd(qn$Observed - qn$Predicted, na.rm = TRUE)
sd(qn$res_std, na.rm = TRUE)
#res
lim_res <- range(qn$res_std, na.rm = TRUE)
tiff(filename = "res_n.jpg", width = 240, height = 180,
     units = "mm", res = 400, 
     compression = "lzw") 
qn$Set <- factor(qn$Set, levels = c(
  "DBH",  "DBH+SHT",  "DBH+SHT+BAH",  "DBH+SHT+BAH+TPH",  "DBH+BAH+TPH"))

ggplot(qn, aes(x = Predicted, y = res_std)) +
  geom_bin2d(bins = 35) +  scale_fill_gradientn(
    colours = c("#c6dbef", "#6baed9", "#2171b9"),
    trans = "sqrt",
    limits = range(1, 55),  
    breaks = c(1, 55),      
    labels = c("1", "50"), 
    name = NULL,             
    guide = guide_colorbar(      direction = "horizontal",
                                 barheight = unit(0.3, "cm"),      barwidth = unit(4, "cm"),
                                 ticks = FALSE,      label.position = "bottom",
                                 title.position = "top"    )  ) +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed", color = "black") +
  facet_grid(Family ~ Set, drop = FALSE) +
  labs(title = "Natural origin stand",
       x = "Predicted HCB (m)",       y = "Standartized residuals (m)") +
  theme_bw(base_size = 11) +
  theme(    panel.grid.minor = element_blank(),
            strip.text = element_text(face = "bold"),
            panel.spacing = unit(0.6, "lines"),
            legend.position = "bottom",  # move legend below
            legend.justification = "left"  )

dev.off()


################### Mixed effect M1 V4 ##########
# Using nlme

names(n)
v4 <- c("DBH","StandHtDomOrCodom","TreeBALiveTotalPHa","TreeDensityLiveTotalPHa")
dat <- na.omit(n[, c("PlotName","HtLiveBranch","HtTot", v4)])
names(dat)

model_formula <- function(b0, b1, b2, b3, b4,DBH, StandHtDomOrCodom,
                          TreeBALiveTotalPHa, TreeDensityLiveTotalPHa, HtTot) {
  HtTot * (1 - exp(-(b0 +b1 * DBH +b2 * StandHtDomOrCodom +
                       b3 * TreeBALiveTotalPHa +b4 * TreeDensityLiveTotalPHa)))}

start_vals <- c(b0 = -1,b1 = -0.03, b2 = 0.07, b3 = 0.02, b4 = -0.00003)

mfix <- nls(  HtLiveBranch ~ HtTot * (1 - exp(-(b0 +
                                                  b1 * DBH +b2 * StandHtDomOrCodom +
                                                  b3 * TreeBALiveTotalPHa +b4 * TreeDensityLiveTotalPHa))),
              data = dat,  start = start_vals,  na.action = na.omit)

{sse <- sum(residuals(mfix)^2)
  css <- sum((dat$HtLiveBranch - mean(dat$HtLiveBranch))^2)
  rsq <- 1 - sse / css
  print(rsq)
  error_var <- mfix$sigma^2
  print(error_var)
  coef<-mfix$coefficients
  print(coef)
  RMSE <- sqrt(mean(residuals(mfix)^2))
  print(RMSE)}

mb0 <- nlme(
  HtLiveBranch ~ model_formula(b0, b1, b2, b3, b4, DBH, StandHtDomOrCodom, TreeBALiveTotalPHa, TreeDensityLiveTotalPHa, HtTot),
  data = dat,  fixed = b0 + b1 + b2 + b3 + b4 ~ 1,
  random = b0 ~ 1 | PlotName, 
  start = start_vals,  na.action = na.omit,
  control = nlmeControl(pnlsTol = 1e-6, maxIter = 400),
  method = "REML")

{sse <- sum(residuals(mb0)^2)
  css <- sum((dat$HtLiveBranch - mean(dat$HtLiveBranch))^2)
  rsq <- 1 - sse / css
  print(rsq)
  error_var <- mb0$sigma^2
  print(error_var)
  coef<-mb0$coefficients
  rand <- coef$random$Plot
  rand_var <- var(rand)  
  print(rand_var)
  RMSE <- sqrt(mean(residuals(mb0)^2))
  print(RMSE)}

summary(mb0)
coefmb0<-mb0$coefficients
coefmb0
# Fixed effects
fixed_df <- data.frame(term = names(coefmb0$fixed),
                       estimate = coefmb0$fixed)

mb1 <- nlme(
  HtLiveBranch ~ model_formula(b0, b1, b2, b3, b4, DBH, StandHtDomOrCodom, TreeBALiveTotalPHa, TreeDensityLiveTotalPHa, HtTot),
  data = dat,  fixed = b0 + b1 + b2 + b3 + b4 ~ 1,
  random = b1 ~ 1 | PlotName, 
  start = start_vals,  na.action = na.omit,
  control = nlmeControl(pnlsTol = 1e-6, maxIter = 400), method = "REML")

{sse <- sum(residuals(mb1)^2)
  css <- sum((dat$HtLiveBranch - mean(dat$HtLiveBranch))^2)
  rsq <- 1 - sse / css
  print(rsq)
  error_var <- mb1$sigma^2
  print(error_var)
  coef<-mb1$coefficients
  rand <- coef$random$Plot
  rand_var <- var(rand)  
  print(rand_var)
  RMSE <- sqrt(mean(residuals(mb1)^2))
  print(RMSE)}

mb2 <- nlme(
  HtLiveBranch ~ model_formula(b0, b1, b2, b3, b4, DBH, StandHtDomOrCodom, TreeBALiveTotalPHa, TreeDensityLiveTotalPHa, HtTot),
  data = dat,  fixed = b0 + b1 + b2 + b3 + b4 ~ 1,
  random = b2 ~ 1 | PlotName, 
  start = start_vals,  na.action = na.omit,
  control = nlmeControl(pnlsTol = 1e-6, maxIter = 400))

{sse <- sum(residuals(mb2)^2)
  css <- sum((dat$HtLiveBranch - mean(dat$HtLiveBranch))^2)
  rsq <- 1 - sse / css
  print(rsq)
  error_var <- mb2$sigma^2
  print(error_var)
  coef<-mb2$coefficients
  rand <- coef$random$Plot
  rand_var <- var(rand)  
  print(rand_var)
  RMSE <- sqrt(mean(residuals(mb2)^2))
  print(RMSE)}

mb3 <- nlme(
  HtLiveBranch ~ model_formula(b0, b1, b2, b3, b4, DBH, StandHtDomOrCodom, TreeBALiveTotalPHa, TreeDensityLiveTotalPHa, HtTot),
  data = dat,  fixed = b0 + b1 + b2 + b3 + b4 ~ 1,
  random =  b3 ~ 1 | PlotName, 
  start = start_vals,  na.action = na.omit,
  control = nlmeControl(pnlsTol = 1e-6, maxIter = 400))

{sse <- sum(residuals(mb3)^2)
  css <- sum((dat$HtLiveBranch - mean(dat$HtLiveBranch))^2)
  rsq <- 1 - sse / css
  print(rsq)
  error_var <- mb3$sigma^2
  print(error_var)
  coef<-mb3$coefficients
  rand <- coef$random$Plot
  rand_var <- var(rand)  
  print(rand_var)
  RMSE <- sqrt(mean(residuals(mb3)^2))
  print(RMSE)}

mb4 <- nlme(
  HtLiveBranch ~ model_formula(b0, b1, b2, b3, b4, DBH, StandHtDomOrCodom, TreeBALiveTotalPHa, TreeDensityLiveTotalPHa, HtTot),
  data = dat,  fixed = b0 + b1 + b2 + b3 + b4 ~ 1,
  random =  b4 ~ 1 | PlotName, 
  start = start_vals,  na.action = na.omit,
  control = nlmeControl(pnlsTol = 1e-6, maxIter = 400))

{sse <- sum(residuals(mb4)^2)
  css <- sum((dat$HtLiveBranch - mean(dat$HtLiveBranch))^2)
  rsq <- 1 - sse / css
  print(rsq)
  error_var <- mb4$sigma^2
  print(error_var)
  coef<-mb4$coefficients
  rand <- coef$random$Plot
  rand_var <- var(rand)  
  print(rand_var)
  RMSE <- sqrt(mean(residuals(mb4)^2))
  print(RMSE)}

mb2b3 <- nlme(
  HtLiveBranch ~ model_formula(b0, b1, b2, b3, b4, DBH, StandHtDomOrCodom, TreeBALiveTotalPHa, TreeDensityLiveTotalPHa, HtTot),
  data = dat,  fixed = b0 + b1 + b2 + b3 + b4 ~ 1,
  random = b2 + b3 ~ 1 | PlotName, 
  start = start_vals,  na.action = na.omit,
  control = nlmeControl(pnlsTol = 1e-6, maxIter = 400))

{sse <- sum(residuals(mb2b3)^2)
  css <- sum((dat$HtLiveBranch - mean(dat$HtLiveBranch))^2)
  rsq <- 1 - sse / css
  print(rsq)
  error_var <- mb2b3$sigma^2
  print(error_var)
  coef<-mb2b3$coefficients
  rand <- coef$random$Plot
  rand_var <- var(rand)  
  print(rand_var)
  RMSE <- sqrt(mean(residuals(mb2b3)^2))
  print(RMSE)}

mb3b4 <- nlme(
  HtLiveBranch ~ model_formula(b0, b1, b2, b3, b4, DBH, StandHtDomOrCodom, TreeBALiveTotalPHa, TreeDensityLiveTotalPHa, HtTot),
  data = dat,  fixed = b0 + b1 + b2 + b3 + b4 ~ 1,
  random =  b3 + b4 ~ 1 | PlotName, 
  start = start_vals,  na.action = na.omit,
  control = nlmeControl(pnlsTol = 1e-6, maxIter = 400))

{sse <- sum(residuals(mb3b4)^2)
  css <- sum((dat$HtLiveBranch - mean(dat$HtLiveBranch))^2)
  rsq <- 1 - sse / css
  print(rsq)
  error_var <- mb3b4$sigma^2
  print(error_var)
  coef<-mb3b4$coefficients
  rand <- coef$random$Plot
  rand_var <- var(rand)  
  print(rand_var)
  RMSE <- sqrt(mean(residuals(mb3b4)^2))
  print(RMSE)}

mb2b4 <- nlme(
  HtLiveBranch ~ model_formula(b0, b1, b2, b3, b4, DBH, StandHtDomOrCodom, TreeBALiveTotalPHa, TreeDensityLiveTotalPHa, HtTot),
  data = dat,  fixed = b0 + b1 + b2 + b3 + b4 ~ 1,
  random = b2 + b4 ~ 1 | PlotName, 
  start = start_vals,  na.action = na.omit,
  control = nlmeControl(pnlsTol = 1e-6, maxIter = 400))

{sse <- sum(residuals(mb2b4)^2)
  css <- sum((dat$HtLiveBranch - mean(dat$HtLiveBranch))^2)
  rsq <- 1 - sse / css
  print(rsq)
  error_var <- mb2b4$sigma^2
  print(error_var)
  coef<-mb2b4$coefficients
  rand <- coef$random$Plot
  rand_var <- var(rand)  
  print(rand_var)
  RMSE <- sqrt(mean(residuals(mb2b4)^2))
  print(RMSE)}

mb2b3b4 <- nlme(
HtLiveBranch ~ model_formula(b0, b1, b2, b3, b4, DBH, StandHtDomOrCodom, TreeBALiveTotalPHa, TreeDensityLiveTotalPHa, HtTot),
data = dat,  fixed = b0 + b1 + b2 + b3 + b4 ~ 1,
random = b2 + b3 + b4 ~ 1 | PlotName, 
start = start_vals,  na.action = na.omit,
control = nlmeControl(pnlsTol = 1e-6, maxIter = 1000))


dat2 <- within(dat, { DBH_s = DBH / 10
SHT_s = StandHtDomOrCodom / 10   
BAH_s = TreeBALiveTotalPHa / 10 
TPH_s = TreeDensityLiveTotalPHa / 1000  })

model_formula <- function(b0, b1, b2, b3, b4,
                          DBH_s, SHT_s, BAH_s, TPH_s, HtTot) {
  HtTot * (1 - exp(-(b0 +b1 * DBH_s +b2 * SHT_s +b3 * BAH_s +b4 * TPH_s)))}

start_vals <- c(  b0 = -1,  b1 = -0.2,  b2 = 0.4, b3 = 0.1,b4 = -0.05)

mb2b3b4 <- nlme(  HtLiveBranch ~ model_formula(b0, b1, b2, b3, b4,
                                               DBH_s, SHT_s, BAH_s, TPH_s, HtTot),
                  data   = dat2,  fixed  = b0 + b1 + b2 + b3 + b4 ~ 1,
                  random = b2 + b3 + b4 ~ 1 | PlotName,
                  start  = start_vals,  control = nlmeControl(
                    pnlsTol    = 1e-6,msMaxIter  = 200,maxIter    = 200,
                    pnlsMaxIter = 20  ))

{sse <- sum(residuals(mb2b3b4)^2)
  css <- sum((dat$HtLiveBranch - mean(dat$HtLiveBranch))^2)
  rsq <- 1 - sse / css
  print(rsq)
  error_var <- mb2b3b4$sigma^2
  print(error_var)
  coef<-mb2b3b4$coefficients
  rand <- coef$random$Plot
  rand_var <- var(rand)  
  print(rand_var)
  RMSE <- sqrt(mean(residuals(mb2b3b4)^2))
  print(RMSE)}



anova(mb0, mb1) #, mb2, mb3, mb2b3, mb3b4, mb2b4, mb2b3b4)
AIC(mfix, mb0)#, mb1,mb2, mb3, mb2b3, mb3b4, mb2b4, mb2b3b4)
BIC(mfix, mb0, mb1)#,mb2, mb3, mb2b3, mb3b4, mb2b4, mb2b3b4)

# List of models
models <- list(mb0, mb1, mb2, mb3, mb2b3, mb3b4, mb2b4, mb2b3b4)
model_names <- c("mb0", "mb1", "mb2", "mb3", "mb2b3", "mb3b4", "mb2b4", "mb2b3b4")

null_logLik <- as.numeric(logLik(mfix))

get_pseudo_r2 <- function(model) {
  tryCatch({    ll <- as.numeric(logLik(model))
  r2 <- 1 - (ll / null_logLik)
  return(r2)  }, error = function(e) {    return(NA)  })}

pseudo_r2_values <- sapply(models, get_pseudo_r2)
names(pseudo_r2_values) <- model_names

print(pseudo_r2_values)

#End
######################## Random effects - trees per plot ##############
dat <- na.omit(dat)
dat$PlotName <- as.factor(dat$PlotName)
dat$RowID    <- seq_len(nrow(dat))

model_formula <- function(b0, b1, b2, b3, b4,
                          DBH, StandHtDomOrCodom,
                          TreeBALiveTotalPHa, TreeDensityLiveTotalPHa, HtTot) {
  HtTot * (1 - exp(-(b0 +b1 * DBH +b2 * StandHtDomOrCodom +
                       b3 * TreeBALiveTotalPHa +b4 * TreeDensityLiveTotalPHa)))}

start_vals <- c(b0 = -1,b1 = -0.03,b2 =0.07,b3 =0.02,b4 = -0.00003)

mb0_full <- nlme(
  HtLiveBranch ~ model_formula(b0, b1, b2, b3, b4,
                               DBH, StandHtDomOrCodom,TreeBALiveTotalPHa, TreeDensityLiveTotalPHa, HtTot),
  data   = dat,
  fixed  = b0 + b1 + b2 + b3 + b4 ~ 1,
  random = b0 ~ 1 | PlotName,           # random intercept
  start  = start_vals,
  control = nlmeControl(
    pnlsTol= 1e-6,msMaxIter= 200,maxIter= 200,pnlsMaxIter = 20), method = "REML")

summary(mb0_full)

{sse <- sum(residuals(mb0_full)^2)
  css <- sum((dat$HtLiveBranch - mean(dat$HtLiveBranch))^2)
  rsq <- 1 - sse / css
  print(rsq)
  error_var <- mb0_full$sigma^2
  print(error_var)
  coef<-mb0_full$coefficients
  rand <- coef$random$Plot
  rand_var <- var(rand)  
  print(rand_var)
  RMSE <- sqrt(mean(residuals(mb0_full)^2))
  print(RMSE)}


fixef(mb0_full)
start_vals <- fixef(mb0_full)

b0 <- start_vals["b0"]
b1 <- start_vals["b1"]
b2 <- start_vals["b2"]
b3 <- start_vals["b3"]
b4 <- start_vals["b4"]

dat$eta0 <- with(dat,b0 +b1 * DBH +b2 * StandHtDomOrCodom +
                   b3 * TreeBALiveTotalPHa +b4 * TreeDensityLiveTotalPHa)

dat$mu_fixed <- with(dat,HtTot * (1 - exp(-eta0)))

dat$z_u <- with(dat,HtTot * exp(-eta0))

vc <- VarCorr(mb0_full)
vc
sigma_u2 <- as.numeric(vc["b0", "Variance"])
sigma_e2 <- as.numeric(vc["Residual",   "Variance"])   
sigma_u2
sigma_e2


u_full <- ranef(mb0_full)[, "b0"]  # named by PlotName

compute_u_hat_for_plot <- function(df_plot_sub, sigma_u2, sigma_e2) {
  r_i <- df_plot_sub$HtLiveBranch - df_plot_sub$mu_fixed
  z_i <- df_plot_sub$z_u
  
  num <- (sigma_u2 / sigma_e2) * sum(z_i * r_i)
  den <- 1 + (sigma_u2 / sigma_e2) * sum(z_i^2)
  
  u_hat <- num / den
  return(u_hat)
}

set.seed(123)
k_values <- 2:10
n_rep    <- 300

plot_sizes <- table(dat$PlotName)

results_list <- list()

for (k in k_values) {
  message("Calibrating with k = ", k, " trees per plot...")
  
  valid_plots <- names(plot_sizes)[plot_sizes >= k]
  dat_k_all   <- dat[dat$PlotName %in% valid_plots, ]
  
  metrics_k <- data.frame(
    RMSE = numeric(0),
    Bias = numeric(0),
    VarBias = numeric(0)
  )
  
  for (rep in 1:n_rep) {
    
    calib_dat <- dat_k_all %>%
      group_by(PlotName) %>%
      slice_sample(n = k) %>%
      ungroup()
    
    u_hat_list <- calib_dat %>%
      group_by(PlotName) %>%
      group_map(~ {
        data.frame(
          PlotName = .y$PlotName,
          u_hat    = compute_u_hat_for_plot(.x, sigma_u2, sigma_e2)
        )
      })
    
    u_hat_df <- do.call(rbind, u_hat_list)
    
    u_hat_df$PlotName <- as.character(u_hat_df$PlotName)
    
    dat_pred <- dat_k_all  # full data for valid plots
    
    dat_pred <- dat_pred %>%
      left_join(u_hat_df, by = "PlotName")
    
    dat_pred <- dat_pred[!is.na(dat_pred$u_hat), ]
    dat_pred$Ht_pred <- with(dat_pred,
                             HtTot * (1 - exp(-(eta0 + u_hat))))
    
    res_y <- dat_pred$HtLiveBranch - dat_pred$Ht_pred
    
    bias_y <- mean(res_y)
    rmse_y <- sqrt(mean(res_y^2))
    
    plot_bias <- tapply(res_y, dat_pred$PlotName, mean)
    var_bias  <- var(plot_bias)
    
    metrics_k <- rbind(
      metrics_k,
      data.frame(RMSE = rmse_y, Bias = bias_y, VarBias = var_bias)
    )
  }
  
  summary_k <- data.frame(
    k           = k,
    RMSE_mean   = mean(metrics_k$RMSE),
    RMSE_sd     = sd(metrics_k$RMSE),
    Bias_mean   = mean(metrics_k$Bias),
    Bias_sd     = sd(metrics_k$Bias),
    VarBias_mean = mean(metrics_k$VarBias),
    VarBias_sd   = sd(metrics_k$VarBias)
  )
  
  results_list[[as.character(k)]] <- summary_k
}

results_df <- do.call(rbind, results_list)
results_df

results_df<-edit(results_df)

names(results_df)
results_df$k
q1<-ggplot(results_df, aes(x = k, y = RMSE_mean)) +
  geom_line(color = "blue") +  geom_point(color = "blue4") +
  geom_errorbar(aes(ymin = RMSE_mean - RMSE_sd, ymax = RMSE_mean + RMSE_sd),
                width = 0.3, color = "blue")  +
  scale_x_continuous(breaks = 2:10) +  labs(x = "",
                                            y = "RSME (m)",title = "a)") +  theme_bw()

q2<- ggplot(results_df, aes(x = k, y = VarBias_mean)) +
  geom_line(color = "blue") +  geom_point(color = "blue4") +
  geom_errorbar(aes(ymin = VarBias_mean - VarBias_sd, ymax = VarBias_mean + VarBias_sd),
                width = 0.3, color = "blue")+ scale_x_continuous(breaks = 2:10) +
  labs(x = "Number of trees per plot",
       y = "Variance of bias",title = "b)") +  theme_bw()
q0<-grid.arrange(q1, q2, ncol = 1,nrow=2)
print(q0)

tiff(filename = "Plot_RSME_n.jpg", width = 100, height = 150,
     units = "mm", res = 400, 
     compression = "lzw") 
q0
grid::grid.draw(q0)
dev.off()

#end

################################# Residual diagnostics across independent variables - natural origin

a<- n
names(a)
str(a)

set_predictors <- list(
  "DBH"             = c("DBH"),
  "DBH+SHT"         = c("DBH", "SHT"),
  "DBH+SHT+BAH"     = c("DBH", "SHT", "BAH"),
  "DBH+SHT+BAH+TPH" = c("DBH", "SHT", "BAH", "TPH"),
  "DBH+BAH+TPH"     = c("DBH", "BAH", "TPH")
)

set_titles <- c(
  "DBH"             = "PC1\n(DBH)",
  "DBH+SHT"         = "PC2\n(DBH + SHT)",
  "DBH+SHT+BAH"     = "PC3\n(DBH + SHT + BAH)",
  "DBH+SHT+BAH+TPH" = "PC4\n(DBH + SHT + BAH + TPH)",
  "DBH+BAH+TPH"     = "PC5\n(DBH + BAH + TPH)"
)

family_order <- c("M1", "M2", "M3")

print(unique(a$Set))

names(a)
dat_long <- a %>%
  dplyr::select(Family, Set, DBH, SHT, BAH, TPH, res_std ) %>%
  tidyr::pivot_longer(
    cols = c(DBH, SHT, BAH, TPH),
    names_to = "Predictor",  values_to = "PredictorValue") %>%
  rowwise() %>%
  filter( Predictor %in% set_predictors[[Set]]) %>%  ungroup()


n_classes <- 10

res_sum <- dat_long %>%  filter( is.finite(PredictorValue),is.finite(res_std)
) %>%
  group_by( Family, Set,Predictor ) %>%  
  mutate(qclass = ntile( PredictorValue,n_classes )) %>%
  group_by(Family,Set,Predictor,qclass) %>%
  summarise(N = n(),
            x_mean =mean( PredictorValue, na.rm = TRUE),
            res_mean =mean( res_std,na.rm = TRUE),
            res_sd =sd( res_std,na.rm = TRUE),
            res_se =res_sd / sqrt(N),
            .groups = "drop")


print(res_sum)

res_sum %>%  count(Family,Set,Predictor) %>%print(n = Inf)

predictor_label <- function(pred) {
  switch( pred,
          "DBH" = "DBH (cm)",
          "SHT" = "SHT (m)",
          "BAH" = expression(
            BAH~(m^2~ha^{-1})),
          "TPH" = expression(
            TPH~(trees~ha^{-1})))}

make_panel <- function(
    data, fam,set_name, pred, show_y = FALSE, show_x = TRUE) {
  
  d <- data %>%
    filter( Family == fam,Set == set_name, Predictor == pred )
  
  ggplot( d,aes(x = x_mean, y = res_mean)) +  
    geom_hline(yintercept = 0,linetype = "dashed",color = "grey75", linewidth = 0.4) +
    
    geom_line(linewidth = 0.35) +
    
    geom_errorbar(aes(ymin = res_mean - res_se, ymax = res_mean + res_se),
                  width = 0,  linewidth = 0.3) +
    
    geom_point( size = 1.5, color="#2171b9") +
    
    annotate( "text", x = -Inf, y = Inf,label = pred, hjust = -0.15,vjust = 1.35,
              fontface = "plain",  size = 3.1) +
    labs( x = NULL, y = NULL) +  theme_bw(base_size = 9) +
    theme(panel.grid.minor =element_blank(),
          panel.grid.major =element_line(linewidth = 0.25,color = "grey90"),
          axis.text.x = if (show_x) element_text(size = 7) else element_blank(),
          axis.ticks.x =if (show_x) element_line() else element_blank(),
          axis.text.y = if (show_y) element_text(size = 7) else element_blank(),
          axis.ticks.y =if (show_y) element_line() else element_blank(),
          plot.margin =  ggplot2::margin( t = 2, r = 3, b = 2, l = 3))}


make_predictor_row <- function( set_name, pred) {
  p1 <- make_panel( data = res_sum, fam = "M1",set_name = set_name,pred = pred,show_y = TRUE,show_x = TRUE  )
  p2 <- make_panel(data = res_sum,fam = "M2",set_name = set_name,pred = pred,show_y = FALSE,show_x = TRUE  )
  p3 <- make_panel(data = res_sum,fam = "M3",set_name = set_name,pred = pred,show_y = FALSE,show_x = TRUE  )
  patchwork::wrap_plots( list( p1, p2, p3 ), nrow = 1, widths = c( 1, 1, 1 ))}

make_set_block <- function(set_name) {
  preds <- set_predictors[[set_name]]
  predictor_rows <- lapply( preds, function(pr) {
    make_predictor_row( set_name = set_name,pred = pr)})
  patchwork::wrap_plots(predictor_rows, ncol = 1)}

make_PC_label <- function(label) {
  ggplot() +annotate( "text", x = 0.5, y = 0.5, label = label, angle= 90 ,fontface = "plain", size = 3.6, lineheight = 1.15) +
    xlim(0, 1 ) + ylim( 0, 1 ) + theme_void() +  
    theme(plot.background = element_rect( fill = "white", color = "black", linewidth = 0.5))}

make_complete_block <- function(set_name) {
  label_plot <-make_PC_label(set_titles[[set_name]])
  diagnostic_block <-make_set_block(set_name)
  patchwork::wrap_plots( list(label_plot, diagnostic_block),
                         nrow = 1, widths = c(0.15, 1))}

make_header <- function(label) {ggplot() +  
    annotate( "text", x = 0.5, y = 0.5, label = label, fontface = "plain", size = 3.6) +
    xlim( 0, 1) + ylim( 0,1) +  theme_void() +
    theme( plot.background = element_rect( fill = "white", color = "black", linewidth = 0.5))}

header_blank <- ggplot() +theme_void()
header_M1 <-make_header("M1")
header_M2 <-make_header("M2")
header_M3 <-make_header("M3")

model_headers <-patchwork::wrap_plots(
  list(header_M1, header_M2, header_M3),nrow = 1)

top_header <-  patchwork::wrap_plots( list( header_blank, model_headers ),
                                      nrow = 1, widths = c( 0.15, 1))

PC1 <-make_complete_block("DBH")
PC2 <-make_complete_block("DBH+SHT")
PC3 <-make_complete_block("DBH+SHT+BAH")
PC4 <-make_complete_block("DBH+SHT+BAH+TPH")
PC5 <-make_complete_block("DBH+BAH+TPH")


f.fig <-patchwork::wrap_plots(
  list( top_header, PC1, PC2, PC3, PC4, PC5),
  ncol = 1, heights = c(0.35, 1, 2, 3, 4, 3 )) +
  patchwork::plot_annotation(title ="Stand Origin: Natural",
                             theme = theme(plot.title =element_text(hjust = 0.5, face = "plain", size = 14)))


tiff(filename = "Res_per_Vars_n02.jpg", width = 150, height = 220,
     units = "mm", res = 400, 
     compression = "lzw")

f.fig

dev.off()


#other save
ggsave(filename = "Res_Natural.png",
       plot =f.fig, width = 13, height = 18,  units = "in", dpi = 600, bg = "white")





###=====================================================================================================
###------------------------- Plantation origin Jack pine -----------------------------------------------

# data
p <- read.csv("Plantaion origin jack pine data file")

############################## Parametric modeling ########################

set.seed(2002)

predictor_sets <- list(
  DBH                 = c("DBH"),
  DBH_SHT             = c("DBH", "StandHtDomOrCodom"),
  DBH_SHT_BAH         = c("DBH", "StandHtDomOrCodom", "TreeBALiveTotalPHa"),
  DBH_SHT_BAH_TPH     = c("DBH", "StandHtDomOrCodom", "TreeBALiveTotalPHa", "TreeDensityLiveTotalPHa"),
  DBH_BAH_TPH         = c("DBH", "TreeBALiveTotalPHa", "TreeDensityLiveTotalPHa")
)

cols_needed_all <- unique(c("PlotName","HtLiveBranch","HtTot","CrownClassCode", unlist(predictor_sets)))
p0 <- na.omit(p[, cols_needed_all])

#### Data summaries
names(p0)
summary_p <- p0 %>%
  summarise(
    n_Plots = n_distinct(PlotName),
    n_Trees = n(),
    HTLCB_mean = mean(HtLiveBranch, na.rm = TRUE),
    HTLCB_min  = min(HtLiveBranch, na.rm = TRUE),
    HTLCB_max  = max(HtLiveBranch, na.rm = TRUE),
    HTLCB_sd = sd (HtLiveBranch, na.rm = TRUE),
    DBH_mean = mean(DBH, na.rm = TRUE),
    DBH_min  = min(DBH, na.rm = TRUE),
    DBH_max  = max(DBH, na.rm = TRUE),
    DBH_sd  = sd(DBH, na.rm = TRUE),
    Ht_mean = mean(HtTot, na.rm = TRUE),
    Ht_min  = min(HtTot, na.rm = TRUE),
    Ht_max  = max(HtTot, na.rm = TRUE),
    Ht_sd  = sd(HtTot, na.rm = TRUE),
    SHT_mean = mean(StandHtDomOrCodom, na.rm = TRUE),
    SHT_min  = min(StandHtDomOrCodom, na.rm = TRUE),
    SHT_max  = max(StandHtDomOrCodom, na.rm = TRUE),
    SHT_sd  = sd(StandHtDomOrCodom, na.rm = TRUE),
    BAH_mean = mean(TreeBALiveTotalPHa, na.rm = TRUE),
    BAH_min  = min(TreeBALiveTotalPHa, na.rm = TRUE),
    BAH_max  = max(TreeBALiveTotalPHa, na.rm = TRUE),
    BAH_sd  = sd(TreeBALiveTotalPHa, na.rm = TRUE),
    TPH_mean = mean(TreeDensityLiveTotalPHa, na.rm = TRUE),
    TPH_min  = min(TreeDensityLiveTotalPHa, na.rm = TRUE),
    TPH_max  = max(TreeDensityLiveTotalPHa, na.rm = TRUE),
    TPH_sd  = sd(TreeDensityLiveTotalPHa, na.rm = TRUE)
  )

summary_p
clamp_ratio <- function(df, eps = 1e-8) {
  r <- pmin(1 - eps, pmax(eps, df$HtLiveBranch / pmax(eps, df$HtTot)))
  df$._ratio_ <- r
  df
}

metrics_fullIC <- function(obs, pred, k_params) {
  n   <- length(obs)
  rss <- sum((obs - pred)^2)
  r2   <- cor(obs, pred)^2
  rmse <- sqrt(rss / n)
  bias <- mean(pred - obs)
  neg2loglik <- n * (log(2*pi) + 1 + log(rss / n))
  aic <- neg2loglik + 2 * k_params
  bic <- neg2loglik + log(n) * k_params
  c(R2 = r2, RMSE = rmse, Bias = bias, AIC = aic, BIC = bic)
}

bind_rows_safe <- function(dflist) {
  if (!length(dflist)) return(data.frame())
  all_cols <- unique(unlist(lapply(dflist, names)))
  padded <- lapply(dflist, function(df) {
    miss <- setdiff(all_cols, names(df))
    if (length(miss)) df[miss] <- NA
    df[all_cols]
  })
  do.call(rbind, padded)
}

.link_eta_start <- function(df, mt, vars) {
  eps <- 1e-8
  r <- pmin(1 - eps, pmax(eps, df$HtLiveBranch / pmax(eps, df$HtTot)))
  X <- as.data.frame(df[, vars, drop = FALSE])
  X <- cbind(Intercept = 1, X)
  
  y <- switch(mt,
              M1 = -log(1 - r),            
              M2 = log(r / (1 - r)),        
              M3 = -log(1/r - 1)            
  )
  ols <- try(lm(y ~ . - 1, data = as.data.frame(X)), silent = TRUE)
  if (inherits(ols, "try-error") || any(!is.finite(coef(ols)))) return(NULL)
  
  co <- as.numeric(coef(ols)); names(co) <- colnames(X)
  b0 <- co["Intercept"]; if (is.na(b0)) b0 <- 0
  bj <- co[setdiff(names(co), "Intercept")]
  names(bj) <- paste0("b", seq_along(vars))
  list(b0 = b0, bj = bj)
}

build_formula_M <- function(model_type, vars) {
  et <- paste(c("b0", paste0("b", seq_along(vars), " * ", vars)), collapse = " + ")
  switch(model_type,
         M1 = as.formula(paste0("HtLiveBranch ~ HtTot * (1 - exp(-(", et, ")))")),
         M2 = as.formula(paste0("HtLiveBranch ~ HtTot / (1 + exp(-(", et, ")))")),
         M3 = as.formula(paste0("HtLiveBranch ~ HtTot / ((1 + exp(-(", et, ")))^(1/m))"))
  )
}

fit_M_once <- function(train, test) {
  out <- list()
  for (mt in c("M1","M2","M3")) {
    for (ps_name in names(predictor_sets)) {
      vars <- predictor_sets[[ps_name]]
      model_id <- paste(mt, ps_name, sep = "_")
      fml  <- build_formula_M(mt, vars)
      
      st_lin  <- .link_eta_start(train, mt, vars)
      k       <- length(vars)
      st_base <- c(b0 = -0.05, setNames(rep(0.05, k), paste0("b", seq_len(k))))
      if (!is.null(st_lin)) {
        st_base["b0"] <- st_lin$b0
        for (j in seq_along(vars)) {
          nm <- paste0("b", j)
          if (!is.na(st_lin$bj[nm])) st_base[nm] <- st_lin$bj[nm]
        }
      }
      if (mt == "M3") st_base <- c(st_base, m = 1.0)
      
      start_sets <- list(
        st_base,
        { tmp <- st_base; tmp["b0"] <- tmp["b0"] * 0.5; tmp },
        { tmp <- st_base; tmp[grep("^b\\d+$", names(tmp))] <- 0; tmp },
        { tmp <- st_base; tmp[grep("^b\\d+$", names(tmp))] <- tmp[grep("^b\\d+$", names(tmp))] * 0.5; tmp }
      )
      
      fit <- NULL
      for (st in start_sets) {
        if (mt == "M3" && !("m" %in% names(st))) st <- c(st, m = 1.0)
        fit <- try(
          nlsLM(fml, data = train, start = st,
                control = nls.lm.control(maxiter = 600, ftol = 1e-10, ptol = 1e-10)),
          silent = TRUE
        )
        if (!inherits(fit, "try-error")) break
      }
      if (inherits(fit, "try-error")) {

          row <- data.frame(
          Model = model_id, Type = "Fixed", Predictors = paste(vars, collapse = "+"),
          R2 = NA_real_, RMSE = NA_real_, Bias = NA_real_, AIC = NA_real_, BIC = NA_real_,
          R2_test = NA_real_, RMSE_test = NA_real_, Bias_test = NA_real_,
          b0 = NA_real_, check.names = FALSE
        )
        for (j in seq_along(vars)) row[[paste0("b", j)]] <- NA_real_
        if (mt == "M3")  row[["m"]] <- NA_real_ 
        out[[length(out)+1]] <- row
        next
      }
      
      pred_tr <- as.numeric(predict(fit, newdata = train))
      met_tr  <- c(R2 = cor(train$HtLiveBranch, pred_tr)^2,
                   RMSE = sqrt(mean((train$HtLiveBranch - pred_tr)^2)),
                   Bias = mean(pred_tr - train$HtLiveBranch),
                   AIC = AIC(fit), BIC = BIC(fit))
      

      cf <- coef(fit)
      eta_fun <- function(d) {
        eta <- cf["b0"]
        for (j in seq_along(vars)) eta <- eta + cf[paste0("b", j)] * d[[vars[j]]]
        eta
      }
      pred_te <- switch(mt,
                        M1 = test$HtTot * (1 - exp(-eta_fun(test))),
                        M2 = test$HtTot / (1 + exp(-eta_fun(test))),
                        M3 = {  mpar <- cf["m"]; test$HtTot / ((1 + exp(-eta_fun(test)))^(1/mpar)) }
      )
      met_te <- c(R2_test = cor(test$HtLiveBranch, pred_te)^2,
                  RMSE_test = sqrt(mean((test$HtLiveBranch - pred_te)^2)),
                  Bias_test = mean(pred_te - test$HtLiveBranch))
      
      row <- data.frame(
        Model = model_id, Type = "Fixed", Predictors = paste(vars, collapse = "+"),
        t(met_tr), t(met_te), check.names = FALSE
      )

      param_names <- c("b0", paste0("b", seq_along(vars)), if (mt=="M3") "m")
      for (nm in param_names) row[[nm]] <- if (nm %in% names(cf)) cf[[nm]] else NA_real_
      out[[length(out)+1]] <- row
    }
  }
  bind_rows_safe(out)
}

grouped_kfold_indices <- function(df, K = 10) {
  pl <- sample(unique(df$PlotName))
  split(pl, cut(seq_along(pl), K, labels = FALSE))
}

run_grouped_kfold <- function(df, K = 10) {
  folds <- grouped_kfold_indices(df, K)
  per_fold <- list()
  for (k in seq_along(folds)) {
    test_plots  <- folds[[k]]
    train_plots <- setdiff(unique(df$PlotName), test_plots)
    train <- clamp_ratio(df[df$PlotName %in% train_plots, ])
    test  <- clamp_ratio(df[df$PlotName %in% test_plots, ])
    res   <- fit_M_once(train, test)         # << ONLY M1/M2/M3
    res$Fold <- k
    per_fold[[k]] <- res
  }
  bind_rows_safe(per_fold)
}

# Run CV
K <- 10
cv_results <- run_grouped_kfold(p0, K)

agg <- do.call(rbind, lapply(split(cv_results, list(cv_results$Model, cv_results$Predictors), drop = TRUE), function(d) {
  data.frame(
    Model = d$Model[1],
    Predictors = d$Predictors[1],
    Folds = nrow(d),
    # Train 
    R2_mean = mean(d$R2, na.rm = TRUE),
    RMSE_mean = mean(d$RMSE, na.rm = TRUE),
    Bias_mean = mean(d$Bias, na.rm = TRUE),
    AIC_mean = mean(d$AIC, na.rm = TRUE),
    BIC_mean = mean(d$BIC, na.rm = TRUE),
    # Test
    R2_test_mean   = mean(d$R2_test, na.rm = TRUE),
    R2_test_sd     = sd(d$R2_test, na.rm = TRUE),
    RMSE_test_mean = mean(d$RMSE_test, na.rm = TRUE),
    RMSE_test_sd   = sd(d$RMSE_test, na.rm = TRUE),
    Bias_test_mean = mean(d$Bias_test, na.rm = TRUE),
    Bias_test_sd   = sd(d$Bias_test, na.rm = TRUE),
    Family = sub("_.*$", "", d$Model[1]),
    stringsAsFactors = FALSE
  )
}))
agg <- agg[order(agg$RMSE_test_mean), ]
print(agg)


predict_row_on_data <- function(row, data) {
  fam  <- sub("_.*$", "", row$Model)   # M1/M2/M3
  vars <- strsplit(row$Predictors, "\\+")[[1]]
  
  b0 <- suppressWarnings(as.numeric(row$b0)); if (is.na(b0)) b0 <- 0
  bj <- sapply(seq_along(vars), function(j) suppressWarnings(as.numeric(row[[paste0("b", j)]])))
  bj[is.na(bj)] <- 0
  
  eta <- rep(b0, nrow(data))
  for (j in seq_along(vars)) eta <- eta + bj[j] * data[[vars[j]]]
  
  switch(fam,
         M1 = data$HtTot * (1 - exp(-eta)),
         M2 = data$HtTot / (1 + exp(-eta)),
         M3 = {
           cpar <- suppressWarnings(as.numeric(row$c)); if (is.na(cpar)) cpar <- 0.3
           mpar <- suppressWarnings(as.numeric(row$m)); if (is.na(mpar)) mpar <- 1.0
           data$HtTot / ((1 + cpar * exp(-eta))^(1/mpar))
         },
         NA_real_
  )
}

set.seed(2002)
folds <- grouped_kfold_indices(p0, K)

pred_plot_df <- list()
for (k in seq_along(folds)) {
  test_plots <- folds[[k]]
  test_df    <- p0[p0$PlotName %in% test_plots, , drop = FALSE]
  rows_k <- cv_results[cv_results$Fold == k, , drop = FALSE]
  if (!nrow(rows_k)) next
  
  for (i in seq_len(nrow(rows_k))) {
    rowi <- rows_k[i, ]
    pred <- try(predict_row_on_data(rowi, test_df), silent = TRUE)
    if (inherits(pred, "try-error") || all(!is.finite(pred))) next
    pred_plot_df[[length(pred_plot_df)+1]] <- data.frame(
      Model     = rowi$Model,
      Observed  = test_df$HtLiveBranch,
      Predicted = as.numeric(pred),
      Height = test_df$HtTot,
      CrownClass = test_df$CrownClassCode,
      plotname = test_df$PlotName,
      dbh = test_df$DBH,
      SHT = test_df$StandHtDomOrCodom,
      BAH = test_df$TreeBALiveTotalPHa,
      TPH = test_df$TreeDensityLiveTotalPHa)}}
pred_plot_df <- if (length(pred_plot_df)) do.call(rbind, pred_plot_df) else data.frame()

#### Split Model into Family and Set _ PLOTS ####
pred_plot_p <- pred_plot_df %>%
  mutate(
    Family = sub("_.*$", "", Model),             # M1/M2/M3
    Set    = sub("^[^_]+_", "", Model)
  )

family_levels <- c("M1","M2","M3")
set_levels    <- c("DBH","DBH+SHT","DBH+SHT+BAH","DBH+SHT+BAH+TPH","DBH+BAH+TPH")
pred_plot_p <- pred_plot_p %>%
  filter(Family %in% family_levels) %>%
  mutate(Family = factor(Family, levels = family_levels),
         Set    = factor(Set,    levels = set_levels))

rng <- range(c(pred_plot_p$Observed), na.rm = TRUE, finite = TRUE)
pad <- 0.05 * diff(rng); if (!is.finite(pad)) pad <- 0
lims <- c(rng[1] - pad, rng[2] + pad)


names(pred_plot_p)

pred_plot_p <- pred_plot_p %>%
  filter(Predicted >= 0)


stats_df00 <- pred_plot_p %>%  group_by(Family, Set) %>%
  summarise(    R2 = round(cor(Observed, Predicted, use = "complete.obs")^2, 3),
                RMSE = round(sqrt(mean((Observed - Predicted)^2, na.rm = TRUE)), 3),
                .groups = "drop"  )

tiff(filename = "ObsPred_P_2.jpg", width = 240, height = 180,
     units = "mm", res = 400, 
     compression = "lzw") 
ggplot(pred_plot_p, aes(x = Observed, y = Predicted)) +
  geom_bin2d(bins = 35) +  scale_fill_gradientn(
    colours = c("#dadaeb", "#9e9ac8", "#54278f"),
    trans = "sqrt",
    limits = range(1, 75),  
    breaks = c(1, 75),      
    labels = c("1", "75"), 
    name = NULL,             
    guide = guide_colorbar(      direction = "horizontal",
                                 barheight = unit(0.3, "cm"),      barwidth = unit(4, "cm"),
                                 ticks = FALSE,      label.position = "bottom",
                                 title.position = "top"    )  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_smooth(
    method = "lm",    formula = y ~0+ x,    se = TRUE,    color = "red",
    fill = "red",    alpha = 0.10,    linewidth = 0.5  ) +
  facet_grid(Family ~ Set, drop = FALSE) +
  coord_equal(xlim = lims, ylim = lims, expand = FALSE) +
  labs(title = "Plantation origin stand",
       x = "Observed HCB (m)",       y = "Predicted HCB (m)") +
  theme_bw(base_size = 11) +
  theme(    panel.grid.minor = element_blank(),
            strip.text = element_text(face = "bold"),
            panel.spacing = unit(0.6, "lines"),
            legend.position = "bottom",  # move legend below
            legend.justification = "left"  )+
  geom_text(    data = stats_df00,
                aes(x = Inf, y = -Inf,
                    label = paste0("R² = ", R2, "\nRMSE = ", RMSE)    ),
                hjust = 1.1, vjust = -0.5,
                size = 3.2,
                color = "black",
                inherit.aes = FALSE
  )

dev.off()


qp<- pred_plot_p
names(qp)
qp$HTLCBr <- (qp$Predicted *100) / qp$Height
qp$HTLCBrOBS <- (qp$Observed *100) / qp$Height
qp$res <- qp$Observed - qp$Predicted
qp$res_std <- (qp$Observed - qp$Predicted) / sd(qp$Observed - qp$Predicted, na.rm = TRUE)
sd(qp$res_std, na.rm = TRUE)
#res
lim_res <- range(qp$res_std, na.rm = TRUE)

qp$Set <- factor(qp$Set, levels = c(
  "DBH",  "DBH+SHT",  "DBH+SHT+BAH",  "DBH+SHT+BAH+TPH",  "DBH+BAH+TPH"))
tiff(filename = "res_p.jpg", width = 240, height = 180,
     units = "mm", res = 400, 
     compression = "lzw") 


ggplot(qp, aes(x = Predicted, y = res_std)) +
  geom_bin2d(bins = 35) +  scale_fill_gradientn(
    colours = c("#dadaeb", "#9e9ac8", "#54278f"),
    trans = "sqrt",
    limits = range(1, 50),  
    breaks = c(1, 50),      
    labels = c("1", "50"), 
    name = NULL,             
    guide = guide_colorbar(      direction = "horizontal",
                                 barheight = unit(0.3, "cm"),      barwidth = unit(4, "cm"),
                                 ticks = FALSE,      label.position = "bottom",
                                 title.position = "top"    )  ) +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed", color = "black") +
  facet_grid(Family ~ Set, drop = FALSE) +
  labs(title = "Plantation origin stand",
       x = "Predicted HCB (m)",       y = "Standartized residuals (m)") +
  theme_bw(base_size = 11) +
  theme(    panel.grid.minor = element_blank(),
            strip.text = element_text(face = "bold"),
            panel.spacing = unit(0.6, "lines"),
            legend.position = "bottom",  # move legend below
            legend.justification = "left"  )


dev.off()


################### Mixed effect M2 V4 ##########
# Using nlme
names(p)

v5 <- c("DBH","StandHtDomOrCodom","TreeBALiveTotalPHa","TreeDensityLiveTotalPHa")
dat <- na.omit(p[, c("PlotName","HtLiveBranch","HtTot", v5)])
names(dat)

model_formula <- function(b0, b1, b2, b3, b4, DBH, StandHtDomOrCodom, TreeBALiveTotalPHa, TreeDensityLiveTotalPHa, HtTot) {
  HtTot / (1 + exp(-(b0 + b1 * DBH + b2 * StandHtDomOrCodom + b3 * TreeBALiveTotalPHa + b4 * TreeDensityLiveTotalPHa)))}

start_vals <- c(b0 = -0.05, b1 = -0.03, b2 = 0.04, b3 = 0.01, b4 = -0.00007)

mfix <- nls(
  HtLiveBranch ~ HtTot / (1 + exp(-(b0 + b1 * DBH + b2 * StandHtDomOrCodom + b3 * TreeBALiveTotalPHa + b4 * TreeDensityLiveTotalPHa))),
  data = dat,  start = start_vals,  na.action = na.omit)

{sse <- sum(residuals(mfix)^2)
  css <- sum((dat$HtLiveBranch - mean(dat$HtLiveBranch))^2)
  rsq <- 1 - sse / css
  print(rsq)
  error_var <- mfix$sigma^2
  print(error_var)
  coef<-mfix$coefficients
  print(coef)
  RMSE <- sqrt(mean(residuals(mfix)^2))
  print(RMSE)}

mb0 <- nlme(
  HtLiveBranch ~ model_formula(b0, b1, b2, b3, b4, DBH, StandHtDomOrCodom, TreeBALiveTotalPHa, TreeDensityLiveTotalPHa, HtTot),
  data = dat,  fixed = b0 + b1 + b2 + b3 + b4 ~ 1,
  random = b0 ~ 1 | PlotName, 
  start = start_vals,  na.action = na.omit,
  control = nlmeControl(pnlsTol = 1e-6, maxIter = 400),method = "REML")

{sse <- sum(residuals(mb0)^2)
  css <- sum((dat$HtLiveBranch - mean(dat$HtLiveBranch))^2)
  rsq <- 1 - sse / css
  print(rsq)
  error_var <- mb0$sigma^2
  print(error_var)
  coef<-mb0$coefficients
  rand <- coef$random$Plot
  rand_var <- var(rand)  
  print(rand_var)
  RMSE <- sqrt(mean(residuals(mb0)^2))
  print(RMSE)}


summary(mb0)
coefmb0p<-mb0$coefficients
coefmb0p
# Fixed effects
fixed_dfp <- data.frame(term = names(coefmb0p$fixed),
                        estimate = coefmb0p$fixed)
write.csv(fixed_dfp, "coefmb0p_fixed_p.csv", row.names = FALSE)


mb1 <- nlme(
  HtLiveBranch ~ model_formula(b0, b1, b2, b3, b4, DBH, StandHtDomOrCodom, TreeBALiveTotalPHa, TreeDensityLiveTotalPHa, HtTot),
  data = dat,  fixed = b0 + b1 + b2 + b3 + b4 ~ 1,
  random = b1 ~ 1 | PlotName, 
  start = start_vals,  na.action = na.omit,
  control = nlmeControl(pnlsTol = 1e-6, maxIter = 400), method = "REML")

{sse <- sum(residuals(mb1)^2)
  css <- sum((dat$HtLiveBranch - mean(dat$HtLiveBranch))^2)
  rsq <- 1 - sse / css
  print(rsq)
  error_var <- mb1$sigma^2
  print(error_var)
  coef<-mb1$coefficients
  rand <- coef$random$Plot
  rand_var <- var(rand)  
  print(rand_var)
  RMSE <- sqrt(mean(residuals(mb1)^2))
  print(RMSE)}

mb2 <- nlme(
  HtLiveBranch ~ model_formula(b0, b1, b2, b3, b4, DBH, StandHtDomOrCodom, TreeBALiveTotalPHa, TreeDensityLiveTotalPHa, HtTot),
  data = dat,  fixed = b0 + b1 + b2 + b3 + b4 ~ 1,
  random = b2 ~ 1 | PlotName, 
  start = start_vals,  na.action = na.omit,
  control = nlmeControl(pnlsTol = 1e-6, maxIter = 400))

{sse <- sum(residuals(mb2)^2)
  css <- sum((dat$HtLiveBranch - mean(dat$HtLiveBranch))^2)
  rsq <- 1 - sse / css
  print(rsq)
  error_var <- mb2$sigma^2
  print(error_var)
  coef<-mb2$coefficients
  rand <- coef$random$Plot
  rand_var <- var(rand)  
  print(rand_var)
  RMSE <- sqrt(mean(residuals(mb2)^2))
  print(RMSE)}

mb3 <- nlme(
  HtLiveBranch ~ model_formula(b0, b1, b2, b3, b4, DBH, StandHtDomOrCodom, TreeBALiveTotalPHa, TreeDensityLiveTotalPHa, HtTot),
  data = dat,  fixed = b0 + b1 + b2 + b3 + b4 ~ 1,
  random =  b3 ~ 1 | PlotName, 
  start = start_vals,  na.action = na.omit,
  control = nlmeControl(pnlsTol = 1e-6, maxIter = 400))

{sse <- sum(residuals(mb3)^2)
  css <- sum((dat$HtLiveBranch - mean(dat$HtLiveBranch))^2)
  rsq <- 1 - sse / css
  print(rsq)
  error_var <- mb3$sigma^2
  print(error_var)
  coef<-mb3$coefficients
  rand <- coef$random$Plot
  rand_var <- var(rand)  
  print(rand_var)
  RMSE <- sqrt(mean(residuals(mb3)^2))
  print(RMSE)}

mb4 <- nlme(
  HtLiveBranch ~ model_formula(b0, b1, b2, b3, b4, DBH, StandHtDomOrCodom, TreeBALiveTotalPHa, TreeDensityLiveTotalPHa, HtTot),
  data = dat,  fixed = b0 + b1 + b2 + b3 + b4 ~ 1,
  random =  b4 ~ 1 | PlotName, 
  start = start_vals,  na.action = na.omit,
  control = nlmeControl(pnlsTol = 1e-6, maxIter = 400))

{sse <- sum(residuals(mb4)^2)
  css <- sum((dat$HtLiveBranch - mean(dat$HtLiveBranch))^2)
  rsq <- 1 - sse / css
  print(rsq)
  error_var <- mb4$sigma^2
  print(error_var)
  coef<-mb4$coefficients
  rand <- coef$random$Plot
  rand_var <- var(rand)  
  print(rand_var)
  RMSE <- sqrt(mean(residuals(mb4)^2))
  print(RMSE)}

mb2b3 <- nlme(
  HtLiveBranch ~ model_formula(b0, b1, b2, b3, b4, DBH, StandHtDomOrCodom, TreeBALiveTotalPHa, TreeDensityLiveTotalPHa, HtTot),
  data = dat,  fixed = b0 + b1 + b2 + b3 + b4 ~ 1,
  random = b2 + b3 ~ 1 | PlotName, 
  start = start_vals,  na.action = na.omit,
  control = nlmeControl(pnlsTol = 1e-6, maxIter = 400))

{sse <- sum(residuals(mb2b3)^2)
  css <- sum((dat$HtLiveBranch - mean(dat$HtLiveBranch))^2)
  rsq <- 1 - sse / css
  print(rsq)
  error_var <- mb2b3$sigma^2
  print(error_var)
  coef<-mb2b3$coefficients
  rand <- coef$random$Plot
  rand_var <- var(rand)  
  print(rand_var)
  RMSE <- sqrt(mean(residuals(mb2b3)^2))
  print(RMSE)}

mb3b4 <- nlme(
  HtLiveBranch ~ model_formula(b0, b1, b2, b3, b4, DBH, StandHtDomOrCodom, TreeBALiveTotalPHa, TreeDensityLiveTotalPHa, HtTot),
  data = dat,  fixed = b0 + b1 + b2 + b3 + b4 ~ 1,
  random =  b3 + b4 ~ 1 | PlotName, 
  start = start_vals,  na.action = na.omit,
  control = nlmeControl(pnlsTol = 1e-6, maxIter = 400))

{sse <- sum(residuals(mb3b4)^2)
  css <- sum((dat$HtLiveBranch - mean(dat$HtLiveBranch))^2)
  rsq <- 1 - sse / css
  print(rsq)
  error_var <- mb3b4$sigma^2
  print(error_var)
  coef<-mb3b4$coefficients
  rand <- coef$random$Plot
  rand_var <- var(rand)  
  print(rand_var)
  RMSE <- sqrt(mean(residuals(mb3b4)^2))
  print(RMSE)}

mb2b4 <- nlme(
  HtLiveBranch ~ model_formula(b0, b1, b2, b3, b4, DBH, StandHtDomOrCodom, TreeBALiveTotalPHa, TreeDensityLiveTotalPHa, HtTot),
  data = dat,  fixed = b0 + b1 + b2 + b3 + b4 ~ 1,
  random = b2 + b4 ~ 1 | PlotName, 
  start = start_vals,  na.action = na.omit,
  control = nlmeControl(pnlsTol = 1e-6, maxIter = 400))

{sse <- sum(residuals(mb2b4)^2)
  css <- sum((dat$HtLiveBranch - mean(dat$HtLiveBranch))^2)
  rsq <- 1 - sse / css
  print(rsq)
  error_var <- mb2b4$sigma^2
  print(error_var)
  coef<-mb2b4$coefficients
  rand <- coef$random$Plot
  rand_var <- var(rand)  
  print(rand_var)
  RMSE <- sqrt(mean(residuals(mb2b4)^2))
  print(RMSE)}

mb2b3b4 <- nlme(
  HtLiveBranch ~ model_formula(b0, b1, b2, b3, b4, DBH, StandHtDomOrCodom, TreeBALiveTotalPHa, TreeDensityLiveTotalPHa, HtTot),
  data = dat,  fixed = b0 + b1 + b2 + b3 + b4 ~ 1,
  random = b2 + b3 + b4 ~ 1 | PlotName, 
  start = start_vals,  na.action = na.omit,
  control = nlmeControl(pnlsTol = 1e-6, maxIter = 1000))

{sse <- sum(residuals(mb2b3b4)^2)
  css <- sum((dat$HtLiveBranch - mean(dat$HtLiveBranch))^2)
  rsq <- 1 - sse / css
  print(rsq)
  error_var <- mb2b3b4$sigma^2
  print(error_var)
  coef<-mb2b3b4$coefficients
  rand <- coef$random$Plot
  rand_var <- var(rand)  
  print(rand_var)
  RMSE <- sqrt(mean(residuals(mb2b3b4)^2))
  print(RMSE)}

coef<-mb2b3b4$coefficients
print(coef)

anova(mb0, mb1)#, mb2, mb3, mb4, mb2b3, mb3b4, mb2b4, mb2b3b4)
AIC(mfix, mb0)#, mb1,mb2, mb3, mb4, mb2b3, mb3b4, mb2b4, mb2b3b4)
BIC(mfix, mb0, mb1)#,mb2, mb3, mb4, mb2b3, mb3b4, mb2b4, mb2b3b4)

# List of models
models <- list(mb0, mb1, mb2, mb3, mb4, mb2b3, mb3b4, mb2b4, mb2b3b4)
model_names <- c("mb0", "mb1", "mb2", "mb4", "mb3", "mb2b3", "mb3b4", "mb2b4", "mb2b3b4")

null_logLik <- as.numeric(logLik(mfix))

get_pseudo_r2 <- function(model) {
  tryCatch({    ll <- as.numeric(logLik(model))
  r2 <- 1 - (ll / null_logLik)
  return(r2)  }, error = function(e) {    return(NA)  })}

pseudo_r2_values <- sapply(models, get_pseudo_r2)
names(pseudo_r2_values) <- model_names

print(pseudo_r2_values)

#End
######################## Random effects - trees per plot ##############

dat <- na.omit(dat)
dat$PlotName <- as.factor(dat$PlotName)
dat$RowID    <- seq_len(nrow(dat))

# Your model:
model_formula <- function(b0, b1, b2, b3, b4,
                          DBH, StandHtDomOrCodom,
                          TreeBALiveTotalPHa, TreeDensityLiveTotalPHa, HtTot) {
  HtTot / (1 + exp(-(b0 + b1 * DBH + b2 * StandHtDomOrCodom + 
                       b3 * TreeBALiveTotalPHa + b4 * TreeDensityLiveTotalPHa)))}

start_vals <- c(b0 = -0.05, b1 = -0.03, b2 = 0.04, b3 = 0.01, b4 = -0.00007)

mb0_full <- nlme(
  HtLiveBranch ~ model_formula(b0, b1, b2, b3, b4,
                               DBH, StandHtDomOrCodom,TreeBALiveTotalPHa, TreeDensityLiveTotalPHa, HtTot),
  data   = dat,
  fixed  = b0 + b1 + b2 + b3 + b4 ~ 1,
  random = b0 ~ 1 | PlotName,           # random intercept
  start  = start_vals,
  control = nlmeControl(
    pnlsTol= 1e-6,msMaxIter= 200,maxIter= 200,pnlsMaxIter = 20), method = "REML")



summary(mb0_full)

{sse <- sum(residuals(mb0_full)^2)
  css <- sum((dat$HtLiveBranch - mean(dat$HtLiveBranch))^2)
  rsq <- 1 - sse / css
  print(rsq)
  error_var <- mb0_full$sigma^2
  print(error_var)
  coef<-mb0_full$coefficients
  rand <- coef$random$Plot
  rand_var <- var(rand)  
  print(rand_var)
  RMSE <- sqrt(mean(residuals(mb0_full)^2))
  print(RMSE)}


fixef(mb0_full)
start_vals <- fixef(mb0_full)

b0 <- start_vals["b0"]
b1 <- start_vals["b1"]
b2 <- start_vals["b2"]
b3 <- start_vals["b3"]
b4 <- start_vals["b4"]

dat$eta0 <- with(dat,b0 +b1 * DBH +b2 * StandHtDomOrCodom +
                   b3 * TreeBALiveTotalPHa +b4 * TreeDensityLiveTotalPHa)

dat$mu_fixed <- with(dat,HtTot / (1 + exp(-eta0)))

dat$z_u <- with(dat,HtTot / (1 + exp(-eta0)))

vc <- VarCorr(mb0_full)
vc
sigma_u2 <- as.numeric(vc["b0", "Variance"])
sigma_e2 <- as.numeric(vc["Residual",   "Variance"])   # residual variance
sigma_u2
sigma_e2

u_full <- ranef(mb0_full)[, "b0"]  # named by PlotName

compute_u_hat_for_plot <- function(df_plot_sub, sigma_u2, sigma_e2) {
  r_i <- df_plot_sub$HtLiveBranch - df_plot_sub$mu_fixed
  z_i <- df_plot_sub$z_u
  
  num <- (sigma_u2 / sigma_e2) * sum(z_i * r_i)
  den <- 1 + (sigma_u2 / sigma_e2) * sum(z_i^2)
  
  u_hat <- num / den
  return(u_hat)
}

set.seed(123)
k_values <- 2:10
n_rep    <- 300

plot_sizes <- table(dat$PlotName)

results_list <- list()

for (k in k_values) {
  message("Calibrating with k = ", k, " trees per plot...")
  
  valid_plots <- names(plot_sizes)[plot_sizes >= k]
  dat_k_all   <- dat[dat$PlotName %in% valid_plots, ]
  
  metrics_k <- data.frame(
    RMSE = numeric(0),
    Bias = numeric(0),
    VarBias = numeric(0))
  
  for (rep in 1:n_rep) {
    calib_dat <- dat_k_all %>%
      group_by(PlotName) %>%
      slice_sample(n = k) %>%
      ungroup()
    
    u_hat_list <- calib_dat %>%
      group_by(PlotName) %>%
      group_map(~ {
        data.frame(
          PlotName = .y$PlotName,
          u_hat    = compute_u_hat_for_plot(.x, sigma_u2, sigma_e2))})
    
    u_hat_df <- do.call(rbind, u_hat_list)
    
    u_hat_df$PlotName <- as.character(u_hat_df$PlotName)
    
    dat_pred <- dat_k_all  
    
    dat_pred <- dat_pred %>%
      left_join(u_hat_df, by = "PlotName")
    
    dat_pred <- dat_pred[!is.na(dat_pred$u_hat), ]
    
    dat_pred$Ht_pred <- with(dat_pred,
                             HtTot / (1 + exp(-(eta0 + u_hat))))
    
    res_y <- dat_pred$HtLiveBranch - dat_pred$Ht_pred
    
    bias_y <- mean(res_y)
    rmse_y <- sqrt(mean(res_y^2))
    
    plot_bias <- tapply(res_y, dat_pred$PlotName, mean)
    var_bias  <- var(plot_bias)
    
    metrics_k <- rbind(
      metrics_k,
      data.frame(RMSE = rmse_y, Bias = bias_y, VarBias = var_bias)
    )
  }
  
  summary_k <- data.frame(
    k           = k,
    RMSE_mean   = mean(metrics_k$RMSE),
    RMSE_sd     = sd(metrics_k$RMSE),
    Bias_mean   = mean(metrics_k$Bias),
    Bias_sd     = sd(metrics_k$Bias),
    VarBias_mean = mean(metrics_k$VarBias),
    VarBias_sd   = sd(metrics_k$VarBias)
  )
  
  results_list[[as.character(k)]] <- summary_k
}

results_df <- do.call(rbind, results_list)
results_df
results_df<-edit(results_df)


names(results_df)
results_df$k
q1<-ggplot(results_df, aes(x = k, y = RMSE_mean)) +
  geom_line(color = "purple") +  geom_point(color = "purple4") +
  geom_errorbar(aes(ymin = RMSE_mean - RMSE_sd, ymax = RMSE_mean + RMSE_sd),
                width = 0.3, color = "purple")  +
  scale_x_continuous(breaks = 2:10) +  labs(x = "",
                                            y = "RSME (m)",title = "a)") +  theme_bw()

q2<- ggplot(results_df, aes(x = k, y = VarBias_mean)) +
  geom_line(color = "purple") +  geom_point(color = "purple4") +
  geom_errorbar(aes(ymin = VarBias_mean - VarBias_sd, ymax = VarBias_mean + VarBias_sd),
                width = 0.3, color = "purple")+ scale_x_continuous(breaks = 2:10) +
  labs(x = "Number of trees per plot",
       y = "Variance of bias",title = "b)") +  theme_bw()
q0<-grid.arrange(q1, q2, ncol = 1,nrow=2)
q0

tiff(filename = "Plot_RSME_p.jpg", width = 100, height = 150,
     units = "mm", res = 400, 
     compression = "lzw") 
q0
dev.off()

#end
################################# Residual diagnostics across independent variables - natural origin

a<- p
names(a)
str(a)

set_predictors <- list(
  "DBH"             = c("DBH"),
  "DBH+SHT"         = c("DBH", "SHT"),
  "DBH+SHT+BAH"     = c("DBH", "SHT", "BAH"),
  "DBH+SHT+BAH+TPH" = c("DBH", "SHT", "BAH", "TPH"),
  "DBH+BAH+TPH"     = c("DBH", "BAH", "TPH")
)

set_titles <- c(
  "DBH"             = "PC1\n(DBH)",
  "DBH+SHT"         = "PC2\n(DBH + SHT)",
  "DBH+SHT+BAH"     = "PC3\n(DBH + SHT + BAH)",
  "DBH+SHT+BAH+TPH" = "PC4\n(DBH + SHT + BAH + TPH)",
  "DBH+BAH+TPH"     = "PC5\n(DBH + BAH + TPH)"
)

family_order <- c("M1", "M2", "M3")

print(unique(a$Set))

names(a)
dat_long <- a %>%
  dplyr::select(Family, Set, DBH, SHT, BAH, TPH, res_std ) %>%
  tidyr::pivot_longer(
    cols = c(DBH, SHT, BAH, TPH),
    names_to = "Predictor",  values_to = "PredictorValue") %>%
  rowwise() %>%
  filter( Predictor %in% set_predictors[[Set]]) %>%  ungroup()


n_classes <- 10

res_sum <- dat_long %>%  filter( is.finite(PredictorValue),is.finite(res_std)
) %>%
  group_by( Family, Set,Predictor ) %>%  
  mutate(qclass = ntile( PredictorValue,n_classes )) %>%
  group_by(Family,Set,Predictor,qclass) %>%
  summarise(N = n(),
            x_mean =mean( PredictorValue, na.rm = TRUE),
            res_mean =mean( res_std,na.rm = TRUE),
            res_sd =sd( res_std,na.rm = TRUE),
            res_se =res_sd / sqrt(N),
            .groups = "drop")


print(res_sum)

res_sum %>%  count(Family,Set,Predictor) %>%print(n = Inf)

predictor_label <- function(pred) {
  switch( pred,
          "DBH" = "DBH (cm)",
          "SHT" = "SHT (m)",
          "BAH" = expression(
            BAH~(m^2~ha^{-1})),
          "TPH" = expression(
            TPH~(trees~ha^{-1})))}

make_panel <- function(
    data, fam,set_name, pred, show_y = FALSE, show_x = TRUE) {
  
  d <- data %>%
    filter( Family == fam,Set == set_name, Predictor == pred )
  
  ggplot( d,aes(x = x_mean, y = res_mean)) +  
    geom_hline(yintercept = 0,linetype = "dashed",color = "grey75", linewidth = 0.4) +
    
    geom_line(linewidth = 0.35) +
    
    geom_errorbar(aes(ymin = res_mean - res_se, ymax = res_mean + res_se),
                  width = 0,  linewidth = 0.3) +
    
    geom_point( size = 1.5, color="#54278f") +
    
    annotate( "text", x = -Inf, y = Inf,label = pred, hjust = -0.15,vjust = 1.35,
              fontface = "plain",  size = 3.1) +
    labs( x = NULL, y = NULL) +  theme_bw(base_size = 9) +
    theme(panel.grid.minor =element_blank(),
          panel.grid.major =element_line(linewidth = 0.25,color = "grey90"),
          axis.text.x = if (show_x) element_text(size = 7) else element_blank(),
          axis.ticks.x =if (show_x) element_line() else element_blank(),
          axis.text.y = if (show_y) element_text(size = 7) else element_blank(),
          axis.ticks.y =if (show_y) element_line() else element_blank(),
          plot.margin =  ggplot2::margin( t = 2, r = 3, b = 2, l = 3))}


make_predictor_row <- function( set_name, pred) {
  p1 <- make_panel( data = res_sum, fam = "M1",set_name = set_name,pred = pred,show_y = TRUE,show_x = TRUE  )
  p2 <- make_panel(data = res_sum,fam = "M2",set_name = set_name,pred = pred,show_y = FALSE,show_x = TRUE  )
  p3 <- make_panel(data = res_sum,fam = "M3",set_name = set_name,pred = pred,show_y = FALSE,show_x = TRUE  )
  patchwork::wrap_plots( list( p1, p2, p3 ), nrow = 1, widths = c( 1, 1, 1 ))}

make_set_block <- function(set_name) {
  preds <- set_predictors[[set_name]]
  predictor_rows <- lapply( preds, function(pr) {
    make_predictor_row( set_name = set_name,pred = pr)})
  patchwork::wrap_plots(predictor_rows, ncol = 1)}

make_PC_label <- function(label) {
  ggplot() +annotate( "text", x = 0.5, y = 0.5, label = label, angle= 90 ,fontface = "plain", size = 3.6, lineheight = 1.15) +
    xlim(0, 1 ) + ylim( 0, 1 ) + theme_void() +  
    theme(plot.background = element_rect( fill = "white", color = "black", linewidth = 0.5))}

make_complete_block <- function(set_name) {
  label_plot <-make_PC_label(set_titles[[set_name]])
  diagnostic_block <-make_set_block(set_name)
  patchwork::wrap_plots( list(label_plot, diagnostic_block),
                         nrow = 1, widths = c(0.15, 1))}

make_header <- function(label) {ggplot() +  
    annotate( "text", x = 0.5, y = 0.5, label = label, fontface = "plain", size = 3.6) +
    xlim( 0, 1) + ylim( 0,1) +  theme_void() +
    theme( plot.background = element_rect( fill = "white", color = "black", linewidth = 0.5))}

header_blank <- ggplot() +theme_void()
header_M1 <-make_header("M1")
header_M2 <-make_header("M2")
header_M3 <-make_header("M3")

model_headers <-patchwork::wrap_plots(
  list(header_M1, header_M2, header_M3),nrow = 1)

top_header <-  patchwork::wrap_plots( list( header_blank, model_headers ),
                                      nrow = 1, widths = c( 0.15, 1))

PC1 <-make_complete_block("DBH")
PC2 <-make_complete_block("DBH+SHT")
PC3 <-make_complete_block("DBH+SHT+BAH")
PC4 <-make_complete_block("DBH+SHT+BAH+TPH")
PC5 <-make_complete_block("DBH+BAH+TPH")


f.fig <-patchwork::wrap_plots(
  list( top_header, PC1, PC2, PC3, PC4, PC5),
  ncol = 1, heights = c(0.35, 1, 2, 3, 4, 3 )) +
  patchwork::plot_annotation(title ="Stand Origin: Plantation",
                             theme = theme(plot.title =element_text(hjust = 0.5, face = "plain", size = 14)))


tiff(filename = "Res_per_Vars_p02.jpg", width = 150, height = 220,
     units = "mm", res = 400, 
     compression = "lzw")

f.fig

dev.off()


#other save
ggsave(filename = "Res_Natural.png",
       plot =f.fig, width = 13, height = 18,  units = "in", dpi = 600, bg = "white")


##########  Plot both origins ############
datn <-  na.omit(n[, c("PlotName","HtLiveBranch","HtTot", "DBH","StandHtDomOrCodom","TreeBALiveTotalPHa","TreeDensityLiveTotalPHa")])
datp <-  na.omit(p[, c("PlotName","HtLiveBranch","HtTot", "DBH","StandHtDomOrCodom","TreeBALiveTotalPHa","TreeDensityLiveTotalPHa")])

names(datn)
names(datp)
dim(datp)
model_formulan <- function(b0, b1, b2, b3, b4,DBH, StandHtDomOrCodom,
                           TreeBALiveTotalPHa, TreeDensityLiveTotalPHa, HtTot) {
  HtTot * (1 - exp(-(b0 +b1 * DBH +b2 * StandHtDomOrCodom +
                       b3 * TreeBALiveTotalPHa +b4 * TreeDensityLiveTotalPHa)))}
start_vals <- c(b0 = -1,b1 = -0.03, b2 = 0.07, b3 = 0.02, b4 = -0.00003)
mb0n <- nlme(
  HtLiveBranch ~ model_formulan(b0, b1, b2, b3, b4, DBH, StandHtDomOrCodom, TreeBALiveTotalPHa, TreeDensityLiveTotalPHa, HtTot),
  data = datn,  fixed = b0 + b1 + b2 + b3 + b4 ~ 1,
  random = b0 ~ 1 | PlotName, 
  start = start_vals,  na.action = na.omit,
  control = nlmeControl(pnlsTol = 1e-6, maxIter = 400))
## Just Coaf: 
summary(mb0n)
coefn<-mb0n$coefficients
print(coefn)

model_formulap <- function(b0, b1, b2, b3, b4, DBH, StandHtDomOrCodom, TreeBALiveTotalPHa, TreeDensityLiveTotalPHa, HtTot) {
  HtTot / (1 + exp(-(b0 + b1 * DBH + b2 * StandHtDomOrCodom + b3 * TreeBALiveTotalPHa + b4 * TreeDensityLiveTotalPHa)))}
start_vals <- c(b0 = -0.05, b1 = -0.03, b2 = 0.04, b3 = 0.01, b4 = -0.00007)
mb0p <- nlme(
  HtLiveBranch ~ model_formulap(b0, b1, b2, b3, b4, DBH, StandHtDomOrCodom, TreeBALiveTotalPHa, TreeDensityLiveTotalPHa, HtTot),
  data = datp,  fixed = b0 + b1 + b2 + b3 + b4 ~ 1,
  random = b0 ~ 1 | PlotName, 
  start = start_vals,  na.action = na.omit,
  control = nlmeControl(pnlsTol = 1e-6, maxIter = 400))
## Just Coaf: 
summary(mb0p)
coefp<-mb0p$coefficients
print(coefp)

datn$Pred <- predict(mb0n, level = 0)
datp$Pred <- predict(mb0p, level = 0)


#  n
HtTot_seq <- seq(5, 40, length.out = 100)
dbh_mean <- mean(datn$DBH, na.rm = TRUE)
sht_mean <- mean(datn$StandHtDomOrCodom, na.rm = TRUE)
bah_mean <- mean(datn$TreeBALiveTotalPHa, na.rm = TRUE)
tph_mean <- 2000 #mean(datn$TreeDensityLiveTotalPHa, na.rm = TRUE)
coefixn <- fixed.effects(mb0n)
print(coefixn)
Plot<- datn$PlotName
uPlot <- as.data.frame(coefn$random$PlotName)
coefn$random$PlotName
str(uPlot)
# f part
model_lfn <- function(HtTot, coefixn) {
  HtTot * (1 - exp(-(coefixn["b0"] + coefixn["b1"] * dbh_mean +coefixn["b2"] * sht_mean +
                       coefixn["b3"] * bah_mean +coefixn["b4"] * tph_mean)))}

pred_f <- sapply(HtTot_seq, function(ht) {
  model_lfn(ht, coefixn)})
fix.linen <- data.frame(HtTot = HtTot_seq, Pred = pred_f)
# m part (b0 + u0)
model_len <- function(HtTot, coefixn, rand_b0) {
  HtTot * (1 - exp(-( (coefixn["b0"] + rand_b0) +  # b0 + u0
                        coefixn["b1"] * dbh_mean +coefixn["b2"] * sht_mean +
                        coefixn["b3"] * bah_mean +coefixn["b4"] * tph_mean )))}

plot_linesn <- list()

for (plot_name in rownames(uPlot)) {
  u0 <- uPlot[plot_name, "b0"]  # random intercept for this plot
  pred_values <- sapply(HtTot_seq, function(ht) {
    model_len(ht, coefixn, u0)  })
  plot_linesn[[plot_name]] <- data.frame(HtTot = HtTot_seq, Pred = pred_values)}

all_linesn <- bind_rows(plot_linesn, .id = "PlotName")
names(all_linesn)
names(fix.linen)

all_linesn <- all_linesn %>%
  filter(PlotName != "1211696PIP")


nfig <- ggplot() +
  geom_line(data = all_linesn, aes(x = HtTot, y = Pred, group = PlotName) ,alpha = 0.3, linewidth = 0.3, color = "gray") + 
  geom_line(data = fix.linen, aes(x= HtTot, y = Pred), alpha = 0.6, linewidth = 0.9 ,color = "blue", linetype = "dashed")+
  labs(title = "",x = "HT (m)",y = "Predicted HCB (m)") + theme_bw() #a) Natural origin stand
nfig
#  p
mean(datn$TreeDensityLiveTotalPHa)
median(datn$TreeDensityLiveTotalPHa)
HtTot_seq <- seq(5, 40, length.out = 100)
dbh_mean <- mean(datp$DBH, na.rm = TRUE)
sht_mean <- mean(datp$StandHtDomOrCodom, na.rm = TRUE)
bah_mean <- mean(datp$TreeBALiveTotalPHa, na.rm = TRUE)
tph_mean <- 2000 #mean(datp$TreeDensityLiveTotalPHa, na.rm = TRUE)
coefixp <- fixed.effects(mb0p)
print(coefixp)
Plot<- datp$PlotName
uPlot <- as.data.frame(coefp$random$PlotName)
coefp$random$PlotName
str(uPlot)
# f part
model_lfp <- function(HtTot, coefixp) {
  HtTot / (1 + exp(-(coefixp["b0"] + coefixp["b1"] * dbh_mean +coefixp["b2"] * sht_mean +
                       coefixp["b3"] * bah_mean +coefixp["b4"] * tph_mean)))}

pred_f <- sapply(HtTot_seq, function(ht) {
  model_lfp(ht, coefixp)})
fix.linep <- data.frame(HtTot = HtTot_seq, Pred = pred_f)
# m part (b0 + u0)
model_lep <- function(HtTot, coefixp, rand_b0) {
  HtTot / (1 + exp(-( (coefixp["b0"] + rand_b0) +  # b0 + u0
                        coefixp["b1"] * dbh_mean +coefixp["b2"] * sht_mean +
                        coefixp["b3"] * bah_mean +coefixp["b4"] * tph_mean )))}

plot_linesp <- list()

for (plot_name in rownames(uPlot)) {
  u0 <- uPlot[plot_name, "b0"]  # random intercept for this plot
  pred_values <- sapply(HtTot_seq, function(ht) {
    model_lep(ht, coefixp, u0)  })
  plot_linesp[[plot_name]] <- data.frame(HtTot = HtTot_seq, Pred = pred_values)}

all_linesp <- bind_rows(plot_linesp, .id = "PlotName")
names(all_linesp)
names(fix.linep)

pfig <- ggplot() +
  geom_line(data = all_linesp, aes(x = HtTot, y = Pred, group = PlotName) ,alpha = 0.3, linewidth = 0.3 ,color = "gray") + 
  geom_line(data = fix.linep, aes(x= HtTot, y = Pred), alpha = 0.6, linewidth = 0.9 ,color = "purple", linetype = "dashed")+
  labs(title = "",x = "HT (m)",y = "Predicted HCB (m)") + theme_bw()

pfig
Z0<-grid.arrange(nfig, pfig, ncol = 1,nrow=2)
Z0

tiff(filename = "Plot_lines_p-t.jpg", width = 90, height = 90,
     units = "mm", res = 400, 
     compression = "lzw") 
pfig

dev.off()



####### Per BAH classes 
BAH_classes <- seq(10, 50, by = 10)
HtTot_seq <- seq(5, 40, length.out = 200)

StandHt_mean <- 20
DBH_mean <- 20
Density_mean <- 5000

coef_plant <- fixed.effects(mb0p)
coef_natural <- fixed.effects(mb0n)

model_linen <- function(HtTot, BAH_classes, coefs) {
  HtTot * (1 - exp(-(coefs["b0"] + coefs["b1"] * DBH_mean +coefs["b2"] * StandHt_mean +
                       coefs["b3"] * BAH_classes +coefs["b4"] * Density_mean)))}

model_linep <- function(HtTot, BAH_classes, coefs) {
  HtTot / (1 + exp(-(coefs["b0"] + coefs["b1"] * DBH_mean +coefs["b2"] * StandHt_mean +
                       coefs["b3"] * BAH_classes +coefs["b4"] * Density_mean)))}

df_lines <- expand.grid(HtTot = HtTot_seq, BAH = BAH_classes) %>%
  mutate(    Plantation = model_linep(HtTot, BAH, coef_plant),
             Natural = model_linen(HtTot, BAH, coef_natural)  ) %>%
  pivot_longer(cols = c("Plantation", "Natural"), names_to = "Origin", values_to = "HtLiveBranch")

# Plot
df_lines$BAH <-as.factor(df_lines$BAH)

w3<- ggplot(df_lines, aes(x = HtTot, y = HtLiveBranch, color = BAH)) +
  geom_line(size = 0.6) +  facet_wrap(~ Origin, ncol = 2) +
  labs(title = "a) TPH = 5000",x = "Total Height (m)",#
       y = "Predicted HTLCB (m)", color = "BAH (m²/ha)") +
  theme_bw(base_size = 14) +
  theme(    legend.position = "bottom",
            legend.box = "horizontal",
            strip.text = element_text(face = "bold"),
            panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold", hjust = 0.5))+
  #annotate("text", x = 12, y = 35, label = "a) TPH = 1000", color = "black", size = 3)+
  theme_bw()

w3

tiff(filename = "df_lines_BAH.jpg", width = 120, height = 240,
     units = "mm", res = 400, 
     compression = "lzw") 

w0<-grid.arrange(w1, w2, w3, ncol = 1, nrow=3)

dev.off()


####### Covariate plots #########

names(datn)
names(datp)

n1<- ggplot(datn, aes(x = HtTot, y = HtLiveBranch)) +
  geom_bin2d(bins = 50) +
  scale_fill_gradientn(colours = c( "#c6dbef", "#6baed9", "#2171b9"),
                       name = "Frequency",trans = "sqrt") +
  labs(title = "(a)",x = "total height (m)",y = "HTLCB (m)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.spacing = unit(0.6, "lines")) #+ theme(legend.position = "none")
n1

n2<- ggplot(datn, aes(x = DBH, y = HtLiveBranch)) +
  geom_bin2d(bins = 50) +
  scale_fill_gradientn(colours = c( "#c6dbef", "#6baed9", "#2171b9"),
                       name = "Frequency",trans = "sqrt") +
  labs(title = "(b)",x = "DBH (cm)",y = "HTLCB (m)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.spacing = unit(0.6, "lines")) #+ theme(legend.position = "none")
n2

n3<- ggplot(datn, aes(x = StandHtDomOrCodom, y = HtLiveBranch)) +
  geom_bin2d(bins = 50) +
  scale_fill_gradientn(colours = c( "#c6dbef", "#6baed9", "#2171b9"),
                       name = "Frequency",trans = "sqrt") +
  labs(title = "(c)",x = "SHT (m)",y = "HTLCB (m)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.spacing = unit(0.6, "lines")) #+ theme(legend.position = "none")
n3

n4<- ggplot(datn, aes(x = TreeBALiveTotalPHa, y = HtLiveBranch)) +
  geom_bin2d(bins = 50) +
  scale_fill_gradientn(colours = c( "#c6dbef", "#6baed9", "#2171b9"),
                       name = "Frequency",trans = "sqrt") +
  labs(title = "(d)",x = "BAH (m²/ha)",y = "HTLCB (m)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.spacing = unit(0.6, "lines")) #+ theme(legend.position = "none")
n4
n5<- ggplot(datn, aes(x = TreeDensityLiveTotalPHa, y = HtLiveBranch)) +
  geom_bin2d(bins = 50) +
  scale_fill_gradientn(colours = c( "#c6dbef", "#6baed9", "#2171b9"),
                       name = "Frequency",trans = "sqrt") +
  labs(title = "(e)",x = "TPH (n/ha)",y = "HTLCB (m)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.spacing = unit(0.6, "lines")) #+ theme(legend.position = "none")
n5

freqs <- bind_rows(
  ggplot_build(n1)$data[[1]],
  ggplot_build(n2)$data[[1]],
  ggplot_build(n3)$data[[1]],
  ggplot_build(n4)$data[[1]],
  ggplot_build(n5)$data[[1]]
)

max_freq <- max(freqs$count)

fill_scale <- scale_fill_gradientn(
  colours = c("#c6dbef", "#6baed9", "#2171b9"),
  name = "Frequency",
  trans = "sqrt",
  limits = c(0, max_freq),  
  breaks = c(1, 5, 10, 20, 40)  
)

n1 <- ggplot(datn, aes(x = HtTot, y = HtLiveBranch)) +
  geom_bin2d(bins = 50) +
  fill_scale +
  labs(title = "a) Natural origin stand", x = "", y = "HCB (m)") +
  theme_bw(base_size = 12)

n2 <- ggplot(datn, aes(x = DBH, y = HtLiveBranch)) +
  geom_bin2d(bins = 50) +
  fill_scale +
  labs(title = "a) Natural origin stand", x = "", y = "HCB (m)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

n3 <- ggplot(datn, aes(x = StandHtDomOrCodom, y = HtLiveBranch)) +
  geom_bin2d(bins = 50) +
  fill_scale +
  labs(title = "", x = "", y = "") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

n4 <- ggplot(datn, aes(x = TreeBALiveTotalPHa, y = HtLiveBranch)) +
  geom_bin2d(bins = 50) +
  fill_scale +
  labs(title = "", x = "", y = "") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

n5 <- ggplot(datn, aes(x = TreeDensityLiveTotalPHa, y = HtLiveBranch)) +
  geom_bin2d(bins = 50) +
  fill_scale +
  labs(title = "", x = "", y = "") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

pfN <- ( n2 | n3 | n4 | n5) + plot_layout(guides = "collect") & theme(legend.position = "right")
pfN



names(datp)

p1<- ggplot(datp, aes(x = HtTot, y = HtLiveBranch)) +
  geom_bin2d(bins = 50) +
  scale_fill_gradientn(colours = c( "#dadaeb", "#9e9ac8", "#54278f"),
                       name = "Frequency",trans = "sqrt") +
  labs(title = "(a)",x = "total height (m)",y = "HTLCB (m)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.spacing = unit(0.6, "lines")) #+ theme(legend.position = "none")
p1

p2<- ggplot(datp, aes(x = DBH, y = HtLiveBranch)) +
  geom_bin2d(bins = 50) +
  scale_fill_gradientn(colours = c( "#dadaeb", "#9e9ac8", "#54278f"),
                       name = "Frequency",trans = "sqrt") +
  labs(title = "(b)",x = "DBH (cm)",y = "HTLCB (m)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.spacing = unit(0.6, "lines")) #+ theme(legend.position = "none")
p2

p3<- ggplot(datp, aes(x = StandHtDomOrCodom, y = HtLiveBranch)) +
  geom_bin2d(bins = 50) +
  scale_fill_gradientn(colours = c( "#dadaeb", "#9e9ac8", "#54278f"),
                       name = "Frequency",trans = "sqrt") +
  labs(title = "(c)",x = "SHT (m)",y = "HTLCB (m)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.spacing = unit(0.6, "lines")) #+ theme(legend.position = "none")
p3

p4<- ggplot(datp, aes(x = TreeBALiveTotalPHa, y = HtLiveBranch)) +
  geom_bin2d(bins = 50) +
  scale_fill_gradientn(colours = c("#dadaeb", "#9e9ac8", "#54278f"),
                       name = "Frequency",trans = "sqrt") +
  labs(title = "(d)",x = "BAH (m²/ha)",y = "HTLCB (m)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.spacing = unit(0.6, "lines")) #+ theme(legend.position = "none")
p4
p5<- ggplot(datp, aes(x = TreeDensityLiveTotalPHa, y = HtLiveBranch)) +
  geom_bin2d(bins = 50) +
  scale_fill_gradientn(colours = c("#dadaeb", "#9e9ac8", "#54278f"),
                       name = "Frequency",trans = "sqrt") +
  labs(title = "(e)",x = "TPH (n/ha)",y = "HTLCB (m)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.spacing = unit(0.6, "lines")) #+ theme(legend.position = "none")
p5

freqs <- bind_rows(
  ggplot_build(p1)$data[[1]],
  ggplot_build(p2)$data[[1]],
  ggplot_build(p3)$data[[1]],
  ggplot_build(p4)$data[[1]],
  ggplot_build(p5)$data[[1]]
)

max_freqp <- max(freqs$count)

fill_scale <- scale_fill_gradientn(
  colours = c("#dadaeb", "#9e9ac8", "#54278f"),
  name = "Frequency",
  trans = "sqrt",
  limits = c(0, max_freqp),  
  breaks = c(1, 5, 10, 20, 40)  
)

p1 <- ggplot(datp, aes(x = HtTot, y = HtLiveBranch)) +
  geom_bin2d(bins = 50) +
  fill_scale +
  labs(title = "b) Plantation origin stand", x = "Total height (m)", y = "HTLCB (m)") +
  theme_bw(base_size = 12)

p2 <- ggplot(datp, aes(x = DBH, y = HtLiveBranch)) +
  geom_bin2d(bins = 50) +
  fill_scale +
  labs(title = "b) Plantation origin stand", x = "DBH (cm)", y = "HCB (m)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

p3 <- ggplot(datp, aes(x = StandHtDomOrCodom, y = HtLiveBranch)) +
  geom_bin2d(bins = 50) +
  fill_scale +
  labs(title = "", x = "SHT (m)", y = "") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

p4 <- ggplot(datp, aes(x = TreeBALiveTotalPHa, y = HtLiveBranch)) +
  geom_bin2d(bins = 50) +
  fill_scale +
  labs(title = "", x = "BAH (m²/ha)", y = "") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

p5 <- ggplot(datp, aes(x = TreeDensityLiveTotalPHa, y = HtLiveBranch)) +
  geom_bin2d(bins = 50) +
  fill_scale +
  labs(title = "", x = "TPH (n/ha)", y = "") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")


pfP <- (p2 | p3 | p4 | p5) + plot_layout(guides = "collect") & theme(legend.position = "right")
pfP


n0 <- pfN / pfP  
n0

tiff(filename = "Covars_NP.jpg", width = 320, height = 240,
     units = "mm", res = 400, 
     compression = "lzw")
n0
dev.off()


#### Test the origins ############
datn$Origin <- "Natural"
datp$Origin <- "Plantation"

df_all <- rbind(datn, datp)
names(df_all)

df_all$HCB<-df_all$HtLiveBranch
df_all$HT<-df_all$HtTot
df_all$SHT<-df_all$StandHtDomOrCodom
df_all$BAH<-df_all$TreeBALiveTotalPHa
df_all$TPH<-df_all$TreeDensityLiveTotalPHa


vars <- c("HCB", "DBH", "HT",
          "SHT","BAH","TPH")

pvals <- sapply(vars, function(v) {
  wilcox.test(df_all[[v]] ~ df_all$Origin)$p.value
})

round(pvals, 4)


df_long <- df_all %>%
  select(Origin,
         HCB, DBH, HT,
         SHT,BAH,TPH) %>%
  pivot_longer(
    cols = -Origin,
    names_to = "Variable",
    values_to = "Value"  )
Variable_order <- c( "DBH", "HT","HCB","SHT","BAH","TPH")
df_long<- df_long %>% mutate(Variable = factor(Variable, levels = Variable_order))

f0<-ggplot(df_long, aes(x = Origin, y = Value, fill = Origin)) +
  geom_boxplot(width = 0.5, outlier.size = 0, outliers = FALSE) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 3) +
  geom_text(
    aes(x = 2, y = Inf, label = "p < 0.001"),
    inherit.aes = FALSE,vjust = 1.8,size = 2.5) +
  scale_fill_manual(values = c("Natural" = "#2171c9",
                               "Plantation" = "#9e9ac8")) +
  labs(x = "Stand origin",y = NULL) +
  theme_bw(base_size = 10) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank()  )



tiff(filename = "varcopm_NP.jpg", width = 100, height = 150,
     units = "mm", res = 400, 
     compression = "lzw")
f0
dev.off()


