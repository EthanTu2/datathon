
# oda_regression_alternatives.R
# Purpose: Theory-motivated alternative models that reduce multicollinearity
#          when predicting ODA per capita with health-need indicators.
#
# Data: joined_output.csv (place in the same folder as this script)
# Outputs:
#   - oda_alt_summary.txt
#   - oda_alt_results.csv
#   - oda_alt_model_stats.csv
#   - oda_alt_analysis_sample.csv
#   - oda_alt_predictor_correlations.csv

# ---------------------------
# Helpers: cluster-robust vcov
# ---------------------------
vcov_cluster_lm <- function(model, cluster) {
  # model: lm()
  # cluster: vector of cluster IDs aligned with model.frame(model)
  X <- model.matrix(model)
  u <- residuals(model)
  cl <- as.factor(cluster)
  G <- nlevels(cl)
  N <- nrow(X)
  K <- ncol(X)

  # (X'X)^-1
  XtX_inv <- solve(crossprod(X))

  # Meat: sum_g (X_g' u_g)(X_g' u_g)'
  meat <- matrix(0, K, K)
  for (g in levels(cl)) {
    idx <- which(cl == g)
    Xg <- X[idx, , drop = FALSE]
    ug <- u[idx]
    Xu <- crossprod(Xg, ug) # K x 1
    meat <- meat + (Xu %*% t(Xu))
  }

  # Finite-sample correction (common CR1)
  # See e.g. Cameron & Miller (2015)
  scale <- (G / (G - 1)) * ((N - 1) / (N - K))
  V <- scale * (XtX_inv %*% meat %*% XtX_inv)
  return(V)
}

tidy_cluster <- function(model, cluster, model_name) {
  V <- vcov_cluster_lm(model, cluster)
  se <- sqrt(diag(V))
  est <- coef(model)
  tval <- est / se

  cl <- as.factor(cluster)
  G <- nlevels(cl)
  df <- G - 1

  pval <- 2 * pt(abs(tval), df = df, lower.tail = FALSE)
  crit <- qt(0.975, df = df)
  ci_low <- est - crit * se
  ci_high <- est + crit * se

  out <- data.frame(
    model = model_name,
    term = names(est),
    estimate = as.numeric(est),
    se_cluster = as.numeric(se),
    t = as.numeric(tval),
    p = as.numeric(pval),
    ci_low = as.numeric(ci_low),
    ci_high = as.numeric(ci_high),
    stringsAsFactors = FALSE
  )
  return(out)
}

model_stats <- function(model, model_name, n_countries) {
  s <- summary(model)
  data.frame(
    model = model_name,
    n_obs = nobs(model),
    n_countries = n_countries,
    r2 = s$r.squared,
    adj_r2 = s$adj.r.squared,
    stringsAsFactors = FALSE
  )
}

# ---------------------------
# Load data
# ---------------------------
df <- read.csv("joined_output.csv", check.names = FALSE)

# Variables
y_name <- "ODA per capita by recipient"
u5_name <- "mort_rate_per1000_under5"
mmr_name <- "est_maternal_mortality_ratio_per100kbirths"
hiv_name <- "hiv_prev_ages15to49"
tb_name <- "TB_death_rate_per100k"

# Keep only true countries (exclude income-group aggregates with missing code)
df <- df[!is.na(df$Code), ]

# Keep only years with overlap for the 4 predictors
df <- df[df$Year %in% c(2004, 2014), ]

# Keep complete cases for these variables
keep <- complete.cases(df[, c("Entity","Code","Year", y_name, u5_name, mmr_name, hiv_name, tb_name)])
dat <- df[keep, c("Entity","Code","Year", y_name, u5_name, mmr_name, hiv_name, tb_name)]

# Rename and transform
names(dat) <- c("country","code","year","oda_pc","mort_u5","mmr","hiv_prev","tb_death")
dat$year <- as.integer(dat$year)

# Transforms
dat$y_asinh <- asinh(as.numeric(dat$oda_pc))
dat$x_u5  <- log1p(as.numeric(dat$mort_u5))
dat$x_mmr <- log1p(as.numeric(dat$mmr))
dat$x_hiv <- log1p(as.numeric(dat$hiv_prev))
dat$x_tb  <- log1p(as.numeric(dat$tb_death))
dat$year2014 <- ifelse(dat$year == 2014, 1, 0)

# Predictor correlations
corr <- cor(dat[, c("x_u5","x_mmr","x_hiv","x_tb")], use = "complete.obs")

# PCA burden index (standardized predictors)
pca <- prcomp(dat[, c("x_u5","x_mmr","x_hiv","x_tb")], center = TRUE, scale. = TRUE)
dat$burden_pc1 <- pca$x[,1]

# Models (all include a year indicator because we pool 2004 & 2014)
m_full <- lm(y_asinh ~ year2014 + x_u5 + x_mmr + x_hiv + x_tb, data = dat)

# If your goal is a single “need proxy” (less collinearity):
m_u5 <- lm(y_asinh ~ year2014 + x_u5, data = dat)

# If you want to keep *all four* indicators but reduce collinearity:
m_pc1 <- lm(y_asinh ~ year2014 + burden_pc1, data = dat)

# Cluster-robust results
G <- length(unique(dat$code))

res_full <- tidy_cluster(m_full, dat$code, "full_4_predictors")
res_u5   <- tidy_cluster(m_u5, dat$code, "u5_only")
res_pc1  <- tidy_cluster(m_pc1, dat$code, "pca_burden_pc1")
res_all  <- rbind(res_full, res_u5, res_pc1)

stats_all <- rbind(
  model_stats(m_full, "full_4_predictors", G),
  model_stats(m_u5, "u5_only", G),
  model_stats(m_pc1, "pca_burden_pc1", G)
)

# Export files
write.csv(res_all, "oda_alt_results.csv", row.names = FALSE)
write.csv(stats_all, "oda_alt_model_stats.csv", row.names = FALSE)
write.csv(dat[, c("country","code","year","y_asinh","x_u5","x_mmr","x_hiv","x_tb","burden_pc1")],
          "oda_alt_analysis_sample.csv", row.names = FALSE)
write.csv(corr, "oda_alt_predictor_correlations.csv")

# Human-readable summary
sink("oda_alt_summary.txt")
cat("Alternative, theory-motivated models to reduce multicollinearity and improve power\n\n")
cat("DATA NOTE: The four predictors overlap only in 2004 and 2014 in this dataset.\n")
cat("Outcome: y = asinh(ODA per capita by recipient)\n")
cat("Predictors: log1p() transforms\n")
cat(sprintf("Sample: N = %d country-years; Countries (clusters) = %d\n\n", nrow(dat), G))

cat("Predictor correlations (log1p scale):\n")
print(round(corr, 3))
cat("\nPCA burden index:\n")
cat(sprintf("  PC1 variance explained: %.1f%%\n\n", 100 * (pca$sdev[1]^2 / sum(pca$sdev^2))))

cat("MODEL: full_4_predictors (clustered SE by country)\n")
print(res_full)
cat("\nMODEL: u5_only (clustered SE by country)\n")
print(res_u5)
cat("\nMODEL: pca_burden_pc1 (clustered SE by country)\n")
print(res_pc1)

cat("\nModel stats:\n")
print(stats_all)
sink()

cat("Done. Wrote: oda_alt_summary.txt, oda_alt_results.csv, oda_alt_model_stats.csv, oda_alt_analysis_sample.csv, oda_alt_predictor_correlations.csv\n")
