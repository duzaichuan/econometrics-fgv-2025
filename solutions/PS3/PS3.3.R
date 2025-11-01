# --- Packages ---
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2)
  library(lubridate); library(stringr); library(purrr); library(tidyr); library(here)
})

# --- (a) Load and plot IPCA --------------------------------------------------
csv_path <-  here("solutions", "PS3", "ipca_igpm.csv")
ipca_df <- readr::read_csv(csv_path, show_col_types = FALSE)
names(ipca_df)= tolower(names(ipca_df))
ipca_df$date = as.Date(paste0(ipca_df$date,"-01")) #Adjust date
ipca_df = na.omit(ipca_df) #remove NA
Tn <- length(y)

# Plot
g_ipca <- ggplot(ipca_df, aes(date, ipca)) +
  geom_line() +
  labs(title = "IPCA inflation (time series)", x = NULL, y = "IPCA") +
  theme_minimal(base_size = 12)
print(g_ipca)

# --- (b) Conditional log-likelihood for ARMA(p,q) -----------------------------

# Helper: check MA invertibility (all roots of 1 + theta_1 z + ... + theta_q z^q outside unit circle)
ma_invertible <- function(theta, tol = 1e-6) {
  q <- length(theta)
  if (q == 0) return(TRUE)
  rt <- tryCatch(polyroot(c(1, theta)), warning = function(w) NA, error = function(e) NA)
  if (any(is.na(rt))) return(FALSE)
  all(Mod(rt) > 1 + tol)
}

# Negative conditional Gaussian log-likelihood 
ncll_arma <- function(pars, y, p, q) {
  rr <- arma_recursion(y, p, q, pars)         
  m <- max(p, q)
  e_eff <- rr$e[(m + 1):length(y)]
  if (length(e_eff) < 1L || !is.finite(rr$s2) || rr$s2 < 1e-12 || any(!is.finite(e_eff)))
    return(1e12)
  0.5 * sum(log(2 * pi) + log(rr$s2) + (e_eff^2) / rr$s2)
}

# --- (c) Estimation by conditional ML + table of results ----------------------

# Robust CML with OLS starts for AR and invertibility-safe MA handling
fit_arma_cml <- function(y, p, q, start = NULL, loglik_fn = ncll_arma, use_bounds = TRUE) {
  y <- as.numeric(y)
  Tn <- length(y)
  m  <- max(p, q)
  neff <- Tn - m
  if (neff <= 5) stop("Series too short for this (p,q).")
  
  # ---------- OLS starts for AR (phi) ----------
  if (is.null(start)) {
    if (p > 0) {
      M <- embed(y, p + 1L)             # columns: y_t, y_{t-1}, ..., y_{t-p}
      Y <- M[, 1L]
      X <- cbind(1, M[, 2:(p + 1L)])
      ols <- tryCatch(stats::lm.fit(X, Y), error = function(e) NULL)
      
      if (!is.null(ols)) {
        co   <- as.numeric(ols$coefficients)
        c0   <- co[1L]
        phi0 <- co[-1L]
        s20  <- stats::var(drop(Y - X %*% co))
      } else {
        phi0 <- rep(0, p)
        c0   <- mean(y, na.rm = TRUE)
        s20  <- stats::var(y - c0, na.rm = TRUE)
      }
      
      # shrink inside the stationary region if needed
      if (any(Mod(polyroot(c(1, -phi0))) <= 1)) {
        tries <- 0L
        while (any(Mod(polyroot(c(1, -phi0))) <= 1) && tries < 30L) {
          phi0 <- 0.98 * phi0
          tries <- tries + 1L
        }
      }
    } else {
      phi0 <- numeric(0)
      c0   <- mean(y, na.rm = TRUE)
      s20  <- stats::var(y - c0, na.rm = TRUE)
    }
    
    # MA starts: zeros are invertible and stable
    theta0  <- rep(0, q)
    log_s20 <- log(max(s20, 1e-8))
    start   <- c(c0, phi0, theta0, log_s20)
  }
  
  # Names in order
  nm <- c(
    "c",
    if (p > 0) paste0("phi", seq_len(p)) else character(0),
    if (q > 0) paste0("theta", seq_len(q)) else character(0),
    "log_sigma2"
  )
  names(start) <- nm
  
  # ---------- Safe objective: finite + invertibility penalty ----------
  BIG <- 1e12  # moderate "big" to avoid overflow/NaN in finite differencing
  safe_fn <- function(pars) {
    # penalize non-invertible MA to keep the search stable
    if (q > 0) {
      theta <- pars[1L + p + seq_len(q)]
      if (!ma_invertible(theta)) return(BIG)
    }
    val <- tryCatch(loglik_fn(pars, y = y, p = p, q = q),
                    error = function(e) NA_real_,
                    warning = function(w) NA_real_)
    if (!is.finite(val)) BIG else val
  }
  
  # ---------- Bounds (gentle) ----------
  method <- if (use_bounds) "L-BFGS-B" else "BFGS"
  lower <- upper <- NULL
  if (use_bounds) {
    k <- length(start)
    lower <- rep(-Inf, k); upper <- rep( Inf, k)
    
    # keep AR and MA coefficients in reasonable boxes (not sufficient for stationarity/invertibility,
    # but it avoids extreme excursions that cause non-finite evals)
    if (p > 0) {
      idx_phi <- 1L + seq_len(p)
      lower[idx_phi] <- -0.999
      upper[idx_phi] <-  0.999
    }
    if (q > 0) {
      idx_theta <- 1L + p + seq_len(q)
      lower[idx_theta] <- -1.5
      upper[idx_theta] <-  1.5
    }
    # variance bounds
    idx_logsig <- k
    lower[idx_logsig] <- log(1e-8)
    upper[idx_logsig] <- log(1e+8)
  }
  
  opt <- optim(
    par = start, fn = safe_fn, method = method, hessian = TRUE,
    lower = lower, upper = upper,
    control = list(maxit = 5000, reltol = 1e-10,
                   ndeps = rep(1e-5, length(start)))  # a bit larger step helps stability
  )
  if (opt$convergence != 0) {
    warning(sprintf("optim convergence code %d for ARMA(%d,%d)", opt$convergence, p, q))
  }
  
  est <- opt$par; names(est) <- nm
  
  # Variance-covariance, SEs, CIs
  vc <- tryCatch(solve(opt$hessian), error = function(e) NULL)
  se <- if (is.null(vc)) rep(NA_real_, length(est)) else sqrt(diag(vc))
  names(se) <- nm
  ci_low  <- est - 1.96 * se
  ci_high <- est + 1.96 * se
  
  sigma_hat <- sqrt(exp(est["log_sigma2"]))
  se_sigma  <- if (is.null(vc)) NA_real_ else 0.5 * sigma_hat * sqrt(vc["log_sigma2","log_sigma2"])
  sigma_ci  <- c(sigma = sigma_hat,
                 low = if (is.na(se_sigma)) NA_real_ else sigma_hat - 1.96 * se_sigma,
                 high= if (is.na(se_sigma)) NA_real_ else sigma_hat + 1.96 * se_sigma)
  
  ll  <- -safe_fn(est)          # finite by construction
  kpar <- p + q + 2L
  AIC <- 2 * kpar - 2 * ll
  BIC <- kpar * log(neff) - 2 * ll
  
  list(
    spec = sprintf("ARMA(%d,%d)", p, q),
    est = est, se = se, ci_low = ci_low, ci_high = ci_high,
    sigma = sigma_ci, logLik = ll, AIC = AIC, BIC = BIC,
    convergence = opt$convergence, vcov = vc, start = start
  )
}

# Fit the 8 models: (p,q) in {0,1,2} excluding (0,0)
grid <- expand.grid(p = 0:2, q = 0:2) %>% 
  filter(!(p == 0 & q == 0))

y <- ipca_df$ipca

fits <- purrr::pmap(grid, ~ fit_arma_cml(y, ..1, ..2))
 
# Build a compact table of estimates and 95% CIs
summ_tbl <- purrr::map_dfr(fits, function(fit) {
  est <- fit$est
  se  <- fit$se
  ciL <- fit$ci_low
  ciH <- fit$ci_high
  
  # Assemble parameter rows (report sigma instead of log_sigma2)
  pars_to_show <- c(setdiff(names(est), "log_sigma2"), "sigma")
  vals <- c(est[setdiff(names(est), "log_sigma2")], sigma = fit$sigma["sigma"])
  lows <- c(ciL[setdiff(names(est), "log_sigma2")], sigma = fit$sigma["low"])
  highs<- c(ciH[setdiff(names(est), "log_sigma2")], sigma = fit$sigma["high"])
  
  tibble::tibble(
    model = fit$spec,
    param = pars_to_show,
    estimate = as.numeric(vals),
    ci95 = sprintf("[%.4f, %.4f]", as.numeric(lows), as.numeric(highs)),
    logLik = fit$logLik, AIC = fit$AIC, BIC = fit$BIC
  )
}) %>%
  group_by(model) %>%
  mutate(
    logLik = unique(logLik),
    AIC = unique(AIC),
    BIC = unique(BIC)
  ) %>%
  ungroup()

# A more compact, one-row-per-model view (parameters concatenated)
compact_tbl <- summ_tbl %>%
  mutate(term = paste0(param, "=", sprintf("%.3f", estimate), " ", ci95)) %>%
  group_by(model) %>%
  summarise(
    Params = paste(term, collapse = "; "),
    logLik = first(logLik),
    AIC = first(AIC),
    BIC = first(BIC),
    .groups = "drop"
  ) %>%
  arrange(BIC)

# --- (d) Model choice (AIC and BIC) -------------------------------------------

best_AIC <- compact_tbl %>% slice_min(AIC, n = 1)
best_BIC <- compact_tbl %>% slice_min(BIC, n = 1)

cat("\nBest by AIC:\n"); print(best_AIC)
cat("\nBest by BIC:\n"); print(best_BIC)

if (best_AIC$model != best_BIC$model) {
  message("AIC and BIC disagree. Prefer BIC for parsimony unless residual diagnostics suggest misspecification.")
}

# --- (e) Plots: fitted vs actual; residuals vs time ---------------------------

#Simulate the estimated process
arma_recursion <- function(y, p, q, pars) {
  # pars = c(c, phi_1..phi_p, theta_1..theta_q, log_sigma2)
  cst <- pars[1]
  phi <- if (p>0) pars[2:(1+p)] else numeric()
  theta <- if (q>0) pars[(2+p):(1+p+q)] else numeric()
  log_s2 <- pars[2+p+q]; s2 <- exp(log_s2)
  
  e <- numeric(length(y))
  yhat <- numeric(length(y))
  for (t in seq_along(y)) {
    ar <- if (p==0 || t==1) 0 else {
      idx <- seq_len(min(p, t-1))
      sum(phi[idx] * y[t-idx])
    }
    ma <- if (q==0 || t==1) 0 else {
      idx <- seq_len(min(q, t-1))
      sum(theta[idx] * e[t-idx])
    }
    yhat[t] <- cst + ar + ma
    e[t] <- y[t] - yhat[t]
  }
  list(e=e, yhat=yhat, s2=s2)
}

#ARMA(1,0) was the best fit:
fit <- fit_arma_cml(y,1,0)
est <- fit$est
rec <- arma_recursion(y, p=1, q=0, pars=est)
m <- 2L
dates_eff <- ipca_df$date[(m+1):Tn]
fitted_eff <- rec$yhat[(m+1):Tn]
resid_eff  <- rec$e[(m+1):Tn]
sigma_hat  <- sqrt(exp(est["log_sigma2"]))

#Make plots
p_fit <- ggplot() +
  geom_line(aes(ipca_df$date, y), color="red") +
  geom_line(aes(dates_eff, fitted_eff),color="blue") +
  labs(title="IPCA: actual (red) vs fitted (ARMA(1,0))", x=NULL, y="IPCA") +
  theme_minimal(base_size=12)
print(p_fit)

p_res <- ggplot(data.frame(date = dates_eff, resid = resid_eff),
                aes(date, resid)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line() +
  labs(title="Residuals from ARMA(1,0)", x=NULL, y="Residual") +
  theme_minimal(base_size=12)
print(p_res)

# --- (f) Residual ACF plot with CI bands --------------------------------------
acf_obj <- acf(resid_eff, plot=FALSE, lag.max = 36)
acf_df <- data.frame(lag = acf_obj$lag[-1]*length(resid_eff),
                     acf = acf_obj$acf[-1])
ci <- 1.96 / sqrt(length(resid_eff))
p_acf <- ggplot(acf_df, aes(x=lag, y=acf)) +
  geom_hline(yintercept = c(-ci, ci), linetype=2) +
  geom_segment(aes(x=lag, xend=lag, y=0, yend=acf)) +
  labs(title="ACF of residuals (ARMA(0,1))", x="Lag", y="ACF") +
  theme_minimal(base_size=12)
print(p_acf)

# --- (g) One simulated 24-step forecast path ----------------------------------
arma_forecast_path <- function(y_last, e_last, pars, H) {
  # y_last: c(y_T, y_{T-1}); e_last: c(e_T, e_{T-1})
  cst <- pars["c"]; phi1 <- pars["phi1"]; phi2 <- pars["phi2"]
  th1 <- pars["theta1"]; th2 <- pars["theta2"]
  s2 <- exp(pars["log_sigma2"])
  y_hist <- y_last         # length 2: (y_T, y_{T-1})
  e_hist <- e_last         # length 2: (e_T, e_{T-1})
  eps <- rnorm(H, 0, sqrt(s2))
  y_future <- numeric(H)
  for (h in 1:H) {
    ma_part <- eps[h] + ifelse(is.na(th1),0,th1)*e_hist[1] + 
      ifelse(is.na(th2),0,th2)*e_hist[2]
    ar_part <- ifelse(is.na(phi1),0,phi1)*y_hist[1] +
      ifelse(is.na(phi2),0,phi2)*y_hist[2]
    y_next <- cst + ar_part + ma_part
    y_future[h] <- y_next
    # update histories (most recent first)
    y_hist <- c(y_next, y_hist[1])
    e_hist <- c(eps[h], e_hist[1])
  }
  list(y_future = y_future)
}

H <- 24
y_last <- c(y[Tn], y[Tn-1])
e_full <- rec$e
e_last <- c(e_full[Tn], e_full[Tn-1])

set.seed(123)
one_path <- arma_forecast_path(y_last, e_last, est, H)$y_future
future_dates <- seq(from = ipca_df$date[Tn] %m+% months(1), by = "1 month", length.out = H)

p_one <- ggplot() +
  geom_line(aes(ipca_df$date, y), color="grey40") +
  geom_line(aes(future_dates, one_path),color="blue") +
  labs(title="Observed series and one simulated 24-month path (ARMA(1,0))",
       x=NULL, y="IPCA") +
  theme_minimal(base_size=12)
print(p_one)

# --- (h) 1000 simulated paths, 5%/95% bands -----------------------------------

H <- 24
y_last <- c(y[Tn], y[Tn-1])
e_full <- rec$e
e_last <- c(e_full[Tn], e_full[Tn-1])

B <- 1000L
set.seed(2025)
paths <- matrix(NA_real_, nrow = H, ncol = B)
for (b in 1:B) {
  paths[, b] <- arma_forecast_path(y_last, e_last, est, H)$y_future
}

# Quantiles with na.rm=TRUE (belt-and-suspenders)
q_lo <- apply(paths, 1, quantile, probs = 0.05, na.rm = TRUE)
q_md <- apply(paths, 1, quantile, probs = 0.50, na.rm = TRUE)
q_hi <- apply(paths, 1, quantile, probs = 0.95, na.rm = TRUE)

df_fore <- tibble(date = future_dates, path_g = one_path, q05 = q_lo, q50 = q_md, q95 = q_hi)
p_band <- ggplot() +
  geom_line(aes(ipca_df$date, y), color="grey40") +
  geom_ribbon(data=df_fore, aes(x=date, ymin=q05, ymax=q95), alpha=0.25) +
  geom_line(data=df_fore, aes(date, path_g)) +
  labs(title="ARMA(1,0) forecast, median and 90% band (parametric bootstrap)",
       x=NULL, y="IPCA") +
  theme_minimal(base_size=12)
print(p_band)



# --- (Latex Help) -------------------------------------------------------------

## ==== Assumes you already have `summ_tbl` with columns:
## model, param, estimate, ci95, logLik, AIC, BIC

## 1) Remove variance rows
tb <- subset(summ_tbl, !(param %in% c("sigma","log_sigma2")))

## 2) Format numbers as strings
tb$Estimate_txt <- ifelse(is.na(tb$estimate), "", sprintf("%.3f", tb$estimate))
tb$logLik_txt   <- ifelse(is.na(tb$logLik),   "", sprintf("%.2f", tb$logLik))
tb$AIC_txt      <- ifelse(is.na(tb$AIC),      "", sprintf("%.2f", tb$AIC))
tb$BIC_txt      <- ifelse(is.na(tb$BIC),      "", sprintf("%.2f", tb$BIC))

## 3) Build LaTeX body rows, grouped by model; model name only on first line
mods <- unique(tb$model)
rows <- character(0)

for (mod in mods) {
  df <- tb[tb$model == mod, , drop = FALSE]
  
  ## Order columns as: c, phi1.., theta1.., then anything else
  ix_c     <- which(df$param == "c")
  ix_phi   <- grep("^phi[0-9]+$", df$param)
  ix_theta <- grep("^theta[0-9]+$", df$param)
  if (length(ix_phi))   ix_phi   <- ix_phi[order(as.integer(sub("^phi",   "", df$param[ix_phi])))]
  if (length(ix_theta)) ix_theta <- ix_theta[order(as.integer(sub("^theta","", df$param[ix_theta])))]
  
  ord <- c(ix_c, ix_phi, ix_theta, setdiff(seq_len(nrow(df)), c(ix_c, ix_phi, ix_theta)))
  df  <- df[ord, , drop = FALSE]
  
  model_col <- c(mod, rep("", nrow(df) - 1L))
  
  for (i in seq_len(nrow(df))) {
    line <- sprintf("%s & %s & %s & %s & %s & %s & %s \\\\",
                    model_col[i],
                    df$param[i],
                    df$Estimate_txt[i],
                    df$ci95[i],
                    df$logLik_txt[i],
                    df$AIC_txt[i],
                    df$BIC_txt[i])
    rows <- c(rows, line)
  }
  rows <- c(rows, "\\addlinespace")
}

body <- paste(rows, collapse = "\n")

## 4) Full LaTeX table code (booktabs + threeparttable)
caption <- "ARMA CML estimates (95\\% Wald intervals)"
label   <- "tab:summ_tbl"

table_tex <- sprintf(
  "%% Requires in preamble: \\usepackage{booktabs,threeparttable}
\\begin{table}[!ht]
\\centering
\\begin{threeparttable}
\\caption{%s}
\\label{%s}
\\setlength{\\tabcolsep}{6pt}
\\renewcommand{\\arraystretch}{1.1}
\\begin{tabular}{l l r l r r r}
\\toprule
Model & Parameter & Estimate & 95\\%% CI & logLik & AIC & BIC \\\\
\\midrule
%s
\\bottomrule
\\end{tabular}
\\begin{tablenotes}[flushleft]
\\footnotesize Notes: Gaussian conditional maximum likelihood (CML). Wald 95\\%% CIs.
Model-level statistics repeat across parameter rows. Innovation variance rows omitted.
\\end{tablenotes}
\\end{threeparttable}
\\end{table}
", caption, label, body)

## 5) Output:
cat(table_tex)

#