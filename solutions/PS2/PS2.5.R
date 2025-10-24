# ================================
# IPCA vs IGP-M: Descriptive Analysis
# ================================

# ---- 0) Packages ----
library(here)
library(readr)
library(tidyverse)

# ---- 1) Load & clean data ----
csv_path <-  here("solutions", "PS2", "ipca_igpm.csv")
df <- readr::read_csv(csv_path, show_col_types = FALSE)
names(df) <- tolower(names(df))
df$date <- as.Date(paste0(df$date, "-01"))
df$date <- as.Date(df$date)
df <- df[!is.na(df$igpm) & !is.na(df$ipca), ]

# ---- 2) Time-series plot (MoM) ----
df_long <- df |>
  tidyr::pivot_longer(ipca:igpm, names_to = "index", values_to = "mom") |>
  mutate(index = toupper(index))


p <- ggplot(df_long, aes(date, mom, color = index)) +
  geom_line() +
  labs(title = "Monthly MoM Inflation: IPCA vs IGP-M",
       x = NULL, y = "MoM (%)", color = NULL) +
  theme_minimal(base_size = 12)
print(p)

# ---- 3) Summary stats & autocovariances ----
summarise_series <- function(x) {
  data.frame(
    mean      = mean(x, na.rm = TRUE),
    variance  = var(x, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

as_cov <- function(x, L = 24) {
  ac <- acf(x, type = "covariance", plot = FALSE, na.action = na.omit, lag.max = L)
  data.frame(lag = as.integer(ac$lag[, , 1]),
             autocov = as.numeric(ac$acf[, , 1]))
}

ip_stats <- summarise_series(df$ipca); rownames(ip_stats) <- "IPCA"
ig_stats <- summarise_series(df$igpm); rownames(ig_stats) <- "IGP-M"
tab_stats <- rbind(ip_stats, ig_stats)

cat("\n=== Sample Mean & Variance (MoM, in % units) ===\n")
print(round(tab_stats, 4))

ip_cov <- as_cov(df$ipca, 24); ip_cov$series <- "IPCA"
ig_cov <- as_cov(df$igpm, 24); ig_cov$series <- "IGP-M"

cat("\n=== IPCA Autocovariances (lags 0..24) ===\n")
print(round(ip_cov[, c("lag","autocov")], 5))
cat("\n=== IGP-M Autocovariances (lags 0..24) ===\n")
print(round(ig_cov[, c("lag","autocov")], 5))


# ---- 4) ACF plots (persistence) ----
par(mfrow = c(1, 2))
acf(na.omit(df$ipca), lag.max = 24, main = "IPCA: ACF (0–24)")
acf(na.omit(df$igpm), lag.max = 24, main = "IGP-M: ACF (0–24)")
par(mfrow = c(1, 1))

