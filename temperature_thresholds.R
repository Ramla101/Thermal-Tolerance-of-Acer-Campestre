# ---- Load libraries ----
library(tidyverse)
library(minpack.lm)
library(broom)

data_folder <- "."  # folder containing CSV files
temp_col <- "Temperature"
fv_col   <- "Fv/m"
id_col   <- "Sample #"
file_names <- c("T1", "T1_24hr", "T2", "T2_24hr", "T3", "T3_24hr", "T4", "T4_24hr")

# ---- combine all CSV files ----
df_list <- list()

for (fname in file_names) {
  file_path <- file.path(data_folder, paste0(fname, ".csv"))
  
  if (file.exists(file_path)) {
    cat("Reading", fname, "...\n")
    
    # Read the CSV
    temp_df <- read_csv(file_path, guess_max = 10000, show_col_types = FALSE)
    
    # Add Test column
    temp_df$Test <- fname
    
    df_list[[fname]] <- temp_df
  } else {
    warning(paste("File not found:", file_path))
  }
}

# Combine all dataframes
df_raw <- bind_rows(df_list)

cat("\nTotal rows read:", nrow(df_raw), "\n")
cat("Tests found:", paste(unique(df_raw$Test), collapse = ", "), "\n")

df <- df_raw %>%
  rename(
    rep = !!rlang::sym(id_col),
    Temp = !!rlang::sym(temp_col),
    Fv.m = !!rlang::sym(fv_col)
  )

# Convert to numeric and filter
df <- df %>%
  mutate(
    Temp = as.numeric(Temp),
    Fv.m = as.numeric(Fv.m),
    rep = as.factor(rep),
    Test = as.factor(Test)
  ) %>%
  filter(!is.na(Temp), !is.na(Fv.m)) %>%
  filter(Temp >= 20 & Temp <= 60)

# Check data per Test
cat("\nData summary by Test:\n")
df %>%
  group_by(Test) %>%
  summarise(
    n_obs = n(),
    n_reps = n_distinct(rep),
    temp_range = paste(min(Temp), "-", max(Temp)),
    .groups = "drop"
  ) %>%
  print()

# ---- Summarise mean & SE per Temp per Test ----
df_avg <- df %>%
  group_by(Test, Temp) %>%
  summarise(
    Avg_Value = mean(Fv.m, na.rm = TRUE),
    SE = sd(Fv.m, na.rm = TRUE) / sqrt(sum(!is.na(Fv.m))),
    .groups = "drop"
  )

# ---- Temperature grid and logistic function ----
x_grid <- seq(20, 60, by = 0.1)

logistic_fun <- function(x, d, b, c) {
  d / (1 + exp(b * (x - c)))
}

# ---- Fit per replicate per Test ----
fit_rep <- function(df_rep) {
  x <- df_rep$Temp
  y <- df_rep$Fv.m
  start_d <- max(y, na.rm = TRUE)
  start_b <- 0.3
  start_c <- mean(range(x, na.rm = TRUE))
  fit <- tryCatch({
    nlsLM(y ~ logistic_fun(x, d, b, c),
          start = list(d = start_d, b = start_b, c = start_c),
          control = nls.lm.control(maxiter = 1000))
  }, error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  coefs <- coef(fit)
  d <- coefs["d"]; b <- coefs["b"]; c <- coefs["c"]
  pred <- tibble(x = x_grid, yhat = logistic_fun(x_grid, d, b, c))
  deriv <- - d * b * exp(b * (x_grid - c)) / (1 + exp(b * (x_grid - c)))^2
  steepest <- - d * b / 4
  threshold_deriv <- 0.15 * steepest
  idx_crit <- which(deriv <= threshold_deriv)
  Tcrit <- if (length(idx_crit)>0) x_grid[min(idx_crit)] else NA_real_
  half_val <- 0.5 * d
  idx_t50 <- which(pred$yhat <= half_val)
  T50 <- if (length(idx_t50)>0) pred$x[min(idx_t50)] else NA_real_
  tibble(
    rep = unique(df_rep$rep),
    Test = unique(df_rep$Test),
    d = d, b = b, c = c,
    Tcrit = Tcrit, T50 = T50,
    pred = list(pred)
  )
}

# Apply fitting grouped by Test and rep
fits_all <- df %>%
  group_by(Test, rep) %>%
  group_split() %>%
  map(fit_rep) %>%
  compact() %>%
  bind_rows()

# Expand predicted curves for plotting
pred_indiv <- fits_all %>% select(Test, rep, pred) %>% unnest(pred)

# Mean predicted curve per Test
mean_pred <- pred_indiv %>%
  group_by(Test, x) %>%
  summarise(mean_yhat = mean(yhat, na.rm=TRUE), .groups = "drop")

# Fit logistic to mean_pred per Test
mean_thresholds <- mean_pred %>%
  group_by(Test) %>%
  nest() %>%
  mutate(mean_fit = map(data, ~tryCatch({
    nlsLM(mean_yhat ~ logistic_fun(x, d, b, c),
          data = .x,
          start = list(d = max(.x$mean_yhat, na.rm=TRUE), b = 0.3, c = mean(range(.x$x))),
          control = nls.lm.control(maxiter = 1000))
  }, error = function(e) NULL))) %>%
  mutate(
    coefs = map(mean_fit, ~ if (!is.null(.x)) coef(.x) else NULL),
    mean_Tcrit = map_dbl(coefs, function(co) {
      if (is.null(co)) return(NA_real_)
      d_m <- co["d"]; b_m <- co["b"]; c_m <- co["c"]
      deriv_m <- - d_m * b_m * exp(b_m * (x_grid - c_m)) / (1 + exp(b_m * (x_grid - c_m)))^2
      threshold_deriv_m <- 0.15 * (-d_m * b_m / 4)
      idx <- which(deriv_m <= threshold_deriv_m)
      if (length(idx)>0) x_grid[min(idx)] else NA_real_
    }),
    mean_T50 = map_dbl(coefs, function(co) {
      if (is.null(co)) return(NA_real_)
      d_m <- co["d"]; b_m <- co["b"]; c_m <- co["c"]
      pred_m <- logistic_fun(x_grid, d_m, b_m, c_m)
      idx <- which(pred_m <= 0.5 * d_m)
      if (length(idx)>0) x_grid[min(idx)] else NA_real_
    })
  ) %>%
  select(Test, mean_Tcrit, mean_T50)

# Merge with mean_pred
mean_pred <- mean_pred %>% left_join(mean_thresholds, by = "Test")

# Summary of replicate thresholds per Test
summary_rep <- fits_all %>%
  group_by(Test) %>%
  summarise(
    mean_Tcrit_rep = mean(Tcrit, na.rm = TRUE),
    sd_Tcrit_rep = sd(Tcrit, na.rm = TRUE),
    mean_T50_rep = mean(T50, na.rm = TRUE),
    sd_T50_rep = sd(T50, na.rm = TRUE),
    n_rep = n(),
    .groups = "drop"
  )

# labels for plotting ----
labels_df <- summary_rep %>%
  rename(Tcrit = mean_Tcrit_rep, T50 = mean_T50_rep) %>%
  mutate(
    y_pos_tcrit = 0.08,
    y_pos_t50 = 0.65,
    Tcrit_label = if_else(!is.na(Tcrit), 
                          paste0("bold(T[crit] == ", round(Tcrit, 1), "*degree*C)"), 
                          NA_character_),
    T50_label = if_else(!is.na(T50), 
                        paste0("bold(T[50] == ", round(T50, 1), "*degree*C)"), 
                        NA_character_)
  )


test_dates <- tibble(
  Test = c("T1", "T1_24hr", "T2", "T2_24hr", "T3", "T3_24hr", "T4", "T4_24hr"),
  Date = c("03 Jul 2023", "04 Jul 2023", "13 Jul 2023", "14 Jul 2023",
           "23 Jul 2023", "24 Jul 2023", "28 Aug 2023", "29 Aug 2023")
)

# Merge to all relevant data frames used in plotting
df_avg       <- df_avg       %>% left_join(test_dates, by = "Test")
pred_indiv   <- pred_indiv   %>% left_join(test_dates, by = "Test")
mean_pred    <- mean_pred    %>% left_join(test_dates, by = "Test")
labels_df    <- labels_df    %>% left_join(test_dates, by = "Test")

df_avg$facet_label     <- paste0(df_avg$Test, ": ", df_avg$Date)
pred_indiv$facet_label <- paste0(pred_indiv$Test, ": ", pred_indiv$Date)
mean_pred$facet_label  <- paste0(mean_pred$Test, ": ", mean_pred$Date)
labels_df$facet_label  <- paste0(labels_df$Test, ": ", labels_df$Date)

p <- ggplot() +
  geom_line(data = pred_indiv, aes(x = x, y = yhat, group = rep),
            color = "darkseagreen4", linetype = "dashed", alpha = 0.6) +
  geom_line(data = mean_pred, aes(x = x, y = mean_yhat),
            color = "purple4", size = 1.1) +
  geom_errorbar(data = df_avg,
                aes(x = Temp, ymin = Avg_Value - SE, ymax = Avg_Value + SE),
                width = 0.6, color = "steelblue4") +
  geom_point(data = df_avg, aes(x = Temp, y = Avg_Value),
             color = "steelblue4") +
  geom_vline(data = labels_df, aes(xintercept = Tcrit),
             color = "red3", linetype = "dashed", size = 0.8, na.rm = TRUE) +
  geom_vline(data = labels_df, aes(xintercept = T50),
             color = "red3", linetype = "dashed", size = 0.8, na.rm = TRUE) +
  geom_text(data = labels_df %>% filter(!is.na(Tcrit)),
            aes(x = Tcrit, y = y_pos_tcrit, label = Tcrit_label),
            color = "red3", hjust = -0.1, size = 3.5, parse = TRUE) +
  geom_text(data = labels_df %>% filter(!is.na(T50)),
            aes(x = T50, y = y_pos_t50, label = T50_label),
            color = "red3", hjust = -0.1, size = 3.5, parse = TRUE) +
  labs(x = "Temperature (°C)", y = "Fv/Fm") +
  scale_x_continuous(limits = c(20, 60), breaks = seq(20, 60, 2)) +
  scale_y_continuous(limits = c(0, 0.85), breaks = seq(0, 0.8, 0.1)) +
  theme_classic(base_size = 14) +
  theme(panel.border = element_rect(color = "black", fill = NA)) +
  facet_wrap(~facet_label, ncol = 2)  # <<— use combined label

p


# Save plot
ggsave("TRC_panel.tif", p, width = 12, height = 15, dpi = 300)



#Confidence Interval
library(dplyr)
library(readr)

# ---- compute replicate-based summaries & empirical 95% percentile CIs ----
replicate_summary <- fits_all %>%
  filter(!is.na(Tcrit) | !is.na(T50)) %>%   # keep rows where at least one threshold exists
  group_by(Test) %>%
  summarise(
    n_reps_total = n(),                                 # total replicate rows for that test (including NA fits)
    n_successful_Tcrit = sum(!is.na(Tcrit)),
    n_successful_T50  = sum(!is.na(T50)),
    # Tcrit summaries (use only non-NA)
    Tcrit_mean = ifelse(n_successful_Tcrit>0, mean(Tcrit, na.rm = TRUE), NA_real_),
    Tcrit_sd   = ifelse(n_successful_Tcrit>1, sd(Tcrit, na.rm = TRUE), NA_real_),
    Tcrit_min  = ifelse(n_successful_Tcrit>0, min(Tcrit, na.rm = TRUE), NA_real_),
    Tcrit_p025 = ifelse(n_successful_Tcrit>0, quantile(Tcrit, probs = 0.025, na.rm = TRUE, type = 7), NA_real_),
    Tcrit_p975 = ifelse(n_successful_Tcrit>0, quantile(Tcrit, probs = 0.975, na.rm = TRUE, type = 7), NA_real_),
    Tcrit_max  = ifelse(n_successful_Tcrit>0, max(Tcrit, na.rm = TRUE), NA_real_),
    # T50 summaries (use only non-NA)
    T50_mean = ifelse(n_successful_T50>0, mean(T50, na.rm = TRUE), NA_real_),
    T50_sd   = ifelse(n_successful_T50>1, sd(T50, na.rm = TRUE), NA_real_),
    T50_min  = ifelse(n_successful_T50>0, min(T50, na.rm = TRUE), NA_real_),
    T50_p025 = ifelse(n_successful_T50>0, quantile(T50, probs = 0.025, na.rm = TRUE, type = 7), NA_real_),
    T50_p975 = ifelse(n_successful_T50>0, quantile(T50, probs = 0.975, na.rm = TRUE, type = 7), NA_real_),
    T50_max  = ifelse(n_successful_T50>0, max(T50, na.rm = TRUE), NA_real_)
  ) %>%
  ungroup()

replicate_summary