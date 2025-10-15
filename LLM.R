library(lme4)
library(lmerTest) 
library(emmeans)
library(broom.mixed)
library(ggplot2)
library(dplyr)

dat <- fits_all %>%
  mutate(
    # Ensure rep is factor
    rep = as.factor(rep),
    Test = as.factor(Test),
     Treatment = case_when(
      grepl("_24hr", Test, ignore.case = TRUE) ~ "After24",
      grepl("^T4($|_)", Test, ignore.case = TRUE) & !grepl("_24hr", Test, ignore.case = TRUE) ~ "Initial_Drought",
      TRUE ~ "Initial"
    )
  ) %>%
  mutate(Treatment = factor(Treatment, levels = c("Initial", "After24", "Initial_Drought")))


# Remove rows with NA thresholds for each model separately
dat_tcrit <- dat %>% filter(!is.na(Tcrit))
dat_t50   <- dat %>% filter(!is.na(T50))

# ---- Fit random-intercept LMMs
# Use REML for parameter estimates; use ML for model comparisons
m_tcrit_re  <- lmer(Tcrit ~ Treatment + (1 | rep), data = dat_tcrit, REML = TRUE)
m_t50_re    <- lmer(T50   ~ Treatment + (1 | rep), data = dat_t50,   REML = TRUE)

# Summaries and ANOVA (Type III)
summary(m_tcrit_re)
anova(m_tcrit_re, type = "III")  
summary(m_t50_re)
anova(m_t50_re, type = "III")

# ---- Estimated marginal means / pairwise comparisons ----
emm_tcrit <- emmeans(m_tcrit_re, ~ Treatment)
pairs(emm_tcrit, adjust = "tukey")
summary(emm_tcrit)

emm_t50 <- emmeans(m_t50_re, ~ Treatment)
pairs(emm_t50, adjust = "tukey")
summary(emm_t50)

# emmeans tables
emm_tcrit_df <- as.data.frame(emm_tcrit)
emm_t50_df   <- as.data.frame(emm_t50)


# Repeat for T50
m_t50_rs <- try(lmer(T50 ~ Treatment + (1 + Treatment | rep), data = dat_t50, REML = FALSE), silent = TRUE)
m_t50_re_ML <- update(m_t50_re, REML = FALSE)
if(!inherits(m_t50_rs, "try-error")) {
  anova(m_t50_re_ML, m_t50_rs)
} else {
  cat("Random-slope model for T50 failed to converge or errored; keep random-intercept model.\n")
}

# ---- Post-hoc effect sizes: difference in emmeans ----
emm_diff_tcrit <- contrast(emm_tcrit, method = "pairwise", adjust = "tukey")
summary(emm_diff_tcrit)  # shows estimates, SE, t.ratio, p.value

