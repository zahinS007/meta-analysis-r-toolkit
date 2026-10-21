# =============================================================================
# 02_single_arm_templates.R — proportions (prevalence/incidence) & single means
# =============================================================================
source("R/ma_toolkit.R")
ma <- ma_import("data.xlsx")

# -----------------------------------------------------------------------------
# A. PROPORTIONS (single-arm: events / n)   -> metaprop
#    Sheet columns: Author, Year, event, n
#    RCTs can be included — use the intervention arm only.
#    Transformation choice:
#      GLMM + PLOGIT  : default & recommended (no study weights shown)
#      Inverse + PLOGIT: proportions < 0.2 or > 0.8, need weights
#      Inverse + PFT  : many 0% or 100% studies (sensitivity-check with PLOGIT)
#      Inverse + PRAW : proportions between 0.2 and 0.8
# -----------------------------------------------------------------------------
m.prop <- metaprop(event, n,
                   data = ma$OUTCOME,
                   studlab = Author,
                   method = "GLMM",          # or "Inverse" if you need weights
                   sm = "PLOGIT",            # "PFT", "PRAW", "PAS"
                   method.ci = "CP",         # Clopper-Pearson study CIs
                   common = FALSE, random = TRUE)
summary(m.prop)
ma_report(m.prop)
ma_forest(m.prop, xlab = "Prevalence of OUTCOME (%)", pscale = 100,
          year_col = "Year")
# NB: with method = "GLMM" the Weight column is empty; drop it with
# leftcols = c("studlab","event","n","effect","ci")

# -----------------------------------------------------------------------------
# B. SINGLE MEANS   -> metamean
#    Sheet columns: Author, n, mean, sd   (sm = "MRAW"; "MLN" for skewed data)
# -----------------------------------------------------------------------------
# B1. Some studies give median + IQR/range? Estimate mean (SD) first.
#     Columns: n, median, q1, q3, min, max (leave blank cells as NA)
#     qe = quantile estimation (recommended default); bc = Box-Cox; mln
med <- ma$OUTCOME_MEDIANS
med <- ma_median_to_mean_df(med, method = "qe")
med[, c("Author", "est.mean", "est.sd")]
# one-off:
ma_median_to_mean(n = 1390, median = 38.5, q1 = 35.7, q3 = 40.5)

# B2. Pool
m.mean <- metamean(n, mean, sd,
                   data = ma$OUTCOME,
                   studlab = Author,
                   sm = "MRAW",
                   method.tau = "REML",
                   common = FALSE, random = TRUE)
summary(m.mean)
ma_forest(m.mean, xlab = "OUTCOME (units)")
