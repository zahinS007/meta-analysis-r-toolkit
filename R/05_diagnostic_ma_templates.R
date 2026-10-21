# =============================================================================
# 05_diagnostic_ma_templates.R — diagnostic test accuracy (DTA) meta-analysis
# Bivariate Reitsma model (mada) + univariate sens/spec forest plots (meta)
# Sheet columns (exact names): study, TP, FP, FN, TN   — one sheet per test
# No 2x2 table? Back-calculate from prevalence, n, sens, spec:
#   diseased = round(n*prev); TP = round(sens*diseased); FN = diseased - TP
#   healthy = n - diseased;   TN = round(spec*healthy);  FP = healthy - TN
# =============================================================================
library(mada); library(meta); library(readxl); library(dplyr)
source("R/ma_toolkit.R")

dta <- ma_import("dta.xlsx")
d1  <- dta$TEST_A                    # e.g. MRI
# d2 <- dta$TEST_B                   # e.g. CECT (for comparisons)

# helper: back-calculate 2x2 from summary values
dta_from_summary <- function(n, prev, sens, spec) {
  dis <- round(n * prev); hea <- n - dis
  TP <- round(sens * dis); TN <- round(spec * hea)
  data.frame(TP = TP, FN = dis - TP, TN = TN, FP = hea - TN)
}

# ---- 1. Univariate pooled sensitivity & specificity (forest plots) ----------
sens1 <- metaprop(event = TP, n = TP + FN, data = d1, studlab = study,
                  method = "GLMM", sm = "PLOGIT", method.ci = "CP",
                  common = FALSE, random = TRUE)
spec1 <- metaprop(event = TN, n = TN + FP, data = d1, studlab = study,
                  method = "GLMM", sm = "PLOGIT", method.ci = "CP",
                  common = FALSE, random = TRUE)
dta_forest <- function(m, lab) meta::forest(m, layout = "RevMan5", xlim = c(0, 1),
  leftcols = c("studlab", "event", "n", "effect", "ci"),
  leftlabs = c("Study", if (lab == "Sensitivity") "TP" else "TN",
               if (lab == "Sensitivity") "TP+FN" else "TN+FP", lab, "95% CI"),
  rightcols = FALSE, just = "center", random = TRUE, common = FALSE,
  col.square = "blue", col.diamond = "black", print.tau.ci = TRUE,
  colgap.forest = "15mm")
ma_save(dta_forest(sens1, "Sensitivity"), "TEST_A_sensitivity", 10, 2.5 + 0.3 * nrow(d1))
ma_save(dta_forest(spec1, "Specificity"), "TEST_A_specificity", 10, 2.5 + 0.3 * nrow(d1))

# ---- 2. Bivariate model (Reitsma) — the primary DTA analysis ----------------
fit1 <- reitsma(d1)                  # needs columns TP FN FP TN
summary(fit1)                        # pooled sens, FPR (1-spec), AUC, I² (Zhou/Dendukuri)
AUC(fit1)

# ---- 3. SROC curve ----------------------------------------------------------
sroc_plot <- function() {
  plot(fit1, xlim = c(0, 1), ylim = c(0, 1),
       main = "SROC curve of TEST_A for TARGET_CONDITION")
  points(fpr(d1), sens(d1), pch = 2, cex = 0.7)
  legend("bottomright", c("Studies", "Summary estimate"), pch = c(2, 1), bty = "n")
}
ma_save(sroc_plot(), "TEST_A_sroc", 7, 7)

# ---- 4. Other summary measures ----------------------------------------------
mr <- madauni(d1, type = "DOR", method = "DSL"); summary(mr)   # pooled DOR
madauni(d1, type = "posLR"); madauni(d1, type = "negLR")       # LR+ / LR-

# ---- 5. Publication bias — Deeks' funnel asymmetry test (k >= 10) -----------
# Correct DOR set-up: metabin(TP, TP+FN, FP, FP+TN, sm = "DOR")
dor <- metabin(TP, TP + FN, FP, FP + TN, data = d1, studlab = study,
               sm = "DOR", method = "Inverse", common = FALSE, random = TRUE)
ma_save(funnel(dor, studlab = TRUE), "TEST_A_funnel_DOR", 8, 7)
if (nrow(d1) >= 10) metabias(dor, method.bias = "Deeks")

# ---- 6. Comparing two tests -------------------------------------------------
# a) Meta-regression on the bivariate model (preferred)
# both <- bind_rows(mutate(d1, test = "TEST_A"), mutate(d2, test = "TEST_B"))
# fitc <- reitsma(both, formula = cbind(tsens, tfpr) ~ test); summary(fitc)
# b) Overlaid SROC
# fit2 <- reitsma(d2)
# plot(fit1, xlim = c(0,1), ylim = c(0,1)); lines(sroc(fit2), lty = 2)
# ROCellipse(fit2, lty = 2, pch = 2, add = TRUE)
# legend("bottomright", c("TEST_A", "TEST_B"), lty = 1:2, pch = 1:2)
# c) Difference in AUC with bootstrap CI and p (dmetatools, GitHub):
# devtools::install_github("nomahi/dmetatools")
# dmetatools::AUC_comparison(d1$TP, d1$FP, d1$FN, d1$TN, d2$TP, d2$FP, d2$FN, d2$TN)
# d) Sensitivity by test as subgroups (p for difference):
# sc <- metaprop(TP, TP + FN, data = both, studlab = study, subgroup = test,
#                method = "GLMM", sm = "PLOGIT"); summary(sc)
