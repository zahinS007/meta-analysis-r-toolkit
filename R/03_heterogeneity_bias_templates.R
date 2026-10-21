# =============================================================================
# 03_heterogeneity_bias_templates.R — explore heterogeneity + publication bias
# Works on any meta object `m` (metabin / metacont / metagen / metaprop ...)
# =============================================================================
source("R/ma_toolkit.R")
library(metafor)
# m <- m.bin   # your main model from 01/02

# I² guide (Higgins 2002): <=25% low · 26–50% moderate · >50% substantial
summary(m)

# --- One call does LOO + Baujat + funnel (+Egger if k>=10) + influence -------
sens <- ma_sensitivity(m, prefix = "OUTCOME")

# --- Leave-one-out ----------------------------------------------------------
loo <- metainf(m, pooled = "random")
meta::forest(loo, rightcols = c("effect", "ci", "I2"),
       col.diamond = "gray",
       xlab = "Favors INTERVENTION    Favors CONTROL", ff.xlab = "bold")

# --- Baujat (x = contribution to Q; y = influence on pooled effect) ----------
par(cex.lab = 1.2, font.lab = 2); baujat(m, bg = "red", cex.studlab = 1.1)

# --- Subgroups (pre-specify in protocol!) ------------------------------------
m.sub <- update(m, subgroup = SUBGROUP_COLUMN)
summary(m.sub); ma_forest(m.sub)

# --- Meta-regression (study-level covariate: year, mean age, % DM, dose ...) -
#     Rule of thumb: >= 10 studies per covariate.
mreg <- metareg(m, ~ COVARIATE)
mreg
bubble(mreg, studlab = TRUE, bg = "red", cex.studlab = 0.6, pos.studlab = 3,
       xlab = "COVARIATE LABEL",
       ylab = "Treatment effect (log scale)")       # add xlim/ylim to fit bubbles

# --- Influence diagnostics (metafor; replaces dmetar::InfluenceAnalysis) -----
rm <- rma(yi = m$TE, sei = m$seTE, method = "REML", slab = m$studlab)
inf <- influence(rm); inf; plot(inf)                 # flagged studies marked *

# --- GOSH (all 2^k - 1 subsets; slow for k > 15) -----------------------------
g <- gosh(rm, subsets = if (m$k > 15) 1e5 else NA)
plot(g, alpha = 0.1, col = "darkcyan")
# cluster diagnostics need dmetar:
# remotes::install_github("MathiasHarrer/dmetar"); dmetar::gosh.diagnostics(g)

# --- L'Abbé plot (binary outcomes) -------------------------------------------
labbe(m, studlab = TRUE, bg = "red")

# =============================================================================
# PUBLICATION BIAS
# =============================================================================
funnel(m, studlab = TRUE, bg = "red", contour = c(0.9, 0.95, 0.99),
       col.contour = c("gray75", "gray85", "gray95"))
legend("topright", c("p < 0.10", "p < 0.05", "p < 0.01"),
       fill = c("gray75", "gray85", "gray95"), bty = "n")

if (m$k >= 10) {
  metabias(m, method.bias = "linreg")        # Egger
  metabias(m, method.bias = "Harbord")       # binary OR alternative
  tf <- trimfill(m); summary(tf); funnel(tf) # trim-and-fill sensitivity
} else message("k < 10: report funnel visually only; tests are underpowered.")

# =============================================================================
# SAVE everything
# =============================================================================
ma_save(funnel(m, studlab = TRUE, bg = "red"), "OUTCOME_funnel", 8, 7)
ma_save(bubble(mreg, studlab = TRUE, bg = "red"), "OUTCOME_metareg", 8, 7)
