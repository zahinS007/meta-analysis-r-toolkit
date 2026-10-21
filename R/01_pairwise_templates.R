# =============================================================================
# 01_pairwise_templates.R — two-arm meta-analyses (binary, continuous, generic)
# Copy the block you need, replace the CAPS placeholders, run.
# =============================================================================
source("R/ma_toolkit.R")
setwd("PATH/TO/YOUR/WORKING/DIRECTORY")      # folder holding data.xlsx
ma <- ma_import("data.xlsx")                 # one sheet per outcome

# -----------------------------------------------------------------------------
# A. BINARY OUTCOME (events / totals per arm)  -> OR or RR
#    Sheet columns: Author, Year, event.e, n.e, event.c, n.c
# -----------------------------------------------------------------------------
m.bin <- metabin(event.e, n.e, event.c, n.c,
                 data = ma$OUTCOME,            # e.g. ma$acm
                 studlab = Author,
                 sm = "RR",                    # "RR" or "OR" (or "RD")
                 method = "MH",                # MH default; "Peto" for rare events + OR
                 method.tau = "PM",            # PM/REML preferred over DL (see guide)
                 method.random.ci = "HK",      # Hartung-Knapp; use "classic" to switch off
                 common = FALSE, random = TRUE)
summary(m.bin)
ma_report(m.bin)

ma_forest(m.bin, lab.e = "INTERVENTION", lab.c = "CONTROL", year_col = "Year")
sz <- ma_forest_size(m.bin)
ma_save(ma_forest(m.bin, lab.e = "INTERVENTION", lab.c = "CONTROL", year_col = "Year"),
        "OUTCOME_forest", width = sz["width"], height = sz["height"])

# -----------------------------------------------------------------------------
# B. CONTINUOUS OUTCOME (mean, SD, n per arm)  -> MD or SMD
#    Sheet columns: Author, n.e, mean.e, sd.e, n.c, mean.c, sd.c
#    Medians/IQR? convert first with ma_median_to_mean() (see 02_single_arm.R)
# -----------------------------------------------------------------------------
m.cont <- metacont(n.e, mean.e, sd.e, n.c, mean.c, sd.c,
                   data = ma$OUTCOME,
                   studlab = Author,
                   sm = "MD",                  # "MD" same scale; "SMD" different scales
                   method.tau = "REML",
                   method.random.ci = "HK",
                   common = FALSE, random = TRUE)
summary(m.cont)
ma_forest(m.cont, lab.e = "INTERVENTION", lab.c = "CONTROL")

# -----------------------------------------------------------------------------
# C. GENERIC INVERSE VARIANCE — studies report only OR/RR/HR (+95% CI or SE)
# -----------------------------------------------------------------------------
# C1. Easiest: columns Author, HR, lower, upper  -> convert in R
d <- ma$OUTCOME
d <- cbind(d, ma_ci_to_se(d$HR, d$lower, d$upper, ratio = TRUE))   # adds TE, seTE
m.gen <- metagen(TE = TE, seTE = seTE,
                 data = d, studlab = Author,
                 sm = "HR",                    # "OR","RR","HR" (log scale) or "MD"
                 method.tau = "PM",
                 method.random.ci = "HK",
                 common = FALSE, random = TRUE)
# (equivalent shortcut: metagen(TE = log(HR), lower = log(lower), upper = log(upper), ...))
summary(m.gen)
ma_forest(m.gen, lab.e = "INTERVENTION", lab.c = "CONTROL")

# C2. Single study from 2x2 counts or OR + CI (quick calculator)
ma_counts_to_logOR(e.event = 17, n.e = 50, c.event = 8, n.c = 52)
ma_counts_to_logRR(e.event = 17, n.e = 50, c.event = 8, n.c = 52)
ma_ci_to_se(est = 1.52, lower = 0.85, upper = 2.71)

# -----------------------------------------------------------------------------
# D. SUBGROUPS — add a column (e.g. drug, design, RoB) and pass subgroup =
# -----------------------------------------------------------------------------
m.sub <- update(m.bin, subgroup = SUBGROUP_COLUMN)   # e.g. subgroup = drug
summary(m.sub)                                       # test for subgroup differences
ma_forest(m.sub, lab.e = "INTERVENTION", lab.c = "CONTROL")

# -----------------------------------------------------------------------------
# E. COMMON-EFFECT version (only if no heterogeneity / strong reason)
# -----------------------------------------------------------------------------
m.fixed <- update(m.bin, common = TRUE, random = FALSE)
ma_forest(m.fixed, random = FALSE)
