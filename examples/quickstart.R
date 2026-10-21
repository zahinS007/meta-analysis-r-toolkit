# =============================================================================
# quickstart.R — a complete binary-outcome meta-analysis in ~20 lines
# Run from the repository root:  Rscript examples/quickstart.R
# Data are simulated (examples/data/example_data.xlsx), not from real trials.
# =============================================================================
source("R/ma_toolkit.R")
dir.create("output", showWarnings = FALSE)

ma <- ma_import("examples/data/example_data.xlsx")   # one sheet per outcome

m <- metabin(event.e, n.e, event.c, n.c, data = ma$acm, studlab = Author,
             sm = "RR", method = "MH", method.tau = "PM",
             method.random.ci = "HK", common = FALSE, random = TRUE)
summary(m)
cat("\n", ma_report(m), "\n\n")

sz <- ma_forest_size(m)
ma_save(ma_forest(m, lab.e = "Intervention", lab.c = "Control", year_col = "Year"),
        "output/forest_plot", sz["width"], sz["height"])

# Subgroups by drug, meta-regression on % with diabetes, sensitivity bundle
m.sub <- update(m, subgroup = drug)
ma_save(ma_forest(m.sub, lab.e = "Intervention", lab.c = "Control"),
        "output/forest_subgroups", 12, ma_forest_size(m.sub)["height"])
print(metareg(m, ~ pct_dm))
ma_sensitivity(m, prefix = "output/acm")
