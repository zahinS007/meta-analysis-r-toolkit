# =============================================================================
# 99_test_all.R — runs every analysis on synthetic data (smoke test)
# Rscript R/99_test_all.R   (run from the repository root)
# =============================================================================
source("R/ma_toolkit.R")
library(metafor); library(netmeta); library(mada); library(writexl)
set.seed(42)
od <- Sys.getenv("MA_TEST_OUT", file.path(tempdir(), "ma_test")); dir.create(od, showWarnings = FALSE)
old <- setwd(od); on.exit(setwd(old))

k <- 12
bin <- data.frame(Author = paste0("Study", 1:k), Year = 2010 + 1:k,
                  n.e = sample(80:400, k), n.c = sample(80:400, k))
bin$event.e <- rbinom(k, bin$n.e, 0.10); bin$event.c <- rbinom(k, bin$n.c, 0.15)
bin$drug <- rep(c("empa", "dapa", "cana"), length.out = k)
bin$pct_dm <- runif(k, 20, 80)
cont <- data.frame(Author = paste0("Study", 1:8), n.e = 50:57, mean.e = rnorm(8, 120, 4),
                   sd.e = runif(8, 8, 12), n.c = 52:59, mean.c = rnorm(8, 126, 4),
                   sd.c = runif(8, 8, 12))
gen <- data.frame(Author = paste0("Study", 1:6), HR = c(.8, .7, .9, .85, .6, 1.05))
gen$lower <- gen$HR * 0.75; gen$upper <- gen$HR / 0.75
prop <- data.frame(Author = paste0("S", 1:10), Year = 2011:2020,
                   n = sample(50:300, 10)); prop$event <- rbinom(10, prop$n, 0.08)
mean_df <- data.frame(Author = paste0("S", 1:6), n = c(100, 200, 150, 90, 300, 120),
                      mean = rnorm(6, 28, 1), sd = runif(6, 3, 5))
med <- data.frame(Author = c("A", "B", "C"), n = c(1390, 520, 200), median = c(38.5, 39.8, 40),
                  q1 = c(35.7, 36.9, NA), q3 = c(40.5, 42.1, NA),
                  min = c(30.5, NA, 31), max = c(43.8, NA, 47))
write_xlsx(list(acm = bin, sbp = cont, hr = gen, prev = prop, bmi = mean_df,
                bmi_med = med), "data.xlsx")

ma <- ma_import("data.xlsx")
ok <- function(label) cat("  [OK]", label, "\n")

# --- pairwise ---
m.bin <- metabin(event.e, n.e, event.c, n.c, data = ma$acm, studlab = Author, sm = "RR",
                 method = "MH", method.tau = "PM", method.random.ci = "HK",
                 common = FALSE, random = TRUE)
cat(ma_report(m.bin), "\n")
sz <- ma_forest_size(m.bin)
ma_save(ma_forest(m.bin, lab.e = "SGLT2i", lab.c = "Placebo", year_col = "Year"),
        "acm_forest", sz["width"], sz["height"]); ok("binary + forest")
m.sub <- update(m.bin, subgroup = drug)
ma_save(ma_forest(m.sub, lab.e = "SGLT2i", lab.c = "Placebo"), "acm_sub",
        12, ma_forest_size(m.sub)["height"]); ok("subgroup forest")
m.cont <- metacont(n.e, mean.e, sd.e, n.c, mean.c, sd.c, data = ma$sbp, studlab = Author,
                   sm = "MD", method.tau = "REML", method.random.ci = "HK",
                   common = FALSE, random = TRUE)
ma_save(ma_forest(m.cont), "sbp_forest", 13, 5); ok("continuous + forest")
d <- cbind(ma$hr, ma_ci_to_se(ma$hr$HR, ma$hr$lower, ma$hr$upper))
m.gen <- metagen(TE, seTE, data = d, studlab = Author, sm = "HR", method.tau = "PM",
                 method.random.ci = "HK", common = FALSE, random = TRUE)
ma_save(ma_forest(m.gen), "hr_forest", 10, 4); cat(ma_report(m.gen), "\n"); ok("generic HR")
print(ma_counts_to_logOR(17, 50, 8, 52)); print(ma_counts_to_logRR(17, 50, 8, 52))
m.fix <- update(m.bin, common = TRUE, random = FALSE)
ma_save(ma_forest(m.fix, random = FALSE), "acm_fixed", 12, 6); ok("common-effect forest")

# --- single arm ---
m.prop <- metaprop(event, n, data = ma$prev, studlab = Author, method = "GLMM",
                   sm = "PLOGIT", method.ci = "CP", common = FALSE, random = TRUE)
ma_save(ma_forest(m.prop, xlab = "Prevalence (%)", pscale = 100, year_col = "Year",
                  leftcols = c("studlab", "Year", "event", "n", "effect", "ci"),
                  leftlabs = c("Study", "Year", "Events", "Total", "Prevalence", "95% CI")),
        "prev_forest", 10, 5); cat(ma_report(m.prop), "\n"); ok("proportion")
print(ma_median_to_mean_df(ma$bmi_med)); ok("median -> mean")
m.mean <- metamean(n, mean, sd, data = ma$bmi, studlab = Author, sm = "MRAW",
                   method.tau = "REML", common = FALSE, random = TRUE)
ma_save(ma_forest(m.mean, xlab = "BMI"), "bmi_forest", 10, 4); ok("single mean")

# --- heterogeneity / bias ---
s <- ma_sensitivity(m.bin, prefix = "acm"); ok("sensitivity bundle")
mreg <- metareg(m.bin, ~ pct_dm); print(mreg)
ma_save(bubble(mreg, studlab = TRUE, bg = "red", xlab = "% with DM"), "acm_bubble", 8, 7)
ok("meta-regression")
rm <- rma(yi = m.bin$TE, sei = m.bin$seTE, method = "REML", slab = m.bin$studlab)
g <- gosh(rm, subsets = 2000); ma_save(plot(g, alpha = 0.1, col = "darkcyan"), "acm_gosh", 8, 8)
ok("GOSH")
ma_save(labbe(m.bin, studlab = TRUE, bg = "red"), "acm_labbe", 7, 7); ok("L'Abbe")
tf <- trimfill(m.bin); ok("trim-and-fill")

# --- NMA ---
arms <- data.frame(study = rep(paste0("T", 1:8), each = 2),
                   treatment = c("A","B","A","C","B","C","A","D","C","D","A","B","B","D","A","C"),
                   sampleSize = 100)
arms$responders <- rbinom(16, 100, c(A = .3, B = .22, C = .2, D = .25)[arms$treatment])
pw <- pairwise(treat = treatment, event = responders, n = sampleSize, studlab = study,
               data = arms, sm = "OR")
nma <- netmeta(pw, sm = "OR", common = FALSE, random = TRUE, reference.group = "A",
               method.tau = "REML", sep.trts = " vs ")
ma_save(netgraph(nma, thickness = "number.of.studies", number.of.studies = TRUE,
                 points = TRUE, col.points = "darkblue", cex.points = 3), "nma_net", 8, 8)
ma_save(meta::forest(nma, reference.group = "A"), "nma_forest", 9, 4)
netleague(nma, common = FALSE, random = TRUE, digits = 2, path = "league.xlsx", overwrite = TRUE)
print(netrank(nma, small.values = "good"))
rk <- rankogram(nma, nsim = 2000, small.values = "good"); ma_save(plot(rk), "nma_rank", 9, 6)
ns <- netsplit(nma); print(ns)
ma_save(meta::forest(ns), "nma_split", 9, 8); ok("network MA")

# --- DTA ---
dis <- sample(30:120, 10); hea <- sample(80:300, 10)
dta <- data.frame(study = paste0("D", 1:10), TP = rbinom(10, dis, .85),
                  FP = rbinom(10, hea, .2))
dta$FN <- dis - dta$TP; dta$TN <- hea - dta$FP
sens1 <- metaprop(TP, TP + FN, data = dta, studlab = study, method = "GLMM", sm = "PLOGIT")
spec1 <- metaprop(TN, TN + FP, data = dta, studlab = study, method = "GLMM", sm = "PLOGIT")
ma_save(meta::forest(sens1, xlim = c(0, 1)), "dta_sens", 10, 5)
fit1 <- reitsma(dta); print(summary(fit1)$coefficients); print(AUC(fit1))
ma_save({plot(fit1, xlim = c(0, 1), ylim = c(0, 1)); points(fpr(dta), sens(dta), pch = 2)},
        "dta_sroc", 7, 7)
dor <- metabin(TP, TP + FN, FP, FP + TN, data = dta, studlab = study, sm = "DOR",
               method = "Inverse", common = FALSE, random = TRUE)
print(metabias(dor, method.bias = "Deeks")); ok("diagnostic MA")

cat("\nAll outputs in:", od, "\n"); print(list.files(od))
