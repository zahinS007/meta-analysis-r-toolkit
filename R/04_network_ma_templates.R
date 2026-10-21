# =============================================================================
# 04_network_ma_templates.R — network meta-analysis (frequentist netmeta;
# Bayesian gemtc optional). Check transitivity (similar populations across
# comparisons) BEFORE pooling; check consistency (netsplit) AFTER.
# =============================================================================
library(netmeta); library(meta); library(readxl); library(writexl)
source("R/ma_toolkit.R")

# Arm-based sheet, one row per arm:
#   study | treatment | responders | sampleSize        (binary)
#   study | treatment | mean | sd | n                  (continuous)
arms <- as.data.frame(read_excel("data.xlsx", sheet = "NMA"))

# 1. Arm-based -> contrast-based
pw <- pairwise(treat = treatment, event = responders, n = sampleSize,
               studlab = study, data = arms, sm = "OR")
# continuous: pairwise(treat = treatment, mean = mean, sd = sd, n = n, studlab = study, data = arms, sm = "MD")
write_xlsx(as.data.frame(pw), "nma_pairwise.xlsx")

# 2. Fit
nma <- netmeta(pw, sm = "OR",
               common = FALSE, random = TRUE,
               reference.group = "REFERENCE_TREATMENT",   # usually placebo/control
               method.tau = "REML",
               details.chkmultiarm = TRUE, sep.trts = " vs ")
summary(nma)

# 3. Network graph (edge thickness = no. of studies)
ma_save(netgraph(nma, thickness = "number.of.studies", number.of.studies = TRUE,
                 points = TRUE, col.points = "darkblue", cex.points = 3,
                 plastic = FALSE, multiarm = TRUE),
        "nma_network", 8, 8)

# 4. Forest vs reference
ma_save(meta::forest(nma, reference.group = "REFERENCE_TREATMENT",
               sortvar = TE, smlab = "Odds ratio vs reference"),
        "nma_forest", 9, 2.5 + 0.35 * nma$n)

# 5. League table (Excel)
netleague(nma, common = FALSE, random = TRUE, backtransf = TRUE, digits = 2,
          path = "nma_league.xlsx", overwrite = TRUE)

# 6. Ranking — P-scores (small.values = "good" if lower OR/event = better)
netrank(nma, small.values = "good")
rk <- rankogram(nma, nsim = 10000, random = TRUE, small.values = "good")
ma_save(plot(rk), "nma_rankogram", 9, 6)

# 7. Consistency: design decomposition + node-splitting
decomp.design(nma)
ns <- netsplit(nma, random = TRUE, common = FALSE)
ns
ma_save(meta::forest(ns), "nma_netsplit", 9, 3 + 0.4 * nrow(ns$random))
netheat(nma, random = TRUE)

# 8. Comparison-adjusted funnel (order = treatments from newest to oldest)
ma_save(funnel(nma, order = c("TRT_NEW", "TRT_OLD", "REFERENCE_TREATMENT"),
               linreg = TRUE, pch = 19),
        "nma_funnel", 8, 7)

# -----------------------------------------------------------------------------
# OPTIONAL Bayesian NMA (gemtc + JAGS; needs JAGS installed on the machine)
# Data columns must be named: study, treatment, responders, sampleSize
# -----------------------------------------------------------------------------
# library(gemtc); library(rjags)
# net  <- mtc.network(data.ab = arms)
# mod  <- mtc.model(net, type = "consistency", likelihood = "binom",
#                   link = "logit", linearModel = "random", n.chain = 4)
# fit  <- mtc.run(mod, n.adapt = 5000, n.iter = 1e5, thin = 10)
# gelman.diag(fit)              # PSRF ~1 = converged
# round(exp(relative.effect.table(fit)), 2)        # league table (OR)
# rank.probability(fit, preferredDirection = -1)   # -1 if lower = better
# ns_b <- mtc.nodesplit(net, likelihood = "binom", link = "logit",
#                       linearModel = "random"); summary(ns_b)
