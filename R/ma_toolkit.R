# =============================================================================
# ma_toolkit.R — helper functions for meta-analysis in R (`meta` package)
# Author: Zahin Shahriar (https://github.com/zahinS007) · MIT License
# Requires meta >= 7 (tested with meta 8.5, metafor 5.2, netmeta 3.7)
# Usage:  source("R/ma_toolkit.R")
# =============================================================================

suppressPackageStartupMessages({
  library(meta)
  library(readxl)
})

# -----------------------------------------------------------------------------
# 1. IMPORT — one Excel workbook, one sheet per outcome (acm, mb, ...)
#    Column names: no spaces; use event.e, n.e, event.c, n.c, mean.e, sd.e ...
# -----------------------------------------------------------------------------
ma_import <- function(path) {
  sheets <- readxl::excel_sheets(path)
  out <- lapply(sheets, function(s) as.data.frame(readxl::read_excel(path, sheet = s)))
  names(out) <- sheets
  for (s in sheets) {
    message("Sheet '", s, "': ", nrow(out[[s]]), " studies; columns: ",
            paste(names(out[[s]]), collapse = ", "))
  }
  out
}

# -----------------------------------------------------------------------------
# 2. SAVE ANY PLOT — wraps png/pdf/svg + dev.off() so nothing gets lost
#    Forest plots are not auto-rescaled: raise width/height if text is cut off.
#    ma_save(forest(m), "acm_forest", width = 12, height = 6, format = "png")
# -----------------------------------------------------------------------------
ma_save <- function(plot_expr, file, width = 12, height = 6,
                    format = c("png", "pdf", "svg", "tiff"), res = 300) {
  format <- match.arg(format)
  fn <- paste0(sub(paste0("\\.", format, "$"), "", file), ".", format)
  switch(format,
    png  = png(fn,  width = width, height = height, units = "in", res = res),
    tiff = tiff(fn, width = width, height = height, units = "in", res = res,
                compression = "lzw"),
    pdf  = pdf(fn,  width = width, height = height),
    svg  = svg(fn,  width = width, height = height)   # svg() is in INCHES
  )
  on.exit(dev.off())
  force(plot_expr)
  message("Saved: ", normalizePath(fn, mustWork = FALSE))
  invisible(fn)
}

# Rough size so a forest plot is not cut off (inches)
ma_forest_size <- function(m, base_h = 2.5, per_row = 0.28, width = 12) {
  k <- m$k.all + if (!is.null(m$subgroup)) 3 * length(unique(m$subgroup)) else 0
  c(width = width, height = base_h + per_row * k)
}

# -----------------------------------------------------------------------------
# 3. ACADEMY-STYLE FOREST PLOT (RevMan 5 look) for any meta object
#    type is detected from the object: metabin / metacont / metagen /
#    metaprop / metamean. Pass extra forest() args via ...
# -----------------------------------------------------------------------------
ma_forest <- function(m,
                      lab.e = "Experimental", lab.c = "Control",
                      label.left = paste("Favors", lab.e),
                      label.right = paste("Favors", lab.c),
                      xlab = NULL,
                      random = TRUE,
                      year_col = NULL,          # e.g. "Year" if in data
                      layout = "RevMan5",
                      sortvar = NULL,
                      ...) {
  w <- if (random) "w.random" else "w.common"
  sm <- m$sm
  cls <- class(m)[1]

  yr <- if (!is.null(year_col)) year_col else character(0)
  yrlab <- if (!is.null(year_col)) "Year" else character(0)

  cols <- switch(cls,
    metabin  = list(c("studlab", yr, "event.e", "n.e", "event.c", "n.c", w, "effect", "ci"),
                    c("Study", yrlab, "Events", "Total", "Events", "Total", "Weight", sm, "95% CI")),
    metacont = list(c("studlab", yr, "mean.e", "sd.e", "n.e", "mean.c", "sd.c", "n.c", w, "effect", "ci"),
                    c("Study", yrlab, "Mean", "SD", "Total", "Mean", "SD", "Total", "Weight", sm, "95% CI")),
    metaprop = list(c("studlab", yr, "event", "n", w, "effect", "ci"),
                    c("Study", yrlab, "Events", "Total", "Weight", "Proportion", "95% CI")),
    metamean = list(c("studlab", yr, "n", "mean", "sd", w, "effect", "ci"),
                    c("Study", yrlab, "Total", "Mean", "SD", "Weight", "Mean", "95% CI")),
    # metagen (generic inverse variance: OR/RR/HR/MD from log effect + SE or CI)
    list(c("studlab", yr, "effect", "ci", w),
         c("Study", yrlab, sm, "95% CI", "Weight"))
  )

  args <- list(
    x = m,
    layout = layout,
    sortvar = if (is.null(sortvar)) m$studlab else sortvar,
    leftcols = cols[[1]], leftlabs = cols[[2]],
    rightcols = FALSE,
    just = "center",
    random = random, common = !random,
    pooled.events = TRUE, pooled.totals = TRUE,
    test.overall = TRUE,
    overall.hetstat = TRUE,
    print.Q = TRUE, print.pval.Q = TRUE, print.tau.ci = TRUE,
    digits = 2, digits.pval = 3, digits.sd = 2,
    col.square = "darkblue", col.square.lines = "black",
    col.diamond = "black", col.diamond.lines = "black",
    fs.heading = 12, fs.study = 11, fs.hetstat = 10,
    colgap = "5mm", colgap.forest = "5mm"
  )
  # two-arm plots get arm headers + favours labels; single-arm get an x label
  if (cls %in% c("metabin", "metacont", "metagen")) {
    # meta >= 8 renamed lab.e/lab.c to label.e/label.c
    new_names <- "label.e" %in% names(formals(meta:::forest.meta))
    arm <- if (new_names) list(label.e = lab.e, label.c = lab.c)
           else list(lab.e = lab.e, lab.c = lab.c)
    args <- c(args, arm, list(label.left = label.left, label.right = label.right,
                              ff.lr = "bold"))
  } else {
    args$col.diamond <- "maroon"
    if (!is.null(xlab)) args <- c(args, list(xlab = xlab, ff.xlab = "bold", fs.xlab = 12))
  }
  # subgroup display if the model has subgroups
  if (!is.null(m$subgroup)) {
    args <- c(args, list(subgroup = TRUE, print.subgroup.name = FALSE,
                         col.subgroup = "black", test.subgroup = TRUE,
                         test.effect.subgroup = TRUE))
  }
  extra <- list(...)
  args[names(extra)] <- extra
  do.call(meta::forest, args)
}

# -----------------------------------------------------------------------------
# 4. EFFECT-SIZE CONVERSIONS (for metagen when studies give only OR/RR/HR + CI)
# -----------------------------------------------------------------------------
# Ratio measure (OR/RR/HR) + 95% CI  ->  log effect (TE) and SE
ma_ci_to_se <- function(est, lower, upper, ratio = TRUE, level = 0.95) {
  z <- qnorm(1 - (1 - level) / 2)
  if (ratio) {
    data.frame(TE = log(est), seTE = (log(upper) - log(lower)) / (2 * z))
  } else {                                   # MD / difference measures
    data.frame(TE = est, seTE = (upper - lower) / (2 * z))
  }
}

# 2x2 counts -> log OR / log RR and SE (a = events exp, c = events ctrl)
ma_counts_to_logOR <- function(e.event, n.e, c.event, n.c) {
  a <- e.event; b <- n.e - e.event; c <- c.event; d <- n.c - c.event
  data.frame(TE = log((a * d) / (b * c)), seTE = sqrt(1/a + 1/b + 1/c + 1/d))
}
ma_counts_to_logRR <- function(e.event, n.e, c.event, n.c) {
  # correct SE for log RR (the course page used the OR formula by mistake)
  data.frame(TE = log((e.event / n.e) / (c.event / n.c)),
             seTE = sqrt(1/e.event - 1/n.e + 1/c.event - 1/n.c))
}

# -----------------------------------------------------------------------------
# 5. MEDIAN / IQR / RANGE -> MEAN (SD)   (estmeansd package)
#    method = "qe" (quantile estimation, default), "bc" (Box-Cox), "mln"
#    Give whichever you have: min/max (S1), q1/q3 (S2), or all five (S3).
# -----------------------------------------------------------------------------
ma_median_to_mean <- function(n, median, q1 = NA, q3 = NA, min = NA, max = NA,
                              method = c("qe", "bc", "mln")) {
  method <- match.arg(method)
  if (!requireNamespace("estmeansd", quietly = TRUE))
    stop("install.packages('estmeansd') first")
  f <- switch(method, qe = estmeansd::qe.mean.sd,
                      bc = estmeansd::bc.mean.sd,
                      mln = estmeansd::mln.mean.sd)
  args <- list(n = n, med.val = median)
  if (!is.na(q1))  args$q1.val  <- q1
  if (!is.na(q3))  args$q3.val  <- q3
  if (!is.na(min)) args$min.val <- min
  if (!is.na(max)) args$max.val <- max
  r <- do.call(f, args)
  c(mean = unname(r$est.mean), sd = unname(r$est.sd))
}

# Vectorised over a data frame with columns n, median, q1, q3, min, max
ma_median_to_mean_df <- function(df, method = "qe") {
  get <- function(col, i) if (col %in% names(df)) df[[col]][i] else NA
  res <- t(sapply(seq_len(nrow(df)), function(i)
    ma_median_to_mean(get("n", i), get("median", i), get("q1", i), get("q3", i),
                      get("min", i), get("max", i), method = method)))
  cbind(df, est.mean = res[, "mean"], est.sd = res[, "sd"])
}

# -----------------------------------------------------------------------------
# 5b. CLEAN BAUJAT AND SROC PLOTS
#     ma_baujat(): labels only the studies that stand out, so labels never pile
#                  up at the origin. Change the cut-offs with x_cut / y_cut.
#     ma_sroc():   bivariate (Reitsma) SROC with study points sized by n,
#                  confidence + prediction regions, legend and pooled estimates.
# -----------------------------------------------------------------------------
ma_baujat <- function(m, x_cut = NULL, y_cut = NULL, bg = "firebrick", ...) {
  op <- par(cex.lab = 1.1, font.lab = 2); on.exit(par(op))
  bj <- meta::baujat(m, bg = bg, cex = 1.3, studlab = FALSE, ...)
  if (is.null(x_cut)) x_cut <- max(2, 0.25 * max(bj$x))
  if (is.null(y_cut)) y_cut <- 0.25 * max(bj$y)
  flag <- bj$x > x_cut | bj$y > y_cut
  if (any(flag)) text(bj$x[flag], bj$y[flag], bj$studlab[flag],
                      pos = ifelse(bj$x[flag] > 0.6 * max(bj$x), 2, 4),
                      cex = 0.9, font = 2)
  if (any(!flag)) legend("topleft", bty = "n", cex = 0.85,
                         legend = "Unlabelled points: studies with little impact")
  invisible(bj)
}

ma_sroc <- function(fit, data, main = "Summary ROC curve (bivariate model)",
                    col = "steelblue") {
  s   <- summary(fit)$coefficients
  se  <- s["sensitivity", c(1, 5, 6)]
  fp  <- s["false pos. rate", c(1, 5, 6)]
  auc <- mada::AUC(fit)$AUC
  n   <- with(data, TP + FN + FP + TN)
  op  <- par(mar = c(4.5, 4.5, 3, 1), las = 1); on.exit(par(op))
  plot(fit, xlim = c(0, 1), ylim = c(0, 1), sroclwd = 2, predict = TRUE, ann = FALSE)
  title(xlab = "1 - Specificity (false positive rate)", ylab = "Sensitivity", main = main)
  abline(0, 1, lty = 3, col = "gray60")
  points(mada::fpr(data), mada::sens(data), pch = 21, bg = col, col = "white",
         cex = 0.6 + 1.6 * sqrt(n / max(n)))
  legend("bottomright", bty = "n", cex = 0.85,
         legend = c("Individual studies (size ~ n)", "Summary estimate",
                    "95% confidence region", "95% prediction region", "SROC curve"),
         pch = c(21, 1, NA, NA, NA), pt.bg = c(col, NA, NA, NA, NA),
         col = c("white", "black", "black", "black", "black"),
         lty = c(NA, NA, 1, 3, 1), lwd = c(NA, NA, 1, 1, 2))
  text(0.97, 0.42, adj = 1, cex = 0.85, sprintf(
    "Sensitivity %.2f (%.2f-%.2f)\nSpecificity %.2f (%.2f-%.2f)\nAUC %.2f",
    se[1], se[2], se[3], 1 - fp[1], 1 - fp[3], 1 - fp[2], auc))
  invisible(fit)
}

# -----------------------------------------------------------------------------
# 6. SENSITIVITY BUNDLE — leave-one-out, Baujat, influence, funnel + Egger
#    ma_sensitivity(m, prefix = "acm")  -> saves PNGs and prints results
# -----------------------------------------------------------------------------
ma_sensitivity <- function(m, prefix = "outcome", save = TRUE, random = TRUE) {
  out <- list()
  # Leave-one-out
  out$loo <- meta::metainf(m, pooled = if (random) "random" else "common")
  print(out$loo)
  loo_plot <- function() meta::forest(out$loo, rightcols = c("effect", "ci", "I2"),
                                      col.diamond = "gray")
  if (save) ma_save(loo_plot(), paste0(prefix, "_leave1out"), width = 10,
                    height = 2 + 0.3 * m$k) else loo_plot()
  # Baujat
  bj <- function() ma_baujat(m)
  if (save) ma_save(bj(), paste0(prefix, "_baujat"), 8, 7) else bj()
  # Funnel + Egger (only meaningful with >= 10 studies)
  fp <- function() meta::funnel(m, studlab = TRUE, bg = "red",
                                common = !random, random = random)
  if (save) ma_save(fp(), paste0(prefix, "_funnel"), 8, 7) else fp()
  if (m$k >= 10) {
    out$egger <- meta::metabias(m, method.bias = "linreg")
    print(out$egger)
  } else {
    message("Egger's test skipped: k = ", m$k, " (< 10 studies, test underpowered).")
  }
  # Influence diagnostics via metafor (no dmetar needed)
  if (requireNamespace("metafor", quietly = TRUE)) {
    rma_m <- metafor::rma(yi = m$TE, sei = m$seTE,
                          method = if (m$method.tau == "DL") "DL" else m$method.tau,
                          slab = m$studlab)
    out$influence <- stats::influence(rma_m)
    print(out$influence)
    ip <- function() plot(out$influence)
    if (save) ma_save(ip(), paste0(prefix, "_influence"), 10, 9) else ip()
  }
  invisible(out)
}

# -----------------------------------------------------------------------------
# 7. ONE-LINE RESULTS SENTENCE for the manuscript
# -----------------------------------------------------------------------------
ma_report <- function(m, random = TRUE, digits = 2) {
  te  <- if (random) m$TE.random else m$TE.common
  lo  <- if (random) m$lower.random else m$lower.common
  up  <- if (random) m$upper.random else m$upper.common
  p   <- if (random) m$pval.random else m$pval.common
  bt  <- function(x) if (m$sm %in% c("OR", "RR", "HR", "IRR", "DOR", "ROM")) exp(x) else x
  if (inherits(m, "metaprop") || inherits(m, "metamean")) {
    bt <- function(x) meta:::backtransf(x, m$sm, n = mean(m$n))
  }
  fmtp <- if (is.null(p) || is.na(p)) "no test vs null" else
          if (p < 0.001) "p<0.001" else paste0("p=", formatC(p, digits = 3, format = "f"))
  lab <- if (inherits(m, "metaprop")) "Pooled proportion" else
         if (inherits(m, "metamean")) "Pooled mean" else m$sm
  sprintf("%s %s (95%% CI %s to %s; %s; I²=%s%%; %d studies, %s model)",
          lab, formatC(bt(te), digits = digits, format = "f"),
          formatC(bt(lo), digits = digits, format = "f"),
          formatC(bt(up), digits = digits, format = "f"),
          fmtp, round(100 * m$I2), m$k, if (random) "random-effects" else "common-effect")
}
