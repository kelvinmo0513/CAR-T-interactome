#!/usr/bin/env Rscript
# Selected manuscript models only. Input preparation performs all eligibility
# decisions; this stage fails rather than silently changing a pooling roster.
suppressPackageStartupMessages({
  library(survival)
  library(metafor)
  library(jsonlite)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript analyze.R INPUT_PREPARED_DIR OUTPUT_DIR")
inp <- normalizePath(args[1], mustWork = TRUE)
out <- normalizePath(args[2], mustWork = FALSE)
if (inp == out) stop("Input and output directories must differ")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
outputs <- c("model_statistics.csv", "pooled_statistics.csv", "pool_weights.csv",
             "KM_statistics.csv", "km_curves.csv", "risk_tables.csv",
             "R_validation.json", "R_sessionInfo.txt")
if (any(file.exists(file.path(out, outputs)))) stop("Refusing to overwrite existing analysis outputs")
read_input <- function(name, columns) {
  d <- read.csv(file.path(inp, name), stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(columns %in% names(d))) stop("Missing required columns in ", name)
  d
}
assert_unique <- function(d) {
  stopifnot(nrow(d) > 0, !anyNA(d$patient_id), !anyNA(d$source_id),
            !anyDuplicated(d$patient_id), !anyDuplicated(d$source_id))
}
check_fraction <- function(x, allow_missing = FALSE) {
  if (!allow_missing) stopifnot(!anyNA(x))
  stopifnot(all(is.finite(x[!is.na(x)])), all(x[!is.na(x)] >= 0),
            all(x[!is.na(x)] <= 1 + 1e-10))
}
capture_fit <- function(expr) {
  warnings <- character()
  error <- ""
  fit <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = function(e) { error <<- conditionMessage(e); NULL })
  list(fit = fit, warnings = warnings, error = error)
}
is_good <- function(r, term, logistic = FALSE) {
  if (is.null(r$fit) || nzchar(r$error)) return(FALSE)
  b <- coef(r$fit)[term]
  se <- sqrt(vcov(r$fit)[term, term])
  is.finite(b) && is.finite(se) && se > 0 &&
    !any(grepl("infinite|converge|numerically 0 or 1", r$warnings)) &&
    (!logistic || isTRUE(r$fit$converged))
}
failure_reason <- function(r) {
  paste(c("Non-estimable model", r$error, r$warnings)[
    nzchar(c("Non-estimable model", r$error, r$warnings))], collapse = "; ")
}
glm_control <- glm.control(epsilon = 1e-14, maxit = 100)
surv_data <- read_input("survival_analysis_data.csv", c("cohort", "endpoint", "marker",
  "source_id", "patient_id", "time_months", "event", "fraction", "high", "group"))
response_data <- read_input("response_analysis_data.csv", c("cohort", "endpoint", "marker",
  "source_id", "patient_id", "fraction", "bad"))
chemo_data <- read_input("chemo_analysis_data.csv", c("cohort", "source_id", "patient_id",
  "OS_Time", "OS_Status", "fraction"))
chemo_roster <- c("monti", "rchop", "chop")
response_roster <- c("moffitt", "zuma7", "stahl", "li_cart")
markers <- c("CCL8", "CCL13", "CCL8_CCL13", "CCL18")
stopifnot(setequal(unique(surv_data$cohort), "li_cart"),
          setequal(unique(surv_data$marker), markers),
          setequal(unique(chemo_data$cohort), chemo_roster),
          setequal(unique(response_data$cohort), response_roster),
          all(response_data$marker == "CCL8_CCL13"))

# Independently repeat all 24 absolute-score normalizations in R and verify the
# exact selected 22 fractions used below. Clinical exclusion decisions remain
# in preparation and are checked separately against the source records.
components <- c("CCL8+CCL13+ Myeloid DC", "Myeloid DC", "B cells", "T cells", "Tregs")
absolute <- read_input("cologne_absolute_scores_LOCAL_ONLY.csv", c("Mixture", components,
  "Absolute score (sig.score)"))
shares <- read_input("cologne_normalized_shares.csv", c("source_id", components))
stopifnot(nrow(absolute) == 24L, nrow(shares) == 24L,
          !anyNA(absolute$Mixture), !anyDuplicated(absolute$Mixture),
          !anyNA(shares$source_id), !anyDuplicated(shares$source_id),
          setequal(absolute$Mixture, shares$source_id))
score_matrix <- as.matrix(absolute[components])
stopifnot(all(is.finite(score_matrix)), all(score_matrix >= 0))
denominator <- rowSums(score_matrix)
stopifnot(all(denominator > 0),
          all(is.finite(absolute[["Absolute score (sig.score)"]])),
          max(abs(denominator - absolute[["Absolute score (sig.score)"]])) < 1e-7)
r_shares <- sweep(score_matrix, 1, denominator, "/")
python_shares <- as.matrix(shares[match(absolute$Mixture, shares$source_id), components])
normalization_difference <- max(abs(r_shares - python_shares))
stopifnot(all(is.finite(python_shares)), normalization_difference < 1e-12,
          max(abs(rowSums(r_shares) - 1)) < 1e-12)
cologne <- response_data[response_data$cohort == "stahl", ]
stopifnot(nrow(cologne) == 22L, !anyDuplicated(cologne$source_id),
          all(cologne$source_id %in% absolute$Mixture),
          max(abs(cologne$fraction - r_shares[match(cologne$source_id, absolute$Mixture), 1])) < 1e-12)
normalization_checks <- list(status = "PASS", available_samples = nrow(absolute),
  response_evaluable_samples = nrow(cologne), components = length(components),
  maximum_R_Python_difference = normalization_difference,
  denominator = "Sum of all retained absolute component scores", selected_fractions_match = TRUE,
  target_all_zero = all(r_shares[, 1] == 0), target_min = min(r_shares[, 1]),
  target_max = max(r_shares[, 1]), target_zero_count = sum(r_shares[, 1] == 0))

model_rows <- list()
scaling_checks <- list()
for (kind in c("survival", "response")) {
  roster <- if (kind == "survival") chemo_roster else response_roster
  for (co in roster) {
    z <- if (kind == "survival") chemo_data[chemo_data$cohort == co, ] else
      response_data[response_data$cohort == co, ]
    assert_unique(z)
    check_fraction(z$fraction, allow_missing = kind == "survival")
    joined_n <- nrow(z)
    # Kelvin's chemotherapy IQR is calculated BEFORE endpoint complete cases.
    # CAR-T inputs have already been restricted to endpoint-evaluable patients.
    q <- quantile(z$fraction, c(.25, .5, .75), type = 7, na.rm = kind == "survival")
    iqr <- unname(q[3] - q[1])
    if (kind == "survival") {
      complete <- complete.cases(z[c("OS_Time", "OS_Status", "fraction")])
      z <- z[complete, ]
      stopifnot(nrow(z) > 0, all(is.finite(z$OS_Time)), all(z$OS_Time >= 0),
                all(z$OS_Status %in% c(0, 1)))
      endpoint <- "OS"
    } else {
      stopifnot(!anyNA(z$bad), all(z$bad %in% c(0, 1)), length(unique(z$endpoint)) == 1L)
      endpoint <- z$endpoint[1]
    }
    b <- se <- pv <- ph_p <- NA_real_
    reason <- if (!is.finite(iqr) || iqr <= 0) "IQR = 0; per-IQR effect not estimable" else ""
    estimable <- FALSE
    warnings <- character()
    converged <- NA
    if (is.finite(iqr) && iqr > 0) {
      z$x <- (z$fraction - q[2]) / iqr
      r <- capture_fit(if (kind == "survival")
        coxph(Surv(OS_Time, OS_Status) ~ x, z, ties = "efron", x = TRUE, y = TRUE) else
        glm(bad ~ x, z, family = binomial(), control = glm_control))
      warnings <- r$warnings
      converged <- if (kind == "response" && !is.null(r$fit)) isTRUE(r$fit$converged) else NA
      estimable <- is_good(r, "x", logistic = kind == "response")
      if (estimable) {
        b <- unname(coef(r$fit)["x"])
        se <- sqrt(vcov(r$fit)["x", "x"])
        unscaled <- capture_fit(if (kind == "survival")
          coxph(Surv(OS_Time, OS_Status) ~ fraction, z, ties = "efron") else
          glm(bad ~ fraction, z, family = binomial(), control = glm_control))
        stopifnot(is_good(unscaled, "fraction", kind == "response"))
        delta_b <- abs(b - unname(coef(unscaled$fit)["fraction"]) * iqr)
        delta_se <- abs(se - sqrt(vcov(unscaled$fit)["fraction", "fraction"]) * iqr)
        stopifnot(delta_b < 1e-7, delta_se < 1e-6)
        z$minus_x <- -z$x
        reversed <- capture_fit(if (kind == "survival")
          coxph(Surv(OS_Time, OS_Status) ~ minus_x, z, ties = "efron") else
          glm(I(1 - bad) ~ x, z, family = binomial(), control = glm_control))
        reverse_term <- if (kind == "survival") "minus_x" else "x"
        stopifnot(is_good(reversed, reverse_term, kind == "response"),
                  abs(b + coef(reversed$fit)[reverse_term]) < 1e-7)
        scaling_checks[[paste(kind, co)]] <- list(cohort = co, beta_difference = delta_b,
          se_difference = delta_se, reverse_direction_verified = TRUE)
        pv <- 2 * pnorm(-abs(b / se))
        if (kind == "survival") {
          ph <- capture_fit(cox.zph(r$fit))
          if (!is.null(ph$fit)) ph_p <- ph$fit$table["x", "p"]
          warnings <- c(warnings, ph$warnings, if (nzchar(ph$error)) paste("PH diagnostic:", ph$error))
        }
      } else reason <- failure_reason(r)
    }
    model_rows[[paste(kind, co)]] <- data.frame(kind = kind, cohort = co, endpoint = endpoint,
      n = nrow(z), events = if (kind == "survival") sum(z$OS_Status) else sum(z$bad),
      q25 = unname(q[1]), median = unname(q[2]), q75 = unname(q[3]), IQR = iqr,
      zeros = sum(z$fraction == 0), estimable = estimable, reason = reason,
      beta = b, se = se, ratio = exp(b), CI_low = exp(b - qnorm(.975) * se),
      CI_high = exp(b + qnorm(.975) * se), p = pv, PH_test_p = ph_p,
      warnings = paste(warnings, collapse = "; "), joined_before_endpoint_omission = joined_n,
      IQR_basis = if (kind == "survival") "joined_before_endpoint_complete_cases" else
        "endpoint_evaluable", converged = converged)
  }
}
models <- do.call(rbind, model_rows)
rownames(models) <- NULL
write.csv(models, file.path(out, "model_statistics.csv"), row.names = FALSE, na = "")
if (any(!models$estimable)) {
  write_json(list(status = "FAIL", reason = "A required cohort was non-estimable; no pool produced",
    non_estimable = models$cohort[!models$estimable]), file.path(out, "R_validation.json"),
    auto_unbox = TRUE, pretty = TRUE)
  stop("Required cohort model is non-estimable; inspect model_statistics.csv")
}
pool_rows <- weight_rows <- list()
for (name in c("chemo_OS", "durability_primary")) {
  roster <- if (name == "chemo_OS") chemo_roster else response_roster
  m <- models[match(roster, models$cohort), ]
  stopifnot(identical(m$cohort, roster), all(m$estimable))
  fit_pool <- function(beta) {
    if (name == "chemo_OS") rma(yi = beta, sei = m$se, method = "REML") else
      rma.uni(yi = beta, sei = m$se, method = "REML", test = "z",
              control = list(threshold = 1e-12, maxiter = 10000))
  }
  fit <- fit_pool(m$beta)
  rev <- fit_pool(-m$beta)
  stopifnot(abs(as.numeric(fit$b) + as.numeric(rev$b)) < 1e-8,
            abs(fit$ci.lb + rev$ci.ub) < 1e-8)
  pool_rows[[name]] <- data.frame(pool = name, k = fit$k, n = sum(m$n),
    beta = as.numeric(fit$b), se = fit$se, ratio = exp(as.numeric(fit$b)),
    CI_low = exp(fit$ci.lb), CI_high = exp(fit$ci.ub), p = fit$pval,
    tau2 = fit$tau2, I2 = fit$I2)
  weight_rows[[name]] <- data.frame(pool = name, cohort = m$cohort, weight_percent = weights(fit))
}
write.csv(do.call(rbind, pool_rows), file.path(out, "pooled_statistics.csv"), row.names = FALSE)
write.csv(do.call(rbind, weight_rows), file.path(out, "pool_weights.csv"), row.names = FALSE)

# Recompute all four CAR-T KMs from the prepared total-cell fractions.
km_rows <- curves <- risks <- list()
for (marker in markers) {
  z <- surv_data[surv_data$marker == marker, ]
  assert_unique(z)
  check_fraction(z$fraction)
  stopifnot(all(z$endpoint == "EFS"), all(is.finite(z$time_months)), all(z$time_months >= 0),
    all(z$event %in% c(0, 1)), all(z$high == as.integer(z$fraction > median(z$fraction))),
    all(z$group == ifelse(z$high == 1, "High", "Low")), length(unique(z$high)) == 2L)
  r <- capture_fit(coxph(Surv(time_months, event) ~ high, z, ties = "efron"))
  ok <- is_good(r, "high")
  b <- if (ok) unname(coef(r$fit)["high"]) else NA_real_
  se <- if (ok) sqrt(vcov(r$fit)["high", "high"]) else NA_real_
  lr <- survdiff(Surv(time_months, event) ~ high, z, rho = 0)
  km_rows[[marker]] <- data.frame(cohort = "li_cart", marker = marker, endpoint = "EFS",
    n = nrow(z), events = sum(z$event), high_n = sum(z$high), low_n = sum(1 - z$high),
    median_fraction = median(z$fraction), ties_at_median = sum(z$fraction == median(z$fraction)),
    beta = b, se = se, HR = exp(b), CI_low = exp(b - qnorm(.975) * se),
    CI_high = exp(b + qnorm(.975) * se), cox_p = 2 * pnorm(-abs(b / se)),
    logrank_chisq = unname(lr$chisq), logrank_p = pchisq(lr$chisq, 1, lower.tail = FALSE),
    estimable = ok, warnings = paste(c(r$warnings, r$error)[nzchar(c(r$warnings, r$error))], collapse = "; "))
  ticks <- seq(0, 36, 12)
  if (max(z$time_months) >= 42) stop("Follow-up exceeds approved display range; inspect inputs")
  for (group in c("High", "Low")) {
    a <- z[z$group == group, ]
    sf <- survfit(Surv(time_months, event) ~ 1, a, conf.type = "log-log", conf.int = .95)
    curves[[paste(marker, group)]] <- data.frame(cohort = "li_cart", marker = marker,
      group = group, time = sf$time, survival = sf$surv, lower = sf$lower, upper = sf$upper,
      n_risk = sf$n.risk, n_event = sf$n.event, n_censor = sf$n.censor)
    rt <- summary(sf, times = ticks, extend = TRUE)
    direct <- vapply(ticks, function(t) sum(a$time_months >= t), integer(1))
    stopifnot(identical(as.integer(rt$n.risk), direct))
    risks[[paste(marker, group)]] <- data.frame(cohort = "li_cart", marker = marker,
      group = group, time = ticks, n_risk = direct)
  }
}
ks <- do.call(rbind, km_rows)
ks$logrank_q_BH4 <- p.adjust(ks$logrank_p, method = "BH")
write.csv(ks, file.path(out, "KM_statistics.csv"), row.names = FALSE, na = "")
write.csv(do.call(rbind, curves), file.path(out, "km_curves.csv"), row.names = FALSE, na = "")
write.csv(do.call(rbind, risks), file.path(out, "risk_tables.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out, "R_sessionInfo.txt"))
write_json(list(status = if (all(ks$estimable)) "PASS" else "REVIEW REQUIRED",
  predictor = "Total-mixture/total-cell fractions; per-IQR forests and median-split KMs",
  quantile_type = 7, chemo_IQR = "Before endpoint complete-case omission",
  CAR_T_IQR = "Endpoint-evaluable patients", Cox_ties = "Efron",
  confidence_intervals = "95% Wald for models; log-log Greenwood for KMs",
  pooling = "REML, normal Wald; original chemotherapy optimizer defaults preserved",
  non_estimable_KMs = ks$marker[!ks$estimable], scaling_checks = scaling_checks,
  KM_risk_count_checks = TRUE, median_tie_rule = "High > median; Low <= median",
  Cologne_reference = "oversampled two-cell AND reference; exploratory",
  Cologne_normalization_checks = normalization_checks,
  reference_uncertainty_in_model_CI = FALSE), file.path(out, "R_validation.json"),
  auto_unbox = TRUE, pretty = TRUE, digits = 16)
print(do.call(rbind, pool_rows), row.names = FALSE)
