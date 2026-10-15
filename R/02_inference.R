# Inference for ONE completed, two-arm trial. Base R only.
# A method that is undefined returns NA plus a note, never a fake zero p-value.

.car_check_trial <- function(y, treatment, prepared) {
  if (length(y) != prepared$n || any(!is.finite(y)))
    stop("y must be numeric, finite, and match the number of patients")
  if (length(treatment) != prepared$n || anyNA(treatment) ||
      !all(treatment %in% c(0, 1))) stop("treatment must contain 0 and 1 only")
  if (!any(treatment == 0) || !any(treatment == 1))
    stop("Both arms must contain at least one patient")
  invisible(TRUE)
}
.car_na_test <- function() c(estimate = NA_real_, se = NA_real_, p_value = NA_real_)
.car_two_sided <- function(x, df = Inf) {
  if (!is.finite(x) || !(df > 0)) return(NA_real_)
  if (is.infinite(df)) 2 * pnorm(-abs(x)) else 2 * pt(-abs(x), df)
}
.car_regression <- function(y, treatment, X) {
  n <- length(y)
  # Residualize y and treatment on the intercept and the analysis features.
  Q <- qr(cbind(`(Intercept)` = 1, X), tol = 1e-10)
  yd <- qr.resid(Q, y); dd <- qr.resid(Q, treatment)
  den <- sum(dd^2)
  if (den <= 1e-12) return(list(beta = NA_real_, se = NA_real_,
                                hc0 = NA_real_, df = NA_integer_))
  beta <- sum(dd * yd) / den
  e <- yd - beta * dd
  df <- n - Q$rank - 1L
  se <- if (df > 0) sqrt(sum(e^2) / (df * den)) else NA_real_
  hc0 <- sqrt(sum(dd^2 * e^2)) / den
  list(beta = beta, se = se, hc0 = hc0, df = df)
}
.car_strata_t <- function(y, d, z) {
  ns <- tabulate(z); k <- length(ns)
  nt <- tabulate(z[d == 1L], nbins = k); nc <- ns - nt
  if (any(nt[ns > 0] < 2L) || any(nc[ns > 0] < 2L)) return(.car_na_test())
  wt <- nt * nc / ns; wt <- wt / sum(wt)
  betas <- vT <- vC <- numeric(k)
  for (s in seq_len(k)) {
    yyT <- y[z == s & d == 1L]; yyC <- y[z == s & d == 0L]
    betas[s] <- mean(yyT) - mean(yyC)
    vT[s] <- var(yyT); vC[s] <- var(yyC)
  }
  est <- sum(wt * betas)
  se <- sqrt(sum(wt^2 * (vT / nt + vC / nc)))
  c(estimate = est, se = se,
    p_value = if (se > 0) .car_two_sided(est / se) else NA_real_)
}
.car_vanelteren <- function(y, d, z) {
  n <- length(y); k <- max(z)
  W <- mu <- varpart <- 0
  for (s in seq_len(k)) {
    ix <- which(z == s); nn <- length(ix)
    if (nn == 0L) next
    nt <- sum(d[ix]); nc <- nn - nt
    rs <- rank(y[ix], ties.method = "average")
    w <- 1 / (nn + 1)
    W <- W + w * sum(rs[d[ix] == 1])
    mu <- mu + nt / 2
    varpart <- varpart + nc * nt * w
  }
  se <- sqrt(varpart / 12)
  c(estimate = NA_real_, se = NA_real_,
    p_value = if (se > 0) .car_two_sided((W - mu) / se) else NA_real_)
}
.car_fe_hc0 <- function(y, d, z) {
  # Stratum fixed effects + HC0; retained as an exploratory legacy diagnostic.
  sums <- rowsum(cbind(n = 1, nt = d, sy = y), group = z, reorder = TRUE)
  occ <- as.integer(rownames(sums))
  ix <- match(z, occ)
  mean_d <- sums[ix, "nt"] / sums[ix, "n"]
  mean_y <- sums[ix, "sy"] / sums[ix, "n"]
  dr <- d - mean_d; yr <- y - mean_y
  den <- sum(dr^2)
  if (den <= 1e-12) return(.car_na_test())
  est <- sum(dr * yr) / den
  e <- yr - est * dr
  se <- sqrt(sum(dr^2 * e^2)) / den
  c(estimate = est, se = se,
    p_value = if (se > 0) .car_two_sided(est / se, length(y) - 1L) else NA_real_)
}
.car_bugni <- function(y, d, z) {
  n <- length(y); k <- max(z)
  ns <- tabulate(z, nbins = k)
  nt <- tabulate(z[d == 1], nbins = k); nc <- ns - nt
  if (any((ns > 0) & (nt == 0 | nc == 0))) return(.car_na_test())
  mT <- mean(y[d == 1]); mC <- mean(y[d == 0]); delta <- mT - mC
  first <- withinT <- withinC <- 0
  for (s in which(ns > 0)) {
    mt <- mean(y[z == s & d == 1]); mc <- mean(y[z == s & d == 0])
    first <- first + ns[s] * (mt - mc - delta)^2 / n
    withinT <- withinT + ns[s] * mt^2 / n
    withinC <- withinC + ns[s] * mc^2 / n
  }
  vT <- mean(y[d == 1]^2) - withinT
  vC <- mean(y[d == 0]^2) - withinC
  v <- (first + 2 * (vT + vC)) / n
  fe <- .car_fe_hc0(y, d, z)
  if (!is.finite(v) || v <= 0 || !is.finite(fe["estimate"])) return(.car_na_test())
  est <- unname(fe["estimate"]); se <- sqrt(v)
  c(estimate = est, se = se, p_value = .car_two_sided(est / se, n - 1L))
}

car_bootstrap <- function(y, treatment, design, B = 200L, seed = NULL,
                          progress_every = 0L) {
  if (!inherits(design, "car_design")) stop("design must be a car_design")
  .car_check_trial(y, treatment, design$prepared)
  if (B < 2L || B != as.integer(B)) stop("B must be an integer >= 2")
  if (!is.null(seed)) set.seed(seed)
  n <- length(y); draws <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    ix <- sample.int(n, n, replace = TRUE)
    boot_d <- car_bootstrap_design(design, ix)
    db <- car_allocate(boot_d)$treatment
    if (sum(db) > 0 && sum(db) < n)
      draws[b] <- mean(y[ix][db == 1]) - mean(y[ix][db == 0])
    if (progress_every > 0L && b %% progress_every == 0L)
      message("bootstrap: ", b, "/", B)
  }
  if (anyNA(draws)) {
    return(list(estimate = mean(y[treatment == 1]) - mean(y[treatment == 0]),
                se = NA_real_, p_value = NA_real_, invalid = sum(is.na(draws)),
                B = B, distribution = draws))
  }
  se <- sd(draws); est <- mean(y[treatment == 1]) - mean(y[treatment == 0])
  list(estimate = est, se = se,
       p_value = if (se > 0) .car_two_sided(est / se, n - 1L) else NA_real_,
       invalid = 0L, B = B, distribution = draws)
}

car_population <- function(y, treatment, prepared, bootstrap = NULL,
                           include_diagnostic_hc0_fe = FALSE) {
  .car_check_trial(y, treatment, prepared)
  d <- as.integer(treatment); z <- prepared$analysis_strata
  rg <- .car_regression(y, d, prepared$analysis_X)
  res <- list(
    t = .car_strata_t(y, d, z),
    W = .car_vanelteren(y, d, z),
    treg = c(estimate = rg$beta, se = rg$se,
             p_value = if (is.finite(rg$se) && rg$se > 0)
               .car_two_sided(rg$beta / rg$se, rg$df) else NA_real_),
    treg_HC0 = c(estimate = rg$beta, se = rg$hc0,
                 p_value = if (is.finite(rg$hc0) && rg$hc0 > 0)
                   .car_two_sided(rg$beta / rg$hc0, rg$df) else NA_real_),
    tadj = .car_bugni(y, d, z)
  )
  if (!is.null(bootstrap)) {
    res$tB <- c(estimate = bootstrap$estimate, se = bootstrap$se,
               p_value = bootstrap$p_value)
  }
  if (include_diagnostic_hc0_fe) res$FE_HC0_diagnostic <- .car_fe_hc0(y, d, z)
  ans <- data.frame(test = names(res), do.call(rbind, res), row.names = NULL)
  ans$interpretation <- ifelse(is.na(ans$p_value), "undefined in this trial", "computed")
  ans
}

.car_prepare_rr_stats <- function(y, z, X) {
  n <- length(y); k <- max(z)
  ns <- tabulate(z, nbins = k)
  # Rank statistic R: within-stratum centered ranks weighted by n_s/n.
  rank_score <- numeric(n)
  stratum_res <- numeric(n)
  for (s in seq_len(k)) {
    ix <- which(z == s)
    if (!length(ix)) next
    rank_score[ix] <- (rank(y[ix], ties.method = "average") - (length(ix) + 1) / 2) * length(ix) / n
    stratum_res[ix] <- y[ix] - mean(y[ix])
  }
  resid <- qr.resid(qr(cbind(`(Intercept)` = 1, X), tol = 1e-10), y)
  rank_res <- rank(resid, ties.method = "average") - (n + 1) / 2
  function(d) {
    nt <- sum(d); nc <- n - nt
    nT <- tabulate(z[d == 1], nbins = k)
    out <- c(D = NA_real_, R = sum(rank_score[d == 1]),
             Dreg = NA_real_, Rreg = sum(rank_res[d == 1]), FWL = NA_real_)
    if (nt > 0 && nc > 0) out["Dreg"] <- mean(resid[d == 1]) - mean(resid[d == 0])
    denom <- sum(nT * (ns - nT) / ns)
    if (denom > 0) out["FWL"] <- sum(stratum_res[d == 1]) / denom
    if (all(nT[ns > 0] > 0) && all(nT[ns > 0] < ns[ns > 0])) {
      total <- rowsum(y, z, reorder = TRUE)
      sumT <- rowsum(y * d, z, reorder = TRUE)
      # Both sums have a row for every observed stratum.
      mean_diff <- as.numeric(sumT) / nT - (as.numeric(total) - as.numeric(sumT)) / (ns - nT)
      out["D"] <- sum(ns * mean_diff) / n
    }
    out
  }
}

car_rerandomization <- function(y, treatment, design, L = 999L,
                                seed = NULL, undefined = c("strict", "conservative"),
                                tie_tolerance = 1e-12, progress_every = 0L,
                                keep_draws = FALSE) {
  if (!inherits(design, "car_design")) stop("Use car_design()")
  .car_check_trial(y, treatment, design$prepared)
  if (L < 2L || L != as.integer(L)) stop("L must be integer >= 2")
  undefined <- match.arg(undefined)
  if (!is.null(seed)) set.seed(seed)
  p <- design$prepared
  stats <- .car_prepare_rr_stats(y, p$analysis_strata, p$analysis_X)
  observed <- stats(as.integer(treatment))
  exceed <- invalid <- setNames(integer(length(observed)), names(observed))
  draws <- if (keep_draws) matrix(NA_real_, L - 1L, length(observed),
                                dimnames = list(NULL, names(observed))) else NULL
  for (b in seq_len(L - 1L)) {
    d <- car_allocate(design)$treatment
    sim <- stats(d)
    if (keep_draws) draws[b, ] <- sim
    bad <- !is.finite(sim)
    invalid <- invalid + as.integer(bad)
    valid <- !bad & is.finite(observed)
    exceed[valid] <- exceed[valid] + as.integer(
      abs(sim[valid]) + tie_tolerance * pmax(1, abs(observed[valid])) >= abs(observed[valid]))
    if (progress_every > 0L && b %% progress_every == 0L)
      message("re-randomization: ", b, "/", L - 1L)
  }
  numerator <- 1 + exceed + if (undefined == "conservative") invalid else 0
  pv <- numerator / L
  pv[!is.finite(observed) | (undefined == "strict" & invalid > 0)] <- NA_real_
  ans <- data.frame(test = names(observed), observed_statistic = as.numeric(observed),
                    extreme = as.integer(exceed), invalid = as.integer(invalid),
                    p_value = as.numeric(pv), L = L,
                    row.names = NULL)
  list(tests = ans, simulated_statistics = draws, undefined = undefined)
}
