# Covariate-adaptive randomization: readable, base-R implementations.
# One call to car_allocate() draws ONE allocation, in recruitment order.
# For scientific descriptions, see README.md and the source Rcpp suite.

CAR_METHODS <- c("CR", "PBD", "RBD", "BSD", "CA-BCD", "HH", "PS",
                 "DA-BCD", "ECADE", "KER", "COV")

car_discretize <- function(x, cutpoints, labels = NULL) {
  if (!is.numeric(x) || anyNA(x)) stop("x must be a complete numeric vector")
  cp <- sort(unique(as.numeric(cutpoints)))
  if (any(!is.finite(cp))) stop("cutpoints must be finite")
  cut(x, breaks = c(-Inf, cp, Inf), labels = labels, include.lowest = TRUE)
}

car_prepare <- function(data, categorical = character(), continuous = character(),
                        feature_cat = categorical, feature_cont = continuous,
                        analysis_covariates = c(feature_cat, feature_cont),
                        analysis_strata = NULL) {
  if (!is.data.frame(data) || nrow(data) < 3L) stop("data must have >= 3 rows")
  cols <- unique(c(categorical, continuous, feature_cat, feature_cont,
                   analysis_covariates))
  if (anyDuplicated(names(data)) || !all(cols %in% names(data)))
    stop("Check covariate names and duplicate columns")
  if (length(intersect(categorical, continuous))) stop("Categorical and continuous names overlap")
  for (nm in unique(c(categorical, feature_cat))) {
    if (anyNA(data[[nm]])) stop("Missing categorical value: ", nm)
  }
  for (nm in unique(c(continuous, feature_cont))) {
    if (!is.numeric(data[[nm]]) || any(!is.finite(data[[nm]])))
      stop("Continuous variables must be numeric and finite: ", nm)
  }
  n <- nrow(data)
  cat_df <- if (length(categorical))
    as.data.frame(lapply(data[categorical], function(v) {
      if (is.factor(v)) factor(v, levels = levels(v)) else factor(v)
    }), optional = TRUE) else
    data.frame(row.names = seq_len(n))
  if (length(categorical)) {
    names(cat_df) <- categorical
    # Store factor levels so bootstrap resamples preserve the original support.
    for (nm in categorical) data[[nm]] <- cat_df[[nm]]
  }
  strata_design <- if (length(categorical)) {
    interaction(cat_df, drop = TRUE, lex.order = FALSE)
  } else factor(rep("all", n))
  if (is.null(analysis_strata)) analysis_strata <- strata_design
  if (length(analysis_strata) != n || anyNA(analysis_strata))
    stop("analysis_strata must be a complete vector with one value per participant")
  analysis_strata <- factor(analysis_strata)
  # Distinct levels in each categorical covariate get distinct marginal codes.
  marginal <- matrix(integer(0), nrow = n, ncol = 0L)
  if (length(categorical)) {
    m <- matrix(0L, n, length(categorical), dimnames = list(NULL, categorical))
    offset <- 0L
    for (j in seq_along(categorical)) {
      fj <- cat_df[[j]]
      m[, j] <- as.integer(fj) + offset
      offset <- offset + nlevels(fj)
    }
    marginal <- m
  }
  # The feature map can intentionally differ from the factors defining strata:
  # e.g. analysis tertiles but ORIGINAL continuous values in ECADE.
  model_frame <- data.frame(row.names = seq_len(n))
  for (nm in feature_cat) {
    v <- data[[nm]]
    model_frame[[nm]] <- if (is.factor(v)) factor(v, levels = levels(v)) else factor(v)
  }
  for (nm in feature_cont) model_frame[[nm]] <- data[[nm]]
  if (ncol(model_frame)) {
    X <- model.matrix(~ ., data = model_frame)[, -1L, drop = FALSE]
  } else X <- matrix(numeric(0), n, 0L)
  X <- as.matrix(X); storage.mode(X) <- "double"
  analysis_frame <- data.frame(row.names = seq_len(n))
  for (nm in analysis_covariates) {
    if (nm %in% categorical || is.factor(data[[nm]]) || is.character(data[[nm]])) {
      analysis_frame[[nm]] <- if (is.factor(data[[nm]]))
        factor(data[[nm]], levels = levels(data[[nm]])) else factor(data[[nm]])
    } else {
      analysis_frame[[nm]] <- data[[nm]]
    }
  }
  if (ncol(analysis_frame)) {
    Xa <- model.matrix(~ ., analysis_frame)[, -1L, drop = FALSE]
  } else Xa <- matrix(numeric(0), n, 0L)
  ker_feature_type <- integer(ncol(X))
  if (ncol(X)) {
    # model.matrix orders feature columns as named factor contrast columns and numerics;
    # derive categorical column flags from separate model.matrix on feature_cat.
    # model.matrix may reorder numerical vs factor features; assign by model term.
    terms_info <- attr(model.matrix(~ ., model_frame), "assign")[-1L]
    term_names <- attr(terms(~ ., data = model_frame), "term.labels")
    ker_feature_type <- as.integer(!(term_names[terms_info] %in% feature_cat))
  }
  list(n = n, data = data, design_strata = as.integer(strata_design),
       analysis_strata = as.integer(analysis_strata),
       marginal = marginal, design_X = cbind(`(Intercept)` = 1, X),
       feature_X = X, analysis_X = Xa,
       ker_type = ker_feature_type, categorical = categorical,
       continuous = continuous, feature_cat = feature_cat,
       feature_cont = feature_cont)
}

.car_prob <- function(score, rho) {
  if (score < 0) rho else if (score > 0) 1 - rho else 0.5
}
.car_sign_prob <- function(score, rho) {
  # The source implementation rounds near-zero scores to 1e-10.
  .car_prob(round(score, 10), rho)
}
.car_pinv_action <- function(G, x) {
  # Moore-Penrose inverse action for positive-semidefinite G.
  if (!length(x)) return(numeric())
  eg <- eigen(G, symmetric = TRUE)
  mx <- max(abs(eg$values))
  if (mx <= 0) return(rep(0, length(x)))
  keep <- eg$values > max(1, length(x)) * mx * .Machine$double.eps * 64
  if (!any(keep)) return(rep(0, length(x)))
  v <- eg$vectors[, keep, drop = FALSE]
  as.vector(v %*% (crossprod(v, x) / eg$values[keep]))
}
.car_prefix_actions <- function(X, method) {
  n <- nrow(X); p <- ncol(X)
  V <- matrix(0, n, p)
  G <- matrix(0, p, p)
  for (i in seq_len(n)) {
    xi <- X[i, ]
    if (method == "ECADE") G <- G + tcrossprod(xi)
    if (i > 1L) V[i, ] <- .car_pinv_action(G, xi)
    if (method == "DA-BCD") G <- G + tcrossprod(xi)
  }
  V
}
car_design <- function(prepared, method = "ECADE", block_size = 6L,
                       random_block_sizes = c(4L, 6L, 8L), mti = 3L,
                       a = 1, rho = 0.85, ps_rho = 0.75,
                       ker_rho = 0.8, cov_rho = 0.85,
                       hh_weights = NULL, cov_weights = c(1, 1, 1),
                       ecade_link = c("Efron", "Norm1", "Norm2"),
                       ker_burnin = 10L,
                       ker_legacy_scale = TRUE) {
  if (!is.list(prepared) || is.null(prepared$design_X)) stop("Use car_prepare() first")
  method <- match.arg(method, CAR_METHODS)
  ecade_link <- match.arg(ecade_link)
  if (any(c(block_size, random_block_sizes) < 2) ||
      any(c(block_size, random_block_sizes) %% 2 != 0)) stop("Block sizes must be even >= 2")
  if (mti < 1 || a < 0 || !is.finite(a)) stop("Invalid mti or a")
  if (any(!is.finite(c(rho, ps_rho, ker_rho, cov_rho))) ||
      any(c(rho, ps_rho, ker_rho, cov_rho) < .5) ||
      any(c(rho, ps_rho, ker_rho, cov_rho) > 1)) stop("rho must be in [0.5, 1]")
  if (length(cov_weights) != 3L || any(cov_weights < 0) || !any(cov_weights > 0))
    stop("Three nonnegative COV weights, at least one positive, required")
  if (ker_burnin %% 2 || ker_burnin < 2) stop("ker_burnin must be even >= 2")
  J <- ncol(prepared$marginal)
  if (is.null(hh_weights)) hh_weights <- rep(1 / (J + 2), J + 2)
  if (length(hh_weights) != J + 2 || any(hh_weights < 0) || sum(hh_weights) <= 0)
    stop("HH weights require J+2 nonnegative components: overall, J marginals, strata")
  hh_weights <- hh_weights / sum(hh_weights)
  if (method %in% c("PS", "HH") && J == 0)
    stop("PS and HH require categorical factors; discretize prospectively if necessary")
  X <- prepared$design_X
  if (method == "COV") {
    f <- prepared$feature_X
    k <- ncol(f); phi <- matrix(0, nrow(f), 1L + k + k * k)
    phi[, 1L] <- sqrt(cov_weights[1])
    if (k) {
      phi[, 1L + seq_len(k)] <- sqrt(cov_weights[2]) * f
      for (i in seq_len(nrow(f))) {
        phi[i, (k + 2L):ncol(phi)] <-
          sqrt(cov_weights[3]) * as.vector(t(tcrossprod(f[i, ])))
      }
    }
    X <- phi
  }
  V <- if (method %in% c("DA-BCD", "ECADE"))
    .car_prefix_actions(X, method) else NULL
  structure(list(prepared = prepared, method = method, X = X, V = V,
                 block_size = as.integer(block_size),
                 random_block_sizes = as.integer(random_block_sizes),
                 cov_weights = cov_weights, mti = as.integer(mti), a = a, rho = rho, ps_rho = ps_rho,
                 ker_rho = ker_rho, cov_rho = cov_rho,
                 hh_weights = hh_weights, ecade_link = ecade_link,
                 ker_burnin = as.integer(ker_burnin),
                 ker_legacy_scale = isTRUE(ker_legacy_scale)), class = "car_design")
}

car_allocate <- function(design, seed = NULL) {
  if (!inherits(design, "car_design")) stop("Use car_design()")
  if (!is.null(seed)) set.seed(seed)
  dat <- design$prepared; n <- dat$n; method <- design$method
  z <- dat$design_strata; K <- max(z)
  marg <- dat$marginal; J <- ncol(marg)
  D <- integer(K); N <- integer(K)
  global <- 0L
  margin_counts <- if (J) integer(max(marg)) else integer()
  remaining <- integer(K); remaining_A <- integer(K)
  delta <- integer(n); prob <- numeric(n)
  KX <- dat$feature_X; kt <- dat$ker_type
  # IMPORTANT: KER updates imbalance on feature_X (without the intercept),
  # while DA-BCD, ECADE, and COV update on design$X. The two dimensions
  # generally differ by one; using the wrong width silently recycles vectors.
  bt <- if (method == "KER") numeric(ncol(KX)) else numeric(ncol(design$X))
  nt <- 0L
  ker_sum <- rep(0, ncol(KX)); ker_sumsq <- ker_sum
  has_cont <- method == "KER" && any(kt == 1L)
  burn_left <- min(design$ker_burnin, n)
  burn_A <- burn_left %/% 2L
  for (i in seq_len(n)) {
    s <- z[i]; pp <- 0.5; score <- 0
    if (method %in% c("PBD", "RBD")) {
      if (remaining[s] == 0L) {
        bs <- if (method == "PBD") design$block_size else
          design$random_block_sizes[sample.int(length(design$random_block_sizes), 1L)]
        remaining[s] <- bs; remaining_A[s] <- bs %/% 2L
      }
      pp <- remaining_A[s] / remaining[s]
    } else if (method == "BSD") {
      if (D[s] >= design$mti) pp <- 0
      else if (D[s] <= -design$mti) pp <- 1
    } else if (method == "CA-BCD") {
      v <- abs(D[s]) ^ design$a
      pp <- if (D[s] > 0) 1 / (1 + v) else
        if (D[s] < 0) v / (1 + v) else 0.5
    } else if (method %in% c("HH", "PS")) {
      if (method == "HH") {
        ww <- design$hh_weights
        score <- ww[1L] * global + ww[J + 2L] * D[s] +
          sum(ww[1L + seq_len(J)] * margin_counts[marg[i, ]])
        if (abs(score) < 32 * .Machine$double.eps *
            (abs(global) + abs(D[s]) + sum(abs(margin_counts[marg[i, ]])))) score <- 0
        pp <- .car_prob(score, design$rho)
      } else {
        score <- sum(margin_counts[marg[i, ]])
        pp <- .car_prob(score, design$ps_rho)
      }
    } else if (method %in% c("DA-BCD", "ECADE") && i > 1L) {
      score <- sum(design$V[i, ] * bt)
      if (method == "DA-BCD") pp <- 0.5 - score / (1 + score^2)
      else {
        score <- 4 * i * score
        if (design$ecade_link == "Efron") pp <- .car_sign_prob(score, design$rho)
        else if (design$ecade_link == "Norm1") pp <- .1 + .8 * pnorm(-score)
        else pp <- .2 + .6 * pnorm(-score)
      }
    } else if (method == "KER") {
      if (has_cont && i <= burn_left) {
        pp <- burn_A / (burn_left - i + 1L)
      } else if (i > 1L) {
        for (j in seq_len(ncol(KX))) {
          if (kt[j] == 0L) {
            val <- KX[i, j]
            score <- score + if (val > 0.5) bt[j] / (i - 1L) else
              (global - bt[j]) / (i - 1L)
          } else if (i > 2L && nt > 0 && nt < i - 1L) {
            nn <- i - 1L
            sdj <- sqrt(max(0, (ker_sumsq[j] - ker_sum[j]^2 / nn) / (nn - 1L)))
            if (sdj > 0) {
              hT <- sdj / nt^0.2; hC <- sdj / (nn - nt)^0.2
              past <- seq_len(nn)
              # Legacy mode reproduces the Rcpp normalization (X/s) before KDE;
              # raw mode uses the usual dimensionally consistent (X_i-X_r)/h.
              diff <- KX[i, j] - KX[past, j]
              if (design$ker_legacy_scale) diff <- diff / sdj
              kernels <- ifelse(delta[past] == 1L,
                                dnorm(diff / hT) / hT, dnorm(diff / hC) / hC)
              score <- score + (sum(kernels[delta[past] == 1L]) -
                                sum(kernels[delta[past] == 0L])) / nn
            }
          }
        }
        pp <- .car_sign_prob(score, design$ker_rho)
      }
    } else if (method == "COV" && i > 1L) {
      pp <- .car_sign_prob(sum(bt * design$X[i, ]), design$cov_rho)
    }
    if (!is.finite(pp) || pp < -1e-10 || pp > 1 + 1e-10)
      stop("Invalid assignment probability at participant ", i)
    pp <- max(0, min(1, pp))
    prob[i] <- pp
    # Same Bernoulli assignment convention as the source engine.
    delta[i] <- as.integer(runif(1) < pp)
    change <- 2L * delta[i] - 1L
    D[s] <- D[s] + change; N[s] <- N[s] + 1L
    global <- global + change; nt <- nt + delta[i]
    if (method %in% c("PBD", "RBD")) {
      remaining[s] <- remaining[s] - 1L
      remaining_A[s] <- remaining_A[s] - delta[i]
    }
    if (method %in% c("PS", "HH"))
      for (j in seq_len(J)) margin_counts[marg[i, j]] <-
        margin_counts[marg[i, j]] + change
    if (method %in% c("DA-BCD", "ECADE", "COV")) bt <- bt + change * design$X[i, ]
    if (method == "KER") {
      if (length(bt) != ncol(KX))
        stop("Internal KER feature/imbalance dimension mismatch")
      bt <- bt + change * KX[i, ]
      if (has_cont && i <= burn_left) burn_A <- burn_A - delta[i]
      for (j in which(kt == 1L)) {
        ker_sum[j] <- ker_sum[j] + KX[i, j]
        ker_sumsq[j] <- ker_sumsq[j] + KX[i, j]^2
      }
    }
  }
  data.frame(patient = seq_len(n), treatment = delta,
             probability_T = prob, row.names = NULL)
}

car_bootstrap_design <- function(design, ids) {
  d <- design$prepared
  # Rebuild features in sampled patient order; precomputed DA/ECADE V changes.
  resampled <- car_prepare(d$data[ids, , drop = FALSE],
            categorical = d$categorical, continuous = d$continuous,
            feature_cat = d$feature_cat, feature_cont = d$feature_cont,
            analysis_covariates = character(),
            analysis_strata = d$analysis_strata[ids])
  car_design(resampled, method = design$method,
             block_size = design$block_size,
             random_block_sizes = design$random_block_sizes,
             mti = design$mti, a = design$a, rho = design$rho,
             ps_rho = design$ps_rho, ker_rho = design$ker_rho,
             cov_rho = design$cov_rho, hh_weights = design$hh_weights,
             ecade_link = design$ecade_link,
             cov_weights = design$cov_weights,
             ker_burnin = design$ker_burnin,
             ker_legacy_scale = design$ker_legacy_scale)
}
