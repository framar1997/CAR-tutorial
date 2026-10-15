# Convenience wrappers and formatted output for one completed trial.

car_analyze_trial <- function(y, treatment, design, L = 999L, B = 200L,
                              seed = 2026L, rerandomize = TRUE,
                              bootstrap = TRUE, include_diagnostic_hc0_fe = FALSE,
                              progress_every = 0L) {
  if (!inherits(design, "car_design")) stop("Use car_design() first")
  .car_check_trial(y, treatment, design$prepared)
  if (!is.null(seed)) set.seed(seed)
  # Rerandomization and bootstrap use separate deterministic streams if seed is set.
  boot <- if (bootstrap) car_bootstrap(
    y, treatment, design, B = B,
    progress_every = progress_every) else NULL
  population <- car_population(y, treatment, design$prepared, boot,
                               include_diagnostic_hc0_fe = include_diagnostic_hc0_fe)
  rr <- if (rerandomize) car_rerandomization(
    y, treatment, design, L = L,
    progress_every = progress_every) else NULL
  structure(list(method = design$method,
                 n = design$prepared$n,
                 group_sizes = c(control = sum(treatment == 0), treatment = sum(treatment == 1)),
                 population = population,
                 rerandomization = if (is.null(rr)) NULL else rr$tests,
                 bootstrap = if (is.null(boot)) NULL else
                   list(B = boot$B, invalid = boot$invalid, se = boot$se),
                 settings = list(L = if (rerandomize) L else 0L,
                                 B = if (bootstrap) B else 0L,
                                 seed = seed)), class = "car_trial_result")
}

print.car_trial_result <- function(x, ...) {
  cat("CAR single-trial analysis |", x$method, "| n =", x$n, "\n")
  cat("Control:", x$group_sizes[1], " Treatment:", x$group_sizes[2], "\n")
  cat("\nPopulation-based results:\n")
  print(x$population, row.names = FALSE)
  if (!is.null(x$rerandomization)) {
    cat("\nRe-randomization results (L =", x$settings$L, "):\n")
    print(x$rerandomization, row.names = FALSE)
  }
  invisible(x)
}

car_example_data <- function(n = 120L, seed = 123L) {
  set.seed(seed)
  data.frame(
    sex = factor(sample(c("F", "M"), n, replace = TRUE)),
    stage = factor(sample(c("low", "medium", "high"), n, replace = TRUE)),
    age_z = rnorm(n)
  )
}
