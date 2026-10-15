# Trial example 2: Pocock-Simon minimization with marginal balancing.
# Run from repository root: source("examples/03_inference_PS.R")
source("R/01_designs.R")
source("R/02_inference.R")
source("R/03_trial.R")

set.seed(4001)
n <- 240L
baseline <- data.frame(
  sex = factor(sample(c("F", "M"), n, replace = TRUE)),
  stage = factor(sample(c("I", "II", "III"), n, replace = TRUE)),
  x1 = rnorm(n)
)

prepared <- car_prepare(
  baseline,
  categorical = c("sex", "stage"),
  continuous = "x1",
  analysis_covariates = c("sex", "stage", "x1")
)

# PS targets MARGINAL categorical balance rather than complete-stratum balance.
plan <- car_design(prepared, method = "PS", ps_rho = 0.75)
allocation <- car_allocate(plan, seed = 4002)

# Synthetic outcomes illustrate a heavier-tailed response distribution.
# In an actual study the recorded outcome vector replaces these draws.
set.seed(4003)
errors <- ifelse(runif(n) < 0.20, rt(n, df = 3), rnorm(n))
y <- 0.35 * allocation$treatment + 0.6 * baseline$x1 +
  0.25 * (baseline$stage == "III") + errors

result <- car_analyze_trial(
  y, allocation$treatment, plan,
  L = 199, B = 49, seed = 4004
)
print(result)

# The design-based analysis above replays PS, NOT stratified blocks:
# the assignment reference distribution must reflect marginal minimization.
# Its D/R statistics still use analysis strata defined by sex and stage.
print(result$population[, c("test", "p_value")], row.names = FALSE)
print(result$rerandomization[, c("test", "p_value", "invalid")],
      row.names = FALSE)
