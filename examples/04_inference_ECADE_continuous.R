# Trial example 3: balancing continuous covariates directly (ECADE), while
# subsequent stratified analyses use their prespecified tertile categories.
# Run from repository root: source("examples/04_inference_ECADE_continuous.R")
source("R/01_designs.R")
source("R/02_inference.R")
source("R/03_trial.R")

set.seed(5001)
n <- 240L
baseline <- data.frame(x1 = rnorm(n), x2 = rnorm(n))

# Prospective cutpoints based on the known standard normal distribution.
# In a real trial, cutpoints must be specified before recruitment begins.
tertiles <- qnorm(c(1/3, 2/3))
baseline$x1_tertile <- car_discretize(baseline$x1, tertiles)
baseline$x2_tertile <- car_discretize(baseline$x2, tertiles)

prepared <- car_prepare(
  baseline,
  categorical = c("x1_tertile", "x2_tertile"),
  continuous = c("x1", "x2"),
  feature_cat = character(),
  feature_cont = c("x1", "x2"),
  analysis_covariates = c("x1", "x2")
)

# ECADE balances x1 and x2 in their original continuous representation.
plan <- car_design(prepared, method = "ECADE", rho = 0.85,
                   ecade_link = "Efron")
allocation <- car_allocate(plan, seed = 5002)

# A linear prognostic signal of reference SNR = 1 (standard normal errors).
set.seed(5003)
y <- 0.4 * allocation$treatment +
  (baseline$x1 + baseline$x2) / sqrt(2) + rnorm(n)

result <- car_analyze_trial(
  y, allocation$treatment, plan,
  L = 199, B = 49, seed = 5004
)
print(result)

# Comparison of inference strategies for ONE observed trial:
# - t, W: nine groups based on x1_tertile x x2_tertile;
# - treg, treg_HC0: regression on continuous x1, x2;
# - Dreg, Rreg: randomization tests using fixed residuals from x1, x2;
# - all randomization tests: ECADE allocation rule based on continuous x1, x2.
print(result$population[, c("test", "estimate", "se", "p_value")],
      row.names = FALSE)
print(result$rerandomization[, c("test", "p_value", "invalid")],
      row.names = FALSE)

# Type I error and statistical power are repeated-trial properties and cannot
# be estimated from this single example. Compare the tutorial's SNR study.
