# Trial example 1: stratified permuted blocks and subsequent inference.
# Run from repository root: source("examples/02_inference_PBD.R")
source("R/01_designs.R")
source("R/02_inference.R")
source("R/03_trial.R")

set.seed(3001)
n <- 240L
baseline <- data.frame(
  sex = factor(sample(c("F", "M"), n, replace = TRUE)),
  stage = factor(sample(c("I", "II", "III"), n, replace = TRUE)),
  age_z = rnorm(n)
)

prepared <- car_prepare(
  baseline,
  categorical = c("sex", "stage"),
  continuous = "age_z",
  analysis_covariates = c("sex", "stage", "age_z")
)

# Randomize within the 2-by-3 strata using fixed blocks of size six.
plan <- car_design(prepared, method = "PBD", block_size = 6)
allocation <- car_allocate(plan, seed = 3002)

# Synthetic outcomes for illustration. In a real trial use the observed data.
set.seed(3003)
y <- 0.4 * allocation$treatment + 0.5 * baseline$age_z +
  0.35 * (baseline$stage == "III") + rnorm(n)

# Analyze this completed trial. B bootstrap resamples and L-1 new
# allocations approximate the reference distributions; they are not trials.
result <- car_analyze_trial(
  y, allocation$treatment, plan,
  L = 199, B = 49, seed = 3004
)
print(result)

# The observed patient-level trial record can be assembled as follows:
trial_record <- cbind(baseline, allocation, outcome = y)
print(head(trial_record))

# Notes: 't', 'W', and 'D' use the 6 categorical strata, while regression
# methods also adjust for age_z. 'tB' uses patient resampling and allocation
# under the original PBD rule. Sharp-null randomization p-values require
# reproducing the original allocation method and block size.
