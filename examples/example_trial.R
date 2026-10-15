# Run this from the GitHub repository root (no packages required).
source("R/01_designs.R")
source("R/02_inference.R")
source("R/03_trial.R")

# 1) Baseline covariates listed in the original recruitment order.
baseline <- car_example_data(n = 120, seed = 1)

# 2) Specify analysis features before randomization.
#    Here both categorical factors and age are prognostic.
prepared <- car_prepare(baseline,
                        categorical = c("sex", "stage"),
                        continuous = "age_z")

# 3) ONE trial randomization. Change "ECADE" to any CAR_METHODS name.
design <- car_design(prepared, method = "ECADE", rho = 0.85)
allocation <- car_allocate(design, seed = 12345)

# 4) Trial responses: demonstration only. For a completed clinical trial,
#    replace this line with the OBSERVED, patient-level outcome vector.
set.seed(44)
y <- 0.4 * allocation$treatment + 0.5 * baseline$age_z +
     0.3 * (baseline$stage == "high") + rnorm(nrow(baseline))

# 5) Analyze ONE completed trial. The L-1 allocations and B bootstrap
#    samples below form reference distributions, NOT new simulated trials.
result <- car_analyze_trial(y, allocation$treatment, design,
                            L = 499, B = 99, seed = 900)
print(result)

# For a production run, increase L (e.g. L=20000) and B (e.g. B=500).
# Time depends strongly on the method. Continuous KER is especially slow in R.

# An optional plain CSV containing the original allocation and responses:
trial_data <- cbind(baseline, allocation[c("treatment", "probability_T")], y = y)
# write.csv(trial_data, "trial_example.csv", row.names = FALSE)
