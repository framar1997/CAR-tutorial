# Allocation examples for the procedures discussed in the tutorial.
# Run from the repository root: source("examples/01_all_designs.R")
source("R/01_designs.R")

# One sequence of enrolled patients, in their actual enrollment order.
set.seed(2001)
n <- 96L
baseline <- data.frame(
  sex = factor(sample(c("F", "M"), n, replace = TRUE)),
  stage = factor(sample(c("I", "II", "III"), n, replace = TRUE)),
  x1 = rnorm(n), x2 = rnorm(n)
)

# Stratified designs use the categorical factors. Quantitative x1/x2 are
# retained for analysis and for feature-based designs where requested.
p_cat <- car_prepare(
  baseline,
  categorical = c("sex", "stage"),
  continuous = c("x1", "x2"),
  analysis_covariates = c("sex", "stage", "x1", "x2")
)

# Here the feature-based rules use the ORIGINAL continuous measurements
# (plus the intercept) rather than categorical stratum indicators.
p_cont <- car_prepare(
  baseline,
  categorical = c("sex", "stage"),
  continuous = c("x1", "x2"),
  feature_cat = character(),
  feature_cont = c("x1", "x2"),
  analysis_covariates = c("sex", "stage", "x1", "x2")
)

# Explicit car_design() examples, including the principal tuning parameters.
plans <- list(
  CR       = car_design(p_cat, method = "CR"),
  PBD      = car_design(p_cat, method = "PBD", block_size = 6),
  RBD      = car_design(p_cat, method = "RBD", random_block_sizes = c(4, 6, 8)),
  BSD      = car_design(p_cat, method = "BSD", mti = 3),
  "CA-BCD" = car_design(p_cat, method = "CA-BCD", a = 1),
  HH       = car_design(p_cat, method = "HH", rho = 0.85,
                        hh_weights = c(0.25, 0.25, 0.25, 0.25)),
  PS       = car_design(p_cat, method = "PS", ps_rho = 0.75),
  "DA-BCD" = car_design(p_cont, method = "DA-BCD"),
  ECADE    = car_design(p_cont, method = "ECADE", rho = 0.85,
                        ecade_link = "Efron"),
  KER      = car_design(p_cont, method = "KER", ker_rho = 0.80),
  COV      = car_design(p_cont, method = "COV", cov_rho = 0.85,
                        cov_weights = c(1, 1, 1))
)

# Draw ONE treatment sequence for each allocation procedure.
allocations <- lapply(seq_along(plans), function(j) {
  car_allocate(plans[[j]], seed = 2100L + j)
})
names(allocations) <- names(plans)

# A concise description of the realized sequences, not Monte Carlo averages.
allocation_summary <- do.call(rbind, lapply(names(plans), function(method) {
  a <- allocations[[method]]
  z <- plans[[method]]$prepared$design_strata
  within_stratum <- vapply(split(a$treatment, z), function(d) {
    sum(2L * d - 1L)
  }, numeric(1))
  data.frame(
    method = method,
    treatment_n = sum(a$treatment == 1L),
    control_n = sum(a$treatment == 0L),
    maximum_absolute_stratum_imbalance = max(abs(within_stratum)),
    row.names = NULL
  )
}))

print(allocation_summary, row.names = FALSE)

# Access individual patient-level allocations and probabilities, e.g.:
print(head(allocations$ECADE, 6))
print(head(allocations$PBD, 6))

# For each design, the assigned treatment (0 or 1) and probability_T are
# indexed by recruitment order. Do not interpret one allocation as an estimate
# of a procedure's operating characteristics across repeated trials.
