# A SINGLE continuous-covariate trial: original X balanced with ECADE,
# tertile strata used only for t and W; regression uses original X.
# Run from the repository root (base R only).
source("R/01_designs.R")
source("R/02_inference.R")
source("R/03_trial.R")

n <- 120L
set.seed(2026)
d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
prespecified_tertiles <- qnorm(c(1/3, 2/3))
d$bin1 <- car_discretize(d$x1, prespecified_tertiles)
d$bin2 <- car_discretize(d$x2, prespecified_tertiles)

prepared <- car_prepare(d,
  categorical = c("bin1", "bin2"), continuous = c("x1", "x2"),
  feature_cat = character(), feature_cont = c("x1", "x2"),
  analysis_covariates = c("x1", "x2"))

plan <- car_design(prepared, method = "ECADE")
a <- car_allocate(plan, seed = 555)

# Demo only: after a real trial use OBSERVED outcomes rather than this formula.
set.seed(556)
y <- .4 * a$treatment + sqrt(1/2) * (d$x1+d$x2) + rnorm(n)

out <- car_analyze_trial(y, a$treatment, plan, L = 299, B = 59, seed = 99)
print(out)

# A single dataset cannot be used to estimate type-I error or power.
# Those are repeated-trial operating characteristics, not individual-trial outputs.
