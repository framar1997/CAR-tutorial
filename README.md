# Covariate-Adaptive Randomization: Design and Inference in R

**Companion code for the tutorial _Theoretical and Practical Guidance for Covariate-Adaptive Randomization_**  
Rosamarie Frieri, Francesco Mariani, Marco Novelli and Sara Cecconi

## Overview

Covariate-adaptive randomization (CAR) incorporates baseline prognostic information into treatment allocation to limit clinically relevant imbalances. The choice of an allocation rule involves a trade-off between covariate balance and allocation unpredictability; importantly, the properties of the subsequent statistical analysis also depend on how baseline covariates were used during randomization.

This repository accompanies the tutorial _Theoretical and Practical Guidance for Covariate-Adaptive Randomization_. It provides reproducible examples of sequential treatment allocation and outcome analysis in two-arm randomized trials. The code illustrates stratified randomization, minimization, and procedures that can balance continuous covariates directly, together with population-based and randomization-based inferential methods.

The central principle is to **prespecify the allocation procedure and a compatible analysis strategy**. In particular, the information balanced during randomization need not be equivalent to the information retained after categorizing continuous covariates, and fully stratified analyses may be unsuitable when strata are numerous and sparsely populated.

## Methods

### Randomization procedures

All procedures target equal allocation to treatment (`1`) and control (`0`). Participants are processed in the order of the rows supplied to the functions.

| Abbreviation | Procedure | Covariate information used for allocation | Main argument(s) |
|:--|:--|:--|:--|
| `CR` | Complete randomization | No covariate restrictions | `method = "CR"` |
| `PBD` | Stratified permuted block design | Joint categorical strata | `block_size = 6` |
| `RBD` | Stratified randomized block design | Joint categorical strata | `random_block_sizes = c(4, 6, 8)` |
| `BSD` | Big Stick Design | Joint categorical strata | `mti = 3` |
| `CA-BCD` | Covariate-adaptive biased coin design | Joint categorical strata | `a = 1` |
| `HH` | Hu and Hu randomization | Overall, marginal and stratum imbalance | `rho = 0.85`, `hh_weights` |
| `PS` | Pocock--Simon minimization | Marginal categorical imbalance | `ps_rho = 0.75` |
| `DA-BCD` | Atkinson's DA-optimum biased coin design | Prespecified design features | `method = "DA-BCD"` |
| `ECADE` | Efficient covariate-adaptive design | Prespecified design features | `rho = 0.85`, `ecade_link` |
| `KER` | Kernel-based covariate-adaptive randomization | Categorical and/or continuous features | `ker_rho = 0.80` |
| `COV` | First- and second-moment feature-map balancing | Prespecified features and cross-products | `cov_rho = 0.85`, `cov_weights` |

Here *design features* refer to the categorical indicators, quantitative covariates or other prespecified transformations supplied by `car_prepare()`. The current interface exposes the ECADE allocation links `"Efron"`, `"Norm1"`, and `"Norm2"`; the weighted ECADE specification discussed in the tutorial is not separately exposed as a `car_design()` method.

For a comparison of the balancing properties and allocation predictability of these procedures, including the distinction between PBD and BSD and the behavior of stratified designs with many sparse strata, see the accompanying tutorial.

### Inference after randomization

The analysis functions consider one completed trial, using its **observed outcomes and recorded treatment assignments**.

| Approach | Test or statistic | Function/output name |
|:--|:--|:--|
| Population-based | Stratified mean-based test | `t` |
| Population-based | Stratified Wilcoxon/van Elteren test | `W` |
| Population-based | Linear regression adjusted for analysis covariates | `treg` |
| Population-based | Regression with HC0 standard errors | `treg_HC0` |
| Population-based | Stratum-adjusted test | `tadj` |
| Population-based | Bootstrap-adjusted mean-based test | `tB` |
| Randomization-based | Weighted stratum-specific mean differences | `D` |
| Randomization-based | Stratified centered-rank statistic | `R` |
| Randomization-based | Difference in means of fixed regression residuals | `Dreg` |
| Randomization-based | Rank statistic on fixed regression residuals | `Rreg` |

The analysis also offers optional exploratory diagnostics (`FWL` for re-randomization and `FE_HC0_diagnostic` for population-based inference). These are not recommended as interchangeable substitutes for the primary procedures, particularly in sparse strata.

For randomization-based tests, the code reproduces the **prespecified sequential allocation mechanism**, conditional on the enrolled participants, their baseline covariates, and recruitment order. With `L - 1` new allocations, the observed assignment is included in the Monte Carlo reference distribution, and the two-sided p-value is

$$
\widehat p = \frac{1 + \sum_{m=1}^{L-1}
\mathbf{1}\{|T_m|\geq |T_{\mathrm{obs}}|\}}{L}.
$$

This calculation has a finite-sample justification under the **sharp null of no individual treatment effect**, when the allocation mechanism is correctly reconstructed and the statistic is defined throughout the reference distribution. It does not automatically yield an exact test of a weak null of zero average treatment effect under heterogeneous effects.

The bootstrap test `tB` estimates its standard error using `B` resamples of patients followed by new allocations under the specified design. Bootstrap calibration is not guaranteed for every covariate structure; the tutorial documents anti-conservative behavior for some designs with many sparse strata.

## Getting started

Download the repository and run R from its root directory:

```r
source("R/01_designs.R")
source("R/02_inference.R")
source("R/03_trial.R")
```

The workflow has four steps: prepare the baseline data, specify an allocation design, obtain the treatment assignments, and—**after outcomes have been observed**—analyze the completed trial.

```r
# Baseline records, in enrollment order
baseline <- car_example_data(n = 180, seed = 100)

prepared <- car_prepare(
  baseline,
  categorical = c("sex", "stage"),
  continuous = "age_z"
)

# Allocation rule and one treatment sequence
plan <- car_design(prepared, method = "PBD", block_size = 6)
allocation <- car_allocate(plan, seed = 101)

# Illustrative outcomes only; replace with observed trial outcomes.
set.seed(102)
y <- 0.4 * allocation$treatment + 0.5 * baseline$age_z +
  0.3 * (baseline$stage == "high") + rnorm(nrow(baseline))

results <- car_analyze_trial(
  y = y,
  treatment = allocation$treatment,
  design = plan,
  L = 199, B = 49, seed = 103
)
print(results)
```

`L = 199` denotes **199 additional allocations plus the observed allocation**, not 199 simulated clinical trials. `B = 49` denotes 49 bootstrap samples used to estimate the bootstrap variance. The small values in the introductory example are for illustration; more accurate Monte Carlo approximation requires larger values. For a completed trial, supply the **actual recorded allocation**, and use the rule and settings that generated it—do not generate a replacement allocation after observing the responses.

### Selecting covariates for design and analysis

`car_prepare()` distinguishes three objects that need not coincide:

- **`categorical`** identifies variables defining categorical strata for the corresponding allocation procedures and default stratified comparisons.
- **`feature_cat` and `feature_cont`** identify the features used by model- and feature-based allocation rules, including ECADE.
- **`analysis_covariates`** defines the regression adjustment. `analysis_strata` can be supplied separately when the inferential strata differ from the design strata.

For example, continuous covariates may be included in ECADE's allocation criterion at their original scale, while prespecified tertiles define the strata for `t` and `W`. Regression-based methods may then adjust for the original continuous variables. This distinction is central to the tutorial's discussion of inference following direct balancing of continuous covariates.

Covariate transformations, cutpoints and design parameters must be determined prospectively. In particular, transformations for a sequential allocation should not use values from patients who have not yet entered the study.

## Worked examples

Each script can be run independently from the repository root using `source()` or `Rscript`. The scripts illustrate **individual trial workflows**; they do not estimate repeated-trial type I error or power.

| File | Illustration |
|:--|:--|
| [`examples/01_all_designs.R`](examples/01_all_designs.R) | Explicit `car_design()` calls for all 11 allocation procedures and one allocation per design. |
| [`examples/02_inference_PBD.R`](examples/02_inference_PBD.R) | Stratified PBD allocation, population-based tests, bootstrap and design-based re-randomization. |
| [`examples/03_inference_PS.R`](examples/03_inference_PS.R) | Pocock--Simon minimization with marginal balancing, followed by inference using the original minimization rule. |
| [`examples/04_inference_ECADE_continuous.R`](examples/04_inference_ECADE_continuous.R) | ECADE balancing continuous covariates directly; discretized strata for `t`/`W`, original covariates for regression and adjusted re-randomization. |
| [`examples/example_trial.R`](examples/example_trial.R) | Compact end-to-end trial example. |
| [`examples/example_continuous.R`](examples/example_continuous.R) | Compact continuous-covariate example. |

For instance:

```r
source("examples/01_all_designs.R")
source("examples/02_inference_PBD.R")
source("examples/03_inference_PS.R")
source("examples/04_inference_ECADE_continuous.R")
```

The first script summarizes the realized allocations; the remaining three print population-based and re-randomization results separately. Because the allocation procedures are sequential, **row order matters** in both the original trial and subsequent re-randomization.

## Functions and output

| Function | Purpose |
|:--|:--|
| `car_discretize()` | Apply prospectively chosen cutpoints to quantitative covariates. |
| `car_prepare()` | Encode stratification factors, balancing features and analysis covariates. |
| `car_design()` | Specify one allocation rule and its tuning parameters. |
| `car_allocate()` | Generate a treatment sequence and record conditional treatment probabilities. |
| `car_population()` | Compute the available population-based analyses. |
| `car_bootstrap()` | Construct a bootstrap distribution for the adjusted `tB` test. |
| `car_rerandomization()` | Approximate design-based randomization distributions using the specified allocation procedure. |
| `car_analyze_trial()` | Run population-based and randomization-based analyses for a completed trial. |

Some stratified statistics require both treatment arms—and, for variance estimation, more than one patient per arm—in each occupied stratum. When a statistic cannot be computed, the corresponding result is reported as `NA`. By default, a re-randomization test is also reported as undefined if its statistic cannot be computed for a reference allocation; this avoids silently changing the reference distribution.

The functions assume complete numeric outcomes, an observed binary assignment, two treatment arms and equal target allocation. The examples do not address missing data, clustered allocation, unequal treatment proportions or time-to-event outcomes. They illustrate the statistical methods discussed in the tutorial and are not a substitute for an independently validated trial randomization or trial-management system.

## Reproducibility

All examples set seeds for illustration. To run the repository's basic checks from its root directory:

```bash
Rscript tests/smoke_tests.R
```

The checks cover allocation reproducibility, all implemented procedures, selected design restrictions and the inference interfaces. Before using the functions in a new study, users should confirm that their design settings, covariate processing, statistical estimand and analysis assumptions match those for which their chosen procedures are justified.

## Reference

Frieri, R., Mariani, F., Novelli, M., and Cecconi, S. **Theoretical and Practical Guidance for Covariate-Adaptive Randomization**. Tutorial manuscript. [Publication information to be added upon acceptance.]
