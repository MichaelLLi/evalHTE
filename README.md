# evalHTE

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/evalHTE)](https://CRAN.R-project.org/package=evalHTE)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/evalHTE)](https://CRAN.R-project.org/package=evalHTE)
<!-- badges: end -->

## Overview

**evalHTE** provides statistical methods for evaluating heterogeneous treatment effects (HTE) in randomized experiments. The package implements the methodology developed in [Imai and Li (2025)](https://doi.org/10.1080/07350015.2024.2421995) for:

- Estimating **Grouped Average Treatment Effects (GATEs)** with uniform confidence bands
- Identifying **exceptional responders** — individuals who benefit the most (or are harmed) by a treatment
- Conducting **hypothesis tests** for treatment effect heterogeneity and consistency
- Supporting both **sample splitting** and **cross-validation** approaches

## Installation

### From CRAN

```r
install.packages("evalHTE")
```

### Development Version from GitHub

```r
# install.packages("devtools")
devtools::install_github("MichaelLLi/evalHTE")
```

## Key Features

| Feature | Description |
|---------|-------------|
| **GATE Estimation** | Estimate grouped average treatment effects sorted by ML-predicted treatment effects |
| **Uniform Confidence Bands** | Construct valid confidence intervals that account for multiple testing |
| **Exceptional Responders** | Identify subgroups with the largest treatment effects using URATE |
| **Hypothesis Testing** | Test for heterogeneity and consistency of treatment effects across groups |
| **Cross-Validation Support** | Proper inference under cross-fitting to avoid overfitting |
| **Flexible ML Integration** | Works with various ML algorithms via `caret` or user-defined models |

## Quick Start

```r
library(evalHTE)

# Prepare your data with outcome (y), treatment (z), and covariates
data <- data.frame(
  y = outcome_variable,
  z = treatment_indicator,
  x1 = covariate_1,
  x2 = covariate_2
)

# Step 1: Estimate heterogeneous treatment effects
fit <- estimate_hte(
  treatment = "z",
  form = y ~ z * (x1 + x2),
  data = data,
  algorithms = c("causal_forest"),
  n_folds = 5,
  ngates = 5
)

# Step 2: Evaluate the estimated HTE
est <- evaluate_hte(fit)

# Step 3: View summary results
summary(est)

# Step 4: Conduct hypothesis tests
test_results <- test_itr(est)
summary(test_results)

# Step 5: Visualize results
plot(est)                # Plot GATE estimates
plot_CI(est, alpha = 0.05)  # Plot uniform confidence intervals
```

## Main Functions

### Estimation

| Function | Description |
|----------|-------------|
| `estimate_hte()` | Estimate heterogeneous treatment effects using ML algorithms |
| `evaluate_hte()` | Evaluate and summarize HTE estimates |
| `test_itr()` | Conduct hypothesis tests for heterogeneity and consistency |

### GATE Functions

| Function | Description |
|----------|-------------|
| `GATE()` | Estimate Grouped Average Treatment Effects (sample splitting) |
| `GATEcv()` | Estimate GATEs under cross-validation |

### Hypothesis Tests

| Function | Description |
|----------|-------------|
| `het.test()` | Test for heterogeneous treatment effects across groups |
| `hetcv.test()` | Heterogeneity test under cross-validation |
| `consist.test()` | Test for consistency of treatment effect ranking |
| `consistcv.test()` | Consistency test under cross-validation |

### Exceptional Responders

| Function | Description |
|----------|-------------|
| `URATE()` | Estimate the Uplift RATE curve for identifying exceptional responders |

### Visualization

| Function | Description |
|----------|-------------|
| `plot.hte()` | Plot GATE estimates with confidence intervals |
| `plot_CI()` | Plot uniform and pointwise confidence bands |

## Example with Simulated Data

```r
library(evalHTE)
library(grf)

# Simulate data
set.seed(123)
n <- 1000
X <- matrix(rnorm(n * 5), n, 5)
colnames(X) <- paste0("x", 1:5)
tau <- X[, 1] + 0.5 * X[, 2]  # True treatment effect heterogeneity
W <- rbinom(n, 1, 0.5)        # Random treatment assignment
Y <- tau * W + rnorm(n)       # Observed outcome

# Create data frame
df <- data.frame(X, y = Y, z = W)

# Estimate HTE using causal forest
fit <- estimate_hte(
  treatment = "z",
  form = y ~ z * (x1 + x2 + x3 + x4 + x5),
  data = df,
  algorithms = c("causal_forest"),
  n_folds = 5,
  ngates = 5
)

# Evaluate results
est <- evaluate_hte(fit)
summary(est)

# Visualize GATE estimates
plot(est)
```

## Supported Algorithms

The package supports various machine learning algorithms for estimating treatment effects:

- **Causal Forest** (`causal_forest`) via the `grf` package
- **S-learner** and **T-learner** meta-learners
- Any algorithm supported by `caret`
- User-defined custom models

## Citation

If you use this package in your research, please cite:

```bibtex
@article{imai2025statistical,
  title={Statistical Inference for Heterogeneous Treatment Effects Discovered by Generic Machine Learning in Randomized Experiments},
  author={Imai, Kosuke and Li, Michael Lingzhi},
  journal={Journal of Business \& Economic Statistics},
  volume={43},
  number={1},
  pages={256--268},
  year={2025},
  doi={10.1080/07350015.2024.2421995}
}
```

Or in R:

```r
citation("evalHTE")
```

## Authors

- **Michael Lingzhi Li** (maintainer) - Harvard Business School - [mili@hbs.edu](mailto:mili@hbs.edu)
- **Kosuke Imai** - Harvard University - [imai@harvard.edu](mailto:imai@harvard.edu)

Contributor:
- **Jialu Li** - Harvard University

## Related Packages

- [evalITR](https://CRAN.R-project.org/package=evalITR) - Evaluation of Individualized Treatment Rules
- [grf](https://CRAN.R-project.org/package=grf) - Generalized Random Forests
- [caret](https://CRAN.R-project.org/package=caret) - Classification and Regression Training

## License
MIT License

## References

Imai, K. and Li, M. L. (2025). "Statistical Inference for Heterogeneous Treatment Effects Discovered by Generic Machine Learning in Randomized Experiments." *Journal of Business & Economic Statistics*, Vol. 43, No. 1, pp. 256-268. https://doi.org/10.1080/07350015.2024.2421995
