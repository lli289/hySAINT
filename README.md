# hySAINT

**Hybrid Genetic and Simulated Annealing Algorithm for High-Dimensional Linear Models with Interaction Effects**

[![R](https://img.shields.io/badge/R-4.0%2B-blue)](https://www.r-project.org/)
[![License: GPL-2](https://img.shields.io/badge/License-GPL%202-green.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

## Overview

`hySAINT` (Hybrid Simulated Annealing and Intelligent Selection) is an R package that provides a stage-wise selection method using genetic and simulated annealing algorithms to efficiently identify main effects and two-way interaction effects in high-dimensional linear regression models. The algorithm can enforce heredity constraints (Strong, Weak, or No heredity) and is particularly useful when the number of predictors is large relative to the sample size.

### Key Features

- **Hybrid Algorithm**: Combines genetic algorithms with simulated annealing for robust variable selection
- **Interaction Detection**: Efficiently identifies two-way interaction effects in high-dimensional settings
- **Heredity Constraints**: Supports Strong, Weak, or No heredity conditions
- **Scalability**: Handles high-dimensional data where p >> n
- **Model Selection**: Uses ABC criterion for model evaluation

## Installation

### Prerequisites

You need **R 4.0+** on your system. Check your R version:

```r
R.version.string
```

### Install from GitHub

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install hySAINT
devtools::install_github("lli289/hySAINT")
```

### Install Dependencies

The package requires the following R packages:
- `utils`
- `Matrix`
- `pracma`
- `stats`
- `SIS`

These will be installed automatically when you install `hySAINT`.

### Install from source

If you have the package source code locally:

```r
# Set working directory to the package parent folder
setwd("path/to/hySAINT")

# Install dependencies
install.packages(c("Matrix", "pracma", "SIS"))

# Install the package
install.packages("hySAINT", repos = NULL, type = "source")
```

## Quick Start

Here's a simple example to get started:

```r
library(hySAINT)

# Simulate data
set.seed(123)
n <- 100
p <- 20
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("X", 1:p)

# True model: y = 1 + X1 + X2 + X3 + X1*X2 + X2*X3 + noise
y <- 1 + X[,1] + X[,2] + X[,3] + X[,1]*X[,2] + X[,2]*X[,3] + rnorm(n, 0, 0.5)

# Estimate sigma (noise standard deviation)
sigma <- mad(residuals(lm(y ~ X)))

# Run hySAINT
result <- hySAINT(
  X = X,
  y = y,
  heredity = "Strong",
  r1 = 5,           # Max number of main effects
  r2 = 3,           # Max number of interactions
  sigma = sigma,
  max.iter = 50     # Reduce for quick testing
)

# View results
print(result)
```

## Example with R Built-in Dataset: `mtcars`

Here's a practical example using the `mtcars` dataset to predict fuel efficiency (mpg):

```r
library(hySAINT)

# Load the mtcars dataset
data(mtcars)

# Prepare the data
# Use multiple predictors to find main effects and interactions
X <- as.matrix(mtcars[, c("cyl", "disp", "hp", "drat", "wt", "qsec", "vs", "am", "gear", "carb")])
y <- mtcars$mpg

# Standardize predictors (recommended for better performance)
X <- scale(X)

# Estimate sigma using MAD estimator
lm_fit <- lm(y ~ X)
sigma <- mad(residuals(lm_fit))

cat("Dataset info:\n")
cat("Sample size (n):", nrow(X), "\n")
cat("Number of predictors (p):", ncol(X), "\n")
cat("Estimated sigma:", round(sigma, 3), "\n\n")

# Run hySAINT to find main effects and interactions
set.seed(42)
result <- hySAINT(
  X = X,
  y = y,
  heredity = "Strong",    # Enforce strong heredity
  r1 = 6,                 # Allow up to 6 main effects
  r2 = 3,                 # Allow up to 3 interactions
  sigma = sigma,
  numElite = 20,          # Number of elite parents
  max.iter = 100,         # Number of iterations
  initial.temp = 1000,    # Initial temperature for SA
  cooling.rate = 0.95,    # Cooling rate
  lambda = 10             # Lambda parameter for ABC
)

# Display results
cat("============================================\n")
cat("hySAINT Results for mtcars Dataset\n")
cat("============================================\n\n")

cat("Selected Variables:\n")
print(result$Final.variable.names)

cat("\nSelected Variable Indices:\n")
print(result$Final.variable.idx)

cat("\nFinal ABC Score:", round(result$Final.model.score, 2), "\n")

cat("\nNumber of Iterations:", result$n.iterations, "\n")

# Plot convergence
if (length(result$All.iter.score) > 1) {
  plot(0:(length(result$All.iter.score)-1), result$All.iter.score,
       type = "l", col = "blue", lwd = 2,
       xlab = "Iteration", ylab = "Best ABC Score",
       main = "hySAINT Convergence on mtcars Data")
  grid()
}

# Fit the final model with selected variables
# Note: Use Extract() to properly build the design matrix with interactions
if (length(result$Final.variable.idx) > 0) {
  X_selected <- Extract(X, varind = result$Final.variable.idx)
  final_model <- lm(y ~ X_selected)
  
  cat("\n============================================\n")
  cat("Final Model Summary\n")
  cat("============================================\n")
  print(summary(final_model))
}
```

### Expected Output

The package will identify important main effects (like weight `wt`, horsepower `hp`, number of cylinders `cyl`) and their interactions that significantly affect fuel efficiency.

## Understanding Parameters

### Required Parameters

- **`X`**: Input data matrix (n × p). Each row is an observation, each column is a predictor.
- **`y`**: Response vector (n × 1).
- **`r1`**: Maximum number of main effects to include.
- **`r2`**: Maximum number of interaction effects to include.
- **`sigma`**: Standard deviation of the noise term. Estimate using `sigma = mad(residuals(lm(y ~ X)))`.

### Optional Parameters

- **`heredity`**: Heredity constraint type. Options:
  - `"Strong"` (default): An interaction X_i * X_j can only be selected if both X_i and X_j are in the model.
  - `"Weak"`: An interaction X_i * X_j can be selected if at least one of X_i or X_j is in the model.
  - `"No"`: No heredity constraints.
  
- **`varind`**: Indices of variables to screen (for high-dimensional data). Leave as `NULL` for automatic screening.

- **`numElite`**: Number of elite parents in the genetic algorithm (default: 40).

- **`max.iter`**: Maximum number of iterations (default: 500).

- **`initial.temp`**: Initial temperature for simulated annealing (default: 1000).

- **`cooling.rate`**: Cooling rate for simulated annealing (default: 0.95).

- **`lambda`**: Penalty parameter for ABC criterion. Must satisfy λ ≥ 5.1/log(2) (default: 10).

## Output

The function returns an object of class `"hySAINT"` with the following components:

- **`Final.variable.names`**: Names of the selected effects (main effects and interactions).
- **`Final.variable.idx`**: Indices of the selected effects.
- **`Final.model.score`**: Final ABC score of the best model.
- **`All.iter.score`**: Vector of best ABC scores from each iteration (useful for convergence analysis).
- **`n.iterations`**: Number of iterations performed.

## Advanced Usage

### High-Dimensional Example (p >> n)

When the number of predictors exceeds the sample size, hySAINT automatically applies screening:

```r
set.seed(456)
n <- 100
p <- 500  # High-dimensional setting

X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("X", 1:p)

# True sparse model
true_main <- c(1, 2, 5, 10, 25)
true_interactions <- matrix(c(1,2, 5,10), ncol=2, byrow=TRUE)

y <- 1 + rowSums(X[, true_main]) + 
     X[,1]*X[,2] + X[,5]*X[,10] + 
     rnorm(n, 0, 1)

sigma <- mad(residuals(lm(y ~ X)))

# Run hySAINT (automatic screening will be applied)
result <- hySAINT(
  X = X,
  y = y,
  heredity = "Strong",
  r1 = 10,
  r2 = 5,
  sigma = sigma,
  max.iter = 200
)

print(result$Final.variable.names)
```

### Tuning the Algorithm

For better performance, consider:

1. **Increasing `max.iter`**: More iterations allow better exploration (default: 500).
2. **Adjusting `numElite`**: More elite parents provide diversity (default: 40).
3. **Tuning `cooling.rate`**: Slower cooling (closer to 1) allows more thorough search (default: 0.95).
4. **Selecting `lambda`**: Higher values impose stronger penalties for model complexity (default: 10).


## Citation

If you use this package in your research, please cite:

```
L. Li, Z. You, and C. Ye (2026). hySAINT: Hybrid Genetic and Simulated Annealing Algorithm for High-Dimensional Linear Models with Interaction Effects. R package version 1.3.0. https://github.com/lli289/hySAINT
```

## Authors

- **Leiyue Li** (Maintainer) - [lli289.git@gmail.com](mailto:lli289.git@gmail.com)
- **Chenglong Ye** - [chenglong.ye@uky.edu](mailto:chenglong.ye@uky.edu)
- **Zhengzhong You** -

## License

This package is licensed under the GPL-2 License. See the `LICENSE` file for details.

## Troubleshooting

### Common Issues

**Problem**: Installation fails with dependency errors.
```r
# Solution: Install dependencies manually
install.packages(c("Matrix", "pracma", "SIS"))
```

**Problem**: Out of memory errors.
```r
# Solution: Reduce max.iter, r1, or r2, or use screening for high-dimensional data
```

**Problem**: Slow performance.
```r
# Solution: Reduce max.iter or numElite for faster (but potentially less accurate) results
result <- hySAINT(X, y, r1=5, r2=2, sigma=sigma, max.iter=50, numElite=20)
```

### Getting Help

- Check the function documentation: `?hySAINT`
- View examples: `example(hySAINT)`
- Report issues: [GitHub Issues](https://github.com/lli289/hySAINT/issues)

## Acknowledgments

This work is based on the paper by Chenglong Ye and Yuhong Yang (2019) on high-dimensional variable selection with interaction effects.
