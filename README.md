# pRRophetic

**Predict clinical chemotherapeutic response from before-treatment tumor gene expression levels**

[![R Version](https://img.shields.io/badge/R-4.0%2B-blue.svg)](https://cran.r-project.org/)
[![License](https://img.shields.io/badge/License-GPL--2-green.svg)](https://www.gnu.org/licenses/gpl-2.0.html)

## Overview

pRRophetic is an R package that predicts clinical chemotherapeutic response from before-treatment tumor gene expression levels. The package uses cell lines from the Cancer Genome Project (CGP) as a training set to build statistical models that can predict drug sensitivity in clinical samples.

**⚠️ Important: R 4.0+ Compatibility**  
This repository contains fixes for R 4.0+ compatibility issues that were present in the original package.

## Key Features

- **Drug Sensitivity Prediction**: Predict IC50 values for over 130 anti-cancer drugs
- **Multiple Tissue Types**: Train models on specific tissue types or all solid tumors
- **Cross-Validation**: Built-in cross-validation for model assessment
- **Batch Effect Correction**: Multiple methods including ComBat, quantile normalization
- **Flexible Models**: Both linear (ridge regression) and logistic regression models
- **Clinical Validation**: Validated on multiple clinical datasets

## Installation

### From GitHub (Recommended - R 4.0+ Compatible)

```r
# Install devtools if you haven't already
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

# Install the R 4.0+ compatible version
devtools::install_github("yourusername/pRRophetic")
```

### Dependencies

```r
# Install required Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("sva", "genefilter", "preprocessCore"))

# Install CRAN packages
install.packages(c("car", "ridge"))
```

## Quick Start

### Basic Drug Sensitivity Prediction

```r
library(pRRophetic)

# Load example data
data("bortezomibData")

# Predict bortezomib sensitivity
predicted_sensitivity <- pRRopheticPredict(
    testMatrix = exprDataBortezomib,
    drug = "Bortezomib",
    tissueType = "blood",
    selection = 1
)

# View predictions
head(predicted_sensitivity)
```

### Cross-Validation Assessment

```r
# Perform 5-fold cross-validation
cv_results <- pRRopheticCV(
    drug = "Doxorubicin",
    cvFold = 5,
    testExprData = exprDataBortezomib
)

# Summary of prediction accuracy
summary(cv_results)
plot(cv_results)
```

### Available Drugs

The package supports prediction for 130+ drugs including:

- **Chemotherapy**: Doxorubicin, Cisplatin, Paclitaxel, Etoposide
- **Targeted Therapy**: Erlotinib, Imatinib, Lapatinib, Sunitinib
- **Novel Agents**: Bortezomib, Temsirolimus, AZD6244

```r
# Check drug distribution normality
pRRopheticQQplot("Doxorubicin")

# Get available tissue types
# "all", "allSolidTumors", "blood", "breast", "CNS", 
# "GI tract", "lung", "skin", "upper aerodigestive"
```

## Advanced Usage

### Tissue-Specific Models

```r
# Train on blood cancers only
blood_prediction <- pRRopheticPredict(
    testMatrix = expression_data,
    drug = "Bortezomib",
    tissueType = "blood"
)

# Train on solid tumors only
solid_prediction <- pRRopheticPredict(
    testMatrix = expression_data,
    drug = "Paclitaxel",
    tissueType = "allSolidTumors"
)
```

### Custom Training Data

```r
# Use your own training data
custom_prediction <- calcPhenotype(
    trainingExprData = training_expression,
    trainingPtype = training_phenotype,
    testExprData = test_expression,
    batchCorrect = "eb",
    powerTransformPhenotype = TRUE
)
```

### Logistic Regression Models

For drugs with highly skewed distributions:

```r
# Use logistic model instead of linear
logistic_prediction <- pRRopheticLogisticPredict(
    testMatrix = expression_data,
    drug = "Erlotinib",
    tissueType = "allSolidTumors"
)
```

## Data Requirements

### Expression Data Format

Your expression matrix should be:
- **Rows**: Genes (with official gene symbols as rownames)
- **Columns**: Samples (with sample IDs as colnames)
- **Values**: Normalized expression values (microarray or RNA-seq)

```r
# Example format
head(expression_matrix[1:5, 1:3])
#         Sample1 Sample2 Sample3
# EGFR       8.2     7.8     9.1
# TP53       6.5     6.2     5.9
# MYC        7.1     8.3     7.8
```

## R 4.0+ Compatibility Fixes

This version includes important fixes for R 4.0+ compatibility:

### Issues Fixed
- **class() function changes**: Fixed "condition has length > 1" errors
- **Matrix class inheritance**: Updated type checking for multiple class inheritance
- **Function parameter validation**: Improved object type detection

### Functions Updated
- `calcPhenotype()`
- `summarizeGenesByMean()`
- All internal type checking functions

### Before (causing errors in R 4.0+):
```r
if (class(exprMatUnique) == "numeric") {
    # This fails in R 4.0+
}
```

### After (R 4.0+ compatible):
```r
if (is.vector(exprMatUnique) || (!is.matrix(exprMatUnique) && !is.data.frame(exprMatUnique))) {
    # This works in all R versions
}
```

## Performance and Validation

### Benchmark Results
- **Correlation with clinical response**: r = 0.4-0.6 for most drugs
- **Cross-validation R²**: 0.2-0.4 depending on drug and tissue type
- **Computational time**: ~2-5 minutes per drug prediction

### Clinical Validation
The package has been validated on multiple clinical datasets:
- Multiple myeloma patients (bortezomib)
- Breast cancer patients (various chemotherapies)  
- Ovarian cancer patients (platinum compounds)

## Citation

If you use pRRophetic in your research, please cite:

```
Geeleher P, Cox NJ, Huang RS. Clinical drug response can be predicted using baseline gene expression levels and in vitro drug sensitivity in cell lines. Genome Biology. 2014;15(3):R47.
```

## Examples and Tutorials

### Example 1: Predicting Doxorubicin Sensitivity

```r
library(pRRophetic)

# Load your expression data
# expression_data <- read.csv("your_expression_data.csv", row.names = 1)

# Predict doxorubicin sensitivity
dox_sensitivity <- pRRopheticPredict(
    testMatrix = expression_data,
    drug = "Doxorubicin",
    tissueType = "allSolidTumors",
    selection = 1,
    printOutput = TRUE
)

# Lower values indicate higher sensitivity
summary(dox_sensitivity)
```

### Example 2: Comparing Multiple Drugs

```r
drugs <- c("Doxorubicin", "Cisplatin", "Paclitaxel", "Etoposide")
predictions <- list()

for (drug in drugs) {
    predictions[[drug]] <- pRRopheticPredict(
        testMatrix = expression_data,
        drug = drug,
        selection = 1,
        printOutput = FALSE
    )
}

# Create results matrix
results_matrix <- do.call(cbind, predictions)
colnames(results_matrix) <- drugs
```

## Troubleshooting

### Common Issues

**1. Gene Symbol Mismatches**
```r
# Ensure your gene symbols match HGNC official symbols
# Use org.Hs.eg.db for symbol conversion if needed
```

**2. Sample Size Warnings**
```r
# Minimum 10 samples recommended
# Increase minNumSamples parameter if needed
pRRopheticPredict(..., minNumSamples = 5)
```

**3. Memory Issues with Large Datasets**
```r
# For large datasets, consider subsetting or using more memory
# Increase R memory limit if necessary
```

### Getting Help

- **Issues**: Report bugs and issues on GitHub
- **Documentation**: Use `?pRRopheticPredict` for detailed help
- **Examples**: See the package vignette for comprehensive examples

## Contributing

We welcome contributions! Please:

1. Fork the repository
2. Create a feature branch
3. Make your changes (ensure R 4.0+ compatibility)
4. Add tests if applicable
5. Submit a pull request

### Development Notes

When contributing, please ensure:
- All functions work with R 4.0+
- Use safe type checking (avoid `class(x) == "type"`)
- Include appropriate error handling
- Update documentation as needed

## License

This package is licensed under GPL-2. See the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Original package developed by Paul Geeleher
- R 4.0+ compatibility fixes by the Abu Saadat

---

**Note**: This repository contains the R 4.0+ compatible version of pRRophetic with important bug fixes. The original package may not work correctly with newer R versions.
