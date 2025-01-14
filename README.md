---
title: "Code for Testing Elliptical Models in High Dimensions"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Level


*elliptical_test.m* -- calculate $T_n$, $\sigma_n^2$ and the p-value.

*setting.m* -- settings used in Section 4.1 to assess the level of the proposed test.

*elliptical_level.m* -- calculate the empirical level for the proposed test with a 5% nominal level.

*beta_prime.m* -- Use the Beta distribution to generate a Beta Prime distribution. The input parameters are a positive integer \(n\), a positive integer \(p\), and a positive number \(\tau\). The function generates an \(n\)-dimensional independent Beta Prime distribution with parameters \(\frac{(1 + p + \tau)p}{\tau}\) and \(\frac{1 + p + 2\tau}{\tau}\).


## Power


*tests.R* -- contains the function 'elliptical_test' which calculates $T_n$, $\sigma_n^2$ and the p-value of the proposed test, as well as the function 'normality_test' which calculates the p-value of the normality test.

*elliptical_power.R* -- compute the empirical power for both the proposed test and the normality test.

*power_plot.R* -- plot the figures in Section 4.2.

*funcs.r* and *funcs2.R* -- functions used in the normality test.

### Stock dataset

The dataset consists of the monthly log returns of 480 stocks from the S\&P 500 that maintained complete trading histories between July 2012 and June 2022. It is organized such that each row represents a specific month, and each column corresponds to a stock.

*stock_data_elliptical.csv* -- The dataset used in Section 4.3; the dataset is structured such that each row represents a specific month, and each column corresponds to a stock.

*stock_dataset.R* -- The code used in Section 4.3 calculates the power of the elliptical test and the normality test.


### Breast cancer dataset
The dataset is available in the R package `cancerdata`, with a detailed description provided in the documentation, `cancerdata.pdf`. For our analysis, we use the VIJVER dataset. The original data is a \(24481 \times 295\) matrix, where each row represents a gene, and each column corresponds to a tumor sample.

*cancerdata.pdf* -- A comprehensive description of the R package `cancerdata` is provided.

*breast_cancer_dataset.R* -- The code used in Section 4.4 calculates the power of the elliptical test and the normality test.
