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



## Power


*tests.R* -- contains the function 'elliptical_test' which calculates $T_n$, $\sigma_n^2$ and the p-value of the proposed test, as well as the function 'normality_test' which calculates the p-value of the normality test.

*elliptical_power.R* -- compute the empirical power for both the proposed test and the normality test.

*power_plot.R* -- plot the figures in Section 4.2.

### Stock dataset
There are two main ways to get stock history data from Yahoo Finance:

1. **Direct Download**:  
   We can directly download historical data for a specific stock from Yahoo's website. For example, the page [3M Company (MMM) Historical Data](https://finance.yahoo.com/quote/MMM/history/) provides the historical data for 3M Company (MMM).

2. **Using Kaggle's Method**:  
   Kaggle provides a Python script for downloading historical data of multiple stocks in the dataset titled ["S&P 500 Stock Data"](https://www.kaggle.com/datasets/camnugent/sandp500/data?select=getSandP.py). The script downloads data for various stocks and combines them into a single dataset, `all_stocks_5yr.csv`.
   
    We implemented the code provided in the Kaggle dataset with the following modifications:

- Changed the `datetime` parameters to extend the historical data period:
  ```python
  now_time = datetime(2022, 6, 30)
  start_time = datetime(now_time.year - 10, now_time.month, now_time.day - 3)
  ```
*stock_data_elliptical.csv* -- preprocess the historic dataset of stocks to get the dataset we will test -- stock_data_elliptical.csv

*stock_data_elliptical.csv* -- the dataset used in Section 4.3.

*stock_dataset.R* -- the code used in Section 4.3.

*breast_cancer.csv* -- the dataset used in Section 4.4.

*breast_cancer_dataset.R* -- the code used in Section 4.4.

*funcs.r* and *funcs2.R* -- functions used in the normality test.
