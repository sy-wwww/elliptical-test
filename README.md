# Code for Testing Elliptical Models in High Dimensions

For details on the methodology, refer to the paper on arXiv: https://arxiv.org/abs/2408.05514.
The contents of the 'level' and 'power' directories of the reproducibility materials are described below. 

## Level 


*elliptical_test.m* -- This script calculates the elliptical test statistic $T_n$, as well as an estimate for its variance, and its p-value.

*setting.m* -- This script contains the settings used in Section 4.1 to assess the level of the proposed test.

*elliptical_level.m* -- This script calculates the empirical level for the proposed test with a 5% nominal level.

*beta_prime.m* -- This script generates Beta Prime random variables. The input parameters are a positive integer n, a positive integer p, and a positive number $\tau$. The function generates n i.i.d. Beta Prime random variables with parameters $\frac{(1 + p + \tau)  p}{\tau}$ and $\frac{1 + p + 2\tau}{\tau}$.


## Power


*tests.R* -- This script contains the function 'elliptical_test', which calculates the elliptical test statistic $T_n$, as well as an estimate for its variance, and its p-value. In addition, this script contains the function 'normality_test', which calculates the p-value of the normality test.

*elliptical_power.R* -- This script supports the experiments in Section 4.2 of the paper, and computes p-valueS for both the elliptical test and the normality test.

*power_plot.R* -- This script plots the figures in Section 4.2.

*funcs.r* and *funcs2.R* -- These functions are called by tests.R.

### Stock dataset

*stock_data_elliptical.csv* -- The stock market dataset discussed in section 4.3 of the paper is `stock_data_elliptical.csv` in the 'power' directory of the reproducibility materials, which consists of 120 rows and 480 columns. The columns correspond to different stocks, which are labeled by their respective ticker symbols, e.g. AMZN for Amazon. The rows correspond to dates, which are labeled in the month/day/year format. Each numerical entry in the .csv file represents the monthly log return of a given stock (column) computed on a given date (row).


*stock_dataset.R* -- This script supports the experiments in Section 4.3 of the paper, and calculates p-values for the elliptical test and normality test. 


### Breast cancer dataset
This is the VIJVER dataset, available in the R package 'cancerdata'. Documentation is provided in the manual 'cancerdata.pdf' located in the 'power' directory of the reproducibility materials. In brief, the VIJVER dataset is structured such that each row represents a gene (variable), and each column corresponds to a tumor sample.

*cancerdata.pdf* -- This is the manual for the R package 'cancerdata'.

*breast_cancer_dataset.R* -- This script supports the experiments in Section 4.4 of the paper, and calculates p-values for the elliptical test and the normality test.
