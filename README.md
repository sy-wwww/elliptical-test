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

*funcs.r* and *funcs2.R* -- functions used in the normality test.

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
- Adjusted the s_and_p list:
  ```python
  s_and_p = ['MMM', 'ABT', 'ABBV', 'ABMD', 'ACN', 'ATVI', 'ADBE', 'AMD', 'AAP', 'AES', 'AFL', 'A', 'APD', 'AKAM', 'ALK', 'ALB', 
		   'ARE', 'ALGN', 'ALLE', 'LNT', 'ALL', 'GOOGL', 'GOOG', 'MO', 'AMZN', 'AMCR', 'AEE', 'AAL', 'AEP', 'AXP', 'AIG', 'AMT',
		   'AWK', 'AMP', 'ABC', 'AME', 'AMGN', 'APH', 'ADI', 'ANSS', 'ANTM', 'AON', 'AOS', 'APA', 'AAPL', 'AMAT', 'APTV', 'ADM',
		   'ANET', 'AJG', 'AIZ', 'T', 'ATO', 'ADSK', 'ADP', 'AZO', 'AVB', 'AVY', 'BKR', 'BALL', 'BLL', 'BAC', 'BBWI', 'BAX', 'BDX', 'BRK.B',
		   'BBY', 'BIO', 'TECH', 'BIIB', 'BLK', 'BK', 'BA', 'BKNG', 'BWA', 'BXP', 'BSX', 'BMY', 'AVGO', 'BR', 'BRO', 'BF.B', 'CHRW',
		   'CDNS', 'CZR', 'CPB', 'COF', 'CAH', 'KMX', 'CCL', 'CARR', 'CTLT', 'CAT', 'CBOE', 'CBRE', 'CDW', 'CE', 'CNC', 'CNP', 'CDAY',
		   'CERN', 'CF', 'CRL', 'SCHW', 'CHTR', 'CVX', 'CMG', 'CB', 'CHD', 'CI', 'CINF', 'CTAS', 'CSCO', 'C', 'CFG', 'CTXS', 'CLX',
		   'CME', 'CMS', 'KO', 'CTSH', 'CL', 'CMCSA', 'CMA', 'CAG', 'COP', 'ED', 'STZ', 'CEG', 'COO', 'CPRT', 'GLW', 'CTVA', 'COST', 'CTRA',
		   'CCI', 'CSX', 'CMI', 'CVS', 'DHI', 'DHR', 'DRI', 'DVA', 'DE', 'DAL', 'XRAY', 'DVN', 'DXCM', 'FANG', 'DLR', 'DFS', 'DISCA',
		   'DISCK', 'DISH', 'DIS', 'DG', 'DLTR', 'D', 'DPZ', 'DOV', 'DOW', 'DTE', 'DUK', 'DRE', 'DD', 'DXC', 'EMN', 'ETN', 'EBAY', 'ECL',
		   'EIX', 'EW', 'EA', 'EMR', 'ENPH', 'ETR', 'EOG', 'EPAM', 'EFX', 'EQIX', 'EQR', 'ESS', 'EL', 'ETSY', 'EVRG', 'ES', 'RE', 'EXC', 'EXPE',
		   'EXPD', 'EXR', 'XOM', 'FFIV', 'FDS', 'FB', 'FAST', 'FRT', 'FDX', 'FIS', 'FITB', 'FE', 'FRC', 'FISV', 'FLT', 'FMC', 'F', 'FTNT', 'FTV',
		   'FBHS', 'FOXA', 'FOX', 'BEN', 'FCX', 'GPS', 'GRMN', 'IT', 'GNRC', 'GD', 'GE', 'GIS', 'GM', 'GPC', 'GILD', 'GL', 'GPN', 'GS',
		   'GWW', 'HAL', 'HBI', 'HIG', 'HAS', 'HCA', 'PEAK', 'HSIC', 'HSY', 'HES', 'HPE', 'HLT', 'HOLX', 'HD', 'HON', 'HRL', 'HST', 'HWM',
		   'HPQ', 'HUM', 'HBAN', 'HII', 'IEX', 'IDXX', 'INFO', 'ITW', 'ILMN', 'INCY', 'IR', 'INTC', 'ICE', 'IBM', 'IP', 'IPG', 'IFF', 'INTU',
		   'ISRG', 'IVZ', 'IPGP', 'IQV', 'IRM', 'JKHY', 'J', 'JBHT', 'SJM', 'JNJ', 'JCI', 'JPM', 'JNPR', 'KSU', 'K', 'KEY', 'KEYS', 'KMB',
		   'KIM', 'KMI', 'KLAC', 'KHC', 'KR', 'LHX', 'LH', 'LRCX', 'LW', 'LVS', 'LEG', 'LDOS', 'LEN', 'LLY', 'LNC', 'LIN', 'LYV', 'LKQ', 
		   'LMT', 'L', 'LOW', 'LUMN', 'LYB', 'MTB', 'MRO', 'MPC', 'MKTX', 'MAR', 'MMC', 'MLM', 'MAS', 'MA', 'MTCH', 'MKC', 'MCD', 'MCK', 'MDT',
		   'MRK', 'MET', 'MTD', 'MGM', 'MCHP', 'MU', 'MSFT', 'MAA', 'MRNA', 'MHK', 'MOH', 'TAP', 'MDLZ', 'MPWR', 'MNST', 'MCO', 'MS', 'MOS', 'MSI', 'MSCI',
		   'NDAQ', 'NTAP', 'NFLX', 'NWL', 'NEM', 'NWSA', 'NWS', 'NEE', 'NLSN', 'NKE', 'NI', 'NDSN', 'NSC', 'NTRS', 'NOC', 'NLOK', 'NCLH', 'NRG', 'NUE',
		   'NVDA', 'NVR', 'NXPI', 'ORLY', 'OXY', 'ODFL', 'OMC', 'OKE', 'ORCL', 'OGN', 'OTIS', 'PCAR', 'PKG', 'PH', 'PAYX', 'PAYC', 'PYPL', 
		   'PENN', 'PNR', 'PBCT', 'PEP', 'PKI', 'PFE', 'PM', 'PSX', 'PNW', 'PXD', 'PNC', 'POOL', 'PPG', 'PPL', 'PFG', 'PG', 'PGR', 'PLD', 
		   'PRU', 'PTC', 'PEG', 'PSA', 'PHM', 'PVH', 'QRVO', 'PWR', 'QCOM', 'DGX', 'RL', 'RJF', 'RTX', 'O', 'REG', 'REGN', 'RF', 'RSG', 
		   'RMD', 'RHI', 'ROK', 'ROL', 'ROP', 'ROST', 'RCL', 'SPGI', 'CRM', 'SBAC', 'SLB', 'STX', 'SEE', 'SRE', 'NOW', 'SHW', 'SBNY', 'SPG', 'SWKS',
		   'SNA', 'SEDG', 'SO', 'LUV', 'SWK', 'SBUX', 'STT', 'STE', 'SYK', 'SIVB', 'SYF', 'SNPS', 'SYY', 'TMUS', 'TROW', 'TTWO', 'TPR', 'TGT', 'TEL',
		   'TDY', 'TFX', 'TER', 'TSLA', 'TXN', 'TXT', 'TMO', 'TJX', 'TSCO', 'TT', 'TDG', 'TRV', 'TRMB', 'TFC', 'TWTR', 'TYL', 'TSN', 'UDR',
		   'ULTA', 'USB', 'UAA', 'UA', 'UNP', 'UAL', 'UNH', 'UPS', 'URI', 'UHS', 'VLO', 'VTR', 'VRSN', 'VRSK', 'VZ', 'VRTX', 'VFC', 'VIAC', 
		   'VTRS', 'VICI', 'V', 'VNO', 'VMC', 'WRB', 'WAB', 'WMT', 'WBD', 'WBA', 'DIS', 'WM', 'WAT', 'WEC', 'WFC', 'WELL', 'WST', 'WDC', 'WU', 'WRK', 'WY', 
		   'WHR', 'WMB', 'WTW', 'WLTW', 'WYNN', 'XEL', 'XLNX', 'XYL', 'YUM', 'ZBRA', 'ZBH', 'ZION', 'ZTS', 'CPT', 'VNT', 'MXIM', 'WCG', 'HFC', 'NKTR',
		   'DWDP', 'Q', 'BHF', 'EVHC', 'ARNC', 'COTY', 'AYI', 'FL', 'CXO', 'CSRA', 'CMCSK', 'SIG', 'CPGX', 'BXLT', 'SLG', 'ENDP', 'LVLT',
		   'MNK', 'AMG', 'XEC', 'NAVI', 'GMCR', 'GGP', 'KORS', 'RIG', 'MAC', 'DLPH', 'PETM', 'KRFT', 'ADT', 'ESV', 'ALXN', 'FOSL', 'WPX',
		   'TRIP', 'PRGO', 'GAS', 'CBE', 'ANR', 'COV', 'JOY', 'CVC', 'NFX', 'TYC', 'QEP', 'HP', 'MJN', 'CLF', 'SAI', 'PCLN', 'ARG', 'FTI'
		   'OI', 'HRS', 'AKS', 'COG', 'MEE', 'LO', 'RRC', 'GME', 'JEC', 'TDC', 'TSO', 'LUK', 'KFT', 'ESRX', 'SBL', 'AYE', 'ABK', 'JDSU',
		   'CVG', 'YHOO']
  ```
  The historical data for each stock is stored in the folder stock_data, and the combined dataset is saved as 'all_stocks_5yr.csv'.

The 'all_stocks_5yr.csv' dataset contains historical stock data with the following columns:
- **Date**: The trading date for the stock data.
- **High**: The highest price of the stock on that date.
- **Low**: The lowest price of the stock on that date.
- **Open**: The open price of the stock on that date.
- **Close**: The close price of the stock on that date adjusted for splits.
- **Volume**: The total number of shares traded on that date.
- **Adj Close**: The adjusted close price adjusted for splits and dividend and/or capital gain distributions.
- **Name**: The ticker symbol of the stock (e.g., `AAPL` for Apple, `MMM` for 3M Company).

  
*stock_preprocessing.R* -- Preprocess the `all_stocks_5yr.csv` dataset to generate the testing dataset, `stock_data_elliptical.csv`.

*stock_data_elliptical.csv* -- The dataset used in Section 4.3; the dataset is structured such that each row represents a specific datetime, and each column corresponds to a stock.

*stock_dataset.R* -- The code used in Section 4.3 calculates the power of the elliptical test and the normality test.


### Breast cancer dataset
The dataset is available in the R package `cancerdata`, with a detailed description provided in the documentation, `cancerdata.pdf`. For our analysis, we use the VIJVER dataset. The original data is a \(24481 \times 295\) matrix, where each row represents a gene, and each column corresponds to a tumor sample.

*cancerdata.pdf* -- A comprehensive description of the R package `cancerdata` is provided.

*breast_cancer_preprocessing.R* -- Load and preprocess the VIJVER dataset to generate the testing dataset, `elliptical_breast_cancer.csv`.

*breast_cancer.csv* -- The dataset used in Section 4.4; the dataset is structured such that each row represents a specific tumor sample, and each column corresponds to a gene.

*breast_cancer_dataset.R* -- The code used in Section 4.4 calculates the power of the elliptical test and the normality test.

