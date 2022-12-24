# Postprocessing with Autoregressive Heteroscedastic EMOS (AR-EMOS)

AR-EMOS scripts for the ESSD benchmark. Provide the AR-EMOS output file (see the method's details below).

This code is provided as supplementary material with:

* ...: The EUPPBench postprocessing benchmark dataset v1.0, ...

**Please cite this article if you use (a part of) this code for a publication.**

## Method

AR-EMOS was introduced by [Möller and Groß (2016)](https://doi.org/10.1002/qj.2741) for postprocessing of 24h ahead forecasts of t2m. The method utilizes the autoregressive information in the forecast errors of the ensemble forecasts in order to estimate the parameters of a Gaussian predictive distribution in an EMOS-like fashion. 

Later on the estimation of the predictive standard deviation was further refined to combine the variance component obtained from the longitudinal time series information with the variance component from the cross-sectional empirical ensemble variance. The refined method is described in [Möller and Groß (2020)](https://doi.org/10.1002/qj.3667), along with an extension that makes it possible to apply the method to arbitrary forecast horizons. 



## Implementation Details

The original AR-EMOS method as described in the references above is implemented in the programming language [R](https://www.r-project.org), within the R package [ensAR](https://github.com/JuGross/ensAR). 


### Additional R Packages

- [ncdf4](https://cran.r-project.org/web/packages/ncdf4/index.html): For importing the `*.nc`-files into R.
- [zoo](https://cran.r-project.org/web/packages/zoo/index.html): For imputing individual missing values. 

### Data

First, if you do not have it, get the ESSD benchmark dataset using [the download script](https://github.com/EUPP-benchmark/ESSD-benchmark-datasets). This will fetch the dataset into NetCDF files on your disk.

### R-Script Information

- To construct data frames from the provided `*.nc`-files, the the R-script `benchmarkdata.R` was used. This results i the data frames `benchmark_t2m_train.Rdata`, which denotes the training data and `benchmark_t2m_test.Rdata`, which contains the test/validation data. 
- The R-script can be found in the D-vine postprocessing repository of David Jobst (https://github.com/EUPP-benchmark/ESSD-DVQR).
- For AR-EMOS there are 3 R-scripts, the file `FitArEmos.R` contains the function to fit the model to given training data, the file `PredictArEmos.R` contains the function applying the fitted model object of the fit function to (new test) data. The script `ArEmos_BenchmarkData.R` executes the fit and predict function on the pre-processed ESSD benchmark data in a loop over all forecast horizons and stations. 


### Details of AR-EMOS Implementation

- The postprocessing is performed **locally**, i.e. the `locations` (stations) are postprocessed separately, for each `leadtime`individually.
- Contrary to the original implementation in the ensAR package, where the model parameters were estimated based on a rolling training period throughout the full data, the implementation for the benchmark estimates the model only **once from a static set of training data**. The estimation is performed for each station and each lead time separately. 
- This change was necessary due to the fact that training and test data have different distances between the time points at which data is observed (every 3-4 days in the training data, each day in the test data), so that a rolling through the training data with first prediction date being the first date in the test data would not have been reasonable. Implementing any other form of semi-rolling training period (seasonal, weekly, monthly) also is not straightforward due to the different time resolutions in the two data sets. This is supposed to be investigated in more detail for the second benchmark part. 
- Furthermore, the number of ensemble members in the training and test data differs, requiring further modifications. While in in the original implementation an AR process was fitted to each ensemble member individually, this is not possible when test and training data have different number of members. 
- To compromise, an **AR process was fitted separately to the control member** (as this is present in both data sets) **and one to the mean of the remaining members**. The estimated parameters from the control member were used to obtain the AR-corrected control member, the estimated parameters from the ensemble mean were then applied to each of the members (instead of applying it to the ensemble mean also) in order to generate an ensemble of Ar-corrected forecasts and not only an AR-corrected mean forecast. 
- As AR-EMOS is based on a time series method a time series with equally spaced time points is required. What the time unit is (1 day, 1 hour, 3 days, ...) is not important, it is more important that the time points are as regular as possible, without (too extensive) gaps in between. Thus, missing values, especially a longer consecutive series of missing values can deteriorate the model estimation. In that regard, longer consecutive series of missing values were removed completely from the training data and model estimation performed only on the remaining data. 
- The really long blocks of missing values are mostly concentrated at the beginning of the training data and sometimes at the end. So in fact the training data was reduced by cutting out these blocks. For individual or only few consecutive missing values in the middle of the data a spline-based interpolation was used for imputation - as in the original implementation of AR-EMOS. 

Author: Annette Möller
