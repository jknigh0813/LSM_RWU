# LSM_RWU
Root water uptake model files

[![DOI](https://zenodo.org/badge/216648866.svg)](https://zenodo.org/badge/latestdoi/216648866)

## General Instructions
* NetCDF files from Chen et al (2008) and Livneh et al (2013) must be downloaded locally.
* Use Preprocess_ForcingData_Sites.m to extract a time series of daily minimum and maximum air temperatures and daily precipitation from Chen et al (2008) (i.e. CPC) and Livneh et al (2013) precipitation datasets. A time series will be generated for each catchment listed in HCDN_Sites.csv. 
* Lumped_VSA_model_Mod.R calls the modified version of JoFlo (discussed below). 
* Use function DDS_Wrapper_LSM_GW.R to execute the Dynamically Dimensioned Search (DDS) routine to calibrate the modified JoFlo model to a USGS gage. Arguments required are the USGS gage ID number and the number of calibration iterations.
* Use Cluster_DDS_Wrapper.R to call DDS_Wrapper_LSM_GW.R from catchments listed in HCDN_Sites.csv
* Cluster_DDS_Wrapper.R is called twice to calibrate Lumped_VSA_model_Mod.R to both Chen et al (2008) and Livneh et al (2013)

## Lumped_VSA_model_Mod.R

* Modified from Archibald et al. (2014) available in Fuka et al (2014)
* includes exponential term for gravity driven baseflow: A1*S^B1
* includes five vertically stacked soil layers
* includes Root Water Uptake (RWU) depth-demand function 

## References
Archibald, J. A., Buchanan, B. P., Fuka, D. R., Georgakakos, C. B., Lyon, S. W., & Walter, M. T. (2014). A simple, regionally parameterized model for predicting nonpoint source areas in the northeastern US. Journal of Hydrology: Regional Studies, 1, 74-91.

Chen, M., Shi, W., Xie, P., Silva, V. B., Kousky, V. E., Wayne Higgins, R., & Janowiak, J. E. (2008). Assessing objective techniques for gauge‐based analyses of global daily precipitation. Journal of Geophysical Research: Atmospheres, 113(D4).

Fuka, D. R., Walter, M. T., Archibald, J. A., Steenhuis, T. S., Easton, Z. M., Fuka, M. D., & KeepSource, T. R. U. E. (2014). Package ‘EcoHydRology’.

Livneh, B., Rosenberg, E. A., Lin, C., Nijssen, B., Mishra, V., Andreadis, K. M., ... & Lettenmaier, D. P. (2013). A long-term hydrologically based dataset of land surface fluxes and states for the conterminous United States: Update and extensions. Journal of Climate, 26(23), 9384-9392.
