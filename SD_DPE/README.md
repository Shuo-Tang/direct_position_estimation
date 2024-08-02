# SD-DPE

In this project, we extend the Direct Position Estimation (DPE) approach to Differential Positioning GNSS (DPGS) by utilizing raw satellite signals from multiple receivers. This method, termed SD-DPE, is designed to retain the high sensitivity of standard DPE while enhancing positioning performance by eliminating common error terms, such as ionospheric and tropospheric errors.

We propose the SD-DPE approach by reconstructing raw signals from distributed base stations. This offers the following advantages:

1. Circumvents noise caused by cross-correlation between signals from different satellites.
2. Facilitates the implementation and popularization of the DGPS DPE algorithm due to the widespread distribution of base stations and easy access to their measurements.

The performance of the SD-DPE approach is influenced by three key factors:

1. Common ionospheric and tropospheric errors.
2. Uncertainty in the measurements from the base station.
3. Signal quality from the rover.

## Scripts

- **`DPE_SD_ARS_varying_CNo.m`**  
  Demonstrates the relationship between positioning precision and the carrier-to-noise ratio (CN0) of the rover signal.  
  <img src="figs/error_vs_CNo.png" alt="Positioning Error vs CN0" width="600"/>

- **`DPE_SD_ARS_varying_iono.m`**  
  Demonstrates how positioning precision is affected by ionospheric error.

- **`DPE_SD_ARS_varying_sigma_range.m`**  
  Demonstrates the impact of measurement uncertainty from the base station on positioning precision.
