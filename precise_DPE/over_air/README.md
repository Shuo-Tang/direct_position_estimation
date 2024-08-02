# Precise DPE -- Over-the-air Signal

This section contains scripts related to the implementation of Direct Position Estimation (DPE) and Precise DPE (PDPE) using over-the-air GNSS signals.

## Scripts

- **`init_1ms.m`**  
  Implementation for GNSS-SDR (traditional 2SP approach).

- **`dpe_limeSDR_4Dxyzdt_coherent_sol_01_02.m`**  
  Implementation for Direct Position Estimation (DPE).

- **`pdpe_limeSDR_3Dxyz_20ms_ars_21_25.m`**  
  Implementation for Precise DPE (PDPE).

## Notes

1. The 20ms coherent integration and common error cancellation by the reference receiver have been implemented to enable precise positioning.
2. Due to time constraints, only a few solutions for specific time instances are provided. For example, `01_02` refers to the first two solutions out of a total of 75, with a solution rate of 1Hz.
3. We would like to thank Dr. Adrià Gusi-Amigó and Albora Technologies for their assistance with data collection. The data is not included in this repository to avoid commercial usage. For more details on the data, please refer to the following paper or contact the corresponding author.
	```
     Liu, X., Ribot, M.Á., Gusi-Amigó, A., Closas, P., Garcia, A.R. and Subirana, J.S.,
	 2020, September. RTK feasibility analysis for GNSS snapshot positioning. In Proceedings of the 33rd International Technical Meeting of the Satellite Division of The Institute of Navigation (ION GNSS+ 2020) (pp. 2911-2921).
     ```
     [INO GNSS+ 2020: https://www.ion.org/publications/abstract.cfm?articleID=17768](https://www.ion.org/publications/abstract.cfm?articleID=17768)