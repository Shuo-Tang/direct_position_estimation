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
