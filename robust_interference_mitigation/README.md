For more details or to cite this work, refer to the paper:  
```
Closas, P. and Gusi-Amigo, A.,
2017. Direct position estimation of GNSS receivers: Analyzing main results, architectures, enhancements, and challenges. IEEE Signal Processing Magazine, 34(5), pp.72-84.  
```  
[IEEE SPM 2017: https://ieeexplore.ieee.org/abstract/document/8026199](https://ieeexplore.ieee.org/abstract/document/8026199)  

This project can be separated into two parts:

1. LoE (Loss of Efficiency) calculation to evaluate the performance of robust DPE in the abscence of interference threats. The calculation restuls are included in the figure below, where we can see that the robust DPE share the same performance as the robust 2SP method and LoE of DD-RIM (Dual Domain Robust Interference Mitigation) techniques is twice as the one of SD-RIM (Single Domain Robust Interference Mitigation) technique:
   ![image](https://github.com/user-attachments/assets/d90a76f4-41a0-445e-bc3c-72c05ee7c7f0)

2. RMSE (Rooted Mean Square Error) calculation to evaluate the performance of robust DPE in the existence of interference threats. The calculation results are included in the figures below. The upper figure shows the RMSE of robust DPE against CW (Continuous Wave) jamming signal, while the lower one shows the RMSE of robust DPE against DME (Distance Measurement Equipment) signal. In the two figures, we can see the efficiency of robust DPE against interference threats compared with the standard DPE:
   ![image](https://github.com/user-attachments/assets/c89f9dc1-1687-47ba-a0ff-009c65ca2c1d)
   ![image](https://github.com/user-attachments/assets/d2e6da6e-be91-4b71-bc8d-0d0ae9ec82e6)
