# Benchmarks  

This benchmark directory contains some simulation resuls of GKV.  

You will be able to reproduce these results by appropriately setting the resolutions in "gkvp_header.f90" and physical parameters in "gkvp_namelist". Note that all benchmark chooses *init_random=.false.* to be independent of numerical random number, conventional simulations of GKV employs *init_random=.true.*   


---
## Benchmark - SmallGridIO  
#### Description  
Check standard I/O by a quick computation with small grid numbers.  
Since the grid numbers are too small, the results have no physical meanings.  
#### Environmental informations on this benchmark  
- Benchmark date:  Oct. 12, 2020  
- GKV+ version: gkvp_f0.58  
- Architecture: Fugaku (RIKEN R-CCS)  
- Enviromental setting: -     


---
## Benchmark - ITGae-lin  
#### Description  
Linear ITG with adiabatic electron (nprocs=1)  
Cyclone-base-case parameters (R0_Ln=2.22,R0_Lt=6.92,eps_r=0.18,q_0=1.4,s_hat=0.8)  
Electrostatic beta=0, equib_type="analytic", LB collision  
Settings from GKV hands-on seminar at NIFS on Dec 15, 2017.  
#### Environmental informations on this benchmark  
- Benchmark date:  Oct. 12, 2020  
- GKV+ version: gkvp_f0.58  
- Architecture: Fugaku (RIKEN R-CCS, Kobe, Japan)  
- Enviromental setting: -   
#### Results: linear dispersion  
k_y [1/rho_i], frequency [v_ti/R], growthrate [v_ti/R]  
   0.05    -0.101     0.030  
   0.10    -0.219     0.079  
   0.15    -0.354     0.144  
   0.20    -0.498     0.200  
   0.25    -0.645     0.236  
   0.30    -0.790     0.249  
   0.35    -0.928     0.238  
   0.40    -1.056     0.208  
   0.45    -1.173     0.162  
   0.50    -1.278     0.102  
   0.55    -1.370     0.311  
   0.60    -1.453    -0.488  


---
## Benchmark - ITGae-nl  
#### Description  
Nonlinear ITG with adiabatic electron (nprocs=1)  
Cyclone-base-case parameters (R0_Ln=2.22,R0_Lt=6.92,eps_r=0.18,q_0=1.4,s_hat=0.8)  
Electrostatic beta=0, equib_type="analytic", LB collision  
Settings from GKV hands-on seminar at NIFS on Dec 15, 2017.  
#### Environmental informations on this benchmark  
- Benchmark date:  Oct. 12, 2020  
- GKV+ version: gkvp_f0.58  
- Architecture: Fugaku (RIKEN R-CCS, Kobe, Japan)  
- Enviromental setting: -    


---
## Benchmark - ITGke-lin   
#### Description  
Linear ITG with kinetic electron (nprocs=2)  
Cyclone-base-case parameters (R0_Ln=2.22,R0_Lt=6.92,eps_r=0.18,q_0=1.4,s_hat=0.8)  
Electromagnetic beta=0.001, equib_type="analytic", LB collision  
#### Environmental informations on this benchmark  
- Benchmark date:  Oct. 12, 2020  
- GKV+ version: gkvp_f0.58  
- Architecture: Fugaku (RIKEN R-CCS)  
- Enviromental setting: -  
#### Results: linear dispersion  
k_y [1/rho_i], frequency [v_ti/R], growthrate [v_ti/R]  
   0.10    -0.282     0.121  
   0.15    -0.459     0.240  
   0.20    -0.641     0.345  
   0.25    -0.817     0.413  
   0.30    -0.985     0.439  
   0.35    -1.140     0.427  
   0.40    -1.279     0.385  
   0.45    -1.400     0.319  
   0.50    -1.504     0.239  


---
## Benchmark - ITGke-nl  
#### Description  
Nonlinear ITG with kinetic electron (nprocs=1)  
Cyclone-base-case parameters (R0_Ln=2.22,R0_Lt=6.92,eps_r=0.18,q_0=1.4,s_hat=0.8)  
Electromagnetic beta=0.001, equib_type="analytic", LB collision  
#### Environmental informations on this benchmark  
- Benchmark date: Oct. 12, 2020  
- GKV+ version: gkvp_f0.58  
- Architecture: Fugaku (RIKEN R-CCS, Kobe, Japan)  
- Enviromental setting: -  


