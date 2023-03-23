Updated for gkvp_f0.62                                  S. Maeyama   March 2023
Updated for gkvp_f0.61                                  S. Maeyama   March 2021
Updated for gkvp_f0.55                                   M. Nakata     Dec 2019
Updated for gkvp_f0.50                                  S. Maeyama     Sep 2017
Updated for gkvp_f0.48                                  S. Maeyama     Dec 2016
Updated for gkvp_f0.47                                  S. Maeyama     Nov 2016
Updated for gkvp_f0.46                                  S. Maeyama     May 2016
Updated for gkvp_f0.45                                   M. Nakata    July 2015
Updated for gkvp_f0.40                                   M. Nakata    June 2014
NOTE for gkvp_f0.30                                     S. Maeyama   March 2013

%%% How to run the code %%%

1. make
2. ./shoot start_num end_num (JOB_ID)

   Examples: Single job submission (First, *.001)       -  ./shoot 1 1
             Single job submission (Second, *.002)      -  ./shoot 2 2
             Step job submission (from Third to Fifth)  -  ./shoot 3 5
             Sequential job submission                  -  ./shoot 6 7 11223
             (It is the case that there is a Fifth job in queue having 
              JOB_ID=11223, and you continue step jobs from Sixth to Seventh)



%%% For multi-platforms %%%

1. Create "shoot", "Makefile" and "sub.q" in "run/", which have already been
   prepared for helios(IFERC), k(RIKEN), ps(NIFS), nu(Nagoya), oakleaf(Tokyo).
2. Set the directory for output "DIR" in "shoot".
3. Set node number, elapsed time limit, and so on in "sub.q".
4. Compile and run the code.


%%% For numerical and physical settings %%%

Grid number and mpi process number are defined in "gkvp_f0.50_header.f90".

calc_type: "linear"    -  for linear runs
           "lin_freq"  -  for linear runs with frequency check
           "nonlinear" -  for nonlinear runs

z_bound: "zerofixed"   -  fixed boundary in zz
         "outflow"     -  outflow boundary in zz
         "mixed"       -  outflow boundary in zz only for ff

z_filt: "on"           -  enable 4th-order filtering in zz  
        "off"          -  disable it 

z_calc: "cf4"          -  4th-order central finite difference for df/dz (nzb=2)
        "up5"          -  5th-order upwind finite difference for df/dz (nzb=3)

art_diff:              -  coefficient of artificial diffusion for z_calc=cf4

init_random:           - Switch on/off random number for initialization.

num_triad_diag:        - Number of triad transfer diagnostics, which should be 
                         consistent with the number of "&triad mxt=**,myt=**/".

&triad mxt=**,myt=**/  - Diagnosed mode number of triad transfer analysis.
                         Add lines of "&triad mxt=**,myt=**/" as desired.

equib_type: "analytic" -  Analytic helical field with the metrics in cylinder
            "s-alpha"  -  s-alpha model with alpha = 0 (cylindrical metrics)
       "s-alpha-shift" -  s-alpha model with Shafranov shift
            "circ-MHD" -  Concentric circular field with the consistent metrics 
            "vmec"     -  Tokamak/stellarator field from the VMEC code
            "eqdsk"    -  Tokamak field (MEUDAS/TOPICS or G-EQDSK) via IGS code
            "slab"     -  Shearless slab geometry
            "ring"     -  Ring dipole geometry

inum: current shot number

ch_res: Change perpendicular resolutions (Settings are somewhat complicated.)

f_log: data directory for log data  
f_hst: data directory for time-series data  
f_phi: data directory for field quantity data   
f_fxv: data directory for distribution function data  
f_cnt: data directory for continue data  

e_limit: Elapsed time limit [sec]

tend: End of simulation time [L_ref/v_ref]

dtout_fxv, dtout_ptn, dtout_eng: Time spacing for data output

dtout_dtc: Time spacing for time-step-size adaption

dt_max: Maximum time step size

adapt_dt: Time-step-size adaption
          (If adapt_dt = .false., time step size is set to be dt_max.)

courant_num: courant number for time-step-size adaption.

time_advnc: "rkg4"      - Explicit time integration by 4th-order Runge-Kutta-Gill method
            "imp_colli" - 2nd-order operator split + 2nd-order implicit collision solver + 4th-order RKG method for collisionless physics
            "auto_init" - If collision restrict linear time step size, time_advnc="imp_colli". Otherwise, time_advnc="rkg4"

!!! NOTE THAT L_ref is set to the major radius on magnetic axis, Raxi, on GKV. !!! 
!!! NOTE THAT m_ref is set to the proton mass, mp, on GKV. !!! 
!!! NOTE THAT e_ref is set to the elementaly charge, e, on GKV. !!! 
!!! NOTE THAT n_ref is set to the local electron density, ne, on GKV. !!! 
!!! NOTE THAT T_ref is set to the first ion temperature, Ti, on GKV. !!! 
!!! NOTE THAT rho_ref is set to the thermal proton gyroradius, rho_tp (with v_ref = vtp = sqrt(Ti/mp)), on GKV. !!! 

R0_Ln: normalized density gradient parameter, L_ref/L_ne, L_ref/L_ni, ...

R0_Lt: normalized temperature gradient parameterm, L_ref/L_te, L_ref/L_ti, ...

nu: bias factor for LB collision model, e.g., 1.d0, 0.5d0, 2.d0, ...
!!! NOTE THAT, after ver. f0.40, the collision frequency is consistently calculated by (Nref, Tref, Lref), 
     and nu is just used as a bias factor only for LB case. Also, nu is not used in multi-species collisions (full). !!!    

Anum: Mass number, m_e/m_ref, m_i/m_ref, ...

Znum: Atomic number, |e_e/e_ref|, |e_i/e_ref|, ...

fcs: charge fraction, |e_e*n_e/(e_ref*n_ref)|, |e_i*n_i/(e_ref*n_ref)|, ...
!!! NOTE THAT fcs for electron shoud be 1.0, and the summation of fcs over all ion species should also be 1.0. !!!

sgn: sign of chaege, e_e/|e_e|, e_i/|e_i|, ...

tau: normalized temperature, T_e/T_ref, T_i/T_ref, ...
!!! NOTE THAT T_i/T_ref = 1 for the first ion species, because T_ref = T_i(first ion) in GKV!!! 

dns1: initial perturbation amplitude, (L_ref/rho_ref)*delta-n_e/n_ref, (L_ref/rho_ref)*delta-n_i/n_ref, ...

tau_ad: Ti/Te for single species ITG-ae(sgn=+1), Te/Ti for single species ETG-ai(sgn=-1) 

lambda_i: ratio of (Debey_length/rho_ref)^2 = epsilon_0*B_ref**2/(m_ref*n_ref)

beta: local beta value evaluated with B_ref, mu_0*n_ref*T_ref/B_ref*2

ibprime:   "1"  -  enable a grad-p (finete beta-prime) contribution on the magnetic drift kvd for equib_type = eqdsk
           "0"  -  ignore it 

vmax: velocity domain size in the unit of each thermal speed vts, max[v/vts]

nx0: the radial mode number assigned for the initial perturbation 
!!! NOTE that if nx0 exceeds nx, nx0 is reset to nx. 
       A sufficiently large value, thus, gives an uniform pertubation for entire kx-modes. !!!    

mach: not used in f0.55
uprime: not used in f0.55

gamma_e: mean ExB shearing rate defined by the 2nd-order radial derivative: (1/B_ref)*d^2(Phi)/dx^2 / (V_ref/L_ref) (at x=0: fluxtube center), 
         where Phi(x) is the equilibrium electrostatic potential.

ntheta: the length of fluxtube, zz-domain = +/-pi times ntheta

kymin: minimum poloidal wave number

m_j: mode connection number in fluxtube model, kxmin = |2*pi*s_hat*kymin/m_j|

del_c: mode connection phase in fluxtube model

eps_r ~~ malpha : geometrical parameters such as safety factor, B-shear, etc.  

&ring: parameters for ring dipole geometry
       !  There is a ring current at R=a. The field line passing through (R,Z)=(R0,0) is picked up as a flux-tube domain.
       !  The reference length is set to be R0 (not the ring current at R=a).
       !  The reference magnetic field strength is B0 at (R,Z)=(R0,0).

ring_a: = a / R0, which specify a flux tube of the ring dipole.

kxmin: Minimum wavenumber in kx, valid only when equib_type == "ring"

&vmecp -- &bozxf : parameters for vmec equilibrium

&igsp -- &igsf : parameters for tokamak (g-eqdsk) equilibrium

s_input: reference radial flux surface, rho

mc_type:   "0"  -  Axisymmetric
           "1"  -  Boozer
           "2"  -  Hamada

q_type:    "1"  -  use consistent q-value on g-eqdsk equilibrium (Recommended)
           "0"  -  use inconsistent, but given q_0 value described above.

nss: the number of radial grids on METRIC data
ntheta: (the number of poloidal grids on METRIC data) + 1 = global_nz*2 + 1

f_igs: file location of METRIC data produced by IGS code

&nu_ref: parameters for collisions  

Nref: local electron density in m^-3
Lref: reference length (= Raxi) in m 
Tref: main ion temperature in keV

col_type: "LB"      -  Lenard-Bernstein type collision operator 
          "lorentz" -  Lorentz model collision operator  
          "full"    -  multi-species linearized collision operator  

iFLR:     "1"     -  enable the FLR terms (for LB and full)
          "0"     -  disable it (DK-limit)

icheck:   "0"     -  for production runs   
          "1"     -  debug test with Maxwellian Annihilation (should be used with IFLR = 0)


