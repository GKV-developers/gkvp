MODULE GKV_geom
!-------------------------------------------------------------------------------
!
!    Calculate geometric constants  
!
!    Update history of gkvp_geom.f90
!    --------------
!      gkvp_f0.62 (S. Maeyama, Mar 2023)
!        - First implementation.
!        - Geometric constants, which had been set in gkvp_set.f90, are moved.
!          This module will be called from gkvp_set.f90 for initialization,
!          and from gkvp_advnc.f90 for update in rotating flux-tube model.
!        - Subroutines geom_* are public, can be called from other module.
!        - Subroutines metric_* are private, treating metric structure.
!        - equib_type = "ring" is added for ring dipole geometry.
!
!-------------------------------------------------------------------------------

  use GKV_header
  use GKV_mpienv
  use GKV_math,   only: math_j0, math_j1, math_j2, math_g0
  use GKV_intgrl, only: intgrl_fsrf, intgrl_v0_moment_ms
! for vmec equilibrium w/ Booz_xform by M. Nakata & M. Nunami  (Aug. 2016)
  use GKV_vmecbzx, only: vmecbzx_boozx_read, vmecbzx_boozx_coeff
! for tokamak(eqdsk) equilibrium
  use GKV_igs,    only: igs_read, igs_coeff
  !sakano_ring-dipole st 202303
  use GKV_ring,   only: ring_coordinates
  !sakano_ring-dipole end 202303

  implicit none

  private

  public   geom_read_nml,      & ! Called once in GKV_set. Read namelist.
           geom_init_kxkyzvm,  & ! Called once in GKV_set. Set time-indep. grids.
           geom_init_metric,   & ! Called once in GKV_set. Set metrics at t=0.
           geom_set_operators, & ! Set operators (e.g., ksq) by using metrics at t. 
           geom_reset_time,    & ! Called once in GKV_set. Reset metrics and operators at a given time.
           geom_increment_time   ! Called in GKV_advnc. Increment metrics and operators by a given time step.
           

  type metric_global
     ! Global metrics at t=0 are stored.
     ! Metrics in GKV coordinates (x,y,z)
     ! Metrics in flux coordinates (r,t,q)=(rho,theta,zeta)
      real(kind=DP), dimension(-global_nz:global_nz-1) :: zz    ! The rotating flux tube coordinate (= z'')
      real(kind=DP), dimension(-global_nz:global_nz-1) :: theta ! The poloidal angle theta_pol, not the flux-coordinate theta
      real(kind=DP), dimension(-global_nz:global_nz-1) :: omg   ! Magnetic field strength
      real(kind=DP), dimension(-global_nz:global_nz-1) :: &
          domgdx, domgdy, domgdz, gxx, gxy, gxz, gyy, gyz, gzz, rootg_xyz
      real(kind=DP), dimension(-global_nz:global_nz-1) :: &
          domgdr, domgdt, domgdq, grr, grt, grq, gtt, gtq, gqq, rootg_rtq
    contains
      procedure :: init => metric_global_init
      procedure :: xyz2rtq => metric_global_xyz2rtq
      procedure :: rtq2xyz => metric_global_rtq2xyz
  end type

  type metric_fourier
     ! Metrics in flux coordinates at t=0, stored in Fourier coefficient
     !   e.g.,   fourier_omg(kz) = \int omg(z) * exp(-i*kz*z) dz / \int dz
     !           omg(z) = \sum_k fourier_omg(kz) * exp(i*kz*z)
     ! Thus, omg(z) at arbitrary z is obtained by Fourier interpolation.
      real(kind=DP), dimension(-global_nz:global_nz-1) :: kz
      complex(kind=DP), dimension(-global_nz:global_nz-1) :: theta_tilde, omg
      complex(kind=DP), dimension(-global_nz:global_nz-1) :: &
          domgdr, domgdt, domgdq, grr, grt, grq, gtt, gtq, gqq, rootg_rtq
    contains
      procedure :: init => metric_fourier_init
      procedure :: dft_rtq2coef => metric_fourier_dft_rtq2coef
  end type

  type metric_local
     ! Local metrics at any time t are stored.
     ! They are updated with time integration, and used for solving
     ! Gyrokinetic equation in the rotating flux tube model.
      real(kind=DP), dimension(-nz:nz-1) :: zz          ! The rotating flux tube coordinate (= z'')
      real(kind=DP), dimension(-nz:nz-1) :: zz_labframe ! The flux-coordinate theta in the lab frame (= z''+t*gamma_e/s_hat)
      real(kind=DP), dimension(-nz:nz-1) :: theta ! The geometrical poloidal angle theta_pol, not the flux-coordinate theta
      real(kind=DP), dimension(-nz:nz-1) :: omg   ! Magnetic field strength
      real(kind=DP), dimension(-nz:nz-1) :: &
          domgdx, domgdy, domgdz, gxx, gxy, gxz, gyy, gyz, gzz, rootg_xyz
      real(kind=DP), dimension(-nz:nz-1) :: &
          domgdr, domgdt, domgdq, grr, grt, grq, gtt, gtq, gqq, rootg_rtq
    contains
      procedure :: copy_global => metric_local_copy_global
      procedure :: init => metric_local_init
      procedure :: update => metric_local_update
      procedure :: dft_coef2rtq => metric_local_dft_coef2rtq
      procedure :: rtq2xyz => metric_local_rtq2xyz
  end type

  type(metric_global),  save :: mtr_global
  type(metric_fourier), save :: mtr_fourier
  type(metric_local),   save :: mtr_local

  real(kind=DP), save :: cx, cy, cb


! for s-alpha model with Shafranov shift
    real(kind=DP) ::  p_total, dp_totaldx, beta_total, alpha_MHD

    real(kind=DP) :: r_major

    integer, parameter :: num_omtr = 13
!   real(kind=DP) :: metric_l(1:num_omtr,-nz:nz-1), metric_g(1:num_omtr,-global_nz:global_nz-1)

    real(kind=DP) :: s_hat

    real(kind=DP) :: eps_r

    real(kind=DP) :: lz, kxmin, kymin, dz, mmax, dm, del_c
    real(kind=DP) :: z0, z0_l
    integer       :: n_tht, m_j

    real(kind=DP) :: rdeps00, eps_hor, lprd, mprd, lmmq, malpha
    real(kind=DP) :: eps_mor, eps_por, lprdm1, lprdp1, lmmqm1, lmmqp1
    real(kind=DP) :: eps_rnew, rdeps1_0, rdeps1_10, rdeps2_10, rdeps3_10

    real(kind=DP) :: s_input, s_0      ! radial label of fluxtube center 
    integer       :: mc_type           ! 0:Axisym., 1:Boozer, 2:Hamada
    integer       :: q_type            ! 0:use q and s_hat value in confp, 1:calclated by IGS
    integer       :: isw, nss, ntheta, nzeta
    real(kind=DP) :: phi_ax            ! axisymmetric toroidal angle 

!sakano_ring-dipole st 202303
    real(kind=DP) :: ring_a
!sakano_ring-dipole end 202303

    real(kind=DP) :: lz_l


CONTAINS


!--------------------------------------
  SUBROUTINE geom_read_nml
!--------------------------------------
    implicit none

    real(kind=DP) :: theta
    real(kind=DP), dimension(0:ns-1) :: eta
    real(kind=DP) :: domgdx, domgdy, domgdz
    real(kind=DP), dimension(1:3,1:3) :: gg
    integer :: iz, is, isw


    namelist /physp/ R0_Ln,  &    ! R0/Lns
                     R0_Lt,  &    ! R0/Lts
                        nu,  &    ! factor for collision freq. in LB model    
                      Anum,  &    ! mass number
                      Znum,  &    ! charge number 
                       fcs,  &    ! charge-density fraction 
                       sgn,  &    ! signs of charge 
                       tau,  &    ! T-ratio Ts/T0, T0=reference ion temp. of ranks=1
                      dns1,  &    ! initial perturbation amplitude
                    tau_ad,  &    ! Ti/Te for ITG-ae, Te/Ti for ETG-ai
                  lambda_i,  &    ! (Debye/rho_tp)^2 
                      beta,  &    ! mu0*ni*Ti/B^2
                   ibprime,  &    ! flag for finite beta-prime effect on kvd
                      vmax,  &    ! maximum v_para in unit of v_ts
                       nx0        ! mode number for the initial perturbation

    namelist /rotat/ mach, uprime, gamma_e

    namelist /nperi/ n_tht, kymin, m_j, del_c
    namelist /confp/ eps_r, eps_rnew,                       &
                     q_0, s_hat,                            &
                     lprd, mprd, eps_hor, eps_mor, eps_por, &
                     rdeps00, rdeps1_0, rdeps1_10,          & 
                     rdeps2_10, rdeps3_10, malpha
!    namelist /vmecp/ q_0, rad_a,                            &
!                     R0_unit, r_edge,                       &
!                     b0b00, alpha_fix
    namelist /vmecp/ s_input, nss, ntheta, nzeta

    namelist /igsp/ s_input, mc_type, q_type, nss, ntheta
!sakano_ring-dipole st 202303
    namelist /ring/ ring_a, kxmin
!sakano_ring-dipole end 202303

      tau(:)   = 1.0_DP
      nu(:)    = 0.002_DP
      R0_Ln(:) = 2.5_DP
      R0_Lt(:) = 7.5_DP


      read(inml,nml=physp)


        do is = 0, ns-1
          if( R0_Ln(is) /= 0._DP ) then
            eta(is) = R0_Lt(is) / R0_Ln(is)
          else
            eta(is) = 1.d+20
          end if
        end do


        write( olog, * ) " # Physical parameters"
        write( olog, * ) ""
        write( olog, * ) " # r_major/L_ns = ", R0_Ln(:)
        write( olog, * ) " # r_major/L_ts = ", R0_Lt(:)
        write( olog, * ) " # eta          = ", eta(:)
        write( olog, * ) " # nu           = ", nu(:)
        write( olog, * ) " # A-number     = ", Anum(:)
        write( olog, * ) " # Z-number     = ", Znum(:)
        write( olog, * ) " # fcs          = ", fcs(:)
        write( olog, * ) " # sgn          = ", sgn(:)
        write( olog, * ) " # tau          = ", tau(:)
        write( olog, * ) " # dns1         = ", dns1(:)
        write( olog, * ) " # tau_ad       = ", tau_ad
        write( olog, * ) " # lambda_i^2   = ", lambda_i
        write( olog, * ) " # beta_i       = ", beta
        write( olog, * ) " # ibprime      = ", ibprime
        write( olog, * ) " # nx0          = ", nx0
        write( olog, * ) ""


      mach = 0._DP
      uprime = 0._DP
      gamma_e = 0._DP

      read(inml,nml=rotat)

        write( olog, * ) " # Mean rotation parameters"
        write( olog, * ) ""
        write( olog, * ) " # Mach number  = ", mach
        write( olog, * ) " # uptime       = ", uprime
        write( olog, * ) " # gamma_ExB    = ", gamma_e
        write( olog, * ) ""


      n_tht = 1

      read(inml,nml=nperi)


      if( trim(equib_type) == "slab") then

        read(inml,nml=confp)

        lprdm1   = 0._DP
        lprdp1   = 0._DP

        lmmq     = 0._DP
        lmmqm1   = 0._DP
        lmmqp1   = 0._DP

        q_0      = 1._DP ! For now, fixed q_0=1. Changing q_0 can extend parallel z-box size.
        s_hat    = 0._DP ! only shear less slab
        eps_r    = 1._DP

        eps_hor  = 0._DP
        lprd     = 0._DP
        mprd     = 0._DP
        malpha   = 0._DP

        rdeps00  = 0._DP
        eps_mor  = 0._DP
        eps_por  = 0._DP

        write( olog, * ) " # Configuration parameters"
        write( olog, * ) ""
        write( olog, * ) " # q_0          = ", q_0
        write( olog, * ) " # s_hat        = ", s_hat
        write( olog, * ) " # eps_r        = ", eps_r
        write( olog, * ) ""

        write( olog, * ) " # eps_hor      = ", eps_hor
        write( olog, * ) " # lprd         = ", lprd
        write( olog, * ) " # mprd         = ", mprd
        write( olog, * ) " # malpha       = ", malpha
        write( olog, * ) " # rdeps00      = ", rdeps00

        write( olog, * ) " # eps_mor      = ", eps_mor
        write( olog, * ) " # lprdm1       = ", lprdm1
        write( olog, * ) " # eps_por      = ", eps_por
        write( olog, * ) " # lprdp1       = ", lprdp1
        write( olog, * ) ""

      else if( trim(equib_type) == "analytic"  .OR.  &
               trim(equib_type) == "s-alpha"   .OR.  &
               trim(equib_type) == "s-alpha-shift"   .OR.  &
               trim(equib_type) == "circ-MHD" ) then


        read(inml,nml=confp)


        lprdm1   = lprd - 1.0_DP
        lprdp1   = lprd + 1.0_DP

        lmmq     = lprd   - mprd * q_0
        lmmqm1   = lprdm1 - mprd * q_0
        lmmqp1   = lprdp1 - mprd * q_0


        write( olog, * ) " # Configuration parameters"
        write( olog, * ) ""
        write( olog, * ) " # q_0          = ", q_0
        write( olog, * ) " # s_hat        = ", s_hat
        write( olog, * ) " # eps_r        = ", eps_r
        write( olog, * ) ""

        write( olog, * ) " # eps_hor      = ", eps_hor
        write( olog, * ) " # lprd         = ", lprd
        write( olog, * ) " # mprd         = ", mprd
        write( olog, * ) " # malpha       = ", malpha
        write( olog, * ) " # rdeps00      = ", rdeps00

        write( olog, * ) " # eps_mor      = ", eps_mor
        write( olog, * ) " # lprdm1       = ", lprdm1
        write( olog, * ) " # eps_por      = ", eps_por
        write( olog, * ) " # lprdp1       = ", lprdp1
        write( olog, * ) ""


     else if( trim(equib_type) == "vmec" ) then


        read(inml,nml=confp)

        read(inml,nml=vmecp)

        call vmecbzx_boozx_read( nss, ntheta, nzeta )

         isw = 0
         iz = 0
         call vmecbzx_boozx_coeff( isw,  nss,  ntheta,  nzeta,  s_input, iz, 0._DP,  lz_l,   &  ! input 
                           s_0,           q_0,     s_hat,    eps_r,  phi_ax,             &  ! output
                           omg(iz), rootg(iz),    domgdx,   domgdz,  domgdy,             &
                           gg(1,1),   gg(1,2),   gg(1,3),  gg(2,2),                      &
                           gg(2,3),   gg(3,3)  )

        write( olog, * ) " # Configuration parameters"
        write( olog, * ) ""
        write( olog, * ) " # r_major/L_ns = ", R0_Ln(:)
        write( olog, * ) " # r_major/L_ts = ", R0_Lt(:)
        write( olog, * ) " # eta          = ", eta(:)
        write( olog, * ) " # q_0          = ", q_0
        write( olog, * ) " # s_hat        = ", s_hat
        write( olog, * ) " # eps_r        = ", eps_r
        write( olog, * ) " # s_input, s_0 = ", s_input, s_0
        write( olog, * ) " # nss, ntheta, nzeta  = ", nss, ntheta, nzeta


     else if( trim(equib_type) == "eqdsk" ) then


        read(inml,nml=confp)

        read(inml,nml=igsp)

        call igs_read( mc_type, nss, ntheta )

        if ( q_type == 1 ) then
         isw = 0
         iz = 0
         call igs_coeff( isw,  mc_type,   nss,    ntheta,  s_input,  0._DP, lz_l,   &  ! input 
                         s_0,       q_0,     s_hat,    eps_r,  theta,               &  ! output
                           omg(iz), rootg(iz),    domgdx,   domgdz, domgdy,         &
                           gg(1,1),   gg(1,2),   gg(1,3),  gg(2,2),                 &
                           gg(2,3),   gg(3,3)  )
        end if

        write( olog, * ) " # Configuration parameters"
        write( olog, * ) ""
        write( olog, * ) " # r_major/L_ns = ", R0_Ln(:)
        write( olog, * ) " # r_major/L_ts = ", R0_Lt(:)
        write( olog, * ) " # eta          = ", eta(:)
        write( olog, * ) " # q_0          = ", q_0
        write( olog, * ) " # s_hat        = ", s_hat
        write( olog, * ) " # eps_r        = ", eps_r
        write( olog, * ) " # s_input, s_0 = ", s_input, s_0
        write( olog, * ) " # nss, ntheta  = ", nss, ntheta

!sakano_ring-dipole st 202303
      else if ( trim(equib_type) == "ring" ) then

        read(inml,nml=confp)

        read(inml,nml=ring)

        s_hat = 0._DP

        write( olog, * ) " # Configuration parameters for ring dipole configuration"
        write( olog, * ) ""
        write( olog, * ) " # s_hat        = ", s_hat
        write( olog, * ) " # kxmin        = ", kxmin
        write( olog, * ) " # ring_a       = ", ring_a
        write( olog, * ) " # eps_r        = ", eps_r
        write( olog, * ) " # q_0          = ", q_0
!sakano_ring-dipole end 202303

      else

        write( olog, * ) " # wrong choice of the equilibrium "
        call flush(olog)
        call MPI_Finalize(ierr_mpi)
        stop

      end if

  END SUBROUTINE geom_read_nml


!--------------------------------------
  SUBROUTINE geom_init_kxkyzvm(lx, ly, eps_r_temp)
!--------------------------------------
    implicit none
    real(kind=DP), intent(out) :: lx, ly, eps_r_temp
    integer       :: global_iv, global_im
    integer       :: mx, my, iz, iv, im

      eps_r_temp = eps_r

      if ( trim(equib_type) /= "ring" ) then
        if (abs(s_hat) < 1.d-10) then ! When s_hat == ZERO
          m_j = 0
          kxmin = kymin
        else if (m_j == 0) then
          kxmin = kymin
        else
          kxmin    = abs(2._DP * pi * s_hat * kymin / real(m_j, kind=DP))
        end if
      end if
      lx       = pi / kxmin
      ly       = pi / kymin
                    ! kymin=pi/ly=pi/[r_minor*pi/(q0*n_alp)]=q0*n_alp/r_minor

      lz       = real( n_tht, kind=DP ) * pi        ! total z-length
      lz_l     = lz / real( nprocz, kind=DP )       ! local z-length

      do mx = -nx, nx
        kx(mx)   = kxmin * real( mx, kind=DP )
      end do

      ky(:) = 0._DP
      do my = ist_y_g, iend_y_g
        ky(my-ist_y_g)   = kymin * real( my, kind=DP )
      end do

      kxmin_g = kxmin
      kymin_g = kymin 

      z0       = - lz                     ! global lower boundary
      z0_l     = 2._DP * lz_l * real( rankz, kind=DP ) + z0
                                          ! local lower boundary

      dz       = lz_l / real( nz, kind=DP )

      do iz = -nz, nz-1
        zz(iz)   = dz * real( iz + nz, kind=DP ) + z0_l
      end do


      dv       = 2._DP * vmax / real( 2 * nv * nprocv -1, kind=DP )

      do iv = 1, 2*nv
        global_iv = 2 * nv * rankv + iv
        vl(iv)    = dv * ( real( global_iv - nv * nprocv - 1, kind=DP ) + 0.5_DP )
      end do
                                          ! --- debug
                                          !   write( olog, * ) " *** iv, vl "
                                          !   do iv = 1, 2*nv
                                          !     global_iv = 2 * nv * rankv + iv
                                          !     write( olog, * ) iv, global_iv, vl(iv)
                                          !   end do
                                          !   write( olog, * ) ""

      mmax     = vmax
      dm       = mmax / real( nprocm * ( nm+1 ) - 1, kind=DP )
                                          ! --- equal spacing in vperp

      do im = 0, nm
        global_im = ( nm+1 ) * rankm + im
        mu(im)    = 0.5_DP * ( dm * real( global_im, kind=DP ) )**2
      end do


      do my = ist_y_g, iend_y_g
        ck(my-ist_y_g)   = exp( ui * 2._DP * pi * del_c &
                                   * real( n_tht * my, kind=DP ) )
        dj(my-ist_y_g)   = - m_j * n_tht * my
                                          !  del_c = q_0*n_alp-int(q_0*n_alp)
                                          !  m_j   = 2*n_alp*q_d
      end do


      write( olog, * ) " # Numerical parameters"
      write( olog, * ) ""
      write( olog, * ) " # n_tht = ", n_tht
      write( olog, * ) " # lx, ly, lz   = ", lx, ly, lz
      write( olog, * ) " # lz,   z0     = ", lz, z0
      write( olog, * ) " # lz_l, z0_l   = ", lz_l, z0_l
      write( olog, * ) " # kxmin, kymin = ", kxmin, kymin
      write( olog, * ) " # kxmax, kymax = ", kxmin*nx, kymin*global_ny
      write( olog, * ) " # kperp_max    = ", sqrt((kxmin*nx)**2+(kymin*global_ny)**2)
      write( olog, * ) " # m_j, del_c   = ", m_j, del_c
      write( olog, * ) " # dz           = ", dz
      write( olog, * ) " # dv, vmax     = ", dv, vmax
      write( olog, * ) " # dm, mmax     = ", dm, mmax
      write( olog, * ) ""

      if (gamma_e == 0._DP) then
        tlim_exb = 999999.d0
      else
        tlim_exb = (kxmin*(nx-nx0))/(kymin*global_ny*abs(gamma_e))
      end if
      write( olog, * ) " # ExB limit time tlim_exb  = ", tlim_exb
      write( olog, * ) " # for (mx=nx0,my=global_ny) initial perturbation: "
      write( olog, * ) " # tlim_exb = kxmin*(nx-nx0)/(kymax*|gamma_e|)"
      write( olog, * ) ""

  END SUBROUTINE geom_init_kxkyzvm


!--------------------------------------
  SUBROUTINE geom_init_metric
!--------------------------------------
    real(kind=DP) :: r_0
    real(kind=DP) :: wzz, theta, gomg
    real(kind=DP) :: gdomgdx, gdomgdy, gdomgdz, &
                     ggxx, ggxy, ggxz, ggyy, ggyz, ggzz, grootg_xyz
    real(kind=DP) :: gdomgdr, gdomgdt, gdomgdq, &
                     ggrr, ggrt, ggrq, ggtt, ggtq, ggqq, grootg_rtq
    integer :: iz, is
!sakano_ring-dipole st 202303
    real(kind=DP) :: ub_dot_grdb, ub_crs_grdb
!sakano_ring-dipole end 202303

      s_hat_g = s_hat

      !- zero clear -
      gdomgdr = 0._DP
      gdomgdt = 0._DP
      gdomgdq = 0._DP
      ggrr = 1._DP
      ggrt = 0._DP
      ggrq = 0._DP
      ggtt = 1._DP
      ggtq = 0._DP
      ggqq = 1._DP
      grootg_rtq = 1._DP

      do iz = -global_nz, global_nz-1

        wzz = dz * iz

        if ( trim(equib_type) == "slab") then
         !- Shearless slab geometry-
         !  Consider translationally symmetric flux surface
         !  (r,t,q)=(x_car,y_car,z_car).
         !  GKV coordinates are
         !    x = x_car,              -lx<=x/gyroradius<lx
         !    y = y_car,              -ly<=y/gyroradius<ly
         !    z = z_car/q_0*r_major,  -n_tht*pi<=z<n_tht*pi
         !  z is dimensionless. 
         !  Length in field-aligned z_car is 2*q_0*lz*r_major.
         !  r_major = 1 in the R0 unit.
         !  Magnetic field is constant B=B0 in z, so omg = 1 in the B0 unit.
         !  Normalized cb = B0/B0 = 1
         !-
          r_major = 1._DP ! in the R0 unit
          cx = 1._DP
          cy = 1._DP
          cb = 1._DP

          q_bar = q_0 
          theta = wzz
          gomg = 1._DP
         !- Metrics in GKV coordinates (x,y,z)
          gdomgdx = 0._DP
          gdomgdy = 0._DP
          gdomgdz = 0._DP
          ggxx = 1._DP
          ggxy = 0._DP
          ggxz = 0._DP
          ggyy = 1._DP
          ggyz = 0._DP
          ggzz = 1._DP/(q_0*r_major)**2
          grootg_xyz = q_0*r_major
         !- Metrics in cartesian (x_car,y_car,z_car)
         !gdomgdr = 0._DP
         !gdomgdt = 0._DP
         !gdomgdq = 0._DP
         !ggrr = 1._DP
         !ggrt = 0._DP
         !ggrq = 0._DP
         !ggtt = 1._DP
         !ggtq = 0._DP
         !ggqq = 1._DP
         !grootg_rtq = 1._DP


        else if( trim(equib_type) == "analytic" ) then
         !- Analytic model of large-aspect-ratio helical system -
         !  [Ref.1] H. Sugama and T.-H. Watanabe, Phys. Plasmas 11, 3068 (2004).
         !  [Ref.2] T.-H.Watanabe, H. Sugama, and S. Ferrando-Margalet,
         !          Nucl. Fusion 47, 1383 (2007).
         !
         !  Consider concentric circular, but helically twisted,
         !  flux surface (r,theta,zeta). 
         !  GKV coordinates are
         !    x = cx*(r-r0)
         !    y = cy*(q(r)*theta-zeta)
         !    z = theta
         !  with cx=1, cy=cx*r0/q0.
         !  In the large-aspect ratio limit, the geometrical length
         !  in the field-aligned direction is dpara=q*r_major*dz.
         !  r_major = 1 in the R0 unit.
         !  Finite aspect ratio eps_r = r0/R0 is retained only in 
         !  magnetic field omg, domgdx, domgdy, domgdz, but not for metrics.
         !  Flux-surface averaged magnetic field is <omg>=B_ax,
         !  where B_ax is the value at the magnetic axis.
         !  cb = (psi_p(r))'/(cx*cy) = B_ax
         !  Normalized <omg> = 1 and cb = 1 in the B_ax unit.
         !-
          r_major = 1._DP       ! Major radius of magnetic axis in the R0 unit
          r_0 = r_major * eps_r ! Minor radius of flux-tube center
          cx = 1._DP
          cy = r_0/q_0
          cb = 1._DP

          q_bar = q_0 
          theta = wzz
          gomg = 1._DP                                          &
               - eps_r * ( cos( wzz )                           &
                       + eps_hor * cos( lmmq   * wzz - malpha ) &
                       + eps_mor * cos( lmmqm1 * wzz - malpha ) &
                       + eps_por * cos( lmmqp1 * wzz - malpha ) )
         !- Metrics in GKV coordinates (x,y,z)
          gdomgdz = eps_r * ( sin(wzz)                                        &
                          + eps_hor * lmmq   * sin( lmmq   * wzz - malpha )   &
                          + eps_mor * lmmqm1 * sin( lmmqm1 * wzz - malpha )   &
                          + eps_por * lmmqp1 * sin( lmmqp1 * wzz - malpha ) )
          gdomgdy = - eps_rnew / r_major * (                              &
                    - ( sin( wzz )                                        &
                      + eps_hor * lprd   * sin( lmmq   * wzz - malpha )   &
                      + eps_mor * lprdm1 * sin( lmmqm1 * wzz - malpha )   &
                      + eps_por * lprdp1 * sin( lmmqp1 * wzz - malpha )   &
                    ) - (-1._DP/eps_r) * gdomgdz )
          gdomgdx = eps_rnew / r_major * (                                &
                    - (                                                   &
                        rdeps00                                           &
                      + rdeps1_0 * cos( wzz )                             &
                      + rdeps2_10 * cos( lmmq   * wzz - malpha )          &
                      + rdeps1_10 * cos( lmmqm1 * wzz - malpha )          &
                      + rdeps3_10 * cos( lmmqp1 * wzz - malpha )          &
                      + s_hat * wzz * ( sin( wzz )                        &
                      + eps_hor * lprd   * sin( lmmq   * wzz - malpha )   &
                      + eps_mor * lprdm1 * sin( lmmqm1 * wzz - malpha )   &
                      + eps_por * lprdp1 * sin( lmmqp1 * wzz - malpha ) ) &
                    ) - (-s_hat*wzz/eps_r) * gdomgdz )
          ggxx = 1._DP
          ggxy = s_hat*wzz
          ggxz = 0._DP
          ggyy = 1._DP + (s_hat*wzz)**2
          ggyz = 1._DP/r_0
          ggzz = 1._DP/r_0**2
          grootg_xyz = q_0*r_major/gomg
         !- Metrics in flux coordinates (r,theta,zeta)
         !ggrr = 1._DP
         !ggrt = 0._DP
         !ggrq = 0._DP
         !ggtt = 1._DP/r_0**2
         !ggtq = 0._DP
         !ggqq = 0._DP ! /=1._DP/r_major**2, because of large-aspect ratio approximation
         !grootg_rtq = r_0*r_major/gomg


        else if( trim(equib_type) == "s-alpha" .or. trim(equib_type) == "s-alpha-shift" ) then
         !- Analytic model of large-aspect-ratio tokamak system -
         !  Consider concentric circular flux surface (r,theta,zeta). 
         !  GKV coordinates are
         !    x = cx*(r-r0)
         !    y = cy*(q(r)*theta-zeta)
         !    z = theta
         !  with cx=1, cy=cx*r0/q0.
         !  In the large-aspect ratio limit, the geometrical length
         !  in the field-aligned direction is dpara=q*r_major*dz.
         !  r_major = 1 in the R0 unit.
         !  Finite aspect ratio eps_r = r0/R0 is retained only in 
         !  magnetic field omg, domgdx, domgdy, domgdz, but not for metrics.
         !  Flux-surface averaged magnetic field is <omg>=B_ax,
         !  where B_ax is the value at the magnetic axis.
         !  cb = (psi_p(r))'/(cx*cy) = B_ax
         !  Normalized <omg> = 1 and cb = 1 in the B_ax unit.
         !-
          r_major = 1._DP       ! Major radius of magnetic axis in the R0 unit
          r_0 = r_major * eps_r ! Minor radius of flux-tube center
          cx = 1._DP
          cy = r_0/q_0
          cb = 1._DP

          if (trim(equib_type) == "s-alpha") then
            !--- s-alpha model without Shafranov shift -
            alpha_MHD = 0._DP
          else if (trim(equib_type) == "s-alpha-shift") then
            !--- s-alpha model with Shafranov shift ----
            p_total = 0._DP
            dp_totaldx = 0._DP
            beta_total = 0._DP
            do is = 0, ns-1
              p_total = p_total + fcs(is) * tau(is) / Znum(is)
              dp_totaldx = dp_totaldx - fcs(is) * tau(is) / Znum(is) * (R0_Ln(is) + R0_Lt(is))
              beta_total = beta_total + 2._DP * beta * fcs(is) * tau(is) / Znum(is)
            end do
            alpha_MHD = - q_0**2 * r_major * beta_total * dp_totaldx / p_total
          end if
          q_bar = q_0
          theta = wzz
          gomg = 1._DP - eps_r * cos( theta ) ! s-alpha with eps-expansion
          !!!!gomg = 1._DP / (1._DP + eps_r * cos( theta )) ! for benchmark
         !- Metrics in GKV coordinates (x,y,z)
          gdomgdz = eps_r * sin( theta )
          !!!!gdomgdz = eps_r * sin( theta ) * omg(iz)**2 ! for benchmark
          gdomgdx = - cos( theta ) / r_major
          gdomgdy = 0._DP
          ggxx = 1._DP
          ggxy = s_hat*wzz - alpha_MHD*sin(wzz) ! with Shafranov shift
          ggxz = 0._DP
          ggyy = 1._DP + (s_hat*wzz - alpha_MHD*sin(wzz))**2 ! with Shafranov shift
          ggyz = 1._DP/r_0
          ggzz = 1._DP/r_0**2
          grootg_xyz = q_0*r_major/gomg
         !- Metrics in flux coordinates (r,theta,zeta)
         !gdomgdr = - cos( theta ) / r_major
         !gdomgdt = eps_r * sin( theta )
         !gdomgdq = 0._DP
         !ggrr = 1._DP
         !ggrt = 0._DP
         !ggrq = 0._DP + sin(wzz)*alpha_MHD*q_0/r_0 ! with Shafranov shift
         !ggtt = 1._DP/r_0**2
         !ggtq = 0._DP
         !ggqq = 0._DP & ! /=1._DP/r_major**2, because of large-aspect ratio approximation
         !     + sin(wzz)*(alpha_MHD*q_0/r_0)**2 ! with Shafranov shift
         !grootg_rtq = r_0*r_major/gomg


        else if( trim(equib_type) == "circ-MHD" ) then
         !- Circular MHD equilibrium -
         !  [Ref.] X. Lapillonne, et al., Phys. Plasmas 16, 032308 (2009).
         !
         !  Consider concentric circular flux surface (r,theta,zeta). 
         !  GKV coordinates are
         !    x = cx*(r-r0)
         !    y = cy*(q(r)*theta-zeta)
         !    z = theta
         !  with cx=1, cy=cx*r0/q0.
         !  In the large-aspect ratio limit, the geometrical length
         !  in the field-aligned direction is dpara=q*r_major*dz.
         !  r_major = 1 in the R0 unit.
         !  In contrast to the s-alpha model, finite aspect ratio eps_r = r0/R0
         !  is retained in both of magnetic field and metrics.
         !  Difference between the flux-surface averaged magnetic field <omg>
         !  and the value at the magnetic axis B_ax also appears.
         !  cb = (psi_p(r))'/(cx*cy) = B_ax
         !  Normalized omg = B(z)/B_ax and cb = 1 in the B_ax unit.
         !-
          r_major = 1._DP       ! in the R0 unit
          r_0 = r_major * eps_r ! Minor radius of flux-tube center
          cx = 1._DP
          cy = r_0/q_0
          cb = 1._DP

          theta = 2._DP*atan( sqrt( (1._DP+eps_r)/(1._DP-eps_r) ) &
                             * tan(wzz/2._DP) )
          q_bar = dsqrt( 1._DP - eps_r**2 )*q_0
          gomg = sqrt( q_bar**2 + eps_r**2 ) &
               / ( 1._DP + eps_r*cos( theta ) ) / q_bar
         !- Metrics in GKV coordinates (x,y,z)
          gdomgdz = eps_r * sin(theta) * sqrt( q_bar**2 + eps_r**2 ) &
                     / ( 1._DP + eps_r * cos( theta ) )**2           &
                     / ( 1._DP - eps_r * cos( wzz) ) / q_0
          gdomgdx = -( cos(theta)                                           &
                       - eps_r*(1._DP-s_hat+eps_r**2*q_0**2/q_bar**2)        &
                              *(1._DP+eps_r*cos(theta))/(q_bar**2+eps_r**2) &
                       - eps_r*sin(theta)**2/(1._DP-eps_r**2)               &
                       ) / ((1._DP + eps_r*cos(theta))**2)                  &
                         * sqrt(q_bar**2+eps_r**2) / q_bar / r_major
          gdomgdy = 0._DP
          ggxx = (q_0/q_bar)**2
          ggxy = ( s_hat*wzz*q_0/q_bar - eps_r*sin(wzz)/(1._DP-eps_r**2) )*q_0/q_bar
          ggxz = - sin(wzz)/(1._DP-eps_r**2)/r_major*q_0/q_bar
          ggyy = (s_hat*wzz*q_0/q_bar)**2                                      &
                 - 2._DP*q_0/q_bar*s_hat*wzz*eps_r*sin(wzz)/(1._DP-eps_r**2)   &
                 + (q_bar**2+eps_r**2)/((1._DP+eps_r*cos(theta))**2)/(q_0**2) &
                 + (eps_r*sin(wzz))**2/(1._DP-eps_r**2)**2
          ggyz = ( -s_hat*wzz*q_0/q_bar*sin(wzz)/(1._DP-eps_r**2)           &
                   + ((q_bar/q_0)**2)/((1._DP+eps_r*cos(theta))**2)/eps_r  &
                   + eps_r*(sin(wzz)**2)/((1._DP-eps_r**2)**2)              &
                 ) / r_major
          ggzz = ( ((q_bar/q_0)**2)/((1._DP+eps_r*cos(theta))**2)/(eps_r**2) &
                   + (sin(wzz)**2)/((1._DP-eps_r**2)**2)                      &
                 ) / (r_major**2)
          grootg_xyz = q_0*r_major*( 1._DP+eps_r*cos(theta) )**2


        else if( trim(equib_type) == "vmec" ) then
         !- VMEC-BoozXform interface for stellarator equilibirum -
         !  References on the previous implementation by VMEC-NEWBOZ is
         !  [Ref.1] M. Nunami, T.-H. Watanabe, H. Sugama, Plasma Fusion Res. 5,
         !          016 (2010).
         !  New interface for VMEC-BoozXform is developed by M. Nakata and 
         !  M. Nunami (Aug. 2016) in the same manner for IGS.
         !
         !  Consider flux coordinates (rho,theta,zeta).
         !  Using the toroidal flux psi_t, the normalized minor radius is
         !  rho= sqrt(psi_t/psi_ta), and the minor radius at the last closed
         !  flux surface is a=sqrt(2*psi_ta/B_ax).
         !  Poloidal and toroidal angles are defined in the Boozer coordinates.
         !  GKV coordinates (x,y,z) are
         !    x = cx*(rho-rho0)
         !    y = cy*(q(r)*theta-zeta)
         !    z = theta
         !  with cx=a, cy=cx*rho0/q0.
         !  In these definitions, the factor on the magnetic field 
         !  B = cb * \nabla x \times \nabla y is
         !  cb = (psi_p(rho))'/(cx*cy) = B_ax.
         !  Normalized omg = B(z)/B_ax and cb = 1 in the B_ax unit.
         !  The reference length is set to be r_major at the magnetic axis R0.
         !-
          r_major = 1._DP ! in the R0 unit
          q_bar = q_0
          isw = 1
          call vmecbzx_boozx_coeff( isw, nss, ntheta, nzeta, s_input, iz, wzz, lz_l, &  ! input 
                                    s_0, q_0, s_hat, eps_r, phi_ax,                  &  ! output
                                    gomg, grootg_xyz, gdomgdx, gdomgdz, gdomgdy,     &
                                    ggxx, ggxy, ggxz, ggyy, ggyz, ggzz )
         ! NOTE: phi_ax axisymmetric toroidal angle is stored for vmec, rather than theta
          theta = phi_ax
          cx = eps_r/s_0 ! = eps_a = a/R0
          cy = cx*s_0/q_0
          cb = 1._DP

        else if( trim(equib_type) == "eqdsk" ) then
         !- EQDSK-IGS interface for tokamak equilibirum -
         !  [Ref.] M. Nakata, A. Matsuyama, N. Aiba, S. Maeyama, M. Nunami,
         !         and T.-H. Watanabe, Plasma Fusion Res. 9, 1403029 (2014).
         !
         !  Consider flux coordinates (rho,theta,zeta).
         !  GKV coordinates (x,y,z) are
         !    x = cx*(rho-rho0)
         !    y = cy*(q(r)*theta-zeta)
         !    z = theta
         !  with cx=a, cy=cx*rho0/q0. All explanation is the same as that in
         !  equib_type == "vmec", except that poloidal and toroidal angles have
         !  a choice of freedom: Hamada, Boozer, or axisymmetric coordinates.
         !-
          r_major = 1._DP ! in the R0 unit
          q_bar = q_0
          isw = 1
          call igs_coeff( isw, mc_type, nss, ntheta, s_input, wzz, lz_l, &  ! input 
                          s_0, q_0, s_hat, eps_r, theta,                 &  ! output
                          gomg, grootg_xyz, gdomgdx, gdomgdz, gdomgdy,   &
                          ggxx, ggxy, ggxz, ggyy, ggyz, ggzz )
          cx = eps_r/s_0 ! = eps_a = a/R0
          cy = cx*s_0/q_0
          cb = 1._DP

!sakano_ring-dipole st 202303
        else if( trim(equib_type) == "ring" ) then
         !- Ring dipole geometry -
         !  [Ref.] J. Sakano, Master thesis, Nagoya University (in Japanese).
         !
         !  Consider flux coordinates (Psi,Theta,phi), where the magnetic
         !  poloidal flux Psi<0, the geometrical poloidal angle Theta = arctan(Z/(R-a)), 
         !  the azimuthal angle of the cylindrical coordinate phi.
         !  There is a ring current in direction of phi at R=a. The field line 
         !  passing through (R,Z)=(R0,0) is picked up as a flux-tube domain.
         !  
         !  GKV coordinates (x,y,z) are (right-handed system)
         !    x = cx*(Psi0 - Psi)/Psi0
         !    y = cy*phi
         !    z = Theta
         !  with cx=R0, cy=R0. Note that Psi0 is the magnetic poloidal flux 
         !  at the center of the considered flux-tube domain.
         !  In these definitions, the factor on the magnetic field 
         !  B = cb * \nabla x \times \nabla y is cb = Psi0/(R0*R0) = B0,
         !  where B0 is the magnetic field strength at (R,Z)=(R0,0).
         !  Normalized omg = B(z)/B0 and cb = 1 in the B0 unit.
         !  The reference length is set to be R0 (not the ring current at R=a).
         !  The normalized parameter to specify the flux-tube is 
         !    ring_a = a / R0
         !-
          r_major = 1._DP ! in the R0 unit
          q_bar   = 0._DP
          theta = wzz
          call ring_coordinates( ring_a, wzz, &                                ! input
                                 gomg, ub_dot_grdb, ub_crs_grdb, ggxx, ggxy, & ! output
                                 ggxz, ggyy, ggyz, ggzz, grootg_xyz, gdomgdx, gdomgdz )
          gdomgdy = 0._DP
          cx = 1._DP
          cy = 1._DP
          cb = 1._DP
!sakano_ring-dipole end 202303

        else

          write( olog, * ) " # wrong choice of the equilibrium "
          call flush(olog)
          call MPI_Finalize(ierr_mpi)
          stop

        end if

        call mtr_global%init(iz, wzz, theta, gomg,                  &
                             gdomgdx, gdomgdy, gdomgdz, ggxx, ggxy, &
                             ggxz, ggyy, ggyz, ggzz, grootg_xyz,    &
                             gdomgdr, gdomgdt, gdomgdq, ggrr, ggrt, &
                             ggrq, ggtt, ggtq, ggqq, grootg_rtq)

      end do   ! iz loop ends

      call mtr_global%xyz2rtq
      call mtr_fourier%init
      call mtr_fourier%dft_rtq2coef(mtr_global)

      call mtr_local%copy_global(mtr_global)

      if ( rankg == 0 ) then
        do iz = -global_nz, global_nz-1
          write( omtr, fmt="(f15.8,SP,256E24.14e3)") &
            mtr_global%zz(iz), mtr_global%theta(iz),      &
            mtr_global%omg(iz), mtr_global%domgdx(iz),    &
            mtr_global%domgdy(iz), mtr_global%domgdz(iz), &
            mtr_global%gxx(iz), mtr_global%gxy(iz),       &
            mtr_global%gxz(iz), mtr_global%gyy(iz),       &
            mtr_global%gyz(iz), mtr_global%gzz(iz),       &
            mtr_global%rootg_xyz(iz)
        end do
        !call flush(omtr)
        do iz = -global_nz, global_nz-1
          write( omtf, fmt="(f15.8,SP,256E24.14e3)") &
            mtr_global%zz(iz), mtr_global%theta(iz),      &
            mtr_global%omg(iz), mtr_global%domgdr(iz),    &
            mtr_global%domgdt(iz), mtr_global%domgdq(iz), &
            mtr_global%grr(iz), mtr_global%grt(iz),       &
            mtr_global%grq(iz), mtr_global%gtt(iz),       &
            mtr_global%gtq(iz), mtr_global%gqq(iz),       &
            mtr_global%rootg_rtq(iz)
        end do
        !call flush(omtr)
      end if

     !%%% For debug %%%
     ! do iz = -nz, nz-1
     !   write( 990000000+rankg, fmt="(f15.8,SP,256E24.14e3)") &
     !     mtr_local%zz(iz), mtr_local%theta(iz),      &
     !     mtr_local%omg(iz), mtr_local%domgdx(iz),    &
     !     mtr_local%domgdy(iz), mtr_local%domgdz(iz), &
     !     mtr_local%gxx(iz), mtr_local%gxy(iz),       &
     !     mtr_local%gxz(iz), mtr_local%gyy(iz),       &
     !     mtr_local%gyz(iz), mtr_local%gzz(iz),       &
     !     mtr_local%rootg_xyz(iz)
     ! end do
     ! call mtr_global%rtq2xyz
     ! call mtr_global%xyz2rtq
     ! call mtr_local%init(mtr_fourier, time_shearflow=0._DP)
     ! if ( rankg == 0 ) then
     !   do iz = -global_nz, global_nz-1
     !     write( 900000011, fmt="(f15.8,SP,256E24.14e3)") &
     !       mtr_global%zz(iz), mtr_global%theta(iz),      &
     !       mtr_global%omg(iz), mtr_global%domgdx(iz),    &
     !       mtr_global%domgdy(iz), mtr_global%domgdz(iz), &
     !       mtr_global%gxx(iz), mtr_global%gxy(iz),       &
     !       mtr_global%gxz(iz), mtr_global%gyy(iz),       &
     !       mtr_global%gyz(iz), mtr_global%gzz(iz),       &
     !       mtr_global%rootg_xyz(iz)
     !   end do
     !   do iz = -global_nz, global_nz-1
     !     write( 900000012, fmt="(f15.8,SP,256E24.14e3)") &
     !       mtr_global%zz(iz), mtr_global%theta(iz),      &
     !       mtr_global%omg(iz), mtr_global%domgdr(iz),    &
     !       mtr_global%domgdt(iz), mtr_global%domgdq(iz), &
     !       mtr_global%grr(iz), mtr_global%grt(iz),       &
     !       mtr_global%grq(iz), mtr_global%gtt(iz),       &
     !       mtr_global%gtq(iz), mtr_global%gqq(iz),       &
     !       mtr_global%rootg_rtq(iz)
     !   end do
     ! end if
     ! do iz = -nz, nz-1
     !   write( 980000000+rankg, fmt="(f15.8,SP,256E24.14e3)") &
     !     mtr_local%zz(iz), mtr_local%theta(iz),      &
     !     mtr_local%omg(iz), mtr_local%domgdx(iz),    &
     !     mtr_local%domgdy(iz), mtr_local%domgdz(iz), &
     !     mtr_local%gxx(iz), mtr_local%gxy(iz),       &
     !     mtr_local%gxz(iz), mtr_local%gyy(iz),       &
     !     mtr_local%gyz(iz), mtr_local%gzz(iz),       &
     !     mtr_local%rootg_xyz(iz)
     ! end do
     ! do iz = -nz, nz-1
     !   write( 970000000+rankg, fmt="(f15.8,SP,256E24.14e3)") &
     !     mtr_local%zz_labframe(iz), mtr_local%theta(iz),      &
     !     mtr_local%omg(iz), mtr_local%domgdr(iz),    &
     !     mtr_local%domgdt(iz), mtr_local%domgdq(iz), &
     !     mtr_local%grr(iz), mtr_local%grt(iz),       &
     !     mtr_local%grq(iz), mtr_local%gtt(iz),       &
     !     mtr_local%gtq(iz), mtr_local%gqq(iz),       &
     !     mtr_local%rootg_rtq(iz)
     ! end do
     !%%%%%%%%%%%%%%%%%%

  END SUBROUTINE geom_init_metric
 

!--------------------------------------
  SUBROUTINE geom_set_operators
!--------------------------------------
    implicit none
    real(kind=DP) :: wzz    ! The rotating flux tube coordinate (= z'')
    real(kind=DP) :: zz_lab ! The flux-coordinate theta in the lab frame (= z''+t*gamma_e/s_hat)
    real(kind=DP) :: kkx, kky, domgdz, domgdx, domgdy
    real(kind=DP) :: bb, kmo
    real(kind=DP) :: gg0

    real(kind=DP) :: cfsrf_l
    real(kind=DP), dimension(1:3,1:3) :: gg
    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: wf
    complex(kind=DP), dimension(:,:,:), allocatable :: nw
    real(kind=DP), dimension(:,:,:), allocatable :: ww

    integer :: mx, my, iz, iv, im, is

        do iz = -nz, nz-1

          wzz = zz(iz)
          zz_lab = mtr_local%zz_labframe(iz)
          omg(iz)   = mtr_local%omg(iz)
          domgdx    = mtr_local%domgdx(iz)
          domgdy    = mtr_local%domgdy(iz)
          domgdz    = mtr_local%domgdz(iz)
          gg(1,1)   = mtr_local%gxx(iz)
          gg(1,2)   = mtr_local%gxy(iz)
          gg(1,3)   = mtr_local%gxz(iz)
          gg(2,2)   = mtr_local%gyy(iz)
          gg(2,3)   = mtr_local%gyz(iz)
          gg(3,3)   = mtr_local%gzz(iz)
          rootg(iz) = mtr_local%rootg_xyz(iz)
          gg(2,1)   = gg(1,2)
          gg(3,1)   = gg(1,3)
          gg(3,2)   = gg(2,3)

!!! for slab model
          if ( trim(equib_type) == "slab") then

            dpara(iz) = dz * q_0 * r_major

            do im = 0, nm
              vp(iz,im)  = sqrt( 2._DP * mu(im) )!* omg(iz) )
              mir(iz,im) = 0._DP
              do iv = 1, 2*nv
                vdx(iz,iv,im) = 0._DP
                vdy(iz,iv,im) = 0._DP
                vsy(iz,iv,im) =                                           &
                  - sgn(ranks) * tau(ranks) / Znum(ranks)                 & 
                  * ( R0_Ln(ranks) + R0_Lt(ranks) * ( 0.5_DP*vl(iv)**2    &
                                           + omg(iz)*mu(im) - 1.5_DP ) )
              end do
            end do   ! im loop ends

            ksq(:,:,iz) = 0._DP
            do my = ist_y, iend_y
              do mx = -nx, nx
                ksq(mx,my,iz) = kx(mx)**2 + ky(my)**2
              end do
            end do

!!! for the concentric and large-aspect-ratio model !!!
          else if( trim(equib_type) == "analytic" ) then

            dpara(iz) = dz * q_0 * r_major

            do im = 0, nm

              vp(iz,im)  = sqrt( 2._DP * mu(im) * omg(iz) )
              mir(iz,im) = mu(im) * eps_r / ( q_0 * r_major )  &
                       * ( sin(zz_lab)                                          &
                         + eps_hor * lmmq   * sin( lmmq   * zz_lab - malpha )   &
                         + eps_mor * lmmqm1 * sin( lmmqm1 * zz_lab - malpha )   &
                         + eps_por * lmmqp1 * sin( lmmqp1 * zz_lab - malpha ) )

              do iv = 1, 2*nv
                    vdx(iz,iv,im)=                                            &
                      - ( vl(iv)**2 + omg(iz)*mu(im) ) * eps_rnew / r_major         &
                      * ( 0._DP * ( rdeps00 + rdeps1_0 * cos( zz_lab )             &
                             + rdeps2_10 * cos( lmmq   * zz_lab - malpha )          &
                             + rdeps1_10 * cos( lmmqm1 * zz_lab - malpha )          &
                             + rdeps3_10 * cos( lmmqp1 * zz_lab - malpha ) )        &
                             + ( 1._DP + s_hat * wzz * 0._DP )                 &
                             * ( sin( zz_lab )                                      &
                             + eps_hor * lprd   * sin( lmmq   * zz_lab - malpha )   &
                             + eps_mor * lprdm1 * sin( lmmqm1 * zz_lab - malpha )   &
                             + eps_por * lprdp1 * sin( lmmqp1 * zz_lab - malpha ) ) &
                         ) * ( sgn(ranks) * tau(ranks) / Znum(ranks) )

                    vdy(iz,iv,im)=                                            &
                      - ( vl(iv)**2 + omg(iz)*mu(im) ) * eps_rnew / r_major         &
                      * ( 1._DP * ( rdeps00 + rdeps1_0 * cos( zz_lab )             &
                             + rdeps2_10 * cos( lmmq   * zz_lab - malpha )          &
                             + rdeps1_10 * cos( lmmqm1 * zz_lab - malpha )          &
                             + rdeps3_10 * cos( lmmqp1 * zz_lab - malpha ) )        &
                             + ( 0._DP + s_hat * wzz * 1._DP )                 &
                             * ( sin( zz_lab )                                      &
                             + eps_hor * lprd   * sin( lmmq   * zz_lab - malpha )   &
                             + eps_mor * lprdm1 * sin( lmmqm1 * zz_lab - malpha )   &
                             + eps_por * lprdp1 * sin( lmmqp1 * zz_lab - malpha ) ) &
                         ) * ( sgn(ranks) * tau(ranks) / Znum(ranks) )

                    vsy(iz,iv,im) =                                           &
                      - sgn(ranks) * tau(ranks) / Znum(ranks)                 & 
                      * ( R0_Ln(ranks) + R0_Lt(ranks) * ( 0.5_DP*vl(iv)**2    &
                                               + omg(iz)*mu(im) - 1.5_DP ) )

              end do

            end do   ! im loop ends

            ksq(:,:,iz) = 0._DP
            do my = ist_y, iend_y
              do mx = -nx, nx
                ksq(mx,my,iz) = ( kx(mx) + s_hat * wzz * ky(my) )**2 + ky(my)**2
              end do
            end do

!!! for s-alpha !!! <--- the current version is the same as "analytic"
          else if( trim(equib_type) == "s-alpha" .or. trim(equib_type) == "s-alpha-shift" ) then

            dpara(iz) = dz* q_0 * r_major

            kkx = -r_major * (q_0/q_bar) &
                           * ( gg(1,1)*gg(2,3) - gg(1,2)*gg(1,3) )*domgdz
            kky =  r_major * (q_bar/q_0) &
                           * ( domgdx - ( gg(1,2)*gg(2,3) - gg(2,2)*gg(1,3) )*domgdz )

            do im = 0, nm

              vp(iz,im)  = sqrt( 2._DP * mu(im) * omg(iz) )

              mir(iz,im) = mu(im) * (q_0/q_bar) * domgdz / ( omg(iz)*rootg(iz) )

              do iv = 1, 2*nv
                    vdx(iz,iv,im) =                                     &
                       ( vl(iv)**2 + omg(iz)*mu(im) ) / r_major         &   
                      * kkx                                             &
                      * ( sgn(ranks) * tau(ranks) / Znum(ranks) )

                    vdy(iz,iv,im) =                                     &
                       ( vl(iv)**2 + omg(iz)*mu(im) ) / r_major         &   
                      * kky                                             &
                      * ( sgn(ranks) * tau(ranks) / Znum(ranks) )

                    vsy(iz,iv,im) =                                           &
                      - sgn(ranks) * tau(ranks) / Znum(ranks)                 &
                      * ( R0_Ln(ranks) + R0_Lt(ranks) * ( 0.5_DP*vl(iv)**2    &
                                               + omg(iz)*mu(im) - 1.5_DP ) )  &
                      * (q_bar/q_0)
              end do

            end do   ! im loop ends


            ksq(:,:,iz) = 0._DP
            do my = ist_y, iend_y
              do mx = -nx, nx
                ksq(mx,my,iz) = ( kx(mx) + ( s_hat * wzz - alpha_MHD*sin(zz_lab) ) &
                                * ky(my) )**2 + ky(my)**2 ! with Shafranov shift
              end do
            end do

!!! for circular MHD equilibrium !!!
          else if( trim(equib_type) == "circ-MHD" ) then

            dpara(iz) = dz * omg(iz) * rootg(iz)

            kkx =  r_major*( -domgdy + (gg(1,3)*gg(1,2) - gg(1,1)*gg(2,3))*domgdz/omg(iz)**2 )
            kky =  r_major*(  domgdx + (gg(1,3)*gg(2,2) - gg(1,2)*gg(2,3))*domgdz/omg(iz)**2 )

            do im = 0, nm

              vp(iz,im)  = sqrt( 2._DP * mu(im) * omg(iz) )

              mir(iz,im) = mu(im) * domgdz / ( omg(iz)*rootg(iz) )

              do iv = 1, 2*nv
                    vdx(iz,iv,im)=                                              &
                       ( vl(iv)**2 + omg(iz)*mu(im) ) / ( r_major*omg(iz) )     &
                      * kkx                                                     &
                      * ( sgn(ranks) * tau(ranks) / Znum(ranks) )

                    vdy(iz,iv,im)=                                              &
                       ( vl(iv)**2 + omg(iz)*mu(im) ) / ( r_major*omg(iz) )     &
                      * kky                                                     &
                      * ( sgn(ranks) * tau(ranks) / Znum(ranks) )

                    vsy(iz,iv,im) =                                           &
                      - sgn(ranks) * tau(ranks) / Znum(ranks)                 &
                      * ( R0_Ln(ranks) + R0_Lt(ranks) * ( 0.5_DP*vl(iv)**2    &
                                               + omg(iz)*mu(im) - 1.5_DP ) )  
              end do

            end do   ! im loop ends

            ksq(:,:,iz) = 0._DP
            do my = ist_y, iend_y
              do mx = -nx, nx
                ksq(mx,my,iz) = (kx(mx)**2)*gg(1,1)         &
                              + 2._DP*kx(mx)*ky(my)*gg(1,2) &
                              + (ky(my)**2)*gg(2,2)
              end do
            end do

!  this is new vmec-BoozXform interface  by M. Nakata & M. Nunami  (Aug. 2016)
          else if( trim(equib_type) == "vmec" ) then

            dpara(iz) = dz * omg(iz) * rootg(iz)

            kkx =  r_major*( -domgdy + (gg(1,3)*gg(1,2) - gg(1,1)*gg(2,3))*domgdz/omg(iz)**2 )
            kky =  r_major*(  domgdx + (gg(1,3)*gg(2,2) - gg(1,2)*gg(2,3))*domgdz/omg(iz)**2 )

            do im = 0, nm

              vp(iz,im)  = sqrt( 2._DP * mu(im) * omg(iz) )

              mir(iz,im) = mu(im) * domgdz / ( omg(iz)*rootg(iz) )

              do iv = 1, 2*nv
                    vdx(iz,iv,im) =                                              &
                       ( vl(iv)**2 + omg(iz)*mu(im) ) / ( r_major*omg(iz) )      &
                      * kkx                                                      &
                      * ( sgn(ranks) * tau(ranks) / Znum(ranks) )


                    vdy(iz,iv,im) =                                              &
                       ( vl(iv)**2 + omg(iz)*mu(im) ) / ( r_major*omg(iz) )      &
                      * kky                                                      &
                      * ( sgn(ranks) * tau(ranks) / Znum(ranks) )                &
                     - real(ibprime,kind=DP) * vl(iv)**2 / r_major / omg(iz)**2  & ! grad-p (beta-prime) term 
                      * ( beta*(R0_Ln(ranks) + R0_Lt(ranks)) )                   &
                      * ( sgn(ranks) * tau(ranks) / Znum(ranks) )

                    vsy(iz,iv,im) =                                           &
                      - sgn(ranks) * tau(ranks) / Znum(ranks)                 &
                      * ( R0_Ln(ranks) + R0_Lt(ranks) * ( 0.5_DP*vl(iv)**2    &
                                               + omg(iz)*mu(im) - 1.5_DP ) )  
              end do

            end do   ! im loop ends

            ksq(:,:,iz) = 0._DP
            do my = ist_y, iend_y
              do mx = -nx, nx
                ksq(mx,my,iz) = (kx(mx)**2)*gg(1,1)         &
                              + 2._DP*kx(mx)*ky(my)*gg(1,2) &
                              + (ky(my)**2)*gg(2,2)
              end do
            end do

          else if( trim(equib_type) == "eqdsk" ) then

            dpara(iz) = dz * omg(iz) * rootg(iz)

            kkx =  r_major*( -domgdy + (gg(1,3)*gg(1,2) - gg(1,1)*gg(2,3))*domgdz/omg(iz)**2 )
            kky =  r_major*(  domgdx + (gg(1,3)*gg(2,2) - gg(1,2)*gg(2,3))*domgdz/omg(iz)**2 )

            do im = 0, nm

              vp(iz,im)  = sqrt( 2._DP * mu(im) * omg(iz) )

              mir(iz,im) = mu(im) * domgdz / ( omg(iz)*rootg(iz) )

              do iv = 1, 2*nv
                    vdx(iz,iv,im) =                                              &
                       ( vl(iv)**2 + omg(iz)*mu(im) ) / ( r_major*omg(iz) )      &
                      * kkx                                                      &
                      * ( sgn(ranks) * tau(ranks) / Znum(ranks) )

                    vdy(iz,iv,im) =                                              &
                       ( vl(iv)**2 + omg(iz)*mu(im) ) / ( r_major*omg(iz) )      &
                      * kky                                                      &
                      * ( sgn(ranks) * tau(ranks) / Znum(ranks) )                &
                     - real(ibprime,kind=DP) * vl(iv)**2 / r_major / omg(iz)**2  & ! grad-p (beta-prime) term 
                      * ( beta*(R0_Ln(ranks) + R0_Lt(ranks)) )                   &
                      * ( sgn(ranks) * tau(ranks) / Znum(ranks) )

                    vsy(iz,iv,im) =                                           &
                      - sgn(ranks) * tau(ranks) / Znum(ranks)                 &
                      * ( R0_Ln(ranks) + R0_Lt(ranks) * ( 0.5_DP*vl(iv)**2    &
                                               + omg(iz)*mu(im) - 1.5_DP ) )  
              end do

            end do   ! im loop ends

            ksq(:,:,iz) = 0._DP
            do my = ist_y, iend_y
              do mx = -nx, nx
                ksq(mx,my,iz) = (kx(mx)**2)*gg(1,1)         &
                              + 2._DP*kx(mx)*ky(my)*gg(1,2) &
                              + (ky(my)**2)*gg(2,2)
              end do
            end do

!sakano_ring-dipole st 202303
          else if( trim(equib_type) == "ring" ) then

            dpara(iz) = dz * omg(iz) * rootg(iz)

            kkx = 0._DP
            kky = r_major*( domgdx + (gg(1,3)*gg(2,2) - gg(1,2)*gg(2,3))*domgdz/omg(iz)**2 )

            do im = 0, nm
! r_major = 1 is assumed as the equilibrium length unit
! B on the equatorial plane is also unity

              vp(iz,im)  = sqrt( 2._DP * mu(im) * omg(iz) )

              !mir(iz,im) = mu(im) * ub_dot_grdb
              mir(iz,im) = mu(im) * domgdz / ( omg(iz)*rootg(iz) )

              do iv = 1, 2*nv
                vdx(iz,iv,im) = 0._DP
                   
                !vdy(iz,iv,im) =                                        &
                !        ( vl(iv)**2 + omg(iz)*mu(im) )                 &
                !      * ( ub_crs_grdb / omg(iz)**2 ) * sqrt( gg(2,2) ) &
                !      * ( sgn(ranks) * tau(ranks) / Znum(ranks) )      ! ion's vdy is negative y direction
                vdy(iz,iv,im) =                                              &
                        ( vl(iv)**2 + omg(iz)*mu(im) ) / ( r_major*omg(iz) ) &
                      * kky                                                  &
                      * ( sgn(ranks) * tau(ranks) / Znum(ranks) ) 
                vsy(iz,iv,im)=                                             &
                      - sgn(ranks) * tau(ranks) / Znum(ranks)              &
                      * ( R0_Ln(ranks) + R0_Lt(ranks) * ( 0.5_DP*vl(iv)**2 &
                                               + omg(iz)*mu(im) - 1.5_DP ) ) ! ion's vsy is negative y directuin
              end do

            end do   ! im loop ends

            ksq(:,:,iz) = 0._DP
            do my = ist_y, iend_y
              do mx = -nx, nx
                ksq(mx,my,iz) = ( kx(mx) * gg(1,1) )**2 + ( ky(my) * gg(2,2) )**2
              end do
            end do
!sakano_ring-dipole end 202303

          else

            write( olog, * ) " # wrong choice of the equilibrium "
            call flush(olog)
            call MPI_Finalize(ierr_mpi)
            stop

          end if


          do im = 0, nm
            do my = ist_y, iend_y
              do mx = -nx, nx
                kmo           = sqrt( 2._DP * ksq(mx,my,iz) * mu(im) / omg(iz) ) &
                               * dsqrt( tau(ranks)*Anum(ranks) ) / Znum(ranks)
                call math_j0( kmo, j0(mx,my,iz,im) )
                call math_j1( kmo, j1(mx,my,iz,im) )
                call math_j2( kmo, j2(mx,my,iz,im) )
              end do
            end do
          end do


          do my = ist_y, iend_y
            do mx = -nx, nx
              bb     = ksq(mx,my,iz) / omg(iz)**2 &
                        * tau(ranks)*Anum(ranks)/(Znum(ranks)**2)
              call math_g0( bb, g0(mx,my,iz) )
            end do
          end do

        end do   ! iz loop ends

        cfsrf   = 0._DP
        cfsrf_l = 0._DP
        do iz = -nz, nz-1
          cfsrf_l   = cfsrf_l + rootg(iz)
                                            ! normalization coefficient for 
                                            ! the surface average
        end do
        call MPI_Allreduce( cfsrf_l, cfsrf, 1, MPI_DOUBLE_PRECISION, &
                            MPI_SUM, zsp_comm_world, ierr_mpi )

        if ( vel_rank == 0 ) then
          do iz = -nz, nz-1
            dvp(iz)  = sqrt( 2._DP * (0.5_DP * dm**2) * omg(iz) )
          end do
        end if
        call MPI_Bcast( dvp, 2*nz, MPI_DOUBLE_PRECISION, 0, &
                        vel_comm_world, ierr_mpi )

        do im = 0, nm
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              fmx(iz,iv,im)   = exp( - 0.5_DP * vl(iv)**2 - omg(iz) * mu(im) ) &
                              / sqrt( twopi**3 )
            end do
          end do
        end do

        allocate( ww(-nx:nx,0:ny,-nz:nz-1) )

! --- GK polarization factor for efield calculation 
        fct_poisson(:,:,:) = 0._DP
        fct_e_energy(:,:,:) = 0._DP

        ww(:,:,:) = 0._DP
        do iz = -nz, nz-1
          do my = ist_y, iend_y
            do mx = -nx, nx

              if ( rankw == 0 .and. mx == 0 .and. my == 0 ) then !- (0,0) mode

                fct_poisson(mx,my,iz) = 0._DP
                fct_e_energy(mx,my,iz) = 0._DP

              else

                ww(mx,my,iz) = lambda_i * ksq(mx,my,iz)
                do is = 0, ns-1
                  bb   = ksq(mx,my,iz) / omg(iz)**2 &
                          * tau(is)*Anum(is)/(Znum(is)**2)
                  call math_g0( bb, gg0 )
                  ww(mx,my,iz) = ww(mx,my,iz)  &
                               + Znum(is) * fcs(is) / tau(is) * ( 1._DP - gg0 )
                end do
                fct_poisson(mx,my,iz) = 1._DP / ww(mx,my,iz)
                fct_e_energy(mx,my,iz) = ww(mx,my,iz)

              end if

            end do
          end do
        end do


! --- ZF-factor for adiabatic model
        if ( ns == 1 ) then

          ww(:,:,:) = 0._DP
          do iz = -nz, nz-1
            my = 0
              do mx = -nx, nx
                ww(mx,my,iz) = ( 1._DP - g0(mx,my,iz) )       &
                             / ( 1._DP - g0(mx,my,iz) + tau(0)*tau_ad )
              end do
          end do

          call intgrl_fsrf ( ww, fctgt )

          if ( rankw == 0 )  then
            fctgt(0)   = ( 1._DP - g0(0,0,0) ) / ( 1._DP - g0(0,0,0) + tau(0)*tau_ad )
                                              ! g0(0,0,iz) has no z dependence
          endif

        endif

        deallocate( ww )

        allocate( wf(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) )
        allocate( nw(-nx:nx,0:ny,-nz:nz-1) )
        wf(:,:,:,:,:) = ( 0._DP, 0._DP )
        nw(:,:,:) = ( 0._DP, 0._DP )

! --- GK polarization factor for mfield calculation 
        fct_ampere(:,:,:) = 0._DP
        fct_m_energy(:,:,:) = 0._DP

        if ( beta .ne. 0._DP ) then
       
          do im = 0, nm
            do iv = 1, 2*nv
              do iz = -nz, nz-1
                do my = ist_y, iend_y
                  do mx = -nx, nx
                    wf(mx,my,iz,iv,im) = Znum(ranks) * fcs(ranks) / Anum(ranks)  &
                                       * vl(iv)**2 * j0(mx,my,iz,im)**2 * fmx(iz,iv,im)
                  end do
                end do
              end do
            end do
          end do
  
          call intgrl_v0_moment_ms ( wf, nw )
  
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                fct_ampere(mx,my,iz) = 1._DP / real( ksq(mx,my,iz) + beta * nw(mx,my,iz), kind=DP )
                fct_m_energy(mx,my,iz) = ksq(mx,my,iz) / beta
              end do
            end do
          end do
  
          if ( rankw == 0 ) then
            do iz = -nz, nz-1
              fct_ampere(0,0,iz) = 0._DP
              fct_m_energy(0,0,iz) = 0._DP
            end do
          end if

        end if

        deallocate( wf )
        deallocate( nw )

  END SUBROUTINE geom_set_operators

!--------------------------------------
  SUBROUTINE geom_reset_time(time_shearflow)
!--------------------------------------
    implicit none
    real(kind=DP), intent(in) :: time_shearflow
      call mtr_local%init(mtr_fourier, time_shearflow)
      call geom_set_operators
      !NOTE: colliimp_set_param in GKV_colliimp should also be updated.
  END SUBROUTINE geom_reset_time

!--------------------------------------
  SUBROUTINE geom_increment_time(dt_shearflow)
!--------------------------------------
    implicit none
    real(kind=DP), intent(in) :: dt_shearflow
      call mtr_local%update(mtr_fourier, dt_shearflow)
      call geom_set_operators
      !NOTE: colliimp_set_param in GKV_colliimp should also be updated.
  END SUBROUTINE geom_increment_time


!--------------------------------------
  SUBROUTINE metric_global_init(self, iz, wzz, theta, gomg,            &
                                gdomgdx, gdomgdy, gdomgdz, ggxx, ggxy, &
                                ggxz, ggyy, ggyz, ggzz, grootg_xyz,    &
                                gdomgdr, gdomgdt, gdomgdq, ggrr, ggrt, &
                                ggrq, ggtt, ggtq, ggqq, grootg_rtq)
!--------------------------------------
    implicit none
    class(metric_global), intent(inout) :: self
    integer, intent(in) :: iz
    real(kind=DP), intent(in) :: wzz, theta, gomg
    real(kind=DP), intent(in) :: gdomgdx, gdomgdy, gdomgdz, &
                     ggxx, ggxy, ggxz, ggyy, ggyz, ggzz, grootg_xyz
    real(kind=DP), intent(in) :: gdomgdr, gdomgdt, gdomgdq, &
                     ggrr, ggrt, ggrq, ggtt, ggtq, ggqq, grootg_rtq

      self%zz(iz)     = wzz
      self%theta(iz)  = theta
      self%omg(iz)    = gomg
      self%domgdx(iz) = gdomgdx
      self%domgdy(iz) = gdomgdy
      self%domgdz(iz) = gdomgdz
      self%gxx(iz) = ggxx
      self%gxy(iz) = ggxy
      self%gxz(iz) = ggxz
      self%gyy(iz) = ggyy
      self%gyz(iz) = ggyz
      self%gzz(iz) = ggzz
      self%rootg_xyz(iz) = grootg_xyz
      self%domgdr(iz) = gdomgdr
      self%domgdt(iz) = gdomgdt
      self%domgdq(iz) = gdomgdq
      self%grr(iz) = ggrr
      self%grt(iz) = ggrt
      self%grq(iz) = ggrq
      self%gtt(iz) = ggtt
      self%gtq(iz) = ggtq
      self%gqq(iz) = ggqq
      self%rootg_rtq(iz) = grootg_rtq

  END SUBROUTINE metric_global_init


!--------------------------------------
  SUBROUTINE metric_global_xyz2rtq(self)
!--------------------------------------
    implicit none
    class(metric_global), intent(inout) :: self
    real(kind=DP) :: wzz
    real(kind=DP) :: gdomgdx, gdomgdy, gdomgdz, &
                     ggxx, ggxy, ggxz, ggyy, ggyz, ggzz, grootg_xyz
    real(kind=DP) :: gdomgdr, gdomgdt, gdomgdq, &
                     ggrr, ggrt, ggrq, ggtt, ggtq, ggqq, grootg_rtq
    integer :: iz

      do iz = -global_nz, global_nz-1
       ! load (x,y,z)
        wzz = self%zz(iz)
        gdomgdx = self%domgdx(iz)
        gdomgdy = self%domgdy(iz)
        gdomgdz = self%domgdz(iz)
        ggxx = self%gxx(iz)
        ggxy = self%gxy(iz)
        ggxz = self%gxz(iz)
        ggyy = self%gyy(iz)
        ggyz = self%gyz(iz)
        ggzz = self%gzz(iz)
        grootg_xyz = self%rootg_xyz(iz)

       ! translate (x,y,z)->(r,t,q)=(rho,theta,zeta)
       !   NOTE: cx*rho0/(cy*q_0=1) is used.
        gdomgdr = cx*gdomgdx + cx*s_hat*wzz*gdomgdy
        gdomgdt = gdomgdz + cy*q_0*gdomgdy
        gdomgdq = - cy*gdomgdy
        ggrr = ggxx/cx**2
        ggrt = ggxz/cx
        ggrq = (s_hat*wzz*ggxx-ggxy)/(cx*cy) + q_0*ggxz/cx
        ggtt = ggzz
        ggtq = (s_hat*wzz*ggxz-ggyz)/cy + q_0*ggzz
        ggqq = (s_hat*wzz/cy)**2*ggxx - 2._DP*(s_hat*wzz/cy**2)*ggxy &
             + 2._DP*(q_0*s_hat*wzz/cy)*ggxz + ggyy/cy**2            &
             - 2._DP*(q_0/cy)*ggyz + q_0**2*ggzz
        grootg_rtq = cx*cy*grootg_xyz

       ! store (r,t,q)
        self%domgdr(iz) = gdomgdr
        self%domgdt(iz) = gdomgdt
        self%domgdq(iz) = gdomgdq
        self%grr(iz) = ggrr
        self%grt(iz) = ggrt
        self%grq(iz) = ggrq
        self%gtt(iz) = ggtt
        self%gtq(iz) = ggtq
        self%gqq(iz) = ggqq
        self%rootg_rtq(iz) = grootg_rtq
      end do

  END SUBROUTINE metric_global_xyz2rtq


!--------------------------------------
  SUBROUTINE metric_global_rtq2xyz(self)
!--------------------------------------
    implicit none
    class(metric_global), intent(inout) :: self
    real(kind=DP) :: wzz
    real(kind=DP) :: gdomgdx, gdomgdy, gdomgdz, &
                     ggxx, ggxy, ggxz, ggyy, ggyz, ggzz, grootg_xyz
    real(kind=DP) :: gdomgdr, gdomgdt, gdomgdq, &
                     ggrr, ggrt, ggrq, ggtt, ggtq, ggqq, grootg_rtq
    integer :: iz

      do iz = -global_nz, global_nz-1
       ! load (r,t,q)=(rho,theta,zeta)
        wzz = self%zz(iz)
        gdomgdr = self%domgdr(iz)
        gdomgdt = self%domgdt(iz)
        gdomgdq = self%domgdq(iz)
        ggrr = self%grr(iz)
        ggrt = self%grt(iz)
        ggrq = self%grq(iz)
        ggtt = self%gtt(iz)
        ggtq = self%gtq(iz)
        ggqq = self%gqq(iz)
        grootg_rtq = self%rootg_rtq(iz)

       ! translate (r,t,q)->(x,y,z)
       !   NOTE: cx*rho0/(cy*q_0=1) is used.
        gdomgdx = gdomgdr/cx + s_hat*wzz*gdomgdq/cy
        gdomgdy = - gdomgdq/cy
        gdomgdz = gdomgdt + q_0*gdomgdq
        ggxx = cx**2*ggrr
        ggxy = cx**2*s_hat*wzz*ggrr + cx*cy*(q_0*ggrt - ggrq)
        ggxz = cx*ggrt
        ggyy = (cx*s_hat*wzz)**2*ggrr + 2._DP*cx*cy*s_hat*wzz*(q_0*ggrt-ggrq) &
             + (cy*q_0)**2*ggtt - 2._DP*cy**2*q_0*ggtq + cy**2*ggqq
        ggyz = cx*s_hat*wzz*ggrt + cy*q_0*ggtt - cy*ggtq
        ggzz = ggtt
        grootg_xyz = grootg_rtq/(cx*cy)

       ! store (x,y,z)
        self%domgdx(iz) = gdomgdx
        self%domgdy(iz) = gdomgdy
        self%domgdz(iz) = gdomgdz
        self%gxx(iz) = ggxx
        self%gxy(iz) = ggxy
        self%gxz(iz) = ggxz
        self%gyy(iz) = ggyy
        self%gyz(iz) = ggyz
        self%gzz(iz) = ggzz
        self%rootg_xyz(iz) = grootg_xyz
      end do

  END SUBROUTINE metric_global_rtq2xyz


!--------------------------------------
  SUBROUTINE metric_fourier_init(self)
!--------------------------------------
    implicit none
    class(metric_fourier), intent(inout) :: self
    real(kind=DP) :: kzmin
    integer :: iz

      kzmin = 2._DP * pi / (2._DP * lz)
      do iz = -global_nz, global_nz-1
        self%kz(iz) = iz * kzmin
      end do

  END SUBROUTINE metric_fourier_init


!--------------------------------------
  SUBROUTINE forward_dft_globalz(zz_global,kz,fz,fk)
!--------------------------------------
    implicit none
    real(kind=DP), intent(in),     &
      dimension(-global_nz:global_nz-1) :: zz_global, kz, fz
    complex(kind=DP), intent(out), & 
      dimension(-global_nz:global_nz-1) :: fk
    integer :: iz, mz

      fk(:) = (0._DP, 0._DP)
      do mz = -global_nz, global_nz-1
        do iz = -global_nz, global_nz-1
          fk(mz) = fk(mz) + fz(iz)*exp(-ui*kz(mz)*zz_global(iz))*dz/(2._DP*lz)
        end do
      end do

  END SUBROUTINE forward_dft_globalz


!!--------------------------------------
!  SUBROUTINE backward_dft_globalz(zz_global,kz,fk,fz)
!!--------------------------------------
!    implicit none
!    real(kind=DP), intent(in),     &
!      dimension(-global_nz:global_nz-1) :: zz_global, kz
!    complex(kind=DP), intent(in), & 
!      dimension(-global_nz:global_nz-1) :: fk
!    real(kind=DP), intent(out),     &
!      dimension(-global_nz:global_nz-1) :: fz
!    integer :: iz, mz
!
!      fz(:) = 0._DP
!      do iz = -global_nz, global_nz-1
!        do mz = -global_nz, global_nz-1
!          fz(iz) = fz(iz) + real(fk(mz)*exp(ui*kz(mz)*zz_global(iz)), kind=DP)
!        end do
!      end do
!
!  END SUBROUTINE backward_dft_globalz


!--------------------------------------
  SUBROUTINE backward_dft_localz(zz_local,kz,fk,fz)
!--------------------------------------
    implicit none
    real(kind=DP), intent(in),     &
      dimension(-nz:nz-1) :: zz_local
    real(kind=DP), intent(in),     &
      dimension(-global_nz:global_nz-1) :: kz
    complex(kind=DP), intent(in), & 
      dimension(-global_nz:global_nz-1) :: fk
    real(kind=DP), intent(out),     &
      dimension(-nz:nz-1) :: fz
    integer :: iz, mz

      fz(:) = 0._DP
      do iz = -nz, nz-1
        do mz = -global_nz, global_nz-1
          fz(iz) = fz(iz) + real(fk(mz)*exp(ui*kz(mz)*zz_local(iz)), kind=DP)
        end do
      end do

  END SUBROUTINE backward_dft_localz


!--------------------------------------
  SUBROUTINE metric_fourier_dft_rtq2coef(self, mtr_g)
!--------------------------------------
    implicit none
    class(metric_fourier), intent(inout) :: self
    class(metric_global), intent(in) :: mtr_g
    real(kind=DP), dimension(-global_nz:global_nz-1) :: theta_tilde

     ! theta = zz + theta_tilde(zz), theta_tilde is a periodic function.
      if (trim(equib_type) == "vmec") then
        theta_tilde = mtr_g%theta - q_0 * mtr_g%zz ! Axisymmetric toroidal angle phi_ax
      else
        theta_tilde = mtr_g%theta - mtr_g%zz
      end if
      call forward_dft_globalz(mtr_g%zz, self%kz, theta_tilde,   self%theta_tilde)
      call forward_dft_globalz(mtr_g%zz, self%kz, mtr_g%omg      , self%omg      )
      call forward_dft_globalz(mtr_g%zz, self%kz, mtr_g%domgdr   , self%domgdr   )
      call forward_dft_globalz(mtr_g%zz, self%kz, mtr_g%domgdt   , self%domgdt   )
      call forward_dft_globalz(mtr_g%zz, self%kz, mtr_g%domgdq   , self%domgdq   )
      call forward_dft_globalz(mtr_g%zz, self%kz, mtr_g%grr      , self%grr      )
      call forward_dft_globalz(mtr_g%zz, self%kz, mtr_g%grt      , self%grt      )
      call forward_dft_globalz(mtr_g%zz, self%kz, mtr_g%grq      , self%grq      )
      call forward_dft_globalz(mtr_g%zz, self%kz, mtr_g%gtt      , self%gtt      )
      call forward_dft_globalz(mtr_g%zz, self%kz, mtr_g%gtq      , self%gtq      )
      call forward_dft_globalz(mtr_g%zz, self%kz, mtr_g%gqq      , self%gqq      )
      call forward_dft_globalz(mtr_g%zz, self%kz, mtr_g%rootg_rtq, self%rootg_rtq)
      ! NOTE:
      ! Arguments are (zz_global(in),kz_global(in),omg_global(in),coef_global(out))

  END SUBROUTINE metric_fourier_dft_rtq2coef


!--------------------------------------
  SUBROUTINE metric_local_dft_coef2rtq(self, mtr_f)
!--------------------------------------
    implicit none
    class(metric_local), intent(inout) :: self
    class(metric_fourier), intent(in) :: mtr_f 
    real(kind=DP), dimension(-nz:nz-1) :: theta_tilde

     ! theta = zz + theta_tilde(zz), theta_tilde is a periodic function.
      call backward_dft_localz(self%zz_labframe, mtr_f%kz, mtr_f%theta_tilde, theta_tilde )
      if (trim(equib_type) == "vmec") then
        self%theta = q_0 * self%zz_labframe + theta_tilde ! Axisymmetric toroidal angle phi_ax = q_0*zz + phi_tilde(zz)
      else
        self%theta = self%zz_labframe + theta_tilde
      end if
      call backward_dft_localz(self%zz_labframe, mtr_f%kz, mtr_f%omg      , self%omg      )
      call backward_dft_localz(self%zz_labframe, mtr_f%kz, mtr_f%domgdr   , self%domgdr   )
      call backward_dft_localz(self%zz_labframe, mtr_f%kz, mtr_f%domgdt   , self%domgdt   )
      call backward_dft_localz(self%zz_labframe, mtr_f%kz, mtr_f%domgdq   , self%domgdq   )
      call backward_dft_localz(self%zz_labframe, mtr_f%kz, mtr_f%grr      , self%grr      )
      call backward_dft_localz(self%zz_labframe, mtr_f%kz, mtr_f%grt      , self%grt      )
      call backward_dft_localz(self%zz_labframe, mtr_f%kz, mtr_f%grq      , self%grq      )
      call backward_dft_localz(self%zz_labframe, mtr_f%kz, mtr_f%gtt      , self%gtt      )
      call backward_dft_localz(self%zz_labframe, mtr_f%kz, mtr_f%gtq      , self%gtq      )
      call backward_dft_localz(self%zz_labframe, mtr_f%kz, mtr_f%gqq      , self%gqq      )
      call backward_dft_localz(self%zz_labframe, mtr_f%kz, mtr_f%rootg_rtq, self%rootg_rtq)
      ! NOTE:
      ! Arguments are (zz_local(in),kz_global(in), coef_global(in), omg_local(out)).
      ! Fourier coefficients have been evaluated in the lab frame at t=0.
      ! self%zz_labframe (= z''+t*gamma_e/s_hat) is the time-dependent flux-coordinate theta in the lab frame.

  END SUBROUTINE metric_local_dft_coef2rtq


!--------------------------------------
  SUBROUTINE metric_local_rtq2xyz(self)
!--------------------------------------
    implicit none
    class(metric_local), intent(inout) :: self
    real(kind=DP) :: wzz
    real(kind=DP) :: gdomgdx, gdomgdy, gdomgdz, &
                     ggxx, ggxy, ggxz, ggyy, ggyz, ggzz, grootg_xyz
    real(kind=DP) :: gdomgdr, gdomgdt, gdomgdq, &
                     ggrr, ggrt, ggrq, ggtt, ggtq, ggqq, grootg_rtq
    integer :: iz

      do iz = -nz, nz-1
       ! load (r,t,q)=(rho,theta,zeta)
        wzz = self%zz(iz)
        gdomgdr = self%domgdr(iz)
        gdomgdt = self%domgdt(iz)
        gdomgdq = self%domgdq(iz)
        ggrr = self%grr(iz)
        ggrt = self%grt(iz)
        ggrq = self%grq(iz)
        ggtt = self%gtt(iz)
        ggtq = self%gtq(iz)
        ggqq = self%gqq(iz)
        grootg_rtq = self%rootg_rtq(iz)

       ! translate (r,t,q)->(x,y,z)
       !   NOTE: cx*rho0/(cy*q_0=1) is used.
        gdomgdx = gdomgdr/cx + s_hat*wzz*gdomgdq/cy
        gdomgdy = - gdomgdq/cy
        gdomgdz = gdomgdt + q_0*gdomgdq
        ggxx = cx**2*ggrr
        ggxy = cx**2*s_hat*wzz*ggrr + cx*cy*(q_0*ggrt - ggrq)
        ggxz = cx*ggrt
        ggyy = (cx*s_hat*wzz)**2*ggrr + 2._DP*cx*cy*s_hat*wzz*(q_0*ggrt-ggrq) &
             + (cy*q_0)**2*ggtt - 2._DP*cy**2*q_0*ggtq + cy**2*ggqq
        ggyz = cx*s_hat*wzz*ggrt + cy*q_0*ggtt - cy*ggtq
        ggzz = ggtt
        grootg_xyz = grootg_rtq/(cx*cy)

       ! store (x,y,z)
        self%domgdx(iz) = gdomgdx
        self%domgdy(iz) = gdomgdy
        self%domgdz(iz) = gdomgdz
        self%gxx(iz) = ggxx
        self%gxy(iz) = ggxy
        self%gxz(iz) = ggxz
        self%gyy(iz) = ggyy
        self%gyz(iz) = ggyz
        self%gzz(iz) = ggzz
        self%rootg_xyz(iz) = grootg_xyz
      end do

  END SUBROUTINE metric_local_rtq2xyz


!--------------------------------------
  SUBROUTINE metric_local_copy_global(self, mtr_g)
!--------------------------------------
    implicit none
    class(metric_local), intent(inout) :: self
    class(metric_global), intent(in) :: mtr_g
    integer :: iz, giz

      do iz = -nz, nz-1
        giz = iz - global_nz + 2*nz * rankz + nz
        self%zz_labframe(iz) = mtr_g%zz(giz)    
        self%zz(iz)        = mtr_g%zz(giz)    
        self%theta(iz)     = mtr_g%theta(giz) 
        self%omg(iz)       = mtr_g%omg(giz)   
        self%domgdx(iz)    = mtr_g%domgdx(giz)
        self%domgdy(iz)    = mtr_g%domgdy(giz)
        self%domgdz(iz)    = mtr_g%domgdz(giz)
        self%gxx(iz)       = mtr_g%gxx(giz)
        self%gxy(iz)       = mtr_g%gxy(giz)
        self%gxz(iz)       = mtr_g%gxz(giz)
        self%gyy(iz)       = mtr_g%gyy(giz)
        self%gyz(iz)       = mtr_g%gyz(giz)
        self%gzz(iz)       = mtr_g%gzz(giz)
        self%rootg_xyz(iz) = mtr_g%rootg_xyz(giz)
      end do

  END SUBROUTINE metric_local_copy_global


!--------------------------------------
  SUBROUTINE metric_local_init(self, mtr_f, time_shearflow)
!--------------------------------------
    implicit none
    class(metric_local), intent(inout) :: self
    class(metric_fourier), intent(in) :: mtr_f
    real(kind=DP), intent(in) :: time_shearflow

      self%zz(:) = zz(:)
      self%zz_labframe(:) = zz(:) + time_shearflow * gamma_e / s_hat
      call self%dft_coef2rtq(mtr_f)
      call self%rtq2xyz
      
  END SUBROUTINE metric_local_init


!--------------------------------------
  SUBROUTINE metric_local_update(self, mtr_f, dt_shearflow)
!--------------------------------------
    implicit none
    class(metric_local), intent(inout) :: self
    class(metric_fourier), intent(in) :: mtr_f
    real(kind=DP), intent(in) :: dt_shearflow

      self%zz_labframe(:) = self%zz_labframe(:) + dt_shearflow * gamma_e / s_hat
      call self%dft_coef2rtq(mtr_f)
      call self%rtq2xyz
      
  END SUBROUTINE metric_local_update

      
END MODULE GKV_geom
