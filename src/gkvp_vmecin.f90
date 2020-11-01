MODULE GKV_vmecin
!-------------------------------------------------------------------------------
!
!    Calculate the magnetic field components and metric coefficients
!      from the VMEC equilibrium
!
!      GKV-plus r0.6 ( T.-H.Watanabe and M. Nunami, Dec 2011 )
!
!      This module utilizes subroutines developed by M. Nunami
!         for the GKV-X code
!
!-------------------------------------------------------------------------------

  use GKV_header
  use GKV_mpienv

  implicit none

  private

  public :: vmecin_fileopen, vmecin_coeff, vmecin_read


!*****************************************************************************************
!****************                                                         ****************
!********                    shot# 088343   t = 1.833 [s]                         ********
!****************                                                         ****************
!*****************************************************************************************

!== Paramters at rho=0.65 ======I N W A R D============================================
!-*** common parameter at rho=0.65 ***------------------------------------------------
!   real(kind=DP), parameter :: Ln_unit =-4.230701_DP                      ! Ln [m]
!   real(kind=DP), parameter :: Lt_unit = 0.3135611_DP                     ! Lt [m]
!-------------------------------------------------------------------------------------

!-- For inward-LHD vacuum" -----------------------------------------------------------
!   real(kind=DP), parameter :: R0_unit = 3.599858_DP                      ! R0 [m]
!   real(kind=DP), parameter :: r_edge  = 0.6362872D0                      ! r_edge [m]
!   real(kind=DP), parameter :: b0b00   = 2.88846853946973647d0/2.940307D0 ! b00mode/B0
!-------------------------------------------------------------------------------------


!========================================================================
!==  MODULES variables FOR VMEC                                                            
!========================================================================

!!!!! === MODULE NEWBOZ
! ... originally contained param.
  integer(kind=4), parameter :: NBZMIN = -20,NBZMAX = 20, MBZMAX = 100, &
                                IBHAR  = ( NBZMAX - NBZMIN + 1 )*MBZMAX 
! ... from pmnsd
  integer(kind=4), parameter :: NSD_MAX=501

! ... from pmims
! integer(kind=4), parameter :: KMSH_MAX=NSD_MAX-1, MDMX_MAX=101
! ... from newboz
  real(kind=8), save, dimension(IBHAR,NSD_MAX) :: BBOZH,RBOZH,ZBOZH,PBOZH
  real(kind=8), save, dimension(NSD_MAX)       :: PSIBZ,EOTBZ,CUIBZ,CUGBZ
  integer(kind=4), save, dimension(IBHAR)      :: MBOZ,NBOZ
  integer(kind=4), save                        :: NFP,NMBOZ

!!!!! === MODULE bcoef
  real(kind=8), save, dimension(:,:), allocatable :: BCO,C1BF,C2BF,C3BF
  real(kind=8), save, dimension(:),   allocatable :: CM,CN,SPOS,PSIB,    &
       EOT,C1ET,C2ET,C3ET,CUI,C1I,C2I,C3I,CUG,C1G,C2G,C3G,         &
       TXI,C1TI,C2TI,C3TI,TXE,C1TE,C2TE,C3TE,DXI,C1DI,C2DI,C3DI,   &
       DLN,C1LN,C2LN,C3LN     
  real(kind=8), save, dimension(:,:), allocatable :: bco0
  real(kind=8), save, dimension(:),   allocatable :: spos0

!!!!! === MODULE xcoef
  real(kind=8), save, dimension(:,:), allocatable :: RCO,ZCO,PCO,        &
                               C1R,C2R,C3R,C1Z,C2Z,C3Z,C1P,C2P,C3P
  real(kind=8), save, dimension(:,:), allocatable :: rco0,zco0,pco0

!!!!! === MODULE pcoef
  real(kind=8), save, dimension(0:10) :: dn0,ti0,te0

!!!!! === MODULE param1
  integer(kind=4), save :: kmsh, kmsh1, nsd, mdmx, itype
!
!!!!! === MODULE bpara
  real(kind=8), save    :: bb0,psia,sa,rmaj0,zi,nmass

  real(kind=DP), save   :: ss, cug1, cui1

  integer, save :: jp0

!========================================================================

  integer :: ndiskc = inbz
!!!  integer :: ndiskc = 14
!!!  integer :: ivmc = 15, ovmc = olog  ! move to GKV_header


CONTAINS


!--------------------------------------
  SUBROUTINE vmecin_fileopen
!--------------------------------------

    implicit none

    character(512) :: f_nbz, f_vmc

    namelist /vmecf/ f_nbz, f_vmc


      read(inml,nml=vmecf)

      write(olog,*) "# newboz and vmec input files : "
      write(olog,*) trim(f_nbz), trim(f_vmc)

      open( inbz, file=trim(f_nbz), form="unformatted", &
                  status="old", action="read" )
      open( ivmc, file=trim(f_vmc), &
                  status="old", action="read" )


  END SUBROUTINE vmecin_fileopen


!--------------------------------------
  SUBROUTINE vmecin_coeff( rad_a, R0_unit, rho2R_0, q_input, theta, alpha_fix, &
                           r_0, r_minor, s_hat,                                &
                           gdwss, gdwtt, gdwzz, gdwst, gdwsz, gdwtz,           &
                           gupss, guptt, gupzz, gupst, gupsz, guptz,           &
                           babs, Bs, Bth, Bzt, dBds, dBdt, dBdz,               &
                           dBdt_mir, rootg, rootgft, rootgbz )
!--------------------------------------
    implicit none

! --- arguments

    real(kind=DP), intent(in)  :: rad_a, R0_unit, rho2R_0, q_input, theta, alpha_fix
    real(kind=DP), intent(out) :: r_0, r_minor, s_hat 
    real(kind=DP), intent(out) ::                          &
                 gdwss, gdwtt, gdwzz, gdwst, gdwsz, gdwtz, &
                 gupss, guptt, gupzz, gupst, gupsz, guptz, &
                 babs, Bs  , Bth , Bzt , dBds, dBdt, dBdz, &
                 dBdt_mir, rootg, rootgft, rootgbz


! -- For read_VMEC routine -------------------------------------
!!!    integer, parameter :: npmax=100
!!!    real(kind=DP)  :: ssi(npmax)
    real(kind=DP), dimension(:), allocatable  :: ssi

    integer, save  :: npsi=0,ntheta,nzeta

    real(kind=DP)  :: zeta,rmaj,ph,cx,sx
    real(kind=DP)  :: dRds,dZds,dPds
    real(kind=DP)  :: dRdt,dZdt,dPdt
    real(kind=DP)  :: dRdz,dZdz,dPdz
    real(kind=DP)  :: bci, rci, zci, pci
    real(kind=DP)  :: dbci,drci,dzci,dpci

    real(kind=DP)  :: rg2inv

!!!12.22    real(kind=DP)  :: Bs, Bth, Bzt

    real(kind=DP)  :: B00mode, dB00mode
    real(kind=DP)  :: Bm1_0mode,   dBm1_0mode
    real(kind=DP)  :: Bm1_10mode,   B0_10mode,  Bp1_10mode
    real(kind=DP)  :: dBm1_10mode,  dB0_10mode, dBp1_10mode

    real(kind=DP)  :: ds, ss, r_a
    real(kind=DP)  :: q_00, dq00, eot0, diff_q, q_vmec
    real(kind=DP)  :: dq_0, eot1, cug1, cui1
    real(kind=DP)  :: dsfix
    integer        :: inm, is, iis, jp, jp0

    real(kind=DP)  :: q_0

!!!12.22    real(kind=DP)  :: rho_unit, L_n, L_t, rho2R_0  ! nunami (10.04.15)


    integer, save :: isw = 0
  

    namelist /zahyoin/npsi,ntheta,nzeta


    if( isw == 0 ) then
      read(inml,nml=zahyoin)
    end if

    allocate ( ssi(npsi+1) )


    r_a      = rad_a * rho2R_0             ! --- rad_a / R_0

    do iis=1,npsi+1
      ssi(iis) = dble(iis-1)/dble(npsi)
      ss       = ssi(iis)
      jp0      = int(ss*dble(kmsh))+1
      ds       = ss - spos(jp0)

! -- For Safety Factor----------------------------------------------------------------------
      eot0  =   eot(jp0)+(c1et(jp0)+(c2et(jp0)+c3et(jp0)*ds)*ds)*ds          ! rotational transform
      q_00  =   1.d0/eot0                                                    ! safety factor
      dq00  = - q_00 * q_00 * (c1et(jp0)+(2.d0*c2et(jp0)+3.d0*c3et(jp0)*ds)*ds)  ! dq/ds
! ------------------------------------------------------------------------------------------

      diff_q  = q_input - q_00

      if(abs(diff_q) < 0.02d0 ) then
        q_0   = q_input
        q_vmec= q_00
        eot1  = eot0
        is    = iis
        dq_0  = dq00
        jp    = jp0
        dsfix = ds
      endif

    enddo

    ss = ssi(is)
    ds = dsfix

! This is to adjust to Nunami's version (Dec 2011)
! It should be confirmed how q is determined in consistent with r_minor
    r_minor  = rad_a * ss                       ! in unit of rho

    s_hat    = (ss / q_0) * dq_0

    r_0  = r_minor * rho2R_0    !! by nunami (10.04.18)

!!!    s_hat0   = (r_0 / q_0) * dq_0

! -- troidal & poloidal current ------------------------------------------------------------
    cug1  = cug(jp) + ( c1g(jp) + ( c2g(jp) + c3g(jp) * ds ) * ds ) * ds   ! G : poloidal
    cui1  = cui(jp) + ( c1i(jp) + ( c2i(jp) + c3i(jp) * ds ) * ds ) * ds   ! I : toroidal
! ------------------------------------------------------------------------------------------



! ======================================================================================

!!!    th_shift = -alpha_fix / q_0

!!!      theta = zz(iz)                   ! theta
      zeta  = theta * q_0 + alpha_fix  ! from alpha = zeta - q*theta = 0

      rmaj = 0.0d0
      babs = 0.0d0

      dBds = 0.0d0
      dRds = 0.0d0
      dZds = 0.0d0
      dPds = 0.0d0

      dBdt = 0.0d0
      dRdt = 0.0d0
      dZdt = 0.0d0
      dPdt = 0.0d0

      dBdt_mir = 0.0d0

      dBdz = 0.0d0
      dRdz = 0.0d0
      dZdz = 0.0d0
      dPdz = 1.0d0   ! !!should be 1.0, not 0.0!!

! --- Summation of Fourier components ---
      do inm=1, mdmx+1
        bci  = bco(jp,inm) + ( c1bf(jp,inm)+ ( c2bf(jp,inm) + c3bf(jp,inm) * ds ) * ds ) * ds
        rci  = rco(jp,inm) + ( c1r(jp,inm) + ( c2r(jp,inm)  + c3r(jp,inm)  * ds ) * ds ) * ds
        zci  = zco(jp,inm) + ( c1z(jp,inm) + ( c2z(jp,inm)  + c3z(jp,inm)  * ds ) * ds ) * ds
        pci  = pco(jp,inm) + ( c1p(jp,inm) + ( c2p(jp,inm)  + c3p(jp,inm)  * ds ) * ds ) * ds
        dbci = c1bf(jp,inm)+ ( 2.d0 * c2bf(jp,inm) + 3.d0 * c3bf(jp,inm) * ds ) * ds 
        drci = c1r(jp,inm) + ( 2.d0 * c2r(jp,inm)  + 3.d0 * c3r(jp,inm)  * ds ) * ds
        dzci = c1z(jp,inm) + ( 2.d0 * c2z(jp,inm)  + 3.d0 * c3z(jp,inm)  * ds ) * ds
        dpci = c1p(jp,inm) + ( 2.d0 * c2p(jp,inm)  + 3.d0 * c3p(jp,inm)  * ds ) * ds

        if(cn(inm)==0 .and. cm(inm)==0) then
          B00mode = bci !bco(jp,inm)               ! (m,n)=(0,0) component of B
          dB00mode = dbci
        endif

        if(cn(inm)==0) then
          if(cm(inm)==-1) then
             Bm1_0mode  =  bci
            dBm1_0mode  = dbci
          endif
        endif


        if(cn(inm)==-10) then
          if(cm(inm)==-2) then                 ! (m,n)=(0,10) mode
             B0_10mode  =  bci
            dB0_10mode  = dbci
          elseif(cm(inm)==-3) then
             Bm1_10mode =  bci
            dBm1_10mode = dbci
          elseif(cm(inm)==-1) then
             Bp1_10mode =  bci
            dBp1_10mode = dbci
          endif
        endif

        ph = cn(inm) * zeta - cm(inm) * theta
        cx = cos(ph)                          ! cos(n*zeta - m*theta)
        sx = sin(ph)                          ! sin(n*zeta - m*theta)

        rmaj =  rci * cx + rmaj               ! R in (R,Phi,Z)
        babs =  bci * cx + babs               ! Absolute of B 

        dBds =  dbci * cx + dBds              ! dB/ds
        dRds =  drci * cx + dRds              ! dR/ds
        dZds =  dzci * sx + dZds              ! dZ/ds
        dPds =  dpci * sx + dPds              ! d(phi)/ds

        dBdt =  bci  * sx * cm(inm) + dBdt    ! dB/d(theta)
        dRdt =  rci  * sx * cm(inm) + dRdt    ! dR/d(theta)
        dZdt = -zci  * cx * cm(inm) + dZdt    ! dZ/d(theta)
        dPdt = -pci  * cx * cm(inm) + dPdt    ! d(phi)/d(theta)

!-- dB/d(theta) for mirror term "mir"  --- 100109 nunami -----------------
        dBdt_mir =  bci  * sx * (cm(inm)-cn(inm)*q_0) + dBdt_mir
!-------------------------------------------------------------------------

        dBdz = -bci  * sx * cn(inm) + dBdz    ! dB/d(zeta)
        dRdz = -rci  * sx * cn(inm) + dRdz    ! dR/d(zeta)
        dZdz =  zci  * cx * cn(inm) + dZdz    ! dZ/d(zeta)
        dPdz =  pci  * cx * cn(inm) + dPdz    ! d(phi)/d(zeta)
      end do
! ------------------------------------------------------------------------
!
! ---- Covariant componets of the metric tensor --------------------------------------
      gdwss  = dRds**2   + dZds**2   + (rmaj*dPds)**2         ! g_s_s
      gdwtt  = dRdt**2   + dZdt**2   + (rmaj*dPdt)**2         ! g_theta_theta
      gdwzz  = dRdz**2   + dZdz**2   + (rmaj*dPdz)**2         ! g_zeta_zeta
      gdwst  = dRds*dRdt + dZds*dZdt + (rmaj**2)*dPds*dPdt    ! g_s_theta
      gdwsz  = dRds*dRdz + dZds*dZdz + (rmaj**2)*dPds*dPdz    ! g_s_zeta
      gdwtz  = dRdt*dRdz + dZdt*dZdz + (rmaj**2)*dPdt*dPdz    ! g_theta_zeta
! ------------------------------------------------------------------------------------
!
! ---- For Jacobian in (s,theta,zeta) coordinates ------------------------------------
      rootg  = (2.d0 * ss * psia) * ( cug1 + eot1 * cui1 ) /  babs / babs  ! sqrt(g)
      rg2inv = 1.0d0 / rootg / rootg                          ! 1/g

!!! THW (Dec 7, 2011)
!!!      jcob(iz) = rootg
!!! THW (Dec 7, 2011)

! ------------------------------------------------------------------------------------
!
! ---- Contravariant componets of the metric tensor ----------------------------------
      gupss  = rg2inv * ( gdwtt * gdwzz - gdwtz * gdwtz )     ! g^s^s
      guptt  = rg2inv * ( gdwss * gdwzz - gdwsz * gdwsz )     ! g^theta^theta
      gupzz  = rg2inv * ( gdwss * gdwtt - gdwst * gdwst )     ! g^zeta^zeta 
      gupst  = rg2inv * ( gdwtz * gdwsz - gdwst * gdwzz )     ! g^s^theta
      gupsz  = rg2inv * ( gdwst * gdwtz - gdwtt * gdwsz )     ! g^s^zeta
      guptz  = rg2inv * ( gdwst * gdwsz - gdwss * gdwtz )     ! g^theta^zeta
! ------------------------------------------------------------------------------------
!
! ---- Covariant componets of B ------------------------------------------------------
      Bs   = ( gdwsz + eot1 * gdwst ) / rootg                 ! B_s
      Bth  = cui1                                             ! B_theta
      Bzt  = cug1                                             ! B_zeta
! ------------------------------------------------------------------------------------
!

! *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*= Edit by nunami 10.04.15 =*=*=*=*=*=*=*=*
! ---- Normalization for length with R_0 ---------------------------------------------
      rmaj  = rmaj  / R0_unit
!     ---------------------------------
      gdwss = gdwss / R0_unit / R0_unit
      gdwtt = gdwtt / R0_unit / R0_unit
      gdwzz = gdwzz / R0_unit / R0_unit
      gdwst = gdwst / R0_unit / R0_unit
      gdwsz = gdwsz / R0_unit / R0_unit
      gdwtz = gdwtz / R0_unit / R0_unit
!     ---------------------------------
      gupss = gupss * R0_unit * R0_unit
      guptt = guptt * R0_unit * R0_unit
      gupzz = gupzz * R0_unit * R0_unit
      gupst = gupst * R0_unit * R0_unit
      gupsz = gupsz * R0_unit * R0_unit
      guptz = guptz * R0_unit * R0_unit
! ------------------------------------------------------------------------------------
!
!
! ---- Normalization for B componets etc. --------------------------------------------
      babs  = babs / B00mode
      Bs    = Bs   / B00mode / R0_unit
      Bth   = Bth  / B00mode / R0_unit
      Bzt   = Bzt  / B00mode / R0_unit
      dBds  = dBds / B00mode
      dBdt  = dBdt / B00mode
      dBdz  = dBdz / B00mode

      dBdt_mir  = dBdt_mir / B00mode
! ------------------------------------------------------------------------------------

! *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*= Edit by nunami 10.04.15 =*=*=*=*=*=*=*=*


! --- Normalization for rootg with R_0 -- by nunami 10.04.15 -------------------------
      rootg = rootg / R0_unit / R0_unit / R0_unit

! ---- sqrt(g) in (x,y,z), i.e. flux tube coordinates --------------------------------
      rootgft   = ( q_0 / (r_0 * r_a) ) * rootg  ! sqrt(g_ft)

! ---- sqrt(g) in (r,theta,zeta), i.e. Boozer coordinates ----------------------------
      rootgbz   = rootg / r_a


! ======================================================================================


      if( isw == 0 ) then

        write( olog, * ) "================= "
        write( olog, * ) " #  B00mode     = ",  B00mode
        write( olog, * ) " #  Bm1_0mode   = ",  Bm1_0mode
        write( olog, * ) " #  B0_10mode   = ",  B0_10mode
        write( olog, * ) " #  Bm1_10mode  = ",  Bm1_10mode
        write( olog, * ) " #  Bp1_10mode  = ",  Bp1_10mode
        write( olog, * ) " # dB00mode     = ", dB00mode
        write( olog, * ) " # dBm1_0mode   = ", dBm1_0mode
        write( olog, * ) " # dB0_10mode   = ", dB0_10mode
        write( olog, * ) " # dBm1_10mode  = ", dBm1_10mode
        write( olog, * ) " # dBp1_10mode  = ", dBp1_10mode
        write( olog, * ) "================="

        write( olog, * ) ""

        write( olog, * ) "================= "
        write( olog, * ) " # eps_t        = ",  - Bm1_0mode  / B00mode
        write( olog, * ) " # eps_h/eps_t  = ",  - B0_10mode  / Bm1_0mode
        write( olog, * ) " # eps_-/eps_t  = ",  - Bp1_10mode / Bm1_0mode
        write( olog, * ) " # eps_+/eps_t  = ",  - Bm1_10mode / Bm1_0mode
        write( olog, * ) " # rdeps00/eps_t= ",  ss * dB00mode    / Bm1_0mode
        write( olog, * ) " # rdeps_t/eps_t= ",  ss * dBm1_0mode  / Bm1_0mode
        write( olog, * ) " # rdeps_h/eps_t= ",  -ss * dB0_10mode  / Bm1_0mode
        write( olog, * ) " # rdeps_-/eps_t= ",  -ss * dBp1_10mode / Bm1_0mode
        write( olog, * ) " # rdeps_+/eps_t= ",  -ss * dBm1_10mode / Bm1_0mode
        write( olog, * ) "================="

      end if

      isw = 1


    deallocate ( ssi )


  END SUBROUTINE vmecin_coeff



!--------------------------------------
  SUBROUTINE vmecin_read
!--------------------------------------

! +------------------------------------------------------------------------+
! |    Read VMEC equiblium                                                 |
! |    (programme for interface vmec(newboz))                              | 
! |                                                                        |
! |                                  by M. Nunami (Nov. 2009)              |
! +------------------------------------------------------------------------+


! --- read vmec-nwboz field ----------
    call iodisk

! --- calculate spline coeff ---------
    call setfld


  END SUBROUTINE vmecin_read


!--------------------------------------
  SUBROUTINE iodisk
!--------------------------------------

! +-------------------------------------------------------------------------+
! |   iodisk ; read input data from newboz                                  |
! |                                                                         |
! |   nsd          number of vmec grid points                               |
! |   mbzmax       number of boozer theta harmonics desired                 |
! |   nbzmin       number of negative boozer zeta harmonics desired         |
! |   nbzmax       number of positive boozer zeta harmonics desired         |
! |                                                                         |
! |   bbozh : mod b spectrum in boozer coordinates                          |
! |                                                                         |
! |         b(i,t,z) =     sum bbozh(m,i) cos(-t*mboz(m) + z*nboz(m) )      |
! |                                                                         |
! |   note   bco(i,m) = bbozh(m,i)                                          |
! |                                                                         |
! |   transformation from boozer to cylindrical coordinates                 |
! |                                                                         |
! |       ( psi, theta, zeta ) -> ( r, phi, z )                             |
! |                                                                         |
! |         r(i,t,z) =     sum rbozh(m,i) cos(-t*mboz(m) + z*nboz(m) )      |
! |         p(i,t,z) = z + sum pbozh(m,i) sin(-t*mboz(m) + z*nboz(m) )      |
! |         z(i,t,z) =     sum zbozh(m,i) sin(-t*mboz(m) + z*nboz(m) )      |
! |                                                                         |
! |   surface quantities ( after renormalized )                             |
! |                                                                         |
! |         psib   : toroidal flux within a flux surface / ( 2*pi )         |
! |         eot    : rotational transform                                   |
! |         cui    : toroidal current within a flux surface                 |
! |         cug    : poloidal current without a flux surface                |
! |                                                                         |
! |   ---------------------------------------------------------------       |
! |                                                                         |
! |  Input Parameters                                                       |
! |  NEWBZ:                                                                 |
! |   mdmx : number of important modes extracted from nmboz modes           |
! |              ( nmboz - 1 >= mdmx is required )                          |
! |            Note : total mode is mdmx+1 (including B_(0,0) mode)         |
! |   nlim : if nlim/=0, |n|>nlim modes are skipped.                        |
! |                                                                         |
! |   kmsh : number of radial mesh of field data for FORTEC-3D              |
! |            Note : total mesh is kmsh+1 (including rho=0.0)              |
! |   dcheck : if TRUE, original data from NEWBOZ is put out.               |
! |   lrchk  : if TRUE, check and adjust the coordinate system to RHS.      |
! |   bmag  : magnification factor of the magnetic field strength.          |
! +-------------------------------------------------------------------------+
!********************************************************************************
  implicit real(KIND=8)(a-h,o-z)

  implicit integer (i-n) ! by THW Dec 22, 2011

  logical dcheck,lrchk
!
  dimension ifsel(ibhar)
  dimension bmax(nsd_max), bmin(nsd_max)
  dimension bsel(ibhar),icont(ibhar), msel(ibhar)
!
!!  COMMON/BPARA/BB0,psia,sa,rmaj0,zi,nmass
!!  common/param1/kmsh,kmsh1,nsd,mdmx,itype
  namelist/newbz/dcheck,lrchk,mdmx,kmsh,nlim,bmag
! ----------------------------------------------------------------------------- !
!
!
  read(inml,nml=newbz)
  if ( rank == 0 ) then
    write(ovmc,nml=newbz)
  endif
!
!
! ... Read input data (the same file used to make input data for VMEC) for density and temp.
! ... The density and temperature profiles -> given by subroutine "ntfunc0" or "ntfunc1"
! ... Choose input type ( Satake's exponential style      : itype=0, 
! ...                     expansion in rho like in GSRAKE : itype=1, 
! ...                     Idomura's style                 : itype=2 )
!
  read(ivmc,*)
  read(ivmc,*) itype
  if(itype==0) len=4
  if(itype==1) len=10
  if(itype==2) len=3
!
  read(ivmc,*)
  read(ivmc,*)dn0(0:len)
  read(ivmc,*)
  read(ivmc,*)ti0(0:len)
  read(ivmc,*)
  read(ivmc,*)te0(0:len)
  read(ivmc,*)
  read(ivmc,*)dummy
  read(ivmc,*)
  read(ivmc,*)dummy,R0,dummy,a0,zi,nmass
!
  if(itype==2) then 
    dn0(2)=dn0(1)*dn0(2)*a0/R0
    ti0(2)=ti0(1)*ti0(2)*a0/R0
    te0(2)=te0(1)*te0(2)*a0/R0
  end if
!
!!    dn0: ni(0) [m^(-3)]
!!    ti0, te0: temperature at the axis [keV]
!!    zi:  ion charge number
!!    nmass : relative mass of the ion species to proton
!
!!    BB0:  toroidal field strength at the magnetic axis  [T]
!!    Rmaj0 and a0 : major and minor radius [m]
!
! --- constant parameters --------
  pi2    = twopi
  cmu0   = 4.0d0*pi*1.0d-7
!
! --- Read newboz data ------------- 
  read(ndiskc) nmboz, nsd, nfp
  read(ndiskc) Rmaj0, BB0
!
  if ( rank == 0 ) then
    write(ovmc,6001) nmboz, nsd, nfp
  endif
!
  nsd=nsd+1
  kmsh1=kmsh+1
  mmx1=mdmx+1
!
  if(nmboz-1.lt.mdmx) then
    if ( rank == 0 ) then
      write(ovmc,9001) nmboz - 1, mdmx
    endif
    stop
  else
    if ( rank == 0 ) then
      write(ovmc,9003) nmboz - 1, mdmx
    endif
  endif
!
  if ( rank == 0 ) then
    write(ovmc,9004)nsd,kmsh
  endif
!
  allocate(CUI(kmsh1), CUG(kmsh1), TXI(kmsh1), TXE(kmsh1), DXI(kmsh1), DLN(kmsh1),    &
           C1I(kmsh1), C2I(kmsh1), C3I(kmsh1), C1G(kmsh1), C2G(kmsh1), C3G(kmsh1),    &
           C1TI(kmsh1),C2TI(kmsh1),C3TI(kmsh1),C1TE(kmsh1),C2TE(kmsh1),C3TE(kmsh1),   &
           C1DI(kmsh1),C2DI(kmsh1),C3DI(kmsh1),SPOS(kmsh1),PSIB(kmsh1),EOT(kmsh1),    &
           C1ET(kmsh1),C2ET(kmsh1),C3ET(kmsh1),C1LN(kmsh1),C2LN(kmsh1),C3LN(kmsh1))
  allocate(BCO(kmsh1,MMX1),C1BF(kmsh1,MMX1),C2BF(kmsh1,MMX1),C3BF(kmsh1,MMX1),        &
           CM(MMX1),CN(MMX1))
  allocate(RCO(kmsh1,MMX1),ZCO(kmsh1,MMX1),PCO(kmsh1,MMX1),                           &
           C1R(kmsh1,MMX1),C2R(kmsh1,MMX1),C3R(kmsh1,MMX1),                           &
           C1Z(kmsh1,MMX1),C2Z(kmsh1,MMX1),C3Z(kmsh1,MMX1),                           &
           C1P(kmsh1,MMX1),C2P(kmsh1,MMX1),C3P(kmsh1,MMX1))
  allocate(bco0(nsd,mmx1),rco0(nsd,mmx1),zco0(nsd,mmx1),pco0(nsd,mmx1),               &
           spos0(nsd))
!
!
! --- Read newboz data ----------------------
  READ(ndiskc) (PSIBZ(i), i = 2, NSD)
  READ(ndiskc) (EOTBZ(i), i = 2, NSD)
!
  READ(ndiskc) (CUIBZ(i), i = 2, NSD)
  READ(ndiskc) (CUGBZ(i), i = 2, NSD)
!
!
  DO M = 1, NMBOZ
    READ(ndiskc)  MBOZ(M), NBOZ(M)
    READ(ndiskc) (BBOZH(M,I), I = 2, NSD)
  END DO
!
  DO M = 1, NMBOZ
    READ(ndiskc)  MBOZ(M), NBOZ(M)
    READ(ndiskc) (RBOZH(M,I), I = 2, NSD)
  END DO
!
  DO M = 1, NMBOZ
    READ(ndiskc)  MBOZ(M), NBOZ(M)
    READ(ndiskc) (ZBOZH(M,I), I = 2, NSD)
  END DO
!
  DO M = 1, NMBOZ
    READ(ndiskc)  MBOZ(M), NBOZ(M)
    READ(ndiskc) (PBOZH(M,I), I = 2, NSD)
  END DO
!
!
! --- magnetic field magnification -------------------------------------
  if(bmag/=1.0d0) then
    if( rank == 0 ) then
      write(ovmc,"(a,f8.3,a)")'** bfield is magnified by ',bmag,' **'
    endif
    BB0=BB0*bmag
    PSIBZ=PSIBZ*bmag
    CUIBZ=CUIBZ*bmag
    CUGBZ=CUGBZ*bmag
    BBOZH=BBOZH*bmag
  end if
!
! --- check of LHS or RHS ----------------------------------------------
  if ( .not.lrchk ) goto  999 !check lhs or rhs if lrchk = .true.
!
  ichk = nsd/2
  chi = pi2/90.0d0
  zsum  =  0.0d0
  do i = 1, nmboz
    zsum  =  zsum + zbozh(i,ichk)*dsin(-mboz(i)*chi)
  end do
!
  if( zsum .gt. 0.0d0 ) then
!  ....   counterclockwise
    drthta =  1.0d0
    if( rank == 0 ) then
      write(ovmc,*)
      write(ovmc,*)' --------------------------------------'
      write(ovmc,*)'     LHS ( theta is counterclockwise ) '
      write(ovmc,*)' --------------------------------------'
      write(ovmc,*)
    endif
  else
!  ....   clockwise
    drthta = -1.0d0
    if( rank == 0 ) then
      write(ovmc,*)
      write(ovmc,*)' --------------------------------------'
      write(ovmc,*)'       RHS ( theta is clockwise )      '
      write(ovmc,*)' --------------------------------------'
      write(ovmc,*)
    endif
  endif
!
!
!+++++ updated on  9/29 94 by N^2 : for left -> right start
!
!      - theta -> theta in VMEC
!
! irtchg = 0
  if( drthta .gt. 0.0d0 ) then
!   irtchg = 1
    do i = 2, nsd
!      change of sign of psi is due to miss-interface
      PSIBZ(i) = - PSIBZ(i)
      EOTBZ(i) = - EOTBZ(i)
      CUIBZ(i) = - CUIBZ(i)
    end do
!
    do m = 1, nmboz
      mboz(m) = - mboz(m)
    end do
!
    if( rank == 0 ) then
      write(ovmc,*)
      write(ovmc,*)' --------------------------------------'
      write(ovmc,*)'     LHS -> RHS ( theta is clockwise ) '
      write(ovmc,*)' --------------------------------------'
      write(ovmc,*)
    endif
    drthta = -1.0d0
  endif
!
!..... psi > 0 must be ensured in the usual Boozer coordinates
!
!      - curpol -> curpol
!      - curtor -> curtor in VMEC
!
! ipschg = 0
  if( psibz(2) .lt. 0.0d0 ) then
!   ipschg = 1
    do i = 2, nsd
      PSIBZ(i) = - PSIBZ(i)
      CUIBZ(i) = - CUIBZ(i)
      CUGBZ(i) = - CUGBZ(i)
    end do
!
    if( rank == 0 ) then
      write(ovmc,*)
      write(ovmc,*)' --------------------------------------'
      write(ovmc,*)'                    psi < 0 -> psi > 0 '
      write(ovmc,*)' --------------------------------------'
      write(ovmc,*)
    endif
  endif
!
!      winding law is inverted
!      - curtor -> curtor
!
  if( eotbz(nsd) .lt. 0.0d0 ) then
    sgniot =-1.0d0
  else
    sgniot = 1.0d0
  endif
! iitchg = 0
! if( litchg .and. sgniot .lt. 0.0d0 ) then
  if(sgniot .lt. 0.0d0 ) then
!   iitchg = 1
    do i = 2, nsd
      EOTBZ(i) = - EOTBZ(i)
      CUIBZ(i) = - CUIBZ(i)
    end do
    do m = 1, nmboz
      nboz(m) = - nboz(m)
      do i = 2, nsd
        pbozh(m,i) = - pbozh(m,i)
      end do
    end do
    sgniot = 1.0d0
!
    if( rank == 0 ) then
      write(ovmc,*)
      write(ovmc,*)' --------------------------------------'
      write(ovmc,*)'                  iota < 0 -> iota > 0 '
      write(ovmc,*)' --------------------------------------'
      write(ovmc,*)
    endif
  endif
!
!...  transpose LHS --> RHS end 
!
999 continue
!     
!.... from half mesh to integer mesh
!
  do i = 2, nsd-1
    psibz(i)   = ( psibz(i) + psibz(i+1) )/2.0d0
    eotbz(i)   = ( eotbz(i) + eotbz(i+1) )/2.0d0
    cuibz(i)   = ( cuibz(i) + cuibz(i+1) )/2.0d0
    cugbz(i)   = ( cugbz(i) + cugbz(i+1) )/2.0d0
    do j = 1, nmboz
      bbozh(j,i) = ( bbozh(j,i) + bbozh(j,i+1) )/2.0d0
      rbozh(j,i) = ( rbozh(j,i) + rbozh(j,i+1) )/2.0d0
      zbozh(j,i) = ( zbozh(j,i) + zbozh(j,i+1) )/2.0d0
      pbozh(j,i) = ( pbozh(j,i) + pbozh(j,i+1) )/2.0d0
    enddo
  enddo

!     
!+++++ extrapolation of the values at the magnetic axis
!
  do j = 1, nmboz
    if( mboz(j) .eq. 0 ) then
      bbozh(j,1) = 3.0d0*bbozh(j,2) - 3.0d0*bbozh(j,3) + bbozh(j,4)
      rbozh(j,1) = 3.0d0*rbozh(j,2) - 3.0d0*rbozh(j,3) + rbozh(j,4)
      zbozh(j,1) = 3.0d0*zbozh(j,2) - 3.0d0*zbozh(j,3) + zbozh(j,4)
      pbozh(j,1) = 3.0d0*pbozh(j,2) - 3.0d0*pbozh(j,3) + pbozh(j,4)
    else
      bbozh(j,1) = 0.0d0
      rbozh(j,1) = 0.0d0
      zbozh(j,1) = 0.0d0
      pbozh(j,1) = 0.0d0
    endif
  enddo

  if( psibz(nsd-1) .gt. 0.0d0 ) then
    psisgn   =  1.0d0
  else
    psisgn   = -1.0d0
  endif
  psibz(  1) = 1.0d-18*psisgn
  eotbz(  1) = 3.0d0*eotbz(  2) - 3.0d0*eotbz(  3) + eotbz(  4)
  cuibz(  1) = 3.0d0*cuibz(  2) - 3.0d0*cuibz(  3) + cuibz(  4)
  cugbz(  1) = 3.0d0*cugbz(  2) - 3.0d0*cugbz(  3) + cugbz(  4)
!
!+++++ extrapolation of the values at the boundary
!
  do j = 1, nmboz
    bbozh(j,nsd) = 3.0d0*bbozh(j,nsd-1)-3.0d0*bbozh(j,nsd-2)+bbozh(j,nsd-3)
    rbozh(j,nsd) = 3.0d0*rbozh(j,nsd-1)-3.0d0*rbozh(j,nsd-2)+rbozh(j,nsd-3)
    zbozh(j,nsd) = 3.0d0*zbozh(j,nsd-1)-3.0d0*zbozh(j,nsd-2)+zbozh(j,nsd-3)
    pbozh(j,nsd) = 3.0d0*pbozh(j,nsd-1)-3.0d0*pbozh(j,nsd-2)+pbozh(j,nsd-3)
  enddo
!
  psibz(nsd) = 3.0d0*psibz(nsd-1)-3.0d0*psibz(nsd-2)+psibz(nsd-3)
  eotbz(nsd) = 3.0d0*eotbz(nsd-1)-3.0d0*eotbz(nsd-2)+eotbz(nsd-3)
  cuibz(nsd) = 3.0d0*cuibz(nsd-1)-3.0d0*cuibz(nsd-2)+cuibz(nsd-3)
  cugbz(nsd) = 3.0d0*cugbz(nsd-1)-3.0d0*cugbz(nsd-2)+cugbz(nsd-3)
!     
!.... normalization to vmec calculation
!
  cnorm  = rmaj0*bb0/cugbz(nsd)
  do i = 1, nsd
    cuibz(i)   = cuibz(i)*cnorm
    cugbz(i)   = cugbz(i)*cnorm
    psibz(i)   = psibz(i)*cnorm
    do j = 1, nmboz
      bbozh(j,i) = bbozh(j,i)*cnorm
    enddo
  enddo

  psia = psibz(nsd)
  sa    = sqrt(psia*2.d0/bb0)
  if( rank == 0 ) then
    write(ovmc,305) rmaj0,bb0,psia,sa
    write(ovmc,306) cnorm
  endif
305 format(' rmaj0= ',d15.7,' , bb0 = ',d15.7,' , psia = ',d15.7,' , sa = ',d15.7)
306 format(' cnom =',d15.7)
!
  if( rank == 0 ) then
    write(ovmc,6002)
    write(ovmc,6003) (i, psibz(i), eotbz(i), cuibz(i), cugbz(i), i = 1, nsd)
  endif
!
  if( dcheck ) then                            !check bbozh
    nout  = 5
    nprt  = nsd/nout
    if( mod(nsd,nout) .ne. 0 ) nprt = nprt + 1
!     
    do ii=1,nprt
      iis  =(ii-1)*nout+1
      iie  = ii   *nout
      if( iie .gt. nsd) iie=nsd 
      if( rank == 0 ) then
        write(ovmc,6011)         ( kk , kk = iis, iie )
      endif
!     
      do j = 1, nmboz
        if( rank == 0 ) then
          write(ovmc,6012) j, mboz(j), nboz(j),( bbozh(j,kk), kk = iis, iie )
        endif
      enddo
    enddo
  endif
!
!        modified by Satake 2004/02/26
!
!..... assuming ijf = 1 correponds to mboz = 0 and nboz = 0.
!
!
!..... initialization of labels of fourier modes for each surface
!
  do ijf = 2, nmboz
    bsel(ijf) = -1.0d0
    msel(ijf) = 0
    icont(ijf)= 0
  end do
  icont(1)=1
  imttl   =1
!
!..... remove |n|>nlim modes 
  if(nlim.ne.0) then
    if( rank == 0 ) then
      write(ovmc,*)
!     write(ovmc,'("***    |n|>",i3," modes are skipped.   ***")'),nlim
    endif
    do ijf = 2, nmboz
      if(abs(nboz(ijf))>nlim) then
        bbozh(ijf,:)=0.0d0
      end if
    end do
  end if
!
!..... selection of modes for each surface
!
  do is = 1, nsd
!
    bmax(is)  = -100.0*bb0
    bmin(is)  =  100.0*bb0
    do ijf = 2, nmboz
      bmax(is)  = dmax1( bbozh(ijf,is), bmax(is) )
      bmin(is)  = dmin1( bbozh(ijf,is), bmin(is) )
    end do
!     
    bamax = dmax1( dabs(bmax(is)), dabs(bmin(is)) )
!     
    do ijf = 2, nmboz
      asel = dabs(bbozh(ijf,is))/bamax
      if (asel.ge.bsel(ijf)) then
        bsel(ijf) = asel
        msel(ijf) = is
      end if
    end do
  end do
!     
!.. selection of modes 
!.. search major modes of mag. field spectrum (total mdmx+1 modes)
!.. b(m=0,n=0) mode is always comes first.
!
  do ijf = 2, nmboz
    igt = 0
    do ijf2 = 2, nmboz
      if (ijf==ijf2) cycle
      if (bsel(ijf2).gt.bsel(ijf)) then
        igt = igt+1
        cycle
      end if
      if (bsel(ijf2)==1.0d0.and.bsel(ijf)==1.0d0) then
        if (dabs( bbozh(ijf2,msel(ijf2))).gt.dabs( bbozh(ijf,msel(ijf)))) then
          igt = igt+1
        end if
      end if
    end do
!     
    if (igt.lt.mdmx) then 
      icont(ijf) = igt+2
      imttl=imttl+1
    end if
  end do
!
  if( rank == 0 ) then
    write(ovmc,6035) 
  endif
!
  do i=1,imttl
    do ijf = 1, nmboz
      if (icont(ijf)==i) then
        bmax1=-100.0*bb0
        do is=1,nsd
          bmax1=dmax1(dsqrt(bbozh(ijf,is)**2),bmax1)
        end do
        if( rank == 0 ) then
          write(ovmc,6036) i,mboz(ijf),nboz(ijf),icont(ijf),msel(ijf),bsel(ijf),bmax1
        endif
        do is = 1, nsd
          bco0(is,i) = bbozh(ijf,is)
          rco0(is,i) = rbozh(ijf,is)
          zco0(is,i) = zbozh(ijf,is)
          pco0(is,i) = pbozh(ijf,is)
        end do
        cm(i) = dble(mboz(ijf))
        cn(i) = dble(nboz(ijf))
        cycle
      end if
    end do
  end do
!
  if (imttl.ne.mmx1) then
    if( rank == 0 ) then
      write(ovmc,*) 'imttl diff. from mmx1. abort! '
    endif
    stop
  end if
!
  if( rank == 0 ) then
    write(ovmc,*)
    write(ovmc,*)' selected fourier spectrum of b '
    write(ovmc,*)
  endif
!
  nout  = 5
  nprt  = nsd/nout
  if( mod(nsd,nout) .ne. 0 ) nprt = nprt + 1
!     
  do ii=1,nprt
    iis  =(ii-1)*nout+1
    iie  = ii   *nout
    if( ii .eq. nprt ) iie  = nsd
    if( rank == 0 ) then
      write(ovmc,6033)  ( kk , kk = iis, iie )
    endif
!
    do j = 1, imttl
      if( rank == 0 ) then
        write(ovmc,6034) j, idint(cm(j)), idint(cn(j)),( bco0(kk,j), kk = iis, iie )
      endif
    end do
  end do
!
6033 format(//4x,' j',2x,'(   m,   n)',5(5x,'b _',i4)/)
6034 format(  2x,i4,2x,'(',i4,',',i4,')',5(1x,1pd11.4))
6035 format(/2x,'selected mode j   (   m,   n)_j','  rank   most significant on', &
     &       '   b_max(rel.)     b_max(abs.)')
6036 format(10x,i4,'      (',i4,',',i4,')',4x,i4,8x,i6,8x,1p,2e13.5)
!
6001 format(//9x,'nmboz  = ',i4,' : nsd - 1 = ',i4,' : nfp    =',i4,//)
6002 format(//6x,'i',11x,'psi',10x,'iota',13x,'i',13x,'g'/)
6003 format(/(2x,i5,4(2x,1pd12.5)))
6011 format(//3x,' j',3x,'(   m,   n)',5(5x,'b _',i4)/)
6012 format(  2x,i4,2x,'(',i4,',',i4,')',5(1x,1pd11.4))
9001 format(//10x,'!err! nmboz - 1 = ',i4,', mdmx  = ',i4//)
9003 format(//10x,'nmboz - 1 = ',i4,' >= mdmx  = ',i4//)
!!c 9002 format(//10x,'!err! nsd   - 1 = ',i4,' =/ kmsh = ',i4//)
9004 format(//10x,'nsd = ',i4,', kmsh = ',i4//)

  return


  END SUBROUTINE iodisk



!
!
!
!--------------------------------------
  SUBROUTINE setfld
!--------------------------------------

! +-------------------------------------------------------------------------+
! |   Set field data                                                        |
! |                                                                         |
! |   Modified 2004/3/9: Spline tables are in s_bar=psi/psia coordinate.    |
! |   Modified 2005/10/24: kmsh can be taken to be different from nsd.      |
! +-------------------------------------------------------------------------+
!********************************************************************************
! ----------------------------------------------------------------------------- !
  implicit real(KIND=8)(a-h,o-z)

  implicit integer (i-n) ! by THW Dec 22, 2011

!
  real(KIND=8), dimension(kmsh1) :: wkbf,wkc1,wkc2,wkc3,ds,dsl
  real(KIND=8), dimension(nsd) :: wkbf0,wkc10,wkc20,wkc30
  integer(kind=4), dimension(kmsh1) :: jpl
  real(KIND=8) dy(2),fu(0:3)
!!  COMMON/BPARA/bb0,psia,sa,rmaj0,zi,nmass
!!  common/param1/kmsh,kmsh1,nsd,mdmx,itype
! ----------------------------------------------------------------------------- !
!
!
!------  radial position -----
!
  spos0(1) = 0.0d0
  spos(1)  = 0.0d0
  psib(1)  = psibz(1)
  do j=2,nsd
    spos0(j) = sqrt( psibz(j)/psia ) ! normalized minor radius
  end do
!
  dsp=1.0d0/dble(kmsh)
  do j=2,kmsh1
    spos(j)  = dsp*dble(j-1)
    psib(j)  = psia*spos(j)**2
    ds(j-1)  = spos(j)-spos(j-1)
  end do
  ds(kmsh1)=ds(kmsh)
!... relation between spos0 and spos
  dsl(1)=0.0d0
  jpl(1)=1
  do j=2,kmsh1
    sp = spos(j)
    jpl(j)=nsd
    do k=2,nsd
      if(sp<spos0(k)) then
        jpl(j)=k-1
        exit
      end if
    end do
    sj     = spos0(jpl(j))
    dsl(j) = sp - sj
  end do
!
!------  make n,t,ln_lambda tables on spos ---!
!
  do j=1,kmsh1
    select case(itype)
    case(0)
      call ntfunc0(spos(j),txe(j),txi(j),dxi(j))
    case(1)
      call ntfunc1(spos(j),txe(j),txi(j),dxi(j))
    case(2)
      call ntfunc2(spos(j),txe(j),txi(j),dxi(j))
    end select
    dln(j)=32.2d0+1.15d0*dlog10(txe(j)**2*1.0d6/dxi(j)/zi)
  end do
!
!--------------------------------------------!
!
!   make the spline table for magnetic field and (r,phi,z)
!   
!--------------------------------------------!
!
  do i=1,mdmx+1
!... first, make spline table on spos0
    do j=1,nsd
      wkbf0(j)=bco0(j,i)
    end do
!
    call fit_x3(spos0(1:4),wkbf0(1:4),1,dy(1),fu)
    call fit_x3(spos0(nsd-3:nsd),wkbf0(nsd-3:nsd),2,dy(2),fu)
    call spline_fit2(spos0,wkbf0,dy,nsd,wkc10,wkc20,wkc30)
!... second, interporate value at spos
    call val_intp(jpl,dsl,wkbf0,wkc10,wkc20,wkc30,kmsh1,nsd,wkbf)
!... last, make spline table on spos
    call fit_x3(spos(1:4),wkbf(1:4),1,dy(1),fu)
    call fit_x3(spos(kmsh-2:kmsh1),wkbf(kmsh-2:kmsh1),2,dy(2),fu)
    call spline_fit2(spos,wkbf,dy,kmsh1,wkc1,wkc2,wkc3)
!
    do j=1, kmsh
      bco (j,i) =  wkbf(j)
      c1bf(j,i) =  wkc1(j)
      c2bf(j,i) =  wkc2(j)
      c3bf(j,i) =  wkc3(j)
    end do
!.. on the outer edge
    bco (kmsh1,i) = fu(0)
    c1bf(kmsh1,i) = fu(1)
    c2bf(kmsh1,i) = fu(2)
    c3bf(kmsh1,i) = fu(3)
!!
    do j=1,nsd
      wkbf0(j)=rco0(j,i)
    end do
!
    call fit_x3(spos0(1:4),wkbf0(1:4),1,dy(1),fu)
    call fit_x3(spos0(nsd-3:nsd),wkbf0(nsd-3:nsd),2,dy(2),fu)
    call spline_fit2(spos0,wkbf0,dy,nsd,wkc10,wkc20,wkc30)
!
    call val_intp(jpl,dsl,wkbf0,wkc10,wkc20,wkc30,kmsh1,nsd,wkbf)
!
    call fit_x3(spos(1:4),wkbf(1:4),1,dy(1),fu)
    call fit_x3(spos(kmsh-2:kmsh1),wkbf(kmsh-2:kmsh1),2,dy(2),fu)
    call spline_fit2(spos,wkbf,dy,kmsh1,wkc1,wkc2,wkc3)
!
    do j=1, kmsh
      rco(j,i) =  wkbf(j)
      c1r(j,i) =  wkc1(j)
      c2r(j,i) =  wkc2(j)
      c3r(j,i) =  wkc3(j)
    end do
    rco(kmsh1,i) = fu(0)
    c1r(kmsh1,i) = fu(1)
    c2r(kmsh1,i) = fu(2)
    c3r(kmsh1,i) = fu(3)
!!
    do j=1,nsd
      wkbf0(j)=zco0(j,i)
    end do
!
    call fit_x3(spos0(1:4),wkbf0(1:4),1,dy(1),fu)
    call fit_x3(spos0(nsd-3:nsd),wkbf0(nsd-3:nsd),2,dy(2),fu)
    call spline_fit2(spos0,wkbf0,dy,nsd,wkc10,wkc20,wkc30)
!
    call val_intp(jpl,dsl,wkbf0,wkc10,wkc20,wkc30,kmsh1,nsd,wkbf)
!
    call fit_x3(spos(1:4),wkbf(1:4),1,dy(1),fu)
    call fit_x3(spos(kmsh-2:kmsh1),wkbf(kmsh-2:kmsh1),2,dy(2),fu)
    call spline_fit2(spos,wkbf,dy,kmsh1,wkc1,wkc2,wkc3)
!
    do j=1, kmsh1
      zco(j,i) =  wkbf(j)
      c1z(j,i) =  wkc1(j)
      c2z(j,i) =  wkc2(j)
      c3z(j,i) =  wkc3(j)
    end do
    zco(kmsh1,i) = fu(0)
    c1z(kmsh1,i) = fu(1)
    c2z(kmsh1,i) = fu(2)
    c3z(kmsh1,i) = fu(3)
!!
    do j=1,nsd
      wkbf0(j)=pco0(j,i)
    end do
!
    call fit_x3(spos0(1:4),wkbf0(1:4),1,dy(1),fu)
    call fit_x3(spos0(nsd-3:nsd),wkbf0(nsd-3:nsd),2,dy(2),fu)
    call spline_fit2(spos0,wkbf0,dy,nsd,wkc10,wkc20,wkc30)
!
    call val_intp(jpl,dsl,wkbf0,wkc10,wkc20,wkc30,kmsh1,nsd,wkbf)
!
    call fit_x3(spos(1:4),wkbf(1:4),1,dy(1),fu)
    call fit_x3(spos(kmsh-2:kmsh1),wkbf(kmsh-2:kmsh1),2,dy(2),fu)
    call spline_fit2(spos,wkbf,dy,kmsh1,wkc1,wkc2,wkc3)
!
    do j=1, kmsh1
      pco(j,i) =  wkbf(j)
      c1p(j,i) =  wkc1(j)
      c2p(j,i) =  wkc2(j)
      c3p(j,i) =  wkc3(j)
    end do
    pco(kmsh1,i) = fu(0)
    c1p(kmsh1,i) = fu(1)
    c2p(kmsh1,i) = fu(2)
    c3p(kmsh1,i) = fu(3)
!
  end do
!
!...
!
!!$if( rank == 0 ) then
!!$  write(8,*)'### spline coefficients for B ###'
!!$  write(8,*)'### kmsh    mdmx ###'
!!$  write(8,7001)kmsh,mdmx
!!$  write(8,*)'### bco          c1bf           c2bf           c3bf            (  j,  i)###'
!!$  do i=1,mdmx+1
!!$    do j=1,kmsh1
!!$      write(8,7002)bco(j,i),c1bf(j,i),c2bf(j,i),c3bf(j,i),j,i
!!$    end do
!!$  end do
!!$!
!!$  write(8,*)'### spline coefficients for R ###'
!!$  write(8,*)'### kmsh    mdmx ###'
!!$  write(8,7001)kmsh,mdmx
!!$  write(8,*)'### rco          c1r            c2r            c3r             (  j,  i)###'
!!$  do i=1,mdmx+1
!!$    write(8,7000)int(cm(i)),int(cn(i))
!!$    do j=1,kmsh1
!!$      write(8,7002)rco(j,i),c1r(j,i),c2r(j,i),c3r(j,i),j,i
!!$    end do
!!$  end do
!!$!
!!$  write(8,*)'### spline coefficients for Z ###'
!!$  write(8,*)'### kmsh    mdmx ###'
!!$  write(8,7001)kmsh,mdmx
!!$  write(8,*)'### zco          c1z            c2z            c3z             (  j,  i)###'
!!$  do i=1,mdmx+1
!!$    do j=1,kmsh1
!!$      write(8,7002)zco(j,i),c1z(j,i),c2z(j,i),c3z(j,i),j,i
!!$    end do
!!$  end do
!!$!
!!$  write(8,*)'### spline coefficients for PHI ###'
!!$  write(8,*)'### kmsh    mdmx ###'
!!$  write(8,7001)kmsh,mdmx
!!$  write(8,*)'### pco          c1p            c2p            c3p             (  j,  i)###'
!!$  do i=1,mdmx+1
!!$    do j=1,kmsh1
!!$      write(8,7002)pco(j,i),c1p(j,i),c2p(j,i),c3p(j,i),j,i
!!$    end do
!!$  end do
!!$!
!!$endif
!!
7000 format(4x,'cm = ',i4,' , cn = ',i4)
7001 format(6x,i4,4x,i4)
7002 format(4(1x,d14.6),1x,2i4)
!
!--------------------------------------------!
!
!     make the spline table for the flux functions;
!     eot, cui, cug, txi, txe, dxi, dln
!     Note that dT/dx, dn/dx (at x=0) = 0 are assumed. 
!--------------------------------------------!
!
!!!   eot
  call fit_x3(spos0(1:4),eotbz(1:4),1,dy(1),fu)
  call fit_x3(spos0(nsd-3:nsd),eotbz(nsd-3:nsd),2,dy(2),fu)
  call spline_fit2(spos0,eotbz,dy,nsd,wkc10,wkc20,wkc30)
!
  call val_intp(jpl,dsl,eotbz,wkc10,wkc20,wkc30,kmsh1,nsd,eot)
!
  call fit_x3(spos(1:4),eot(1:4),1,dy(1),fu)
  call fit_x3(spos(kmsh-2:kmsh1),eot(kmsh-2:kmsh1) ,2,dy(2),fu)
  call spline_fit2(spos,eot,dy,kmsh1,c1et,c2et,c3et)
!.. on the outer edge
  eot (kmsh1) = fu(0)
  c1et(kmsh1) = fu(1)
  c2et(kmsh1) = fu(2)
  c3et(kmsh1) = fu(3)
!!!   cui
  call fit_x3(spos0(1:4),cuibz(1:4),1,dy(1),fu)
  call fit_x3(spos0(nsd-3:nsd),cuibz(nsd-3:nsd),2,dy(2),fu)
  call spline_fit2(spos0,cuibz,dy,nsd,wkc10,wkc20,wkc30)
!
  call val_intp(jpl,dsl,cuibz,wkc10,wkc20,wkc30,kmsh1,nsd,cui)
!
  call fit_x3(spos(1:4),cui(1:4),1,dy(1),fu)
  call fit_x3(spos(kmsh-2:kmsh1),cui(kmsh-2:kmsh1),2,dy(2),fu)
  call spline_fit2(spos,cui,dy,kmsh1,c1i,c2i,c3i)
  cui(kmsh1) = fu(0)
  c1i(kmsh1) = fu(1)
  c2i(kmsh1) = fu(2)
  c3i(kmsh1) = fu(3)
!!!  cug
  call fit_x3(spos0(1:4),cugbz(1:4),1,dy(1),fu)
  call fit_x3(spos0(nsd-3:nsd),cugbz(nsd-3:nsd),2,dy(2),fu)
  call spline_fit2(spos0,cugbz,dy,nsd,wkc10,wkc20,wkc30)
!
  call val_intp(jpl,dsl,cugbz,wkc10,wkc20,wkc30,kmsh1,nsd,cug)
!
  call fit_x3(spos(1:4),cug(1:4),1,dy(1),fu)
  call fit_x3(spos(kmsh-2:kmsh1),cug(kmsh-2:kmsh1),2,dy(2),fu)
  call spline_fit2(spos,cug,dy,kmsh1,c1g,c2g,c3g)
  cug(kmsh1) = fu(0)
  c1g(kmsh1) = fu(1)
  c2g(kmsh1) = fu(2)
  c3g(kmsh1) = fu(3)
!!!  txi
  dy(1)=0.0d0
!!$      call fit_x3(spos(1:4),txi(1:4),1,dy(1),fu)
  call fit_x3(spos(kmsh-2:kmsh1),txi(kmsh-2:kmsh1),2,dy(2),fu)
  call spline_fit2(spos,txi,dy,kmsh1,c1ti,c2ti,c3ti)
  txi(kmsh1) = fu(0)
  c1ti(kmsh1) = fu(1)
  c2ti(kmsh1) = fu(2)
  c3ti(kmsh1) = fu(3)
!!!  txe
  dy(1)=0.0d0
!!$      call fit_x3(spos(1:4),txe(1:4),1,dy(1),fu)
  call fit_x3(spos(kmsh-2:kmsh1),txe(kmsh-2:kmsh1),2,dy(2),fu)
  call spline_fit2(spos,txe,dy,kmsh1,c1te,c2te,c3te)
  txe(kmsh1) = fu(0)
  c1te(kmsh1) = fu(1)
  c2te(kmsh1) = fu(2)
  c3te(kmsh1) = fu(3)
!!!  dxi 
  dy(1)=0.0d0
!!$      call fit_x3(spos(1:4),dxi(1:4),1,dy(1),fu)
  call fit_x3(spos(kmsh-2:kmsh1),dxi(kmsh-2:kmsh1),2,dy(2),fu)
  call spline_fit2(spos,dxi,dy,kmsh1,c1di,c2di,c3di)
  dxi(kmsh1) = fu(0)
  c1di(kmsh1) = fu(1)
  c2di(kmsh1) = fu(2)
  c3di(kmsh1) = fu(3)
!!!  dln
  dy(1)=0.0d0
!!$      call fit_x3(spos(1:4),dln(1:4),1,dy(1),fu)
  call fit_x3(spos(kmsh-2:kmsh1),dln(kmsh-2:kmsh1),2,dy(2),fu)
  call spline_fit2(spos,dln,dy,kmsh1,c1ln,c2ln,c3ln)
  dln(kmsh1) = fu(0)
  c1ln(kmsh1) = fu(1)
  c2ln(kmsh1) = fu(2)
  c3ln(kmsh1) = fu(3)
!
!
!!$  write(8,*)'### spline coefficients for flux-surface funcs. ###'
!!$  write(8,*)'### eot         c1et           c2et            c3et'
!!$  do j=1,kmsh1
!!$    write(8,7003)eot(j),c1et(j),c2et(j),c3et(j),j
!!$  end do
!!$  write(8,*)'### cui         c1i            c2i             c3i '
!!$  do j=1,kmsh1
!!$    write(8,7003)cui(j),c1i(j),c2i(j),c3i(j),j
!!$  end do
!!$  write(8,*)'### cug         c1g            c2g             c3g '
!!$  do j=1,kmsh1
!!$    write(8,7003)cug(j),c1g(j),c2g(j),c3g(j),j
!!$  end do
!!$  write(8,*)'### txi         c1ti           c2ti            c3ti'
!!$  do j=1,kmsh1
!!$    write(8,7003)txi(j),c1ti(j),c2ti(j),c3ti(j),j
!!$  end do
!!$  write(8,*)'### txe         c1te           c2te            c3te'
!!$  do j=1,kmsh1
!!$    write(8,7003)txe(j),c1te(j),c2te(j),c3te(j),j
!!$  end do
!!$  write(8,*)'### dxi         c1di           c2di            c3di'
!!$  do j=1,kmsh1
!!$    write(8,7003)dxi(j),c1di(j),c2di(j),c3di(j),j
!!$  end do
!!$  write(8,*)'### dln         c1ln           c2ln            c3ln'
!!$  do j=1,kmsh1
!!$    write(8,7003)dln(j),c1ln(j),c2ln(j),c3ln(j),j
!!$  end do
!!$!     
!!$  write(8,*)'### psib        spos           spos0           ds &
!!$     &           jp  jpl    dsl   ###'
!!$  write(8,7004)psib(1),spos(1),spos0(1),0.0d0,1,jpl(1),dsl(1)
!!$  do j=2,kmsh1
!!$    jp=jpl(j)
!!$    write(8,7004)psib(j),spos(j),spos0(jp),spos(j)-spos(j-1),j,jp,dsl(j)
!!$  end do
!
7003 format(4(1x,d14.6),i4)
7004 format(4(1x,d14.6),2i4,d14.6)
!
  deallocate(bco0,rco0,zco0,pco0,spos0)
!
  return
!

  END SUBROUTINE setfld



!
!
!********************************************************************************
SUBROUTINE val_intp(jpl,dsl,c0,c1,c2,c3,kmsh1,nsd,cval)
!********************************************************************************
  IMPLICIT REAL(KIND=8) (A-H,O-Z)

  implicit integer (i-n) ! by THW Dec 22, 2011

  integer(kind=4) :: kmsh1,nsd,jpl(kmsh1)
  real(kind=8), dimension(nsd) :: c0,c1,c2,c3
  real(kind=8), dimension(kmsh1) ::dsl,cval
! ----------------------------------------------------------------------------- !
!
  do j=1,kmsh1
    jp=jpl(j)
    ds=dsl(j)
    cval(j)=c0(jp)+(c1(jp) +(c2(jp) +c3(jp)*ds)*ds)*ds
  end do

  return

END SUBROUTINE val_intp


!
!********************************************************************************
SUBROUTINE  spline_fit2(X,Y,DY,N,C,D,E)
!********************************************************************************
  IMPLICIT REAL(KIND=8) (A-H,O-Z)

  implicit integer (i-n) ! by THW Dec 22, 2011

  real(kind=8) ::  X(N),Y(N),DY(2),C(N),D(N),E(N)
! ----------------------------------------------------------------------------- !
!
  N1 = N - 1
!!! dy/dx on edges are specified by dy
!
  d(1)=6.0d0*((y(2)-y(1))/(x(2)-x(1))-dy(1))
  d(n)=6.0d0*(dy(2)-(y(n)-y(n1))/(x(n)-x(n1)))
  c(1)=2.0d0*(x(2)-x(1)) !u(1,1)!
  do i=2,n1
    d(i)=6.0d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))
  end do
!
  do i=1,n-2
    e(i)=(x(i+1)-x(i))/c(i)   !l(i+1,i)!
    c(i+1)=2.0d0*(x(i+2)-x(i))-e(i)*(x(i+1)-x(i))   !u(i+1,i+1)!
  end do
  e(n1)=(x(n)-x(n1))/c(n1)
  c(n)=2.0d0*(x(n)-x(n1))-e(n1)*(x(n)-x(n1))   
!
  do i=2,n
    d(i)=d(i)-e(i-1)*d(i-1)
  end do
!
  d(n)=d(n)/c(n)
  do i=n-1,1,-1
    d(i)=(d(i)-(x(i+1)-x(i))*d(i+1))/c(i)
  end do
!
  do i=1,n1
    xh=(x(i+1)-x(i))
    c(i)=(y(i+1)-y(i))/xh-xh*(2.0d0*d(i)+d(i+1))/6.0d0
    e(i)=(d(i+1)-d(i))/6.0d0/xh
    d(i)=0.5d0*d(i)
  end do
  xh=(x(n)-x(n1))
  c(n)=(y(n)-y(n1))/xh+xh*(2.0d0*d(n)+d(n1))/6.0d0
  e(n)=(d(n)-d(n1))/6.0d0/xh
  d(n)=0.5d0*d(n)
!
  return

END SUBROUTINE spline_fit2



!
!********************************************************************************
SUBROUTINE fit_x3(x,y,icon,dy,fu)
!-------------------------------------------
  implicit real(kind=8) (a-h,o-z)

  implicit integer (i-n) ! by THW Dec 22, 2011

  real(kind=8) :: x(4),y(4),fu(0:3)
  integer(kind=4) :: icon
! ----------------------------------------------------------------------------- !
!
! evaluate dy/dx on x=x(1) (icon==1) or x=x(4) (icon/=1)
!
  x1=x(1)
  x2=x(2)
  x3=x(3)
  x4=x(4)
  y1=y(1)
  y2=y(2)
  y3=y(3)
  y4=y(4)
!
  a=y4/(x4-x1)/(x4-x2)/(x4-x3)
  b=y1/(x1-x2)/(x1-x3)/(x1-x4)
  c=y2/(x2-x3)/(x2-x4)/(x2-x1)
  d=y3/(x3-x4)/(x3-x1)/(x3-x2)
!
  if(icon==1) then      ! factor at x=x1
    dy=a*(x1-x2)*(x1-x3)+c*(x1-x3)*(x1-x4)+d*(x1-x4)*(x1-x2)
    dy=dy+y1*(1.0d0/(x1-x2)+1.0d0/(x1-x3)+1.0d0/(x1-x4))
    fu(0)=y1
    fu(1)=dy
    fu(2)=a*(2.0d0*x1-x2-x3)+b*(3.0d0*x1-x2-x3-x4)+c*(2.0d0*x1-x3-x4)+d*(2.0d0*x1-x4-x2)
    fu(3)=a+b+c+d
  else                  ! factor at x=x4
    dy=b*(x4-x2)*(x4-x3)+c*(x4-x3)*(x4-x1)+d*(x4-x1)*(x4-x2)
    dy=dy+y4*(1.0d0/(x4-x1)+1.0d0/(x4-x2)+1.0d0/(x4-x3))
    fu(0)=y4
    fu(1)=dy
    fu(2)=a*(3.0d0*x4-x1-x2-x3)+b*(2.0d0*x4-x2-x3)+c*(2.0d0*x4-x3-x1)+d*(2.0d0*x4-x1-x2)
    fu(3)=a+b+c+d
  end if
!... Here, fu(i) is the coefficient if y(x) is written in the forrowings
!    y(x)=fu(0)+fu(1)(x-xa)+fu(2)(x-xa)**2+fu(3)(x-xa)**3 where xa= x1 or x4.

  return

END SUBROUTINE fit_x3



!
!********************************************************************************
SUBROUTINE ntfunc0(xx,te,ti,dni)
! +---------------------------------------------------------------------------+
! |  xx = sqrt (psib/psia) : label of flux surface, psib = toroidal flux      |
! |  answer ion temp ti, electron temp. te, and ion density dni               |
! |         as functions of xx.                                               |
! |                                                                           |
! |   == satake's way ==                                                      |
! +---------------------------------------------------------------------------+
!********************************************************************************
  implicit real(KIND=8)(a-h,o-z)

  implicit integer (i-n) ! by THW Dec 22, 2011

! ----------------------------------------------------------------------------- !
!
  dni=dn0(0)*(dn0(1)+dn0(2)*dexp(-dn0(3)*xx**dn0(4)))
  ti =ti0(0)*(ti0(1)+ti0(2)*dexp(-ti0(3)*xx**ti0(4)))
  te =te0(0)*(te0(1)+te0(2)*dexp(-te0(3)*xx**te0(4)))
!
  return

END SUBROUTINE ntfunc0



!
!********************************************************************************
SUBROUTINE ntfunc1(xx,te,ti,dni)
! +---------------------------------------------------------------------------+
! |   == GSRAKE input type ==                                                 |
! +---------------------------------------------------------------------------+
!********************************************************************************
  implicit real(KIND=8)(a-h,o-z)

  implicit integer (i-n) ! by THW Dec 22, 2011

! ----------------------------------------------------------------------------- !
!
  te =0.0d0
  ti =0.0d0
  dni=0.0d0
  do i=10,1,-1
    dni=(dni+dn0(i))*xx
    ti =(ti +ti0(i))*xx
    te =(te +te0(i))*xx
  end do
!
  dni=dni+dn0(0)
  ti =ti +ti0(0)
  te =te +te0(0)
!
  return

END SUBROUTINE ntfunc1



!
!********************************************************************************
SUBROUTINE ntfunc2(xx,te,ti,dni)
! +---------------------------------------------------------------------------+
! |   == ??? input type ==                                                    |
! +---------------------------------------------------------------------------+
!********************************************************************************
  implicit real(KIND=8)(a-h,o-z)

  implicit integer (i-n) ! by THW Dec 22, 2011

  integer, parameter :: n=10
! ----------------------------------------------------------------------------- !
!
  dni=dn0(0)*exp(-dn0(2)*tanh((xx-dn0(3))/dn0(1)))
  ti =ti0(0)*exp(-ti0(2)*tanh((xx-ti0(3))/ti0(1)))
  te =te0(0)*exp(-te0(2)*tanh((xx-te0(3))/te0(1)))
!
  return

END SUBROUTINE ntfunc2


END MODULE GKV_vmecin
