MODULE GKV_vmecbzx
!-------------------------------------------------------------------------------
!
!    Calculate the magnetic field components and metric coefficients
!      from the VMEC equilibrium by using BZX code. 
!
!    Update history of gkvp_vmecbxz.f90
!    --------------
!      gkvp_f0.62 (S. Maeyama, Mar 2023)
!        - Input of vmecbzx_boozxcoef is modified from local iz for each rankz
!          to global index giz.
!      gkvp_f0.57 (S. Maeyama, Oct 2020)
!        - Version number f0.57 is removed from filename.
!
!-------------------------------------------------------------------------------

  use GKV_header
  use GKV_mpienv

  implicit none

  private

  public :: vmecbzx_boozx_read, vmecbzx_boozx_coeff

  real(kind=DP), dimension(:),         allocatable   :: ss_bz, q_bz, shat_bz, eps_bz
  real(kind=DP), dimension(:),         allocatable   :: theta_bz, zeta_bz
  real(kind=DP), dimension(:,:,:,:,:), allocatable   :: ggup_bz
  real(kind=DP), dimension(:,:,:),     allocatable   :: B_bz, rootg_bz, rootg_bz0
  real(kind=DP), dimension(:,:,:),     allocatable   :: dBds_bz, dBdt_bz, dBdz_bz
  real(kind=DP), dimension(:,:),       allocatable   :: bbozc_bz, bbozs_bz
  real(kind=DP), dimension(:,:,:),     allocatable   :: rr_bz, zz_bz, ph_bz
  integer,       dimension(:),         allocatable   :: ixn_bz, ixm_bz

  integer        :: nss_bz, ntheta_bz, nzeta_bz
  integer        :: nfp_bz, mnboz_bz, mboz_bz, nboz_bz, asym_flg
  real(kind=DP)  :: Rax_bz, Bax_bz, aa_bz, volume_bz, alpha_fix

CONTAINS


!--------------------------------------
  SUBROUTINE vmecbzx_boozx_read( nss, ntheta, nzeta )
!--------------------------------------

    implicit none

    integer, intent(in)        :: nss, ntheta, nzeta
    integer                    :: ibzx
    character(512) :: f_bozx
!   character(512) :: env_string       !fj

    namelist /bozxf/ f_bozx

    allocate(q_bz(1:nss),shat_bz(1:nss),eps_bz(1:nss))
    allocate(ss_bz(1:nss),theta_bz(0:ntheta),zeta_bz(0:nzeta))
    allocate(ggup_bz(1:nss,0:ntheta,0:nzeta,1:3,1:3))
    allocate(B_bz(1:nss,0:ntheta,0:nzeta),rootg_bz(1:nss,0:ntheta,0:nzeta),rootg_bz0(1:nss,0:ntheta,0:nzeta))
    allocate(dBds_bz(1:nss,0:ntheta,0:nzeta),dBdt_bz(1:nss,0:ntheta,0:nzeta),dBdz_bz(1:nss,0:ntheta,0:nzeta))

      read(inml,nml=bozxf)
      ibzx = 5000

!      call getenv ( 'fu09',env_string )  !fj
!      open(ibzx, file=env_string, status="old", action="read", form="unformatted", CONVERT='LITTLE_ENDIAN' )          !fj


      open( ibzx, file=trim(f_bozx)//"metric_boozer.bin.dat", status="old", &
                  action="read", form="unformatted", CONVERT='LITTLE_ENDIAN' )
      write(olog,*) "# mag.coord.(Booz_xform) input file : "
      write(olog,*) trim(f_bozx)//"metric_boozer.bin.dat"

! --- read B-field and metric components
      ! binary
      read(unit=ibzx) nfp_bz, nss_bz, ntheta_bz, nzeta_bz, mnboz_bz, mboz_bz, nboz_bz, & 
                      Rax_bz, Bax_bz, aa_bz, volume_bz, asym_flg, alpha_fix
      if (nss_bz /= nss .OR. ntheta_bz /= ntheta .OR. nzeta_bz /= nzeta ) then 
        print*, "nss_bz, ntheta_bz, nzeta_bz = ", nss_bz, ntheta_bz, nzeta_bz
        call MPI_Finalize
        stop "# nss or ntheta or nzeta is inconsistent between gkv_namelist and metric_boozer.dat! --> stop"
      end if 

      read(unit=ibzx) ss_bz, theta_bz, zeta_bz, q_bz, shat_bz, eps_bz, B_bz, rootg_bz, rootg_bz0, & 
                      ggup_bz, dBds_bz, dBdt_bz, dBdz_bz

      allocate(rr_bz(1:nss,0:ntheta,0:nzeta),zz_bz(1:nss,0:ntheta,0:nzeta),ph_bz(1:nss,0:ntheta,0:nzeta))
      allocate(bbozc_bz(1:mnboz_bz,1:nss),bbozs_bz(1:mnboz_bz,1:nss))
      allocate(ixn_bz(1:mnboz_bz),ixm_bz(1:mnboz_bz))

      read(unit=ibzx) rr_bz, zz_bz, ph_bz, bbozc_bz, ixn_bz, ixm_bz


                  !!! ! ascii for debug: NOT Used usually
                  !!! read(unit=ibzx, fmt="(7I5,4ES24.15e3)") nfp_bz, nss, ntheta, nzeta, mnboz_bz, mboz_bz, nboz_bz, &
                  !!!                                         Rax_bz, Bax_bz, aa_bz, volume_bz
                  !!! do js = 1, nss
                  !!!   do iz = 0, nzeta
                  !!!     do it = 0, ntheta
                  !!!       write(unit=ibzx, fmt="(256ES24.15e3)")                       &
                  !!!               ss_bz(js),            theta_bz(it),    zeta_bz(iz),  &
                  !!!                q_bz(js),             shat_bz(js),     eps_bz(js),  &
                  !!!          B_bz(js,it,iz),      rootg_bz(js,it,iz),                  &
                  !!!         ggup_bz(js,it,iz,1,1), ggup_bz(js,it,iz,1,2), ggup_bz(js,it,iz,1,3),    &
                  !!!         ggup_bz(js,it,iz,2,2), ggup_bz(js,it,iz,2,3), ggup_bz(js,it,iz,3,3),    & 
                  !!!         Rax_bz, Bax_bz, aa_bz, volume_bz, nss, ntheta, nzeta, nfp_bz
                  !!!     end do
                  !!!     write(unit=ibzx, fmt=*)
                  !!!   end do
                  !!!     write(unit=ibzx, fmt=*)
                  !!!     write(unit=ibzx, fmt=*)
                  !!! end do

      write(olog,*) "# Rax_bz[m], Bax_bz[T], aa_bz[m], volume_bz[m3] = ", Rax_bz, Bax_bz, aa_bz, volume_bz
      write(olog,*) "# nss(=ns_b), ntheta(=ntht), nzeta, asym_flg    = ", nss, ntheta, nzeta, asym_flg
      write(olog,*) "# nfp_bz, mboz_bz, nboz_bz, mnboz_bz            = ", nfp_bz, mboz_bz, nboz_bz, mnboz_bz
      if(nzeta == 0) write(olog,*) "# fixed alpha (= zeta - q*theta) = ", alpha_fix
      write(olog,*) 


! --- normalizartion with Bax, Rax 
      
      B_bz       =  B_bz/Bax_bz
      rootg_bz   =  rootg_bz/Rax_bz**3
      ggup_bz    =  ggup_bz*Rax_bz**2
      dBds_bz    =  dBds_bz/Bax_bz
      dBdt_bz    =  dBdt_bz/Bax_bz
      dBdz_bz    =  dBdz_bz/Bax_bz
      bbozc_bz   =  bbozc_bz/Bax_bz


      return 

  END SUBROUTINE vmecbzx_boozx_read



!----------------------------------------------------------------------------------
!smae start 202303
!  SUBROUTINE vmecbzx_boozx_coeff( isw,  nss,  ntheta,  nzeta,  s_input,  iz, zz,  lz_l,    &  ! input 
  SUBROUTINE vmecbzx_boozx_coeff( isw,  nss,  ntheta,  nzeta,  s_input,  giz, zz,  lz_l,    &  ! input 
!smae end 202303
                                   s_0,       q_0,    s_hat,   eps_r,   phi_ax,           &  ! output
                                   omg,     rootg,   domgdx,  domgdz,  domgdy,            &
                                  gg11,      gg12,     gg13,    gg22,                     &
                                  gg23,      gg33  )
!----------------------------------------------------------------------------------

!smae start 202303
!    integer, intent(in)        :: isw, nss, ntheta, nzeta, iz
    integer, intent(in)        :: isw, nss, ntheta, nzeta, giz
!smae end 202303
    real(kind=DP), intent(in)  :: s_input, zz, lz_l

    real(kind=DP), intent(inout) :: s_0, q_0, s_hat, eps_r, phi_ax
    real(kind=DP), intent(out) :: omg, rootg, domgdx, domgdz, domgdy
    real(kind=DP), intent(out) :: gg11, gg12, gg13, gg22, gg23, gg33

! --- local variables 
    integer                    :: is0, jj0, zt0
    real(kind=DP)              :: eps_a


   ! is0 = nint(s_input*(nss+1))
    is0 = nint(s_input*(nss-1))+1
    
    if ( isw == 0 ) then 

      s_0   =   ss_bz(is0)
      q_0   =    q_bz(is0)
      s_hat = shat_bz(is0)
      eps_r =  eps_bz(is0)

    else if ( isw == 1 ) then

      s_0   =   ss_bz(is0)
!smae start 202303
!      giz = (-global_nz + 2*nz*rankz + iz + nz )          
!smae end 202303
      jj0 = giz + ntheta/2
  
     
! --- for debug
!      write(2000+rankg,*) giz, jj0

      eps_a =  eps_bz(nss)

      zt0    =  0  ! for fixed alpha (i.e., single fluxtube)
      phi_ax =  ph_bz(is0,jj0,zt0)
      omg    =  B_bz(is0,jj0,zt0)
      rootg  =  (eps_a**(-2))*q_0/s_0*rootg_bz(is0,jj0,zt0)
      domgdx =  ( dBds_bz(is0,jj0,zt0) + q_0*s_hat*zz/s_0*dBdz_bz(is0,jj0,zt0) )/eps_a
      domgdz =    dBdt_bz(is0,jj0,zt0) + q_0*dBdz_bz(is0,jj0,zt0)
      domgdy =  - q_0/s_0*dBdz_bz(is0,jj0,zt0)/eps_a

      gg11   =  eps_a**2 * ggup_bz(is0,jj0,zt0,1,1)

      gg12   =  eps_a**2 * (  s_hat*zz*ggup_bz(is0,jj0,zt0,1,1)   &
                             +     s_0*ggup_bz(is0,jj0,zt0,1,2)   &
                             - s_0/q_0*ggup_bz(is0,jj0,zt0,1,3) )

      gg13   =  eps_a*ggup_bz(is0,jj0,zt0,1,2)

      gg22   =  eps_a**2 * (            (s_hat*zz)**2*ggup_bz(is0,jj0,zt0,1,1)   &
                             +                 s_0**2*ggup_bz(is0,jj0,zt0,2,2)   &
                             +           (s_0/q_0)**2*ggup_bz(is0,jj0,zt0,3,3)   & 
                             +     2._DP*s_0*s_hat*zz*ggup_bz(is0,jj0,zt0,1,2)   & 
                             -       2._DP*s_0**2/q_0*ggup_bz(is0,jj0,zt0,2,3)   & 
                             - 2._DP*s_0/q_0*s_hat*zz*ggup_bz(is0,jj0,zt0,1,3) )

      gg23   =     eps_a * (  s_hat*zz*ggup_bz(is0,jj0,zt0,1,2)   &
                             +     s_0*ggup_bz(is0,jj0,zt0,2,2)   &
                             - s_0/q_0*ggup_bz(is0,jj0,zt0,2,3) )

      gg33   =  ggup_bz(is0,jj0,zt0,2,2)

    end if

    return
    
  END SUBROUTINE vmecbzx_boozx_coeff


END MODULE GKV_vmecbzx
