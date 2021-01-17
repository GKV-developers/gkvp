MODULE GKV_igs
!-------------------------------------------------------------------------------
!
!    Calculate the magnetic field components and metric coefficients
!      from the MEUDAS or G-EQDSK equilibrium by using IGS code. 
!
!    Notes
!    -----
!      1D and 2D equilibrium profiles, e.g., q(Psi), Psi(R,Z), 
!        from a free-boundary 2D Grad-shafranov solver MEUDAS or from G-EQDSK 
!          are converted to SFL magnetic coordinates (Axisym/Boozer/Hamada) 
!            by the interface code IGS developed by A. Matsuyama and M. Nakata. 
!
!    Update history of gkvp_igs.f90
!    --------------
!      gkvp_f0.57 (S. Maeyama, Oct 2020)
!        - Version number f0.57 is removed from filename.
!      gkvp_r1.3 (M. Nakata, May 2013)
!        - First implementation.
!
!-------------------------------------------------------------------------------

  use GKV_header
  use GKV_mpienv

  implicit none

  private

  public :: igs_read, igs_coeff

  real(kind=DP), dimension(:),       allocatable   :: ss_mc, q_mc, shat_mc, eps_mc, bsq_mc
  real(kind=DP), dimension(:),       allocatable   :: theta_mc
  real(kind=DP), dimension(:,:,:,:), allocatable   :: ggup_mc, ggdn_mc
  real(kind=DP), dimension(:,:),     allocatable   :: Bupt_mc, Bupz_mc, Bdns_mc, Bdnt_mc, Bdnz_mc
  real(kind=DP), dimension(:,:),     allocatable   :: dBdt_mc, dBds_mc
  real(kind=DP), dimension(:,:),     allocatable   :: B_mc, rootg_mc
  real(kind=DP), dimension(:,:),     allocatable   :: real2axi_mc, axi2mag_mc

CONTAINS


!--------------------------------------
  SUBROUTINE igs_read(mc_type, nss, ntheta)
!--------------------------------------

    implicit none

    integer, intent(in)        :: mc_type, nss, ntheta
    integer                    :: is, jj, ierr, imds
    character(512) :: f_igs

    namelist /igsf/ f_igs

    allocate(ss_mc(1:nss),q_mc(1:nss),shat_mc(1:nss),eps_mc(1:nss),bsq_mc(1:nss))
    allocate(theta_mc(1:ntheta))
    allocate(ggup_mc(1:ntheta,1:nss,1:3,1:3),ggdn_mc(1:ntheta,1:nss,1:3,1:3))
    allocate(Bupt_mc(1:ntheta,1:nss),Bupz_mc(1:ntheta,1:nss))
    allocate(Bdns_mc(1:ntheta,1:nss),Bdnt_mc(1:ntheta,1:nss),Bdnz_mc(1:ntheta,1:nss))
    allocate(dBdt_mc(1:ntheta,1:nss),dBds_mc(1:ntheta,1:nss))
    allocate(B_mc(1:ntheta,1:nss),rootg_mc(1:ntheta,1:nss))
    allocate(real2axi_mc(1:ntheta,1:nss),axi2mag_mc(1:ntheta,1:nss))

      read(inml,nml=igsf)
      imds = 5000
    
      if ( mc_type == 0 ) then

        open( imds, file=trim(f_igs)//"METRIC_axi.OUT", status="old", action="read" )
        write(olog,*) "# mag.coord.(IGS) input files : "
        write(olog,*) trim(f_igs)//"METRIC_axi.OUT"
        write( olog, * ) 

      else if ( mc_type == 1 ) then

        open( imds, file=trim(f_igs)//"METRIC_boz.OUT", status="old", action="read" )
        write(olog,*) " # mag.coord.(IGS) input files : "
        write(olog,*) trim(f_igs)//"METRIC_boz.OUT"
        write( olog, * ) 

      else if ( mc_type == 2 ) then

        open( imds, file=trim(f_igs)//"METRIC_ham.OUT", status="old", action="read" )
        write(olog,*) " # mag.coord.(IGS) input files : "
        write(olog,*) trim(f_igs)//"METRIC_ham.OUT"
        write( olog, * ) 

      else 

        write(*,*) "## Invalid mc_type setting!!" 
        call MPI_Finalize(ierr)
        stop 

      end if  

! --- read B-field and metric components
      do is = 1, nss
        do jj = 1, ntheta
          read(imds,fmt="(2f15.8, SP, 27ES24.15e3)")   ss_mc(is), theta_mc(jj),   &
                    q_mc(is),           shat_mc(is),            eps_mc(is),   &
                 B_mc(jj,is),            bsq_mc(is),       rootg_mc(jj,is),   &
          ggdn_mc(jj,is,1,1),    ggdn_mc(jj,is,1,2),    ggdn_mc(jj,is,1,3),   &
          ggdn_mc(jj,is,2,2),    ggdn_mc(jj,is,2,3),    ggdn_mc(jj,is,3,3),   &
          ggup_mc(jj,is,1,1),    ggup_mc(jj,is,1,2),    ggup_mc(jj,is,1,3),   &
          ggup_mc(jj,is,2,2),    ggup_mc(jj,is,2,3),    ggup_mc(jj,is,3,3),   &
              Bupt_mc(jj,is),        Bupz_mc(jj,is),        Bdnt_mc(jj,is),   &  
              Bdnz_mc(jj,is),        Bdns_mc(jj,is),        dBds_mc(jj,is),   &
              dBdt_mc(jj,is),    real2axi_mc(jj,is),     axi2mag_mc(jj,is)
        end do
        read(imds,fmt=*) 
      end do

      return 

  END SUBROUTINE igs_read


!----------------------------------------------------------------------------------
  SUBROUTINE igs_coeff( isw,  mc_type,  nss,  ntheta,  s_input,  zz,  lz_l, &  ! input 
                                   s_0,       q_0,    s_hat,   eps_r,   theta,   &  ! output
                                   omg,     rootg,   domgdx,  domgdz,  domgdy,   &
                                  gg11,      gg12,     gg13,    gg22,            &
                                  gg23,      gg33  )
!----------------------------------------------------------------------------------

    integer, intent(in)        :: isw, mc_type, nss, ntheta
    real(kind=DP), intent(in)  :: s_input, zz, lz_l

    real(kind=DP), intent(inout) :: s_0, q_0, s_hat, eps_r, theta
    real(kind=DP), intent(out) :: omg, rootg, domgdx, domgdz, domgdy
    real(kind=DP), intent(out) :: gg11, gg12, gg13, gg22, gg23, gg33

! --- local variables 
    integer                    :: is0, iz0, nz0, jj0
    real(kind=DP)              :: eps_a


    is0 = nint(s_input*(nss-1))+1
    
    if (mc_type == 0 ) then 
      axi2mag_mc(:,:) = 0._DP 
    end if

    if ( isw == 0 ) then 

      s_0   =   ss_mc(is0)
      q_0   =    q_mc(is0)
      s_hat = shat_mc(is0)
      eps_r =  eps_mc(is0)

    else if ( isw == 1 ) then

      s_0   =   ss_mc(is0)

      iz0 = nint(nz*zz/lz_l)  ! get global_iz
      nz0 = nint(pi*nz/lz_l)        ! get grid number on z=[0,pi)

      if ( mod(iz0,nz0) /= 0 .AND. mod(iz0/nz0,2) == 0 ) then 
   
        jj0 = mod(iz0,nz0) + (nz0+1)

      else if ( mod(iz0,nz0) /= 0 .AND. mod(iz0/nz0,2) /= 0 ) then    

        jj0 = mod(iz0,nz0) - (mod(iz0,nz0)/abs(mod(iz0,nz0)))*nz0 + (nz0+1)

      else if ( mod(iz0,nz0) == 0 .AND. mod(iz0/nz0,2) == 0 ) then    

        jj0 = 0 + (nz0+1)

      else if ( mod(iz0,nz0) == 0 .AND. mod(iz0/nz0,2) /= 0 ) then    

        jj0 = 1

      end if
              
! --- for debug
!      write(2000+rankg,*) iz0, jj0, nz0, zz

      eps_a =  eps_mc(nss)
     
      theta  =  zz - real2axi_mc(jj0,is0) - axi2mag_mc(jj0,is0)
      omg    =  B_mc(jj0,is0)
      rootg  =  (eps_a**(-2))*q_0/s_0*rootg_mc(jj0,is0)
      domgdx =  dBds_mc(jj0,is0)/eps_a
      domgdz =  dBdt_mc(jj0,is0)
      domgdy =  0._DP

      gg11   =  eps_a**2 * ggup_mc(jj0,is0,1,1)

      gg12   =  eps_a**2 * (  s_hat*zz*ggup_mc(jj0,is0,1,1)   &
                             +     s_0*ggup_mc(jj0,is0,1,2)   &
                             - s_0/q_0*ggup_mc(jj0,is0,1,3) )

      gg13   =  eps_a*ggup_mc(jj0,is0,1,2)

      gg22   =  eps_a**2 * (            (s_hat*zz)**2*ggup_mc(jj0,is0,1,1)   &
                             +                 s_0**2*ggup_mc(jj0,is0,2,2)   &
                             +           (s_0/q_0)**2*ggup_mc(jj0,is0,3,3)   & 
                             +     2._DP*s_0*s_hat*zz*ggup_mc(jj0,is0,1,2)   & 
                             -       2._DP*s_0**2/q_0*ggup_mc(jj0,is0,2,3)   & 
                             - 2._DP*s_0/q_0*s_hat*zz*ggup_mc(jj0,is0,1,3) )

      gg23   =     eps_a * (  s_hat*zz*ggup_mc(jj0,is0,1,2)   &
                             +     s_0*ggup_mc(jj0,is0,2,2)   &
                             - s_0/q_0*ggup_mc(jj0,is0,2,3) )

      gg33   =  ggup_mc(jj0,is0,2,2)

    end if

    return
    
  END SUBROUTINE igs_coeff


END MODULE GKV_igs
