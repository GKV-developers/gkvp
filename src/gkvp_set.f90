MODULE GKV_set
!-------------------------------------------------------------------------------
!
!    Set file I/O, and read parameters from namelist
!
!    Update history of gkvp_set.f90
!    --------------
!      gkvp_f0.62 (S. Maeyama, Mar 2023)
!        - Contents of subroutine set_cnfig are moved to GKV_geom, to implement
!          time-dependent metrics and operators in rotating flux-tube model.
!      gkvp_f0.61 (S. Maeyama, Mar 2021)
!        - equib_type = "s-alpha-shift" is added.
!        - Initial random phase rr is set by global (mx,gmy) indices.
!      gkvp_f0.60 (S. Maeyama, Jan 2021)
!        - Use fileio module to switch Fortran/NetCDF binary output.
!      gkvp_f0.58 (S. Maeyama, Oct 2020)
!        - init_random is added to switch random number for initialization.
!      gkvp_f0.57 (S. Maeyama, Oct 2020)
!        - Version number f0.57 is removed from filename.
!        - equib_type="slab" is added for shearless slab geometry.
!        - Set ky=0, ksq=0 for padding iend_y<my.
!
!-------------------------------------------------------------------------------

  use GKV_header
  use GKV_mpienv
  use GKV_math,   only: math_random
  use GKV_fld,    only: fld_esfield, fld_emfield_ff, fld_ff2hh
  use GKV_bndry,  only: bndry_zvm_bound_f
  use GKV_advnc,  only: caldlt_rev
  use GKV_dtc,    only: dtc_init
  use GKV_colli,  only: colli_set_param
  use GKV_colliimp,  only: colliimp_set_param
  use GKV_tips,   only: tips_reality
  !fj start 202010
  use GKV_fileio
  !fj end 202010
  use GKV_geom, only : geom_read_nml, geom_init_kxkyzvm,      &
                       geom_init_metric, geom_set_operators,  &
                       geom_reset_time

  implicit none

  private

  public   set_init, set_close

CONTAINS

!--------------------------------------
  SUBROUTINE set_init( ff, phi, Al, hh, time )
!--------------------------------------

    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1)                       :: phi, Al
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm)           :: hh
    real(kind=DP), intent(out) :: time


      call set_start
      call set_param


      if ( trim(equib_type) == "slab"     .OR. &
           trim(equib_type) == "analytic" .OR. &
           trim(equib_type) == "s-alpha"  .OR. &
           trim(equib_type) == "s-alpha-shift"  .OR. &
           trim(equib_type) == "circ-MHD" .OR. &
           trim(equib_type) == "vmec"     .OR. &
           trim(equib_type) == "eqdsk"    .OR. &
           trim(equib_type) == "ring" ) then

        call set_cnfig

      else

        if ( rankg == 0 ) then
          write(*,*) "set_cnfig_error!! on namelist: equib"
        end if
        call MPI_Finalize (ierr_mpi)
        stop

      end if

      call set_value( ff, phi, Al, hh, time )

    return


  END SUBROUTINE set_init


!--------------------------------------
  SUBROUTINE set_start
!--------------------------------------

    character(128) :: memo

    character(6)   :: crank
    character(1)   :: srank
    character(3)   :: cold, cnew

    character(10)   :: cdate, ctime

    namelist /cmemo/ memo
    namelist /calct/ calc_type, z_bound, z_filt, z_calc, art_diff, &
                     init_random, num_triad_diag
    namelist /equib/ equib_type
    namelist /run_n/ inum, ch_res
    namelist /files/ f_log, f_hst, f_phi, f_fxv, f_cnt

    character(256)   :: env_string       !fj


      call getenv ( 'fu05',env_string )  !fj
      open(inml, file=env_string )       !fj


      call date_and_time( cdate, ctime )

      read(inml,nml=cmemo)


      read(inml,nml=calct)
      read(inml,nml=equib)
      if (trim(z_calc) == "up5") art_diff = 0._DP


      inum = 1
      ch_res = .false.
      read(inml,nml=run_n)


      read(inml,nml=files)

      write( crank, fmt="(i6.6)" ) rankg
      write( srank, fmt="(i1.1)" ) ranks
      write( cold,  fmt="(i3.3)" ) inum-1
      write( cnew,  fmt="(i3.3)" ) inum


      open( olog, file=trim(f_log)//crank//"."//srank//".log."//cnew )

      if ( inum > 1 ) then
        !fj start 202010
        !open( icnt, file=trim(f_cnt)//crank//".cnt."//cold, &
        !      form="unformatted", status="old", action="read" )
        call fileio_open_icnt( trim(f_cnt) )
        !fj end 202010
      end if

      !fj start 202010
      !open( ofxv, file=trim(f_fxv)//crank//"."//srank//".fxv."//cnew, form="unformatted" )
      !open( ocnt, file=trim(f_cnt)//crank//".cnt."//cnew, form="unformatted" )
      call fileio_open_fxv( trim(f_fxv) )
      call fileio_open_cnt( trim(f_cnt) )
      !fj end 202010

      !fj start 202011
      !if ( vel_rank == 0 ) then
      !  open( omom, file=trim(f_phi)//crank//"."//srank//".mom."//cnew, form="unformatted" )
      !end if
      call fileio_open_mom( trim(f_phi) )
      !fj end 202011

      !fj start 202011
      !if ( ranks == 0 .AND. vel_rank == 0 ) then
      !  open( ophi, file=trim(f_phi)//crank//"."//srank//".phi."//cnew, form="unformatted" )
      !  open(  oAl, file=trim(f_phi)//crank//"."//srank//".Al."//cnew, form="unformatted" )
      !end if
      call fileio_open_phi( trim(f_phi) )
      call fileio_open_Al( trim(f_phi) )
      !fj end 202011

      !fj start 202011
      !if ( zsp_rank == 0 .AND. vel_rank == 0 ) then
      !  open( otrn, file=trim(f_phi)//crank//"."//srank//".trn."//cnew, form="unformatted" )
      !end if
      call fileio_open_trn( trim(f_phi) )
      !fj end 202011

      if( rankg == 0 ) then
        open( omtr, file=trim(f_hst)//"mtr."//cnew )
        open( omtf, file=trim(f_hst)//"mtf."//cnew )
        open( odtc, file=trim(f_hst)//"dtc."//cnew )
        open( oeng, file=trim(f_hst)//"eng."//cnew )
        open( omen, file=trim(f_hst)//"men."//cnew )
        open( owes, file=trim(f_hst)//"wes."//cnew )
        open( owem, file=trim(f_hst)//"wem."//cnew )
        if ( trim(calc_type) == "lin_freq" ) then
          open( ofrq, file=trim(f_hst)//"frq."//cnew )
          open( odsp, file=trim(f_hst)//"dsp."//cnew )
        end if
      end if
      if( rankg == nprocz/2 ) then
        open( ocst, file=trim(f_hst)//"cst."//cnew )
      end if

      if( rank == 0 ) then
        open( obln, file=trim(f_hst)//"bln."//srank//"."//cnew )
        open( oges, file=trim(f_hst)//"ges."//srank//"."//cnew )
        open( ogem, file=trim(f_hst)//"gem."//srank//"."//cnew )
        open( oqes, file=trim(f_hst)//"qes."//srank//"."//cnew )
        open( oqem, file=trim(f_hst)//"qem."//srank//"."//cnew )
      end if

      write( olog, * ) "##### ", trim(memo), " #####"
      write( olog, * ) ""
      write( olog, * ) "# Date : ", cdate
      write( olog, * ) "# Time : ", ctime
      write( olog, * ) ""
      write( olog, * ) "# Type of calc. : ", trim(calc_type)
      write( olog, * ) "# Boundary condition in zz       : ", trim(z_bound)
      write( olog, * ) "# 4th-order filter in zz         : ", trim(z_filt)
      write( olog, * ) "# Finite difference scheme in zz : ", trim(z_calc)
      write( olog, * ) "# Artificial diffusion in zz     : ", art_diff
      write( olog, * ) "# Number of triad transfer diag. : ", num_triad_diag
      write( olog, * ) "# Type of equib. : ", trim(equib_type)
      write( olog, * ) ""
      write( olog, * ) "# Run number = ", inum
      write( olog, * ) "# Resolution change = ", ch_res
      write( olog, * ) ""


    return
  

  END SUBROUTINE set_start


!--------------------------------------
  SUBROUTINE set_close
!--------------------------------------

     close( olog )

     !fj start 202010
     !close( icnt )
     !close( ofxv )
     !close( ocnt )
     call fileio_close_icnt
     call fileio_close_fxv
     call fileio_close_cnt
     !fj end 202010

     !fj start 202011
     !if ( vel_rank == 0 ) then
     !  close( omom )
     !end if
     call fileio_close_mom
     !fj end 202011

     !fj start 202011
     !if ( ranks == 0 .AND. vel_rank == 0 ) then
     !  close( ophi )
     !  close( oAl  )
     !end if
     call fileio_close_phi
     call fileio_close_Al
     !fj end 202011

     !fj start 202011
     !if ( zsp_rank == 0 .AND. vel_rank == 0 ) then
       !close( otrn )
     !end if
     call fileio_close_trn
     !fj end 202011

     if( rankg == 0 ) then
       close( omtr )
       close( omtf )
       close( odtc )
       close( oeng )
       close( omen )
       close( owes )
       close( owem )
       if ( trim(calc_type) == "lin_freq" ) then
         close( ofrq )
         close( odsp )
       end if
     end if
      if( rankg == nprocz/2 ) then
       close( ocst )
     end if

     if( rank == 0 ) then
       close( obln )
       close( oges )
       close( ogem )
       close( oqes )
       close( oqem )
     end if

  END SUBROUTINE set_close


!--------------------------------------
  SUBROUTINE set_param
!--------------------------------------

    namelist /runlm/ e_limit
    namelist /times/ tend, dtout_fxv, dtout_ptn, dtout_eng, dtout_dtc
    namelist /deltt/ dt_max, adapt_dt, courant_num, time_advnc


      e_limit   = 5._DP * 3600._DP - 300._DP

      read(inml,nml=runlm)


      tend      = 10.00_DP
      dtout_fxv = 5._DP
      dtout_ptn = 0.5_DP
      dtout_eng = 0.005_DP
      dtout_dtc = 0.005_DP

      read(inml,nml=times)


      dt_max    = 0.005_DP
      adapt_dt  = .false.
      courant_num  = 0.6_DP

      read(inml,nml=deltt)


        write( olog, * ) " # Numerical parameters "
        write( olog, * ) ""
        write( olog, * ) " # nxw, nyw  = ", nxw, nyw
        write( olog, * ) " # global_ny = ", global_ny
        write( olog, * ) " # global_nz = ", global_nz
        write( olog, * ) " # global_nv, global_nm = ", global_nv, global_nm
        write( olog, * ) ""
        write( olog, * ) " # nx, ny, nz   = ", nx, ny, nz
        write( olog, * ) " # nv, nm       = ", nv, nm
        write( olog, * ) " # nzb, nvb     = ", nzb, nvb
        write( olog, * ) " # nxw_size (local xw allocation size) = ", nxw_size
        write( olog, * ) " # number of species  = ", nprocs
        write( olog, * ) " # ranks=0: Electron"
        write( olog, * ) " # ranks=1: main ion"
        write( olog, * ) " # ranks>1: other ions"
        write( olog, * ) " # Note that proton mass mp and main ion tmep. Ti is used for normalizations"
        write( olog, * ) " # kx, ky are normalized with rho_tp = mp*vtp/e/Baxi, where vtp = sqrt(Ti/mp)"
        write( olog, * ) " # time t is normalized with Raxi/vtp"
        write( olog, * ) ""
        write( olog, * ) " # e_limit      = ", e_limit
        write( olog, * ) " # tend         = ", tend
        write( olog, * ) " # dtout_fxv, dtout_ptn = ", dtout_fxv, dtout_ptn
        write( olog, * ) " # dtout_eng, dtout_dtc = ", dtout_eng, dtout_dtc
        write( olog, * ) ""
        write( olog, * ) " # maximum time step dt_max = ", dt_max
        write( olog, * ) " # adaptive time step = ", adapt_dt
        write( olog, * ) ""

        write( olog, * ) " # MPI parallelization parameters "
        write( olog, * ) ""
        write( olog, * ) " # nproc , rankg = ", nproc , rankg
        write( olog, * ) " # nprocw, rankw = ", nprocw, rankw
        write( olog, * ) " # nprocz, rankz = ", nprocz, rankz
        write( olog, * ) " # nprocv, rankv = ", nprocv, rankv
        write( olog, * ) " # nprocm, rankm = ", nprocm, rankm
        write( olog, * ) " # nprocs, rank  = ", nprocs , rank
        write( olog, * ) " # izup, izdn    = ", izup, izdn
        write( olog, * ) " # ivup, ivdn    = ", ivup, ivdn
        write( olog, * ) " # imup, imdn    = ", imup, imdn
        write( olog, * ) ""
        write( olog, * ) " # fft_nproc , fft_rank  = ", fft_nproc , fft_rank
        write( olog, * ) " # zsp_nproc , zsp_rank  = ", zsp_nproc , zsp_rank
        write( olog, * ) " # vel_nproc , vel_rank  = ", vel_nproc , vel_rank
        write( olog, * ) " # ranks                 = ", ranks
        write( olog, * ) ""

        write( olog, * ) " # ist_y     = ", ist_y
        write( olog, * ) " # iend_y    = ", iend_y
        write( olog, * ) " # nsize_y   = ", nsize_y
        write( olog, * ) " # ist1_y    = ", ist1_y
        write( olog, * ) " # ist_y_g   = ", ist_y_g
        write( olog, * ) " # iend_y_g  = ", iend_y_g
        write( olog, * ) ""
        write( olog, * ) " # ist_xw    = ", ist_xw
        write( olog, * ) " # iend_xw   = ", iend_xw
        write( olog, * ) " # nsize_xw  = ", nsize_xw
        write( olog, * ) " # ist_xw_g  = ", ist_xw_g
        write( olog, * ) " # iend_xw_g = ", iend_xw_g
        write( olog, * ) ""


  END SUBROUTINE set_param


!--------------------------------------
  SUBROUTINE set_cnfig
!--------------------------------------

    real(kind=DP), dimension(0:ns-1,0:ns-1) :: nust
    real(kind=DP) :: lx, ly, eps_r
    integer       :: is1, is2
   
    namelist /nu_ref/ Nref,     & ! reference (electron) density in m^(-3)
                      Lref,     & ! reference length (=R_axis) in m
                      Tref,     & ! reference main-ion (ranks=1) temperature in keV 
                      col_type, & ! flag for collision type: LB or full
                      iFLR,     & ! flag for GK- or DK-limit in collision
                      icheck      ! flag for Maxwellain anihilation test (w/ iFLR=0)

! --- read GKV namelist relating to configurations ---
        call geom_read_nml

! --- coordinate settings (time-indep.) ---
        call geom_init_kxkyzvm(lx, ly, eps_r)

! --- set collision frequencies and v-space functions for multi-species GK collision
        read(inml,nml=nu_ref)
        call colli_set_param(q_0, eps_r, nust)
        write( olog, * ) " # Collision parameters"
        write( olog, * ) ""
        write( olog, * ) " # Nref [m^-3]  = ", Nref
        write( olog, * ) " # Lref [m]     = ", Lref
        write( olog, * ) " # Tref [keV]   = ", Tref
        write( olog, * ) " # col_type     = ", col_type
        write( olog, * ) " # iFLR         = ", iFLR
        write( olog, * ) " # icheck       = ", icheck
        write( olog, * ) 
        write( olog, * ) " # Normalized collisionality: nu*"
        do is1 = 0, ns-1
        do is2 = 0, ns-1
        write( olog, * ) " # a, b, nu*_ab = ", is1, is2, nust(is1,is2)
        end do
        end do
        write( olog, * ) 
        if ( trim(col_type) == "LB" ) then
          write( olog, * ) " # col.-freq. bias factor for LB, nu = ", nu(:)
          write( olog, * ) 
        end if

        Zeff = 0._DP
        do is1 = 1, ns-1
          Zeff = Zeff + fcs(is1)*Znum(is1)
        end do
        write( olog, * ) " # Zeff         = ", Zeff
        write( olog, * ) 


! --- coordinate settings (explicitly time-dependent metrics) ---
        call geom_init_metric

! --- operator settings (time-dependent through metrics) ---
        call geom_set_operators
        if (trim(col_type) == "full" .or. trim(col_type) == "lorentz" .or. trim(time_advnc) == "imp_colli") then
          call colliimp_set_param
        end if

! --- initial estimate of time steps ---
        call dtc_init( lx, ly, vmax )

  END SUBROUTINE set_cnfig


!--------------------------------------
  SUBROUTINE set_value( ff, phi, Al, hh, time )
!--------------------------------------

    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1)                       :: phi, Al
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm)           :: hh

    real(kind=DP), intent(out) :: time


    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: wf, dh, cf, ef
    real(kind=DP),    dimension(:),         allocatable :: rr
    character(15) :: colliflag
    integer :: input_status
    integer :: mx, my, iz, iv, im, nx_init, gmy


      ff(:,:,:,:,:) = ( 0._DP, 0._DP )
      phi(:,:,:)    = ( 0._DP, 0._DP )
      Al(:,:,:)     = ( 0._DP, 0._DP )
      hh(:,:,:,:,:) = ( 0._DP, 0._DP )


      if( inum == 1 ) then


        time     = 0._DP

        allocate( rr((2*nx+1)*(global_ny+1)) )

        if (init_random) then
          call math_random ( rr )
        else
          rr(:) = 0._DP
        end if

        if ( nx0 > nx ) then 
          nx_init = nx
        else        
          nx_init = nx0
        end if

!!!        my = 0 for R-H test
        do im = 0, nm
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              do my = ist1_y, iend_y
            ! do my = ist_y, iend_y
                gmy = my + (ny+1)*rankw
                do mx = -nx_init, nx_init
                  ff(mx,my,iz,iv,im)   = dns1(ranks) * fmx(iz,iv,im)  &
                          * ( 1._DP + vl(iv) + zz(iz) )**2            &
                          * exp( -zz(iz)**2 / (0.2_DP*pi)**2 ) &
                          * exp( ui * twopi * rr(mx+nx+1+(2*nx+1)*gmy) )
                end do
              end do
            end do
          end do
        end do

        if ( rankw == 0 ) then
          do im = 0, nm
            do iv = 1, 2*nv
              do iz = -nz, nz-1
                my = 0
                  do mx = 1, nx
                    ff(-mx,my,iz,iv,im) = conjg( ff(mx,my,iz,iv,im) )
                  end do
                ff(0,0,iz,iv,im) = ( 0._DP, 0._DP )
              end do
            end do
          end do
        end if

        deallocate( rr )


        if ( ch_res ) call set_ch_resolution ( ff, time )


      else


        allocate( wf(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) )

        time   = - 1._DP

        do 
          !fj start 202010
          !read( unit=icnt, iostat=input_status ) time, wf
          call fileio_read_cnt( wf, time, input_status )
          !fj end 202010

          if ( input_status < 0 ) then
            write( olog, * ) &
               " # end of file of unit=30 is detected --> stop"
            call flush(olog)
            call MPI_Abort(MPI_COMM_WORLD, ierr_mpi)
            stop
          end if

          if ( input_status > 0 ) then
            write( olog, * ) &
               " # input error of unit=30 is detected --> stop"
            call flush(olog)
            call MPI_Abort(MPI_COMM_WORLD, ierr_mpi)
            stop
          end if

          if ( time > eps ) exit
        end do


        do im = 0, nm
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  ff(mx,my,iz,iv,im) = wf(mx,my,iz,iv,im)
                end do
              end do
            end do
          end do
        end do

        deallocate( wf )


          write( olog, * ) ""
          write( olog, * ) " # simulation is re-started at t = ", time
          write( olog, * ) ""


      end if

     !%%% For shearflow rotating flux tube model %%%
      if (gamma_e /= 0._DP .and. trim(flag_shearflow) =="rotating") then
        call geom_reset_time(time)
        if (trim(col_type) == "full" .or. trim(col_type) == "lorentz" .or. trim(time_advnc) == "imp_colli") then
          call colliimp_set_param
        end if
      end if
     !%%%

      call bndry_zvm_bound_f( ff )

      call fld_esfield( ff, phi )
      if ( beta .ne. 0._DP ) then
        call fld_emfield_ff( ff, Al )
      end if
      call fld_ff2hh( ff, Al, hh )

      call tips_reality( hh )

      allocate( dh(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) )
      allocate( cf(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) )
      allocate( ef(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) )

      colliflag = "collisional"
      call caldlt_rev( colliflag, ff, phi, Al, hh, dh, cf, ef )

      deallocate( dh )
      deallocate( cf )
      deallocate( ef )

!! --- for debug
!      call MPI_Finalize(ierr_mpi)
!      stop


  END SUBROUTINE set_value


!--------------------------------------
  SUBROUTINE set_ch_resolution ( ff, time )
!--------------------------------------

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    real(kind=DP), intent(out) :: time


    !--- Set perpendicular resolution employed in the input file. ---
    !!! NOTE !!!
    !    Resolutions in (z,v,m,s) should be kept the same.
    !    Since lx and ly should be larger than or equal to lx0, ly0,
    !    xfold and yfold are set to be integers.
    !!!!!!!!!!!!
      integer, parameter :: xfold = 1   ! kxmin0 = xfold * kxmin
      integer, parameter :: yfold = 1   ! kymin0 = yfold * kymin
      integer, parameter :: nx0 = 95, global_ny0 = 1, nprocw0 = 1
      integer, parameter :: ny0 = global_ny0 / nprocw0
      real(kind=DP) :: amplify = 1._DP
    !----------------------------------------------------------------


    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: wf
    integer :: input_status
    integer :: ny_size0, nwk0, irw, ist_y_g0, iend_y_g0
    integer :: mx, my, iz, iv, im, mxw, myw
    character(6)   :: crank
    character(1)   :: srank



      allocate( wf(-nx0:nx0,0:ny0,-nz:nz-1,1:2*nv,0:nm) )

      do irw = 0, nprocw0-1

        ny_size0 = global_ny0 + 1 
        if( mod(ny_size0,nprocw0) == 0 )  then
          nwk0    = ny_size0 / nprocw0
        else
          nwk0    = ny_size0 / nprocw0 + 1
        endif
        !--- global index range ---------------- 
        ist_y_g0  = nwk0*irw
        iend_y_g0 = min( nwk0*(irw+1)-1, (ny_size0-1) )

        write( crank, fmt="(i6.6)" ) irw + nprocw0*rankz + nprocw0*nprocz*rankv  &
            + nprocw0*nprocz*nprocv*rankm + nprocw0*nprocz*nprocv*nprocm*ranks
        write( srank, fmt="(i1.1)" ) ranks

        !fj start 202010
        !open( icnt, file=trim(f_cnt)//crank//".cnt.000", &
        !      form="unformatted", status="old", action="read" )
        !read( unit=icnt, iostat=input_status ) time, wf
        !rewind( icnt )
        !close( icnt )
        call fileio_open_icnt( trim(f_cnt) )
        call fileio_read_cnt( wf, time, input_status )
        call fileio_close_icnt
        !fj end 202010

        if ( ist_y_g <= iend_y_g0 * yfold .and. ist_y_g0 * yfold <= iend_y_g ) then
          do im = 0, nm
            do iv = 1, 2*nv
              do iz = -nz, nz-1
                do my = ist_y_g0, iend_y_g0
                  myw = my * yfold
                  if ( ist_y_g <= myw .and. myw <= iend_y_g ) then
                    do mx = -nx0, nx0
                      mxw = mx * xfold
                      if ( -nx <= mxw .and. mxw <= nx ) then
                        ff(mxw,myw-ist_y_g,iz,iv,im) = wf(mx,my-ist_y_g0,iz,iv,im) * amplify
                      end if
                    end do
                  end if
                end do
              end do
            end do
          end do
        end if

      end do

      deallocate( wf )

      write( olog, * ) ""
      write( olog, * ) " # simulation is re-started at t = ", time
      write( olog, * ) " # perpendicular resolutions are changed:"
      write( olog, * ) " # lx is ", xfold, "times larger" 
      write( olog, * ) " # ly is ", yfold, "times larger" 
      write( olog, * ) " # nx =", nx0, "to", nx
      write( olog, * ) " # global_ny =", global_ny0, "to", global_ny
      write( olog, * ) " # nprocw from", nprocw0, "to", nprocw
      write( olog, * ) ""


  END SUBROUTINE set_ch_resolution


END MODULE GKV_set
