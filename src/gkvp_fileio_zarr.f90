MODULE GKV_fileio
!-------------------------------------------------------------------------------
!
!    File I/O interface for Fortran binary output
!
!    Update history of gkvp_fileio_zarr.f90
!    --------------
!      gkvp_f0.64 (S. Maeyama, June 2025)
!        - Zarr store I/O interface is implemented.
!
!-------------------------------------------------------------------------------

  use GKV_header
  use GKV_mpienv

  implicit  none

  private

  character(1024), save :: path_icnt, path_cnt, path_fxv, &
                           path_phi, path_Al, path_mom,   &
                           path_trn, path_tri
  integer, save :: nt_out_fxv, nt_out_phi, nt_out_Al, nt_out_mom, &
                   nt_out_trn, nt_out_tri
  integer, parameter :: nmom = 6   ! Number of output moments
  integer, parameter :: ntrn = 12  ! Number of output total transfer diagnostics
  integer, parameter :: ntri = 6   ! Number of output triad transfer diagnostics

  public    fileio_open_icnt,fileio_close_icnt, &
            fileio_open_cnt, fileio_close_cnt, &
            fileio_open_fxv, fileio_close_fxv, &
            fileio_open_phi, fileio_close_phi, &
            fileio_open_Al,  fileio_close_Al, &
            fileio_open_mom, fileio_close_mom, &
            fileio_open_trn, fileio_close_trn, &
            fileio_open_tri, fileio_close_tri, &
            
            fileio_read_cnt,  fileio_write_cnt, &
            fileio_write_fxv, fileio_write_phi, fileio_write_Al, &
            fileio_write_mom, fileio_write_trn, fileio_write_tri

CONTAINS

!--------------------------------------
  SUBROUTINE fileio_open_icnt ( path )
!--------------------------------------

    character(*), intent(in) :: path

    character(3) :: cold
    character(6) :: cranks, crankm, crankv, crankz, crankw
    character(1024) :: rank_dir_path

    write( cold, fmt="(i3.3)" ) inum-1
    write( cranks, fmt="(i0)" ) ranks
    write( crankm, fmt="(i0)" ) rankm
    write( crankv, fmt="(i0)" ) rankv
    write( crankz, fmt="(i0)" ) rankz
    write( crankw, fmt="(i0)" ) rankw

    path_icnt = path//'cnt.'//cold//'.zarr/'

    if ( rankg == 0 ) then
      open( icnt_t, file=trim(path_icnt)//'t/c/0',                   &
                    form='unformatted', status='old', action='read', &
                    access='stream', convert='LITTLE_ENDIAN' )
    end if

    rank_dir_path = trim(path_icnt)//'cnt/c/0/'//trim(cranks)//'/'// &
                    trim(crankm)//'/'//trim(crankv)//'/'//           &
                    trim(crankz)//'/'//trim(crankw)//'/'
    call system("mkdir -p "//trim(rank_dir_path))

    open( icnt, file=trim(rank_dir_path)//'0',                   &
                form='unformatted', status='old', action='read', &
                access='stream', convert='LITTLE_ENDIAN' )

  END SUBROUTINE fileio_open_icnt

!--------------------------------------
  SUBROUTINE fileio_close_icnt
!--------------------------------------

    if ( rankg == 0 ) then
      close( icnt_t )
    end if
    close( icnt )

  END SUBROUTINE fileio_close_icnt

 
!--------------------------------------
  SUBROUTINE fileio_open_cnt ( path )
!--------------------------------------

    character(*), intent(in) :: path

    character(3) :: cnew
    character(6) :: cranks, crankm, crankv, crankz, crankw
    character(1024) :: rank_dir_path

    write( cnew, fmt="(i3.3)" ) inum
    write( cranks, fmt="(i0)" ) ranks
    write( crankm, fmt="(i0)" ) rankm
    write( crankv, fmt="(i0)" ) rankv
    write( crankz, fmt="(i0)" ) rankz
    write( crankw, fmt="(i0)" ) rankw

    path_cnt = path//'cnt.'//cnew//'.zarr/'

    if ( rankg == 0 ) then
      call system("mkdir -p "//trim(path_cnt)//'t/c/')
      open( ocnt_t, file=trim(path_cnt)//'t/c/0',                         &
                    form='unformatted', status='replace', action='write', &
                    access='stream', convert='LITTLE_ENDIAN' )
    end if

    rank_dir_path = trim(path_cnt)//'cnt/c/0/'//trim(cranks)//'/'// &
                    trim(crankm)//'/'//trim(crankv)//'/'//          &
                    trim(crankz)//'/'//trim(crankw)//'/'
    call system("mkdir -p "//trim(rank_dir_path))
    open( ocnt, file=trim(rank_dir_path)//'0',                        &
                form='unformatted', status='replace', action='write', &
                access='stream', convert='LITTLE_ENDIAN' )

  END SUBROUTINE fileio_open_cnt

!--------------------------------------
  SUBROUTINE fileio_close_cnt
!--------------------------------------

    character(6) :: cranks, crankm, crankv, crankz, crankw

    write( cranks, fmt="(i0)" ) ranks
    write( crankm, fmt="(i0)" ) rankm
    write( crankv, fmt="(i0)" ) rankv
    write( crankz, fmt="(i0)" ) rankz
    write( crankw, fmt="(i0)" ) rankw

    if ( rankg == 0 ) then
      close( ocnt_t )
    end if
    close( ocnt )

    ! --- output coordinate binary files ---
    if ( rankm == 0 .and. rankv == 0 .and. rankz == 0 .and. rankw == 0 ) then
      call system("mkdir -p "//trim(path_cnt)//'is/c/')
      open( ocnt_t, file=trim(path_cnt)//'is/c/'//trim(cranks),           &
                    form='unformatted', status='replace', action='write', &
                    access='stream', convert='LITTLE_ENDIAN' )
      write( ocnt_t ) int(ranks, kind=4)
      close( ocnt_t )
    end if
    if ( ranks == 0 .and. rankv == 0 .and. rankz == 0 .and. rankw == 0 ) then
      if( vp_coord == 1 ) then
        call system("mkdir -p "//trim(path_cnt)//'vp/c/')
        open( ocnt_t, file=trim(path_cnt)//'vp/c/'//trim(crankm),           &
                      form='unformatted', status='replace', action='write', &
                      access='stream', convert='LITTLE_ENDIAN' )
        write( ocnt_t ) vp(0,:)
        close( ocnt_t )
      else
        call system("mkdir -p "//trim(path_cnt)//'mu/c/')
        open( ocnt_t, file=trim(path_cnt)//'mu/c/'//trim(crankm),           &
                      form='unformatted', status='replace', action='write', &
                      access='stream', convert='LITTLE_ENDIAN' )
        write( ocnt_t ) mu(0,:)
        close( ocnt_t )
      end if
    end if
    if ( ranks == 0 .and. rankm == 0 .and. rankz == 0 .and. rankw == 0 ) then
      call system("mkdir -p "//trim(path_cnt)//'vl/c/')
      open( ocnt_t, file=trim(path_cnt)//'vl/c/'//trim(crankv),           &
                    form='unformatted', status='replace', action='write', &
                    access='stream', convert='LITTLE_ENDIAN' )
      write( ocnt_t ) vl
      close( ocnt_t )
    end if
    if ( ranks == 0 .and. rankm == 0 .and. rankv == 0 .and. rankw == 0 ) then
      call system("mkdir -p "//trim(path_cnt)//'zz/c/')
      open( ocnt_t, file=trim(path_cnt)//'zz/c/'//trim(crankz),           &
                    form='unformatted', status='replace', action='write', &
                    access='stream', convert='LITTLE_ENDIAN' )
      write( ocnt_t ) zz
      close( ocnt_t )
    end if
    if ( ranks == 0 .and. rankm == 0 .and. rankv == 0 .and. rankz == 0 ) then
      call system("mkdir -p "//trim(path_cnt)//'ky/c/')
      open( ocnt_t, file=trim(path_cnt)//'ky/c/'//trim(crankw),           &
                    form='unformatted', status='replace', action='write', &
                    access='stream', convert='LITTLE_ENDIAN' )
      write( ocnt_t ) ky
      close( ocnt_t )
    end if
    if ( rankg == 0 ) then
      call system("mkdir -p "//trim(path_cnt)//'kx/c/')
      open( ocnt_t, file=trim(path_cnt)//'kx/c/0',                        &
                    form='unformatted', status='replace', action='write', &
                    access='stream', convert='LITTLE_ENDIAN' )
      write( ocnt_t ) kx
      close( ocnt_t )
    end if

    ! --- output json files ---
    if ( rankg == 0 ) then
      call write_zarr_json_head(trim(path_cnt)//'zarr.json')
      if( vp_coord == 1 ) then
        call write_zarr_json(trim(path_cnt)//'cnt/zarr.json',          & ! path
                             7,                                        & ! ndims
                             (/1, nprocs, nprocm*(nm+1), nprocv*(2*nv), nprocz*(2*nz), nprocw*(ny+1), 2*nx+1/), & ! shape
                             (/1, 1, nm+1, 2*nv, 2*nz, ny+1, 2*nx+1/), & ! chunk_shape
                             'complex128',                             & ! dtype
                             (/'t ','is','vp','vl','zz','ky','kx'/))     ! dimnames
      else
        call write_zarr_json(trim(path_cnt)//'cnt/zarr.json',          & ! path
                             7,                                        & ! ndims
                             (/1, nprocs, nprocm*(nm+1), nprocv*(2*nv), nprocz*(2*nz), nprocw*(ny+1), 2*nx+1/), & ! shape
                             (/1, 1, nm+1, 2*nv, 2*nz, ny+1, 2*nx+1/), & ! chunk_shape
                             'complex128',                             & ! dtype
                             (/'t ','is','mu','vl','zz','ky','kx'/))     ! dimnames
      end if
      call write_zarr_json(trim(path_cnt)//'t/zarr.json', & ! path
                           1,                             & ! ndims
                           (/1/),                         & ! shape
                           (/1/),                         & ! chunk_shape
                           'float64',                     & ! dtype
                           (/'t'/))                         ! dimnames
      call write_zarr_json(trim(path_cnt)//'is/zarr.json', & ! path
                           1,                              & ! ndims
                           (/nprocs/),                     & ! shape
                           (/1/),                          & ! chunk_shape
                           'int32',                        & ! dtype
                           (/'is'/))                         ! dimnames
      if ( vp_coord == 1 ) then
        call write_zarr_json(trim(path_cnt)//'vp/zarr.json', & ! path
                             1,                              & ! ndims
                             (/nprocm*(nm+1)/),              & ! shape
                             (/nm+1/),                       & ! chunk_shape
                             'float64',                      & ! dtype
                             (/'vp'/))                         ! dimnames
      else
        call write_zarr_json(trim(path_cnt)//'mu/zarr.json', & ! path
                             1,                              & ! ndims
                             (/nprocm*(nm+1)/),              & ! shape
                             (/nm+1/),                       & ! chunk_shape
                             'float64',                      & ! dtype
                             (/'mu'/))                         ! dimnames
      end if
      call write_zarr_json(trim(path_cnt)//'vl/zarr.json', & ! path
                           1,                              & ! ndims
                           (/nprocv*(2*nv)/),              & ! shape
                           (/2*nv/),                       & ! chunk_shape
                           'float64',                      & ! dtype
                           (/'vl'/))                         ! dimnames
      call write_zarr_json(trim(path_cnt)//'zz/zarr.json', & ! path
                           1,                              & ! ndims
                           (/nprocz*(2*nz)/),              & ! shape
                           (/2*nz/),                       & ! chunk_shape
                           'float64',                      & ! dtype
                           (/'zz'/))                         ! dimnames
      call write_zarr_json(trim(path_cnt)//'ky/zarr.json', & ! path
                           1,                              & ! ndims
                           (/nprocw*(ny+1)/),              & ! shape
                           (/ny+1/),                       & ! chunk_shape
                           'float64',                      & ! dtype
                           (/'ky'/))                         ! dimnames
      call write_zarr_json(trim(path_cnt)//'kx/zarr.json', & ! path
                           1,                              & ! ndims
                           (/2*nx+1/),                     & ! shape
                           (/2*nx+1/),                     & ! chunk_shape
                           'float64',                      & ! dtype
                           (/'kx'/))                         ! dimnames
    end if

  END SUBROUTINE fileio_close_cnt


!--------------------------------------
  SUBROUTINE fileio_open_fxv ( path )
!--------------------------------------

    character(*), intent(in) :: path

    character(3) :: cnew
    character(6) :: cranks, crankm, crankv, crankz, crankw
    character(1024) :: rank_dir_path

    write( cnew, fmt="(i3.3)" ) inum
    write( cranks, fmt="(i0)" ) ranks
    write( crankm, fmt="(i0)" ) rankm
    write( crankv, fmt="(i0)" ) rankv
    write( crankz, fmt="(i0)" ) rankz
    write( crankw, fmt="(i0)" ) rankw

    path_fxv = path//'fxv.'//cnew//'.zarr/'

    if ( rankg == 0 ) then
      call system("mkdir -p "//trim(path_fxv)//'t/c/')
      open( ofxv_t, file=trim(path_fxv)//'t/c/0',                         &
                    form='unformatted', status='replace', action='write', &
                    access='stream', convert='LITTLE_ENDIAN' )
    end if
    nt_out_fxv = 0

    rank_dir_path = trim(path_fxv)//'fxv/c/0/'//trim(cranks)//'/'// &
                    trim(crankm)//'/'//trim(crankv)//'/'//          &
                    trim(crankz)//'/'//trim(crankw)//'/'
    call system("mkdir -p "//trim(rank_dir_path))
    open( ofxv, file=trim(rank_dir_path)//'0',                        &
                form='unformatted', status='replace', action='write', &
                access='stream', convert='LITTLE_ENDIAN' )

  END SUBROUTINE fileio_open_fxv

!--------------------------------------
  SUBROUTINE fileio_close_fxv
!--------------------------------------

    character(6) :: cranks, crankm, crankv, crankz, crankw

    write( cranks, fmt="(i0)" ) ranks
    write( crankm, fmt="(i0)" ) rankm
    write( crankv, fmt="(i0)" ) rankv
    write( crankz, fmt="(i0)" ) rankz
    write( crankw, fmt="(i0)" ) rankw

    if ( rankg == 0 ) then
      close( ofxv_t )
    end if
    close( ofxv )

    ! --- output coordinate binary files ---
    if ( rankm == 0 .and. rankv == 0 .and. rankz == 0 .and. rankw == 0 ) then
      call system("mkdir -p "//trim(path_fxv)//'is/c/')
      open( ofxv_t, file=trim(path_fxv)//'is/c/'//trim(cranks),           &
                    form='unformatted', status='replace', action='write', &
                    access='stream', convert='LITTLE_ENDIAN' )
      write( ofxv_t ) int(ranks, kind=4)
      close( ofxv_t )
    end if
    if ( ranks == 0 .and. rankv == 0 .and. rankz == 0 .and. rankw == 0 ) then
      if( vp_coord == 1 ) then
        call system("mkdir -p "//trim(path_fxv)//'vp/c/')
        open( ofxv_t, file=trim(path_fxv)//'vp/c/'//trim(crankm),           &
                      form='unformatted', status='replace', action='write', &
                      access='stream', convert='LITTLE_ENDIAN' )
        write( ofxv_t ) vp(0,:)
        close( ofxv_t )
      else
        call system("mkdir -p "//trim(path_fxv)//'mu/c/')
        open( ofxv_t, file=trim(path_fxv)//'mu/c/'//trim(crankm),           &
                      form='unformatted', status='replace', action='write', &
                      access='stream', convert='LITTLE_ENDIAN' )
        write( ofxv_t ) mu(0,:)
        close( ofxv_t )
      end if
    end if
    if ( ranks == 0 .and. rankm == 0 .and. rankz == 0 .and. rankw == 0 ) then
      call system("mkdir -p "//trim(path_fxv)//'vl/c/')
      open( ofxv_t, file=trim(path_fxv)//'vl/c/'//trim(crankv),           &
                    form='unformatted', status='replace', action='write', &
                    access='stream', convert='LITTLE_ENDIAN' )
      write( ofxv_t ) vl
      close( ofxv_t )
    end if
    if ( ranks == 0 .and. rankm == 0 .and. rankv == 0 .and. rankw == 0 ) then
      call system("mkdir -p "//trim(path_fxv)//'zz/c/')
      open( ofxv_t, file=trim(path_fxv)//'zz/c/'//trim(crankz),           &
                    form='unformatted', status='replace', action='write', &
                    access='stream', convert='LITTLE_ENDIAN' )
      write( ofxv_t ) zz(-nz)
      close( ofxv_t )
    end if
    if ( ranks == 0 .and. rankm == 0 .and. rankv == 0 .and. rankz == 0 ) then
      call system("mkdir -p "//trim(path_fxv)//'ky/c/')
      open( ofxv_t, file=trim(path_fxv)//'ky/c/'//trim(crankw),           &
                    form='unformatted', status='replace', action='write', &
                    access='stream', convert='LITTLE_ENDIAN' )
      write( ofxv_t ) ky
      close( ofxv_t )
    end if
    if ( rankg == 0 ) then
      call system("mkdir -p "//trim(path_fxv)//'kx/c/')
      open( ofxv_t, file=trim(path_fxv)//'kx/c/0',                        &
                    form='unformatted', status='replace', action='write', &
                    access='stream', convert='LITTLE_ENDIAN' )
      write( ofxv_t ) kx
      close( ofxv_t )
    end if

    ! --- output json files ---
    if ( rankg == 0 ) then
      call write_zarr_json_head(trim(path_fxv)//'zarr.json')
      if( vp_coord == 1 ) then
        call write_zarr_json(trim(path_fxv)//'fxv/zarr.json',          & ! path
                             7,                                        & ! ndims
                             (/nt_out_fxv, nprocs, nprocm*(nm+1), nprocv*(2*nv), nprocz, nprocw*(ny+1), 2*nx+1/), & ! shape
                             (/nt_out_fxv, 1, nm+1, 2*nv, 1, ny+1, 2*nx+1/), & ! chunk_shape
                             'complex128',                             & ! dtype
                             (/'t ','is','vp','vl','zz','ky','kx'/))     ! dimnames
      else
        call write_zarr_json(trim(path_fxv)//'fxv/zarr.json',          & ! path
                             7,                                        & ! ndims
                             (/nt_out_fxv, nprocs, nprocm*(nm+1), nprocv*(2*nv), nprocz, nprocw*(ny+1), 2*nx+1/), & ! shape
                             (/nt_out_fxv, 1, nm+1, 2*nv, 1, ny+1, 2*nx+1/), & ! chunk_shape
                             'complex128',                             & ! dtype
                             (/'t ','is','mu','vl','zz','ky','kx'/))     ! dimnames
      end if
      call write_zarr_json(trim(path_fxv)//'t/zarr.json', & ! path
                           1,                             & ! ndims
                           (/nt_out_fxv/),                & ! shape
                           (/nt_out_fxv/),                & ! chunk_shape
                           'float64',                     & ! dtype
                           (/'t'/))                         ! dimnames
      call write_zarr_json(trim(path_fxv)//'is/zarr.json', & ! path
                           1,                              & ! ndims
                           (/nprocs/),                     & ! shape
                           (/1/),                          & ! chunk_shape
                           'int32',                        & ! dtype
                           (/'is'/))                         ! dimnames
      if ( vp_coord == 1 ) then
        call write_zarr_json(trim(path_fxv)//'vp/zarr.json', & ! path
                             1,                              & ! ndims
                             (/nprocm*(nm+1)/),              & ! shape
                             (/nm+1/),                       & ! chunk_shape
                             'float64',                      & ! dtype
                             (/'vp'/))                         ! dimnames
      else
        call write_zarr_json(trim(path_fxv)//'mu/zarr.json', & ! path
                             1,                              & ! ndims
                             (/nprocm*(nm+1)/),              & ! shape
                             (/nm+1/),                       & ! chunk_shape
                             'float64',                      & ! dtype
                             (/'mu'/))                         ! dimnames
      end if
      call write_zarr_json(trim(path_fxv)//'vl/zarr.json', & ! path
                           1,                              & ! ndims
                           (/nprocv*(2*nv)/),              & ! shape
                           (/2*nv/),                       & ! chunk_shape
                           'float64',                      & ! dtype
                           (/'vl'/))                         ! dimnames
      call write_zarr_json(trim(path_fxv)//'zz/zarr.json', & ! path
                           1,                              & ! ndims
                           (/nprocz/),                     & ! shape
                           (/1/),                          & ! chunk_shape
                           'float64',                      & ! dtype
                           (/'zz'/))                         ! dimnames
      call write_zarr_json(trim(path_fxv)//'ky/zarr.json', & ! path
                           1,                              & ! ndims
                           (/nprocw*(ny+1)/),              & ! shape
                           (/ny+1/),                       & ! chunk_shape
                           'float64',                      & ! dtype
                           (/'ky'/))                         ! dimnames
      call write_zarr_json(trim(path_fxv)//'kx/zarr.json', & ! path
                           1,                              & ! ndims
                           (/2*nx+1/),                     & ! shape
                           (/2*nx+1/),                     & ! chunk_shape
                           'float64',                      & ! dtype
                           (/'kx'/))                         ! dimnames
    end if

  END SUBROUTINE fileio_close_fxv


!--------------------------------------
  SUBROUTINE fileio_open_phi ( path )
!--------------------------------------

    character(*), intent(in) :: path

    character(3) :: cnew
    character(6) :: crankz, crankw
    character(1024) :: rank_dir_path

    if ( (ranks /= 0) .OR. (vel_rank /= 0) ) return

    write( cnew, fmt="(i3.3)" ) inum
    write( crankz, fmt="(i0)" ) rankz
    write( crankw, fmt="(i0)" ) rankw

    path_phi = path//'phi.'//cnew//'.zarr/'

    if ( rankg == 0 ) then
      call system("mkdir -p "//trim(path_phi)//'t/c/')
      open( ophi_t, file=trim(path_phi)//'t/c/0',                         &
                    form='unformatted', status='replace', action='write', &
                    access='stream', convert='LITTLE_ENDIAN' )
    end if
    nt_out_phi = 0

    rank_dir_path = trim(path_phi)//'phi/c/0/'//trim(crankz)//'/'// &
                    trim(crankw)//'/'
    call system("mkdir -p "//trim(rank_dir_path))
    open( ophi, file=trim(rank_dir_path)//'0',                        &
                form='unformatted', status='replace', action='write', &
                access='stream', convert='LITTLE_ENDIAN' )

  END SUBROUTINE fileio_open_phi

!--------------------------------------
  SUBROUTINE fileio_close_phi
!--------------------------------------

    character(6) :: crankz, crankw

    if ( (ranks /= 0) .OR. (vel_rank /= 0) ) return

    write( crankz, fmt="(i0)" ) rankz
    write( crankw, fmt="(i0)" ) rankw

    if ( rankg == 0 ) then
      close( ophi_t )
    end if
    close( ophi )

    ! --- output coordinate binary files ---
    if ( ranks == 0 .and. rankm == 0 .and. rankv == 0 .and. rankw == 0 ) then
      call system("mkdir -p "//trim(path_phi)//'zz/c/')
      open( ophi_t, file=trim(path_phi)//'zz/c/'//trim(crankz),           &
                    form='unformatted', status='replace', action='write', &
                    access='stream', convert='LITTLE_ENDIAN' )
      write( ophi_t ) zz
      close( ophi_t )
    end if
    if ( ranks == 0 .and. rankm == 0 .and. rankv == 0 .and. rankz == 0 ) then
      call system("mkdir -p "//trim(path_phi)//'ky/c/')
      open( ophi_t, file=trim(path_phi)//'ky/c/'//trim(crankw),           &
                    form='unformatted', status='replace', action='write', &
                    access='stream', convert='LITTLE_ENDIAN' )
      write( ophi_t ) ky
      close( ophi_t )
    end if
    if ( rankg == 0 ) then
      call system("mkdir -p "//trim(path_phi)//'kx/c/')
      open( ophi_t, file=trim(path_phi)//'kx/c/0',                        &
                    form='unformatted', status='replace', action='write', &
                    access='stream', convert='LITTLE_ENDIAN' )
      write( ophi_t ) kx
      close( ophi_t )
    end if

    ! --- output json files ---
    if ( rankg == 0 ) then
      call write_zarr_json_head(trim(path_phi)//'zarr.json')
      call write_zarr_json(trim(path_phi)//'phi/zarr.json',    & ! path
                           4,                                  & ! ndims
                           (/nt_out_phi, nprocz*(2*nz), nprocw*(ny+1), 2*nx+1/), & ! shape
                           (/nt_out_phi, 2*nz, ny+1, 2*nx+1/), & ! chunk_shape
                           'complex128',                       & ! dtype
                           (/'t ','zz','ky','kx'/))              ! dimnames
      call write_zarr_json(trim(path_phi)//'t/zarr.json', & ! path
                           1,                             & ! ndims
                           (/nt_out_phi/),                & ! shape
                           (/nt_out_phi/),                & ! chunk_shape
                           'float64',                     & ! dtype
                           (/'t'/))                         ! dimnames
      call write_zarr_json(trim(path_phi)//'zz/zarr.json', & ! path
                           1,                              & ! ndims
                           (/nprocz*(2*nz)/),              & ! shape
                           (/2*nz/),                       & ! chunk_shape
                           'float64',                      & ! dtype
                           (/'zz'/))                         ! dimnames
      call write_zarr_json(trim(path_phi)//'ky/zarr.json', & ! path
                           1,                              & ! ndims
                           (/nprocw*(ny+1)/),              & ! shape
                           (/ny+1/),                       & ! chunk_shape
                           'float64',                      & ! dtype
                           (/'ky'/))                         ! dimnames
      call write_zarr_json(trim(path_phi)//'kx/zarr.json', & ! path
                           1,                              & ! ndims
                           (/2*nx+1/),                     & ! shape
                           (/2*nx+1/),                     & ! chunk_shape
                           'float64',                      & ! dtype
                           (/'kx'/))                         ! dimnames
    end if

  END SUBROUTINE fileio_close_phi


!--------------------------------------
  SUBROUTINE fileio_open_Al ( path )
!--------------------------------------

    character(*), intent(in) :: path

    character(3) :: cnew
    character(6) :: crankz, crankw
    character(1024) :: rank_dir_path

    if ( (ranks /= 0) .OR. (vel_rank /= 0) ) return

    write( cnew, fmt="(i3.3)" ) inum
    write( crankz, fmt="(i0)" ) rankz
    write( crankw, fmt="(i0)" ) rankw

    path_Al = path//'Al.'//cnew//'.zarr/'

    if ( rankg == 0 ) then
      call system("mkdir -p "//trim(path_Al)//'t/c/')
      open( oAl_t, file=trim(path_Al)//'t/c/0',                           &
                    form='unformatted', status='replace', action='write', &
                    access='stream', convert='LITTLE_ENDIAN' )
    end if
    nt_out_Al = 0

    rank_dir_path = trim(path_Al)//'Al/c/0/'//trim(crankz)//'/'// &
                    trim(crankw)//'/'
    call system("mkdir -p "//trim(rank_dir_path))
    open( oAl, file=trim(rank_dir_path)//'0',                         &
                form='unformatted', status='replace', action='write', &
                access='stream', convert='LITTLE_ENDIAN' )

  END SUBROUTINE fileio_open_Al

!--------------------------------------
  SUBROUTINE fileio_close_Al
!--------------------------------------

    character(6) :: crankz, crankw

    if ( (ranks /= 0) .OR. (vel_rank /= 0) ) return

    write( crankz, fmt="(i0)" ) rankz
    write( crankw, fmt="(i0)" ) rankw

    if ( rankg == 0 ) then
      close( oAl_t )
    end if
    close( oAl )

    ! --- output coordinate binary files ---
    if ( ranks == 0 .and. rankm == 0 .and. rankv == 0 .and. rankw == 0 ) then
      call system("mkdir -p "//trim(path_Al)//'zz/c/')
      open( oAl_t, file=trim(path_Al)//'zz/c/'//trim(crankz),             &
                    form='unformatted', status='replace', action='write', &
                    access='stream', convert='LITTLE_ENDIAN' )
      write( oAl_t ) zz
      close( oAl_t )
    end if
    if ( ranks == 0 .and. rankm == 0 .and. rankv == 0 .and. rankz == 0 ) then
      call system("mkdir -p "//trim(path_Al)//'ky/c/')
      open( oAl_t, file=trim(path_Al)//'ky/c/'//trim(crankw),             &
                    form='unformatted', status='replace', action='write', &
                    access='stream', convert='LITTLE_ENDIAN' )
      write( oAl_t ) ky
      close( oAl_t )
    end if
    if ( rankg == 0 ) then
      call system("mkdir -p "//trim(path_Al)//'kx/c/')
      open( oAl_t, file=trim(path_Al)//'kx/c/0',                          &
                    form='unformatted', status='replace', action='write', &
                    access='stream', convert='LITTLE_ENDIAN' )
      write( oAl_t ) kx
      close( oAl_t )
    end if

    ! --- output json files ---
    if ( rankg == 0 ) then
      call write_zarr_json_head(trim(path_Al)//'zarr.json')
      call write_zarr_json(trim(path_Al)//'Al/zarr.json',     & ! path
                           4,                                 & ! ndims
                           (/nt_out_Al, nprocz*(2*nz), nprocw*(ny+1), 2*nx+1/), & ! shape
                           (/nt_out_Al, 2*nz, ny+1, 2*nx+1/), & ! chunk_shape
                           'complex128',                      & ! dtype
                           (/'t ','zz','ky','kx'/))             ! dimnames
      call write_zarr_json(trim(path_Al)//'t/zarr.json', & ! path
                           1,                            & ! ndims
                           (/nt_out_Al/),                & ! shape
                           (/nt_out_Al/),                & ! chunk_shape
                           'float64',                    & ! dtype
                           (/'t'/))                        ! dimnames
      call write_zarr_json(trim(path_Al)//'zz/zarr.json', & ! path
                           1,                             & ! ndims
                           (/nprocz*(2*nz)/),             & ! shape
                           (/2*nz/),                      & ! chunk_shape
                           'float64',                     & ! dtype
                           (/'zz'/))                        ! dimnames
      call write_zarr_json(trim(path_Al)//'ky/zarr.json', & ! path
                           1,                             & ! ndims
                           (/nprocw*(ny+1)/),             & ! shape
                           (/ny+1/),                      & ! chunk_shape
                           'float64',                     & ! dtype
                           (/'ky'/))                        ! dimnames
      call write_zarr_json(trim(path_Al)//'kx/zarr.json', & ! path
                           1,                             & ! ndims
                           (/2*nx+1/),                    & ! shape
                           (/2*nx+1/),                    & ! chunk_shape
                           'float64',                     & ! dtype
                           (/'kx'/))                        ! dimnames
    end if

  END SUBROUTINE fileio_close_Al


!--------------------------------------
  SUBROUTINE fileio_open_mom ( path )
!--------------------------------------

    character(*), intent(in) :: path

    character(3)   :: cnew
    character(6)   :: cranks, crankz, crankw
    character(1024) :: rank_dir_path

    if ( vel_rank /= 0 ) return

    write( cnew, fmt="(i3.3)" ) inum
    write( cranks, fmt="(i0)" ) ranks
    write( crankz, fmt="(i0)" ) rankz
    write( crankw, fmt="(i0)" ) rankw

    path_mom = path//'mom.'//cnew//'.zarr/'

    if ( rankg == 0 ) then
      call system("mkdir -p "//trim(path_mom)//'t/c/')
      open( omom_t, file=trim(path_mom)//'t/c/0',                         &
                    form='unformatted', status='replace', action='write', &
                    access='stream', convert='LITTLE_ENDIAN' )
    end if
    nt_out_mom = 0

    rank_dir_path = trim(path_mom)//'mom/c/0/'//trim(cranks)//'/0/'// &
                    trim(crankz)//'/'//trim(crankw)//'/'
    call system("mkdir -p "//trim(rank_dir_path))
    open( omom, file=trim(rank_dir_path)//'0',                        &
                form='unformatted', status='replace', action='write', &
                access='stream', convert='LITTLE_ENDIAN' )

  END SUBROUTINE fileio_open_mom

!--------------------------------------
  SUBROUTINE fileio_close_mom
!--------------------------------------

    character(6) :: cranks, crankz, crankw
    integer :: imom

    if ( vel_rank /= 0 ) return

    write( cranks, fmt="(i0)" ) ranks
    write( crankz, fmt="(i0)" ) rankz
    write( crankw, fmt="(i0)" ) rankw

    if ( rankg == 0 ) then
      close( omom_t )
    end if
    close( omom )

    ! --- output coordinate binary files ---
    if ( rankm == 0 .and. rankv == 0 .and. rankz == 0 .and. rankw == 0 ) then
      call system("mkdir -p "//trim(path_mom)//'is/c/')
      open( omom_t, file=trim(path_mom)//'is/c/'//trim(cranks),           &
                    form='unformatted', status='replace', action='write', &
                    access='stream', convert='LITTLE_ENDIAN' )
      write( omom_t ) int(ranks, kind=4)
      close( omom_t )
    end if
    if ( rankg == 0 ) then
      call system("mkdir -p "//trim(path_mom)//'imom/c/')
      open( omom_t, file=trim(path_mom)//'imom/c/0',                        &
                    form='unformatted', status='replace', action='write', &
                    access='stream', convert='LITTLE_ENDIAN' )
      write( omom_t ) (/ (imom, imom=0,nmom-1) /)
      close( omom_t )
    end if
    if ( ranks == 0 .and. rankm == 0 .and. rankv == 0 .and. rankw == 0 ) then
      call system("mkdir -p "//trim(path_mom)//'zz/c/')
      open( omom_t, file=trim(path_mom)//'zz/c/'//trim(crankz),           &
                    form='unformatted', status='replace', action='write', &
                    access='stream', convert='LITTLE_ENDIAN' )
      write( omom_t ) zz
      close( omom_t )
    end if
    if ( ranks == 0 .and. rankm == 0 .and. rankv == 0 .and. rankz == 0 ) then
      call system("mkdir -p "//trim(path_mom)//'ky/c/')
      open( omom_t, file=trim(path_mom)//'ky/c/'//trim(crankw),           &
                    form='unformatted', status='replace', action='write', &
                    access='stream', convert='LITTLE_ENDIAN' )
      write( omom_t ) ky
      close( omom_t )
    end if
    if ( rankg == 0 ) then
      call system("mkdir -p "//trim(path_mom)//'kx/c/')
      open( omom_t, file=trim(path_mom)//'kx/c/0',                        &
                    form='unformatted', status='replace', action='write', &
                    access='stream', convert='LITTLE_ENDIAN' )
      write( omom_t ) kx
      close( omom_t )
    end if

    ! --- output json files ---
    if ( rankg == 0 ) then
      call write_zarr_json_head(trim(path_mom)//'zarr.json')
      call write_zarr_json(trim(path_mom)//'mom/zarr.json',             & ! path
                           6,                                           & ! ndims
                           (/nt_out_mom, nprocs, nmom, nprocz*(2*nz), nprocw*(ny+1), 2*nx+1/), & ! shape
                           (/nt_out_mom, 1, nmom, 2*nz, ny+1, 2*nx+1/), & ! chunk_shape
                           'complex128',                                & ! dtype
                           (/'t   ','is  ','imom','zz  ','ky  ','kx  '/)) ! dimnames
      call write_zarr_json(trim(path_mom)//'t/zarr.json', & ! path
                           1,                             & ! ndims
                           (/nt_out_mom/),                & ! shape
                           (/nt_out_mom/),                & ! chunk_shape
                           'float64',                     & ! dtype
                           (/'t'/))                         ! dimnames
      call write_zarr_json(trim(path_mom)//'is/zarr.json', & ! path
                           1,                              & ! ndims
                           (/nprocs/),                     & ! shape
                           (/1/),                          & ! chunk_shape
                           'int32',                        & ! dtype
                           (/'is'/))                         ! dimnames
      call write_zarr_json(trim(path_mom)//'imom/zarr.json', & ! path
                           1,                                & ! ndims
                           (/nmom/),                         & ! shape
                           (/nmom/),                         & ! chunk_shape
                           'int32',                          & ! dtype
                           (/'imom'/))                         ! dimnames
      call write_zarr_json(trim(path_mom)//'zz/zarr.json', & ! path
                           1,                              & ! ndims
                           (/nprocz*(2*nz)/),              & ! shape
                           (/2*nz/),                       & ! chunk_shape
                           'float64',                      & ! dtype
                           (/'zz'/))                         ! dimnames
      call write_zarr_json(trim(path_mom)//'ky/zarr.json', & ! path
                           1,                              & ! ndims
                           (/nprocw*(ny+1)/),              & ! shape
                           (/ny+1/),                       & ! chunk_shape
                           'float64',                      & ! dtype
                           (/'ky'/))                         ! dimnames
      call write_zarr_json(trim(path_mom)//'kx/zarr.json', & ! path
                           1,                              & ! ndims
                           (/2*nx+1/),                     & ! shape
                           (/2*nx+1/),                     & ! chunk_shape
                           'float64',                      & ! dtype
                           (/'kx'/))                         ! dimnames
    end if

  END SUBROUTINE fileio_close_mom


!--------------------------------------
  SUBROUTINE fileio_open_trn ( path )
!--------------------------------------

    character(*), intent(in) :: path

    character(3)   :: cnew
    character(6)   :: cranks, crankw
    character(1024) :: rank_dir_path

    if ( (zsp_rank /= 0) .OR. (vel_rank /= 0) ) return

    write( cnew, fmt="(i3.3)" ) inum
    write( cranks, fmt="(i0)" ) ranks
    write( crankw, fmt="(i0)" ) rankw

    path_trn = path//'trn.'//cnew//'.zarr/'

    if ( rankg == 0 ) then
      call system("mkdir -p "//trim(path_trn)//'t/c/')
      open( otrn_t, file=trim(path_trn)//'t/c/0',                         &
                    form='unformatted', status='replace', action='write', &
                    access='stream', convert='LITTLE_ENDIAN' )
    end if
    nt_out_trn = 0

    rank_dir_path = trim(path_trn)//'trn/c/0/'//trim(cranks)//'/0/'// &
                    trim(crankw)//'/'
    call system("mkdir -p "//trim(rank_dir_path))
    open( otrn, file=trim(rank_dir_path)//'0',                        &
                form='unformatted', status='replace', action='write', &
                access='stream', convert='LITTLE_ENDIAN' )

  END SUBROUTINE fileio_open_trn

!--------------------------------------
  SUBROUTINE fileio_close_trn
!--------------------------------------

    character(6) :: cranks, crankw
    integer :: itrn

    if ( (zsp_rank /= 0) .OR. (vel_rank /= 0) ) return

    write( cranks, fmt="(i0)" ) ranks
    write( crankw, fmt="(i0)" ) rankw

    if ( rankg == 0 ) then
      close( otrn_t )
    end if
    close( otrn )

    ! --- output coordinate binary files ---
    if ( rankm == 0 .and. rankv == 0 .and. rankz == 0 .and. rankw == 0 ) then
      call system("mkdir -p "//trim(path_trn)//'is/c/')
      open( otrn_t, file=trim(path_trn)//'is/c/'//trim(cranks),           &
                    form='unformatted', status='replace', action='write', &
                    access='stream', convert='LITTLE_ENDIAN' )
      write( otrn_t ) int(ranks, kind=4)
      close( otrn_t )
    end if
    if ( rankg == 0 ) then
      call system("mkdir -p "//trim(path_trn)//'itrn/c/')
      open( otrn_t, file=trim(path_trn)//'itrn/c/0',                        &
                    form='unformatted', status='replace', action='write', &
                    access='stream', convert='LITTLE_ENDIAN' )
      write( otrn_t ) (/ (itrn, itrn=0,ntrn-1) /)
      close( otrn_t )
    end if
    if ( ranks == 0 .and. rankm == 0 .and. rankv == 0 .and. rankz == 0 ) then
      call system("mkdir -p "//trim(path_trn)//'ky/c/')
      open( otrn_t, file=trim(path_trn)//'ky/c/'//trim(crankw),           &
                    form='unformatted', status='replace', action='write', &
                    access='stream', convert='LITTLE_ENDIAN' )
      write( otrn_t ) ky
      close( otrn_t )
    end if
    if ( rankg == 0 ) then
      call system("mkdir -p "//trim(path_trn)//'kx/c/')
      open( otrn_t, file=trim(path_trn)//'kx/c/0',                        &
                    form='unformatted', status='replace', action='write', &
                    access='stream', convert='LITTLE_ENDIAN' )
      write( otrn_t ) kx
      close( otrn_t )
    end if

    ! --- output json files ---
    if ( rankg == 0 ) then
      call write_zarr_json_head(trim(path_trn)//'zarr.json')
      call write_zarr_json(trim(path_trn)//'trn/zarr.json',       & ! path
                           5,                                     & ! ndims
                           (/nt_out_trn, nprocs, ntrn, nprocw*(ny+1), 2*nx+1/), & ! shape
                           (/nt_out_trn, 1, ntrn, ny+1, 2*nx+1/), & ! chunk_shape
                           'float64',                             & ! dtype
                           (/'t   ','is  ','itrn','ky  ','kx  '/))  ! dimnames
      call write_zarr_json(trim(path_trn)//'t/zarr.json', & ! path
                           1,                             & ! ndims
                           (/nt_out_trn/),                & ! shape
                           (/nt_out_trn/),                & ! chunk_shape
                           'float64',                     & ! dtype
                           (/'t'/))                         ! dimnames
      call write_zarr_json(trim(path_trn)//'is/zarr.json', & ! path
                           1,                              & ! ndims
                           (/nprocs/),                     & ! shape
                           (/1/),                          & ! chunk_shape
                           'int32',                        & ! dtype
                           (/'is'/))                         ! dimnames
      call write_zarr_json(trim(path_trn)//'itrn/zarr.json', & ! path
                           1,                                & ! ndims
                           (/ntrn/),                         & ! shape
                           (/ntrn/),                         & ! chunk_shape
                           'int32',                          & ! dtype
                           (/'itrn'/))                         ! dimnames
      call write_zarr_json(trim(path_trn)//'ky/zarr.json', & ! path
                           1,                              & ! ndims
                           (/nprocw*(ny+1)/),              & ! shape
                           (/ny+1/),                       & ! chunk_shape
                           'float64',                      & ! dtype
                           (/'ky'/))                         ! dimnames
      call write_zarr_json(trim(path_trn)//'kx/zarr.json', & ! path
                           1,                              & ! ndims
                           (/2*nx+1/),                     & ! shape
                           (/2*nx+1/),                     & ! chunk_shape
                           'float64',                      & ! dtype
                           (/'kx'/))                         ! dimnames
    end if

  END SUBROUTINE fileio_close_trn


!--------------------------------------
  SUBROUTINE fileio_open_tri ( path, cmx, cmy, replace )
!--------------------------------------

    character(*), intent(in) :: path
    character(*), intent(in) :: cmx, cmy
    logical, intent(in) :: replace

    character(3)   :: cnew
    character(6)   :: cranks
    character(1024) :: rank_dir_path
    real(kind=DP), dimension(-global_ny:global_ny) :: gky
    integer :: my, itri

    if ( rank /= 0 ) return

    write( cnew, fmt="(i3.3)" ) inum
    write( cranks, fmt="(i0)" ) ranks

    path_tri = path//'tri.mx'//cmx//'my'//cmy//'.'//cnew//'.zarr/'

    if ( replace ) then ! Open a new zarr store

      ! --- output coordinate binary files ---
      if ( rankg == 0 ) then
        call system("mkdir -p "//trim(path_tri)//'t/c/')
        open( otri_t, file=trim(path_tri)//'t/c/0',                         &
                      form='unformatted', status='replace', action='write', &
                      access='stream', convert='LITTLE_ENDIAN' )
      end if
      if ( rankm == 0 .and. rankv == 0 .and. rankz == 0 .and. rankw == 0 ) then
        call system("mkdir -p "//trim(path_tri)//'is/c/')
        open( otri_t, file=trim(path_tri)//'is/c/'//trim(cranks),           &
                      form='unformatted', status='replace', action='write', &
                      access='stream', convert='LITTLE_ENDIAN' )
        write( otri_t ) int(ranks, kind=4)
        close( otri_t )
      end if
      if ( rankg == 0 ) then
        call system("mkdir -p "//trim(path_tri)//'itri/c/')
        open( otri_t, file=trim(path_tri)//'itri/c/0',                      &
                      form='unformatted', status='replace', action='write', &
                      access='stream', convert='LITTLE_ENDIAN' )
        write( otri_t ) (/ (itri, itri=0,ntri-1) /)
        close( otri_t )
      end if
      if ( rankg == 0 ) then
        do my = -global_ny, global_ny
          gky(my) = kymin_g * real( my, kind=DP )
        end do
        call system("mkdir -p "//trim(path_tri)//'ky/c/')
        open( otri_t, file=trim(path_tri)//'ky/c/0',                        &
                      form='unformatted', status='replace', action='write', &
                      access='stream', convert='LITTLE_ENDIAN' )
        write( otri_t ) gky
        close( otri_t )
      end if
      if ( rankg == 0 ) then
        call system("mkdir -p "//trim(path_tri)//'kx/c/')
        open( otri_t, file=trim(path_tri)//'kx/c/0',                        &
                      form='unformatted', status='replace', action='write', &
                      access='stream', convert='LITTLE_ENDIAN' )
        write( otri_t ) kx
        close( otri_t )
      end if
      nt_out_tri = nt_out_mom ! Assume the same call count for *mom* and *tri*.

      ! --- output json files ---
      if ( rankg == 0 ) then
        call write_zarr_json_head(trim(path_tri)//'zarr.json')
        call write_zarr_json(trim(path_tri)//'is/zarr.json', & ! path
                             1,                              & ! ndims
                             (/nprocs/),                     & ! shape
                             (/1/),                          & ! chunk_shape
                             'int32',                        & ! dtype
                             (/'is'/))                         ! dimnames
        call write_zarr_json(trim(path_tri)//'itri/zarr.json', & ! path
                             1,                                & ! ndims
                             (/ntri/),                         & ! shape
                             (/ntri/),                         & ! chunk_shape
                             'int32',                          & ! dtype
                             (/'itri'/))                         ! dimnames
        call write_zarr_json(trim(path_tri)//'ky/zarr.json', & ! path
                             1,                              & ! ndims
                             (/2*global_ny+1/),              & ! shape
                             (/2*global_ny+1/),              & ! chunk_shape
                             'float64',                      & ! dtype
                             (/'ky'/))                         ! dimnames
        call write_zarr_json(trim(path_tri)//'kx/zarr.json', & ! path
                             1,                              & ! ndims
                             (/2*nx+1/),                     & ! shape
                             (/2*nx+1/),                     & ! chunk_shape
                             'float64',                      & ! dtype
                             (/'kx'/))                         ! dimnames
      end if

      rank_dir_path = trim(path_tri)//'tri/c/0/'//trim(cranks)//'/0/0/'
      call system("mkdir -p "//trim(rank_dir_path))
      open( otri, file=trim(rank_dir_path)//'0',                        &
                  form='unformatted', status='replace', action='write', &
                  access='stream', convert='LITTLE_ENDIAN' )

    else ! append data to an existing zarr store

      if ( rankg == 0 ) then
        call system("mkdir -p "//trim(path_tri)//'t/c/')
        open( otri_t, file=trim(path_tri)//'t/c/0',                     &
                      form='unformatted', status='old', action='write', &
                      access='stream', convert='LITTLE_ENDIAN' )
      end if
      rank_dir_path = trim(path_tri)//'tri/c/0/'//trim(cranks)//'/0/0/'
      call system("mkdir -p "//trim(rank_dir_path))
      open( otri, file=trim(rank_dir_path)//'0',                    &
                  form='unformatted', status='old', action='write', &
                  access='stream', convert='LITTLE_ENDIAN' )

    end if

  END SUBROUTINE fileio_open_tri

!--------------------------------------
  SUBROUTINE fileio_close_tri
!--------------------------------------

    character(6) :: cranks, crankw

    if ( rank /= 0 ) return

    write( cranks, fmt="(i0)" ) ranks
    write( crankw, fmt="(i0)" ) rankw

    if ( rankg == 0 ) then
      close( otri_t )
    end if
    close( otri )

    nt_out_tri = nt_out_mom ! Assume the same call count for *mom* and *tri*.
    if ( rankg == 0 ) then
      call write_zarr_json(trim(path_tri)//'tri/zarr.json',       & ! path
                           5,                                     & ! ndims
                           (/nt_out_tri, nprocs, ntri, 2*global_ny+1, 2*nx+1/), & ! shape
                           (/nt_out_tri, 1, ntri, 2*global_ny+1, 2*nx+1/), & ! chunk_shape
                           'float64',                             & ! dtype
                           (/'t   ','is  ','itri','ky  ','kx  '/))  ! dimnames
      call write_zarr_json(trim(path_tri)//'t/zarr.json', & ! path
                           1,                             & ! ndims
                           (/nt_out_tri/),                & ! shape
                           (/nt_out_tri/),                & ! chunk_shape
                           'float64',                     & ! dtype
                           (/'t'/))                         ! dimnames
    end if

  END SUBROUTINE fileio_close_tri





!--------------------------------------
  SUBROUTINE fileio_read_cnt ( wf, time, istatus )
!--------------------------------------

    complex(kind=DP), intent(out), &
         dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: wf
    real(kind=DP), intent(out) :: time
    integer, optional, intent(out) :: istatus

    integer :: input_status

    !read( unit=icnt, iostat=input_status ) time, wf
    if ( rankg == 0 ) then
      read( unit=icnt_t ) time
    end if
    call MPI_Bcast( time, 1, MPI_DOUBLE_PRECISION, 0, &
                    MPI_COMM_WORLD, ierr_mpi )

    read( unit=icnt, iostat=input_status ) wf
    if ( present(istatus) ) then
       istatus = input_status
    endif

  END SUBROUTINE fileio_read_cnt


!--------------------------------------
  SUBROUTINE fileio_write_cnt ( wf, time )
!--------------------------------------

    complex(kind=DP), intent(in), &
         dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: wf
    real(kind=DP), intent(in) :: time

    !rewind ocnt
    !write( unit=ocnt ) time, wf
    if ( rankg == 0 ) then
      write( unit=ocnt_t, pos=1 ) time
    end if
    write( unit=ocnt, pos=1 ) wf

  END SUBROUTINE fileio_write_cnt


!--------------------------------------
  SUBROUTINE fileio_write_fxv ( fout, time )
!--------------------------------------

    complex(kind=DP), intent(in), &
         dimension(-nx:nx,0:ny,1:2*nv,0:nm) :: fout
    real(kind=DP), intent(in) :: time

    !- OUTPUT binary data fxv/"*.fxv.*"
    if ( rankg == 0 ) then
      write( unit=ofxv_t ) time
    end if
    write( unit=ofxv ) fout
    nt_out_fxv = nt_out_fxv + 1

  END SUBROUTINE fileio_write_fxv


!--------------------------------------
  SUBROUTINE fileio_write_phi ( phi, time )
!--------------------------------------

    complex(kind=DP), intent(in), &
         dimension(-nx:nx,0:ny,-nz:nz-1) :: phi
    real(kind=DP), intent(in) :: time

    if ( (ranks /= 0) .OR. (vel_rank /= 0) ) return

    !- OUTPUT binary data phi/"*.phi.*"
    if ( rankg == 0 ) then
      write( unit=ophi_t ) time
    end if
    write( unit=ophi ) phi
    nt_out_phi = nt_out_phi + 1

  END SUBROUTINE fileio_write_phi


!--------------------------------------
  SUBROUTINE fileio_write_Al ( Al, time )
!--------------------------------------

    complex(kind=DP), intent(in), &
         dimension(-nx:nx,0:ny,-nz:nz-1) :: Al
    real(kind=DP), intent(in) :: time

    if ( (ranks /= 0) .OR. (vel_rank /= 0) ) return

    !- OUTPUT binary data phi/"*.Al.*"
    if ( rankg == 0 ) then
      write( unit=oAl_t ) time
    end if
    write( unit=oAl ) Al
    nt_out_Al = nt_out_Al + 1

  END SUBROUTINE fileio_write_Al


!--------------------------------------
  SUBROUTINE fileio_write_mom ( dens, upara, ppara, pperp, qlpara, qlperp, time )
!--------------------------------------

    complex(kind=DP), intent(in), &
         dimension(-nx:nx,0:ny,-nz:nz-1) ::  dens, upara, ppara, pperp, qlpara, qlperp
    real(kind=DP), intent(in) :: time

    if ( vel_rank /= 0 ) return

    !- OUTPUT binary data phi/"*.mom.*"
    if ( rankg == 0 ) then
      write( unit=omom_t ) time
    end if
    write( unit=omom ) dens, upara, ppara, pperp, qlpara, qlperp
    nt_out_mom = nt_out_mom + 1

  END SUBROUTINE fileio_write_mom


!--------------------------------------
  SUBROUTINE fileio_write_trn ( entrpy, fenegy, menegy, peint, pmint, &
                                neint, nmint, dcd, pflux_es, pflux_em, &
                                eflux_es, eflux_em, time )
!--------------------------------------

    real(kind=DP), intent(in), &
         dimension(-nx:nx,0:ny) :: entrpy, fenegy, menegy, peint, pmint, neint, nmint, dcd, &
         pflux_es, pflux_em, eflux_es, eflux_em
    real(kind=DP), intent(in) :: time

    if ( (zsp_rank /= 0) .OR. (vel_rank /= 0) ) return

    !- OUTPUT binary data phi/"*.trn.*"
    if ( rankg == 0 ) then
      write( unit=otrn_t ) time
    end if
    write( unit=otrn ) entrpy, fenegy, menegy,          &
                       peint, pmint, neint, nmint, dcd, &
                       pflux_es, pflux_em, eflux_es, eflux_em
    nt_out_trn = nt_out_trn + 1

  END SUBROUTINE fileio_write_trn


!--------------------------------------
  SUBROUTINE fileio_write_tri ( jkpq_es, jpqk_es, jqkp_es, &
                                jkpq_em, jpqk_em, jqkp_em, time )
!--------------------------------------

    real(kind=DP), intent(in), &
         dimension(-nx:nx,-global_ny:global_ny) :: jkpq_es, jpqk_es, jqkp_es, &
         jkpq_em, jpqk_em, jqkp_em
    real(kind=DP), intent(in) :: time
    integer(kind=DP) :: filesize

    if ( rank /= 0 ) return

    if ( rankg == 0 ) then
       inquire( unit=otri_t, size=filesize ) ! to append
       write( unit=otri_t, pos=filesize+1 ) time
    end if
    inquire( unit=otri, size=filesize ) ! to append
    write( unit=otri, pos=filesize+1  ) jkpq_es, jpqk_es, jqkp_es, jkpq_em, jpqk_em, jqkp_em

  END SUBROUTINE fileio_write_tri


!--------------------------------------
  SUBROUTINE write_zarr_json_head(path)
!--------------------------------------
    ! path: output path for json file (e.g., 'phi.zarr/zz/zarr.json')
    character(*), intent(in) :: path

    open(ojsn, file=path, action="write", status="replace")
    write(ojsn, '(A)') '{'
    write(ojsn, '(A)') '  "attributes": {},'
    write(ojsn, '(A)') '  "zarr_format": 3,'
    write(ojsn, '(A)') '  "consolidated_metadata": null,'
    write(ojsn, '(A)') '  "node_type": "group"'
    write(ojsn, '(A)') '}'
    close(ojsn)

  END SUBROUTINE write_zarr_json_head

!--------------------------------------
  SUBROUTINE write_zarr_json(path, ndims, shape, chunk_shape, dtype, dimnames)
!--------------------------------------
    ! path: output path for json file (e.g., 'phi.zarr/zz/zarr.json')
    ! ndims: number of dimensions
    ! shape(ndims): multi-dimensional shape, C-order (not Fortran-order!!)
    ! chunk_shape(ndims): multi-dimensional chunk, C-order (not Fortran-order!!)
    ! dtype: strings of data type (e.g., 'float64', 'complex128')
    ! dimnames(ndims): names of dimensional axes (e.g., ['zz'], ['t','zz','ky','kx'])
    character(*), intent(in) :: path, dtype
    integer, intent(in) :: ndims, shape(ndims), chunk_shape(ndims)
    character(*), intent(in) :: dimnames(ndims)
    integer :: i

    open(ojsn, file=path, action="write", status="replace")

    write(ojsn, '(A)') '{'

    ! shape
    write(ojsn, '(A)', advance='no') '  "shape": ['
    do i=1, ndims
        if (i > 1) write(ojsn, '(A)', advance='no') ', '
        write(ojsn, '(I0)', advance='no') shape(i)
    end do
    write(ojsn, '(A)') '],'

    ! data_type
    write(ojsn, '(A,A,A)') '  "data_type": "', trim(dtype), '",'

    ! chunk_grid
    write(ojsn, '(A)') '  "chunk_grid": {'
    write(ojsn, '(A)') '    "name": "regular",'
    write(ojsn, '(A)', advance='no') '    "configuration": { "chunk_shape": ['
    do i=1, ndims
        if (i > 1) write(ojsn, '(A)', advance='no') ', '
        write(ojsn, '(I0)', advance='no') chunk_shape(i)
    end do
    write(ojsn, '(A)') '] }'
    write(ojsn, '(A)') '  },'

    ! chunk_key_encoding
    write(ojsn, '(A)') '  "chunk_key_encoding": {'
    write(ojsn, '(A)') '    "name": "default",'
    write(ojsn, '(A)') '    "configuration": { "separator": "/" }'
    write(ojsn, '(A)') '  },'

    ! fill_value
    select case(trim(dtype))
    case('float64', 'float32')
        write(ojsn, '(A)') '  "fill_value": 0.0,'
    case('complex128', 'complex64')
        write(ojsn, '(A)') '  "fill_value": [0.0, 0.0],'
    case('int32', 'int64', 'uint32', 'uint64', 'int16', 'int8')
        write(ojsn, '(A)') '  "fill_value": 0,'
    case default
        write(ojsn, '(A)') '  "fill_value": 0,  ! (default: 0)'
    end select

    ! codecs
    write(ojsn, '(A)') '  "codecs": [{'
    write(ojsn, '(A)') '    "name": "bytes",'
    write(ojsn, '(A)') '    "configuration": { "endian": "little" }'
    write(ojsn, '(A)') '  }],'

    ! attributes
    write(ojsn, '(A)') '  "attributes": {},'

    ! dimension_names
    write(ojsn, '(A)', advance='no') '  "dimension_names": ['
    do i=1, ndims
        if (i > 1) write(ojsn, '(A)', advance='no') ', '
        write(ojsn, '(A)', advance='no') '"' // trim(dimnames(i)) // '"'
    end do
    write(ojsn, '(A)') '],'

    ! 固定項
    write(ojsn, '(A)') '  "zarr_format": 3,'
    write(ojsn, '(A)') '  "node_type": "array",'
    write(ojsn, '(A)') '  "storage_transformers": []'
    write(ojsn, '(A)') '}'

    close(ojsn)

  END SUBROUTINE write_zarr_json

END MODULE GKV_fileio
