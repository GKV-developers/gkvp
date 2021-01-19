MODULE GKV_fileio
!-------------------------------------------------------------------------------
!
!    File I/O interface for NetCDF binary output
!
!    Update history of gkvp_fileio_netcdf.f90
!    --------------
!      gkvp_f0.60 (S. Maeyama, Jan 2021)
!        - NetCDF binary I/O interface by Fujitsu.
!
!-------------------------------------------------------------------------------

  use GKV_header
  use GKV_mpienv

  use netcdf

  implicit  none

  private

  public    fileio_open_icnt, fileio_close_icnt, &
            fileio_open_cnt,  fileio_close_cnt, &
            fileio_open_fxv,  fileio_close_fxv, &
            fileio_open_phi,  fileio_close_phi, &
            fileio_open_Al,   fileio_close_Al, &
            fileio_open_mom,  fileio_close_mom, &
            fileio_open_trn,  fileio_close_trn, &
            fileio_open_tri,  fileio_close_tri, &
            
            fileio_read_cnt,  fileio_write_cnt, &
            fileio_write_fxv, fileio_write_phi, &
            fileio_write_Al,  fileio_write_mom, &
            fileio_write_trn, fileio_write_tri

  integer :: icnt_nc, &
             ocnt_nc, &
             ofxv_nc, &
             ophi_nc, &
             oAl_nc, &
             omom_nc, &
             otrn_nc, &
             otri_nc

  integer :: phi_comm, &
             Al_comm, &
             mom_comm, &
             trn_comm, &
             tri_comm

CONTAINS

!--------------------------------------
  SUBROUTINE fileio_open_icnt ( path )
!--------------------------------------

    character(*), intent(in) :: path

    integer :: ierr_nf90

    character(3)   :: cold

    write( cold,  fmt="(i3.3)" ) inum-1

    ierr_nf90=nf90_open( path=path//"cnt."//cold//".nc",  &
                         mode=NF90_NOWRITE, ncid=icnt_nc, &
                         comm=MPI_COMM_WORLD, info=MPI_INFO_NULL )
    call check_nf90err( ierr_nf90, "nf90_open" )

  END SUBROUTINE fileio_open_icnt

!--------------------------------------
  SUBROUTINE fileio_close_icnt
!--------------------------------------

    integer :: ierr_nf90

    if ( inum > 1 ) then
      ierr_nf90=nf90_close( icnt_nc )
      call check_nf90err( ierr_nf90, "nf90_close" )
    end if

  END SUBROUTINE fileio_close_icnt

 
!--------------------------------------
  SUBROUTINE fileio_open_cnt ( path )
!--------------------------------------

    character(*), intent(in) :: path

    integer :: ierr_nf90

    character(3)   :: cnew

    write( cnew,  fmt="(i3.3)" ) inum

    ierr_nf90=nf90_create( path=path//"cnt."//cnew//".nc", &
                           cmode=IOR(NF90_NETCDF4,NF90_MPIIO), &
                           ncid=ocnt_nc, &
                           comm=MPI_COMM_WORLD, info=MPI_INFO_NULL ) 
    call check_nf90err( ierr_nf90, "nf90_create" )

  END SUBROUTINE fileio_open_cnt

!--------------------------------------
  SUBROUTINE fileio_close_cnt
!--------------------------------------

    integer :: ierr_nf90

    ierr_nf90=nf90_close( ocnt_nc )
    call check_nf90err( ierr_nf90, "nf90_close" )

  END SUBROUTINE fileio_close_cnt


!--------------------------------------
  SUBROUTINE fileio_open_fxv ( path )
!--------------------------------------

    character(*), intent(in) :: path

    integer :: ierr_nf90

    character(3)   :: cnew

    write( cnew,  fmt="(i3.3)" ) inum

    ierr_nf90=nf90_create( path=path//"fxv."//cnew//".nc", &
                           cmode=IOR(NF90_NETCDF4,NF90_MPIIO), &
                           ncid=ofxv_nc, &
                           comm=MPI_COMM_WORLD, info=MPI_INFO_NULL ) 
    call check_nf90err( ierr_nf90, "nf90_create" )

  END SUBROUTINE fileio_open_fxv

!--------------------------------------
  SUBROUTINE fileio_close_fxv
!--------------------------------------

    integer :: ierr_nf90

    ierr_nf90=nf90_close( ofxv_nc )
    call check_nf90err( ierr_nf90, "nf90_close" )

  END SUBROUTINE fileio_close_fxv


!--------------------------------------
  SUBROUTINE fileio_open_phi ( path )
!--------------------------------------

    character(*), intent(in) :: path

    character(3)   :: cnew

    integer :: ierr_nf90
    integer :: counter, i
    integer :: MPI_GROUP_WORLD, new_group, new_comm, ierr_mpi
    integer :: phi_tf(1)
    integer, allocatable :: rank_list(:), phi_tf_list(:)

    allocate(rank_list(nproc))
    allocate(phi_tf_list(nproc))

    if ( ranks == 0 .AND. vel_rank == 0 ) then
      phi_tf(1) = 1
    else
      phi_tf(1) = 0
    end if

    call MPI_Allgather( phi_tf, 1, MPI_INTEGER, phi_tf_list, 1, MPI_INTEGER, &
                        MPI_COMM_WORLD, ierr_mpi)

    counter = 1
    do i = 1, nproc
      if ( phi_tf_list(i) == 1 ) then
        rank_list(counter) = i-1
        counter = counter+1
      end if
    end do

    call MPI_Comm_group( MPI_COMM_WORLD, MPI_GROUP_WORLD, ierr_mpi )
    call MPI_Group_incl( MPI_GROUP_WORLD, counter-1, rank_list, new_group, &
                         ierr_mpi )
    call MPI_Comm_create( MPI_COMM_WORLD, new_group, phi_comm, ierr_mpi )

    write( cnew,  fmt="(i3.3)" ) inum

    if ( ranks == 0 .AND. vel_rank == 0 ) then
      ierr_nf90=nf90_create( path=path//"phi."//cnew//".nc", &
                             cmode=IOR(NF90_NETCDF4,NF90_MPIIO), &
                             ncid=ophi_nc, &
                             comm=phi_comm, info=MPI_INFO_NULL ) 
      call check_nf90err( ierr_nf90, "nf90_create" )
    end if

    deallocate(phi_tf_list, rank_list)

  END SUBROUTINE fileio_open_phi


!--------------------------------------
  SUBROUTINE fileio_close_phi
!--------------------------------------

    integer :: ierr_nf90

    if ( ranks == 0 .AND. vel_rank == 0 ) then
      ierr_nf90=nf90_close( ophi_nc )
      call check_nf90err( ierr_nf90, "nf90_close" )
    end if

  END SUBROUTINE fileio_close_phi


!--------------------------------------
  SUBROUTINE fileio_open_Al ( path )
!--------------------------------------

    character(*), intent(in) :: path

    character(3)   :: cnew

    integer :: ierr_nf90
    integer :: counter, i
    integer :: MPI_GROUP_WORLD, new_group, new_comm, ierr_mpi
    integer :: Al_tf(1)
    integer, allocatable :: rank_list(:), Al_tf_list(:)

    allocate(rank_list(nproc))
    allocate(Al_tf_list(nproc))

    if ( ranks == 0 .AND. vel_rank == 0 ) then
      Al_tf(1) = 1
    else
      Al_tf(1) = 0
    end if

    call MPI_Allgather( Al_tf, 1, MPI_INTEGER, Al_tf_list, 1, MPI_INTEGER, &
                        MPI_COMM_WORLD, ierr_mpi)

    counter = 1
    do i = 1, nproc
      if ( Al_tf_list(i) == 1 ) then
        rank_list(counter) = i-1
        counter = counter+1
      end if
    end do

    call MPI_Comm_group( MPI_COMM_WORLD, MPI_GROUP_WORLD, ierr_mpi )
    call MPI_Group_incl( MPI_GROUP_WORLD, counter-1, rank_list, new_group, &
                         ierr_mpi )
    call MPI_Comm_create( MPI_COMM_WORLD, new_group, Al_comm, ierr_mpi )

    write( cnew,  fmt="(i3.3)" ) inum

    if ( ranks == 0 .AND. vel_rank == 0 ) then
      ierr_nf90=nf90_create( path=path//"Al."//cnew//".nc", &
                             cmode=IOR(NF90_NETCDF4,NF90_MPIIO), &
                             ncid=oAl_nc, &
                             comm=Al_comm, info=MPI_INFO_NULL )
      call check_nf90err( ierr_nf90, "nf90_create" )
    end if

    deallocate(Al_tf_list, rank_list)

  END SUBROUTINE fileio_open_Al


!--------------------------------------
  SUBROUTINE fileio_close_Al
!--------------------------------------

    integer :: ierr_nf90

    if ( ranks == 0 .AND. vel_rank == 0 ) then
      ierr_nf90=nf90_close( oAl_nc )
      call check_nf90err( ierr_nf90, "nf90_close" )
    end if

  END SUBROUTINE fileio_close_Al


!--------------------------------------
  SUBROUTINE fileio_open_mom ( path )
!--------------------------------------

    character(*), intent(in) :: path

    character(3)   :: cnew

    integer :: ierr_nf90
    integer :: counter, i
    integer :: MPI_GROUP_WORLD, new_group, new_comm, ierr_mpi
    integer :: mom_tf(1)
    integer, allocatable :: rank_list(:), mom_tf_list(:)

    allocate(rank_list(nproc))
    allocate(mom_tf_list(nproc))
    rank_list = 0
    mom_tf_list = 0

    if ( vel_rank == 0 ) then
      mom_tf(1) = 1
    else
      mom_tf(1) = 0
    end if

    call MPI_Allgather( mom_tf, 1, MPI_INTEGER, mom_tf_list, 1, MPI_INTEGER, &
                        MPI_COMM_WORLD, ierr_mpi)

    counter = 1
    do i = 1, nproc
      if ( mom_tf_list(i) == 1 ) then
        rank_list(counter) = i-1
        counter = counter+1
      end if
    end do

    call MPI_Comm_group( MPI_COMM_WORLD, MPI_GROUP_WORLD, ierr_mpi )
    call MPI_Group_incl( MPI_GROUP_WORLD, counter-1, rank_list, new_group, &
                         ierr_mpi )
    call MPI_Comm_create( MPI_COMM_WORLD, new_group, mom_comm, ierr_mpi )

    write( cnew,  fmt="(i3.3)" ) inum

    if ( vel_rank == 0 ) then
      ierr_nf90=nf90_create( path=path//"mom."//cnew//".nc", &
                             cmode=IOR(NF90_NETCDF4,NF90_MPIIO), &
                             ncid=omom_nc, &
                             comm=mom_comm, info=MPI_INFO_NULL )
      call check_nf90err( ierr_nf90, "nf90_create" )
    end if

    deallocate(mom_tf_list, rank_list)

  END SUBROUTINE fileio_open_mom


!--------------------------------------
  SUBROUTINE fileio_close_mom
!--------------------------------------

    integer :: ierr_nf90

    if ( vel_rank == 0 ) then
      ierr_nf90=nf90_close( omom_nc )
      call check_nf90err( ierr_nf90, "nf90_close" )
    end if

  END SUBROUTINE fileio_close_mom


!--------------------------------------
  SUBROUTINE fileio_open_trn( path )
!--------------------------------------

    character(*), intent(in) :: path

    character(3)   :: cnew

    integer :: ierr_nf90
    integer :: counter, i
    integer :: MPI_GROUP_WORLD, new_group, new_comm, ierr_mpi
    integer :: trn_tf(1)
    integer, allocatable :: rank_list(:), trn_tf_list(:)

    allocate(rank_list(nproc))
    allocate(trn_tf_list(nproc))

    if ( zsp_rank == 0 .AND. vel_rank == 0 ) then
      trn_tf(1) = 1
    else
      trn_tf(1) = 0
    end if

    call MPI_Allgather( trn_tf, 1, MPI_INTEGER, trn_tf_list, 1, MPI_INTEGER, &
                        MPI_COMM_WORLD, ierr_mpi)

    counter = 1
    do i = 1, nproc
      if ( trn_tf_list(i) == 1 ) then
        rank_list(counter) = i-1
        counter = counter+1
      end if
    end do

    call MPI_Comm_group( MPI_COMM_WORLD, MPI_GROUP_WORLD, ierr_mpi )
    call MPI_Group_incl( MPI_GROUP_WORLD, counter-1, rank_list, new_group, &
                         ierr_mpi )
    call MPI_Comm_create( MPI_COMM_WORLD, new_group, trn_comm, ierr_mpi )

    write( cnew,  fmt="(i3.3)" ) inum

    if ( zsp_rank == 0 .AND. vel_rank == 0 ) then
      ierr_nf90=nf90_create( path=path//"trn."//cnew//".nc", &
                             cmode=IOR(NF90_NETCDF4,NF90_MPIIO), &
                             ncid=otrn_nc, &
                             comm=trn_comm, info=MPI_INFO_NULL ) 
      call check_nf90err( ierr_nf90, "nf90_create" )
    end if

    deallocate(trn_tf_list, rank_list)

  END SUBROUTINE fileio_open_trn


!--------------------------------------
  SUBROUTINE fileio_close_trn
!--------------------------------------

    integer :: ierr_nf90

    if ( zsp_rank == 0 .AND. vel_rank == 0 ) then
      ierr_nf90=nf90_close( otrn_nc )
      call check_nf90err( ierr_nf90, "nf90_close" )
    end if

  END SUBROUTINE fileio_close_trn


!!--------------------------------------
!  SUBROUTINE fileio_open_tri ( path, cmx, cmy, replace )
!!--------------------------------------
!
!    character(*), intent(in) :: path
!    character(*), intent(in) :: cmx, cmy
!    logical, intent(in) :: replace
!
!    character(3)   :: cnew
!
!    integer :: ierr_nf90
!    integer :: counter, i
!    integer :: MPI_GROUP_WORLD, new_group, new_comm, ierr_mpi
!    integer :: tri_tf(1)
!    integer, allocatable :: rank_list(:), tri_tf_list(:)
!
!    allocate(rank_list(nproc))
!    allocate(tri_tf_list(nproc))
!
!    if ( rank == 0 ) then
!      tri_tf(1) = 1
!    else
!      tri_tf(1) = 0
!    end if
!
!    call MPI_Allgather( tri_tf, 1, MPI_INTEGER, tri_tf_list, 1, MPI_INTEGER, &
!                        MPI_COMM_WORLD, ierr_mpi)
!
!    counter = 1
!    do i = 1, nproc
!      if ( tri_tf_list(i) == 1 ) then
!        rank_list(counter) = i-1
!        counter = counter+1
!      end if
!    end do
!
!    call MPI_Comm_group( MPI_COMM_WORLD, MPI_GROUP_WORLD, ierr_mpi )
!    call MPI_Group_incl( MPI_GROUP_WORLD, counter-1, rank_list, new_group, &
!                         ierr_mpi )
!    call MPI_Comm_create( MPI_COMM_WORLD, new_group, tri_comm, ierr_mpi )
!
!    write( cnew,  fmt="(i3.3)" ) inum
!
!    if ( rank == 0 ) then
!      if ( replace ) then
!        ierr_nf90=nf90_create( &
!                path=path//"mx"//cmx//"my"//cmy//".tri."//cnew//".nc", &
!                cmode=IOR(NF90_NETCDF4,NF90_MPIIO), ncid=otri_nc, &
!                comm=tri_comm, info=MPI_INFO_NULL ) 
!        call check_nf90err( ierr_nf90, "nf90_create" )
!      else
!        ierr_nf90=nf90_open( &
!                path=path//"mx"//cmx//"my"//cmy//".tri."//cnew//".nc", &
!                mode=NF90_WRITE, ncid=otri_nc, &
!                comm=tri_comm, info=MPI_INFO_NULL ) 
!        call check_nf90err( ierr_nf90, "nf90_open" )
!      end if
!    end if
!
!    deallocate(tri_tf_list, rank_list)
!
!  END SUBROUTINE fileio_open_tri

!--------------------------------------
  SUBROUTINE fileio_open_tri ( path, cmx, cmy, replace )
!--------------------------------------

    character(*), intent(in) :: path
    character(*), intent(in) :: cmx, cmy
    logical, intent(in) :: replace

    character(3)   :: cnew

    integer :: ierr_nf90
    integer :: counter, i
    integer :: MPI_GROUP_WORLD, new_group, new_comm, ierr_mpi
    integer :: tri_tf(1)
    integer, allocatable :: rank_list(:), tri_tf_list(:)

    write( cnew,  fmt="(i3.3)" ) inum

    if ( replace ) then ! Initial creation of NetCDF files

      allocate(rank_list(nproc))
      allocate(tri_tf_list(nproc))

      if ( rank == 0 ) then
        tri_tf(1) = 1
      else
        tri_tf(1) = 0
      end if
      call MPI_Allgather( tri_tf, 1, MPI_INTEGER, tri_tf_list, 1, MPI_INTEGER, &
                          MPI_COMM_WORLD, ierr_mpi)
      counter = 1
      do i = 1, nproc
        if ( tri_tf_list(i) == 1 ) then
          rank_list(counter) = i-1
          counter = counter+1
        end if
      end do
      call MPI_Comm_group( MPI_COMM_WORLD, MPI_GROUP_WORLD, ierr_mpi )
      call MPI_Group_incl( MPI_GROUP_WORLD, counter-1, rank_list, new_group, &
                           ierr_mpi )
      call MPI_Comm_create( MPI_COMM_WORLD, new_group, tri_comm, ierr_mpi )

      if ( rank == 0 ) then
        ierr_nf90=nf90_create( &
                path=path//"mx"//cmx//"my"//cmy//".tri."//cnew//".nc", &
                cmode=IOR(NF90_NETCDF4,NF90_MPIIO), ncid=otri_nc, &
                comm=tri_comm, info=MPI_INFO_NULL )
        call check_nf90err( ierr_nf90, "nf90_create" )
      end if

      deallocate(tri_tf_list, rank_list)

    else ! For appending data to NetCDF files

      if ( rank == 0 ) then
        ierr_nf90=nf90_open( &
                path=path//"mx"//cmx//"my"//cmy//".tri."//cnew//".nc", &
                mode=NF90_WRITE, ncid=otri_nc, &
                comm=tri_comm, info=MPI_INFO_NULL )
        call check_nf90err( ierr_nf90, "nf90_open" )
      end if

    end if

  END SUBROUTINE fileio_open_tri


!--------------------------------------
  SUBROUTINE fileio_close_tri
!--------------------------------------

    integer :: ierr_nf90

    if ( rank == 0 ) then
      ierr_nf90=nf90_close( otri_nc )
      call check_nf90err( ierr_nf90, "nf90_close" )
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
    integer :: ierr_nf90
    integer :: varid_tt, varid_recnt, varid_imcnt
    integer(kind=4) :: start_time(1),  count_time(1), &
                       start_cnt(1:7), count_cnt(1:7)
    integer :: start_ny(1), count_ny(1)
    integer :: start_nz(1), count_nz(1)
    integer :: start_nv(1), count_nv(1)
    integer :: start_nm(1), count_nm(1)
    integer :: is
    integer :: ny_st, ny_end
    integer :: t_count
    real(kind=DP), dimension(1) :: t_value
    real(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: recnt
    real(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: imcnt

    integer :: ndims, nvars, ngatts, unlimdimid

    != Information about an open netCDF dataset
    ierr_nf90=nf90_inquire( icnt_nc, ndims, nvars, ngatts, unlimdimid )
    call check_nf90err( ierr_nf90, "nf90_inquire" )

    ierr_nf90=nf90_inq_varid(icnt_nc, "t", varid_tt)
    ierr_nf90=nf90_inq_varid(icnt_nc, "recnt", varid_recnt)
    ierr_nf90=nf90_inq_varid(icnt_nc, "imcnt", varid_imcnt)
    call check_nf90err(ierr_nf90, "nf90_inq_varid")

    ierr_nf90=nf90_inquire_dimension(icnt_nc, unlimdimid, len=t_count)
    call check_nf90err(ierr_nf90, "nf90_inquire_dimension")

    count_ny=nsize_y
    start_ny=ist_y_g+1
    count_nz=nz*2
    start_nz=count_nz*rankz+1
    count_nv=nv*2
    start_nv=count_nv*rankv+1
    count_nm=nm+1
    start_nm=count_nm*rankm+1

    count_time(:) = 1
    start_time(:) = t_count
    ierr_nf90=nf90_get_var( ncid=icnt_nc, varid=varid_tt, values=t_value, &
                            start=start_time, count=count_time )
    time = t_value(1)

    ny_st = 0
    ny_end = count_ny(1)-1
    count_cnt(:) = int((/ 2*nx+1,count_ny,count_nz,count_nv,count_nm,1,1 /), kind=4)
    start_cnt(:) = (/ 1,start_ny,start_nz,start_nv,start_nm,1+ranks,t_count /)
    ierr_nf90=nf90_get_var(icnt_nc, varid_recnt, &
                           values=recnt(:,ny_st:ny_end,:,:,:), &
                           start=start_cnt, count=count_cnt )
    ierr_nf90=nf90_get_var(icnt_nc, varid_imcnt, &
                           values=imcnt(:,ny_st:ny_end,:,:,:), &
                           start=start_cnt, count=count_cnt )
    call check_nf90err(ierr_nf90, "nf90_get_var")

    wf = cmplx(recnt,imcnt)

    istatus = ierr_nf90

  END SUBROUTINE fileio_read_cnt


!--------------------------------------
  SUBROUTINE fileio_write_cnt ( wf, time )
!--------------------------------------

    complex(kind=DP), intent(in), &
         dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: wf
    real(kind=DP), intent(in) :: time

    integer :: dimids(1:7), ierr_nf90
    integer :: varid_kx, varid_ky, varid_zz, varid_vl, varid_mu, varid_is, &
               varid_tt, varid_recnt, varid_imcnt
    integer(kind=4) :: start_time(1),  count_time(1), &
                       start_cnt(1:7), count_cnt(1:7)
    integer :: start_ny(1), count_ny(1)
    integer :: start_nz(1), count_nz(1)
    integer :: start_nv(1), count_nv(1)
    integer :: start_nm(1), count_nm(1)
    integer :: is
    integer :: ny_st, ny_end

    integer :: ndims, nvars, ngatts, unlimdimid
    integer :: ocnt_t_count

    != Information about an open netCDF dataset
    ierr_nf90=nf90_inquire( ocnt_nc, ndims, nvars, ngatts, unlimdimid )
    call check_nf90err( ierr_nf90, "nf90_inquire" )

    if ( ndims == 0 ) then
      != Define dimensions in file
      ierr_nf90=nf90_def_dim( ncid=ocnt_nc, name="kx", &
                              len=int(2*nx+1,kind=4),      dimid=dimids(1) )
      ierr_nf90=nf90_def_dim( ncid=ocnt_nc, name="ky", &
                              len=int(global_ny+1,kind=4), dimid=dimids(2) )
      ierr_nf90=nf90_def_dim( ncid=ocnt_nc, name="zz", &
                              len=int(2*global_nz,kind=4), dimid=dimids(3) )
      ierr_nf90=nf90_def_dim( ncid=ocnt_nc, name="vl", &
                              len=int(2*global_nv,kind=4), dimid=dimids(4) )
      ierr_nf90=nf90_def_dim( ncid=ocnt_nc, name="mu", &
                              len=int(global_nm+1,kind=4), dimid=dimids(5) )
      ierr_nf90=nf90_def_dim( ncid=ocnt_nc, name="is", &
                              len=int(ns,kind=4),          dimid=dimids(6) )
      ierr_nf90=nf90_def_dim( ncid=ocnt_nc, name="t",  &
                              len=NF90_UNLIMITED,          dimid=dimids(7) )
      call check_nf90err( ierr_nf90, "nf90_def_dim" )

      != Define variables in file
      ierr_nf90=nf90_def_var( ncid=ocnt_nc, name="kx",    xtype=NF90_DOUBLE, &
                              dimids=dimids(1),   varid=varid_kx )
      ierr_nf90=nf90_def_var( ncid=ocnt_nc, name="ky",    xtype=NF90_DOUBLE, &
                              dimids=dimids(2),   varid=varid_ky )
      ierr_nf90=nf90_def_var( ncid=ocnt_nc, name="zz",    xtype=NF90_DOUBLE, &
                              dimids=dimids(3),   varid=varid_zz )
      ierr_nf90=nf90_def_var( ncid=ocnt_nc, name="vl",    xtype=NF90_DOUBLE, &
                              dimids=dimids(4),   varid=varid_vl )
      ierr_nf90=nf90_def_var( ncid=ocnt_nc, name="mu",    xtype=NF90_DOUBLE, &
                              dimids=dimids(5),   varid=varid_mu )
      ierr_nf90=nf90_def_var( ncid=ocnt_nc, name="is",    xtype=NF90_INT, &
                              dimids=dimids(6),   varid=varid_is )
      ierr_nf90=nf90_def_var( ncid=ocnt_nc, name="t",     xtype=NF90_DOUBLE, &
                              dimids=dimids(7),   varid=varid_tt )
      ierr_nf90=nf90_def_var( ncid=ocnt_nc, name="recnt", xtype=NF90_DOUBLE, &
                              dimids=dimids(1:7), varid=varid_recnt )
      ierr_nf90=nf90_def_var( ncid=ocnt_nc, name="imcnt", xtype=NF90_DOUBLE, &
                              dimids=dimids(1:7), varid=varid_imcnt )
      call check_nf90err( ierr_nf90, "nf90_def_var" )

      != End of definition of file
      ierr_nf90=nf90_enddef( ncid=ocnt_nc )
      call check_nf90err( ierr_nf90, "nf90_enddef" )

      ocnt_t_count = 0
    else
      ierr_nf90=nf90_inq_varid( ocnt_nc, "t",     varid_tt )
      ierr_nf90=nf90_inq_varid( ocnt_nc, "recnt", varid_recnt )
      ierr_nf90=nf90_inq_varid( ocnt_nc, "imcnt", varid_imcnt )
      call check_nf90err( ierr_nf90, "nf90_inq_varid" )

      ierr_nf90=nf90_inquire_dimension(ocnt_nc, unlimdimid, len=ocnt_t_count)
      call check_nf90err(ierr_nf90, "nf90_inquire_dimension")
    end if

    != Parallel data access type
    ierr_nf90=nf90_var_par_access( ncid=ocnt_nc, varid=varid_kx, &
                                   access=NF90_INDEPENDENT )
    ierr_nf90=nf90_var_par_access( ncid=ocnt_nc, varid=varid_ky, &
                                   access=NF90_INDEPENDENT )
    ierr_nf90=nf90_var_par_access( ncid=ocnt_nc, varid=varid_zz, &
                                   access=NF90_INDEPENDENT )
    ierr_nf90=nf90_var_par_access( ncid=ocnt_nc, varid=varid_vl, &
                                   access=NF90_INDEPENDENT )
    ierr_nf90=nf90_var_par_access( ncid=ocnt_nc, varid=varid_mu, &
                                   access=NF90_INDEPENDENT )
    ierr_nf90=nf90_var_par_access( ncid=ocnt_nc, varid=varid_is, &
                                   access=NF90_INDEPENDENT )
    ierr_nf90=nf90_var_par_access( ncid=ocnt_nc, varid=varid_tt, &
                                   access=NF90_COLLECTIVE )
    ierr_nf90=nf90_var_par_access( ncid=ocnt_nc, varid=varid_recnt, &
                                   access=NF90_COLLECTIVE )
    ierr_nf90=nf90_var_par_access( ncid=ocnt_nc, varid=varid_imcnt, &
                                   access=NF90_COLLECTIVE )
    call check_nf90err( ierr_nf90, "nf90_var_par_access" )

    count_ny=nsize_y
    start_ny=ist_y_g+1
    count_nz=nz*2
    start_nz=count_nz*rankz+1
    count_nv=nv*2
    start_nv=count_nv*rankv+1
    count_nm=nm+1
    start_nm=count_nm*rankm+1

    if ( ndims == 0 ) then
      != Write variables: static coordinates x and y
      ierr_nf90=nf90_put_var( ncid=ocnt_nc, varid=varid_kx, &
                              values=kx(-nx:nx) )
      ierr_nf90=nf90_put_var( ncid=ocnt_nc, varid=varid_ky, &
                              values=ky(0:iend_y), &
                              start=start_ny, count=count_ny )
      ierr_nf90=nf90_put_var( ncid=ocnt_nc, varid=varid_zz, &
                              values=zz(-nz:nz-1), &
                              start=start_nz, count=count_nz )
      ierr_nf90=nf90_put_var( ncid=ocnt_nc, varid=varid_vl, &
                              values=vl(1:2*nv), &
                              start=start_nv, count=count_nv )
      ierr_nf90=nf90_put_var( ncid=ocnt_nc, varid=varid_mu, &
                              values=mu(0:nm), &
                              start=start_nm, count=count_nm )
      ierr_nf90=nf90_put_var( ncid=ocnt_nc, varid=varid_is, &
                              values=(/ (is, is=0,ns-1) /) )
      call check_nf90err( ierr_nf90, "nf90_putvar" )
    end if

    ocnt_t_count = ocnt_t_count+1
    count_time(:) = 1
    start_time(:) = ocnt_t_count
    ierr_nf90=nf90_put_var( ncid=ocnt_nc, varid=varid_tt, values=(/time/), &
                            start=start_time, count=count_time )

    ny_st = 0
    ny_end = count_ny(1)-1
    count_cnt(:) = int((/ 2*nx+1,count_ny,count_nz,count_nv,count_nm,1,1 /), kind=4)
    start_cnt(:) = (/ 1,start_ny,start_nz,start_nv,start_nm,1+ranks,ocnt_t_count /)
    ierr_nf90=nf90_put_var( ncid=ocnt_nc, varid=varid_recnt, &
                            values=dble(wf(:,ny_st:ny_end,:,:,:)), &
                            start=start_cnt, count=count_cnt )
    ierr_nf90=nf90_put_var( ncid=ocnt_nc, varid=varid_imcnt, &
                            values=aimag(wf(:,ny_st:ny_end,:,:,:)), &
                            start=start_cnt, count=count_cnt )
    call check_nf90err( ierr_nf90, "nf90_putvar" )

  END SUBROUTINE fileio_write_cnt


!--------------------------------------
  SUBROUTINE fileio_write_fxv ( fout, time )
!--------------------------------------

    complex(kind=DP), intent(in), &
         dimension(-nx:nx,0:ny,1:2*nv,0:nm) :: fout
    real(kind=DP), intent(in) :: time

    integer :: dimids(1:7), ierr_nf90
    integer :: varid_kx, varid_ky, varid_zz, varid_vl, varid_mu, varid_is, &
               varid_tt, varid_refxv, varid_imfxv
    integer(kind=4) :: start_time(1),  count_time(1), &
                       start_fxv(1:7), count_fxv(1:7)
    integer :: start_ny(1), count_ny(1)
    integer :: start_nz(1), count_nz(1)
    integer :: start_nv(1), count_nv(1)
    integer :: start_nm(1), count_nm(1)
    integer :: is
    integer :: ny_st, ny_end
    integer :: ndims, nvars, ngatts, unlimdimid
    integer :: ofxv_t_count

    != Information about an open netCDF dataset
    ierr_nf90=nf90_inquire( ofxv_nc, ndims, nvars, ngatts, unlimdimid )
    call check_nf90err( ierr_nf90, "nf90_inquire" )

    if ( ndims == 0 ) then
      != Define dimensions in file
      ierr_nf90=nf90_def_dim( ncid=ofxv_nc, name="kx", &
                              len=int(2*nx+1,kind=4),      dimid=dimids(1) )
      ierr_nf90=nf90_def_dim( ncid=ofxv_nc, name="ky", &
                              len=int(global_ny+1,kind=4), dimid=dimids(2) )
      ierr_nf90=nf90_def_dim( ncid=ofxv_nc, name="zz", &
                              len=int(nprocz,kind=4),      dimid=dimids(3) )
      ierr_nf90=nf90_def_dim( ncid=ofxv_nc, name="vl", &
                              len=int(2*global_nv,kind=4), dimid=dimids(4) )
      ierr_nf90=nf90_def_dim( ncid=ofxv_nc, name="mu", &
                              len=int(global_nm+1,kind=4), dimid=dimids(5) )
      ierr_nf90=nf90_def_dim( ncid=ofxv_nc, name="is", &
                              len=int(ns,kind=4),          dimid=dimids(6) )
      ierr_nf90=nf90_def_dim( ncid=ofxv_nc, name="t",  &
                              len=NF90_UNLIMITED,          dimid=dimids(7) )
      call check_nf90err( ierr_nf90, "nf90_def_dim" )

      != Define variables in file
      ierr_nf90=nf90_def_var( ncid=ofxv_nc, name="kx",    xtype=NF90_DOUBLE, &
                              dimids=dimids(1),   varid=varid_kx )
      ierr_nf90=nf90_def_var( ncid=ofxv_nc, name="ky",    xtype=NF90_DOUBLE, &
                              dimids=dimids(2),   varid=varid_ky )
      ierr_nf90=nf90_def_var( ncid=ofxv_nc, name="zz",    xtype=NF90_DOUBLE, &
                              dimids=dimids(3),   varid=varid_zz )
      ierr_nf90=nf90_def_var( ncid=ofxv_nc, name="vl",    xtype=NF90_DOUBLE, &
                              dimids=dimids(4),   varid=varid_vl )
      ierr_nf90=nf90_def_var( ncid=ofxv_nc, name="mu",    xtype=NF90_DOUBLE, &
                              dimids=dimids(5),   varid=varid_mu )
      ierr_nf90=nf90_def_var( ncid=ofxv_nc, name="is",    xtype=NF90_INT, &
                              dimids=dimids(6),   varid=varid_is )
      ierr_nf90=nf90_def_var( ncid=ofxv_nc, name="t",     xtype=NF90_DOUBLE, &
                              dimids=dimids(7),   varid=varid_tt )
      ierr_nf90=nf90_def_var( ncid=ofxv_nc, name="refxv", xtype=NF90_DOUBLE, &
                              dimids=dimids(1:7), varid=varid_refxv )
      ierr_nf90=nf90_def_var( ncid=ofxv_nc, name="imfxv", xtype=NF90_DOUBLE, &
                              dimids=dimids(1:7), varid=varid_imfxv )
      call check_nf90err( ierr_nf90, "nf90_def_var" )

      != End of definition of file
      ierr_nf90=nf90_enddef( ncid=ofxv_nc )
      call check_nf90err( ierr_nf90, "nf90_enddef" )

      ofxv_t_count = 0
    else
      ierr_nf90=nf90_inq_varid( ofxv_nc, "t",     varid_tt )
      ierr_nf90=nf90_inq_varid( ofxv_nc, "refxv", varid_refxv )
      ierr_nf90=nf90_inq_varid( ofxv_nc, "imfxv", varid_imfxv )
      call check_nf90err( ierr_nf90, "nf90_inq_varid" )

      ierr_nf90=nf90_inquire_dimension(ofxv_nc, unlimdimid, len=ofxv_t_count)
      call check_nf90err(ierr_nf90, "nf90_inquire_dimension")
    end if

    != Parallel data access type
    ierr_nf90=nf90_var_par_access( ncid=ofxv_nc, varid=varid_kx, &
                                   access=NF90_INDEPENDENT )
    ierr_nf90=nf90_var_par_access( ncid=ofxv_nc, varid=varid_ky, &
                                   access=NF90_INDEPENDENT )
    ierr_nf90=nf90_var_par_access( ncid=ofxv_nc, varid=varid_zz, &
                                   access=NF90_INDEPENDENT )
    ierr_nf90=nf90_var_par_access( ncid=ofxv_nc, varid=varid_vl, &
                                   access=NF90_INDEPENDENT )
    ierr_nf90=nf90_var_par_access( ncid=ofxv_nc, varid=varid_mu, &
                                   access=NF90_INDEPENDENT )
    ierr_nf90=nf90_var_par_access( ncid=ofxv_nc, varid=varid_is, &
                                   access=NF90_INDEPENDENT )
    ierr_nf90=nf90_var_par_access( ncid=ofxv_nc, varid=varid_tt, &
                                   access=NF90_COLLECTIVE )
    ierr_nf90=nf90_var_par_access( ncid=ofxv_nc, varid=varid_refxv, &
                                   access=NF90_COLLECTIVE )
    ierr_nf90=nf90_var_par_access( ncid=ofxv_nc, varid=varid_imfxv, &
                                   access=NF90_COLLECTIVE )
    call check_nf90err( ierr_nf90, "nf90_var_par_access" )

    count_ny=nsize_y
    start_ny=ist_y_g+1
    count_nz=1
    start_nz=rankz+1
    count_nv=nv*2
    start_nv=count_nv*rankv+1
    count_nm=nm+1
    start_nm=count_nm*rankm+1

    != Write variables: static coordinates x and y
    if ( ndims == 0 ) then
      ierr_nf90=nf90_put_var( ncid=ofxv_nc, varid=varid_kx, &
                              values=kx(-nx:nx) )
      ierr_nf90=nf90_put_var( ncid=ofxv_nc, varid=varid_ky, &
                              values=ky(0:iend_y), &
                              start=start_ny, count=count_ny )
      ierr_nf90=nf90_put_var( ncid=ofxv_nc, varid=varid_zz, &
                              values=zz(-nz:nz-1:2*nz), &
                              start=start_nz, count=count_nz )
      ierr_nf90=nf90_put_var( ncid=ofxv_nc, varid=varid_vl, &
                              values=vl(1:2*nv), &
                              start=start_nv, count=count_nv )
      ierr_nf90=nf90_put_var( ncid=ofxv_nc, varid=varid_mu, &
                              values=mu(0:nm), &
                              start=start_nm, count=count_nm )
      ierr_nf90=nf90_put_var( ncid=ofxv_nc, varid=varid_is, &
                              values=(/ (is, is=0,ns-1) /) )
      call check_nf90err( ierr_nf90, "nf90_putvar" )
    end if

    ofxv_t_count = ofxv_t_count+1
    count_time(:) = 1
    start_time(:) = ofxv_t_count
    ierr_nf90=nf90_put_var( ncid=ofxv_nc, varid=varid_tt, values=(/time/), &
                            start=start_time, count=count_time )

    ny_st = 0
    ny_end = count_ny(1)-1
    count_fxv(:) = int((/ 2*nx+1,count_ny,count_nz,count_nv,count_nm,1,1 /), kind=4)
    start_fxv(:) = (/ 1,start_ny,start_nz,start_nv,start_nm,1+ranks,ofxv_t_count /)
    ierr_nf90=nf90_put_var( ncid=ofxv_nc, varid=varid_refxv, &
                            values=dble(fout(:,ny_st:ny_end,:,:)), &
                            start=start_fxv, count=count_fxv )
    ierr_nf90=nf90_put_var( ncid=ofxv_nc, varid=varid_imfxv, &
                            values=aimag(fout(:,ny_st:ny_end,:,:)), &
                            start=start_fxv, count=count_fxv )
    call check_nf90err( ierr_nf90, "nf90_putvar" )

  END SUBROUTINE fileio_write_fxv


!--------------------------------------
  SUBROUTINE fileio_write_phi ( phi, time )
!--------------------------------------

    complex(kind=DP), intent(in), &
         dimension(-nx:nx,0:ny,-nz:nz-1) :: phi
    real(kind=DP), intent(in) :: time

    integer :: dimids(1:4), ierr_nf90
    integer :: varid_kx, varid_ky, varid_zz, varid_is, &
               varid_tt, varid_rephi, varid_imphi
    integer(kind=4) :: start_time(1), count_time(1), start_phi(1:4), count_phi(1:4)
    integer :: start_ny(1), count_ny(1)
    integer :: start_nz(1), count_nz(1)
    integer :: is
    integer :: ny_st, ny_end
    integer :: ndims, nvars, ngatts, unlimdimid
    integer :: ophi_t_count

    != Write conditions check
    if ( (ranks /= 0) .OR. (vel_rank /= 0) ) then
      return
    end if

    != Information about an open netCDF dataset
    ierr_nf90=nf90_inquire( ophi_nc, ndims, nvars, ngatts, unlimdimid )
    call check_nf90err( ierr_nf90, "nf90_inquire" )

    if ( ndims == 0 ) then
      != Define dimensions in file
      ierr_nf90=nf90_def_dim( ncid=ophi_nc, name="kx", &
                              len=int(2*nx+1,kind=4),      dimid=dimids(1) )
      ierr_nf90=nf90_def_dim( ncid=ophi_nc, name="ky", &
                              len=int(global_ny+1,kind=4), dimid=dimids(2) )
      ierr_nf90=nf90_def_dim( ncid=ophi_nc, name="zz", &
                              len=int(2*global_nz,kind=4), dimid=dimids(3) )
      ierr_nf90=nf90_def_dim( ncid=ophi_nc, name="t",  &
                              len=NF90_UNLIMITED,          dimid=dimids(4) )
      call check_nf90err( ierr_nf90, "nf90_def_dim" )

      != Define variables in file
      ierr_nf90=nf90_def_var( ncid=ophi_nc, name="kx",    xtype=NF90_DOUBLE, &
                              dimids=dimids(1),   varid=varid_kx )
      ierr_nf90=nf90_def_var( ncid=ophi_nc, name="ky",    xtype=NF90_DOUBLE, &
                              dimids=dimids(2),   varid=varid_ky )
      ierr_nf90=nf90_def_var( ncid=ophi_nc, name="zz",    xtype=NF90_DOUBLE, &
                              dimids=dimids(3),   varid=varid_zz )
      ierr_nf90=nf90_def_var( ncid=ophi_nc, name="t",     xtype=NF90_DOUBLE, &
                              dimids=dimids(4),   varid=varid_tt )
      ierr_nf90=nf90_def_var( ncid=ophi_nc, name="rephi", xtype=NF90_DOUBLE, &
                              dimids=dimids(1:4), varid=varid_rephi )
      ierr_nf90=nf90_def_var( ncid=ophi_nc, name="imphi", xtype=NF90_DOUBLE, &
                              dimids=dimids(1:4), varid=varid_imphi )
      call check_nf90err( ierr_nf90, "nf90_def_var" )

      != End of definition of file
      ierr_nf90=nf90_enddef( ncid=ophi_nc )
      call check_nf90err( ierr_nf90, "nf90_enddef" )

      ophi_t_count = 0
    else
      ierr_nf90=nf90_inq_varid( ophi_nc, "t",     varid_tt )
      ierr_nf90=nf90_inq_varid( ophi_nc, "rephi", varid_rephi )
      ierr_nf90=nf90_inq_varid( ophi_nc, "imphi", varid_imphi )
      call check_nf90err( ierr_nf90, "nf90_inq_varid" )

      ierr_nf90=nf90_inquire_dimension(ophi_nc, unlimdimid, len=ophi_t_count)
      call check_nf90err(ierr_nf90, "nf90_inquire_dimension")
    end if

    != Parallel data access type
    ierr_nf90=nf90_var_par_access( ncid=ophi_nc, varid=varid_kx, &
                                   access=NF90_INDEPENDENT )
    ierr_nf90=nf90_var_par_access( ncid=ophi_nc, varid=varid_ky, &
                                   access=NF90_INDEPENDENT )
    ierr_nf90=nf90_var_par_access( ncid=ophi_nc, varid=varid_zz, &
                                   access=NF90_INDEPENDENT )
    ierr_nf90=nf90_var_par_access( ncid=ophi_nc, varid=varid_tt, &
                                   access=NF90_COLLECTIVE )
    ierr_nf90=nf90_var_par_access( ncid=ophi_nc, varid=varid_rephi, &
                                   access=NF90_COLLECTIVE )
    ierr_nf90=nf90_var_par_access( ncid=ophi_nc, varid=varid_imphi, &
                                   access=NF90_COLLECTIVE )
    call check_nf90err( ierr_nf90, "nf90_var_par_access" )

    count_ny=nsize_y
    start_ny=ist_y_g+1
    count_nz=nz*2
    start_nz=count_nz*rankz+1

    != Write variables: static coordinates x and y
    if ( ndims == 0 ) then
      ierr_nf90=nf90_put_var( ncid=ophi_nc, varid=varid_kx, &
                              values=kx(-nx:nx) )
      ierr_nf90=nf90_put_var( ncid=ophi_nc, varid=varid_ky, &
                              values=ky(0:iend_y), &
                              start=start_ny, count=count_ny )
      ierr_nf90=nf90_put_var( ncid=ophi_nc, varid=varid_zz, &
                              values=zz(-nz:nz-1), &
                              start=start_nz, count=count_nz )
      call check_nf90err( ierr_nf90, "nf90_putvar" )
    end if

    ophi_t_count = ophi_t_count+1
    count_time(:) = 1
    start_time(:) = ophi_t_count
    ierr_nf90=nf90_put_var( ncid=ophi_nc, varid=varid_tt, values=(/time/), &
                            start=start_time, count=count_time )

    ny_st = 0
    ny_end = count_ny(1)-1
    count_phi(:) = int((/ 2*nx+1,count_ny,count_nz,1 /), kind=4)
    start_phi(:) = (/ 1,start_ny,start_nz,ophi_t_count /)
    ierr_nf90=nf90_put_var( ncid=ophi_nc, varid=varid_rephi, &
                            values=dble(phi(:,ny_st:ny_end,:)), &
                            start=start_phi, count=count_phi )
    ierr_nf90=nf90_put_var( ncid=ophi_nc, varid=varid_imphi, &
                            values=aimag(phi(:,ny_st:ny_end,:)), &
                            start=start_phi, count=count_phi )
    call check_nf90err( ierr_nf90, "nf90_putvar" )

  END SUBROUTINE fileio_write_phi


!--------------------------------------
  SUBROUTINE fileio_write_Al ( Al, time )
!--------------------------------------

    complex(kind=DP), intent(in), &
         dimension(-nx:nx,0:ny,-nz:nz-1) :: Al
    real(kind=DP), intent(in) :: time

    integer :: dimids(1:4), ierr_nf90
    integer :: varid_kx, varid_ky, varid_zz, varid_is, &
               varid_tt, varid_reAl, varid_imAl
    integer(kind=4) :: start_time(1), count_time(1), &
                       start_Al(1:4), count_Al(1:4)
    integer :: start_ny(1), count_ny(1)
    integer :: start_nz(1), count_nz(1)
    integer :: is
    integer :: ny_st, ny_end
    integer :: ndims, nvars, ngatts, unlimdimid
    integer :: oAl_t_count

    != Write conditions check
    if ( (ranks /= 0) .OR. (vel_rank /= 0) ) then
      return
    end if

    != Information about an open netCDF dataset
    ierr_nf90=nf90_inquire( oAl_nc, ndims, nvars, ngatts, unlimdimid )
    call check_nf90err( ierr_nf90, "nf90_inquire" )

    if ( ndims == 0 ) then
      != Define dimensions in file
      ierr_nf90=nf90_def_dim( ncid=oAl_nc, name="kx", &
                              len=int(2*nx+1,kind=4),      dimid=dimids(1) )
      ierr_nf90=nf90_def_dim( ncid=oAl_nc, name="ky", &
                              len=int(global_ny+1,kind=4), dimid=dimids(2) )
      ierr_nf90=nf90_def_dim( ncid=oAl_nc, name="zz", &
                              len=int(2*global_nz,kind=4), dimid=dimids(3) )
      ierr_nf90=nf90_def_dim( ncid=oAl_nc, name="t",  &
                              len=NF90_UNLIMITED,          dimid=dimids(4) )
      call check_nf90err( ierr_nf90, "nf90_def_dim" )

      != Define variables in file
      ierr_nf90=nf90_def_var( ncid=oAl_nc, name="kx",   xtype=NF90_DOUBLE, &
                              dimids=dimids(1),   varid=varid_kx )
      ierr_nf90=nf90_def_var( ncid=oAl_nc, name="ky",   xtype=NF90_DOUBLE, &
                              dimids=dimids(2),   varid=varid_ky )
      ierr_nf90=nf90_def_var( ncid=oAl_nc, name="zz",   xtype=NF90_DOUBLE, &
                              dimids=dimids(3),   varid=varid_zz )
      ierr_nf90=nf90_def_var( ncid=oAl_nc, name="t",    xtype=NF90_DOUBLE, &
                              dimids=dimids(4),   varid=varid_tt )
      ierr_nf90=nf90_def_var( ncid=oAl_nc, name="reAl", xtype=NF90_DOUBLE, &
                              dimids=dimids(1:4), varid=varid_reAl )
      ierr_nf90=nf90_def_var( ncid=oAl_nc, name="imAl", xtype=NF90_DOUBLE, &
                              dimids=dimids(1:4), varid=varid_imAl )
      call check_nf90err( ierr_nf90, "nf90_def_var" )

      != End of definition of file
      ierr_nf90=nf90_enddef( ncid=oAl_nc )
      call check_nf90err( ierr_nf90, "nf90_enddef" )

      oAl_t_count = 0
    else
      ierr_nf90=nf90_inq_varid( oAl_nc, "t",    varid_tt )
      ierr_nf90=nf90_inq_varid( oAl_nc, "reAl", varid_reAl )
      ierr_nf90=nf90_inq_varid( oAl_nc, "imAl", varid_imAl )
      call check_nf90err( ierr_nf90, "nf90_inq_varid" )

      ierr_nf90=nf90_inquire_dimension(oAl_nc, unlimdimid, len=oAl_t_count)
      call check_nf90err(ierr_nf90, "nf90_inquire_dimension")
    end if

    != Parallel data access type
    ierr_nf90=nf90_var_par_access( ncid=oAl_nc, varid=varid_kx, &
                                   access=NF90_INDEPENDENT )
    ierr_nf90=nf90_var_par_access( ncid=oAl_nc, varid=varid_ky, &
                                   access=NF90_INDEPENDENT )
    ierr_nf90=nf90_var_par_access( ncid=oAl_nc, varid=varid_zz, &
                                   access=NF90_INDEPENDENT )
    ierr_nf90=nf90_var_par_access( ncid=oAl_nc, varid=varid_tt, &
                                   access=NF90_COLLECTIVE )
    ierr_nf90=nf90_var_par_access( ncid=oAl_nc, varid=varid_reAl, &
                                   access=NF90_COLLECTIVE )
    ierr_nf90=nf90_var_par_access( ncid=oAl_nc, varid=varid_imAl, &
                                   access=NF90_COLLECTIVE )
    call check_nf90err( ierr_nf90, "nf90_var_par_access" )

    count_ny=nsize_y
    start_ny=ist_y_g+1
    count_nz=nz*2
    start_nz=count_nz*rankz+1

    != Write variables: static coordinates x and y
    if ( ndims == 0 ) then
      ierr_nf90=nf90_put_var( ncid=oAl_nc, varid=varid_kx, &
                              values=kx(-nx:nx) )
      ierr_nf90=nf90_put_var( ncid=oAl_nc, varid=varid_ky, &
                              values=ky(0:iend_y), &
                              start=start_ny, count=count_ny )
      ierr_nf90=nf90_put_var( ncid=oAl_nc, varid=varid_zz, &
                              values=zz(-nz:nz-1), &
                              start=start_nz, count=count_nz )
      call check_nf90err( ierr_nf90, "nf90_putvar" )
    end if

    oAl_t_count = oAl_t_count+1
    count_time(:) = 1
    start_time(:) = oAl_t_count
    ierr_nf90=nf90_put_var( ncid=oAl_nc, varid=varid_tt, values=(/time/), &
                            start=start_time, count=count_time )

    ny_st = 0
    ny_end = count_ny(1)-1
    count_Al(:) = int((/ 2*nx+1,count_ny,count_nz,1 /), kind=4)
    start_Al(:) = (/ 1,start_ny,start_nz,oAl_t_count /)
    ierr_nf90=nf90_put_var( ncid=oAl_nc, varid=varid_reAl, &
                            values=dble(Al(:,ny_st:ny_end,:)), &
                            start=start_Al, count=count_Al )
    ierr_nf90=nf90_put_var( ncid=oAl_nc, varid=varid_imAl, &
                            values=aimag(Al(:,ny_st:ny_end,:)), &
                            start=start_Al, count=count_Al )
    call check_nf90err( ierr_nf90, "nf90_putvar" )

  END SUBROUTINE fileio_write_Al


!--------------------------------------
  SUBROUTINE fileio_write_mom ( dens, upara, ppara, pperp, qlpara, qlperp, time )
!--------------------------------------

    complex(kind=DP), intent(in), &
         dimension(-nx:nx,0:ny,-nz:nz-1) ::  dens, upara, ppara, pperp, qlpara, qlperp
    real(kind=DP), intent(in) :: time

    integer, parameter :: nmom = 6   ! Number of output moments

    integer :: dimids(1:6), ierr_nf90
    integer :: varid_kx, varid_ky, varid_zz, varid_imom, varid_is, &
               varid_tt, varid_remom, varid_immom
    integer(kind=4) :: start_time(1),  count_time(1), &
                       start_mom(1:6), count_mom(1:6)
    integer :: start_ny(1), count_ny(1)
    integer :: start_nz(1), count_nz(1)
    integer :: start_nm(1), count_nm(1)
    integer :: imom, is
    integer :: ny_st, ny_end
    integer :: ndims, nvars, ngatts, unlimdimid
    integer :: omom_t_count

    != Write conditions check
    if ( vel_rank /= 0 ) return

    != Information about an open netCDF dataset
    ierr_nf90=nf90_inquire( omom_nc, ndims, nvars, ngatts, unlimdimid )
    call check_nf90err( ierr_nf90, "nf90_inquire" )

    if ( ndims == 0 ) then
      != Define dimensions in file
      ierr_nf90=nf90_def_dim( ncid=omom_nc, name="kx", &
                              len=int(2*nx+1,kind=4),      dimid=dimids(1) )
      ierr_nf90=nf90_def_dim( ncid=omom_nc, name="ky", &
                              len=int(global_ny+1,kind=4), dimid=dimids(2) )
      ierr_nf90=nf90_def_dim( ncid=omom_nc, name="zz", &
                              len=int(2*global_nz,kind=4), dimid=dimids(3) )
      ierr_nf90=nf90_def_dim( ncid=omom_nc, name="imom", &
                              len=int(nmom,kind=4),        dimid=dimids(4) )
      ierr_nf90=nf90_def_dim( ncid=omom_nc, name="is", &
                              len=int(ns,kind=4),          dimid=dimids(5) )
      ierr_nf90=nf90_def_dim( ncid=omom_nc, name="t",  &
                              len=NF90_UNLIMITED,          dimid=dimids(6) )
      call check_nf90err( ierr_nf90, "nf90_def_dim" )

      != Define variables in file
      ierr_nf90=nf90_def_var( ncid=omom_nc, name="kx",    xtype=NF90_DOUBLE, &
                              dimids=dimids(1),   varid=varid_kx )
      ierr_nf90=nf90_def_var( ncid=omom_nc, name="ky",    xtype=NF90_DOUBLE, &
                              dimids=dimids(2),   varid=varid_ky )
      ierr_nf90=nf90_def_var( ncid=omom_nc, name="zz",    xtype=NF90_DOUBLE, &
                              dimids=dimids(3),   varid=varid_zz )
      ierr_nf90=nf90_def_var( ncid=omom_nc, name="imom",  xtype=NF90_INT, &
                              dimids=dimids(4),   varid=varid_imom )
      ierr_nf90=nf90_def_var( ncid=omom_nc, name="is",    xtype=NF90_INT, &
                              dimids=dimids(5),   varid=varid_is )
      ierr_nf90=nf90_def_var( ncid=omom_nc, name="t",     xtype=NF90_DOUBLE, &
                              dimids=dimids(6),   varid=varid_tt )
      ierr_nf90=nf90_def_var( ncid=omom_nc, name="remom", xtype=NF90_DOUBLE, &
                              dimids=dimids(1:6), varid=varid_remom )
      ierr_nf90=nf90_def_var( ncid=omom_nc, name="immom", xtype=NF90_DOUBLE, &
                              dimids=dimids(1:6), varid=varid_immom )
      call check_nf90err( ierr_nf90, "nf90_def_var" )

      != End of definition of file
      ierr_nf90=nf90_enddef( ncid=omom_nc )
      call check_nf90err( ierr_nf90, "nf90_enddef" )

      omom_t_count = 0
    else
      ierr_nf90=nf90_inq_varid( omom_nc, "t",     varid_tt )
      ierr_nf90=nf90_inq_varid( omom_nc, "remom", varid_remom )
      ierr_nf90=nf90_inq_varid( omom_nc, "immom", varid_immom )
      call check_nf90err( ierr_nf90, "nf90_inq_varid" )

      ierr_nf90=nf90_inquire_dimension(omom_nc, unlimdimid, len=omom_t_count)
      call check_nf90err(ierr_nf90, "nf90_inquire_dimension")
    end if

    != Parallel data access type
    ierr_nf90=nf90_var_par_access( ncid=omom_nc, varid=varid_kx, &
                                   access=NF90_INDEPENDENT )
    ierr_nf90=nf90_var_par_access( ncid=omom_nc, varid=varid_ky, &
                                   access=NF90_INDEPENDENT )
    ierr_nf90=nf90_var_par_access( ncid=omom_nc, varid=varid_zz, &
                                   access=NF90_INDEPENDENT )
    ierr_nf90=nf90_var_par_access( ncid=omom_nc, varid=varid_imom, &
                                   access=NF90_INDEPENDENT )
    ierr_nf90=nf90_var_par_access( ncid=omom_nc, varid=varid_is, &
                                   access=NF90_INDEPENDENT )
    ierr_nf90=nf90_var_par_access( ncid=omom_nc, varid=varid_tt, &
                                   access=NF90_COLLECTIVE )
    ierr_nf90=nf90_var_par_access( ncid=omom_nc, varid=varid_remom, &
                                   access=NF90_COLLECTIVE )
    ierr_nf90=nf90_var_par_access( ncid=omom_nc, varid=varid_immom, &
                                   access=NF90_COLLECTIVE )
    call check_nf90err( ierr_nf90, "nf90_var_par_access" )

    count_ny=nsize_y
    start_ny=ist_y_g+1
    count_nz=nz*2
    start_nz=count_nz*rankz+1
    count_nm=1
    start_nm=1

    != Write variables: static coordinates x and y
    if ( ndims == 0 ) then
      ierr_nf90=nf90_put_var( ncid=omom_nc, varid=varid_kx, &
                              values=kx(-nx:nx) )
      ierr_nf90=nf90_put_var( ncid=omom_nc, varid=varid_ky, &
                              values=ky(0:iend_y), &
                              start=start_ny, count=count_ny )
      ierr_nf90=nf90_put_var( ncid=omom_nc, varid=varid_zz, &
                              values=zz(-nz:nz-1), &
                              start=start_nz, count=count_nz )
      ierr_nf90=nf90_put_var( ncid=omom_nc, varid=varid_imom, &
                              values=(/ (imom, imom=0,nmom-1) /) )
      ierr_nf90=nf90_put_var( ncid=omom_nc, varid=varid_is, &
                              values=(/ (is, is=0,ns-1) /) )
      call check_nf90err( ierr_nf90, "nf90_putvar" )
    end if

    omom_t_count = omom_t_count+1
    count_time(:) = 1
    start_time(:) = omom_t_count
    ierr_nf90=nf90_put_var( ncid=omom_nc, varid=varid_tt, values=(/time/), &
                            start=start_time, count=count_time )

    ny_st = 0
    ny_end = count_ny(1)-1
    count_mom(:) = int((/ 2*nx+1,count_ny,count_nz,count_nm,1,1 /), kind=4)
    start_mom(:) = (/ 1,start_ny,start_nz,start_nm,1+ranks,omom_t_count /)
    ierr_nf90=nf90_put_var( ncid=omom_nc, varid=varid_remom, &
                            values=dble(dens(:,ny_st:ny_end,:)), &
                            start=start_mom, count=count_mom )
    ierr_nf90=nf90_put_var( ncid=omom_nc, varid=varid_immom, &
                            values=aimag(dens(:,ny_st:ny_end,:)), &
                            start=start_mom, count=count_mom )
    start_nm=2
    start_mom(:) = (/ 1,start_ny,start_nz,start_nm,1+ranks,omom_t_count /)
    ierr_nf90=nf90_put_var( ncid=omom_nc, varid=varid_remom, &
                            values=dble(upara(:,ny_st:ny_end,:)), &
                            start=start_mom, count=count_mom )
    ierr_nf90=nf90_put_var( ncid=omom_nc, varid=varid_immom, &
                            values=aimag(upara(:,ny_st:ny_end,:)), &
                            start=start_mom, count=count_mom )
    start_nm=3
    start_mom(:) = (/ 1,start_ny,start_nz,start_nm,1+ranks,omom_t_count /)
    ierr_nf90=nf90_put_var( ncid=omom_nc, varid=varid_remom, &
                            values=dble(ppara(:,ny_st:ny_end,:)), &
                            start=start_mom, count=count_mom )
    ierr_nf90=nf90_put_var( ncid=omom_nc, varid=varid_immom, &
                            values=aimag(ppara(:,ny_st:ny_end,:)), &
                            start=start_mom, count=count_mom )
    start_nm=4
    start_mom(:) = (/ 1,start_ny,start_nz,start_nm,1+ranks,omom_t_count /)
    ierr_nf90=nf90_put_var( ncid=omom_nc, varid=varid_remom, &
                            values=dble(pperp(:,ny_st:ny_end,:)), &
                            start=start_mom, count=count_mom )
    ierr_nf90=nf90_put_var( ncid=omom_nc, varid=varid_immom, &
                            values=aimag(pperp(:,ny_st:ny_end,:)), &
                            start=start_mom, count=count_mom )
    start_nm=5
    start_mom(:) = (/ 1,start_ny,start_nz,start_nm,1+ranks,omom_t_count /)
    ierr_nf90=nf90_put_var( ncid=omom_nc, varid=varid_remom, &
                            values=dble(qlpara(:,ny_st:ny_end,:)), &
                            start=start_mom, count=count_mom )
    ierr_nf90=nf90_put_var( ncid=omom_nc, varid=varid_immom, &
                            values=aimag(qlpara(:,ny_st:ny_end,:)), &
                            start=start_mom, count=count_mom )
    start_nm=6
    start_mom(:) = (/ 1,start_ny,start_nz,start_nm,1+ranks,omom_t_count /)
    ierr_nf90=nf90_put_var( ncid=omom_nc, varid=varid_remom, &
                            values=dble(qlperp(:,ny_st:ny_end,:)), &
                            start=start_mom, count=count_mom )
    ierr_nf90=nf90_put_var( ncid=omom_nc, varid=varid_immom, &
                            values=aimag(qlperp(:,ny_st:ny_end,:)), &
                            start=start_mom, count=count_mom )
    call check_nf90err( ierr_nf90, "nf90_putvar" )

  END SUBROUTINE fileio_write_mom


!--------------------------------------
  SUBROUTINE fileio_write_trn ( entrpy, fenegy, menegy, peint, pmint, &
                                neint, nmint, dcd, pflux_es, pflux_em, &
                                eflux_es, eflux_em, time )
!--------------------------------------

    real(kind=DP), intent(in), &
         dimension(-nx:nx,0:ny) :: entrpy, fenegy, menegy, peint, &
                                   pmint, neint, nmint, dcd, &
                                   pflux_es, pflux_em, eflux_es, eflux_em
    real(kind=DP), intent(in) :: time

    integer, parameter :: ntrn = 12  ! Number of output total transfer diagnostics

    integer :: dimids(1:5), ierr_nf90
    integer :: varid_kx, varid_ky, varid_itrn, varid_is, varid_tt, varid_trn
    integer(kind=4) :: start_time(1),  count_time(1), &
                       start_trn(1:5), count_trn(1:5)
    integer :: start_ny(1), count_ny(1)
    integer :: start_nt(1), count_nt(1)
    integer :: itrn, is
    integer :: ny_st, ny_end
    integer :: ndims, nvars, ngatts, unlimdimid
    integer :: otrn_t_count

    != Write conditions check
    if ( (zsp_rank /= 0) .OR. (vel_rank /= 0) ) then
      return
    end if

    != Information about an open netCDF dataset
    ierr_nf90=nf90_inquire( otrn_nc, ndims, nvars, ngatts, unlimdimid )
    call check_nf90err( ierr_nf90, "nf90_inquire" )

    if ( ndims == 0 ) then
      != Define dimensions in file
      ierr_nf90=nf90_def_dim( ncid=otrn_nc, name="kx", &
                              len=int(2*nx+1,kind=4),      dimid=dimids(1) )
      ierr_nf90=nf90_def_dim( ncid=otrn_nc, name="ky", &
                              len=int(global_ny+1,kind=4), dimid=dimids(2) )
      ierr_nf90=nf90_def_dim( ncid=otrn_nc, name="itrn", &
                              len=int(ntrn,kind=4),        dimid=dimids(3) )
      ierr_nf90=nf90_def_dim( ncid=otrn_nc, name="is", &
                              len=int(ns,kind=4),          dimid=dimids(4) )
      ierr_nf90=nf90_def_dim( ncid=otrn_nc, name="t",  &
                              len=NF90_UNLIMITED,          dimid=dimids(5) )
      call check_nf90err( ierr_nf90, "nf90_def_dim" )

      != Define variables in file
      ierr_nf90=nf90_def_var( ncid=otrn_nc, name="kx",    xtype=NF90_DOUBLE, &
                              dimids=dimids(1),   varid=varid_kx )
      ierr_nf90=nf90_def_var( ncid=otrn_nc, name="ky",    xtype=NF90_DOUBLE, &
                              dimids=dimids(2),   varid=varid_ky )
      ierr_nf90=nf90_def_var( ncid=otrn_nc, name="itrn",  xtype=NF90_INT, &
                              dimids=dimids(3),   varid=varid_itrn )
      ierr_nf90=nf90_def_var( ncid=otrn_nc, name="is",    xtype=NF90_INT, &
                              dimids=dimids(4),   varid=varid_is )
      ierr_nf90=nf90_def_var( ncid=otrn_nc, name="t",     xtype=NF90_DOUBLE, &
                              dimids=dimids(5),   varid=varid_tt )
      ierr_nf90=nf90_def_var( ncid=otrn_nc, name="trn", xtype=NF90_DOUBLE, &
                              dimids=dimids(1:5), varid=varid_trn )
      call check_nf90err( ierr_nf90, "nf90_def_var" )

      != End of definition of file
      ierr_nf90=nf90_enddef( ncid=otrn_nc )
      call check_nf90err( ierr_nf90, "nf90_enddef" )

      otrn_t_count = 0
    else
      ierr_nf90=nf90_inq_varid( otrn_nc, "t",   varid_tt )
      ierr_nf90=nf90_inq_varid( otrn_nc, "trn", varid_trn )
      call check_nf90err( ierr_nf90, "nf90_inq_varid" )

      ierr_nf90=nf90_inquire_dimension(otrn_nc, unlimdimid, len=otrn_t_count)
      call check_nf90err(ierr_nf90, "nf90_inquire_dimension")
    end if

    != Parallel data access type
    ierr_nf90=nf90_var_par_access( ncid=otrn_nc, varid=varid_kx, &
                                   access=NF90_INDEPENDENT )
    ierr_nf90=nf90_var_par_access( ncid=otrn_nc, varid=varid_ky, &
                                   access=NF90_INDEPENDENT )
    ierr_nf90=nf90_var_par_access( ncid=otrn_nc, varid=varid_itrn, &
                                   access=NF90_INDEPENDENT )
    ierr_nf90=nf90_var_par_access( ncid=otrn_nc, varid=varid_is, &
                                   access=NF90_INDEPENDENT )
    ierr_nf90=nf90_var_par_access( ncid=otrn_nc, varid=varid_tt, &
                                   access=NF90_COLLECTIVE )
    ierr_nf90=nf90_var_par_access( ncid=otrn_nc, varid=varid_trn, &
                                   access=NF90_COLLECTIVE )
    call check_nf90err( ierr_nf90, "nf90_var_par_access" )

    count_ny=nsize_y
    start_ny=ist_y_g+1
    count_nt=1
    start_nt=1

    != Write variables: static coordinates x and y
    if ( ndims == 0 ) then
      ierr_nf90=nf90_put_var( ncid=otrn_nc, varid=varid_kx, &
                              values=kx(-nx:nx) )
      ierr_nf90=nf90_put_var( ncid=otrn_nc, varid=varid_ky, &
                              values=ky(0:iend_y), &
                              start=start_ny, count=count_ny )
      ierr_nf90=nf90_put_var( ncid=otrn_nc, varid=varid_itrn, &
                              values=(/ (itrn, itrn=0,ntrn-1) /) )
      ierr_nf90=nf90_put_var( ncid=otrn_nc, varid=varid_is, &
                              values=(/ (is, is=0,ns-1) /) )
      call check_nf90err( ierr_nf90, "nf90_putvar(kx,ky,itrn,is)" )
    end if

    otrn_t_count = otrn_t_count + 1
    count_time(:) = 1
    start_time(:) = otrn_t_count
    ierr_nf90=nf90_put_var( ncid=otrn_nc, varid=varid_tt, values=(/time/), &
                            start=start_time, count=count_time )
    call check_nf90err( ierr_nf90, "nf90_putvar(time)" )

    ny_st = 0
    ny_end = count_ny(1)-1
    count_trn(:) = int((/ 2*nx+1,count_ny,count_nt,1,1 /), kind=4)
    start_trn(:) = (/ 1,start_ny,start_nt,1+ranks,otrn_t_count /)
    ierr_nf90=nf90_put_var( ncid=otrn_nc, varid=varid_trn, &
                            values=entrpy(:,ny_st:ny_end), &
                            start=start_trn, count=count_trn )
    start_nt=2
    start_trn(:) = (/ 1,start_ny,start_nt,1+ranks,otrn_t_count /)
    ierr_nf90=nf90_put_var( ncid=otrn_nc, varid=varid_trn, &
                            values=fenegy(:,ny_st:ny_end), &
                            start=start_trn, count=count_trn )
    start_nt=3
    start_trn(:) = (/ 1,start_ny,start_nt,1+ranks,otrn_t_count /)
    ierr_nf90=nf90_put_var( ncid=otrn_nc, varid=varid_trn, &
                            values=menegy(:,ny_st:ny_end), &
                            start=start_trn, count=count_trn )
    start_nt=4
    start_trn(:) = (/ 1,start_ny,start_nt,1+ranks,otrn_t_count /)
    ierr_nf90=nf90_put_var( ncid=otrn_nc, varid=varid_trn, &
                            values=peint(:,ny_st:ny_end), &
                            start=start_trn, count=count_trn )
    start_nt=5
    start_trn(:) = (/ 1,start_ny,start_nt,1+ranks,otrn_t_count /)
    ierr_nf90=nf90_put_var( ncid=otrn_nc, varid=varid_trn, &
                            values=pmint(:,ny_st:ny_end), &
                            start=start_trn, count=count_trn )
    start_nt=6
    start_trn(:) = (/ 1,start_ny,start_nt,1+ranks,otrn_t_count /)
    ierr_nf90=nf90_put_var( ncid=otrn_nc, varid=varid_trn, &
                            values=neint(:,ny_st:ny_end), &
                            start=start_trn, count=count_trn )
    start_nt=7
    start_trn(:) = (/ 1,start_ny,start_nt,1+ranks,otrn_t_count /)
    ierr_nf90=nf90_put_var( ncid=otrn_nc, varid=varid_trn, &
                            values=nmint(:,ny_st:ny_end), &
                            start=start_trn, count=count_trn )
    start_nt=8
    start_trn(:) = (/ 1,start_ny,start_nt,1+ranks,otrn_t_count /)
    ierr_nf90=nf90_put_var( ncid=otrn_nc, varid=varid_trn, &
                            values=dcd(:,ny_st:ny_end), &
                            start=start_trn, count=count_trn )
    start_nt=9
    start_trn(:) = (/ 1,start_ny,start_nt,1+ranks,otrn_t_count /)
    ierr_nf90=nf90_put_var( ncid=otrn_nc, varid=varid_trn, &
                            values=pflux_es(:,ny_st:ny_end), &
                            start=start_trn, count=count_trn )
    start_nt=10
    start_trn(:) = (/ 1,start_ny,start_nt,1+ranks,otrn_t_count /)
    ierr_nf90=nf90_put_var( ncid=otrn_nc, varid=varid_trn, &
                            values=pflux_em(:,ny_st:ny_end), &
                            start=start_trn, count=count_trn )
    start_nt=11
    start_trn(:) = (/ 1,start_ny,start_nt,1+ranks,otrn_t_count /)
    ierr_nf90=nf90_put_var( ncid=otrn_nc, varid=varid_trn, &
                            values=eflux_es(:,ny_st:ny_end), &
                            start=start_trn, count=count_trn )
    start_nt=12
    start_trn(:) = (/ 1,start_ny,start_nt,1+ranks,otrn_t_count /)
    ierr_nf90=nf90_put_var( ncid=otrn_nc, varid=varid_trn, &
                            values=eflux_em(:,ny_st:ny_end), &
                            start=start_trn, count=count_trn )
    call check_nf90err( ierr_nf90, "nf90_putvar(trn)" )

  END SUBROUTINE fileio_write_trn


!--------------------------------------
  SUBROUTINE fileio_write_tri ( jkpq_es, jpqk_es, jqkp_es, &
                                jkpq_em, jpqk_em, jqkp_em, time )
!--------------------------------------

    real(kind=DP), intent(in), &
         dimension(-nx:nx,-global_ny:global_ny) :: jkpq_es, jpqk_es, jqkp_es, &
         jkpq_em, jpqk_em, jqkp_em
    real(kind=DP), intent(in) :: time

    integer, parameter :: ntri = 6   ! Number of output triad transfer diagnostics

    integer :: dimids(1:5), ierr_nf90
    integer :: varid_kx, varid_ky, varid_itri, varid_is, varid_tt, varid_tri
    integer(kind=4) :: start_time(1),  count_time(1), &
                       start_tri(1:5), count_tri(1:5)
    integer :: start_ny(1), count_ny(1)
    integer :: start_nt(1), count_nt(1)
    integer :: itri, is
    integer :: ny_st, ny_end
    integer :: ndims, nvars, ngatts, unlimdimid
    integer :: otri_t_count

    real(kind=DP), dimension(0:global_ny)            :: gky
    integer :: my

    character(256)   :: env_string

    integer       :: n_tht, m_j
    real(kind=DP) :: kymin, del_c

    namelist /nperi/ n_tht, kymin, m_j, del_c

    if ( rank == 0 ) then
      != Information about an open netCDF dataset
      ierr_nf90=nf90_inquire( otri_nc, ndims, nvars, ngatts, unlimdimid )
      call check_nf90err( ierr_nf90, "nf90_inquire" )

      if ( ndims == 0 ) then
        != Define dimensions in file
        ierr_nf90=nf90_def_dim( ncid=otri_nc, name="kx", &
                                len=int(2*nx+1,kind=4),        dimid=dimids(1) )
        ierr_nf90=nf90_def_dim( ncid=otri_nc, name="ky", &
                                len=int(2*global_ny+1,kind=4), dimid=dimids(2) )
        ierr_nf90=nf90_def_dim( ncid=otri_nc, name="itri", &
                                len=int(ntri,kind=4),          dimid=dimids(3) )
        ierr_nf90=nf90_def_dim( ncid=otri_nc, name="is", &
                                len=int(ns,kind=4),            dimid=dimids(4) )
        ierr_nf90=nf90_def_dim( ncid=otri_nc, name="t",  &
                                len=NF90_UNLIMITED,            dimid=dimids(5) )
        call check_nf90err( ierr_nf90, "nf90_def_dim" )

        != Define variables in file
        ierr_nf90=nf90_def_var( ncid=otri_nc, name="kx",   xtype=NF90_DOUBLE, &
                                dimids=dimids(1),   varid=varid_kx )
        ierr_nf90=nf90_def_var( ncid=otri_nc, name="ky",   xtype=NF90_DOUBLE, &
                                dimids=dimids(2),   varid=varid_ky )
        ierr_nf90=nf90_def_var( ncid=otri_nc, name="itri", xtype=NF90_INT, &
                                dimids=dimids(3),   varid=varid_itri )
        ierr_nf90=nf90_def_var( ncid=otri_nc, name="is",   xtype=NF90_INT, &
                                dimids=dimids(4),   varid=varid_is )
        ierr_nf90=nf90_def_var( ncid=otri_nc, name="t",    xtype=NF90_DOUBLE, &
                                dimids=dimids(5),   varid=varid_tt )
        ierr_nf90=nf90_def_var( ncid=otri_nc, name="tri",  xtype=NF90_DOUBLE, &
                                dimids=dimids(1:5), varid=varid_tri )
        call check_nf90err( ierr_nf90, "nf90_def_var" )

        != End of definition of file
        ierr_nf90=nf90_enddef( ncid=otri_nc )
        call check_nf90err( ierr_nf90, "nf90_enddef" )

        otri_t_count = 0
      else
        ierr_nf90=nf90_inq_varid( otri_nc, "t",   varid_tt )
        ierr_nf90=nf90_inq_varid( otri_nc, "tri", varid_tri )
        call check_nf90err( ierr_nf90, "nf90_inq_varid" )

        ierr_nf90=nf90_inquire_dimension(otri_nc, unlimdimid, len=otri_t_count)
        call check_nf90err(ierr_nf90, "nf90_inquire_dimension")
      end if

      != Parallel data access type
      ierr_nf90=nf90_var_par_access( ncid=otri_nc, varid=varid_kx, &
                                     access=NF90_INDEPENDENT )
      ierr_nf90=nf90_var_par_access( ncid=otri_nc, varid=varid_ky, &
                                     access=NF90_INDEPENDENT )
      ierr_nf90=nf90_var_par_access( ncid=otri_nc, varid=varid_itri, &
                                     access=NF90_INDEPENDENT )
      ierr_nf90=nf90_var_par_access( ncid=otri_nc, varid=varid_is, &
                                     access=NF90_INDEPENDENT )
      ierr_nf90=nf90_var_par_access( ncid=otri_nc, varid=varid_tt, &
                                     access=NF90_COLLECTIVE )
      ierr_nf90=nf90_var_par_access( ncid=otri_nc, varid=varid_tri, &
                                     access=NF90_COLLECTIVE )
      call check_nf90err( ierr_nf90, "nf90_var_par_access" )

      count_ny=2*global_ny+1
      start_ny=1
      count_nt=1
      start_nt=1

      if ( ndims == 0 ) then
        close(inml)
        call getenv ( 'fu05',env_string )
        open(inml, file=env_string )

        read(inml, nml=nperi)

        do my = 0, global_ny
          gky(my) = kymin * real( my, kind=DP )
        end do
      end if

      != Write variables: static coordinates x and y
      if ( ndims == 0 ) then
        ierr_nf90=nf90_put_var( ncid=otri_nc, varid=varid_kx, &
                                values=kx(-nx:nx) )
        ierr_nf90=nf90_put_var( ncid=otri_nc, varid=varid_ky, &
                                values=(/ -gky(global_ny:1:-1),gky(0:global_ny) /) )
        ierr_nf90=nf90_put_var( ncid=otri_nc, varid=varid_itri, &
                                values=(/ (itri, itri=0,ntri-1) /) )
        ierr_nf90=nf90_put_var( ncid=otri_nc, varid=varid_is, &
                                values=(/ (is, is=0,ns-1) /) )
        call check_nf90err( ierr_nf90, "nf90_putvar" )
      end if

      otri_t_count = otri_t_count+1
      count_time(:) = 1
      start_time(:) = otri_t_count
      ierr_nf90=nf90_put_var( ncid=otri_nc, varid=varid_tt, values=(/time/), &
                              start=start_time, count=count_time )

      ny_st = -global_ny
      ny_end = global_ny
      count_tri(:) = int((/ 2*nx+1,count_ny,count_nt,1,1 /), kind=4)
      start_tri(:) = (/ 1,start_ny,start_nt,1+ranks,otri_t_count /)
      ierr_nf90=nf90_put_var( ncid=otri_nc, varid=varid_tri, &
                              values=jkpq_es(:,ny_st:ny_end), &
                              start=start_tri, count=count_tri )
      start_nt=2
      start_tri(:) = (/ 1,start_ny,start_nt,1+ranks,otri_t_count /)
      ierr_nf90=nf90_put_var( ncid=otri_nc, varid=varid_tri, &
                              values=jpqk_es(:,ny_st:ny_end), &
                              start=start_tri, count=count_tri )
      start_nt=3
      start_tri(:) = (/ 1,start_ny,start_nt,1+ranks,otri_t_count /)
      ierr_nf90=nf90_put_var( ncid=otri_nc, varid=varid_tri, &
                              values=jqkp_es(:,ny_st:ny_end), &
                              start=start_tri, count=count_tri )
      start_nt=4
      start_tri(:) = (/ 1,start_ny,start_nt,1+ranks,otri_t_count /)
      ierr_nf90=nf90_put_var( ncid=otri_nc, varid=varid_tri, &
                              values=jkpq_em(:,ny_st:ny_end), &
                              start=start_tri, count=count_tri )
      start_nt=5
      start_tri(:) = (/ 1,start_ny,start_nt,1+ranks,otri_t_count /)
      ierr_nf90=nf90_put_var( ncid=otri_nc, varid=varid_tri, &
                              values=jpqk_em(:,ny_st:ny_end), &
                              start=start_tri, count=count_tri )
      start_nt=6
      start_tri(:) = (/ 1,start_ny,start_nt,1+ranks,otri_t_count /)
      ierr_nf90=nf90_put_var( ncid=otri_nc, varid=varid_tri, &
                              values=jqkp_em(:,ny_st:ny_end), &
                              start=start_tri, count=count_tri )

      call check_nf90err( ierr_nf90, "nf90_putvar" )
    end if

  END SUBROUTINE fileio_write_tri


  !--------------------------------------
  SUBROUTINE check_nf90err(werr, comment)
  !--------------------------------------
  !  Check error message of nf90
    integer(kind=4), intent(in) :: werr
    character(len=*), intent(in) :: comment

    if(werr /= nf90_noerr) then
      write(*,*) comment//" "//trim(nf90_strerror(werr))
      stop
    end if

  END SUBROUTINE check_nf90err


END MODULE GKV_fileio
