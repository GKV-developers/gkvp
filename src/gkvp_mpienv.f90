MODULE GKV_mpienv
!-------------------------------------------------------------------------------
!
!    Header and settings for using MPI
!
!    Update history of gkvp_mpienv.f90
!    --------------
!      gkvp_f0.57 (S. Maeyama, Oct 2020)
!        - Version number f0.57 is removed from filename.
!
!-------------------------------------------------------------------------------

!--------------------------------------
!  Variables for MPI parallilization
!--------------------------------------
  use GKV_header
  implicit none

  public

  include "mpif.h"

  integer :: rankg, nproc
  integer :: sizedouble_c, ierr_mpi
  integer, dimension(MPI_STATUS_SIZE)  :: status

!fj> add
  integer :: rankw
!fj<
  integer :: rankz, rankv, rankm, ranks
  integer :: izup, izdn, ivup, ivdn, imup, imdn

  integer ::  sendzdn, sendzup, recvzdn, recvzup
  integer ::  sendvdn, sendvup, recvvdn, recvvup
  integer ::  sendmdn, sendmup, recvmdn, recvmup

  integer ::  vel_comm_world, zsp_comm_world, vcolor, zcolor
  integer ::  spc_comm_world, sub_comm_world, scolor
  integer ::  col_comm_world, ccolor
  integer ::  vel_rank, vel_nproc
  integer ::  zsp_rank, zsp_nproc
  integer ::  spc_rank, spc_nproc
  integer ::  col_rank, col_nproc
  integer ::      rank, sub_nproc
!fj> add
  integer ::  fft_comm_world, wcolor
  integer ::  fft_rank, fft_nproc
!fj<
  
  integer :: ornk 

CONTAINS

!--------------------------------------
!fj>
!cc  SUBROUTINE mpienv_init( nprocz, nprocv, nprocm )
  SUBROUTINE mpienv_init( nprocw, nprocz, nprocv, nprocm, nprocs )
!fj<
!--------------------------------------

!fj>
!cc  integer, intent(in) :: nprocz, nprocv, nprocm
  integer, intent(in) :: nprocw, nprocz, nprocv, nprocm, nprocs
!fj<

  integer :: ny_size, nxw_sz, nwk


!--- begin MPI settings

    call MPI_Init ( ierr_mpi )

    call MPI_Comm_rank ( MPI_COMM_WORLD, rankg, ierr_mpi )

    call MPI_Comm_size ( MPI_COMM_WORLD, nproc, ierr_mpi )

    call MPI_Type_size ( MPI_DOUBLE_COMPLEX, sizedouble_c, ierr_mpi )

! --- for multispecies 
    ranks  = rankg / ( nprocz * nprocv * nprocm * nprocw ) 

    call MPI_Comm_split ( MPI_COMM_WORLD, ranks, rankg, sub_comm_world, ierr_mpi )
    call MPI_Comm_rank ( sub_comm_world, rank, ierr_mpi )
    call MPI_Comm_size ( sub_comm_world, sub_nproc, ierr_mpi )

! --- additional comunicator for integration over species
    scolor  = mod(rankg, nprocz*nprocw)
    call MPI_Comm_split ( MPI_COMM_WORLD, scolor, rankg, spc_comm_world, ierr_mpi )
    call MPI_Comm_rank ( spc_comm_world, spc_rank, ierr_mpi )
    call MPI_Comm_size ( spc_comm_world, spc_nproc, ierr_mpi )


! --- process allocation to domain
!fj> mod
!cc    rankz = mod( rank,          nprocz )
!cc    rankv = mod( rank / nprocz, nprocv )
!cc    rankm = rank / ( nprocz * nprocv )
    rankw = mod( rank,                      nprocw )
    rankz = mod( rank / nprocw,             nprocz )
    rankv = mod( rank / ( nprocw * nprocz), nprocv )
    rankm = rank / ( nprocw * nprocz * nprocv )


                             ! rank of the targets
!cc    izup  = rank + 1
!cc    izdn  = rank - 1
!cc    ivup  = rank + nprocz
!cc    ivdn  = rank - nprocz
!cc    imup  = rank + nprocz * nprocv
!cc    imdn  = rank - nprocz * nprocv
    izup  = rank + nprocw
    izdn  = rank - nprocw
    ivup  = rank + nprocw * nprocz
    ivdn  = rank - nprocw * nprocz
    imup  = rank + nprocw * nprocz * nprocv
    imdn  = rank - nprocw * nprocz * nprocv

!cc    if ( rankz == nprocz-1 ) izup  = rank - nprocz + 1
!cc    if ( rankz == 0        ) izdn  = rank + nprocz - 1
    if ( rankz == nprocz-1 ) izup  = rank - nprocw * (nprocz - 1)
    if ( rankz == 0        ) izdn  = rank + nprocw * (nprocz - 1)
!fj<
    if ( rankv == nprocv-1 ) ivup  = MPI_PROC_NULL
    if ( rankv == 0        ) ivdn  = MPI_PROC_NULL
    if ( rankm == nprocm-1 ) imup  = MPI_PROC_NULL
    if ( rankm == 0        ) imdn  = MPI_PROC_NULL

! --- generating sub-worlds

!fj> mod
!cc    vcolor   = rankz
!cc    zcolor   = rank / nprocz
!cc
!cc    call MPI_Comm_split ( sub_comm_world, vcolor, rank, vel_comm_world, ierr_mpi )
!cc    call MPI_Comm_rank ( vel_comm_world, vel_rank, ierr_mpi )
!cc    call MPI_Comm_size ( vel_comm_world, vel_nproc, ierr_mpi )
!cc
!cc    call MPI_Comm_split ( sub_comm_world, zcolor, rank, zsp_comm_world, ierr_mpi )
!cc    call MPI_Comm_rank ( zsp_comm_world, zsp_rank, ierr_mpi )
!cc    call MPI_Comm_size ( zsp_comm_world, zsp_nproc, ierr_mpi )
!cc
!cc
!cc      if ( nproc /= nprocz * nprocv * nprocm ) then

    wcolor   = rank / nprocw
    zcolor   = ( rank / (nprocw*nprocz) ) * nprocw + rankw
    vcolor   = rankz*nprocw  +  rankw

    call MPI_Comm_split ( sub_comm_world, wcolor, rankw, fft_comm_world, ierr_mpi )
    call MPI_Comm_rank ( fft_comm_world, fft_rank, ierr_mpi )
    call MPI_Comm_size ( fft_comm_world, fft_nproc, ierr_mpi )

    call MPI_Comm_split ( sub_comm_world, zcolor, rankz, zsp_comm_world, ierr_mpi )
    call MPI_Comm_rank ( zsp_comm_world, zsp_rank, ierr_mpi )
    call MPI_Comm_size ( zsp_comm_world, zsp_nproc, ierr_mpi )

    call MPI_Comm_split ( sub_comm_world, vcolor, (rankv+nprocv*rankm), vel_comm_world, ierr_mpi )
    call MPI_Comm_rank ( vel_comm_world, vel_rank, ierr_mpi )
    call MPI_Comm_size ( vel_comm_world, vel_nproc, ierr_mpi )

! --- additional communicator for inter-species comm. for field particle part collision
    if ( rank == scolor ) then
      ccolor = mod(rankg, sub_nproc)
    else
      ccolor = 99   ! dummy color
    end if

    call MPI_Comm_split ( MPI_COMM_WORLD, ccolor, rankg, col_comm_world, ierr_mpi )
    call MPI_Comm_rank ( col_comm_world, col_rank, ierr_mpi )
    call MPI_Comm_size ( col_comm_world, col_nproc, ierr_mpi )

    if ( vel_rank > 0 ) then 
      col_rank = MPI_PROC_NULL
    end if
! ---


!cc    write( 6, * ) '##### Debug rank=',rank, &
!cc           ' wcolor,rankw=',wcolor,rankw, ' zcolor, rankz=',zcolor, rankz, &
!cc           ' vcolor,(rankv+nprocv*rankm)=',vcolor,(rankv+nprocv*rankm)

      if ( nproc /= nprocw * nprocz * nprocv * nprocm * nprocs ) then
!fj<
        write( 6, * ) &
           " # proccesor assigment is invalid, nproc = ", nproc
        call MPI_Finalize ( ierr_mpi )
        stop
      end if

!! --- for debug
!      ornk = 1000+rankg
!      write( ornk, "(A25)" ) "# rank config., nproc = ", nproc
!      write( ornk, "(A3)" ) "#  "
!      write( ornk, "(A3)" ) "#  "
!      write( ornk, "(A9,I3)" ) "rankg", rankg
!      write( ornk, "(A9,I3)" ) "ranks", ranks
!      write( ornk, "(A9,I3)" ) "rank", rank
!      write( ornk, "(A9,I3)" ) "rankw", rankm
!      write( ornk, "(A9,I3)" ) "rankz", rankz
!      write( ornk, "(A9,I3)" ) "rankv", rankv
!      write( ornk, "(A9,I3)" ) "rankm", rankm
!      write( ornk, "(A9,I3)" ) "vel_rank", vel_rank
!      write( ornk, "(A9,I3)" ) "zsp_rank", zsp_rank
!      write( ornk, "(A9,I3)" ) "spc_rank", spc_rank
!      write( ornk, "(A9,I3)" ) "fft_rank", fft_rank
!      write( ornk, "(A9,I3)" ) "col_rank", col_rank
!      write( ornk, "(A9,I3)" ) "zcolor", zcolor
!      write( ornk, "(A9,I3)" ) "vcolor", vcolor
!      write( ornk, "(A9,I3)" ) "scolor", scolor
!      write( ornk, "(A9,I3)" ) "ccolor", ccolor


! ---- set y range --------------------------
    ny_size = global_ny + 1 
    if( mod(ny_size,nprocw) == 0 )  then
      nwk    = ny_size / nprocw
    else
      nwk    = ny_size / nprocw + 1
    endif
    !--- global index range ---------------- 
    ist_y_g  = nwk*rankw
    iend_y_g = min( nwk*(rankw+1)-1, (ny_size-1) )
    nsize_y  = iend_y_g - ist_y_g + 1
    !--- local index range ---------------- 
    ist_y    = 0
    iend_y   = iend_y_g - ist_y_g

    if( rankw == 0 )   then
       ist1_y    = 1
    else 
       ist1_y    = 0
    endif

! ---- set xw range ---------------------
    nxw_sz = 2*nxw
    if( mod(nxw_sz,nprocw) == 0 )  then
      nwk    = nxw_sz / nprocw
    else
      nwk    = nxw_sz / nprocw + 1
    endif
    !--- global index range ----------------
    ist_xw_g  = nwk*rankw
    iend_xw_g = min( nwk*(rankw+1)-1, (nxw_sz-1) )
    nsize_xw  = iend_xw_g - ist_xw_g + 1
    !--- local index range ----------------
    ist_xw    = 0
    iend_xw   = iend_xw_g - ist_xw_g


  END SUBROUTINE mpienv_init


END MODULE GKV_mpienv
