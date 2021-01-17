MODULE GKV_fileio
!-------------------------------------------------------------------------------
!
!    File I/O interface for Fortran binary output
!
!    Update history of gkvp_fileio_fortran.f90
!    --------------
!      gkvp_f0.60 (S. Maeyama, Jan 2021)
!        - Fortran binary I/O interface by Fujitsu.
!
!-------------------------------------------------------------------------------

  use GKV_header
  use GKV_mpienv

  implicit  none

  private

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

    character(6)   :: crank
    character(3)   :: cold

    write( crank, fmt="(i6.6)" ) rankg
    write( cold,  fmt="(i3.3)" ) inum-1

    open( icnt, file=path//crank//".cnt."//cold, &
          form="unformatted", status="old", action="read" )

  END SUBROUTINE fileio_open_icnt

!--------------------------------------
  SUBROUTINE fileio_close_icnt
!--------------------------------------

     close( icnt )

  END SUBROUTINE fileio_close_icnt

 
!--------------------------------------
  SUBROUTINE fileio_open_cnt ( path )
!--------------------------------------

    character(*), intent(in) :: path

    character(6)   :: crank
    character(3)   :: cnew

    write( crank, fmt="(i6.6)" ) rankg
    write( cnew,  fmt="(i3.3)" ) inum

    open( ocnt, file=path//crank//".cnt."//cnew, &
          form="unformatted" )

  END SUBROUTINE fileio_open_cnt

!--------------------------------------
  SUBROUTINE fileio_close_cnt
!--------------------------------------

     close( ocnt )

  END SUBROUTINE fileio_close_cnt


!--------------------------------------
  SUBROUTINE fileio_open_fxv ( path )
!--------------------------------------

    character(*), intent(in) :: path

    character(6)   :: crank
    character(1)   :: srank
    character(3)   :: cnew

    write( crank, fmt="(i6.6)" ) rankg
    write( srank, fmt="(i1.1)" ) ranks
    write( cnew,  fmt="(i3.3)" ) inum

    open( ofxv, file=path//crank//"."//srank//".fxv."//cnew, &
          form="unformatted" )

  END SUBROUTINE fileio_open_fxv

!--------------------------------------
  SUBROUTINE fileio_close_fxv
!--------------------------------------

     close( ofxv )

  END SUBROUTINE fileio_close_fxv


!--------------------------------------
  SUBROUTINE fileio_open_phi ( path )
!--------------------------------------

    character(*), intent(in) :: path

    character(6)   :: crank
    character(1)   :: srank
    character(3)   :: cnew

    if ( (ranks /= 0) .OR. (vel_rank /= 0) ) return

    write( crank, fmt="(i6.6)" ) rankg
    write( srank, fmt="(i1.1)" ) ranks
    write( cnew,  fmt="(i3.3)" ) inum

    open( ophi, file=path//crank//"."//srank//".phi."//cnew, &
          form="unformatted" )

  END SUBROUTINE fileio_open_phi

!--------------------------------------
  SUBROUTINE fileio_close_phi
!--------------------------------------

    if ( (ranks /= 0) .OR. (vel_rank /= 0) ) return

    close( ophi )

  END SUBROUTINE fileio_close_phi


!--------------------------------------
  SUBROUTINE fileio_open_Al ( path )
!--------------------------------------

    character(*), intent(in) :: path

    character(6)   :: crank
    character(1)   :: srank
    character(3)   :: cnew

    if ( (ranks /= 0) .OR. (vel_rank /= 0) ) return

    write( crank, fmt="(i6.6)" ) rankg
    write( srank, fmt="(i1.1)" ) ranks
    write( cnew,  fmt="(i3.3)" ) inum

    open( oAl, file=path//crank//"."//srank//".Al."//cnew, &
          form="unformatted" )

  END SUBROUTINE fileio_open_Al

!--------------------------------------
  SUBROUTINE fileio_close_Al
!--------------------------------------

    if ( (ranks /= 0) .OR. (vel_rank /= 0) ) return

    close( oAl )

  END SUBROUTINE fileio_close_Al


!--------------------------------------
  SUBROUTINE fileio_open_mom ( path )
!--------------------------------------

    character(*), intent(in) :: path

    character(6)   :: crank
    character(1)   :: srank
    character(3)   :: cnew

    if ( vel_rank /= 0 ) return

    write( crank, fmt="(i6.6)" ) rankg
    write( srank, fmt="(i1.1)" ) ranks
    write( cnew,  fmt="(i3.3)" ) inum

    open( omom, file=path//crank//"."//srank//".mom."//cnew, &
          form="unformatted" )

  END SUBROUTINE fileio_open_mom

!--------------------------------------
  SUBROUTINE fileio_close_mom
!--------------------------------------

    if ( vel_rank /= 0 ) return

    close( omom )

  END SUBROUTINE fileio_close_mom


!--------------------------------------
  SUBROUTINE fileio_open_trn ( path )
!--------------------------------------

    character(*), intent(in) :: path

    character(6)   :: crank
    character(1)   :: srank
    character(3)   :: cnew

    if ( (zsp_rank /= 0) .OR. (vel_rank /= 0) ) return

    write( crank, fmt="(i6.6)" ) rankg
    write( srank, fmt="(i1.1)" ) ranks
    write( cnew,  fmt="(i3.3)" ) inum

    open( otrn, file=path//crank//"."//srank//".trn."//cnew, &
          form="unformatted" )

  END SUBROUTINE fileio_open_trn

!--------------------------------------
  SUBROUTINE fileio_close_trn
!--------------------------------------

    if ( (zsp_rank /= 0) .OR. (vel_rank /= 0) ) return

    close( otrn )

  END SUBROUTINE fileio_close_trn


!--------------------------------------
  SUBROUTINE fileio_open_tri ( path, cmx, cmy, replace )
!--------------------------------------

    character(*), intent(in) :: path
    character(*), intent(in) :: cmx, cmy
    logical, intent(in) :: replace

    character(1)   :: srank
    character(3)   :: cnew

    if ( rank /= 0 ) return

    write( srank, fmt="(i1.1)" ) ranks
    write( cnew,  fmt="(i3.3)" ) inum

    if ( replace ) then
       open( otri, file=path//"s"//srank//"mx"//cmx//"my"//cmy//".tri."//cnew, &
             form="unformatted", status="replace" )
    else
       open( otri, file=path//"s"//srank//"mx"//cmx//"my"//cmy//".tri."//cnew, &
             form="unformatted", status="unknown", position="append" )
    end if

  END SUBROUTINE fileio_open_tri

!--------------------------------------
  SUBROUTINE fileio_close_tri
!--------------------------------------

    if ( rank /= 0 ) return

    close( otri )

  END SUBROUTINE fileio_close_tri



!--------------------------------------
  SUBROUTINE fileio_read_cnt ( wf, time, istatus )
!--------------------------------------

    complex(kind=DP), intent(out), &
         dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: wf
    real(kind=DP), intent(out) :: time
    integer, optional, intent(out) :: istatus

    integer :: input_status

    read( unit=icnt, iostat=input_status ) time, wf
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

    rewind ocnt
    write( unit=ocnt ) time, wf

    call flush(ocnt)

  END SUBROUTINE fileio_write_cnt


!--------------------------------------
  SUBROUTINE fileio_write_fxv ( fout, time )
!--------------------------------------

    complex(kind=DP), intent(in), &
         dimension(-nx:nx,0:ny,1:2*nv,0:nm) :: fout
    real(kind=DP), intent(in) :: time

    write( unit=ofxv ) time, fout

    call flush(ofxv)

  END SUBROUTINE fileio_write_fxv


!--------------------------------------
  SUBROUTINE fileio_write_phi ( phi, time )
!--------------------------------------

    complex(kind=DP), intent(in), &
         dimension(-nx:nx,0:ny,-nz:nz-1) :: phi
    real(kind=DP), intent(in) :: time

    !- OUTPUT binary data phi/"*.phi.*"
    if ( ranks == 0 .AND. vel_rank == 0 ) then
       write( unit=ophi ) time, phi

       call flush(ophi)
    end if


  END SUBROUTINE fileio_write_phi


!--------------------------------------
  SUBROUTINE fileio_write_Al ( Al, time )
!--------------------------------------

    complex(kind=DP), intent(in), &
         dimension(-nx:nx,0:ny,-nz:nz-1) :: Al
    real(kind=DP), intent(in) :: time

    !- OUTPUT binary data phi/"*.Al.*"
    if ( ranks == 0 .AND. vel_rank == 0 ) then
       write( unit=oAl  ) time, Al

       call flush(oAl)
    end if

  END SUBROUTINE fileio_write_Al


!--------------------------------------
  SUBROUTINE fileio_write_mom ( dens, upara, ppara, pperp, qlpara, qlperp, time )
!--------------------------------------

    complex(kind=DP), intent(in), &
         dimension(-nx:nx,0:ny,-nz:nz-1) ::  dens, upara, ppara, pperp, qlpara, qlperp
    real(kind=DP), intent(in) :: time

    if ( vel_rank /= 0 ) return

    write( unit=omom ) time, dens, upara, ppara, pperp, qlpara, qlperp

    call flush(omom)

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

    write( unit=otrn ) time, entrpy, fenegy, menegy,    &
         peint, pmint, neint, nmint, dcd, &
         pflux_es, pflux_em, eflux_es, eflux_em

    call flush(otrn)

  END SUBROUTINE fileio_write_trn


!--------------------------------------
  SUBROUTINE fileio_write_tri ( jkpq_es, jpqk_es, jqkp_es, &
                                jkpq_em, jpqk_em, jqkp_em, time )
!--------------------------------------

    real(kind=DP), intent(in), &
         dimension(-nx:nx,-global_ny:global_ny) :: jkpq_es, jpqk_es, jqkp_es, &
         jkpq_em, jpqk_em, jqkp_em
    real(kind=DP), intent(in) :: time

    if ( rank /= 0 ) return

    write( unit=otri ) time, jkpq_es, jpqk_es, jqkp_es, jkpq_em, jpqk_em, jqkp_em

  END SUBROUTINE fileio_write_tri

END MODULE GKV_fileio
