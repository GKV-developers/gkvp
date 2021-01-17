MODULE GKV_zfilter
!-------------------------------------------------------------------------------
!
!    Filtering in zz to reduce high-kz numerical oscillations
!
!    Update history of gkvp_zfilter.f90
!    --------------
!      gkvp_f0.57 (S. Maeyama, Oct 2020)
!        - Version number f0.57 is removed from filename.
!
!-------------------------------------------------------------------------------

  use GKV_header
  use GKV_mpienv
  use GKV_clock, only: clock_sta, clock_end

  implicit none

  private

  public   zfilter


CONTAINS


!--------------------------------------
  SUBROUTINE zfilter ( vv )
!--------------------------------------
!     z-derivative of f

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: vv

    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: ww
    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: zb1_bottom, zb1_top
    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: zb2_bottom, zb2_top
    integer :: im

      allocate( ww(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1:2*nv,0:nm) )
      allocate( zb1_bottom(-nx:nx,0:ny, 0:nzb-1, 1:2*nv,0:nm) )
      allocate( zb1_top(-nx:nx,0:ny, 0:nzb-1, 1:2*nv,0:nm) )
      allocate( zb2_bottom(-nx:nx,0:ny, 0:nzb-1, 1:2*nv,0:nm) )
      allocate( zb2_top(-nx:nx,0:ny, 0:nzb-1, 1:2*nv,0:nm) )

!$OMP parallel default (none) &
!$OMP shared (vv,ww,zb1_bottom,zb1_top,zb2_bottom,zb2_top) &
!$OMP private (im)
      call zfilter_copy (         vv(:,:,:,:,0),      ww(:,:,:,:,0),  &
                          zb1_bottom(:,:,:,:,0), zb1_top(:,:,:,:,0),  &
                          zb2_bottom(:,:,:,:,0), zb2_top(:,:,:,:,0) )
!$OMP barrier

      im = 0
!$OMP master
        call zfilter_sendrecv ( zb1_bottom(:,:,:,:,im), zb1_top(:,:,:,:,im),  &
                                zb2_bottom(:,:,:,:,im), zb2_top(:,:,:,:,im) )
!$OMP end master
!!!!! barrier ! for test without overlaps
        call zfilter_copy (         vv(:,:,:,:,im+1),      ww(:,:,:,:,im+1),  &
                            zb1_bottom(:,:,:,:,im+1), zb1_top(:,:,:,:,im+1),  &
                            zb2_bottom(:,:,:,:,im+1), zb2_top(:,:,:,:,im+1) )
!$OMP barrier

      im = 1
!$OMP master
        call zfilter_sendrecv ( zb1_bottom(:,:,:,:,im), zb1_top(:,:,:,:,im),  &
                                zb2_bottom(:,:,:,:,im), zb2_top(:,:,:,:,im) )
!$OMP end master
!!!!! barrier ! for test without overlaps
        call zfilter_copy (         vv(:,:,:,:,im+1),      ww(:,:,:,:,im+1),  &
                            zb1_bottom(:,:,:,:,im+1), zb1_top(:,:,:,:,im+1),  &
                            zb2_bottom(:,:,:,:,im+1), zb2_top(:,:,:,:,im+1) )
        call zfilter_buffout ( zb2_bottom(:,:,:,:,im-1), zb2_top(:,:,:,:,im-1), &
                               ww(:,:,:,:,im-1) )
!$OMP barrier

      do im = 2, nm-1
!$OMP master
        call zfilter_sendrecv ( zb1_bottom(:,:,:,:,im), zb1_top(:,:,:,:,im),  &
                                zb2_bottom(:,:,:,:,im), zb2_top(:,:,:,:,im) )
!$OMP end master
!!!!! barrier ! for test without overlaps
        call zfilter_copy (         vv(:,:,:,:,im+1),      ww(:,:,:,:,im+1),  &
                            zb1_bottom(:,:,:,:,im+1), zb1_top(:,:,:,:,im+1),  &
                            zb2_bottom(:,:,:,:,im+1), zb2_top(:,:,:,:,im+1) )
        call zfilter_buffout ( zb2_bottom(:,:,:,:,im-1), zb2_top(:,:,:,:,im-1), &
                               ww(:,:,:,:,im-1) )
        call zfilter_filtering ( ww(:,:,:,:,im-2), vv(:,:,:,:,im-2) )
!$OMP barrier
      end do

      im = nm
!$OMP master
        call zfilter_sendrecv ( zb1_bottom(:,:,:,:,im), zb1_top(:,:,:,:,im),  &
                                zb2_bottom(:,:,:,:,im), zb2_top(:,:,:,:,im) )
!$OMP end master
!!!!! barrier ! for test without overlaps
        call zfilter_buffout ( zb2_bottom(:,:,:,:,im-1), zb2_top(:,:,:,:,im-1), &
                               ww(:,:,:,:,im-1) )
        call zfilter_filtering ( ww(:,:,:,:,im-2), vv(:,:,:,:,im-2) )
!$OMP barrier

        call zfilter_buffout ( zb2_bottom(:,:,:,:,nm), zb2_top(:,:,:,:,nm), &
                               ww(:,:,:,:,nm) )
        call zfilter_filtering ( ww(:,:,:,:,nm-1), vv(:,:,:,:,nm-1) )
!$OMP barrier

        call zfilter_filtering ( ww(:,:,:,:,nm), vv(:,:,:,:,nm) )
!$OMP end parallel

      deallocate( ww )
      deallocate( zb1_bottom )
      deallocate( zb1_top )
      deallocate( zb2_bottom )
      deallocate( zb2_top )

  END SUBROUTINE zfilter


!--------------------------------------
  SUBROUTINE zfilter_copy( vv, ww, zb1_bottom, zb1_top, zb2_bottom, zb2_top )
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv) :: vv
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1:2*nv) :: ww
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,0:nzb-1,1:2*nv) :: zb1_bottom, zb1_top
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,0:nzb-1,1:2*nv) :: zb2_bottom, zb2_top

    complex(kind=DP) :: wk
    integer  ::  mx, my, iz, iv


!$OMP master
                                           call clock_sta(1521)
                                         ! call fapp_start("zfilter_comm_bufferin",1521,1)
!$OMP end master
!$OMP do schedule (dynamic)
      do iv = 1, 2*nv
        do iz = 0, nzb-1
          do my = ist_y, iend_y
            do mx = -nx, nx
              wk                         = vv(mx,my,-nz+iz  ,iv)
              zb1_bottom(mx,my,iz,iv) = wk
                  ww(mx,my,-nz+iz,iv) = wk
            end do
          end do
        end do
        do iz = -nz+nzb, nz-1-nzb
          do my = ist_y, iend_y
            do mx = -nx, nx
                      ww(mx,my,iz,iv) = vv(mx,my,iz,iv)
            end do
          end do
        end do
        do iz = 0, nzb-1
          do my = ist_y, iend_y
            do mx = -nx, nx
              wk                         = vv(mx,my, nz-nzb+iz,iv)
              zb1_top   (mx,my,iz,iv) = wk
              ww(mx,my, nz-nzb+iz,iv) = wk
            end do
          end do
        end do
        do iz = 0, nzb-1
          do my = ist_y, iend_y
            do mx = -nx, nx
              zb2_bottom(mx,my,iz,iv)  = ( 0._DP, 0._DP )
              zb2_top   (mx,my,iz,iv)  = ( 0._DP, 0._DP )
            end do
          end do
        end do
      end do
!$OMP end do nowait
!$OMP master
                                         ! call fapp_stop("zfilter_comm_bufferin",1521,1)
                                           call clock_end(1521)
!$OMP end master


  END SUBROUTINE zfilter_copy


!--------------------------------------
  SUBROUTINE zfilter_sendrecv( zb1_bottom, zb1_top, zb2_bottom, zb2_top )
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,0:nzb-1,1:2*nv) :: zb1_bottom, zb1_top
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,0:nzb-1,1:2*nv) :: zb2_bottom, zb2_top

    integer  ::  slngz
    integer, dimension(4) :: ireq
    integer, dimension(MPI_STATUS_SIZE,4) :: istatus


      slngz  = (2*nx+1)*(ny+1)*(2*nv) * nzb

                                           call clock_sta(1522)
                                         ! call fapp_start("zfilter_comm_sendrecv",1522,1)
!      call MPI_sendrecv( zb1_bottom, slngz, MPI_DOUBLE_COMPLEX, izdn, 1, &
!                         zb2_top,    slngz, MPI_DOUBLE_COMPLEX, izup, 1, &
!                         sub_comm_world, status, ierr_mpi )
!
!      call MPI_sendrecv( zb1_top,    slngz, MPI_DOUBLE_COMPLEX, izup, 2, &
!                         zb2_bottom, slngz, MPI_DOUBLE_COMPLEX, izdn, 2, &
!                         sub_comm_world, status, ierr_mpi )

      call MPI_irecv( zb2_top,    slngz, MPI_DOUBLE_COMPLEX, izup, 1, &
                      sub_comm_world, ireq(1), ierr_mpi )
      call MPI_irecv( zb2_bottom, slngz, MPI_DOUBLE_COMPLEX, izdn, 2, &
                      sub_comm_world, ireq(2), ierr_mpi )
      call MPI_isend( zb1_bottom, slngz, MPI_DOUBLE_COMPLEX, izdn, 1, &
                      sub_comm_world, ireq(3), ierr_mpi )
      call MPI_isend( zb1_top,    slngz, MPI_DOUBLE_COMPLEX, izup, 2, &
                      sub_comm_world, ireq(4), ierr_mpi )
      call MPI_waitall( 4, ireq, istatus, ierr_mpi )
                                         ! call fapp_stop("zfilter_comm_sendrecv",1522,1)
                                           call clock_end(1522)


   END SUBROUTINE zfilter_sendrecv


!--------------------------------------
  SUBROUTINE zfilter_buffout( zb2_bottom, zb2_top, ww )
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,0:nzb-1,1:2*nv) :: zb2_bottom, zb2_top
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1:2*nv) :: ww

    integer  ::  mx, my, iz, iv, mwn, mwp


!$OMP master
                                           call clock_sta(1523)
                                         ! call fapp_start("zfilter_comm_bufferout",1523,1)
!$OMP end master
      if( rankz /= 0 ) then

!$OMP do schedule (dynamic)
        do iv = 1, 2*nv
          do iz = 0, nzb-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                ww(mx,my,-nz-nzb+iz,iv) = zb2_bottom(mx,my,iz,iv)
              end do
            end do
          end do
        end do
!$OMP end do nowait

      else  ! rankz==0

!$OMP do schedule (dynamic)
        do iv = 1, 2*nv
          do iz = 0, nzb-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                mwn   = mx + dj(my)          ! --- mw = mx + dj for the negative-z 

                if( abs(mwn) > nx ) then
                  ww(mx,my,-nz-nzb+iz,iv) = ( 0._DP, 0._DP )
                else
                  ww(mx,my,-nz-nzb+iz,iv) = ck(my) * zb2_bottom(mwn,my,iz,iv)
                end if

              end do
            end do
          end do
        end do
!$OMP end do nowait

      end if

      if( rankz /= nprocz-1 ) then

!$OMP do schedule (dynamic)
        do iv = 1, 2*nv
          do iz = 0, nzb-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                ww(mx,my,nz+iz,iv) = zb2_top(mx,my,iz,iv)
              end do
            end do
          end do
        end do
!$OMP end do nowait

      else ! rankz==nprocz-1

!$OMP do schedule (dynamic)
        do iv = 1, 2*nv
          do iz = 0, nzb-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                mwp   = mx - dj(my)          ! --- mw = mx - dj for the positive-z 

                if( abs(mwp) > nx ) then
                  ww(mx,my,nz+iz,iv) = ( 0._DP, 0._DP )
                else
                  ww(mx,my,nz+iz,iv) = conjg( ck(my) ) * zb2_top(mwp,my,iz,iv)
                end if

              end do
            end do
          end do
        end do
!$OMP end do nowait

      end if
!$OMP master
                                         ! call fapp_stop("zfilter_comm_bufferout",1523,1)
                                           call clock_end(1523)
!$OMP end master


  END SUBROUTINE zfilter_buffout


!--------------------------------------
  SUBROUTINE zfilter_filtering( ww, vv )
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1:2*nv) :: ww
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv)         :: vv

    real(kind=DP) :: alph, ceff
    integer  ::  mx, my, iz, iv


      alph   = 1._DP
      ceff   = alph / 16._DP

!$OMP master
                                             call clock_sta(1510)
                                           ! call fapp_start("zfilter_calc",1510,1)
!$OMP end master
!$OMP do schedule (dynamic)
      do iv = 1, 2*nv
        do iz = -nz, nz-1
          do my = ist_y, iend_y
            do mx = -nx, nx
              vv(mx,my,iz,iv) =                               &
                      ( 1._DP - alph )  * ww(mx,my,iz  ,iv)   &
                    + ceff * ( -          ww(mx,my,iz+2,iv)   &
                               +  4._DP * ww(mx,my,iz+1,iv)   &
                               + 10._DP * ww(mx,my,iz  ,iv)   &
                               +  4._DP * ww(mx,my,iz-1,iv)   &
                               -          ww(mx,my,iz-2,iv) )
            end do
          end do
        end do
      end do
!$OMP end do nowait
!$OMP master
                                           ! call fapp_stop("zfilter_calc",1510,1)
                                             call clock_end(1510)
!$OMP end master

  END SUBROUTINE zfilter_filtering


END MODULE GKV_zfilter
