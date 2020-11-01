MODULE GKV_bndry
!-------------------------------------------------------------------------------
!
!    Some useful tools and tips
!
!      GKV-plus r0.3 ( T.-H.Watanabe, Jun 2011)
!
!-------------------------------------------------------------------------------

  use GKV_header
  use GKV_mpienv
  use GKV_clock, only: clock_sta, clock_end

  implicit none

  private

  public   bndry_zvm_bound_f, bndry_bound_e,  &
      bndry_bound_f_buffin, bndry_bound_f_sendrecv, bndry_bound_f_buffout,  &
      bndry_shifts_v_buffin, bndry_shifts_v_sendrecv, bndry_shifts_v_buffout,  &
      bndry_zv_buffin, bndry_zv_sendrecv, bndry_zv_buffout, &
      bndry_vm_buffin, bndry_vm_sendrecv, bndry_vm_buffout, &
      bndry_shifts_m_buffin, bndry_shifts_m_sendrecv, bndry_shifts_m_buffout, &
      bndry_vm_sendrecv_v2, bndry_zv_buffin_v2, bndry_zv_sendrecv_v2, bndry_zv_buffout_v2


CONTAINS


!--------------------------------------
  SUBROUTINE bndry_zvm_bound_f( ff )
!--------------------------------------

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff

    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: zb1_bottom, zb1_top
    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: zb2_bottom, zb2_top
    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: vb1, vb2
    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: mb1, mb2
    integer :: im

      allocate( zb1_bottom(-nx:nx,0:ny,0:nzb-1,1:2*nv,0:nm) )
      allocate( zb1_top(-nx:nx,0:ny,0:nzb-1,1:2*nv,0:nm) )
      allocate( zb2_bottom(-nx:nx,0:ny,0:nzb-1,1:2*nv,0:nm) )
      allocate( zb2_top(-nx:nx,0:ny,0:nzb-1,1:2*nv,0:nm) )
      allocate( vb1(-nx:nx,0:ny,-nz:nz-1,1:2*nvb,0:nm) )
      allocate( vb2(-nx:nx,0:ny,-nz:nz-1,1:2*nvb,0:nm) )
      allocate( mb1(-nx:nx,0:ny,-nz:nz-1,1:2*nv,1:2*nvb) )
      allocate( mb2(-nx:nx,0:ny,-nz:nz-1,1:2*nv,1:2*nvb) )

!$OMP parallel default (none) &
!$OMP shared(ff,zb1_bottom,zb1_top,zb2_bottom,zb2_top,vb1,vb2,mb1,mb2) &
!$OMP private(im)
      do im = 0, nm
        call bndry_bound_f_buffin ( ff(:,:,:,:,im), zb1_bottom(:,:,:,:,im), zb1_top(:,:,:,:,im) )
!$OMP barrier
!$OMP master
        call bndry_bound_f_sendrecv ( zb1_bottom(:,:,:,:,im), zb1_top(:,:,:,:,im),  &
                                      zb2_bottom(:,:,:,:,im), zb2_top(:,:,:,:,im) )
!$OMP end master
!$OMP barrier
        call bndry_bound_f_buffout ( zb2_bottom(:,:,:,:,im), zb2_top(:,:,:,:,im), ff(:,:,:,:,im) )
      end do
!$OMP barrier

      do im = 0, nm
        call bndry_shifts_v_buffin ( ff(:,:,:,:,im), vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
!$OMP barrier
!$OMP master
        call bndry_shifts_v_sendrecv ( vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
!$OMP end master
!$OMP barrier
        call bndry_shifts_v_buffout ( vb2(:,:,:,:,im), ff(:,:,:,:,im) )
      end do
!$OMP barrier

      call bndry_shifts_m_buffin ( ff, mb1, mb2 )
!$OMP barrier
!$OMP master
      call bndry_shifts_m_sendrecv ( mb1, mb2 )
!$OMP end master
!$OMP barrier
      call bndry_shifts_m_buffout ( mb2, ff )
!$OMP end parallel

      deallocate( zb1_bottom )
      deallocate( zb1_top )
      deallocate( zb2_bottom )
      deallocate( zb2_top )
      deallocate( vb1 )
      deallocate( vb2 )
      deallocate( mb1 )
      deallocate( mb2 )

  END SUBROUTINE bndry_zvm_bound_f


!--------------------------------------
  SUBROUTINE bndry_bound_f_buffin( ff, zb1_bottom, zb1_top )
!--------------------------------------
!   Impose the modified periodic boundary condition 
!     in the z-direction for the distribution function

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb) :: ff
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,0:nzb-1,1:2*nv) :: zb1_bottom, zb1_top

    integer  ::  mx, my, iz, iv


!$OMP master
                                           call clock_sta(1351)
                                         ! call fapp_start("literm_boundf_bufferin",1351,1)
!$OMP end master

!$OMP do schedule (dynamic)
        do iv = 1, 2*nv
          do iz = 0, nzb-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                zb1_bottom(mx,my,iz,iv) = ff(mx,my,-nz+iz  ,iv)
                zb1_top   (mx,my,iz,iv) = ff(mx,my, nz-nzb+iz,iv)
              end do
            end do
          end do
        end do
!$OMP end do nowait

!$OMP master
                                         ! call fapp_stop("literm_boundf_bufferin",1351,1)
                                           call clock_end(1351)
!$OMP end master


  END SUBROUTINE bndry_bound_f_buffin


!--------------------------------------
  SUBROUTINE bndry_bound_f_sendrecv ( zb1_bottom, zb1_top, zb2_bottom, zb2_top )
!--------------------------------------
!   Impose the modified periodic boundary condition 
!     in the z-direction for the distribution function

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,0:nzb-1,1:2*nv) :: zb1_bottom, zb1_top
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,0:nzb-1,1:2*nv) :: zb2_bottom, zb2_top
    integer  ::  slngz
    integer, dimension(4) :: ireq
    integer, dimension(MPI_STATUS_SIZE,4) :: istatus


      slngz  = (2*nx+1)*(ny+1)*(2*nv) * nzb

                                           call clock_sta(1352)
                                         ! call fapp_start("literm_boundf_sendrecv",1352,1)
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
                                         ! call fapp_stop("literm_boundf_sendrecv",1352,1)
                                           call clock_end(1352)


  END SUBROUTINE bndry_bound_f_sendrecv


!--------------------------------------
  SUBROUTINE bndry_bound_f_buffout ( zb2_bottom, zb2_top, ff )
!--------------------------------------
!   Impose the modified periodic boundary condition 
!     in the z-direction for the distribution function

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,0:nzb-1,1:2*nv) :: zb2_bottom, zb2_top
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb) :: ff

    integer  ::  mx, my, iz, iv, mwn, mwp


! --- substitution
!$OMP master
                                           call clock_sta(1353)
                                         ! call fapp_start("literm_boundf_bufferout",1353,1)
!$OMP end master

      if( rankz /= 0 ) then

!$OMP do schedule (dynamic)
          do iv = 1, 2*nv
            do iz = 0, nzb-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  ff(mx,my,-nz-nzb+iz,iv) = zb2_bottom(mx,my,iz,iv)
                end do
              end do
            end do
          end do
!$OMP end do nowait

      else  ! rankz==0

        if ( trim(z_bound) == "outflow" .OR. trim(z_bound) == "mixed") then

!$OMP do schedule (dynamic)
          do iv = 1, 2*nv
              do my = ist_y, iend_y
                do mx = -nx, nx
                  mwn   = mx + dj(my)          ! --- mw = mx + dj for the negative-z 

                  if( abs(mwn) > nx ) then
                    if ( vl(iv) > 0._DP ) then ! inflow
                      do iz = 0, nzb-1
                        ff(mx,my,-nz-nzb+iz,iv) = ( 0._DP, 0._DP )
                      end do
                    else                       ! outflow
                      ff(mx,my,-nz-1,iv) =   ff(mx,my,-nz  ,iv)
                      ff(mx,my,-nz-2,iv) = - ff(mx,my,-nz+1,iv) + 2._DP * ff(mx,my,-nz  ,iv)
                    end if
                  else
                    do iz = 0, nzb-1
                      ff(mx,my,-nz-nzb+iz,iv) = ck(my) * zb2_bottom(mwn,my,iz,iv)
                    end do
                  end if

                end do
              end do
          end do
!$OMP end do nowait

        else if ( trim(z_bound) == "zerofixed" ) then

!$OMP do schedule (dynamic)
          do iv = 1, 2*nv
            do iz = 0, nzb-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  mwn   = mx + dj(my)          ! --- mw = mx + dj for the negative-z 

                  if( abs(mwn) > nx ) then
                    ff(mx,my,-nz-nzb+iz,iv) = ( 0._DP, 0._DP )
                  else
                    ff(mx,my,-nz-nzb+iz,iv) = ck(my) * zb2_bottom(mwn,my,iz,iv)
                  end if

                end do
              end do
            end do
          end do
!$OMP end do nowait

        else

          write( olog, * ) " # z_bound is to be  outflow  or  zerofixed"
          call flush(olog)
          stop

        end if

      end if

      if( rankz /= nprocz-1 ) then

!$OMP do schedule (dynamic)
          do iv = 1, 2*nv
            do iz = 0, nzb-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  ff(mx,my,nz+iz,iv) = zb2_top(mx,my,iz,iv)
                end do
              end do
            end do
          end do
!$OMP end do nowait

      else ! rankz==nprocz-1

        if ( trim(z_bound) == "outflow" .OR. trim(z_bound) == "mixed") then

!$OMP do schedule (dynamic)
          do iv = 1, 2*nv
              do my = ist_y, iend_y
                do mx = -nx, nx
                  mwp   = mx - dj(my)          ! --- mw = mx - dj for the positive-z 

                  if( abs(mwp) > nx ) then
                    if ( vl(iv) > 0._DP ) then ! outflow
                      ff(mx,my,nz  ,iv) =   ff(mx,my,nz-1,iv)
                      ff(mx,my,nz+1,iv) = - ff(mx,my,nz-2,iv) + 2._DP * ff(mx,my,nz-1,iv)
                    else                       ! inflow
                      do iz = 0, nzb-1
                        ff(mx,my,nz+iz,iv) = ( 0._DP, 0._DP )
                      end do
                    end if
                  else
                    do iz = 0, nzb-1
                      ff(mx,my,nz+iz,iv) = conjg( ck(my) ) * zb2_top(mwp,my,iz,iv)
                    end do
                  end if

                end do
              end do
          end do
!$OMP end do nowait

        else if ( trim(z_bound) == "zerofixed" ) then

!$OMP do schedule (dynamic)
          do iv = 1, 2*nv
            do iz = 0, nzb-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  mwp   = mx - dj(my)          ! --- mw = mx - dj for the positive-z 

                  if( abs(mwp) > nx ) then
                    ff(mx,my,nz+iz,iv) = ( 0._DP, 0._DP )
                  else
                    ff(mx,my,nz+iz,iv) = conjg( ck(my) ) * zb2_top(mwp,my,iz,iv)
                  end if

                end do
              end do
            end do
          end do
!$OMP end do nowait

        else

          write( olog, * ) " # z_bound is to be  outflow  or  zerofixed"
          call flush(olog)
          stop

        end if

      end if

!$OMP master
                                         ! call fapp_stop("literm_boundf_bufferout",1353,1)
                                           call clock_end(1353)
!$OMP end master


  END SUBROUTINE bndry_bound_f_buffout



!--------------------------------------
  SUBROUTINE bndry_shifts_v_buffin( ff, vb1, vb2 )
!--------------------------------------
!     Shift communications in v and m directions

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb) :: ff
    complex(kind=DP), intent(out),  &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nvb) :: vb1, vb2

    integer  ::  mx, my, iz, iv


!$OMP master
                                           call clock_sta(1361)
                                         ! call fapp_start("literm_shifts_bufferin",1361,1)
!$OMP end master
! --- zero clear is required for rankv = 0, nprocv-1 and rankm = 0, nprocm-1
      do iv = 1, 2*nvb
!$OMP do schedule (dynamic)
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                vb2(mx,my,iz,iv) = ( 0._DP, 0._DP )
              end do
            end do
          end do
!$OMP end do nowait
      end do

!$OMP do schedule (dynamic)
        do iz = -nz, nz-1
          do my = ist_y, iend_y
            do mx = -nx, nx
              do iv = 1, nvb
                vb1(mx,my,iz,iv    ) = ff(mx,my,iz,         iv)
                vb1(mx,my,iz,iv+nvb) = ff(mx,my,iz,2*nv-nvb+iv)
!              vb1(mx,my,iz,1) = ff(mx,my,iz,     1)
!              vb1(mx,my,iz,2) = ff(mx,my,iz,     2)
!              vb1(mx,my,iz,3) = ff(mx,my,iz,     3)
!              vb1(mx,my,iz,4) = ff(mx,my,iz,2*nv-2)
!              vb1(mx,my,iz,5) = ff(mx,my,iz,2*nv-1)
!              vb1(mx,my,iz,6) = ff(mx,my,iz,2*nv  )
              end do
            end do
          end do
        end do
!$OMP end do nowait


!$OMP master
                                         ! call fapp_stop("literm_shifts_bufferin",1361,1)
                                           call clock_end(1361)
!$OMP end master


  END SUBROUTINE bndry_shifts_v_buffin


!--------------------------------------
  SUBROUTINE bndry_shifts_v_sendrecv( vb1, vb2 )
!--------------------------------------
!     Shift communications in v and m directions

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nvb) :: vb1
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nvb) :: vb2
    integer  ::  slngv
    integer, dimension(4) :: ireq
    integer, dimension(MPI_STATUS_SIZE,4) :: istatus


      slngv = (2*nx+1)*(ny+1)*(2*nz) * nvb

                                           call clock_sta(1362)
                                         ! call fapp_start("literm_shifts_sendrecv",1362,1)
!      call MPI_sendrecv( vb1(-nx,0,-nz,1    ), slngv, MPI_DOUBLE_COMPLEX, ivdn, 1, &
!                         vb2(-nx,0,-nz,nvb+1), slngv, MPI_DOUBLE_COMPLEX, ivup, 1, &
!                         sub_comm_world, status, ierr_mpi )
!
!      call MPI_sendrecv( vb1(-nx,0,-nz,nvb+1), slngv, MPI_DOUBLE_COMPLEX, ivup, 2, &
!                         vb2(-nx,0,-nz,1    ), slngv, MPI_DOUBLE_COMPLEX, ivdn, 2, &
!                         sub_comm_world, status, ierr_mpi )

      call MPI_irecv( vb2(-nx,0,-nz,nvb+1), slngv, MPI_DOUBLE_COMPLEX, ivup, 1, &
                      sub_comm_world, ireq(1), ierr_mpi )
      call MPI_irecv( vb2(-nx,0,-nz,    1), slngv, MPI_DOUBLE_COMPLEX, ivdn, 2, &
                      sub_comm_world, ireq(2), ierr_mpi )
      call MPI_isend( vb1(-nx,0,-nz,    1), slngv, MPI_DOUBLE_COMPLEX, ivdn, 1, &
                      sub_comm_world, ireq(3), ierr_mpi )
      call MPI_isend( vb1(-nx,0,-nz,nvb+1), slngv, MPI_DOUBLE_COMPLEX, ivup, 2, &
                      sub_comm_world, ireq(4), ierr_mpi )
      call MPI_waitall( 4, ireq, istatus, ierr_mpi )
                                         ! call fapp_stop("literm_shifts_sendrecv",1362,1)
                                           call clock_end(1362)


  END SUBROUTINE bndry_shifts_v_sendrecv


!--------------------------------------
  SUBROUTINE bndry_shifts_v_buffout( vb2, ff )
!--------------------------------------
!     Shift communications in v and m directions

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nvb) :: vb2
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb) :: ff

    integer  ::  mx, my, iz, iv


!$OMP master
                                           call clock_sta(1363)
                                         ! call fapp_start("literm_shifts_bufferout",1363,1)
!$OMP end master

!$OMP do schedule (dynamic)
        do iz = -nz, nz-1
          do my = ist_y, iend_y
            do mx = -nx, nx
              do iv = 1, nvb
                ff(mx,my,iz,-nvb+iv) = vb2(mx,my,iz,iv    )
                ff(mx,my,iz,2*nv+iv) = vb2(mx,my,iz,iv+nvb)
!              ff(mx,my,iz,    -2) = vb2(mx,my,iz,1)
!              ff(mx,my,iz,    -1) = vb2(mx,my,iz,2)
!              ff(mx,my,iz,     0) = vb2(mx,my,iz,3)
!              ff(mx,my,iz,2*nv+1) = vb2(mx,my,iz,4)
!              ff(mx,my,iz,2*nv+2) = vb2(mx,my,iz,5)
!              ff(mx,my,iz,2*nv+3) = vb2(mx,my,iz,6)
              end do
            end do
          end do
        end do
!$OMP end do nowait


!$OMP master
                                         ! call fapp_stop("literm_shifts_bufferout",1363,1)
                                           call clock_end(1363)
!$OMP end master


  END SUBROUTINE bndry_shifts_v_buffout


!--------------------------------------
  SUBROUTINE bndry_shifts_m_buffin( ff, mb1, mb2 )
!--------------------------------------
!     Shift communications in v and m directions

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,1:2*nvb) :: mb1, mb2

    integer  ::  mx, my, iz, iv, im


!$OMP master
                                           call clock_sta(1371)
                                         ! call fapp_start("literm_shifts_bufferin",1371,1)
!$OMP end master

! --- zero clear is required for rankv = 0, nprocv-1 and rankm = 0, nprocm-1
      do im = 1, 2*nvb
!$OMP do schedule (dynamic)
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                mb2(mx,my,iz,iv,im) = ( 0._DP, 0._DP )
              end do
            end do
          end do
        end do
!$OMP end do nowait
      end do

!$OMP do schedule (dynamic)
      do iv = 1, 2*nv
        do iz = -nz, nz-1
          do my = ist_y, iend_y
            do mx = -nx, nx
              do im = 1, nvb
                mb1(mx,my,iz,iv,im    ) = ff(mx,my,iz,iv,     im-1)
                mb1(mx,my,iz,iv,im+nvb) = ff(mx,my,iz,iv,nm-nvb+im)
!              mb1(mx,my,iz,iv,1) = ff(mx,my,iz,iv,   0)
!              mb1(mx,my,iz,iv,2) = ff(mx,my,iz,iv,   1)
!              mb1(mx,my,iz,iv,3) = ff(mx,my,iz,iv,   2)
!              mb1(mx,my,iz,iv,4) = ff(mx,my,iz,iv,nm-2)
!              mb1(mx,my,iz,iv,5) = ff(mx,my,iz,iv,nm-1)
!              mb1(mx,my,iz,iv,6) = ff(mx,my,iz,iv,nm  )
              end do
            end do
          end do
        end do
      end do
!$OMP end do nowait


!$OMP master
                                         ! call fapp_stop("literm_shifts_bufferin",1371,1)
                                           call clock_end(1371)
!$OMP end master


  END SUBROUTINE bndry_shifts_m_buffin


!--------------------------------------
  SUBROUTINE bndry_shifts_m_sendrecv( mb1, mb2 )
!--------------------------------------
!     Shift communications in v and m directions

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,1:2*nvb) :: mb1
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,1:2*nvb) :: mb2
    integer  ::  slngm
    integer, dimension(4) :: ireq
    integer, dimension(MPI_STATUS_SIZE,4) :: istatus


      slngm = (2*nx+1)*(ny+1)*(2*nz)*(2*nv) * nvb

                                           call clock_sta(1372)
                                         ! call fapp_start("literm_shifts_sendrecv",1372,1)
!      call MPI_sendrecv( mb1(-nx,0,-nz,1,1    ), slngm, MPI_DOUBLE_COMPLEX, imdn, 3, &
!                         mb2(-nx,0,-nz,1,nvb+1), slngm, MPI_DOUBLE_COMPLEX, imup, 3, &
!                         sub_comm_world, status, ierr_mpi )
!
!      call MPI_sendrecv( mb1(-nx,0,-nz,1,nvb+1), slngm, MPI_DOUBLE_COMPLEX, imup, 4, &
!                         mb2(-nx,0,-nz,1,1    ), slngm, MPI_DOUBLE_COMPLEX, imdn, 4, &
!                         sub_comm_world, status, ierr_mpi )

      call MPI_irecv( mb2(-nx,0,-nz,1,nvb+1), slngm, MPI_DOUBLE_COMPLEX, imup, 3, &
                      sub_comm_world, ireq(1), ierr_mpi )
      call MPI_irecv( mb2(-nx,0,-nz,1,    1), slngm, MPI_DOUBLE_COMPLEX, imdn, 4, &
                      sub_comm_world, ireq(2), ierr_mpi )
      call MPI_isend( mb1(-nx,0,-nz,1,    1), slngm, MPI_DOUBLE_COMPLEX, imdn, 3, &
                      sub_comm_world, ireq(3), ierr_mpi )
      call MPI_isend( mb1(-nx,0,-nz,1,nvb+1), slngm, MPI_DOUBLE_COMPLEX, imup, 4, &
                      sub_comm_world, ireq(4), ierr_mpi )
      call MPI_waitall( 4, ireq, istatus, ierr_mpi )
                                         ! call fapp_stop("literm_shifts_sendrecv",1372,1)
                                           call clock_end(1372)


  END SUBROUTINE bndry_shifts_m_sendrecv


!--------------------------------------
  SUBROUTINE bndry_shifts_m_buffout( mb2, ff )
!--------------------------------------
!     Shift communications in v and m directions

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,1:2*nvb) :: mb2
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff

    integer  ::  mx, my, iz, iv, im


!$OMP master
                                           call clock_sta(1373)
                                         ! call fapp_start("literm_shifts_bufferout",1373,1)
!$OMP end master

!$OMP do schedule (dynamic)
      do iv = 1, 2*nv
        do iz = -nz, nz-1
          do my = ist_y, iend_y
            do mx = -nx, nx
              do im = 1, nvb
                ff(mx,my,iz,iv,-nvb-1+im) = mb2(mx,my,iz,iv,im    )
                ff(mx,my,iz,iv,nm+im    ) = mb2(mx,my,iz,iv,im+nvb)
!              ff(mx,my,iz,iv,  -3) = mb2(mx,my,iz,iv,1)
!              ff(mx,my,iz,iv,  -2) = mb2(mx,my,iz,iv,2)
!              ff(mx,my,iz,iv,  -1) = mb2(mx,my,iz,iv,3)
!              ff(mx,my,iz,iv,nm+1) = mb2(mx,my,iz,iv,4)
!              ff(mx,my,iz,iv,nm+2) = mb2(mx,my,iz,iv,5)
!              ff(mx,my,iz,iv,nm+3) = mb2(mx,my,iz,iv,6)
              end do
            end do
          end do
        end do
      end do
!$OMP end do nowait

!$OMP master
                                         ! call fapp_stop("literm_shifts_bufferout",1373,1)
                                           call clock_end(1373)
!$OMP end master


  END SUBROUTINE bndry_shifts_m_buffout


!--------------------------------------
  SUBROUTINE bndry_bound_e ( ew )
!--------------------------------------
!  Impose the modified periodic boundary condition 
!    in the z-direction for the electric field

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,0:nm)   :: ew

    complex(kind=DP), dimension(:,:,:,:), allocatable :: zb1e_bottom, zb1e_top
    complex(kind=DP), dimension(:,:,:,:), allocatable :: zb2e_bottom, zb2e_top
    integer  ::  mx, my, iz, im, mwn, mwp
    integer  ::  slngze
    integer, dimension(4) :: ireq
    integer, dimension(MPI_STATUS_SIZE,4) :: istatus


      allocate( zb1e_bottom(-nx:nx,0:ny,0:nzb-1,0:nm) )
      allocate( zb1e_top(-nx:nx,0:ny,0:nzb-1,0:nm) )
      allocate( zb2e_bottom(-nx:nx,0:ny,0:nzb-1,0:nm) )
      allocate( zb2e_top(-nx:nx,0:ny,0:nzb-1,0:nm) )

      slngze  = (2*nx+1)*(ny+1)*(nm+1) * nzb

!$OMP parallel default(none) &
!$OMP shared(zb2e_bottom,zb2e_top,zb1e_bottom,zb1e_top,ist_y,iend_y,ew) &
!$OMP private(mx,my,iz,im)
!$OMP master
                                           call clock_sta(1381)
                                         ! call fapp_start("literm_bounde_bufferin",1381,1)
!$OMP end master

!$OMP do schedule (dynamic)
      do im = 0, nm
        do iz = 0, nzb-1
          do my = ist_y, iend_y
            do mx = -nx, nx
              zb2e_bottom(mx,my,iz,im) = ( 0._DP, 0._DP )
              zb2e_top   (mx,my,iz,im) = ( 0._DP, 0._DP )
            end do
          end do
        end do
      end do
!$OMP end do nowait

!$OMP do schedule (dynamic)
      do im = 0, nm
        do iz = 0, nzb-1
          do my = ist_y, iend_y
            do mx = -nx, nx
              zb1e_bottom(mx,my,iz,im) = ew(mx,my,-nz+iz  ,im)
              zb1e_top   (mx,my,iz,im) = ew(mx,my, nz-nzb+iz,im)
            end do
          end do
        end do
      end do
!$OMP end do

!$OMP master
                                         ! call fapp_stop("literm_bounde_bufferin",1381,1)
                                           call clock_end(1381)
!$OMP end master
!$OMP end parallel

                                           call clock_sta(1382)
                                         ! call fapp_start("literm_bounde_sendrecv",1382,1)
!      call MPI_sendrecv( zb1e_bottom, slngze, MPI_DOUBLE_COMPLEX, izdn, 1, &
!                         zb2e_top,    slngze, MPI_DOUBLE_COMPLEX, izup, 1, &
!                         sub_comm_world, status, ierr_mpi )
!
!      call MPI_sendrecv( zb1e_top,    slngze, MPI_DOUBLE_COMPLEX, izup, 2, &
!                         zb2e_bottom, slngze, MPI_DOUBLE_COMPLEX, izdn, 2, &
!                         sub_comm_world, status, ierr_mpi )

      call MPI_irecv( zb2e_top,    slngze, MPI_DOUBLE_COMPLEX, izup, 1, &
                      sub_comm_world, ireq(1), ierr_mpi )
      call MPI_irecv( zb2e_bottom, slngze, MPI_DOUBLE_COMPLEX, izdn, 2, &
                      sub_comm_world, ireq(2), ierr_mpi )
      call MPI_isend( zb1e_bottom, slngze, MPI_DOUBLE_COMPLEX, izdn, 1, &
                      sub_comm_world, ireq(3), ierr_mpi )
      call MPI_isend( zb1e_top,    slngze, MPI_DOUBLE_COMPLEX, izup, 2, &
                      sub_comm_world, ireq(4), ierr_mpi )
      call MPI_waitall( 4, ireq, istatus, ierr_mpi )
                                         ! call fapp_stop("literm_bounde_sendrecv",1382,1)
                                           call clock_end(1382)

! --- substitution
!$OMP parallel default(none) &
!$OMP shared(zb2e_bottom,zb2e_top,zb1e_bottom,zb1e_top,ist_y,iend_y,ew) &
!$OMP shared(rankz,z_bound,ck,dj) &
!$OMP private(mx,my,iz,im,mwp,mwn)
!$OMP master
                                           call clock_sta(1383)
                                         ! call fapp_start("literm_bounde_bufferout",1383,1)
!$OMP end master

      if( rankz /= 0 ) then

!$OMP do schedule (dynamic)
        do im = 0, nm
          do iz = 0, nzb-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                ew(mx,my,-nz-nzb+iz,im) = zb2e_bottom(mx,my,iz,im)
              end do
            end do
          end do
        end do
!$OMP end do

      else   ! rankz==0

        if ( trim(z_bound) == "outflow" ) then

!$OMP do schedule (dynamic)
          do im = 0, nm
              do my = ist_y, iend_y
                do mx = -nx, nx
                  mwn   = mx + dj(my)          ! --- mw = mx + dj for the negative-z 

                  if( abs(mwn) > nx ) then
                    ew(mx,my,-nz-1,im)   =   ew(mx,my,-nz  ,im)
                    ew(mx,my,-nz-2,im)   = - ew(mx,my,-nz+1,im) + 2._DP * ew(mx,my,-nz  ,im)
                  else
                    do iz = 0, nzb-1
                      ew(mx,my,-nz-nzb+iz,im) = ck(my) * zb2e_bottom(mwn,my,iz,im)
                    end do
                  end if

                end do
              end do
          end do
!$OMP end do

        else if ( trim(z_bound) == "zerofixed" .OR. trim(z_bound) == "mixed" ) then

!$OMP do schedule (dynamic)
          do im = 0, nm
            do iz = 0, nzb-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  mwn   = mx + dj(my)          ! --- mw = mx + dj for the negative-z 

                  if( abs(mwn) > nx ) then
                    ew(mx,my,-nz-nzb+iz,im)   = ( 0._DP, 0._DP )
                  else
                    ew(mx,my,-nz-nzb+iz,im) = ck(my) * zb2e_bottom(mwn,my,iz,im)
                  end if

                end do
              end do
            end do
          end do
!$OMP end do

        else

          write( olog, * ) " # z_bound is to be  outflow  or  zerofixed"
          call flush(olog)
          stop

        end if

      end if

      if( rankz /= nprocz-1 ) then

!$OMP do schedule (dynamic)
        do im = 0, nm
          do iz = 0, nzb-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                ew(mx,my,nz+iz,im) = zb2e_top(mx,my,iz,im)
              end do
            end do
          end do
        end do
!$OMP end do

      else   ! rankz==nprocz-1

        if ( trim(z_bound) == "outflow" ) then

!$OMP do schedule (dynamic)
          do im = 0, nm
              do my = ist_y, iend_y
                do mx = -nx, nx
                  mwp   = mx - dj(my)          ! --- mw = mx - dj for the positive-z 

                  if( abs(mwp) > nx ) then
                    ew(mx,my,nz  ,im)   =   ew(mx,my,nz-1,im)
                    ew(mx,my,nz+1,im)   = - ew(mx,my,nz-2,im) + 2._DP * ew(mx,my,nz-1,im) 
                  else
                    do iz = 0, nzb-1
                      ew(mx,my,nz+iz,im) = conjg( ck(my) ) * zb2e_top(mwp,my,iz,im)
                    end do
                  end if

                end do
              end do
          end do
!$OMP end do

        else if ( trim(z_bound) == "zerofixed" .OR. trim(z_bound) == "mixed" ) then

!$OMP do schedule (dynamic)
          do im = 0, nm
            do iz = 0, nzb-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  mwp   = mx - dj(my)          ! --- mw = mx - dj for the positive-z 

                  if( abs(mwp) > nx ) then
                    ew(mx,my,nz+iz,im)   = ( 0._DP, 0._DP )
                  else
                    ew(mx,my,nz+iz,im) = conjg( ck(my) ) * zb2e_top(mwp,my,iz,im)
                  end if

                end do
              end do
            end do
          end do
!$OMP end do

        else

          write( olog, * ) " # z_bound is to be  outflow  or  zerofixed"
          call flush(olog)
          stop

        end if

      end if

!$OMP master
                                         ! call fapp_stop("literm_bounde_bufferout",1383,1)
                                           call clock_end(1383)
!$OMP end master
!$OMP end parallel

      deallocate( zb1e_bottom )
      deallocate( zb1e_top )
      deallocate( zb2e_bottom )
      deallocate( zb2e_top )


  END SUBROUTINE bndry_bound_e


!--------------------------------------
  SUBROUTINE bndry_zv_buffin( ff, zb1_bottom, zb1_top, vb1 )
!--------------------------------------
!   Impose the modified periodic boundary condition 
!     in the z-direction for the distribution function

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb) :: ff
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,0:nzb-1,1:2*nv) :: zb1_bottom, zb1_top
    complex(kind=DP), intent(out),  &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nvb) :: vb1

    integer :: iz, iv


!$OMP master
                                           call clock_sta(1351)
                                         ! call fapp_start("literm_boundf_bufferin",1351,1)
!$OMP end master

!$OMP do collapse(2) schedule(dynamic)
        do iv = 1, 2*nv
          do iz = 0, nzb-1
            zb1_bottom(:,:,iz,iv) = ff(:,:,-nz+iz  ,iv)
            zb1_top   (:,:,iz,iv) = ff(:,:, nz-nzb+iz,iv)
          end do
        end do
!$OMP end do nowait


!! --- zero clear is required for rankv = 0, nprocv-1 and rankm = 0, nprocm-1
!      do iv = 1, 2*nvb
!!$OMP do schedule(dynamic)
!          do iz = -nz, nz-1
!            do my = ist_y, iend_y
!              do mx = -nx, nx
!                vb2(mx,my,iz,iv) = ( 0._DP, 0._DP )
!              end do
!            end do
!          end do
!!$OMP end do nowait
!      end do

!$OMP do collapse(2) schedule(dynamic)
        do iv = 1, nvb
          do iz = -nz, nz-1
            vb1(:,:,iz,iv    ) = ff(:,:,iz,         iv)
            vb1(:,:,iz,iv+nvb) = ff(:,:,iz,2*nv-nvb+iv)
          end do
        end do
!$OMP end do nowait

!$OMP master
                                         ! call fapp_stop("literm_boundf_bufferin",1351,1)
                                           call clock_end(1351)
!$OMP end master


  END SUBROUTINE bndry_zv_buffin


!--------------------------------------
  SUBROUTINE bndry_zv_sendrecv ( zb1_bottom, zb1_top, zb2_bottom, zb2_top, vb1, vb2 )
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,0:nzb-1,1:2*nv) :: zb1_bottom, zb1_top
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,0:nzb-1,1:2*nv) :: zb2_bottom, zb2_top
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nvb) :: vb1
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nvb) :: vb2

    integer :: slngz, slngv
    integer, dimension(8) :: ireq
    integer, dimension(MPI_STATUS_SIZE,8) :: istatus


      slngz  = (2*nx+1)*(ny+1)*(2*nv) * nzb
      slngv = (2*nx+1)*(ny+1)*(2*nz) * nvb

                                           call clock_sta(1352)
                                         ! call fapp_start("literm_boundf_sendrecv",1352,1)
     !call MPI_sendrecv( zb1_bottom, slngz, MPI_DOUBLE_COMPLEX, izdn, 1, &
     !                   zb2_top,    slngz, MPI_DOUBLE_COMPLEX, izup, 1, &
     !                   sub_comm_world, status, ierr_mpi )
     !call MPI_sendrecv( zb1_top,    slngz, MPI_DOUBLE_COMPLEX, izup, 2, &
     !                   zb2_bottom, slngz, MPI_DOUBLE_COMPLEX, izdn, 2, &
     !                   sub_comm_world, status, ierr_mpi )
     !call MPI_sendrecv( vb1(-nx,0,-nz,1    ), slngv, MPI_DOUBLE_COMPLEX, ivdn, 3, &
     !                   vb2(-nx,0,-nz,nvb+1), slngv, MPI_DOUBLE_COMPLEX, ivup, 3, &
     !                   sub_comm_world, status, ierr_mpi )
     !call MPI_sendrecv( vb1(-nx,0,-nz,nvb+1), slngv, MPI_DOUBLE_COMPLEX, ivup, 4, &
     !                   vb2(-nx,0,-nz,1    ), slngv, MPI_DOUBLE_COMPLEX, ivdn, 4, &
     !                   sub_comm_world, status, ierr_mpi )

      call MPI_irecv( zb2_top,    slngz, MPI_DOUBLE_COMPLEX, izup, 1, &
                      sub_comm_world, ireq(1), ierr_mpi )
      call MPI_irecv( zb2_bottom, slngz, MPI_DOUBLE_COMPLEX, izdn, 2, &
                      sub_comm_world, ireq(2), ierr_mpi )
      call MPI_irecv( vb2(-nx,0,-nz,nvb+1), slngv, MPI_DOUBLE_COMPLEX, ivup, 3, &
                      sub_comm_world, ireq(3), ierr_mpi )
      call MPI_irecv( vb2(-nx,0,-nz,    1), slngv, MPI_DOUBLE_COMPLEX, ivdn, 4, &
                      sub_comm_world, ireq(4), ierr_mpi )
      call MPI_isend( zb1_bottom, slngz, MPI_DOUBLE_COMPLEX, izdn, 1, &
                      sub_comm_world, ireq(5), ierr_mpi )
      call MPI_isend( zb1_top,    slngz, MPI_DOUBLE_COMPLEX, izup, 2, &
                      sub_comm_world, ireq(6), ierr_mpi )
      call MPI_isend( vb1(-nx,0,-nz,    1), slngv, MPI_DOUBLE_COMPLEX, ivdn, 3, &
                      sub_comm_world, ireq(7), ierr_mpi )
      call MPI_isend( vb1(-nx,0,-nz,nvb+1), slngv, MPI_DOUBLE_COMPLEX, ivup, 4, &
                      sub_comm_world, ireq(8), ierr_mpi )
      call MPI_waitall( 8, ireq, istatus, ierr_mpi )
                                         ! call fapp_stop("literm_boundf_sendrecv",1352,1)
                                           call clock_end(1352)


  END SUBROUTINE bndry_zv_sendrecv


!--------------------------------------
  SUBROUTINE bndry_zv_buffout ( zb2_bottom, zb2_top, vb2, ff )
!--------------------------------------
!   Impose the modified periodic boundary condition 
!     in the z-direction for the distribution function

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,0:nzb-1,1:2*nv) :: zb2_bottom, zb2_top
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nvb) :: vb2
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb) :: ff

    integer :: mx, my, iz, iv, mwn, mwp


! --- substitution
!$OMP master
                                           call clock_sta(1353)
                                         ! call fapp_start("literm_boundf_bufferout",1353,1)
!$OMP end master

      if( rankz /= 0 ) then

!$OMP do collapse(2) schedule(dynamic)
          do iv = 1, 2*nv
            do iz = 0, nzb-1
              ff(:,:,-nz-nzb+iz,iv) = zb2_bottom(:,:,iz,iv)
            end do
          end do
!$OMP end do nowait

      else  ! rankz==0

        if ( trim(z_bound) == "outflow" .OR. trim(z_bound) == "mixed") then

!$OMP do collapse(2) schedule(dynamic)
          do iv = 1, 2*nv
              do my = ist_y, iend_y
                do mx = -nx, nx
                  mwn   = mx + dj(my)          ! --- mw = mx + dj for the negative-z 

                  if( abs(mwn) > nx ) then
                    if ( vl(iv) > 0._DP ) then ! inflow
                      do iz = 0, nzb-1
                        ff(mx,my,-nz-nzb+iz,iv) = ( 0._DP, 0._DP )
                      end do
                    else                       ! outflow
                      ff(mx,my,-nz-1,iv) =   ff(mx,my,-nz  ,iv)
                      ff(mx,my,-nz-2,iv) = - ff(mx,my,-nz+1,iv) + 2._DP * ff(mx,my,-nz  ,iv)
                    end if
                  else
                    do iz = 0, nzb-1
                      ff(mx,my,-nz-nzb+iz,iv) = ck(my) * zb2_bottom(mwn,my,iz,iv)
                    end do
                  end if

                end do
              end do
          end do
!$OMP end do nowait

        else if ( trim(z_bound) == "zerofixed" ) then

!$OMP do collapse(2) schedule(dynamic)
          do iv = 1, 2*nv
            do iz = 0, nzb-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  mwn   = mx + dj(my)          ! --- mw = mx + dj for the negative-z 

                  if( abs(mwn) > nx ) then
                    ff(mx,my,-nz-nzb+iz,iv) = ( 0._DP, 0._DP )
                  else
                    ff(mx,my,-nz-nzb+iz,iv) = ck(my) * zb2_bottom(mwn,my,iz,iv)
                  end if

                end do
              end do
            end do
          end do
!$OMP end do nowait

        else

          write( olog, * ) " # z_bound is to be  outflow  or  zerofixed"
          call flush(olog)
          stop

        end if

      end if

      if( rankz /= nprocz-1 ) then

!$OMP do collapse(2) schedule(dynamic)
          do iv = 1, 2*nv
            do iz = 0, nzb-1
              ff(:,:,nz+iz,iv) = zb2_top(:,:,iz,iv)
            end do
          end do
!$OMP end do nowait

      else ! rankz==nprocz-1

        if ( trim(z_bound) == "outflow" .OR. trim(z_bound) == "mixed") then

!$OMP do collapse(2) schedule(dynamic)
          do iv = 1, 2*nv
              do my = ist_y, iend_y
                do mx = -nx, nx
                  mwp   = mx - dj(my)          ! --- mw = mx - dj for the positive-z 

                  if( abs(mwp) > nx ) then
                    if ( vl(iv) > 0._DP ) then ! outflow
                      ff(mx,my,nz  ,iv) =   ff(mx,my,nz-1,iv)
                      ff(mx,my,nz+1,iv) = - ff(mx,my,nz-2,iv) + 2._DP * ff(mx,my,nz-1,iv)
                    else                       ! inflow
                      do iz = 0, nzb-1
                        ff(mx,my,nz+iz,iv) = ( 0._DP, 0._DP )
                      end do
                    end if
                  else
                    do iz = 0, nzb-1
                      ff(mx,my,nz+iz,iv) = conjg( ck(my) ) * zb2_top(mwp,my,iz,iv)
                    end do
                  end if

                end do
              end do
          end do
!$OMP end do nowait

        else if ( trim(z_bound) == "zerofixed" ) then

!$OMP do collapse(2) schedule(dynamic)
          do iv = 1, 2*nv
            do iz = 0, nzb-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  mwp   = mx - dj(my)          ! --- mw = mx - dj for the positive-z 

                  if( abs(mwp) > nx ) then
                    ff(mx,my,nz+iz,iv) = ( 0._DP, 0._DP )
                  else
                    ff(mx,my,nz+iz,iv) = conjg( ck(my) ) * zb2_top(mwp,my,iz,iv)
                  end if

                end do
              end do
            end do
          end do
!$OMP end do nowait

        else

          write( olog, * ) " # z_bound is to be  outflow  or  zerofixed"
          call flush(olog)
          stop

        end if

      end if


        if ( rankv == 0 ) then
!$OMP do collapse(2) schedule(dynamic)
          do iv = 1, nvb
            do iz = -nz, nz-1
              ff(:,:,iz,-nvb+iv) = (0._DP, 0._DP)
              ff(:,:,iz,2*nv+iv) = vb2(:,:,iz,iv+nvb)
            end do
          end do
!$OMP end do nowait
        else if ( rankv == nprocv-1 ) then
!$OMP do collapse(2) schedule(dynamic)
          do iv = 1, nvb
            do iz = -nz, nz-1
              ff(:,:,iz,-nvb+iv) = vb2(:,:,iz,iv    )
              ff(:,:,iz,2*nv+iv) = (0._DP, 0._DP)
            end do
          end do
!$OMP end do nowait
        else
!$OMP do collapse(2) schedule(dynamic)
          do iv = 1, nvb
            do iz = -nz, nz-1
              ff(:,:,iz,-nvb+iv) = vb2(:,:,iz,iv    )
              ff(:,:,iz,2*nv+iv) = vb2(:,:,iz,iv+nvb)
            end do
          end do
!$OMP end do nowait
        end if

!$OMP master
                                         ! call fapp_stop("literm_boundf_bufferout",1353,1)
                                           call clock_end(1353)
!$OMP end master


  END SUBROUTINE bndry_zv_buffout


!--------------------------------------
  SUBROUTINE bndry_vm_buffin( iz, ff, vb1, mb1 )
!--------------------------------------
!     Shift communications in v and m directions

    integer, intent(in) :: iz
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(out),  &
      dimension(-nx:nx,0:ny,0:nm,1:2*nvb) :: vb1
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,1:2*nv,1:2*nvb) :: mb1

    integer :: mx, my, iv, im


!$OMP master
                                           call clock_sta(1371)
                                         ! call fapp_start("literm_shifts_bufferin",1371,1)
!$OMP end master

!! --- zero clear is required for rankv = 0, nprocv-1 and rankm = 0, nprocm-1
!      do iv = 1, 2*nvb
!        do im = 0, nm
!            do my = ist_y, iend_y
!              do mx = -nx, nx
!                vb2(mx,my,im,iv) = ( 0._DP, 0._DP )
!              end do
!            end do
!        end do
!      end do

!$OMP do collapse(2) schedule(dynamic)
        do iv = 1, nvb
          do im = 0, nm
            do my = ist_y, iend_y
              do mx = -nx, nx
                vb1(mx,my,im,iv    ) = ff(mx,my,iz,         iv,im)
                vb1(mx,my,im,iv+nvb) = ff(mx,my,iz,2*nv-nvb+iv,im)
              end do
            end do
          end do
        end do
!$OMP end do nowait

!! --- zero clear is required for rankv = 0, nprocv-1 and rankm = 0, nprocm-1
!      do im = 1, 2*nvb
!        do iv = 1, 2*nv
!            do my = ist_y, iend_y
!              do mx = -nx, nx
!                mb2(mx,my,iv,im) = ( 0._DP, 0._DP )
!              end do
!            end do
!        end do
!      end do

!$OMP do collapse(2) schedule(dynamic)
      do im = 1, nvb
        do iv = 1, 2*nv
          do my = ist_y, iend_y
            do mx = -nx, nx
              mb1(mx,my,iv,im    ) = ff(mx,my,iz,iv,     im-1)
              mb1(mx,my,iv,im+nvb) = ff(mx,my,iz,iv,nm-nvb+im)
            end do
          end do
        end do
      end do
!$OMP end do nowait


!$OMP master
                                         ! call fapp_stop("literm_shifts_bufferin",1371,1)
                                           call clock_end(1371)
!$OMP end master


  END SUBROUTINE bndry_vm_buffin


!--------------------------------------
  SUBROUTINE bndry_vm_sendrecv ( vb1, vb2, mb1, mb2 )
!--------------------------------------
!     Shift communications in v and m directions

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,0:nm,1:2*nvb) :: vb1
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,0:nm,1:2*nvb) :: vb2
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,1:2*nv,1:2*nvb) :: mb1
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,1:2*nv,1:2*nvb) :: mb2

    integer :: slngv, slngm
    integer, dimension(8) :: ireq
    integer, dimension(MPI_STATUS_SIZE,8) :: istatus


      slngv = (2*nx+1)*(ny+1)*(nm+1) * nvb
      slngm = (2*nx+1)*(ny+1)*(2*nv) * nvb

                                           call clock_sta(1372)
                                         ! call fapp_start("literm_shifts_sendrecv",1372,1)
     !call MPI_sendrecv( vb1(-nx,0,0,1    ), slngv, MPI_DOUBLE_COMPLEX, ivdn, 1, &
     !                   vb2(-nx,0,0,nvb+1), slngv, MPI_DOUBLE_COMPLEX, ivup, 1, &
     !                   sub_comm_world, status, ierr_mpi )
     !call MPI_sendrecv( vb1(-nx,0,0,nvb+1), slngv, MPI_DOUBLE_COMPLEX, ivup, 2, &
     !                   vb2(-nx,0,0,1    ), slngv, MPI_DOUBLE_COMPLEX, ivdn, 2, &
     !                   sub_comm_world, status, ierr_mpi )
     !call MPI_sendrecv( mb1(-nx,0,1,1    ), slngm, MPI_DOUBLE_COMPLEX, imdn, 3, &
     !                   mb2(-nx,0,1,nvb+1), slngm, MPI_DOUBLE_COMPLEX, imup, 3, &
     !                   sub_comm_world, status, ierr_mpi )
     !call MPI_sendrecv( mb1(-nx,0,1,nvb+1), slngm, MPI_DOUBLE_COMPLEX, imup, 4, &
     !                   mb2(-nx,0,1,1    ), slngm, MPI_DOUBLE_COMPLEX, imdn, 4, &
     !                   sub_comm_world, status, ierr_mpi )

      call MPI_irecv( vb2(-nx,0,0,nvb+1), slngv, MPI_DOUBLE_COMPLEX, ivup, 1, &
                      sub_comm_world, ireq(1), ierr_mpi )
      call MPI_irecv( vb2(-nx,0,0,    1), slngv, MPI_DOUBLE_COMPLEX, ivdn, 2, &
                      sub_comm_world, ireq(2), ierr_mpi )
      call MPI_irecv( mb2(-nx,0,1,nvb+1), slngm, MPI_DOUBLE_COMPLEX, imup, 3, &
                      sub_comm_world, ireq(3), ierr_mpi )
      call MPI_irecv( mb2(-nx,0,1,    1), slngm, MPI_DOUBLE_COMPLEX, imdn, 4, &
                      sub_comm_world, ireq(4), ierr_mpi )
      call MPI_isend( vb1(-nx,0,0,    1), slngv, MPI_DOUBLE_COMPLEX, ivdn, 1, &
                      sub_comm_world, ireq(5), ierr_mpi )
      call MPI_isend( vb1(-nx,0,0,nvb+1), slngv, MPI_DOUBLE_COMPLEX, ivup, 2, &
                      sub_comm_world, ireq(6), ierr_mpi )
      call MPI_isend( mb1(-nx,0,1,    1), slngm, MPI_DOUBLE_COMPLEX, imdn, 3, &
                      sub_comm_world, ireq(7), ierr_mpi )
      call MPI_isend( mb1(-nx,0,1,nvb+1), slngm, MPI_DOUBLE_COMPLEX, imup, 4, &
                      sub_comm_world, ireq(8), ierr_mpi )
      call MPI_waitall( 8, ireq, istatus, ierr_mpi )
                                         ! call fapp_stop("literm_shifts_sendrecv",1372,1)
                                           call clock_end(1372)

  END SUBROUTINE bndry_vm_sendrecv


!--------------------------------------
  SUBROUTINE bndry_vm_buffout ( iz, vb2, mb2, ff )
!--------------------------------------
!     Shift communications in v and m directions

    integer, intent(in) :: iz
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,0:nm,1:2*nvb) :: vb2
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,1:2*nv,1:2*nvb) :: mb2
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff

    integer :: mx, my, iv, im


! --- substitution
!$OMP master
                                           call clock_sta(1373)
                                         ! call fapp_start("literm_shifts_bufferout",1373,1)
!$OMP end master

      if ( rankv == 0 ) then
!$OMP do collapse(2) schedule(dynamic)
        do iv = 1, nvb
          do im = 0, nm
            do my = ist_y, iend_y
              do mx = -nx, nx
                ff(mx,my,iz,-nvb+iv,im) = (0._DP, 0._DP)
                ff(mx,my,iz,2*nv+iv,im) = vb2(mx,my,im,iv+nvb)
              end do
            end do
          end do
        end do
!$OMP end do nowait
      else if ( rankv == nprocv-1 ) then
!$OMP do collapse(2) schedule(dynamic)
        do iv = 1, nvb
          do im = 0, nm
            do my = ist_y, iend_y
              do mx = -nx, nx
                ff(mx,my,iz,-nvb+iv,im) = vb2(mx,my,im,iv    )
                ff(mx,my,iz,2*nv+iv,im) = (0._DP, 0._DP)
              end do
            end do
          end do
        end do
!$OMP end do nowait
      else
!$OMP do collapse(2) schedule(dynamic)
        do iv = 1, nvb
          do im = 0, nm
            do my = ist_y, iend_y
              do mx = -nx, nx
                ff(mx,my,iz,-nvb+iv,im) = vb2(mx,my,im,iv    )
                ff(mx,my,iz,2*nv+iv,im) = vb2(mx,my,im,iv+nvb)
              end do
            end do
          end do
        end do
!$OMP end do nowait
      end if


      if ( rankm == 0 ) then
!$OMP do collapse(2) schedule(dynamic)
        do im = 1, nvb
          do iv = 1, 2*nv
            do my = ist_y, iend_y
              do mx = -nx, nx
                ff(mx,my,iz,iv,-nvb-1+im) = (0._DP, 0._DP)
                ff(mx,my,iz,iv,nm+im    ) = mb2(mx,my,iv,im+nvb)
              end do
            end do
          end do
        end do
!$OMP end do nowait
      else if ( rankm == nprocm-1 ) then
!$OMP do collapse(2) schedule(dynamic)
        do im = 1, nvb
          do iv = 1, 2*nv
            do my = ist_y, iend_y
              do mx = -nx, nx
                ff(mx,my,iz,iv,-nvb-1+im) = mb2(mx,my,iv,im    )
                ff(mx,my,iz,iv,nm+im    ) = (0._DP, 0._DP)
              end do
            end do
          end do
        end do
!$OMP end do nowait
      else
!$OMP do collapse(2) schedule(dynamic)
        do im = 1, nvb
          do iv = 1, 2*nv
            do my = ist_y, iend_y
              do mx = -nx, nx
                ff(mx,my,iz,iv,-nvb-1+im) = mb2(mx,my,iv,im    )
                ff(mx,my,iz,iv,nm+im    ) = mb2(mx,my,iv,im+nvb)
              end do
            end do
          end do
        end do
!$OMP end do nowait
      end if

!$OMP master
                                         ! call fapp_stop("literm_shifts_bufferout",1373,1)
                                           call clock_end(1373)
!$OMP end master


  END SUBROUTINE bndry_vm_buffout


!--------------------------------------
  SUBROUTINE bndry_vm_sendrecv_v2 ( vb1, vb2, mb1, mb2 )
!--------------------------------------
!     Shift communications in v and m directions

!mae> modify
!    complex(kind=DP), intent(in), &
!      dimension(-nx:nx,0:ny,-nz:nz-1,0:nm,1:2*nvb,0:nm) :: vb1
!    complex(kind=DP), intent(out), &
!      dimension(-nx:nx,0:ny,-nz:nz-1,0:nm,1:2*nvb,0:nm) :: vb2
!    complex(kind=DP), intent(in), &
!      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,1:2*nvb,0:nm) :: mb1
!    complex(kind=DP), intent(out), &
!      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,1:2*nvb,0:nm) :: mb2
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,0:nm,1:2*nvb) :: vb1
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,0:nm,1:2*nvb) :: vb2
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,1:2*nvb) :: mb1
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,1:2*nvb) :: mb2
!<

    integer :: slngv, slngm
    integer, dimension(8) :: ireq
    integer, dimension(MPI_STATUS_SIZE,8) :: istatus


      slngv = (2*nx+1)*(ny+1)*(2*nz)*(nm+1) * nvb
      slngm = (2*nx+1)*(ny+1)*(2*nz)*(2*nv) * nvb

                                           call clock_sta(1372)
                                         ! call fapp_start("literm_shifts_sendrecv",1372,1)
     !call MPI_sendrecv( vb1(-nx,0,0,1    ), slngv, MPI_DOUBLE_COMPLEX, ivdn, 1, &
     !                   vb2(-nx,0,0,nvb+1), slngv, MPI_DOUBLE_COMPLEX, ivup, 1, &
     !                   sub_comm_world, status, ierr_mpi )
     !call MPI_sendrecv( vb1(-nx,0,0,nvb+1), slngv, MPI_DOUBLE_COMPLEX, ivup, 2, &
     !                   vb2(-nx,0,0,1    ), slngv, MPI_DOUBLE_COMPLEX, ivdn, 2, &
     !                   sub_comm_world, status, ierr_mpi )
     !call MPI_sendrecv( mb1(-nx,0,1,1    ), slngm, MPI_DOUBLE_COMPLEX, imdn, 3, &
     !                   mb2(-nx,0,1,nvb+1), slngm, MPI_DOUBLE_COMPLEX, imup, 3, &
     !                   sub_comm_world, status, ierr_mpi )
     !call MPI_sendrecv( mb1(-nx,0,1,nvb+1), slngm, MPI_DOUBLE_COMPLEX, imup, 4, &
     !                   mb2(-nx,0,1,1    ), slngm, MPI_DOUBLE_COMPLEX, imdn, 4, &
     !                   sub_comm_world, status, ierr_mpi )

      call MPI_irecv( vb2(-nx,0,-nz,0,nvb+1), slngv, MPI_DOUBLE_COMPLEX, ivup, 1, &
                      sub_comm_world, ireq(1), ierr_mpi )
      call MPI_irecv( vb2(-nx,0,-nz,0,    1), slngv, MPI_DOUBLE_COMPLEX, ivdn, 2, &
                      sub_comm_world, ireq(2), ierr_mpi )
      call MPI_irecv( mb2(-nx,0,-nz,1,nvb+1), slngm, MPI_DOUBLE_COMPLEX, imup, 3, &
                      sub_comm_world, ireq(3), ierr_mpi )
      call MPI_irecv( mb2(-nx,0,-nz,1,    1), slngm, MPI_DOUBLE_COMPLEX, imdn, 4, &
                      sub_comm_world, ireq(4), ierr_mpi )
      call MPI_isend( vb1(-nx,0,-nz,0,    1), slngv, MPI_DOUBLE_COMPLEX, ivdn, 1, &
                      sub_comm_world, ireq(5), ierr_mpi )
      call MPI_isend( vb1(-nx,0,-nz,0,nvb+1), slngv, MPI_DOUBLE_COMPLEX, ivup, 2, &
                      sub_comm_world, ireq(6), ierr_mpi )
      call MPI_isend( mb1(-nx,0,-nz,1,    1), slngm, MPI_DOUBLE_COMPLEX, imdn, 3, &
                      sub_comm_world, ireq(7), ierr_mpi )
      call MPI_isend( mb1(-nx,0,-nz,1,nvb+1), slngm, MPI_DOUBLE_COMPLEX, imup, 4, &
                      sub_comm_world, ireq(8), ierr_mpi )
      call MPI_waitall( 8, ireq, istatus, ierr_mpi )
                                         ! call fapp_stop("literm_shifts_sendrecv",1372,1)
                                           call clock_end(1372)

  END SUBROUTINE bndry_vm_sendrecv_v2


!--------------------------------------
  SUBROUTINE bndry_zv_buffin_v2( ff, zb1_bottom, zb1_top, vb1 )
!--------------------------------------
!   Impose the modified periodic boundary condition 
!     in the z-direction for the distribution function

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,0:nzb-1,1:2*nv,0:nm) :: zb1_bottom, zb1_top
    complex(kind=DP), intent(out),  &
      dimension(-nx:nx,0:ny,-nz:nz-1,0:nm,1:2*nvb) :: vb1

    integer :: iz, iv, im


!$OMP master
                                           call clock_sta(1351)
                                         ! call fapp_start("literm_boundf_bufferin",1351,1)
!$OMP end master

!!TBI!! !$OMP do collapse(2) schedule(dynamic)
        do im = 0, nm
        do iv = 1, 2*nv
          do iz = 0, nzb-1
            zb1_bottom(:,:,iz,iv,im) = ff(:,:,-nz+iz    ,iv,im)
            zb1_top   (:,:,iz,iv,im) = ff(:,:, nz-nzb+iz,iv,im)
          end do
        end do
        end do
!!TBI!! !$OMP end do nowait


!! --- zero clear is required for rankv = 0, nprocv-1 and rankm = 0, nprocm-1
!      do iv = 1, 2*nvb
!!$OMP do schedule(dynamic)
!          do iz = -nz, nz-1
!            do my = ist_y, iend_y
!              do mx = -nx, nx
!                vb2(mx,my,iz,iv) = ( 0._DP, 0._DP )
!              end do
!            end do
!          end do
!!$OMP end do nowait
!      end do

!!TBI!! !$OMP do collapse(2) schedule(dynamic)
        do iv = 1, nvb
        do im = 0, nm
          do iz = -nz, nz-1
            vb1(:,:,iz,im,iv    ) = ff(:,:,iz,         iv,im)
            vb1(:,:,iz,im,iv+nvb) = ff(:,:,iz,2*nv-nvb+iv,im)
          end do
        end do
        end do
!!TBI!! !$OMP end do nowait

!$OMP master
                                         ! call fapp_stop("literm_boundf_bufferin",1351,1)
                                           call clock_end(1351)
!$OMP end master


  END SUBROUTINE bndry_zv_buffin_v2


!--------------------------------------
  SUBROUTINE bndry_zv_sendrecv_v2 ( zb1_bottom, zb1_top, zb2_bottom, zb2_top, vb1, vb2 )
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,0:nzb-1,1:2*nv,0:nm) :: zb1_bottom, zb1_top
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,0:nzb-1,1:2*nv,0:nm) :: zb2_bottom, zb2_top
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,0:nm,1:2*nvb) :: vb1
    complex(kind=DP), intent(out), &     
      dimension(-nx:nx,0:ny,-nz:nz-1,0:nm,1:2*nvb) :: vb2

    integer :: slngz, slngv
    integer, dimension(8) :: ireq
    integer, dimension(MPI_STATUS_SIZE,8) :: istatus


      slngz = (2*nx+1)*(ny+1)*(2*nv)*(nm+1) * nzb
      slngv = (2*nx+1)*(ny+1)*(2*nz)*(nm+1) * nvb

                                           call clock_sta(1352)
                                         ! call fapp_start("literm_boundf_sendrecv",1352,1)
     !call MPI_sendrecv( zb1_bottom, slngz, MPI_DOUBLE_COMPLEX, izdn, 1, &
     !                   zb2_top,    slngz, MPI_DOUBLE_COMPLEX, izup, 1, &
     !                   sub_comm_world, status, ierr_mpi )
     !call MPI_sendrecv( zb1_top,    slngz, MPI_DOUBLE_COMPLEX, izup, 2, &
     !                   zb2_bottom, slngz, MPI_DOUBLE_COMPLEX, izdn, 2, &
     !                   sub_comm_world, status, ierr_mpi )
     !call MPI_sendrecv( vb1(-nx,0,-nz,1    ), slngv, MPI_DOUBLE_COMPLEX, ivdn, 3, &
     !                   vb2(-nx,0,-nz,nvb+1), slngv, MPI_DOUBLE_COMPLEX, ivup, 3, &
     !                   sub_comm_world, status, ierr_mpi )
     !call MPI_sendrecv( vb1(-nx,0,-nz,nvb+1), slngv, MPI_DOUBLE_COMPLEX, ivup, 4, &
     !                   vb2(-nx,0,-nz,1    ), slngv, MPI_DOUBLE_COMPLEX, ivdn, 4, &
     !                   sub_comm_world, status, ierr_mpi )

      call MPI_irecv( zb2_top,    slngz, MPI_DOUBLE_COMPLEX, izup, 1, &
                      sub_comm_world, ireq(1), ierr_mpi )
      call MPI_irecv( zb2_bottom, slngz, MPI_DOUBLE_COMPLEX, izdn, 2, &
                      sub_comm_world, ireq(2), ierr_mpi )
      call MPI_irecv( vb2(-nx,0,-nz,0,nvb+1), slngv, MPI_DOUBLE_COMPLEX, ivup, 3, &
                      sub_comm_world, ireq(3), ierr_mpi )
      call MPI_irecv( vb2(-nx,0,-nz,0,    1), slngv, MPI_DOUBLE_COMPLEX, ivdn, 4, &
                      sub_comm_world, ireq(4), ierr_mpi )
      call MPI_isend( zb1_bottom, slngz, MPI_DOUBLE_COMPLEX, izdn, 1, &
                      sub_comm_world, ireq(5), ierr_mpi )
      call MPI_isend( zb1_top,    slngz, MPI_DOUBLE_COMPLEX, izup, 2, &
                      sub_comm_world, ireq(6), ierr_mpi )
      call MPI_isend( vb1(-nx,0,-nz,0,    1), slngv, MPI_DOUBLE_COMPLEX, ivdn, 3, &
                      sub_comm_world, ireq(7), ierr_mpi )
      call MPI_isend( vb1(-nx,0,-nz,0,nvb+1), slngv, MPI_DOUBLE_COMPLEX, ivup, 4, &
                      sub_comm_world, ireq(8), ierr_mpi )
      call MPI_waitall( 8, ireq, istatus, ierr_mpi )
                                         ! call fapp_stop("literm_boundf_sendrecv",1352,1)
                                           call clock_end(1352)


  END SUBROUTINE bndry_zv_sendrecv_v2


!--------------------------------------
  SUBROUTINE bndry_zv_buffout_v2 ( zb2_bottom, zb2_top, vb2, ff )
!--------------------------------------
!   Impose the modified periodic boundary condition 
!     in the z-direction for the distribution function

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,0:nzb-1,1:2*nv,0:nm) :: zb2_bottom, zb2_top
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,0:nm,1:2*nvb) :: vb2
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff

    integer :: mx, my, iz, iv, im, mwn, mwp


! --- substitution
!$OMP master
                                           call clock_sta(1353)
                                         ! call fapp_start("literm_boundf_bufferout",1353,1)
!$OMP end master

      if( rankz /= 0 ) then

!!TBI!! !$OMP do collapse(2) schedule(dynamic)
          do im = 0, nm
          do iv = 1, 2*nv
            do iz = 0, nzb-1
              ff(:,:,-nz-nzb+iz,iv,im) = zb2_bottom(:,:,iz,iv,im)
            end do
          end do
          end do
!!TBI!! !$OMP end do nowait

      else  ! rankz==0

        if ( trim(z_bound) == "outflow" .OR. trim(z_bound) == "mixed") then

!!TBI!! !$OMP do collapse(2) schedule(dynamic)
          do im = 0, nm
          do iv = 1, 2*nv
              do my = ist_y, iend_y
                do mx = -nx, nx
                  mwn   = mx + dj(my)          ! --- mw = mx + dj for the negative-z 

                  if( abs(mwn) > nx ) then
                    if ( vl(iv) > 0._DP ) then ! inflow
                      do iz = 0, nzb-1
                        ff(mx,my,-nz-nzb+iz,iv,im) = ( 0._DP, 0._DP )
                      end do
                    else                       ! outflow
                      ff(mx,my,-nz-1,iv,im) =   ff(mx,my,-nz  ,iv,im)
                      ff(mx,my,-nz-2,iv,im) = - ff(mx,my,-nz+1,iv,im) + 2._DP * ff(mx,my,-nz  ,iv,im)
                    end if
                  else
                    do iz = 0, nzb-1
                      ff(mx,my,-nz-nzb+iz,iv,im) = ck(my) * zb2_bottom(mwn,my,iz,iv,im)
                    end do
                  end if

                end do
              end do
          end do
          end do
!!TBI!! !$OMP end do nowait

        else if ( trim(z_bound) == "zerofixed" ) then

!!TBI!! !$OMP do collapse(2) schedule(dynamic)
          do im = 0, nm
          do iv = 1, 2*nv
            do iz = 0, nzb-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  mwn   = mx + dj(my)          ! --- mw = mx + dj for the negative-z 

                  if( abs(mwn) > nx ) then
                    ff(mx,my,-nz-nzb+iz,iv,im) = ( 0._DP, 0._DP )
                  else
                    ff(mx,my,-nz-nzb+iz,iv,im) = ck(my) * zb2_bottom(mwn,my,iz,iv,im)
                  end if

                end do
              end do
            end do
          end do
          end do
!!TBI!! !$OMP end do nowait

        else

          write( olog, * ) " # z_bound is to be  outflow  or  zerofixed"
          call flush(olog)
          stop

        end if

      end if

      if( rankz /= nprocz-1 ) then

!!TBI!! !$OMP do collapse(2) schedule(dynamic)
          do im = 0, nm
          do iv = 1, 2*nv
            do iz = 0, nzb-1
              ff(:,:,nz+iz,iv,im) = zb2_top(:,:,iz,iv,im)
            end do
          end do
          end do
!!TBI!! !$OMP end do nowait

      else ! rankz==nprocz-1

        if ( trim(z_bound) == "outflow" .OR. trim(z_bound) == "mixed") then

!!TBI!! !$OMP do collapse(2) schedule(dynamic)
          do im = 0, nm
          do iv = 1, 2*nv
              do my = ist_y, iend_y
                do mx = -nx, nx
                  mwp   = mx - dj(my)          ! --- mw = mx - dj for the positive-z 

                  if( abs(mwp) > nx ) then
                    if ( vl(iv) > 0._DP ) then ! outflow
                      ff(mx,my,nz  ,iv,im) =   ff(mx,my,nz-1,iv,im)
                      ff(mx,my,nz+1,iv,im) = - ff(mx,my,nz-2,iv,im) + 2._DP * ff(mx,my,nz-1,iv,im)
                    else                       ! inflow
                      do iz = 0, nzb-1
                        ff(mx,my,nz+iz,iv,im) = ( 0._DP, 0._DP )
                      end do
                    end if
                  else
                    do iz = 0, nzb-1
                      ff(mx,my,nz+iz,iv,im) = conjg( ck(my) ) * zb2_top(mwp,my,iz,iv,im)
                    end do
                  end if

                end do
              end do
          end do
          end do
!!TBI!! !$OMP end do nowait

        else if ( trim(z_bound) == "zerofixed" ) then

!!TBI!! !$OMP do collapse(2) schedule(dynamic)
          do im = 0, nm
          do iv = 1, 2*nv
            do iz = 0, nzb-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  mwp   = mx - dj(my)          ! --- mw = mx - dj for the positive-z 

                  if( abs(mwp) > nx ) then
                    ff(mx,my,nz+iz,iv,im) = ( 0._DP, 0._DP )
                  else
                    ff(mx,my,nz+iz,iv,im) = conjg( ck(my) ) * zb2_top(mwp,my,iz,iv,im)
                  end if

                end do
              end do
            end do
          end do
          end do
!!TBI!! !$OMP end do nowait

        else

          write( olog, * ) " # z_bound is to be  outflow  or  zerofixed"
          call flush(olog)
          stop

        end if

      end if


        if ( rankv == 0 ) then
!!TBI!! !$OMP do collapse(2) schedule(dynamic)
          do iv = 1, nvb
          do im = 0, nm
            do iz = -nz, nz-1
              ff(:,:,iz,-nvb+iv,im) = (0._DP, 0._DP)
              ff(:,:,iz,2*nv+iv,im) = vb2(:,:,iz,im,iv+nvb)
            end do
          end do
          end do
!!TBI!! !$OMP end do nowait
        else if ( rankv == nprocv-1 ) then
!!TBI!! !$OMP do collapse(2) schedule(dynamic)
          do iv = 1, nvb
          do im = 0, nm
            do iz = -nz, nz-1
              ff(:,:,iz,-nvb+iv,im) = vb2(:,:,iz,im,iv    )
              ff(:,:,iz,2*nv+iv,im) = (0._DP, 0._DP)
            end do
          end do
          end do
!!TBI!! !$OMP end do nowait
        else
!!TBI!! !$OMP do collapse(2) schedule(dynamic)
          do iv = 1, nvb
          do im = 0, nm
            do iz = -nz, nz-1
              ff(:,:,iz,-nvb+iv,im) = vb2(:,:,iz,im,iv    )
              ff(:,:,iz,2*nv+iv,im) = vb2(:,:,iz,im,iv+nvb)
            end do
          end do
          end do
!!TBI!! !$OMP end do nowait
        end if

!$OMP master
                                         ! call fapp_stop("literm_boundf_bufferout",1353,1)
                                           call clock_end(1353)
!$OMP end master


  END SUBROUTINE bndry_zv_buffout_v2


END MODULE GKV_bndry
