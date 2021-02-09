MODULE GKV_tips
!-------------------------------------------------------------------------------
!
!    Some useful tools and tips
!
!    Update history of gkvp_tips.f90
!    --------------
!      gkvp_f0.60 (S. Maeyama, Jan 2021)
!        - flush for binary is removed, for the usage of NetCDF.
!      gkvp_f0.57 (S. Maeyama, Oct 2020)
!        - Version number f0.57 is removed from filename.
!
!-------------------------------------------------------------------------------

  use GKV_header
  use GKV_mpienv

  implicit none

  private

  public   tips_reality, tips_flush, tips_rescale_for_linear_runs


CONTAINS


!--------------------------------------
  SUBROUTINE tips_reality( wrk )
!--------------------------------------

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: wrk

    integer :: mx

      if( rankw == 0 )  then
        do mx = 0, nx
          wrk(-mx,0,:,:,:) = conjg( wrk(mx,0,:,:,:) )
        end do
      endif


  END SUBROUTINE tips_reality


!--------------------------------------
  SUBROUTINE tips_flush
!--------------------------------------

        call flush(olog)
        !fj start 202010
        !call flush(ocnt)
        !call flush(ofxv)
        !if ( vel_rank == 0 ) then
        !  call flush(omom)
        !end if
        !if ( ranks == 0 .AND. vel_rank == 0 ) then
        !  call flush(ophi)
        !  call flush(oAl)
        !end if
        !if ( zsp_rank == 0 .AND. vel_rank == 0 ) then
        !  call flush(otrn)
        !end if
        !fj end 202010
        if( rankg == 0 ) then
          call flush(odtc)
          call flush(oeng)
          call flush(omen)
          call flush(owes)
          call flush(owem)
          if ( trim(calc_type) == "lin_freq" ) then
            call flush(ofrq)
          end if
        end if
        if( rank == 0 ) then
          call flush(obln)
          call flush(oges)
          call flush(ogem)
          call flush(oqes)
          call flush(oqem)
        end if

  END SUBROUTINE tips_flush


!--------------------------------------
  SUBROUTINE tips_rescale_for_linear_runs(ff, phi, Al, hh, time)
!--------------------------------------

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: phi, Al
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: hh
    real(kind=DP), intent(in) :: time

    integer, parameter :: rescale_max_num = 10
    real(kind=DP) :: ff_max_l, ff_max, wr1
    real(kind=DP), dimension(-nx:nx,0:ny) :: phi_max_l, phi_max, scale_factor
    integer :: mx, my, mxw
    integer, save :: iflg
    data iflg / 0 /


    ff_max_l = maxval(abs(ff))
    call MPI_Allreduce(ff_max_l, ff_max, 1, MPI_DOUBLE_PRECISION, &
                       MPI_MAX, MPI_COMM_WORLD, ierr_mpi )

    if (ff_max > 1.d50) then ! If ff is too large, rescale ff, phi, Al, hh

      write( olog, * ) " # Rescale linear run at time = ", time
      call tips_flush

      !== Determine rescale factor for each linear mode ==
      phi_max_l(:,:) = 0._DP
      do my = ist_y, iend_y
        do mx = -nx, nx
          phi_max_l(mx,my) = maxval(abs(phi(mx,my,:)))
        end do
      end do
      call MPI_Allreduce(phi_max_l, phi_max, (2*nx+1)*(ny+1), MPI_DOUBLE_PRECISION, &
                         MPI_MAX, zsp_comm_world, ierr_mpi )

      scale_factor(:,:) = 1._DP
      do my = ist_y, iend_y

        if (dj(my) == 0) then  ! dj==0, no mode connection in mx

          do mx = -nx, nx
            if (phi_max(mx,my) .ne. 0.d0) then
              scale_factor(mx,my) = 1.d-10 / phi_max(mx,my)
            end if
          end do

        else                   ! dj/=0, mode connection in mxw=mx+-dj

          do mx = -nx, -nx+abs(dj(my))-1
            wr1 = 0.d0
            do mxw = mx, nx, abs(dj(my))
              if (wr1 < phi_max(mxw,my)) wr1 = phi_max(mxw,my)
            end do
            if (wr1 .ne. 0.d0) then
              do mxw = mx, nx, abs(dj(my))
                scale_factor(mxw,my) = 1.d-10 / wr1
              end do
            end if
          end do

        end if

      end do

      !== Rescale ==
      do my = ist_y, iend_y
        do mx = -nx, nx
          ff(mx,my,:,:,:) = scale_factor(mx,my) * ff(mx,my,:,:,:)
          hh(mx,my,:,:,:) = scale_factor(mx,my) * hh(mx,my,:,:,:)
          phi(mx,my,:) = scale_factor(mx,my) * phi(mx,my,:)
          Al(mx,my,:) = scale_factor(mx,my) * Al(mx,my,:)
        end do
      end do

      iflg = iflg + 1

    end if

    if (iflg > rescale_max_num) then
      write( olog, * ) " # Reach max number of rescale:", iflg
      call tips_flush
      call MPI_Finalize(ierr_mpi)
      stop
    end if

  END SUBROUTINE tips_rescale_for_linear_runs


END MODULE GKV_tips
