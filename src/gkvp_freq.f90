MODULE GKV_freq
!-------------------------------------------------------------------------------
!
!    Module for evaluating linear growth rate and real frequency
!    (without shearflows)
!
!    Update history of gkvp_freq.f90
!    --------------
!      gkvp_f0.57 (S. Maeyama, Oct 2020)
!        - Version number f0.57 is removed from filename.
!        - nxfrq=min(5,nx,nx0) to avoid zero division when nx0=0.
!
!-------------------------------------------------------------------------------

  use GKV_header
  use GKV_mpienv

  implicit none

  private

  real(kind=DP), parameter :: eps_omega = 1.d-2, &
                              eps_gamma = 1.d-2, &
                              eps_ineq  = 1.d-2

  integer, save :: nxfrq = nx
  real(kind=DP), save :: time0 = 0._DP
  complex(kind=DP), allocatable, save, dimension(:,:,:) :: phi0
  complex(kind=DP), allocatable, save, dimension(:,:) :: omega0
  real(kind=DP), allocatable, save, dimension(:,:) :: phi0_norm2

  complex(kind=DP), allocatable, save, dimension(:,:) :: omega_g
  complex(kind=DP), allocatable, save, dimension(:,:) :: diff_g
  real(kind=DP), allocatable, save, dimension(:,:) :: ineq_g
  logical, allocatable, save, dimension(:,:) :: freq_conv

  public :: freq_set, freq_write_frq, freq_write_dsp, freq_conv


CONTAINS


!--------------------------------------
  SUBROUTINE freq_set ( time )
!--------------------------------------

    real(kind=DP), intent(in)           :: time

    integer :: mx, my


      nxfrq = min( 5, nx, nx0 )
      time0 = time
  
      allocate( phi0(-nxfrq:nxfrq,0:ny,-nz:nz-1) )
      allocate( omega0(-nxfrq:nxfrq,0:ny) )
      allocate( phi0_norm2(-nxfrq:nxfrq,0:ny) )
      allocate( omega_g(-nxfrq:nxfrq,0:(ny+1)*nprocw-1) )
      allocate( diff_g(-nxfrq:nxfrq,0:(ny+1)*nprocw-1) )
      allocate( ineq_g(-nxfrq:nxfrq,0:(ny+1)*nprocw-1) )
      allocate( freq_conv(-nxfrq:nxfrq,0:(ny+1)*nprocw-1) )
  
      phi0(:,:,:) = ( 0._DP, 0._DP )
      omega0(:,:) = ( 0._DP, 0._DP )
      phi0_norm2(:,:) = 0._DP
      omega_g(:,:) = ( 0._DP, 0._DP )
      diff_g(:,:) = ( 0._DP, 0._DP )
      ineq_g(:,:) = 0._DP
      freq_conv(:,:) = .false.
      
!- write hst/*.frq.* -
      if ( rankg == 0 ) then
        write(ofrq, fmt="(a)") "#  Re[omega], Im[omega] for (kx,ky)"
        write(ofrq, fmt="(99a17)", advance="no") "#            time"
        do mx = 0, nxfrq
          do my = 1, global_ny
            write(ofrq, "(a,f7.3,a,f7.3,a)", advance="no")  &
                "(", kx(mx), ",", ky(1) * real(my, kind=DP), ")"
            write(ofrq, "(a17)", advance="no") " "
          end do
        end do
        if (nxfrq > 0) then
          do mx = -nxfrq, -1
            do my = 1, global_ny
              write(ofrq, "(a,f7.3,a,f7.3,a)", advance="no")  &
                  "(", kx(mx), ",", ky(1) * real(my, kind=DP), ")"
              write(ofrq, "(a17)", advance="no") " "
            end do
          end do
        end if
        write(ofrq, *)
      end if

  END SUBROUTINE freq_set


!--------------------------------------
  SUBROUTINE freq_reset
!--------------------------------------


      deallocate( phi0 )
      deallocate( omega0 )
      deallocate( phi0_norm2 )
      deallocate( omega_g )
      deallocate( diff_g )
      deallocate( ineq_g )
      deallocate( freq_conv )
                                      !%%% For debug %%%
                                      !write(olog,*) "freq_reset"
                                      !%%%%%%%%%%%%%%%%%

  END SUBROUTINE freq_reset


!--------------------------------------
  SUBROUTINE freq_write_frq ( time, phi )
!--------------------------------------

    real(kind=DP), intent(in)           :: time
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1)   :: phi

    complex(kind=DP), dimension(-nxfrq:nxfrq,0:ny) :: phi0phi, omega_l, diff_l
    real(kind=DP), dimension(-nxfrq:nxfrq,0:ny) :: phi_norm2, ineq_l
    complex(kind=DP), dimension(-nxfrq:nxfrq,0:ny) :: wc3
    real(kind=DP), dimension(-nxfrq:nxfrq,0:ny) :: wr3
    integer :: mx, my, iz

    integer, save ::  iflg
    data iflg / 0 /


      if( iflg == 0 ) then
        iflg = 1
        phi0(-nxfrq:nxfrq,:,:) = phi(-nxfrq:nxfrq,:,:)
        phi0_norm2(:,:) = 0._DP
        wr3(:,:) = 0._DP
        do iz = -nz, nz-1
          do my = ist1_y, iend_y
            do mx = -nxfrq, nxfrq
              wr3(mx,my) = wr3(mx,my) + abs( phi(mx,my,iz) )**2
            end do
          end do
        end do
        call MPI_Allreduce( wr3, phi0_norm2, (2*nxfrq+1)*(ny+1), MPI_DOUBLE_PRECISION, &
                            MPI_SUM, zsp_comm_world, ierr_mpi )
        return
      end if


!- calculate interior products -
      phi0phi(:,:) = (0._DP, 0._DP)
      wc3(:,:) = (0._DP, 0._DP)
      do iz = -nz, nz-1
        do my = ist1_y, iend_y
          do mx = -nxfrq, nxfrq
            wc3(mx,my) = wc3(mx,my) + conjg( phi0(mx,my,iz) ) * phi(mx,my,iz)
          end do
        end do
      end do
      call MPI_Allreduce( wc3, phi0phi, (2*nxfrq+1)*(ny+1), MPI_DOUBLE_COMPLEX, &
                          MPI_SUM, zsp_comm_world, ierr_mpi )

      phi_norm2(:,:) = 0._DP
      wr3(:,:) = 0._DP
      do iz = -nz, nz-1
        do my = ist1_y, iend_y
          do mx = -nxfrq, nxfrq
            wr3(mx,my) = wr3(mx,my) + abs( phi(mx,my,iz) )**2
          end do
        end do
      end do
      call MPI_Allreduce( wr3, phi_norm2, (2*nxfrq+1)*(ny+1), MPI_DOUBLE_PRECISION, &
                          MPI_SUM, zsp_comm_world, ierr_mpi )

!- calculate frequency -
      omega_l(:,:) = (0._DP, 0._DP)
      do my = ist1_y, iend_y
        do mx = -nxfrq, nxfrq
          omega_l(mx,my) = log( phi0phi(mx,my) / phi0_norm2(mx,my) )  &
                                      / ( ui * ( time0 - time ) )
        end do
      end do
      call MPI_Allgather( omega_l, (2*nxfrq+1)*(ny+1), MPI_DOUBLE_COMPLEX, &
                          omega_g, (2*nxfrq+1)*(ny+1), MPI_DOUBLE_COMPLEX, &
                          fft_comm_world, ierr_mpi )

!- convergence check -
      diff_l(:,:) = (0._DP, 0._DP)
      do my = ist1_y, iend_y
        do mx = -nxfrq, nxfrq
          diff_l(mx,my) = abs(real(omega_l(mx,my) - omega0(mx,my), kind=DP)  &
                                      / real(omega_l(mx,my), kind=DP)) / (time - time0) &
                        + ui * abs(aimag(omega_l(mx,my) - omega0(mx,my))     &
                                      / aimag(omega_l(mx,my)) ) / (time - time0)
        end do
      end do
      call MPI_Allgather( diff_l, (2*nxfrq+1)*(ny+1), MPI_DOUBLE_COMPLEX, &
                          diff_g, (2*nxfrq+1)*(ny+1), MPI_DOUBLE_COMPLEX, &
                          fft_comm_world, ierr_mpi )

      ineq_l(:,:) = 0._DP
      do my = ist1_y, iend_y
        do mx = -nxfrq, nxfrq
          ineq_l(mx,my) = sqrt(abs(phi0phi(mx,my))**2 / (phi0_norm2(mx,my) * phi_norm2(mx,my)))
        end do
      end do
      call MPI_Allgather( ineq_l, (2*nxfrq+1)*(ny+1), MPI_DOUBLE_PRECISION, &
                          ineq_g, (2*nxfrq+1)*(ny+1), MPI_DOUBLE_PRECISION, &
                          fft_comm_world, ierr_mpi )

      do my = 1, global_ny
        do mx = -nxfrq, nxfrq
          if ( real( diff_g(mx,my), kind=DP ) < eps_omega .and.  &
               aimag( diff_g(mx,my) ) < eps_gamma .and.          &
               abs(1._DP - ineq_g(mx,my)) < eps_ineq ) then
            freq_conv(mx,my) = .true.
          else
            freq_conv(mx,my) = .false.
          end if
        end do
      end do

      freq_conv(:,0) = .true.
      if ( global_ny < (ny+1)*nprocw-1 ) then
        do my = global_ny+1, (ny+1)*nprocw-1
          diff_g(:,my) = ( 0._DP, 0._DP )
          omega_g(:,my) = ( 0._DP, 0._DP )
          ineq_g(:,my) = 0._DP
          freq_conv(:,my) = .true.
        end do
      end if

!- remember the values -
      time0 = time
      phi0(-nxfrq:nxfrq,:,:) = phi(-nxfrq:nxfrq,:,:)
      omega0(:,:) = omega_l(:,:)
      phi0_norm2(:,:) = phi_norm2(:,:)

!- write hst/*.frq.* -
      if ( rankg == 0 ) then
        !write( ofrq, '(9999G17.7e3)' ) time, (omega_g(0,my), my=1,global_ny)
         write( ofrq, '(9999G17.7e3)' ) time, ((omega_g(mx,my), my=1,global_ny), mx=0,nxfrq)
      end if
                                      !%%% For debug %%%
                                      !write(olog,*) "freq_write_frq at t=",time
                                      !%%%%%%%%%%%%%%%%%

END SUBROUTINE freq_write_frq


!--------------------------------------
  SUBROUTINE freq_write_dsp
!--------------------------------------

    integer :: mx, my


      if ( rankg == 0 ) then
        write( odsp, '(99A17)' ) "#              kx","ky","frequency","growthrate",&
                                           "diff(freq)","diff(grow)","1-ineq"
        do mx = -nxfrq, nxfrq
          do my = 1, global_ny
            if ( freq_conv(mx,my) ) then
              write( odsp, '(9999G17.7e3)' ) kx(mx), ky(1) * real( my, kind=DP ),       &
                         real( omega_g(mx,my), kind=DP ), aimag( omega_g(mx,my) ),    &
                         real( diff_g(mx,my), kind=DP ), aimag( diff_g(mx,my) ),      &
                         abs(1._DP - ineq_g(mx,my))
            else
              write( odsp, '(A2,9999G17.7e3)' ) "# ", kx(mx), ky(1) * real( my, kind=DP ),&
                         real( omega_g(mx,my), kind=DP ), aimag( omega_g(mx,my) ),      &
                         real( diff_g(mx,my), kind=DP ), aimag( diff_g(mx,my) ),        &
                         abs(1._DP - ineq_g(mx,my))
            end if
          end do
          write( odsp, * )
        end do

      end if
  
      call freq_reset
                                      !%%% For debug %%%
                                      !write(olog,*) "freq_write_dsp"
                                      !%%%%%%%%%%%%%%%%%


  END SUBROUTINE freq_write_dsp


END MODULE GKV_freq
