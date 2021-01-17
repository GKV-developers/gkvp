MODULE GKV_fft
!-------------------------------------------------------------------------------
!
!    FFT module for E x B term calculation using FFTW
!
!    Update history of gkvp_fft_fftw.f90
!    --------------
!      gkvp_f0.57 (S. Maeyama, Oct 2020)
!        - Version number f0.57 is removed from filename.
!
!-------------------------------------------------------------------------------

  use GKV_header
  use GKV_mpienv
  use GKV_clock, only: clock_sta, clock_end

  implicit none

  include "fftw3.f"

  private

  integer(kind=DP), save      :: plan_x_forward, plan_x_backward
  integer(kind=DP), save      :: plan_y_forward, plan_y_backward

  public   fft_pre,  &
           !fft_backward_Xfft, fft_backward_chXY, fft_backward_Yfft, &
           !fft_forward_Yfft, fft_forward_chYX, fft_forward_Xfft,    &
           plan_x_forward, plan_x_backward,                         &
           plan_y_forward, plan_y_backward


CONTAINS


!--------------------------------------
  SUBROUTINE fft_pre( )
!--------------------------------------
!  Initialization of FFT

    complex(kind=DP) :: wk1_x_z(0:2*nxw-1)
    complex(kind=DP) :: wk2_x_z(0:2*nxw-1)
    complex(kind=DP) :: wk1_y_z(0:nyw)
    real(kind=DP)    :: wk2_y_r(0:2*nyw-1)


     wk1_x_z(:) = (0._DP, 0._DP)
     wk2_x_z(:) = (0._DP, 0._DP)
     wk1_y_z(:) = (0._DP, 0._DP)
     wk2_y_r(:) = 0._DP

     call dfftw_plan_dft_1d( plan_x_backward,     &
                             (2*nxw),             &
                             wk1_x_z,             &   ! complex in
                             wk2_x_z,             &   ! complex out
                             FFTW_BACKWARD,       &
                             FFTW_MEASURE )
                             !! FFTW_ESTIMATE )

     call dfftw_plan_dft_c2r_1d( plan_y_backward, &
                                 (2*nyw),         &
                                 wk1_y_z,         &   ! complex in
                                 wk2_y_r,         &   ! real    out
                                 FFTW_MEASURE )
                                 !! FFTW_ESTIMATE )

     call dfftw_plan_dft_r2c_1d( plan_y_forward,  &
                                 (2*nyw),         &
                                 wk2_y_r,         &   ! real    in
                                 wk1_y_z,         &   ! complex out
                                 FFTW_MEASURE )
                                 !! FFTW_ESTIMATE )

     call dfftw_plan_dft_1d( plan_x_forward,      &
                             (2*nxw),             &
                             wk2_x_z,             &   ! complex in
                             wk1_x_z,             &   ! complex out
                             FFTW_FORWARD,        &
                             FFTW_MEASURE )
                             !! FFTW_ESTIMATE )


  END SUBROUTINE fft_pre


!!--------------------------------------
!  SUBROUTINE fft_backward_Xfft ( exbdf, send_buff, num_trans )
!!--------------------------------------
!
!    complex(kind=DP), intent(in), &
!      dimension(0:2*nxw-1,0:ny,num_trans)   :: exbdf
!    complex(kind=DP), intent(out), &
!      dimension(0:ny, 0:nxw_size, num_trans, 0:nprocw-1) :: send_buff
!    integer, intent(in)  :: num_trans
!
!    complex(kind=DP), dimension(0:2*nxw-1) :: wk_x_out ! for outplace
!    integer :: ist_xw_g_rank
!    integer :: mx, my, i, irank
!
!
!! ---- Backward FFT in X ----------------------------------
!!$OMP master
!                                           call clock_sta(1421)
!                                         ! call fapp_start("nlterm_backward_Xfft",1421,1)
!!$OMP end master
!!$OMP do schedule (dynamic)
!      do i = 1, num_trans
!        do my = ist_y, iend_y
!
!          call dfftw_execute_dft ( &
!                          plan_x_backward,  &
!                          exbdf(0,my,i),    &     ! complex in
!                        !  exbdf(0,my,i)        )  ! complex out
!                          wk_x_out             )  ! complex out ! for outplace
!
!        ! --- set send buffer ---
!          do irank = 0, nprocw-1
!            ist_xw_g_rank  = (nxw_size+1)*irank
!            do mx = ist_xw, iend_xw
!             !  send_buff(my,mx,i,irank) = exbdf(mx+ist_xw_g_rank,my,i)
!               send_buff(my,mx,i,irank) = wk_x_out(mx+ist_xw_g_rank) ! for outplace
!            enddo
!          enddo
!
!        enddo
!      end do
!!$OMP end do nowait
!!$OMP master
!                                         ! call fapp_stop("nlterm_backward_Xfft",1421,1)
!                                           call clock_end(1421)
!!$OMP end master
!
!
!  END SUBROUTINE fft_backward_Xfft
!
!
!!--------------------------------------
!  SUBROUTINE fft_backward_chXY ( send_buff, recv_buff, num_trans )
!!--------------------------------------
!
!    complex(kind=DP), intent(in), &
!      dimension(0:ny, 0:nxw_size, num_trans, 0:nprocw-1) :: send_buff
!    complex(kind=DP), intent(out), &
!      dimension(0:ny, 0:nxw_size, num_trans, 0:nprocw-1) :: recv_buff
!    integer, intent(in)  :: num_trans
!
!    integer :: nsize_com
!    integer :: irc
!
!
!! ---- Data Exchange ----------------------------------
!!         Y divide -> X divide
!      nsize_com = (ny+1)*(nxw_size+1)*num_trans
!
!                                           call clock_sta(1422)
!                                         ! call fapp_start("nlterm_backward_shiftXY",1422,1)
!      call mpi_alltoall( send_buff, &
!                         nsize_com,              &
!                         MPI_DOUBLE_COMPLEX,     &
!                         recv_buff,              &
!                         nsize_com,              &
!                         MPI_DOUBLE_COMPLEX,     &
!                         fft_comm_world,         &
!                         irc        )
!                                         ! call fapp_stop("nlterm_backward_shiftXY",1422,1)
!                                           call clock_end(1422)
!    
!
!  END SUBROUTINE fft_backward_chXY
!
!
!!--------------------------------------
!  SUBROUTINE fft_backward_Yfft ( recv_buff, exbdf_xw, num_trans )
!!--------------------------------------
!
!    complex(kind=DP), intent(in), &
!      dimension(0:ny, 0:nxw_size, num_trans, 0:nprocw-1) :: recv_buff
!    complex(kind=DP), intent(out), &
!      dimension(0:nyw,0:nxw_size,num_trans)   :: exbdf_xw
!    integer, intent(in)  :: num_trans
!
!    complex(kind=DP) :: wk_y_in(0:nyw) ! for outplace
!    integer :: ist_y_g_rank, iend_y_g_rank
!    integer :: mx, my, i, irank
!
!
!! ---- Backward FFT in Y ----------------------------------
!!$OMP master
!                                           call clock_sta(1423)
!                                         ! call fapp_start("nlterm_backward_Yfft",1423,1)
!!$OMP end master
!!$OMP do schedule (dynamic)
!      do i = 1, num_trans
!        do mx = ist_xw, iend_xw
!
!        ! --- restore receive buffer ---
!          do irank = 0, nprocw-1
!            ist_y_g_rank  = (ny+1)*irank
!            iend_y_g_rank = min ( (ny+1)*(irank+1)-1, global_ny )
!            do my = ist_y_g_rank, iend_y_g_rank
!            !  exbdf_xw(my,mx,i) = recv_buff(my-ist_y_g_rank,mx,i,irank)
!              wk_y_in(my) = recv_buff(my-ist_y_g_rank,mx,i,irank) ! for outplace
!            end do
!          end do
!        ! --- set filler ---
!          do my = global_ny+1, nyw
!          !  exbdf_xw(my,mx,i) = ( 0._DP, 0._DP )
!            wk_y_in(my) = ( 0._DP, 0._DP ) ! for outplace
!          end do
!
!          call dfftw_execute_dft_c2r( &
!                          plan_y_backward,   &
!                        !  exbdf_xw(0,mx,i),  &    ! complex in
!                          wk_y_in,            &    ! complex in ! for outplace
!                          exbdf_xw(0,mx,i)     )  ! real    out
!
!        end do
!      end do
!!$OMP end do nowait
!!$OMP master
!                                         ! call fapp_stop("nlterm_backward_Yfft",1423,1)
!                                           call clock_end(1423)
!!$OMP end master
!
!
!  END SUBROUTINE fft_backward_Yfft
!
!
!!--------------------------------------
!  SUBROUTINE fft_forward_Yfft ( exbdf_xw, send_buff, num_trans )
!!--------------------------------------
!
!    complex(kind=DP), intent(in), &
!      dimension(0:nyw,0:nxw_size,num_trans)   :: exbdf_xw
!    complex(kind=DP), intent(out), &
!      dimension(0:ny, 0:nxw_size, num_trans, 0:nprocw-1) :: send_buff
!    integer, intent(in)  :: num_trans
!
!    complex(kind=DP) :: wk_y_out(0:nyw) ! for outplace
!    integer :: ist_y_g_rank, iend_y_g_rank
!    integer :: mx, my, i, irank
!
!
!! ---- Forward FFT in Y ----------------------------------
!!$OMP master
!                                           call clock_sta(1441)
!                                         ! call fapp_start("nlterm_forward_Yfft",1441,1)
!!$OMP end master
!!$OMP do schedule (dynamic)
!      do i = 1, num_trans
!        do mx = ist_xw, iend_xw
!
!          call dfftw_execute_dft_r2c( &
!                          plan_y_forward,    &
!                          exbdf_xw(0,mx,i),  &    ! real in
!                        !  exbdf_xw(0,mx,i)     )  ! complex out
!                          wk_y_out             )  ! complex out ! for outplace
!
!        ! --- set send buffer ---
!          do irank = 0, nprocw-1
!            ist_y_g_rank  = (ny+1)*irank
!            iend_y_g_rank = min ( (ny+1)*(irank+1)-1, global_ny )
!            do my = ist_y_g_rank, iend_y_g_rank
!            !  send_buff(my-ist_y_g_rank,mx,i,irank) = exbdf_xw(my,mx,i)
!              send_buff(my-ist_y_g_rank,mx,i,irank) = wk_y_out(my) ! for outplace
!            end do
!          end do
!
!      end do
!    end do
!!$OMP end do nowait
!!$OMP master
!                                         ! call fapp_stop("nlterm_forward_Yfft",1441,1)
!                                           call clock_end(1441)
!!$OMP end master
!
!
!  END SUBROUTINE fft_forward_Yfft
!
!
!!--------------------------------------
!  SUBROUTINE fft_forward_chYX ( send_buff, recv_buff, num_trans )
!!--------------------------------------
!
!    complex(kind=DP), intent(in), &
!      dimension(0:ny, 0:nxw_size, num_trans, 0:nprocw-1) :: send_buff
!    complex(kind=DP), intent(out), &
!      dimension(0:ny, 0:nxw_size, num_trans, 0:nprocw-1) :: recv_buff
!    integer, intent(in)  :: num_trans
!
!    integer :: nsize_com
!    integer :: irc
!
!
!! ---- Data Exchange ----------------------------------
!!         X divide -> Y divide
!      nsize_com = (ny+1)*(nxw_size+1)*num_trans
!
!                                           call clock_sta(1442)
!                                         ! call fapp_start("nlterm_forward_shiftYX",1442,1)
!      call mpi_alltoall( send_buff, &
!                         nsize_com,              &
!                         MPI_DOUBLE_COMPLEX,     &
!                         recv_buff,              &
!                         nsize_com,              &
!                         MPI_DOUBLE_COMPLEX,     &
!                         fft_comm_world,         &
!                         irc        )
!                                         ! call fapp_stop("nlterm_forward_shiftYX",1442,1)
!                                           call clock_end(1442)
!
!
!  END SUBROUTINE fft_forward_chYX
!
!
!!--------------------------------------
!  SUBROUTINE fft_forward_Xfft ( recv_buff, exbdf, num_trans )
!!--------------------------------------
!
!    complex(kind=DP), intent(in), &
!      dimension(0:ny, 0:nxw_size, num_trans, 0:nprocw-1) :: recv_buff
!    complex(kind=DP), intent(out), &
!      dimension(0:2*nxw-1,0:ny,num_trans)   :: exbdf
!    integer, intent(in)  :: num_trans
!
!    complex(kind=DP), dimension(0:2*nxw-1) :: wk_x_in ! for outplace
!    integer :: ist_xw_g_rank, iend_xw_g_rank
!    integer :: mx, my, i, irank
!
!
!! ---- Forward FFT in X ----------------------------------
!!$OMP master
!                                           call clock_sta(1443)
!                                         ! call fapp_start("nlterm_forward_Xfft",1443,1)
!!$OMP end master
!!$OMP do schedule (dynamic)
!      do i = 1, num_trans
!        do my = ist_y, iend_y
!
!        ! ---  restore receive buffer ---
!          do irank = 0, nprocw-1
!            ist_xw_g_rank  = (nxw_size+1)*irank
!            iend_xw_g_rank = min ( (nxw_size+1)*(irank+1)-1, (2*nxw-1) )
!            do mx = ist_xw_g_rank, iend_xw_g_rank
!            !  exbdf(mx,my,i) = recv_buff(my,mx-ist_xw_g_rank,i,irank)
!              wk_x_in(mx) = recv_buff(my,mx-ist_xw_g_rank,i,irank) ! for outplace
!            enddo
!          enddo
!
!          call dfftw_execute_dft ( &
!                          plan_x_forward,   &
!                        !  exbdf(0,my,i),    &    ! complex in
!                          wk_x_in,          &    ! complex in ! for outplace
!                          exbdf(0,my,i)       )  ! complex out
!
!        end do
!      end do
!!$OMP end do nowait
!!$OMP master
!                                         ! call fapp_stop("nlterm_forward_Xfft",1443,1)
!                                           call clock_end(1443)
!!$OMP end master
!
!
!  END SUBROUTINE fft_forward_Xfft


END MODULE GKV_fft
