MODULE GKV_exb
!-------------------------------------------------------------------------------
!
!    E x B term
!
!      GKV-plus r0.3 ( T.-H.Watanabe, Jun 2011)
!
!-------------------------------------------------------------------------------

  use GKV_header
  use GKV_mpienv
  use GKV_fft, only: &
           plan_xf_y2zm, plan_xf_y2x, &
           plan_xb_y2zm, plan_xb_y2x, &
           plan_yf_y2zm, plan_yf_y2x, &
           plan_yb_y2zm, plan_yb_y2x, &
           planr_xf_y2zm, planr_xf_y2x, &
           planr_xb_y2zm, planr_xb_y2x, &
           planr_yf_y2zm, planr_yf_y2x, &
           planr_yb_y2zm, planr_yb_y2x, &
           plan_x_forward, plan_x_backward, &
           plan_y_forward, plan_y_backward
  use GKV_clock, only: clock_sta, clock_end

  implicit none

  private

  real(kind=DP), save :: exb_maxvx_eachrank, exb_maxvy_eachrank

  integer, parameter :: nbuff = ((2*nz)*(nm+1)-1)/nprocw + 1
                            !if ( mod(2*nz*(nm+1),nprocw) == 0 ) then
                            !  nbuff = 2*nz*(nm+1)/nprocw
                            !else
                            !  nbuff = 2*nz*(nm+1)/nprocw + 1
                            !end if
  real(kind=DP), dimension(0:global_ny), save :: gky
  integer, save :: nchunk_zm = 1, nchunk_yb = 1, nchunk_xb = 1

  integer, save :: nchunk_yzm = 1, nchunk_xzm = 1

  complex(kind=DP), save, dimension(0:ny,0:2*nxw-1) :: uikx_y2x, uiky_y2x


  public   exb_NL_term, exb_maxvx_eachrank, exb_maxvy_eachrank


CONTAINS


!--------------------------------------
  SUBROUTINE exb_NL_term(hh, psi, chi, pb)
!--------------------------------------
!  Nonlinear term calculation interface
    implicit none
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: hh
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,0:nm) :: psi, chi
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: pb

    real(kind=DP) :: dky
    integer, save :: iflg
    integer :: my
!$  integer :: nthreads, omp_get_num_threads
    data iflg / 0 /
                                               !%%% For debug %%%
                                               !complex(kind=DP) ::                        &
                                               !  whh(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm),   &
                                               !  wpsi(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,0:nm), &
                                               !  wchi(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,0:nm)
                                               !integer :: mx, iz, iv, im
                                               !%%%%%%%%%%%%%%%%%

    if( iflg == 0 ) then
      iflg = 1
      dky = ky(1) - ky(0)
      do my = 0, global_ny
        gky(my) = dky * real(my, kind=DP)
      end do
      exb_maxvx_eachrank = eps
      exb_maxvy_eachrank = eps
!$OMP parallel default(shared)
!$OMP master
!$    nthreads = omp_get_num_threads()
!$    if (nthreads > 1) then
!$      nchunk_zm = ((2*nz)*(nm+1)-1) / (nthreads-1) + 1
!$      nchunk_yb = ((global_ny+1)*nbuff-1) / (nthreads-1) + 1
!$      nchunk_xb = ((2*nxw)*nbuff-1) / (nthreads-1) + 1
!$      nchunk_yzm = ((iend_y-ist_y+1)*(2*nz)*(nm+1)-1) / (nthreads-1) + 1
!$      nchunk_xzm = ((iend_xw-ist_xw+1)*(2*nz)*(nm+1)-1) / (nthreads-1) + 1
!$    end if
!$OMP end master
!$OMP end parallel
    end if
                                               !%%% For debug %%%
                                               !whh(:,:,:,:,:) = (0._DP, 0._DP)
                                               !wpsi(:,:,:,:) = (0._DP, 0._DP)
                                               !wchi(:,:,:,:) = (0._DP, 0._DP)
                                               !if (rankw == 0) then
                                               !  whh(0,1,:,:,:) = (0.5_DP, 0._DP)
                                               !  wpsi(2,0,:,:) = (0._DP, 0.5_DP)
                                               !  wpsi(-2,0,:,:) = (0._DP, -0.5_DP)
                                               !end if
                                               !%%%%%%%%%%%%%%%%%

    if (trim(calc_type) == "nonlinear") then

#ifdef USE_TERM_Y2ZM
        call exb_NL_term_y2zm(hh, psi, chi, pb)
#else
        call exb_NL_term_y2x(hh, psi, chi, pb)
#endif

    else

!$OMP parallel workshare
      pb(:,:,:,:,:) = ( 0._DP, 0._DP )
!$OMP end parallel workshare

    end if

                                               !%%% For debug %%%
                                               !if (rankz == 0 .and. rankv == 0 .and. rankm == 0 .and. ranks == 0) then
                                               !  im = 0; iv = 1; iz = 0
                                               !  do my = 0, ny
                                               !    do mx = -nx, nx
                                               !      write(80000+rankg,*) kx(mx), ky(my),  &
                                               !  dble(whh(mx,my,iz,iv,im)), aimag(whh(mx,my,iz,iv,im)), &
                                               !  dble(wpsi(mx,my,iz,im)), aimag(wpsi(mx,my,iz,im)), &
                                               !  dble(pb(mx,my,iz,iv,im)), aimag(pb(mx,my,iz,iv,im))
                                               !    end do
                                               !    write(80000+rankg,*)
                                               !  end do
                                               !end if
                                               !call MPI_Finalize(ierr_mpi)
                                               !stop
                                               !%%%%%%%%%%%%%%%%%


  END SUBROUTINE exb_NL_term


!--------------------------------------
  SUBROUTINE exb_NL_term_y2zm( hh, psi, chi, ef )
!--------------------------------------
!  ExB nonlinear term calculation 

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: hh
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,0:nm) :: psi, chi
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: ef

    real(kind=DP), dimension(:,:,:), allocatable :: dpdx, dpdy, dadx, dady
    complex(kind=DP), dimension(:,:,:,:), allocatable :: &
                      wc1o, wc2o, wc3o, wc4o,   wc1e, wc2e, wc3e, wc4e
    complex(kind=DP), dimension(:,:,:), allocatable ::   &
                         wwdxo, wwdyo, wwefo,      wwdxe, wwdye, wwefe

    integer :: ix,iy,iz,im,ib,ip
    integer :: iv

      allocate(dpdx(0:nbuff-1,0:2*nxw-1,0:2*nyw-1))
      allocate(dpdy(0:nbuff-1,0:2*nxw-1,0:2*nyw-1))
      allocate(dadx(0:nbuff-1,0:2*nxw-1,0:2*nyw-1))
      allocate(dady(0:nbuff-1,0:2*nxw-1,0:2*nyw-1))
      allocate(wc1o(0:nbuff-1,-nx:nx,0:ny,0:nprocw-1))
      allocate(wc2o(0:nbuff-1,-nx:nx,0:ny,0:nprocw-1))
      allocate(wc3o(0:nbuff-1,-nx:nx,0:ny,0:nprocw-1))
      allocate(wc4o(0:nbuff-1,-nx:nx,0:ny,0:nprocw-1))
      allocate(wc1e(0:nbuff-1,-nx:nx,0:ny,0:nprocw-1))
      allocate(wc2e(0:nbuff-1,-nx:nx,0:ny,0:nprocw-1))
      allocate(wc3e(0:nbuff-1,-nx:nx,0:ny,0:nprocw-1))
      allocate(wc4e(0:nbuff-1,-nx:nx,0:ny,0:nprocw-1))
      allocate(wwdxo(0:nbuff-1,0:2*nxw-1,0:global_ny))
      allocate(wwdyo(0:nbuff-1,0:2*nxw-1,0:global_ny))
      allocate(wwefo(0:nbuff-1,0:2*nxw-1,0:global_ny))
      allocate(wwdxe(0:nbuff-1,0:2*nxw-1,0:global_ny))
      allocate(wwdye(0:nbuff-1,0:2*nxw-1,0:global_ny))
      allocate(wwefe(0:nbuff-1,0:2*nxw-1,0:global_ny))

!$OMP parallel default(none)                          &
!$OMP shared(hh,psi,chi,ef,dpdx,dpdy,dadx,dady)       &
!$OMP shared(wc1o,wc2o,wc3o,wc4o,wc1e,wc2e,wc3e,wc4e) &
!$OMP shared(wwdxo,wwdyo,wwefo,wwdxe,wwdye,wwefe)     &
!$OMP private(iv)

!$OMP workshare
      wc1o(:,:,:,:) = (0._DP, 0._DP)
      wc3o(:,:,:,:) = (0._DP, 0._DP)
      wc1e(:,:,:,:) = (0._DP, 0._DP)
      wc3e(:,:,:,:) = (0._DP, 0._DP)
!$OMP end workshare


!!%%% Without overlap %%%
!      call exb_pack_psi_y2zm(psi,wc3o)
!!$OMP barrier
!      call exb_transpose_y2zm(wc3o,wc4o)
!!$OMP barrier
!      call exb_unpack_y2zm(wc4o,wwdxo,wwdyo)
!!$OMP barrier
!      call exb_backwardfft_y2zm(wwdxo,wwdyo,dpdx,dpdy)
!!$OMP barrier
!      call exb_pack_psi_y2zm(chi,wc3e)
!!$OMP barrier
!      call exb_transpose_y2zm(wc3e,wc4e)
!!$OMP barrier
!      call exb_unpack_y2zm(wc4e,wwdxe,wwdye)
!!$OMP barrier
!      call exb_backwardfft_y2zm(wwdxe,wwdye,dadx,dady)
!!$OMP barrier
!      do iv = 1, 2*nv
!        call exb_pack_hh_y2zm(iv,hh,wc1o)
!!$OMP barrier
!        call exb_transpose_y2zm(wc1o,wc2o)
!!$OMP barrier
!        call exb_unpack_y2zm(wc2o,wwdxo,wwdyo)
!!$OMP barrier
!        call exb_realspcal_y2zm(iv,dpdx,dpdy,dadx,dady,wwdxo,wwdyo,wwefo)
!!$OMP barrier
!        call exb_pack_zm2y(wwefo,wc3o)
!!$OMP barrier
!        call exb_transpose_zm2y(wc3o,wc4o)
!!$OMP barrier
!        call exb_unpack_zm2y(iv,wc4o,ef)
!!$OMP barrier
!      end do
!!%%%%%%%%%%%%%%%%%%%%%%%


!%%% With overlap %%%
      call exb_pack_psi_y2zm(psi,wc3o)
!$OMP barrier
      call exb_transpose_y2zm(wc3o,wc4o)
      call exb_pack_psi_y2zm(chi,wc3e)
!$OMP barrier
      call exb_transpose_y2zm(wc3e,wc4e)
      call exb_unpack_y2zm(wc4o,wwdxo,wwdyo)
!$OMP barrier
      call exb_backwardfft_y2zm(wwdxo,wwdyo,dpdx,dpdy)
      call exb_unpack_y2zm(wc4e,wwdxe,wwdye)
!$OMP barrier
      call exb_backwardfft_y2zm(wwdxe,wwdye,dadx,dady)
      do iv = 1, 2*nv+6
        if (mod(iv,2) == 1) then ! odd
          if (1+1<=iv .and. iv<=2*nv+1) call exb_transpose_y2zm(wc1e,wc2e)
          if (1+5<=iv .and. iv<=2*nv+5) call exb_transpose_zm2y(wc3e,wc4e)
          if (1  <=iv .and. iv<=2*nv  ) call exb_pack_hh_y2zm(iv,hh,wc1o)
          if (1+2<=iv .and. iv<=2*nv+2) call exb_unpack_y2zm(wc2o,wwdxo,wwdyo)
          if (1+3<=iv .and. iv<=2*nv+3) call exb_realspcal_y2zm(iv-3,dpdx,dpdy,dadx,dady,wwdxe,wwdye,wwefe)
          if (1+4<=iv .and. iv<=2*nv+4) call exb_pack_zm2y(wwefo,wc3o)
          if (1+6<=iv .and. iv<=2*nv+6) call exb_unpack_zm2y(iv-6,wc4o,ef)
        else                     ! even
          if (1+1<=iv .and. iv<=2*nv+1) call exb_transpose_y2zm(wc1o,wc2o)
          if (1+5<=iv .and. iv<=2*nv+5) call exb_transpose_zm2y(wc3o,wc4o)
          if (1  <=iv .and. iv<=2*nv  ) call exb_pack_hh_y2zm(iv,hh,wc1e)
          if (1+2<=iv .and. iv<=2*nv+2) call exb_unpack_y2zm(wc2e,wwdxe,wwdye)
          if (1+3<=iv .and. iv<=2*nv+3) call exb_realspcal_y2zm(iv-3,dpdx,dpdy,dadx,dady,wwdxo,wwdyo,wwefo)
          if (1+4<=iv .and. iv<=2*nv+4) call exb_pack_zm2y(wwefe,wc3e)
          if (1+6<=iv .and. iv<=2*nv+6) call exb_unpack_zm2y(iv-6,wc4e,ef)
        end if
!$OMP barrier
      end do
!%%%%%%%%%%%%%%%%%%%%

!$OMP end parallel

      call exb_estimate_maxvel_y2zm(dpdx,dpdy,dadx,dady)

      deallocate(dpdx)
      deallocate(dpdy)
      deallocate(dadx)
      deallocate(dady)
      deallocate(wc1o)
      deallocate(wc2o)
      deallocate(wc3o)
      deallocate(wc4o)
      deallocate(wc1e)
      deallocate(wc2e)
      deallocate(wc3e)
      deallocate(wc4e)
      deallocate(wwdxo)
      deallocate(wwdyo)
      deallocate(wwefo)
      deallocate(wwdxe)
      deallocate(wwdye)
      deallocate(wwefe)

  END SUBROUTINE exb_NL_term_y2zm


!!--------------------------------------
!  SUBROUTINE exb_pack_y2zm ( psi, wc4 )
!!--------------------------------------
!!     Data pack for E x B term calculation (y2zm)
!
!    complex(kind=DP), intent(in), &
!      dimension(-nx:nx,0:ny,-nz:nz-1,0:nm) :: psi
!    complex(kind=DP), intent(inout), &
!      dimension(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1) :: wc4
!
!    integer :: mx, my, iz, im, izm, ibuff, iprocw
!
!!$OMP master
!                                           call clock_sta(1410)
!                                         ! call fapp_start("nlterm_pack",1410,1)
!!$OMP end master
!!$OMP do collapse(2) schedule(dynamic)
!      do im = 0, nm
!        do iz = -nz, nz-1
!
!         !%%% PACK: (kx,ky*,z*,m*)->(kx,ky,(z*,m*)*) %%%
!          izm = (2*nz)*im + (iz + nz)
!          ibuff = mod(izm, nbuff)
!          iprocw = izm / nbuff
!          do my = ist_y, iend_y
!            do mx = -nx, nx
!              wc4(mx,my,ibuff,iprocw) = psi(mx,my,iz,im)
!            end do
!          end do
!         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!        end do
!      end do
!!$OMP end do nowait
!!$OMP master
!                                         ! call fapp_stop("nlterm_pack",1410,1)
!                                           call clock_end(1410)
!!$OMP end master
!
!  END SUBROUTINE exb_pack_y2zm


!--------------------------------------
  SUBROUTINE exb_pack_psi_y2zm ( psi, wc4 )
!--------------------------------------
!     Data pack for E x B term calculation (y2zm)

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,0:nm) :: psi
    complex(kind=DP), intent(inout), &
      dimension(0:nbuff-1,-nx:nx,0:ny,0:nprocw-1) :: wc4

    integer :: mx, my, iz, im, izm, ibuff, iprocw

!$OMP master
                                           call clock_sta(1410)
                                         ! call fapp_start("nlterm_pack",1410,1)
!$OMP end master
!!TBI!! !$OMP do collapse(2) schedule(dynamic,nchunk_zm)
      do im = 0, nm
        do iz = -nz, nz-1

         !%%% PACK: (kx,ky*,z*,m*)->(kx,ky,(z*,m*)*) %%%
          izm = (2*nz)*im + (iz + nz)
          ibuff = mod(izm, nbuff)
          iprocw = izm / nbuff
          do my = ist_y, iend_y
            do mx = -nx, nx
              wc4(ibuff,mx,my,iprocw) = psi(mx,my,iz,im)
            end do
          end do
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end do
      end do
!!TBI!! !$OMP end do nowait
!$OMP master
                                         ! call fapp_stop("nlterm_pack",1410,1)
                                           call clock_end(1410)
!$OMP end master

  END SUBROUTINE exb_pack_psi_y2zm


!--------------------------------------
  SUBROUTINE exb_pack_hh_y2zm ( iv, hh, wc4 )
!--------------------------------------
!     Data pack for E x B term calculation (y2zm)

    integer, intent(in) :: iv
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: hh
    complex(kind=DP), intent(inout), &
      dimension(0:nbuff-1,-nx:nx,0:ny,0:nprocw-1) :: wc4

    integer :: mx, my, iz, im, izm, ibuff, iprocw

!$OMP master
                                           call clock_sta(1410)
                                         ! call fapp_start("nlterm_pack",1410,1)
!$OMP end master
!!TBI!! !$OMP do collapse(2) schedule(dynamic,nchunk_zm)
      do im = 0, nm
        do iz = -nz, nz-1

         !%%% PACK: (kx,ky*,z*,m*)->(kx,ky,(z*,m*)*) %%%
          izm = (2*nz)*im + (iz + nz)
          ibuff = mod(izm, nbuff)
          iprocw = izm / nbuff
          do my = ist_y, iend_y
            do mx = -nx, nx
              wc4(ibuff,mx,my,iprocw) = hh(mx,my,iz,iv,im)
            end do
          end do
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end do
      end do
!!TBI!! !$OMP end do nowait
!$OMP master
                                         ! call fapp_stop("nlterm_pack",1410,1)
                                           call clock_end(1410)
!$OMP end master

  END SUBROUTINE exb_pack_hh_y2zm


!--------------------------------------
  SUBROUTINE exb_transpose_y2zm ( wc4in, wc4out )
!--------------------------------------
!     Data transpose for E x B term calculation (y2zm)

    complex(kind=DP), intent(in), &
      dimension(0:nbuff-1,-nx:nx,0:ny,0:nprocw-1) :: wc4in
    complex(kind=DP), intent(out), &
      dimension(0:nbuff-1,-nx:nx,0:ny,0:nprocw-1) :: wc4out

!$OMP master
                                           call clock_sta(1420)
                                         ! call fapp_start("nlterm_alltoall1",1420,1)
      call MPI_Alltoall( wc4in,                 &
                         (2*nx+1)*(ny+1)*nbuff, &
                         MPI_DOUBLE_COMPLEX,    &
                         wc4out,                &
                         (2*nx+1)*(ny+1)*nbuff, &
                         MPI_DOUBLE_COMPLEX,    &
                         fft_comm_world,        &
                         ierr_mpi )
                                         ! call fapp_stop("nlterm_alltoall1",1420,1)
                                           call clock_end(1420)
!$OMP end master

  END SUBROUTINE exb_transpose_y2zm


!--------------------------------------
  SUBROUTINE exb_unpack_y2zm ( wc4, wwdx, wwdy )
!--------------------------------------
!     Data unpack for E x B term calculation (y2zm)

    complex(kind=DP), intent(in), &
      dimension(0:nbuff-1,-nx:nx,0:ny,0:nprocw-1) :: wc4
    complex(kind=DP), intent(out), &
      dimension(0:nbuff-1,0:2*nxw-1,0:global_ny) :: wwdx, wwdy

    complex(kind=DP), dimension(-nx:nx) :: psi
    complex(kind=DP), dimension(0:nbuff-1,0:global_ny,0:2*nxw-1) :: w1, w2
    integer :: mx, my, ibuff, iprocw, global_my

    integer :: ithd
!$  integer :: omp_get_thread_num

     ithd = 0
#ifndef OMP_INSIDE_FFTW
!$   ithd = omp_get_thread_num()
#endif

!$OMP master
                                           call clock_sta(1430)
                                         ! call fapp_start("nlterm_realspcal",1430,1)
!$OMP end master
!!TBI!! !$OMP do schedule(dynamic)
      do global_my = 0, global_ny
        iprocw = global_my / (ny+1)
        my = mod(global_my, ny+1)
        do mx = 0, nx
          do ibuff = 0, nbuff-1
            w1(ibuff,global_my,mx) = ui * kx(mx) * wc4(ibuff,mx,my,iprocw)
          end do
        end do
        do mx = nx+1, 2*nxw-nx-1
          do ibuff = 0, nbuff-1
            w1(ibuff,global_my,mx) = (0._DP, 0._DP) ! FFTW may destroy input array!
          end do
        end do
        do mx = -nx,-1
          do ibuff = 0, nbuff-1
            w1(ibuff,global_my,2*nxw+mx) = ui * kx(mx) * wc4(ibuff,mx,my,iprocw)
          end do
        end do
      end do

      call dfftw_execute_dft(planr_xb_y2zm(ithd), w1, w2)

      do global_my = 0, global_ny
        do mx = 0, 2*nxw-1
          do ibuff = 0, nbuff-1
            wwdx(ibuff,mx,global_my) = w2(ibuff,global_my,mx)
          end do
        end do
      end do

      do global_my = 0, global_ny
        iprocw = global_my / (ny+1)
        my = mod(global_my, ny+1)
        do mx = 0, nx
          do ibuff = 0, nbuff-1
            w1(ibuff,global_my,mx) = ui * gky(global_my) * wc4(ibuff,mx,my,iprocw)
          end do
        end do
        do mx = nx+1, 2*nxw-nx-1
          do ibuff = 0, nbuff-1
            w1(ibuff,global_my,mx) = (0._DP, 0._DP) ! FFTW may destroy input array!
          end do
        end do
        do mx = -nx,-1
          do ibuff = 0, nbuff-1
            w1(ibuff,global_my,2*nxw+mx) = ui * gky(global_my) * wc4(ibuff,mx,my,iprocw)
          end do
        end do
      end do

      call dfftw_execute_dft(planr_xb_y2zm(ithd), w1, w2)

      do global_my = 0, global_ny
        do mx = 0, 2*nxw-1
          do ibuff = 0, nbuff-1
            wwdy(ibuff,mx,global_my) = w2(ibuff,global_my,mx)
          end do
        end do
      end do
!!TBI!! !$OMP end do nowait
!$OMP master
                                         ! call fapp_stop("nlterm_realspcal",1430,1)
                                           call clock_end(1430)
!$OMP end master

  END SUBROUTINE exb_unpack_y2zm


!--------------------------------------
  SUBROUTINE exb_backwardfft_y2zm ( wwdx, wwdy, dpdx, dpdy )
!--------------------------------------
!     Backward FFT of field(psi,chi) for E x B term calculation (y2zm)

    complex(kind=DP), intent(in), &
      dimension(0:nbuff-1,0:2*nxw-1,0:global_ny) :: wwdx, wwdy
    real(kind=DP), intent(out), &
      dimension(0:nbuff-1,0:2*nxw-1,0:2*nyw-1) :: dpdx, dpdy

    complex(kind=DP), dimension(0:nbuff-1,0:2*nxw-1,0:nyw) :: w3
    integer :: ix, ibuff

    integer :: ithd
!$  integer :: omp_get_thread_num

     ithd = 0
#ifndef OMP_INSIDE_FFTW
!$   ithd = omp_get_thread_num()
#endif

!$OMP master
                                           call clock_sta(1430)
                                         ! call fapp_start("nlterm_realspcal",1430,1)
!$OMP end master
!!TBI!! !$OMP do schedule(dynamic)
      do ibuff = 0, nbuff-1
        do ix = 0, 2*nxw-1
          w3(ibuff,ix,0:global_ny) = wwdx(ibuff,ix,0:global_ny)
          w3(ibuff,ix,global_ny+1:nyw) = (0._DP, 0._DP) ! FFTW may destroy input array!
        end do
      end do

      call dfftw_execute_dft_c2r(planr_yb_y2zm(ithd), w3, dpdx)

      do ibuff = 0, nbuff-1
        do ix = 0, 2*nxw-1
          w3(ibuff,ix,0:global_ny) = wwdy(ibuff,ix,0:global_ny)
          w3(ibuff,ix,global_ny+1:nyw) = (0._DP, 0._DP) ! FFTW may destroy input array!
        end do
      end do

      call dfftw_execute_dft_c2r(planr_yb_y2zm(ithd), w3, dpdy)

!!TBI!! !$OMP end do nowait
!$OMP master
                                         ! call fapp_stop("nlterm_realspcal",1430,1)
                                           call clock_end(1430)
!$OMP end master

  END SUBROUTINE exb_backwardfft_y2zm


!--------------------------------------
  SUBROUTINE exb_realspcal_y2zm ( iv, dpdx, dpdy, dadx, dady, wwdx, wwdy, wwef )
!--------------------------------------
!     Calculate Poisson brackets for E x B term calculation (y2zm)

    integer, intent(in) :: iv
    real(kind=DP), intent(in), &
      dimension(0:nbuff-1,0:2*nxw-1,0:2*nyw-1) :: dpdx, dpdy, dadx, dady
    complex(kind=DP), intent(in), &
      dimension(0:nbuff-1,0:2*nxw-1,0:global_ny) :: wwdx, wwdy
    complex(kind=DP), intent(out), &
      dimension(0:nbuff-1,0:2*nxw-1,0:global_ny) :: wwef

    complex(kind=DP), dimension(0:nbuff-1,0:2*nxw-1,0:nyw) :: w3
    real(kind=DP), dimension(0:nbuff-1,0:2*nxw-1,0:2*nyw-1) :: dhdx, dhdy, pbxy
    real(kind=DP) :: cef, cs1
    integer :: ix, iy, ibuff, global_my

    integer :: ithd
!$  integer :: omp_get_thread_num

     ithd = 0
#ifndef OMP_INSIDE_FFTW
!$   ithd = omp_get_thread_num()
#endif

      cef = 1._DP / real(2*nxw*2*nyw, kind=DP)
      cs1 = sqrt(tau(ranks) / Anum(ranks))

!$OMP master
                                           call clock_sta(1430)
                                         ! call fapp_start("nlterm_realspcal",1430,1)
!$OMP end master
!!TBI!! !$OMP do schedule(dynamic)
      do global_my = 0, global_ny
        do ix = 0, 2*nxw-1
          do ibuff = 0, nbuff-1
            w3(ibuff,ix,global_my) = wwdx(ibuff,ix,global_my)
          end do
        end do
      end do
      do global_my = global_ny+1, nyw
        do ix = 0, 2*nxw-1
          do ibuff = 0, nbuff-1
            w3(ibuff,ix,global_my) = (0._DP, 0._DP) ! FFTW may destroy input array!
          end do
        end do
      end do

      call dfftw_execute_dft_c2r(planr_yb_y2zm(ithd), w3, dhdx)

      do global_my = 0, global_ny
        do ix = 0, 2*nxw-1
          do ibuff = 0, nbuff-1
            w3(ibuff,ix,global_my) = wwdy(ibuff,ix,global_my)
          end do
        end do
      end do
      do global_my = global_ny+1, nyw
        do ix = 0, 2*nxw-1
          do ibuff = 0, nbuff-1
            w3(ibuff,ix,global_my) = (0._DP, 0._DP) ! FFTW may destroy input array!
          end do
        end do
      end do

      call dfftw_execute_dft_c2r(planr_yb_y2zm(ithd), w3, dhdy)

      do iy = 0, 2*nyw-1
        do ix = 0, 2*nxw-1
          do ibuff = 0, nbuff-1
            pbxy(ibuff,ix,iy) = cef * ( & ! Normalization for 2D Forward FFT
                       (dpdx(ibuff,ix,iy) - cs1 * vl(iv) * dadx(ibuff,ix,iy)) * dhdy(ibuff,ix,iy) &
                     - (dpdy(ibuff,ix,iy) - cs1 * vl(iv) * dady(ibuff,ix,iy)) * dhdx(ibuff,ix,iy))
          end do
        end do
      end do

      call dfftw_execute_dft_r2c(planr_yf_y2zm(ithd), pbxy, w3)

      do global_my = 0, global_ny
        do ix = 0, 2*nxw-1
          do ibuff = 0, nbuff-1
            wwef(ibuff,ix,global_my) = w3(ibuff,ix,global_my)
          end do
        end do
      end do
!!TBI!! !$OMP end do nowait
!$OMP master
                                         ! call fapp_stop("nlterm_realspcal",1430,1)
                                           call clock_end(1430)
!$OMP end master

  END SUBROUTINE exb_realspcal_y2zm


!--------------------------------------
  SUBROUTINE exb_pack_zm2y ( wwef, wc4 )
!--------------------------------------
!     Data pack for E x B term calculation (y2zm)

    complex(kind=DP), intent(in), &
      dimension(0:nbuff-1,0:2*nxw-1,0:global_ny) :: wwef
    complex(kind=DP), intent(inout), &
      dimension(0:nbuff-1,-nx:nx,0:ny,0:nprocw-1) :: wc4

    complex(kind=DP), dimension(0:nbuff-1,0:global_ny,0:2*nxw-1) :: w1, w2
    complex(kind=DP), dimension(-nx:nx) :: ef
    integer :: mx, my, ibuff, iprocw, global_my

    integer :: ithd
!$  integer :: omp_get_thread_num

     ithd = 0
#ifndef OMP_INSIDE_FFTW
!$   ithd = omp_get_thread_num()
#endif

!$OMP master
                                           call clock_sta(1430)
                                         ! call fapp_start("nlterm_realspcal",1430,1)
!$OMP end master
!!TBI!! !$OMP do schedule(dynamic)
      do global_my = 0, global_ny
        do mx = 0, 2*nxw-1
          do ibuff = 0, nbuff-1
            w2(ibuff,global_my,mx) = wwef(ibuff,mx,global_my) ! FFTW may destroy input array!
          end do
        end do
      end do

      call dfftw_execute_dft(planr_xf_y2zm(ithd), w2, w1)

      do global_my = 0, global_ny
        iprocw = global_my / (ny+1)
        my = mod(global_my, ny+1)
        do mx = 0, nx
          do ibuff = 0, nbuff-1
            wc4(ibuff,mx,my,iprocw) = w1(ibuff,global_my,mx)
          end do
        end do
        do mx = -nx, -1
          do ibuff = 0, nbuff-1
            wc4(ibuff,mx,my,iprocw) = w1(ibuff,global_my,2*nxw+mx)
          end do
        end do
      end do
!!TBI!! !$OMP end do nowait
!$OMP master
                                         ! call fapp_stop("nlterm_realspcal",1430,1)
                                           call clock_end(1430)
!$OMP end master

  END SUBROUTINE exb_pack_zm2y


!--------------------------------------
  SUBROUTINE exb_transpose_zm2y ( wc4in, wc4out )
!--------------------------------------
!     Data transpose for E x B term calculation (y2zm)

    complex(kind=DP), intent(in), &
      dimension(0:nbuff-1,-nx:nx,0:ny,0:nprocw-1) :: wc4in
    complex(kind=DP), intent(out), &
      dimension(0:nbuff-1,-nx:nx,0:ny,0:nprocw-1) :: wc4out

!$OMP master
                                           call clock_sta(1440)
                                         ! call fapp_start("nlterm_alltoall2",1440,1)
      call MPI_Alltoall( wc4in,                 &
                         (2*nx+1)*(ny+1)*nbuff, &
                         MPI_DOUBLE_COMPLEX,    &
                         wc4out,                &
                         (2*nx+1)*(ny+1)*nbuff, &
                         MPI_DOUBLE_COMPLEX,    &
                         fft_comm_world,        &
                         ierr_mpi )
                                         ! call fapp_stop("nlterm_alltoall2",1440,1)
                                           call clock_end(1440)
!$OMP end master

  END SUBROUTINE exb_transpose_zm2y


!--------------------------------------
  SUBROUTINE exb_unpack_zm2y ( iv, wc4, ef )
!--------------------------------------
!     Data unpack for E x B term calculation (y2zm)

    integer, intent(in) :: iv
    complex(kind=DP), intent(in), &
      dimension(0:nbuff-1,-nx:nx,0:ny,0:nprocw-1) :: wc4
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: ef
    !complex(kind=DP), intent(inout), &
    !  dimension(-nx:nx,0:ny,-nz:nz-1,0:nm) :: ef
    !NOTE: A noncontiguous subarray as an argument of a subroutine induces a memory copy.
    !      When the subroutine is called in a OpenMP parallel region, 
    !      the copied subarray may be treated as a thread-private variable.

    integer :: mx, my, iz, im, izm, ibuff, iprocw

!$OMP master
                                           call clock_sta(1450)
                                         ! call fapp_start("nlterm_unpack",1450,1)
!$OMP end master
!!TBI!! !$OMP do collapse(2) schedule(dynamic,nchunk_zm)
      do im = 0, nm
        do iz = -nz, nz-1

         !%%% UNPACK: (kx,ky,(z*,m*)*)->(kx,ky*,z*,m*) %%%
          izm = (2*nz)*im + (iz + nz)
          ibuff = mod(izm, nbuff)
          iprocw = izm / nbuff
          do my = ist_y, iend_y
            do mx = -nx, nx
              ef(mx,my,iz,iv,im) = wc4(ibuff,mx,my,iprocw)
            end do
          end do
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end do
      end do
!!TBI!! !$OMP end do nowait
!$OMP master
                                         ! call fapp_stop("nlterm_unpack",1450,1)
                                           call clock_end(1450)
!$OMP end master

  END SUBROUTINE exb_unpack_zm2y


!--------------------------------------
  SUBROUTINE exb_estimate_maxvel_y2zm ( dpdx, dpdy, dadx, dady )
!--------------------------------------
!     Estimate time step restriction in each MPI processes

    real(kind=DP), intent(in), &
      dimension(0:nbuff-1,0:2*nxw-1,0:2*nyw-1) :: dpdx, dpdy, dadx, dady

    real(kind=DP) :: cs1, wv_nl
    integer :: ix, iy, ibuff, iv

      exb_maxvx_eachrank = eps
      exb_maxvy_eachrank = eps

      cs1 = sqrt(tau(ranks) / Anum(ranks))
      iv = 2*nv
!$OMP parallel default(none)                               &
!$OMP shared(exb_maxvx_eachrank,exb_maxvy_eachrank)        &
!$OMP shared(iv,ist_xw,iend_xw,dpdx,dpdy,dadx,dady,cs1,vl) &
!$OMP private(iy,ix,ibuff,wv_nl)

!$OMP do collapse(2) reduction(max:exb_maxvx_eachrank)
        do iy = 0, 2*nyw-1
          do ix = 0, 2*nxw-1
            do ibuff = 0, nbuff-1
              wv_nl = abs(dpdy(ibuff,ix,iy) - cs1 * vl(iv) * dady(ibuff,ix,iy))
              if (exb_maxvx_eachrank < wv_nl) exb_maxvx_eachrank = wv_nl
            end do
          end do
        end do
!$OMP end do

!$OMP do collapse(2) reduction(max:exb_maxvy_eachrank)
        do iy = 0, 2*nyw-1
          do ix = 0, 2*nxw-1
            do ibuff = 0, nbuff-1
              wv_nl = abs(dpdx(ibuff,ix,iy) - cs1 * vl(iv) * dadx(ibuff,ix,iy))
              if (exb_maxvy_eachrank < wv_nl) exb_maxvy_eachrank = wv_nl
            end do
          end do
        end do
!$OMP end do

!$OMP end parallel

  END SUBROUTINE exb_estimate_maxvel_y2zm


!--------------------------------------
  SUBROUTINE exb_NL_term_y2x( hh, psi, chi, ef )
!--------------------------------------
!  ExB nonlinear term calculation 

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: hh
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,0:nm) :: psi, chi
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: ef

    real(kind=DP), dimension(:,:,:,:), allocatable :: dpdx, dpdy, dadx, dady
    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: &
                                   wdx1o, wdy1o, wdx2o, wdy2o, wef3o, wef4o, &
                                   wdx1e, wdy1e, wdx2e, wdy2e, wef3e, wef4e
    integer :: iv, iprocw
    integer       :: mx, my
    integer, save :: iflg = 0

      allocate(dpdx(-nz:nz-1,0:nm,0:nxw_size,0:2*nyw-1))
      allocate(dpdy(-nz:nz-1,0:nm,0:nxw_size,0:2*nyw-1))
      allocate(dadx(-nz:nz-1,0:nm,0:nxw_size,0:2*nyw-1))
      allocate(dady(-nz:nz-1,0:nm,0:nxw_size,0:2*nyw-1))
      allocate(wdx1o(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1))
      allocate(wdy1o(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1))
      allocate(wdx2o(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1))
      allocate(wdy2o(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1))
      allocate(wef3o(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1))
      allocate(wef4o(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1))
      allocate(wdx1e(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1))
      allocate(wdy1e(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1))
      allocate(wdx2e(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1))
      allocate(wdy2e(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1))
      allocate(wef3e(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1))
      allocate(wef4e(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1))

      if( iflg == 0 ) then
        iflg = 1

!       uikx_y2x(:,:) = ( 0._DP, 0._DP )
!       uiky_y2x(:,:) = ( 0._DP, 0._DP )
        do mx = 0, nx
          do my = ist_y, iend_y
            uikx_y2x(my,mx) = kx(mx) * ui
            uiky_y2x(my,mx) = ky(my) * ui
          end do
        end do
        do mx = -nx, -1
          do my = ist_y, iend_y
            uikx_y2x(my,2*nxw+mx) = kx(mx) * ui
            uiky_y2x(my,2*nxw+mx) = ky(my) * ui
          end do
        end do
      end if

!$OMP parallel default(none)                      &
!$OMP shared(hh,psi,chi,ef,dpdx,dpdy,dadx,dady)   &
!$OMP shared(wdx1o,wdy1o,wdx2o,wdy2o,wef3o,wef4o) &
!$OMP shared(wdx1e,wdy1e,wdx2e,wdy2e,wef3e,wef4e) &
!$OMP private(iv,iprocw)

!$OMP workshare
      ef(:,:,:,:,:) = (0._DP, 0._DP)
      wdx1o(:,:,:,:,:) = (0._DP, 0._DP)
      wdy1o(:,:,:,:,:) = (0._DP, 0._DP)
      wdx1e(:,:,:,:,:) = (0._DP, 0._DP)
      wdy1e(:,:,:,:,:) = (0._DP, 0._DP)
      wef3o(:,:,:,:,:) = (0._DP, 0._DP)
      wef3e(:,:,:,:,:) = (0._DP, 0._DP)
      dpdx(:,:,:,:) = 0._DP
      dpdy(:,:,:,:) = 0._DP
      dadx(:,:,:,:) = 0._DP
      dady(:,:,:,:) = 0._DP
!$OMP end workshare


!!%%% Without overlap %%%
!      call exb_pack_psi_y2x(psi,wdx1o,wdy1o)
!!$OMP barrier
!      call exb_transpose_y2x(wdx1o,wdx2o)
!      call exb_transpose_y2x(wdy1o,wdy2o)
!!$OMP barrier
!      call exb_backwardfft_y2x(wdx2o,wdy2o,dpdx,dpdy)
!!$OMP barrier
!      call exb_pack_psi_y2x(chi,wdx1e,wdy1e)
!!$OMP barrier
!      call exb_transpose_y2x(wdx1e,wdx2e)
!      call exb_transpose_y2x(wdy1e,wdy2e)
!!$OMP barrier
!      call exb_backwardfft_y2x(wdx2e,wdy2e,dadx,dady)
!!$OMP barrier
!      do iv = 1, 2*nv
!        call exb_pack_hh_y2x(iv,hh,wdx1o,wdy1o)
!!$OMP barrier
!        call exb_transpose_y2x(wdx1o,wdx2o)
!        call exb_transpose_y2x(wdy1o,wdy2o)
!!$OMP barrier
!        call exb_realspcal_y2x(iv,dpdx,dpdy,dadx,dady,wdx2o,wdy2o,wef3o)
!!$OMP barrier
!        call exb_transpose_x2y(wef3o,wef4o)
!!$OMP barrier
!        call exb_unpack_x2y(iv,wef4o,ef)
!!$OMP barrier
!      end do
!!%%%%%%%%%%%%%%%%%%%%%%%


!%%% With overlap %%%
      call exb_pack_psi_y2x(psi,wdx1o,wdy1o)
!$OMP barrier
      call exb_transpose_y2x(wdx1o,wdx2o)
      call exb_transpose_y2x(wdy1o,wdy2o)
      call exb_pack_psi_y2x(chi,wdx1e,wdy1e)
!$OMP barrier
      call exb_transpose_y2x(wdx1e,wdx2e)
      call exb_transpose_y2x(wdy1e,wdy2e)
      call exb_backwardfft_y2x(wdx2o,wdy2o,dpdx,dpdy)
!$OMP barrier
      call exb_backwardfft_y2x(wdx2e,wdy2e,dadx,dady)
      do iv = 1, 2*nv+4
        if (mod(iv,2) == 1) then ! odd
          if (1+1<=iv .and. iv<=2*nv+1) call exb_transpose_y2x(wdx1e,wdx2e)
          if (1+1<=iv .and. iv<=2*nv+1) call exb_transpose_y2x(wdy1e,wdy2e)
          if (1+3<=iv .and. iv<=2*nv+3) call exb_transpose_x2y(wef3e,wef4e)
          if (1  <=iv .and. iv<=2*nv  ) call exb_pack_hh_y2x(iv,hh,wdx1o,wdy1o)
          if (1+2<=iv .and. iv<=2*nv+2) call exb_realspcal_y2x(iv-2,dpdx,dpdy,dadx,dady,wdx2o,wdy2o,wef3o)
          if (1+4<=iv .and. iv<=2*nv+4) call exb_unpack_x2y(iv-4,wef4o,ef)
        else                     ! even
          if (1+1<=iv .and. iv<=2*nv+1) call exb_transpose_y2x(wdx1o,wdx2o)
          if (1+1<=iv .and. iv<=2*nv+1) call exb_transpose_y2x(wdy1o,wdy2o)
          if (1+3<=iv .and. iv<=2*nv+3) call exb_transpose_x2y(wef3o,wef4o)
          if (1  <=iv .and. iv<=2*nv  ) call exb_pack_hh_y2x(iv,hh,wdx1e,wdy1e)
          if (1+2<=iv .and. iv<=2*nv+2) call exb_realspcal_y2x(iv-2,dpdx,dpdy,dadx,dady,wdx2e,wdy2e,wef3e)
          if (1+4<=iv .and. iv<=2*nv+4) call exb_unpack_x2y(iv-4,wef4e,ef)
        end if
!$OMP barrier
      end do
!%%%%%%%%%%%%%%%%%%%%

!$OMP end parallel

      call exb_estimate_maxvel_y2x(dpdx,dpdy,dadx,dady)

      deallocate(dpdx)
      deallocate(dpdy)
      deallocate(dadx)
      deallocate(dady)
      deallocate(wdx1o)
      deallocate(wdy1o)
      deallocate(wdx2o)
      deallocate(wdy2o)
      deallocate(wef3o)
      deallocate(wef4o)
      deallocate(wdx1e)
      deallocate(wdy1e)
      deallocate(wdx2e)
      deallocate(wdy2e)
      deallocate(wef3e)
      deallocate(wef4e)

  END SUBROUTINE exb_NL_term_y2x


!!--------------------------------------
!  SUBROUTINE exb_pack_y2x ( psi, wwdx, wwdy )
!!--------------------------------------
!!     Data pack for E x B term calculation (y2x)
!
!    complex(kind=DP), intent(in), &
!      dimension(-nx:nx,0:ny,-nz:nz-1,0:nm) :: psi
!    complex(kind=DP), intent(inout), &
!      dimension(0:ny,0:nxw_size,-nz:nz-1,0:nm,0:nprocw-1) :: wwdx, wwdy
!
!    complex(kind=DP), dimension(0:2*nxw-1) :: w1, w2
!    integer :: ix, my, iz, im, irank, ist_xw_g_rank, iend_xw_g_rank
!
!!$OMP master
!                                           call clock_sta(1410)
!                                         ! call fapp_start("nlterm_pack",1410,1)
!!$OMP end master
!!$OMP do collapse(3) schedule(dynamic)
!      do im = 0, nm
!        do iz = -nz, nz-1
!          do my = ist_y, iend_y
!
!           !%%% Backward x-FFT (kx,ky)->(x,ky) %%%
!            w1(0:nx) = ui * kx(0:nx) * psi(0:nx,my,iz,im)
!            w1(nx+1:2*nxw-nx-1) = (0._DP, 0._DP) ! FFTW may destroy input array!
!            w1(2*nxw-nx:2*nxw-1) = ui * kx(-nx:-1) * psi(-nx:-1,my,iz,im)
!            call dfftw_execute_dft(plan_x_backward, w1, w2)
!           !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!           !%%% PACK: (x,ky*,z*,m*)->(ky,x*,z*,m*) %%%
!            do irank = 0, nprocw-1
!              ist_xw_g_rank  = (nxw_size+1)*irank
!              iend_xw_g_rank = min( (nxw_size+1)*(irank+1)-1, (2*nxw-1) )
!              do ix = ist_xw_g_rank, iend_xw_g_rank
!                wwdx(my,ix-ist_xw_g_rank,iz,im,irank) = w2(ix)
!              enddo
!            enddo
!           !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  
!           !%%% Backward x-FFT (kx,ky)->(x,ky) %%%
!            w1(0:nx) = ui * ky(my) * psi(0:nx,my,iz,im)
!            w1(nx+1:2*nxw-nx-1) = (0._DP, 0._DP) ! FFTW may destroy input array!
!            w1(2*nxw-nx:2*nxw-1) = ui * ky(my) * psi(-nx:-1,my,iz,im)
!            call dfftw_execute_dft(plan_x_backward, w1, w2)
!           !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!           !%%% PACK: (x,ky*,z*,m*)->(ky,x*,z*,m*) %%%
!            do irank = 0, nprocw-1
!              ist_xw_g_rank  = (nxw_size+1)*irank
!              iend_xw_g_rank = min( (nxw_size+1)*(irank+1)-1, (2*nxw-1) )
!              do ix = ist_xw_g_rank, iend_xw_g_rank
!                wwdy(my,ix-ist_xw_g_rank,iz,im,irank) = w2(ix)
!              enddo
!            enddo
!           !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  
!          enddo
!        enddo
!      enddo
!!$OMP end do nowait
!!$OMP master
!                                         ! call fapp_stop("nlterm_pack",1410,1)
!                                           call clock_end(1410)
!!$OMP end master
!
!  END SUBROUTINE exb_pack_y2x


!--------------------------------------
  SUBROUTINE exb_pack_psi_y2x ( psi, wwdx, wwdy ) !done
!--------------------------------------
!     Data pack for E x B term calculation (y2x)

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,0:nm) :: psi
    complex(kind=DP), intent(inout), &
      dimension(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1) :: wwdx, wwdy

    complex(kind=DP), dimension(-nz:nz-1,0:nm,ist_y:iend_y,0:2*nxw-1) :: w1a, w2a
    complex(kind=DP), dimension(-nz:nz-1,0:nm,ist_y:iend_y,0:2*nxw-1) :: w1b, w2b
    integer :: ix, my, iz, im, irank, ist_xw_g_rank, iend_xw_g_rank
    complex(kind=DP) :: ww

    integer :: ithd
!$  integer :: omp_get_thread_num

     ithd = 0
#ifndef OMP_INSIDE_FFTW
!$   ithd = omp_get_thread_num()
#endif

!$OMP master
                                           call clock_sta(1410)
                                         ! call fapp_start("nlterm_pack",1410,1)
!$OMP end master
!!TBI!! !$OMP do collapse(2) schedule(dynamic,nchunk_zm)
      do my = ist_y, iend_y
        do ix = 0, nx
          do im = 0, nm
          do iz = -nz, nz-1
            w1a(iz,im,my,ix) = psi(ix,my,iz,im)
          enddo
          enddo
        enddo
      enddo
      do my = ist_y, iend_y
        do ix = nx+1, 2*nxw-nx-1
          do im = 0, nm
          do iz = -nz, nz-1
            w1a(iz,im,my,ix) = (0._DP, 0._DP) ! FFTW may destroy input array!
          enddo
          enddo
        enddo
      enddo
      do my = ist_y, iend_y
        do ix = -nx, -1
          do im = 0, nm
          do iz = -nz, nz-1
            w1a(iz,im,my,2*nxw+ix) = psi(ix,my,iz,im)
          enddo
          enddo
        enddo
      enddo

      do my = ist_y, iend_y
        do ix = 0, 2*nxw-1
          do im = 0, nm
          do iz = -nz, nz-1
            ww = w1a(iz,im,my,ix)
            w1a(iz,im,my,ix) = uikx_y2x(my,ix) * ww
            w1b(iz,im,my,ix) = uiky_y2x(my,ix) * ww
          enddo
          enddo
        enddo
      enddo

      call dfftw_execute_dft(planr_xb_y2x(ithd), w1a, w2a)
      call dfftw_execute_dft(planr_xb_y2x(ithd), w1b, w2b)

      do irank = 0, nprocw-1
        ist_xw_g_rank  = (nxw_size+1)*irank
        iend_xw_g_rank = min( (nxw_size+1)*(irank+1)-1, (2*nxw-1) )
        do my = ist_y, iend_y
          do ix = ist_xw_g_rank, iend_xw_g_rank
            do im = 0, nm
              do iz = -nz, nz-1
                wwdx(iz,im,ix-ist_xw_g_rank,my,irank) = w2a(iz,im,my,ix)
                wwdy(iz,im,ix-ist_xw_g_rank,my,irank) = w2b(iz,im,my,ix)
              enddo
            enddo
          enddo
        enddo
      enddo
!!TBI!! !$OMP end do nowait
!$OMP master
                                         ! call fapp_stop("nlterm_pack",1410,1)
                                           call clock_end(1410)
!$OMP end master

  END SUBROUTINE exb_pack_psi_y2x


!--------------------------------------
  SUBROUTINE exb_pack_hh_y2x ( iv, hh, wwdx, wwdy ) !done
!--------------------------------------
!     Data pack for E x B term calculation (y2x)

    integer, intent(in) :: iv
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: hh
    complex(kind=DP), intent(inout), &
      dimension(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1) :: wwdx, wwdy

    complex(kind=DP), dimension(-nz:nz-1,0:nm,ist_y:iend_y,0:2*nxw-1) :: w1a, w2a
    complex(kind=DP), dimension(-nz:nz-1,0:nm,ist_y:iend_y,0:2*nxw-1) :: w1b, w2b
    integer :: ix, my, iz, im, irank, ist_xw_g_rank, iend_xw_g_rank
    complex(kind=DP) :: ww

    integer :: ithd
!$  integer :: omp_get_thread_num

     ithd = 0
#ifndef OMP_INSIDE_FFTW
!$   ithd = omp_get_thread_num()
#endif

!$OMP master
                                           call clock_sta(1410)
                                         ! call fapp_start("nlterm_pack",1410,1)
!$OMP end master
!!TBI!! !$OMP do collapse(2) schedule(dynamic,nchunk_zm)
      do my = ist_y, iend_y
       do ix = 0, nx
        do im = 0, nm
          do iz = -nz, nz-1
            w1a(iz,im,my,ix) = hh(ix,my,iz,iv,im)
          enddo
        enddo
       enddo
      enddo
      do my = ist_y, iend_y
       do ix = nx+1, 2*nxw-nx-1
        do im = 0, nm
          do iz = -nz, nz-1
            w1a(iz,im,my,ix) = (0._DP, 0._DP) ! FFTW may destroy input array!
          enddo
        enddo
       enddo
      enddo
      do my = ist_y, iend_y
       do ix = -nx, -1
        do im = 0, nm
          do iz = -nz, nz-1
            w1a(iz,im,my,2*nxw+ix) = hh(ix,my,iz,iv,im)
          enddo
        enddo
       enddo
      enddo

      do my = ist_y, iend_y
       do ix = 0, 2*nxw-1
        do im = 0, nm
          do iz = -nz, nz-1
            ww = w1a(iz,im,my,ix)
            w1a(iz,im,my,ix) = uikx_y2x(my,ix) * ww
            w1b(iz,im,my,ix) = uiky_y2x(my,ix) * ww
          enddo
        enddo
       enddo
      enddo

      call dfftw_execute_dft(planr_xb_y2x(ithd), w1a, w2a)
      call dfftw_execute_dft(planr_xb_y2x(ithd), w1b, w2b)

      do irank = 0, nprocw-1
        ist_xw_g_rank  = (nxw_size+1)*irank
        iend_xw_g_rank = min( (nxw_size+1)*(irank+1)-1, (2*nxw-1) )
        do my = ist_y, iend_y
          do ix = ist_xw_g_rank, iend_xw_g_rank
            do im = 0, nm
              do iz = -nz, nz-1
                wwdx(iz,im,ix-ist_xw_g_rank,my,irank) = w2a(iz,im,my,ix)
                wwdy(iz,im,ix-ist_xw_g_rank,my,irank) = w2b(iz,im,my,ix)
              enddo
            enddo
          enddo
        enddo
      enddo
!!TBI!! !$OMP end do nowait
!$OMP master
                                         ! call fapp_stop("nlterm_pack",1410,1)
                                           call clock_end(1410)
!$OMP end master

  END SUBROUTINE exb_pack_hh_y2x


!--------------------------------------
  SUBROUTINE exb_transpose_y2x ( wwin, wwout ) !done
!--------------------------------------
!     Data transpose for E x B term calculation (y2x)

    complex(kind=DP), intent(in), &
      dimension(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1) :: wwin
    complex(kind=DP), intent(out), &
      dimension(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1) :: wwout

!$OMP master
                                           call clock_sta(1420)
                                         ! call fapp_start("nlterm_alltoall1",1420,1)
      call MPI_Alltoall( wwin,                              &
                         (ny+1)*(nxw_size+1)*(2*nz)*(nm+1), &
                         MPI_DOUBLE_COMPLEX,                &
                         wwout,                             &
                         (ny+1)*(nxw_size+1)*(2*nz)*(nm+1), &
                         MPI_DOUBLE_COMPLEX,                &
                         fft_comm_world,                    &
                         ierr_mpi )
                                         ! call fapp_stop("nlterm_alltoall1",1420,1)
                                           call clock_end(1420)
!$OMP end master

  END SUBROUTINE exb_transpose_y2x


!--------------------------------------
  SUBROUTINE exb_backwardfft_y2x ( wwdx, wwdy, dpdx, dpdy ) !done
!--------------------------------------
!     Backward FFT of field(psi,chi) for E x B term calculation (y2x)

    complex(kind=DP), intent(in), &
      dimension(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1) :: wwdx, wwdy
    real(kind=DP), intent(inout), &
      dimension(-nz:nz-1,0:nm,0:nxw_size,0:2*nyw-1) :: dpdx, dpdy

    complex(kind=DP), dimension(-nz:nz-1,0:nm,ist_xw:iend_xw,0:nyw) :: w3
    real(kind=DP), dimension(-nz:nz-1,0:nm,ist_xw:iend_xw,0:2*nyw-1) :: wtmp
    integer :: ix, my, iz, im, irank, ist_y_g_rank, iend_y_g_rank

    integer :: ithd
!$  integer :: omp_get_thread_num

     ithd = 0
#ifndef OMP_INSIDE_FFTW
!$   ithd = omp_get_thread_num()
#endif

!$OMP master
                                           call clock_sta(1430)
                                         ! call fapp_start("nlterm_realspcal",1430,1)
!$OMP end master
!!TBI!! !$OMP do collapse(2) schedule(dynamic,nchunk_zm)
      do irank = 0, nprocw-1
        ist_y_g_rank  = (ny+1)*irank
        iend_y_g_rank = min( (ny+1)*(irank+1)-1, global_ny )
        do my = ist_y_g_rank, iend_y_g_rank
          do ix = ist_xw, iend_xw
            do im = 0, nm
              do iz = -nz, nz-1
                w3(iz,im,ix,my) = wwdx(iz,im,ix,my-ist_y_g_rank,irank)
              end do
            end do
          end do
        end do
      end do
      do my = global_ny+1, nyw
        do ix = ist_xw, iend_xw
          do im = 0, nm
          do iz = -nz, nz-1
            w3(iz,im,ix,my) = ( 0._DP, 0._DP ) ! FFTW may destroy input array!
          end do
          end do
        end do
      end do

!!TBI!! call dfftw_execute_dft_c2r(plan_yb_y2x(ithd), w3, dpdx(0,ist_xw,iz,im))
      call dfftw_execute_dft_c2r(planr_yb_y2x(ithd), w3, wtmp)
      do my = 0, 2*nyw-1
       do ix = ist_xw, iend_xw
        do im = 0, nm
         do iz = -nz, nz-1
          dpdx(iz,im,ix,my) = wtmp(iz,im,ix,my)
         end do
        end do
       end do
      end do

      do irank = 0, nprocw-1
        ist_y_g_rank  = (ny+1)*irank
        iend_y_g_rank = min( (ny+1)*(irank+1)-1, global_ny )
        do my = ist_y_g_rank, iend_y_g_rank
          do ix = ist_xw, iend_xw
            do im = 0, nm
              do iz = -nz, nz-1
                w3(iz,im,ix,my) = wwdy(iz,im,ix,my-ist_y_g_rank,irank)
              end do
            end do
          end do
        end do
      end do
      do my = global_ny+1, nyw
        do ix = ist_xw, iend_xw
          do im = 0, nm
          do iz = -nz, nz-1
            w3(iz,im,ix,my) = ( 0._DP, 0._DP ) ! FFTW may destroy input array!
          end do
          end do
        end do
      end do

!!TBI!! call dfftw_execute_dft_c2r(plan_yb_y2x(ithd), w3, dpdy(0,ist_xw,iz,im))
      call dfftw_execute_dft_c2r(planr_yb_y2x(ithd), w3, wtmp)
      do my = 0, 2*nyw-1
       do ix = ist_xw, iend_xw
        do im = 0, nm
         do iz = -nz, nz-1
          dpdy(iz,im,ix,my) = wtmp(iz,im,ix,my)
         end do
        end do
       end do
      end do

!!TBI!! !$OMP end do nowait
!$OMP master
                                         ! call fapp_stop("nlterm_realspcal",1430,1)
                                           call clock_end(1430)
!$OMP end master

  END SUBROUTINE exb_backwardfft_y2x


!--------------------------------------
  SUBROUTINE exb_realspcal_y2x ( iv, dpdx, dpdy, dadx, dady, wwdx, wwdy, wwef ) !done
!--------------------------------------
!     Calculate Poisson brackets for E x B term calculation (y2x)

    integer, intent(in) :: iv
    real(kind=DP), intent(in), &
      dimension(-nz:nz-1,0:nm,0:nxw_size,0:2*nyw-1) :: dpdx, dpdy, dadx, dady
    complex(kind=DP), intent(in), &
      dimension(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1) :: wwdx, wwdy
    complex(kind=DP), intent(inout), &
      dimension(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1) :: wwef

    complex(kind=DP), dimension(-nz:nz-1,0:nm,ist_xw:iend_xw,0:nyw) :: w3
    real(kind=DP), dimension(-nz:nz-1,0:nm,ist_xw:iend_xw,0:2*nyw-1) :: dhdx, dhdy, pbxy
    real(kind=DP) :: cef, cs1
    integer :: my, ix, iy, iz, im, irank, ist_y_g_rank, iend_y_g_rank

    integer :: ithd
!$  integer :: omp_get_thread_num

     ithd = 0
#ifndef OMP_INSIDE_FFTW
!$   ithd = omp_get_thread_num()
#endif

      cef = 1._DP / real(2*nxw*2*nyw, kind=DP)
      cs1 = sqrt(tau(ranks) / Anum(ranks))

!$OMP master
                                           call clock_sta(1430)
                                         ! call fapp_start("nlterm_realspcal",1430,1)
!$OMP end master
!!TBI!! !$OMP do collapse(2) schedule(dynamic,nchunk_zm)
      do irank = 0, nprocw-1
        ist_y_g_rank  = (ny+1)*irank
        iend_y_g_rank = min( (ny+1)*(irank+1)-1, global_ny )
        do my = ist_y_g_rank, iend_y_g_rank
          do ix = ist_xw, iend_xw
            do im = 0, nm
              do iz = -nz, nz-1
                w3(iz,im,ix,my) = wwdx(iz,im,ix,my-ist_y_g_rank,irank)
              end do
            end do
          end do
        end do
      end do
      do my = global_ny+1, nyw
       do ix = ist_xw, iend_xw
        do im = 0, nm
          do iz = -nz, nz-1
            w3(iz,im,ix,my) = ( 0._DP, 0._DP ) ! FFTW may destroy input array!
          end do
        end do
       end do
      end do

      call dfftw_execute_dft_c2r(planr_yb_y2x(ithd), w3, dhdx)

      do irank = 0, nprocw-1
        ist_y_g_rank  = (ny+1)*irank
        iend_y_g_rank = min( (ny+1)*(irank+1)-1, global_ny )
        do my = ist_y_g_rank, iend_y_g_rank
          do ix = ist_xw, iend_xw
            do im = 0, nm
              do iz = -nz, nz-1
                w3(iz,im,ix,my) = wwdy(iz,im,ix,my-ist_y_g_rank,irank)
              end do
            end do
          end do
        end do
      end do
      do my = global_ny+1, nyw
       do ix = ist_xw, iend_xw
        do im = 0, nm
          do iz = -nz, nz-1
            w3(iz,im,ix,my) = ( 0._DP, 0._DP ) ! FFTW may destroy input array!
          end do
        end do
       end do
      end do

      call dfftw_execute_dft_c2r(planr_yb_y2x(ithd), w3, dhdy)

!NEC$ NOINTERCHANGE
      do iy = 0, 2*nyw-1
!NEC$ NOINTERCHANGE
        do ix = ist_xw, iend_xw
!NEC$ NOINTERCHANGE
          do im = 0, nm
!NEC$ NOINTERCHANGE
            do iz = -nz, nz-1
              pbxy(iz,im,ix,iy) = cef * ( & ! Normalization for 2D Forward FFT
                         (dpdx(iz,im,ix,iy) - cs1 * vl(iv) * dadx(iz,im,ix,iy)) * dhdy(iz,im,ix,iy) &
                       - (dpdy(iz,im,ix,iy) - cs1 * vl(iv) * dady(iz,im,ix,iy)) * dhdx(iz,im,ix,iy))
            end do
          end do
        end do
      end do

      call dfftw_execute_dft_r2c(planr_yf_y2x(ithd), pbxy, w3)

      do irank = 0, nprocw-1
        ist_y_g_rank  = (ny+1)*irank
        iend_y_g_rank = min( (ny+1)*(irank+1)-1, global_ny )
        do my = ist_y_g_rank, iend_y_g_rank
          do ix = ist_xw, iend_xw
            do im = 0, nm
              do iz = -nz, nz-1
                wwef(iz,im,ix,my-ist_y_g_rank,irank) = w3(iz,im,ix,my)
              end do
            end do
          end do
        end do
      end do
!!TBI!! !$OMP end do nowait
!$OMP master
                                         ! call fapp_stop("nlterm_realspcal",1430,1)
                                           call clock_end(1430)
!$OMP end master

  END SUBROUTINE exb_realspcal_y2x


!--------------------------------------
  SUBROUTINE exb_transpose_x2y ( wwin, wwout ) !done
!--------------------------------------
!     Data transpose for E x B term calculation (y2x)

    complex(kind=DP), intent(in), &
      dimension(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1) :: wwin
    complex(kind=DP), intent(out), &
      dimension(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1) :: wwout

!$OMP master
                                           call clock_sta(1440)
                                         ! call fapp_start("nlterm_alltoall2",1440,1)
      call MPI_Alltoall( wwin,                              &
                         (ny+1)*(nxw_size+1)*(2*nz)*(nm+1), &
                         MPI_DOUBLE_COMPLEX,                &
                         wwout,                             &
                         (ny+1)*(nxw_size+1)*(2*nz)*(nm+1), &
                         MPI_DOUBLE_COMPLEX,                &
                         fft_comm_world,                    &
                         ierr_mpi )
                                         ! call fapp_stop("nlterm_alltoall2",1440,1)
                                           call clock_end(1440)
!$OMP end master

  END SUBROUTINE exb_transpose_x2y


!--------------------------------------
  SUBROUTINE exb_unpack_x2y ( iv, wwef, ef ) !done
!--------------------------------------
!     Data unpack for E x B term calculation (y2x)

    integer, intent(in) :: iv
    complex(kind=DP), intent(in), &
      dimension(-nz:nz-1,0:nm,0:nxw_size,0:ny,0:nprocw-1) :: wwef
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: ef
    !complex(kind=DP), intent(inout), &
    !  dimension(-nx:nx,0:ny,-nz:nz-1,0:nm) :: ef
    !NOTE: A noncontiguous subarray as an argument of a subroutine induces a memory copy.
    !      When the subroutine is called in a OpenMP parallel region, 
    !      the copied subarray may be treated as a thread-private variable.

    complex(kind=DP), dimension(-nz:nz-1,0:nm,ist_y:iend_y,0:2*nxw-1) :: w1, w2
    integer :: ix, my, iz, im, irank, ist_xw_g_rank, iend_xw_g_rank

    integer :: ithd
!$  integer :: omp_get_thread_num

     ithd = 0
#ifndef OMP_INSIDE_FFTW
!$   ithd = omp_get_thread_num()
#endif

!$OMP master
                                           call clock_sta(1450)
                                         ! call fapp_start("nlterm_unpack",1450,1)
!$OMP end master
!!TBI!! !$OMP do collapse(2) schedule(dynamic,nchunk_zm)
      do irank = 0, nprocw-1
        ist_xw_g_rank  = (nxw_size+1)*irank
        iend_xw_g_rank = min( (nxw_size+1)*(irank+1)-1, (2*nxw-1) )
        do my = ist_y, iend_y
          do ix = ist_xw_g_rank, iend_xw_g_rank
            do im = 0, nm
              do iz = -nz, nz-1
                w2(iz,im,my,ix) = wwef(iz,im,ix-ist_xw_g_rank,my,irank)
              enddo
            end do
          end do
        end do
      end do

      call dfftw_execute_dft(planr_xf_y2x(ithd), w2, w1)

      do my = ist_y, iend_y
       do ix = 0, nx
        do im = 0, nm
          do iz = -nz, nz-1
            ef(ix,my,iz,iv,im) = w1(iz,im,my,ix)
          end do
        end do
       end do
      end do
      do my = ist_y, iend_y
       do ix = -nx, -1
        do im = 0, nm
          do iz = -nz, nz-1
            ef(ix,my,iz,iv,im) = w1(iz,im,my,2*nxw+ix)
          end do
        end do
       end do
      end do
!!TBI!! !$OMP end do nowait
!$OMP master
                                         ! call fapp_stop("nlterm_unpack",1450,1)
                                           call clock_end(1450)
!$OMP end master

  END SUBROUTINE exb_unpack_x2y


!--------------------------------------
  SUBROUTINE exb_estimate_maxvel_y2x ( dpdx, dpdy, dadx, dady ) !done
!--------------------------------------
!     Estimate time step restriction in each MPI processes

    real(kind=DP), intent(in), &
      dimension(-nz:nz-1,0:nm,0:nxw_size,0:2*nyw-1) :: dpdx, dpdy, dadx, dady

    real(kind=DP) :: cs1, wv_nl
    integer :: ix, iy, iz, iv, im

      exb_maxvx_eachrank = eps
      exb_maxvy_eachrank = eps

      cs1 = sqrt(tau(ranks) / Anum(ranks))
      im = 0
        iv = 2*nv
!$OMP parallel default(none)                                  &
!$OMP shared(exb_maxvx_eachrank,exb_maxvy_eachrank)           &
!$OMP shared(im,iv,ist_xw,iend_xw,dpdx,dpdy,dadx,dady,cs1,vl) &
!$OMP private(iy,ix,iz,wv_nl)

!$OMP do collapse(2) reduction(max:exb_maxvx_eachrank)
          do iy = 0, 2*nyw-1
            do ix = ist_xw, iend_xw
              do iz = -nz, nz-1
                wv_nl = abs(dpdy(iz,im,ix,iy) - cs1 * vl(iv) * dady(iz,im,ix,iy))
                if (exb_maxvx_eachrank < wv_nl) exb_maxvx_eachrank = wv_nl
              end do
            end do
          end do
!$OMP end do

!$OMP do collapse(2) reduction(max:exb_maxvy_eachrank)
          do iy = 0, 2*nyw-1
            do ix = ist_xw, iend_xw
              do iz = -nz, nz-1
                wv_nl = abs(dpdx(iz,im,ix,iy) - cs1 * vl(iv) * dadx(iz,im,ix,iy))
                if (exb_maxvy_eachrank < wv_nl) exb_maxvy_eachrank = wv_nl
              end do
            end do
          end do
!$OMP end do

!$OMP end parallel

  END SUBROUTINE exb_estimate_maxvel_y2x


END MODULE GKV_exb
