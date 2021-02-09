MODULE GKV_out
!-------------------------------------------------------------------------------
!
!    Data writing
!
!    Update history of gkvp_out.f90
!    --------------
!      gkvp_f0.60 (S. Maeyama, Jan 2021)
!        - Use fileio module to switch Fortran/NetCDF binary output.
!      gkvp_f0.57 (S. Maeyama, Oct 2020)
!        - Version number f0.57 is removed from filename.
!        - Frequency analysis is applied only when calc_type="lin_freq".
!        - menegy=0 when beta=0.
!
!-------------------------------------------------------------------------------

  use GKV_header
  use GKV_mpienv
  use GKV_intgrl, only: intgrl_thet, intgrl_fsrf, &
                        intgrl_v0_moment, intgrl_v2_moment
  use GKV_fld, only: fld_emfield_hh
  use GKV_trans, only: trans_sum, trans_triad
  use GKV_freq, only: freq_write_frq, freq_write_dsp
  use GKV_advnc, only: caldlt_rev
  use GKV_colliimp, only: colliimp_colli
  use GKV_dtc,   only: flag_time_advnc, flag_time_split
  use GKV_tips,  only: tips_flush, tips_rescale_for_linear_runs
  !fj start 202010
  use GKV_fileio
  !fj end 202010

  implicit none

  private

  public   out_cntrl, out_contnu


CONTAINS


!--------------------------------------
  SUBROUTINE out_cntrl( ff, phi, Al, hh, time, id )
!--------------------------------------

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: phi, Al
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: hh
    real(kind=DP), intent(in) :: time
    integer, intent(in) :: id

    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: dh, cf, ef
    real(kind=DP), save :: tout_fxv, tout_ptn, tout_eng
    integer :: flag_updated

      allocate( dh(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) )
      allocate( cf(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) )
      allocate( ef(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) )

      flag_updated = 0

      if( id == 0 ) then

        if ( time == 0._DP ) then
          if (flag_updated == 0) call update_dh(ff, phi, al, hh, dh, cf, ef, flag_updated)
          call wrt ( ff, phi, Al, dh, cf, time, 0 )
          call wrt ( ff, phi, Al, dh, cf, time, 1 )
          call wrt ( ff, phi, Al, dh, cf, time, 2 )
        end if

        tout_fxv  = ( int( ( time + eps )/dtout_fxv ) + 1 ) * dtout_fxv
        tout_ptn  = ( int( ( time + eps )/dtout_ptn ) + 1 ) * dtout_ptn
        tout_eng  = ( int( ( time + eps )/dtout_eng ) + 1 ) * dtout_eng
 
      else if( id == 1 ) then

        if ( time >= tout_fxv - eps ) then
          tout_fxv   = tout_fxv + dtout_fxv
          if (flag_updated == 0) call update_dh(ff, phi, al, hh, dh, cf, ef, flag_updated)
          call wrt ( ff, phi, Al, dh, cf, time, 0 )
          write( olog, * ) " # delta-f data output at time = ", time
        end if

        if ( time >= tout_ptn - eps ) then
          tout_ptn   = tout_ptn + dtout_ptn
          if (flag_updated == 0) call update_dh(ff, phi, al, hh, dh, cf, ef, flag_updated)
          call wrt ( ff, phi, Al, dh, cf, time, 1 )
          write( olog, * ) " # field data output at time = ", time
        end if

        if ( time >= tout_eng - eps ) then
          tout_eng   = tout_eng + dtout_eng
          if (flag_updated == 0) call update_dh(ff, phi, al, hh, dh, cf, ef, flag_updated)
          call wrt ( ff, phi, Al, dh, cf, time, 2 )
        end if

      else if( id == 2 ) then

        if (flag_updated == 0) call update_dh(ff, phi, al, hh, dh, cf, ef, flag_updated)
        call out_contnu ( ff, time )
        !- OUTPUT ascii data hst/*.dsp.* for linear runs
        if ( trim(calc_type) == "lin_freq" ) then
           call freq_write_dsp
        end if

      end if

      if (flag_time_split == 0) then
        if ( trim(calc_type) == "linear" .or. &
             trim(calc_type) == "lin_freq" ) then
          call tips_rescale_for_linear_runs(ff, phi, Al, hh, time)
        end if
      end if

      deallocate( dh )
      deallocate( cf )
      deallocate( ef )


  END SUBROUTINE out_cntrl


!--------------------------------------
SUBROUTINE update_dh( ff, phi, Al, hh, dh, cf, ef, flag_updated )
!--------------------------------------

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: phi, Al
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: hh
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: dh, cf, ef
    integer, intent(out) :: flag_updated

    character(15) :: colliflag


      if (flag_time_advnc == 1) then ! Operator split + implicit collision
        if (flag_time_split == 1) then ! dt/2 integration for 2nd-order split
          call colliimp_colli(0.5_DP*dt, ff, phi, al, hh)
                                               !%%% For debug %%%
                                               ! if (rankg==0) write(olog,*) &
                                               !       "half-step for output"
                                               !%%%%%%%%%%%%%%%%%
        end if
        flag_time_split = 0 ! flag_time_split==0 means you have physical quantities at time
      end if

      colliflag = "collisional"
      call caldlt_rev( colliflag, ff, phi, Al, hh, dh, cf, ef )

      flag_updated = 1

END SUBROUTINE update_dh


!--------------------------------------
  SUBROUTINE out_contnu ( ff, time )
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    real(kind=DP), intent(in) :: time

    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: wf
    integer ::mx, my, iz, iv, im


      allocate( wf(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) )

!$OMP parallel workshare
      wf(:,:,:,:,:) = ( 0._DP, 0._DP )
!$OMP end parallel workshare

!$OMP parallel
      do im = 0, nm
!$OMP do
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                wf(mx,my,iz,iv,im) = ff(mx,my,iz,iv,im)
              end do
            end do
          end do
        end do
!$OMP end do nowait
      end do
!$OMP end parallel

      !fj start 202010
      !rewind ocnt
      !write( unit=ocnt ) time, wf
      call fileio_write_cnt( wf, time )
      !fj end time202010

      deallocate( wf )


  END SUBROUTINE out_contnu


!--------------------------------------
  SUBROUTINE wrt ( ff, phi, Al, dh, cf, time, id )
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: phi, Al
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: dh, cf
    real(kind=DP), intent(in) :: time
    integer, intent(in) :: id

    complex(kind=DP), dimension(:,:,:,:), allocatable :: fout
    real(kind=DP), dimension(-nx:nx,0:ny) :: entrpy, fenegy, menegy,          &
                                             peint, pmint, neint, nmint, dcd, &
                                             pflux_es, pflux_em, eflux_es, eflux_em
    real(kind=DP), dimension(0:(ny+1)*nprocw-1) :: mode_y
    real(kind=DP) :: totl
    integer :: mx, my, iv, im


      if( id == 0 ) then

        allocate( fout(-nx:nx,0:ny,1:2*nv,0:nm) )

        !- OUTPUT binary data fxv/*.fxv.* -
!$OMP parallel workshare
        fout(:,:,:,:) = ( 0._DP, 0._DP )
!$OMP end parallel workshare
!$OMP parallel do collapse(2) private(mx,my,iv,im)
        do im = 0, nm
          do iv = 1, 2*nv
            do my = ist_y, iend_y
              do mx = -nx, nx
                fout(mx,my,iv,im) = ff(mx,my,-nz  ,iv,im)
              end do
            end do
          end do
        end do
        !fj start 202010
        !write( unit=ofxv ) time, fout
        call fileio_write_fxv( fout, time )
        !fj end 202010

        deallocate( fout )

      else if( id == 1 ) then

        !- OUTPUT binary data phi/*.phi.* and phi/*.Al.*
        !fj start 202011
        !if ( ranks == 0 .AND. vel_rank == 0 ) then
          !write( unit=ophi ) time, phi
          !write( unit=oAl  ) time, Al
        !end if
        call fileio_write_phi( phi, time )
        call fileio_write_Al ( Al,  time )
        !fj end 202011

        !- OUTPUT binary data phi/*.mom.*
        call write_moments ( ff, time )

        !- OUTPUT binary data phi/*.tri.*
        if( trim(calc_type) == "nonlinear" ) then
          if ( num_triad_diag > 0 ) then
            call trans_triad ( time, ff, phi, Al )
          end if
        end if

      else if( id == 2 ) then

        !- OUTPUT ascii data hst/*.bln.*
        call balance ( ff, phi, Al, dh, cf, time,                          &
                       entrpy, fenegy, menegy, peint, pmint, neint, nmint, &
                       dcd, pflux_es, pflux_em, eflux_es, eflux_em )

        !- OUTPUT binary data phi/*.trn.*
        !fj start 202011
        !if ( zsp_rank == 0 .and. vel_rank == 0 ) then
          !write( unit=otrn ) time, entrpy, fenegy, menegy,    &
          !                   peint, pmint, neint, nmint, dcd, &
          !                   pflux_es, pflux_em, eflux_es, eflux_em
        !end if
        call fileio_write_trn( entrpy, fenegy, menegy,    &
             peint, pmint, neint, nmint, dcd, &
             pflux_es, pflux_em, eflux_es, eflux_em, time )
        !fj end 202011

        !- OUTPUT ascii data hst/*.eng.*, *.men.*, *.wes.*, *.wem.*,
        !                        *.ges.*, *.gem.*, *.qes.*, *.qem.*
        call mode_energy ( phi, totl, mode_y )
        if( rankg == 0 ) then
          write( unit=oeng, fmt="(f15.8, SP, 9999ES24.15e3)" ) &
                               time, totl, mode_y(0:global_ny)
        end if

        call mode_energy ( Al, totl, mode_y )
        if ( rankg == 0 ) then
          write( unit=omen, fmt="(f15.8, SP, 9999ES24.15e3)" ) &
                               time, totl, mode_y(0:global_ny)
        end if

        call calc_kyspectrum ( fenegy, totl, mode_y )
        if ( rankg == 0 ) then
          write( unit=owes, fmt="(f15.8, SP, 9999ES24.15e3)" ) &
                               time, totl, mode_y(0:global_ny)
        end if

        call calc_kyspectrum ( menegy, totl, mode_y )
        if ( rankg == 0 ) then
          write( unit=owem, fmt="(f15.8, SP, 9999ES24.15e3)" ) &
                               time, totl, mode_y(0:global_ny)
        end if

        call calc_kyspectrum ( pflux_es, totl, mode_y )
        if ( rank == 0 ) then
          write( unit=oges, fmt="(f15.8, SP, 9999ES24.15e3)" ) &
                               time, totl, mode_y(0:global_ny)
        end if

        call calc_kyspectrum ( pflux_em, totl, mode_y )
        if ( rank == 0 ) then
          write( unit=ogem, fmt="(f15.8, SP, 9999ES24.15e3)" ) &
                               time, totl, mode_y(0:global_ny)
        end if

        call calc_kyspectrum ( eflux_es, totl, mode_y )
        if ( rank == 0 ) then
          write( unit=oqes, fmt="(f15.8, SP, 9999ES24.15e3)" ) &
                               time, totl, mode_y(0:global_ny)
        end if
        call calc_kyspectrum ( eflux_em, totl, mode_y )

        if ( rank == 0 ) then
          write( unit=oqem, fmt="(f15.8, SP, 9999ES24.15e3)" ) &
                               time, totl, mode_y(0:global_ny)
        end if

        !- OUTPUT ascii data hst/*.frq.* for linear runs
        if ( trim(calc_type) == "lin_freq" ) then
          call freq_write_frq ( time, phi )
        end if

      end if


  END SUBROUTINE wrt


!--------------------------------------
  SUBROUTINE mode_energy ( phi, totl, mode_y )
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: phi
    real(kind=DP), intent(out), &
      dimension(0:(ny+1)*nprocw-1)    :: mode_y
    real(kind=DP), intent(out)        :: totl

    real(kind=DP), dimension(:,:,:), allocatable :: wr3
    real(kind=DP), dimension(-nx:nx,0:ny)        :: wr2
    real(kind=DP), dimension(0:ny)               :: mode_wk
    integer  ::  mx, my, iz

      allocate( wr3(-nx:nx,0:ny,-nz:nz-1) )

      totl    = 0._DP
      mode_y  = 0._DP
      mode_wk = 0._DP

!$OMP parallel do
      do iz = -nz, nz-1
        do my = ist_y, iend_y
          do mx = -nx, nx
            wr3(mx,my,iz) = real( phi(mx,my,iz) * conjg( phi(mx,my,iz) )  &
                                 , kind=DP )
          end do
        end do
      end do

      call intgrl_thet ( wr3, wr2 )

!$OMP parallel do reduction(+:mode_wk)
      do my = ist1_y, iend_y
        do mx = -nx, nx
          mode_wk(my) = mode_wk(my) + wr2(mx,my)
        end do
      end do

      if( rankw == 0 ) then
        my = 0
          do mx = 0, nx
            mode_wk(my) = mode_wk(my) + wr2(mx,my)
          end do
      endif

      call MPI_Gather( mode_wk, ny+1, MPI_DOUBLE_PRECISION, &
                       mode_y,  ny+1, MPI_DOUBLE_PRECISION, &
                       0, fft_comm_world, ierr_mpi )

      do my = 0, global_ny
        totl = totl + mode_y(my)
      end do

      deallocate( wr3 )

  END SUBROUTINE mode_energy


!--------------------------------------
  SUBROUTINE calc_kyspectrum ( fenegy, totl, mode_y )
!--------------------------------------

    real(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny)          :: fenegy
    real(kind=DP), intent(out), &
      dimension(0:(ny+1)*nprocw-1)    :: mode_y
    real(kind=DP), intent(out)        :: totl

    real(kind=DP), dimension(0:ny) :: mode_wk
    integer :: mx, my


      totl    = 0._DP
      mode_y  = 0._DP
      mode_wk = 0._DP

!$OMP parallel do reduction(+:mode_wk)
      do my = ist1_y, iend_y
        do mx = -nx, nx
          mode_wk(my) = mode_wk(my) + fenegy(mx,my)
        end do
      end do

      if( rankw == 0 ) then
        my = 0
          do mx = 0, nx
            mode_wk(my) = mode_wk(my) + fenegy(mx,my)
          end do
      endif

      call MPI_Gather( mode_wk, ny+1, MPI_DOUBLE_PRECISION, &
                       mode_y,  ny+1, MPI_DOUBLE_PRECISION, &
                       0, fft_comm_world, ierr_mpi )

      mode_y(:) = 2._DP * mode_y(:)
      do my = 0, global_ny
        totl = totl + mode_y(my)
      end do


  END SUBROUTINE calc_kyspectrum


!--------------------------------------
  SUBROUTINE balance ( ff, phi, Al, dh, cf, time,                    &
                       entrpy, fenegy, menegy, peint, pmint, neint, nmint, &
                       dcd, pflux_es, pflux_em, eflux_es, eflux_em )
!--------------------------------------
!     Check the entropy balance equation

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: phi, Al
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: dh, cf
    real(kind=DP), intent(in) :: time
    real(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny) :: entrpy, fenegy, menegy, peint, pmint, neint, nmint, dcd, &
                                pflux_es, pflux_em, eflux_es, eflux_em

    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: wf
    complex(kind=DP), dimension(:,:,:), allocatable :: dens, upara, pres, qpara, ni, wc3
    complex(kind=DP), dimension(-nx:nx,0:ny)            :: wc2
    complex(kind=DP), dimension(-nx:nx)                 :: zf
    real(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1)      :: wr3
    real(kind=DP) :: entrpy_nz, entrpy_zf, fenegy_nz, fenegy_zf, menegy_nz, menegy_zf, &
                     peint_nz, peint_zf, pmint_nz, pmint_zf, neint_nz, neint_zf, nmint_nz, nmint_zf, &
                     dcd_nz, dcd_zf, dgp_es, dgp_em, dqp_es, dqp_em
    real(kind=DP) :: entrpy_wk, fenegy_wk, menegy_wk, peint_wk, pmint_wk, neint_wk, nmint_wk, dcd_wk, dgp_wk, dqp_wk
    integer :: mx, my, iz, iv, im


      allocate( wf(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) )
      allocate(  dens(-nx:nx,0:ny,-nz:nz-1) )
      allocate( upara(-nx:nx,0:ny,-nz:nz-1) )
      allocate(  pres(-nx:nx,0:ny,-nz:nz-1) )
      allocate( qpara(-nx:nx,0:ny,-nz:nz-1) )
      allocate(    ni(-nx:nx,0:ny,-nz:nz-1) )
      allocate(   wc3(-nx:nx,0:ny,-nz:nz-1) )

!$OMP parallel workshare
      wc3(:,:,:) = ( 0._DP, 0._DP )
      wr3(:,:,:) =   0._DP
      wc2(:,:)   = ( 0._DP, 0._DP )
!$OMP end parallel workshare

!$OMP parallel
      do im = 0, nm
!$OMP do
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                wf(mx,my,iz,iv,im) = ff(mx,my,iz,iv,im) * j0(mx,my,iz,im)
              end do
            end do
          end do
        end do
!$OMP end do nowait
      end do
!$OMP end parallel

      call intgrl_v0_moment( wf, dens )
      call intgrl_v2_moment( wf, pres )
     
!$OMP parallel
      do im = 0, nm
!$OMP do
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                wf(mx,my,iz,iv,im) = vl(iv) * ff(mx,my,iz,iv,im) * j0(mx,my,iz,im)
              end do
            end do
          end do
        end do
!$OMP end do nowait
      end do
!$OMP end parallel

      call intgrl_v0_moment( wf, upara )
      call intgrl_v2_moment( wf, qpara )

!$OMP parallel do
      do iz = -nz, nz-1
        do my = ist_y, iend_y
          do mx = -nx, nx
             dens(mx,my,iz) = fcs(ranks) / Znum(ranks) * dens(mx,my,iz)
            upara(mx,my,iz) = fcs(ranks) / Znum(ranks)  &
                            * sqrt( tau(ranks) / Anum(ranks) ) * upara(mx,my,iz)
             pres(mx,my,iz) = fcs(ranks) / Znum(ranks) * tau(ranks) * pres(mx,my,iz)
            qpara(mx,my,iz) = fcs(ranks) / Znum(ranks) * tau(ranks)  &
                            * sqrt( tau(ranks) / Anum(ranks) ) * qpara(mx,my,iz)
          end do
        end do
      end do


! --- \begin{entropy_calculation}

      entrpy(:,:) = 0._DP
      entrpy_nz   = 0._DP
      entrpy_zf   = 0._DP

!$OMP parallel
      do im = 0, nm
!$OMP do
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                wf(mx,my,iz,iv,im) =                                     &
                        ff(mx,my,iz,iv,im) * conjg( ff(mx,my,iz,iv,im) ) &
                        / fmx(iz,iv,im)
              end do
            end do
          end do
        end do
!$OMP end do nowait
      end do
!$OMP end parallel

      call intgrl_v0_moment( wf, wc3 )

      call intgrl_thet ( wc3, wc2 )

!$OMP parallel do
      do my = ist_y, iend_y
        do mx = -nx, nx
          entrpy(mx,my) = 0.5_DP * fcs(ranks) * tau(ranks) / Znum(ranks) * real( wc2(mx,my), kind=DP )
        end do
      end do

!$OMP parallel do reduction(+:entrpy_nz)
      do my = ist1_y, iend_y
        do mx = -nx, nx
          entrpy_nz = entrpy_nz + entrpy(mx,my)
        end do
      end do
      entrpy_wk  =   0._DP
      call MPI_Allreduce( entrpy_nz, entrpy_wk, 1, MPI_DOUBLE_PRECISION, &
                          MPI_SUM, fft_comm_world, ierr_mpi )
      entrpy_nz = 2._DP * entrpy_wk

      if( rankw == 0 )  then
        my   = 0
          do mx = 0, nx
            entrpy_zf = entrpy_zf + entrpy(mx,my)
          end do
      endif
      entrpy_zf = 2._DP * entrpy_zf

! --- \end{entropy_calculation}


! --- \begin{electric_energy_calculation}

      fenegy_nz   = 0._DP
      fenegy_zf   = 0._DP

      if ( ns == 1 ) then    ! --- for ITG-ae, ETG-ai
!$OMP parallel do
        do iz = -nz, nz-1
          do my = ist_y, iend_y
            do mx = -nx, nx
              wr3(mx,my,iz) = 0.5_DP * real( phi(mx,my,iz) * conjg( phi(mx,my,iz) )  &
                             , kind=DP ) * ( 1._DP - g0(mx,my,iz) + tau(0)*tau_ad )
            end do
          end do
        end do
      else                   ! --- for multi-species
!$OMP parallel do
        do iz = -nz, nz-1
          do my = ist_y, iend_y
            do mx = -nx, nx
              wr3(mx,my,iz) = 0.5_DP * fct_e_energy(mx,my,iz)  &
                            * real( phi(mx,my,iz) * conjg( phi(mx,my,iz) ), kind=DP )
            end do
          end do
        end do
      end if

      call intgrl_thet ( wr3, fenegy )

!$OMP parallel do reduction(+:fenegy_nz)
      do my = ist1_y, iend_y
        do mx = -nx, nx
          fenegy_nz   = fenegy_nz + fenegy(mx,my)
        end do
      end do
      fenegy_wk  =   0._DP
      call MPI_Allreduce( fenegy_nz, fenegy_wk, 1, MPI_DOUBLE_PRECISION, &
                          MPI_SUM, fft_comm_world, ierr_mpi )
      fenegy_nz = 2._DP * fenegy_wk

      if( rankw == 0 )  then
        my   = 0
          do mx = 0, nx
            fenegy_zf = fenegy_zf + fenegy(mx,my)
          end do
      endif

      if ( ns == 1 .and. sgn(0) > 0._DP ) then ! --- for ITG-ae

!$OMP parallel
        do im = 0, nm
!$OMP do
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  wf(mx,my,iz,iv,im) = ff(mx,my,iz,iv,im) * j0(mx,my,iz,im)
                end do
              end do
            end do
          end do
!$OMP end do nowait
        end do
!$OMP end parallel
  
        call intgrl_v0_moment( wf, ni )
  
        if( rankw == 0 )  then
          do iz = -nz, nz-1
            ni(0,0,iz) = ( 0._DP, 0._DP )
          end do
        endif
  
        zf = ( 0._DP, 0._DP )

        if( rankw == 0 ) then

          my   = 0
!$OMP parallel do
            do iz = -nz, nz-1
              do mx = -nx, -1
                wc3(mx,my,iz)   = ni(mx,my,iz)                                    &
                                 / ( (  1._DP - g0(mx,my,iz) + tau(0)*tau_ad ) * fctgt(mx) )
              end do
              mx = 0
                wc3(mx,my,iz) = (0._DP, 0._DP)
              do mx = 1, nx
                wc3(mx,my,iz)   = ni(mx,my,iz)                                    &
                                 / ( (  1._DP - g0(mx,my,iz) + tau(0)*tau_ad ) * fctgt(mx) )
              end do
            end do
  
          call intgrl_fsrf ( wc3, zf )
  
          zf(0)   = ( 0._DP, 0._DP )

          do mx = 0, nx
            fenegy_zf   = fenegy_zf - 0.5_DP * real( zf(mx) * conjg( zf(mx) ), kind=DP) * tau(0)*tau_ad
          end do

        end if

      end if
      fenegy_zf = 2._DP * fenegy_zf

! --- \end{electric_energy_calculation}


! --- \begin{magnetic_energy_calculation}

      menegy_nz   = 0._DP
      menegy_zf   = 0._DP

      if ( beta .ne. 0._DP ) then

!$OMP parallel do
        do iz = -nz, nz-1
          do my = ist_y, iend_y
            do mx = -nx, nx
              wr3(mx,my,iz) = 0.5_DP * fct_m_energy(mx,my,iz)  &
                            * real( Al(mx,my,iz) * conjg( Al(mx,my,iz) ), kind=DP )
            end do
          end do
        end do
  
        call intgrl_thet ( wr3, menegy )
  
!$OMP parallel do reduction(+:menegy_nz)
        do my = ist1_y, iend_y
          do mx = -nx, nx
            menegy_nz = menegy_nz + menegy(mx,my)
          end do
        end do
        menegy_wk  =   0._DP
        call MPI_Allreduce( menegy_nz, menegy_wk, 1, MPI_DOUBLE_PRECISION, &
                            MPI_SUM, fft_comm_world, ierr_mpi )
        menegy_nz = 2._DP * menegy_wk
  
        if( rankw == 0 )  then
          my   = 0
            do mx = 0, nx
              menegy_zf = menegy_zf + menegy(mx,my)
            end do
        endif
        menegy_zf = 2._DP * menegy_zf

      else
        menegy(:,:) = 0._DP
      end if

! --- \end{magnetic_energy_calculation}


! --- \begin{particle_field_interaction}

      peint(:,:) = 0._DP
      peint_nz   = 0._DP
      peint_zf   = 0._DP

!$OMP parallel
      do im = 0, nm
!$OMP do
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                wf(mx,my,iz,iv,im) = dh(mx,my,iz,iv,im) * j0(mx,my,iz,im) * sgn(ranks) * fcs(ranks)
              end do
            end do
          end do
        end do
!$OMP end do nowait
      end do
!$OMP end parallel

      call intgrl_v0_moment( wf, wc3 )

!$OMP parallel do
      do iz = -nz, nz-1
        do my = ist_y, iend_y
          do mx = -nx, nx
            wr3(mx,my,iz) = - real( conjg( phi(mx,my,iz) ) * wc3(mx,my,iz), kind=DP )
          end do
        end do
      end do

      call intgrl_thet ( wr3, peint )

!$OMP parallel do reduction(+:peint_nz)
      do my = ist1_y, iend_y
        do mx = -nx, nx
          peint_nz = peint_nz + peint(mx,my)
        end do
      end do
      peint_wk  =   0._DP
      call MPI_Allreduce( peint_nz, peint_wk, 1, MPI_DOUBLE_PRECISION, &
                          MPI_SUM, fft_comm_world, ierr_mpi )
      peint_nz = 2._DP * peint_wk

      if( rankw == 0 )  then
        my   = 0
          do mx = 0, nx
            peint_zf = peint_zf + peint(mx,my)
          end do
      endif
      peint_zf = 2._DP * peint_zf


      pmint(:,:) = 0._DP
      pmint_nz   = 0._DP
      pmint_zf   = 0._DP

      if ( beta .ne. 0._DP ) then

        call fld_emfield_hh ( dh, wc3 )
  
!$OMP parallel do
        do iz = -nz, nz-1
          do my = ist_y, iend_y
            do mx = -nx, nx
              wr3(mx,my,iz) = - real( conjg( upara(mx,my,iz) )  &
                    * sgn(ranks) * Znum(ranks) * wc3(mx,my,iz), kind=DP )
            end do
          end do
        end do
  
        call intgrl_thet ( wr3, pmint )
  
!$OMP parallel do reduction(+:pmint_nz)
        do my = ist1_y, iend_y
          do mx = -nx, nx
            pmint_nz = pmint_nz + pmint(mx,my)
          end do
        end do
        pmint_wk  =   0._DP
        call MPI_Allreduce( pmint_nz, pmint_wk, 1, MPI_DOUBLE_PRECISION, &
                            MPI_SUM, fft_comm_world, ierr_mpi )
        pmint_nz = 2._DP * pmint_wk
  
        if( rankw == 0 )  then
          my   = 0
            do mx = 0, nx
              pmint_zf = pmint_zf + pmint(mx,my)
            end do
        endif
        pmint_zf = 2._DP * pmint_zf

      end if

! --- \begin{particle_field_interaction}


! --- \begin{nonlinear_interaction}

      neint(:,:) = 0._DP
      nmint(:,:) = 0._DP
      neint_nz   = 0._DP
      neint_zf   = 0._DP
      nmint_nz   = 0._DP
      nmint_zf   = 0._DP

      if( trim(calc_type) == "nonlinear" ) then

        call trans_sum ( ff, phi, Al, neint, nmint )

!$OMP parallel do reduction(+:neint_nz)
        do my = ist1_y, iend_y
          do mx = -nx, nx
            neint_nz = neint_nz + neint(mx,my)
          end do
        end do
        neint_wk = 0._DP
        call MPI_Allreduce( neint_nz, neint_wk, 1, MPI_DOUBLE_PRECISION, &
                              MPI_SUM, fft_comm_world, ierr_mpi )
        neint_nz = 2._DP * neint_wk
  
        if( rankw == 0 )  then
          my   = 0
!$OMP parallel do reduction(+:neint_zf)
            do mx = 0, nx
              neint_zf = neint_zf + neint(mx,my)
            end do
        endif
        neint_zf = 2._DP * neint_zf
  
!$OMP parallel do reduction(+:nmint_nz)
        do my = ist1_y, iend_y
          do mx = -nx, nx
            nmint_nz = nmint_nz + nmint(mx,my)
          end do
        end do
        nmint_wk = 0._DP
        call MPI_Allreduce( nmint_nz, nmint_wk, 1, MPI_DOUBLE_PRECISION, &
                              MPI_SUM, fft_comm_world, ierr_mpi )
        nmint_nz = 2._DP * nmint_wk
  
        if( rankw == 0 )  then
          my   = 0
            do mx = 0, nx
              nmint_zf = nmint_zf + nmint(mx,my)
            end do
        endif
        nmint_zf = 2._DP * nmint_zf

      end if

! --- \end{nonlinear_interaction}


! --- \begin{flux_calculation}

      call calc_flux ( dens, phi, pflux_es )
      dgp_es = 0._DP
!$OMP parallel do reduction(+:dgp_es)
      do my = ist1_y, iend_y
        do mx = -nx, nx
          dgp_es = dgp_es + tau(ranks) * (R0_Ln(ranks) + R0_Lt(ranks)) * pflux_es(mx,my)
        end do
      end do
      if( rankw == 0 )   then
        my = 0
          do mx = 0, nx
            dgp_es = dgp_es + tau(ranks) * (R0_Ln(ranks) + R0_Lt(ranks)) * pflux_es(mx,my)
          end do
      endif
      dgp_wk = 0._DP
      call MPI_Allreduce( dgp_es, dgp_wk, 1, MPI_DOUBLE_PRECISION, &
                          MPI_SUM, fft_comm_world, ierr_mpi )
      dgp_es = 2._DP * dgp_wk

      call calc_flux ( upara, -Al, pflux_em )
      dgp_em = 0._DP
!$OMP parallel do reduction(+:dgp_em)
      do my = ist1_y, iend_y
        do mx = -nx, nx
          dgp_em = dgp_em + tau(ranks) * (R0_Ln(ranks) + R0_Lt(ranks)) * pflux_em(mx,my)
        end do
      end do
      if( rankw == 0 )   then
        my = 0
          do mx = 0, nx
            dgp_em = dgp_em + tau(ranks) * (R0_Ln(ranks) + R0_Lt(ranks)) * pflux_em(mx,my)
          end do
      endif
      dgp_wk = 0._DP
      call MPI_Allreduce( dgp_em, dgp_wk, 1, MPI_DOUBLE_PRECISION, &
                          MPI_SUM, fft_comm_world, ierr_mpi )
      dgp_em = 2._DP * dgp_wk

      call calc_flux ( pres, phi, eflux_es )
      dqp_es = 0._DP
!$OMP parallel do reduction(+:dqp_es)
      do my = ist1_y, iend_y
        do mx = -nx, nx
          dqp_es = dqp_es + R0_Lt(ranks) * ( eflux_es(mx,my) - 2.5_DP * tau(ranks) * pflux_es(mx,my) )
        end do
      end do
      if( rankw == 0 )   then
        my = 0
          do mx = 0, nx
            dqp_es = dqp_es + R0_Lt(ranks) * ( eflux_es(mx,my) - 2.5_DP * tau(ranks) * pflux_es(mx,my) )
          end do
      endif
      dqp_wk = 0._DP
      call MPI_Allreduce( dqp_es, dqp_wk, 1, MPI_DOUBLE_PRECISION, &
                          MPI_SUM, fft_comm_world, ierr_mpi )
      dqp_es = 2._DP * dqp_wk

      call calc_flux ( qpara, -Al, eflux_em )
      dqp_em = 0._DP
!$OMP parallel do reduction(+:dqp_em)
      do my = ist1_y, iend_y
        do mx = -nx, nx
          dqp_em = dqp_em + R0_Lt(ranks) * ( eflux_em(mx,my) - 2.5_DP * tau(ranks) * pflux_em(mx,my) )
        end do
      end do
      if( rankw == 0 )   then
        my = 0
          do mx = 0, nx
            dqp_em = dqp_em + R0_Lt(ranks) * ( eflux_em(mx,my) - 2.5_DP * tau(ranks) * pflux_em(mx,my) )
          end do
      endif
      dqp_wk = 0._DP
      call MPI_Allreduce( dqp_em, dqp_wk, 1, MPI_DOUBLE_PRECISION, &
                          MPI_SUM, fft_comm_world, ierr_mpi )
      dqp_em = 2._DP * dqp_wk

! --- \end{flux_calculation}


! --- \begin{collisionaldissipation_calculation}

      dcd(:,:) = 0._DP
      dcd_nz   = 0._DP
      dcd_zf   = 0._DP

!$OMP parallel
      do im = 0, nm
!$OMP do
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                wf(mx,my,iz,iv,im)   = fcs(ranks) / Znum(ranks) * cf(mx,my,iz,iv,im)  &
                      * ( tau(ranks) * conjg( ff(mx,my,iz,iv,im) ) / fmx(iz,iv,im)    &
                        + sgn(ranks) * Znum(ranks) * j0(mx,my,iz,im) * conjg( phi(mx,my,iz) ) )
              end do
            end do
          end do
        end do
!$OMP end do nowait
      end do
!$OMP end parallel

      call intgrl_v0_moment ( wf, wc3 )

      call intgrl_thet ( wc3, wc2 )

!$OMP parallel do
      do my = ist_y, iend_y
        do mx = -nx, nx
          dcd(mx,my) = real( wc2(mx,my), kind=DP )
        end do
      end do

!$OMP parallel do reduction(+:dcd_nz)
      do my = ist1_y, iend_y
        do mx = -nx, nx
          dcd_nz = dcd_nz + dcd(mx,my)
        end do
      end do
      dcd_wk = 0._DP
      call MPI_Allreduce( dcd_nz, dcd_wk, 1, MPI_DOUBLE_PRECISION, &
                            MPI_SUM, fft_comm_world, ierr_mpi )
      dcd_nz = 2._DP * dcd_wk

      if( rankw == 0 )  then
        my   = 0
          do mx = 0, nx
            dcd_zf = dcd_zf + dcd(mx,my)
          end do
      endif
      dcd_zf = 2._DP * dcd_zf

! --- \end{collisionaldissipation_calculation}


      if ( rank == 0 ) then
        write( unit=obln, fmt="(f15.8, SP, 256ES24.15e3)" ) &
             time,      & ! [ 1] Time
             entrpy_nz, & ! [ 2] Entropy S_s (ky/=0)
             entrpy_zf, & ! [ 3] Entropy S_s (ky==0)
             fenegy_nz, & ! [ 4] Electrostatic field energy W_E (ky/=0)
             fenegy_zf, & ! [ 5] Electrostatic field energy W_E (ky/=0)
             menegy_nz, & ! [ 6] Magnetic field energy W_M (ky/=0)
             menegy_zf, & ! [ 7] Magnetic field energy W_M (ky/=0)
             peint_nz,  & ! [ 8] W_E to S_s interaction R_sE (ky/=0)
             peint_zf,  & ! [ 9] W_E to S_s interaction R_sE (ky==0)
             pmint_nz,  & ! [10] W_M to S_s interaction R_sM (ky/=0)
             pmint_zf,  & ! [11] W_M to S_s interaction R_sM (ky==0)
             neint_nz,  & ! [12] S_s(ky==0) to S_s(ky/=0) entropy transfer via ExB nonlinearity -I_sE^(z) (ky/=0)
             neint_zf,  & ! [13] S_s(ky/=0) to S_s(ky==0) entropy transfer via ExB nonlinearity  I_sE^(z) (ky==0)
             nmint_nz,  & ! [14] S_s(ky==0) to S_s(ky/=0) entropy transfer via magnetic nonlinearity -I_sM^(z) (ky/=0)
             nmint_zf,  & ! [15] S_s(ky==0) to S_s(ky/=0) entropy transfer via magnetic nonlinearity  I_sM^(z) (ky==0)
             dcd_nz,    & ! [16] Collisional dissipation D_s (ky/=0)
             dcd_zf,    & ! [17] Collisional dissipation D_s (ky==0)
             dgp_es,    & ! [18] Particle flux term by ExB flows T_s*G_sE/L_ps
             dgp_em,    & ! [19] Particle flux term by magnetic flutters T_s*G_sM/L_ps
             dqp_es,    & ! [20] Heat flux term by ExB flows Theta_sE/L_Ts
             dqp_em       ! [21] Heat flux term by magnetic flutters Theta_sM/L_Ts
      end if

      deallocate( wf )
      deallocate(  dens )
      deallocate( upara )
      deallocate(  pres )
      deallocate( qpara )
      deallocate(    ni )
      deallocate(   wc3 )

! --- divergence trap
        if( trim(calc_type) == "nonlinear" ) then
          if ( entrpy_nz > 1.d100 ) then 
            write(olog,*) "DIVERGE!! at ", time, entrpy_nz
            call tips_flush
            call MPI_Abort(MPI_COMM_WORLD, ierr_mpi)
            stop "STOPPED at DIVERGENCE trap!"
          end if
        end if
! --- 


  END SUBROUTINE balance


!--------------------------------------
  SUBROUTINE calc_flux ( mom, phi, flux )
!--------------------------------------
!     Calculate turbulent flux

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: mom, phi
    real(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny) :: flux

    real(kind=DP), dimension(:,:,:), allocatable :: wr3
    integer  ::  mx, my, iz

      allocate( wr3(-nx:nx,0:ny,-nz:nz-1) )

!$OMP parallel do
      do iz = -nz, nz-1
        do my = ist_y, iend_y
          do mx = -nx, nx
            wr3(mx,my,iz) = real( - ui * ky(my) * phi(mx,my,iz)  &
                                * conjg( mom(mx,my,iz) ), kind=DP )
          end do
        end do
      end do

      call intgrl_thet ( wr3, flux )

      deallocate( wr3 )

  END SUBROUTINE calc_flux


!--------------------------------------
  SUBROUTINE write_moments ( ff, time )
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    real(kind=DP), intent(in) :: time

    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: wf
    complex(kind=DP), dimension(:,:,:), allocatable :: dens, upara, ppara, pperp, qlpara, qlperp
    integer :: mx, my, iz, iv, im

      allocate( wf(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) )
      allocate(   dens(-nx:nx,0:ny,-nz:nz-1) )
      allocate(  upara(-nx:nx,0:ny,-nz:nz-1) )
      allocate(  ppara(-nx:nx,0:ny,-nz:nz-1) )
      allocate(  pperp(-nx:nx,0:ny,-nz:nz-1) )
      allocate( qlpara(-nx:nx,0:ny,-nz:nz-1) )
      allocate( qlperp(-nx:nx,0:ny,-nz:nz-1) )

!$OMP parallel do collapse(2) private(mx,my,iz,iv,im)
      do im = 0, nm
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                wf(mx,my,iz,iv,im) = ff(mx,my,iz,iv,im) * j0(mx,my,iz,im)
              end do
            end do
          end do
        end do
      end do

      call intgrl_v0_moment( wf, dens )

!$OMP parallel do collapse(2) private(mx,my,iz,iv,im)
      do im = 0, nm
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                wf(mx,my,iz,iv,im) = vl(iv) * ff(mx,my,iz,iv,im) * j0(mx,my,iz,im)
              end do
            end do
          end do
        end do
      end do

      call intgrl_v0_moment( wf, upara )

!$OMP parallel do collapse(2) private(mx,my,iz,iv,im)
      do im = 0, nm
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                wf(mx,my,iz,iv,im) = 0.5_DP * vl(iv)**2 * ff(mx,my,iz,iv,im) * j0(mx,my,iz,im)
              end do
            end do
          end do
        end do
      end do

      call intgrl_v0_moment( wf, ppara )

!$OMP parallel do collapse(2) private(mx,my,iz,iv,im)
      do im = 0, nm
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                wf(mx,my,iz,iv,im) = mu(im) * omg(iz) * ff(mx,my,iz,iv,im) * j0(mx,my,iz,im)
              end do
            end do
          end do
        end do
      end do

      call intgrl_v0_moment( wf, pperp )

!$OMP parallel do collapse(2) private(mx,my,iz,iv,im)
      do im = 0, nm
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                wf(mx,my,iz,iv,im) = 0.5_DP * vl(iv)**3 * ff(mx,my,iz,iv,im) * j0(mx,my,iz,im)
              end do
            end do
          end do
        end do
      end do

      call intgrl_v0_moment( wf, qlpara )

!$OMP parallel do collapse(2) private(mx,my,iz,iv,im)
      do im = 0, nm
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                wf(mx,my,iz,iv,im) = vl(iv) * mu(im) * omg(iz) * ff(mx,my,iz,iv,im) * j0(mx,my,iz,im)
              end do
            end do
          end do
        end do
      end do

      call intgrl_v0_moment( wf, qlperp )

!$OMP parallel do
      do iz = -nz, nz-1
        do my = ist_y, iend_y
          do mx = -nx, nx
              dens(mx,my,iz) = fcs(ranks) / Znum(ranks) * dens(mx,my,iz)
             upara(mx,my,iz) = fcs(ranks) / Znum(ranks)  &
                             * sqrt( tau(ranks) / Anum(ranks) ) * upara(mx,my,iz)
             ppara(mx,my,iz) = fcs(ranks) / Znum(ranks) * tau(ranks) * ppara(mx,my,iz)
             pperp(mx,my,iz) = fcs(ranks) / Znum(ranks) * tau(ranks) * pperp(mx,my,iz)
            qlpara(mx,my,iz) = fcs(ranks) / Znum(ranks) * tau(ranks)  &
                             * sqrt( tau(ranks) / Anum(ranks) ) * qlpara(mx,my,iz)
            qlperp(mx,my,iz) = fcs(ranks) / Znum(ranks) * tau(ranks)  &
                             * sqrt( tau(ranks) / Anum(ranks) ) * qlperp(mx,my,iz)
          end do
        end do
      end do

      !fj start 202011
      !if ( vel_rank == 0 ) then
        !write( unit=omom ) time, dens, upara, ppara, pperp, qlpara, qlperp
      !end if
      call fileio_write_mom( dens, upara, ppara, &
                             pperp, qlpara, qlperp, time )
      !fj end 202011

      deallocate( wf )
      deallocate(   dens )
      deallocate(  upara )
      deallocate(  ppara )
      deallocate(  pperp )
      deallocate( qlpara )
      deallocate( qlperp )

  END SUBROUTINE write_moments


END MODULE GKV_out
