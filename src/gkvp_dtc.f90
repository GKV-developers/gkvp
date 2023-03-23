MODULE GKV_dtc
!-------------------------------------------------------------------------------
!
!    Time step size control
!
!    Update history of gkvp_dtc.f90
!    --------------
!      gkvp_f0.57 (S. Maeyama, Oct 2020)
!        - Version number f0.57 is removed from filename.
!        - Minimum value of time step size is set to avoid zero-division
!          e.g., nu_max2=max(nu_max2,1.d-20) for nu=0.
!
!-------------------------------------------------------------------------------

  use GKV_header
  use GKV_mpienv
  use GKV_exb, only: exb_maxvx_eachrank, exb_maxvy_eachrank
  use GKV_colliimp, only: colliimp_colli, gvl, gvp, gnu_ds, gnu_ps, gnu_hs, gnu_gs

  implicit none

  private

  real(kind=DP), save :: dt_linear, dt_nl, dt_limit, dt_col

  real(kind=DP), save :: dx_inv, dy_inv

  integer, save :: flag_time_advnc ! 0: rkg4
                                   ! 1: imp_colli
  integer, save :: flag_time_split ! 0: colli(dt/2)->RK(dt) ! for imp_colli
                                   ! 1: colli(dt  )->RK(dt) ! for imp_colli

  public   dtc_init, dtc_cntrl, flag_time_advnc, flag_time_split


 CONTAINS


!--------------------------------------
  SUBROUTINE dtc_init( lx, ly, vmax )
!--------------------------------------

    real(kind=DP), intent(in) :: lx, ly, vmax

    real(kind=DP) :: dt_perp, dt_zz, dt_vl
    real(kind=DP) :: kvd_max, kvd_max2, vl_max, vl_max2, mir_max, mir_max2
    real(kind=DP) :: ksq_max0, ksq_max, nu_max, nu_max2, nu_temp
    real(kind=DP) :: cs, dx, dy, kvd
    integer :: mx, my, iz, iv, im, is


      ksq_max0 = 0._DP
      ksq_max  = 0._DP
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
!              do mx = 0, nx
                if ( ksq_max0 < ksq(mx,my,iz) ) ksq_max0 = ksq(mx,my,iz)
              end do
            end do
          end do
      call MPI_Allreduce( ksq_max0, ksq_max, 1, MPI_DOUBLE_PRECISION, &
                          MPI_MAX, MPI_COMM_WORLD, ierr_mpi )

      kvd_max = 0._DP
      do im = 0, nm
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
!              do mx = 0, nx
                kvd = kx(mx) * vdx(iz,iv,im) + ky(my) * vdy(iz,iv,im)
                if ( kvd_max < abs(kvd) ) then 
                  kvd_max = abs(kvd)
                end if 
              end do
            end do
          end do
        end do
      end do
      call MPI_Allreduce( kvd_max, kvd_max2, 1, MPI_DOUBLE_PRECISION, &
                          MPI_MAX, MPI_COMM_WORLD, ierr_mpi )
      kvd_max2 = max(kvd_max2, 1.d-20)
      dt_perp = courant_num * pi / kvd_max2

      cs = sqrt( tau(ranks) / Anum(ranks) )
      vl_max = 0._DP
      do iz = -nz, nz-1
        if ( vl_max < cs * vmax / dpara(iz) ) vl_max = cs * vmax / dpara(iz)
      end do
      call MPI_Allreduce( vl_max, vl_max2, 1, MPI_DOUBLE_PRECISION, &
                          MPI_MAX, MPI_COMM_WORLD, ierr_mpi )
      vl_max2 = max(vl_max2, 1.d-20)
      dt_zz = courant_num / vl_max2

      mir_max = 0._DP
      do im = 0, nm
        do iz = -nz, nz-1
          if ( mir_max < cs * mir(iz,im) ) mir_max = cs * mir(iz,im)
        end do
      end do
      call MPI_Allreduce( mir_max, mir_max2, 1, MPI_DOUBLE_PRECISION, &
                          MPI_MAX, MPI_COMM_WORLD, ierr_mpi )
      mir_max2 = max(mir_max2, 1.d-20)
      dt_vl = courant_num * dv / mir_max2

      if ( trim(col_type) == "LB" ) then

        nu_max = nu(ranks)*3._DP*dsqrt(pi)*ctauiv(ranks,ranks)/4._DP &
               * ( 2._DP / dv**2 )
        do iz = -nz, nz-1
          nu_temp = nu(ranks)*3._DP*dsqrt(pi)*ctauiv(ranks,ranks)/4._DP &
                  * ( 2._DP / dvp(iz)**2 )
          if ( nu_max < nu_temp ) nu_max = nu_temp
        end do
        call MPI_Allreduce( nu_max, nu_max2, 1, MPI_DOUBLE_PRECISION, &
                            MPI_MAX, MPI_COMM_WORLD, ierr_mpi )
        nu_max2 = max(nu_max2, 1.d-20)
        dt_col = courant_num / nu_max2
        
      else if ( trim(col_type) == "full" .or. &
                trim(col_type) == "lorentz" ) then

        !nu_max = 0._DP
        !do im = 0, nm
        !  do iv = 1, 2*nv
        !    do iz = -nz, nz-1
        !      nu_temp = ( nu_ps(iz,iv,im) * vl(iv)**2    & 
        !                + nu_ds(iz,iv,im) * vp(iz,im)**2 &
        !                ) * 0.5_DP * ( 2._DP / dv**2 )
        !      if ( nu_max < nu_temp ) nu_max = nu_temp
        !    end do
        !  end do
        !end do
        !do im = 0, nm
        !  do iv = 1, 2*nv
        !    do iz = -nz, nz-1
        !      nu_temp = ( nu_ds(iz,iv,im) * vl(iv)**2    & 
        !                + nu_ps(iz,iv,im) * vp(iz,im)**2 &
        !                ) * 0.5_DP * ( 2._DP / dvp(iz)**2 )
        !      if ( nu_max < nu_temp ) nu_max = nu_temp
        !    end do
        !  end do
        !end do
        nu_max = 0._DP
        do iz = -nz, nz-1
          do is = 0, ns-1
            do im = 0, global_nm
              do iv = 1, 2*global_nv
                nu_temp = ( gnu_ps(iv,im,is,iz) * gvl(iv)**2    & 
                          + gnu_ds(iv,im,is,iz) * gvp(im,iz)**2 &
                          ) * 0.5_DP * ( 2._DP / dv**2 )
                if ( nu_max < nu_temp ) nu_max = nu_temp
              end do
            end do
          end do
        end do
        do iz = -nz, nz-1
          do is = 0, ns-1
            do im = 0, global_nm
              do iv = 1, 2*global_nv
                nu_temp = ( gnu_ds(iv,im,is,iz) * gvl(iv)**2    & 
                          + gnu_ps(iv,im,is,iz) * gvp(im,iz)**2 &
                          ) * 0.5_DP * ( 2._DP / dvp(iz)**2 )
                if ( nu_max < nu_temp ) nu_max = nu_temp
              end do
            end do
          end do
        end do
        call MPI_Allreduce( nu_max, nu_max2, 1, MPI_DOUBLE_PRECISION, &
                            MPI_MAX, MPI_COMM_WORLD, ierr_mpi )
        nu_max2 = max(nu_max2, 1.d-20)
        dt_col = courant_num / nu_max2

      else

        write( olog, * ) "This col_type is not supported by dtc:", trim(col_type)
        dt_col = 99._DP

      end if


      dt_linear = min( dt_perp, dt_zz, dt_vl )

      if (trim(time_advnc) == "rkg4") then

        dt_limit = min( dt_max, dt_linear, dt_col )
        dt = dt_max
        if ( adapt_dt ) then
          if ( dt < 0.8_DP * dt_limit .or. 1.2_DP * dt_limit < dt ) dt = dt_limit
        end if
        flag_time_advnc = 0

      else if (trim(time_advnc) == "imp_colli") then

        dt_limit = min( dt_max, dt_linear )
        dt = dt_max
        if ( adapt_dt ) then
          if ( dt < 0.8_DP * dt_limit .or. 1.2_DP * dt_limit < dt ) dt = dt_limit
        end if
        flag_time_advnc = 1

      else if (trim(time_advnc) == "auto_init") then

        dt_limit = min( dt_max, dt_linear )
        dt = dt_max
        if ( adapt_dt ) then
          if ( dt < 0.8_DP * dt_limit .or. 1.2_DP * dt_limit < dt ) dt = dt_limit
        end if
        if (dt < dt_col) then
           flag_time_advnc = 0
        else
           flag_time_advnc = 1
        end if

      else

        write( olog, * ) " # Wrong type of time_advnc in namelist: ", time_advnc
        call flush(olog)
        call MPI_Finalize ( ierr_mpi )
        stop

      end if

      flag_time_split = 0

        write( olog, * ) " # Time step size control"
        write( olog, * ) ""
        write( olog, * ) " # time_advnc   = ", time_advnc
        write( olog, * ) " # flag_time_adv= ", flag_time_advnc
        write( olog, * ) " # courant num. = ", courant_num
        write( olog, * ) " # ksq_max      = ", ksq_max
        write( olog, * ) " # dt_perp      = ", dt_perp
        write( olog, * ) " # dt_zz        = ", dt_zz
        write( olog, * ) " # dt_vl        = ", dt_vl
        write( olog, * ) " # dt_col       = ", dt_col
        write( olog, * ) " # dt_linear    = ", dt_linear
        write( olog, * ) " # dt_max       = ", dt_max
        write( olog, * ) " # dt           = ", dt
        write( olog, * ) ""

      dx = lx / real( nxw, kind=DP )
      dy = ly / real( nyw, kind=DP )
      dx_inv = 1._DP / dx
      dy_inv = 1._DP / dy

  
  END SUBROUTINE dtc_init


!--------------------------------------
  SUBROUTINE dtc_cntrl( ff, phi, Al, hh, time, id )
!--------------------------------------

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: phi, Al
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: hh
    real(kind=DP), intent(in) :: time
    integer, intent(in) :: id

    real(kind=DP), save :: tout_dtc


      if( id == 0 ) then

        tout_dtc  = ( int( ( time + eps )/dtout_dtc ) + 1 ) * dtout_dtc
 
        if ( trim(calc_type) == "nonlinear" ) call dtc_estimate

        if ( dt < 0.8_DP * dt_limit .or. 1.2_DP * dt_limit < dt ) then
          dt = dt_limit
          write( olog, * ) &
            " # dt is changed at time = ", time, ", dt = ", dt
        end if
     
        if ( rankg == 0 ) then
          write( unit=odtc, fmt="(f10.5, 1p, 3e15.7)" )  &
            time, dt, dt_limit, dt_nl
        end if
  
      else if( id == 1 ) then

        if ( time >= tout_dtc - eps ) then

          tout_dtc   = tout_dtc + dtout_dtc
   
          if ( trim(calc_type) == "nonlinear" ) call dtc_estimate
  
          if ( dt < 0.8_DP * dt_limit .or. 1.2_DP * dt_limit < dt ) then

            if (flag_time_advnc == 1) then ! Operator split + implicit collision
              if (flag_time_split == 1) then ! dt/2 integration for 2nd-order split
                call colliimp_colli(0.5_DP*dt, ff, phi, al, hh)
                                               !%%% For debug %%%
                                               ! if (rankg==0) write(olog,*) &
                                               ! time, "half-step for dtc_cntrl"
                                               !%%%%%%%%%%%%%%%%%
              end if
              flag_time_split = 0 ! flag_time_split==0 means you have physical quantities at time
            end if

            dt = dt_limit
            write( olog, * ) &
              " # dt is changed at time = ", time, ", dt = ", dt
          end if
       
          if ( rankg == 0 ) then
            write( unit=odtc, fmt="(f10.5, 1p, 3e15.7)" )  &
              time, dt, dt_limit, dt_nl
          end if

        end if

      end if


  END SUBROUTINE dtc_cntrl


!!--------------------------------------
!  SUBROUTINE dtc_estimate
!!--------------------------------------
!
!    real(kind=DP) :: wx_nl1, wx_nl2, wy_nl1, wy_nl2, w_nl_max, w_nl_max2
!    real(kind=DP) :: cs
!    integer  ::  mx, my, iz, iv, im
!  
!  
!      w_nl_max = 0._DP
!  
!      cs = sqrt( tau(ranks) / Anum(ranks) )
!
!      do im = 0, nm
!        do iv = 1, 2*nv
!          do iz = -nz, nz-1
!            do mx = ist_xw, iend_xw
!              do my = 0, nyw-1
!                wx_nl1 = abs( dx_inv * real( eyw_xw(my,mx,iz,im) - cs * vl(iv) &
!                                             * byw_xw(my,mx,iz,im), kind=DP ) )
!                wx_nl2 = abs( dx_inv * aimag( eyw_xw(my,mx,iz,im) - cs * vl(iv) &
!                                             * byw_xw(my,mx,iz,im)          ) )
!                wy_nl1 = abs( dy_inv * real( exw_xw(my,mx,iz,im) - cs * vl(iv) &
!                                             * bxw_xw(my,mx,iz,im), kind=DP ) )
!                wy_nl2 = abs( dy_inv * aimag( exw_xw(my,mx,iz,im) - cs * vl(iv) &
!                                             * bxw_xw(my,mx,iz,im)          ) )
!                if ( w_nl_max < wx_nl1 ) w_nl_max = wx_nl1
!                if ( w_nl_max < wx_nl2 ) w_nl_max = wx_nl2
!                if ( w_nl_max < wy_nl1 ) w_nl_max = wy_nl1
!                if ( w_nl_max < wy_nl2 ) w_nl_max = wy_nl2
!              end do
!            end do
!          end do
!        end do
!      end do
!  
!      call MPI_Allreduce( w_nl_max, w_nl_max2, 1, MPI_DOUBLE_PRECISION, &
!                          MPI_MAX, MPI_COMM_WORLD, ierr_mpi )
!  
!      dt_nl = courant_num / w_nl_max2
!  
!      dt_limit = min( dt_max, dt_linear, dt_nl )
!
!  
!  END SUBROUTINE dtc_estimate

!--------------------------------------
  SUBROUTINE dtc_estimate
!--------------------------------------

    real(kind=DP) :: w_nl_max, w_nl_max2

      w_nl_max = max(eps, dx_inv*exb_maxvx_eachrank, dy_inv*exb_maxvy_eachrank)
  
      call MPI_Allreduce( w_nl_max, w_nl_max2, 1, MPI_DOUBLE_PRECISION, &
                          MPI_MAX, MPI_COMM_WORLD, ierr_mpi )
  
      dt_nl = courant_num / w_nl_max2
  
      dt_limit = min( dt_max, dt_linear, dt_nl )

  
  END SUBROUTINE dtc_estimate


END MODULE GKV_dtc
