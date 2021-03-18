MODULE GKV_colli
!-------------------------------------------------------------------------------
!
!    Collision term
!
!      GKV-plus r0.3 ( T.-H.Watanabe, Jun 2011)
!
!-------------------------------------------------------------------------------

  use GKV_header
  use GKV_mpienv
  use GKV_clock, only : clock_sta, clock_end
  use GKV_bndry, only : bndry_vm_buffin, bndry_vm_sendrecv, bndry_vm_buffout, &
           bndry_shifts_v_buffin, bndry_shifts_v_sendrecv, bndry_shifts_v_buffout, &
           bndry_shifts_m_buffin, bndry_shifts_m_sendrecv, bndry_shifts_m_buffout, &
           bndry_vm_sendrecv_v2

  implicit none

  private

  integer, save :: nchunk_xy = 1, nchunk_yvb = 1, nchunk_ymb = 1

  public   colli_set_param, colli_LB!, colli_full


CONTAINS


!--------------------------------------
  SUBROUTINE colli_set_param (q0, eps_r, nust)
!-------------------------------------------------------------------------------
!
!    Set parameters for GK collision term
!
!    by M. Nakata and M. Nunami, April 2014
!
!-------------------------------------------------------------------------------

    real(kind=DP), parameter :: mp      = 1.67262178d-24, & ! proton mass in g
                                ee      = 4.80320425d-10, & ! elementary charge in esu
                                ev2erg  = 1.60217657d-12    ! erg/eV  (cf. 1J = 10^7 erg)
    
    real(kind=DP),                    intent(in)  :: q0, eps_r
    real(kind=DP), dimension(0:ns-1,0:ns-1), intent(out) :: nust

    real(kind=DP), dimension(0:ns-1)        :: tmpr, dens, freq_factor
    real(kind=DP), dimension(0:ns-1,0:ns-1) :: log_lambda
!    real(kind=DP)                           :: cph, dph, cgg

    integer :: is, is1, is2
  
! --- temperature [in eV] and density [in cm^(-3)]
    do is = 0, ns-1
      tmpr(is) = tau(is) * Tref*1.d3
      dens(is) = Nref*1.d-6 * fcs(is)/Znum(is)
    end do

! --- factor for collision frequencies
    do is = 0, ns-1
      freq_factor(is)  = (dens(is) * ee**4 * Lref*1.d2) / (Tref*1.d3*ev2erg)**2 
    end do


! --- Coulomb logarithm in cm^(-3) and eV units (see NRL plasma Formulary)  
    do is1 = 0, ns-1
      if (sgn(is1) < 0.d0) then  !! For is1 = electron

        do is2 = 0, ns-1 
          if (sgn(is2) < 0.d0) then   !! e-e case
            if (dens(is1) == 0.d0) then !-care for tracer particle(dens=0)-
              log_lambda(is1,is2) = 0._DP
            else
              log_lambda(is1,is2) = 23.5_DP - dlog( dsqrt( dens(is1) ) * tmpr(is1)**(-1.25_DP) )  &
                                            - dsqrt( 1.d-5 + (( dlog(tmpr(is1)) - 2._DP )**2 )/16._DP )
            end if
          else                        !! e-i case
            if (dens(is1) == 0.d0) then !-care for tracer particle(dens=0)-
              log_lambda(is1,is2) = 0._DP
            else
              log_lambda(is1,is2) = 24._DP - dlog( dsqrt( dens(is1) ) / tmpr(is1) )
            end if
          end if
        end do

      else                     !! For is1 = ions

        do is2 = 0, ns-1 
          if (sgn(is2) < 0.d0) then   !! i-e case
            if (dens(is2) == 0.d0) then !-care for tracer particle(dens=0)-
              log_lambda(is1,is2) = 0._DP
            else
              log_lambda(is1,is2) = 24._DP - dlog( dsqrt( dens(is2) ) / tmpr(is2) )
            end if
          else                       !! i-i case
            if (dens(is1) == 0.d0 .and. dens(is2) == 0.d0) then !-care for tracer particle(dens=0)-
              log_lambda(is1,is2) = 0._DP
            else
              log_lambda(is1,is2) = 23._DP &
                - dlog( Znum(is1)*Znum(is2)*(Anum(is1)+Anum(is2))/(Anum(is1)*tmpr(is2)+Anum(is2)*tmpr(is1)) &
                        * dsqrt( (dens(is1) * Znum(is1)**2)/tmpr(is1)                                       &
                               + (dens(is2) * Znum(is2)**2)/tmpr(is2) ) )
            end if
          end if
        end do

      end if
    end do

! --- Constant parameters
    do is1 = 0, ns-1
      do is2 = 0, ns-1 

         ctauiv(is1,is2) = freq_factor(is2) * (8._DP*dsqrt(pi)/3._DP/dsqrt(2._DP))*log_lambda(is1,is2)  & 
                                   * ( Znum(is1)**2*Znum(is2)**2/dsqrt(Anum(is1))/tau(is1)**1.5 )

         calpha(is1,is2) = dsqrt( tau(is1) * Anum(is2) / ( tau(is2) * Anum(is1) ) )
         ctheta(is1,is2) = dsqrt( tau(is1) * ( Anum(is1) + Anum(is2) ) / ( tau(is1) * Anum(is2) + tau(is2) * Anum(is1) ) )

         cgamma(is1,is2) = - Anum(is1) * calpha(is1,is2)                                      &
                            * (tau(is1)/tau(is2) + calpha(is1,is2)**2) * ctauiv(is1,is2)      &
                             / (1._DP + calpha(is1,is2)**2)**1.5_DP 

         ceta(is1,is2)   = - tau(is1) * 3._DP * calpha(is1,is2)                               &
                            * (tau(is1)/tau(is2) + calpha(is1,is2)**2) * ctauiv(is1,is2)      &
                             / (1._DP + calpha(is1,is2)**2)**2.5_DP

          cxi(is1,is2)   =  calpha(is1,is2) * ( ctheta(is1,is2) - 1._DP ) * ctauiv(is1,is2)   &
                             / dsqrt(1._DP + calpha(is1,is2)**2) 

         nust(is1,is2)   = q0*(ctauiv(is1,is2)/dsqrt(2._DP))/(eps_r**1.5*dsqrt(tau(is1)/Anum(is1)))

      end do
    end do

!!!%%% Parameters for colli_full %%%
!!!! --- xxa = v/vta/sqrt(2), where vta = sqrt(Ta/ma)
!!!    do im = 0, nm
!!!      do iv = 1, 2*nv
!!!        do iz = -nz, nz-1
!!!          xxa(iz,iv,im) = dsqrt(vl(iv)**2 + vp(iz,im)**2)/dsqrt(2._DP) 
!!!        end do 
!!!      end do
!!!    end do
!!!
!!!! --- collision frequencies and v-space functions
!!!    do im = 0, nm
!!!      do iv = 1, 2*nv
!!!        do iz = -nz, nz-1 
!!!          do is1 = 0, ns-1
!!!            do is2 = 0, ns-1 
!!!
!!!              cph = derf(calpha(is1,is2)*xxa(iz,iv,im))
!!!              dph = 2._DP/dsqrt(pi)*dexp(-calpha(is1,is2)**2*xxa(iz,iv,im)**2)
!!!              cgg = (cph - calpha(is1,is2)*xxa(iz,iv,im)*dph)/(calpha(is1,is2)**2*xxa(iz,iv,im)**2)*0.5_DP
!!!
!!!              nu_d(iz,iv,im,is1,is2) = 0.75_DP*dsqrt(pi)*ctauiv(is1,is2)*(cph-cgg)/xxa(iz,iv,im)**3
!!!              nu_p(iz,iv,im,is1,is2) = 1.50_DP*dsqrt(pi)*ctauiv(is1,is2)*(  cgg  )/xxa(iz,iv,im)**3
!!!              nu_h(iz,iv,im,is1,is2) = 0.75_DP*dsqrt(pi)*ctauiv(is1,is2)*calpha(is1,is2)*dph/xxa(iz,iv,im)**2
!!!              nu_g(iz,iv,im,is1,is2) = nu_p(iz,iv,im,is1,is2)*xxa(iz,iv,im)**2*(1._DP-calpha(is1,is2)**2)
!!!
!!!              c_t0(iz,iv,im,is1,is2,1)  = - (1._DP + calpha(is1,is2)**2)*fmx(iz,iv,im)*nu_p(iz,iv,im,is1,is2)               &
!!!                                             * xxa(iz,iv,im)**2*vl(iv)
!!!              c_t0(iz,iv,im,is1,is2,2)  = - 1.5_DP*dsqrt(pi)*ctauiv(is1,is2)*fmx(iz,iv,im)                                  & 
!!!                                             * ( cph - calpha(is1,is2)*xxa(iz,iv,im)*(1._DP + calpha(is1,is2)**2)*dph )     & 
!!!                                             / calpha(is1,is2)**2 / xxa(iz,iv,im)
!!!
!!!              x_tst(1,iz,iv,im,is1,is2) = (ctheta(is1,is2) - 1._DP)*fmx(iz,iv,im)*vl(iv)
!!!              x_tst(2,iz,iv,im,is1,is2) = x_tst(1,iz,iv,im,is1,is2)*vp(iz,im)/vl(iv) 
!!!              x_tst(3,iz,iv,im,is1,is2) = x_tst(1,iz,iv,im,is1,is2)*(xxa(iz,iv,im)**2/1.5_DP - 1._DP)/vl(iv)
!!!              x_tst(4,iz,iv,im,is1,is2) = (ctheta(is1,is2) - 1._DP)*( c_t0(iz,iv,im,is1,is2,1)                              &
!!!                                                      - (ctheta(is1,is2) - 1._DP)*calpha(is1,is2)*ctauiv(is1,is2)           &
!!!                                                              * fmx(iz,iv,im)*vl(iv)/dsqrt(1._DP + calpha(is1,is2)**2) )
!!!              x_tst(5,iz,iv,im,is1,is2) =  x_tst(4,iz,iv,im,is1,is2)*vp(iz,im)/vl(iv)  
!!!              x_tst(6,iz,iv,im,is1,is2) = (ctheta(is1,is2) - 1._DP)*( c_t0(iz,iv,im,is1,is2,2)*2._DP/3._DP                &
!!!                                                      - (ctheta(is1,is2) - 1._DP)*calpha(is1,is2)*ctauiv(is1,is2)         &
!!!                                                              * fmx(iz,iv,im)*(xxa(iz,iv,im)**2/1.5_DP - 1._DP)*2._DP     &
!!!                                                               / (1._DP + calpha(is1,is2)**2)**1.5 )
!!!
!!!              y_fld(1,iz,iv,im,is1,is2) = - (fcs(is2)/Znum(is2))/(fcs(is1)/Znum(is1))*calpha(is1,is2)*Anum(is1)             & 
!!!                                                    * tau(is2)*ctheta(is1,is2)*ctheta(is2,is1)/tau(is1)/cgamma(is1,is2)     &
!!!                                                    * ( c_t0(iz,iv,im,is1,is2,1) - cxi(is1,is2)*fmx(iz,iv,im)*vl(iv) ) 
!!!              y_fld(2,iz,iv,im,is1,is2) = y_fld(1,iz,iv,im,is1,is2)*vp(iz,im)/vl(iv)  
!!!              y_fld(3,iz,iv,im,is1,is2) = - (fcs(is2)/Znum(is2))/(fcs(is1)/Znum(is1))                                     & 
!!!                                                    * tau(is2)*ctheta(is1,is2)*ctheta(is2,is1)/ceta(is1,is2)              &
!!!                                                    * ( c_t0(iz,iv,im,is1,is2,2)                                          &
!!!                                                           - cxi(is1,is2)/(1._DP+calpha(is1,is2)**2)*fmx(iz,iv,im)        &
!!!                                                               *(2._DP*xxa(iz,iv,im)**2 - 3._DP) ) 
!!!              y_fld(4,iz,iv,im,is1,is2) = - y_fld(1,iz,iv,im,is1,is2)*cxi(is2,is1) 
!!!              y_fld(5,iz,iv,im,is1,is2) = - y_fld(2,iz,iv,im,is1,is2)*cxi(is2,is1) 
!!!              y_fld(6,iz,iv,im,is1,is2) = - y_fld(3,iz,iv,im,is1,is2)*2._DP*cxi(is2,is1)/(1._DP+calpha(is2,is1)**2) 
!!!
!!!            end do 
!!!          end do 
!!!        end do 
!!!      end do
!!!    end do 
!!!
!!!! --- summation of collision frequencies with respect to is2, and adiabatic term (used in colli_GK_CT)
!!!    nu_hs = 0._DP 
!!!    nu_gs = 0._DP 
!!!    nu_ds = 0._DP 
!!!    nu_ps = 0._DP 
!!!    is1 = ranks
!!!      do im = 0, nm
!!!        do iv = 1, 2*nv
!!!          do iz = -nz, nz-1 
!!!            do is2 = 0, ns-1 
!!!                nu_hs(iz,iv,im) = nu_hs(iz,iv,im) + nu_h(iz,iv,im,is1,is2)
!!!                nu_gs(iz,iv,im) = nu_gs(iz,iv,im) + nu_g(iz,iv,im,is1,is2)
!!!                nu_ds(iz,iv,im) = nu_ds(iz,iv,im) + nu_d(iz,iv,im,is1,is2)
!!!                nu_ps(iz,iv,im) = nu_ps(iz,iv,im) + nu_p(iz,iv,im,is1,is2)
!!!            end do 
!!!          end do 
!!!        end do
!!!      end do 
!!!
!!!    if (trim(col_type) == "lorentz") then
!!!      nu_hs(:,:,:) = 0._DP
!!!      nu_ps(:,:,:) = 0._DP
!!!      x_tst(:,:,:,:,:,:) = 0._DP
!!!      y_fld(:,:,:,:,:,:) = 0._DP
!!!      nu_gs(:,:,:) = - nu_ds(:,:,:)
!!!    end if
!!!
!!!! --- adiabatic part (used in colli_GK_CT)
!!!    is1 = ranks
!!!      do im = 0, nm
!!!
!!!        if ( rankm == 0 .AND. im == 0 ) then
!!!
!!!          do iv = 1, 2*nv
!!!            do iz = -nz, nz-1 
!!!              do my = ist_y, iend_y 
!!!                do mx = 0, nx
!!!                  adbtc(mx,my,iz,iv,im) =                                              &
!!!                               - ( nu_ds(iz,iv,im)*vl(iv)**2                           &
!!!                                     * ksq(mx,my,iz)*Anum(is1)/Znum(is1)/omg(iz)**2    &
!!!                                 ) * fmx(iz,iv,im)*sgn(is1)*real(iFLR, kind=DP)
!!!                end do 
!!!              end do 
!!!            end do 
!!!          end do 
!!!
!!!        else
!!!
!!!          do iv = 1, 2*nv
!!!            do iz = -nz, nz-1 
!!!              do my = ist_y, iend_y 
!!!                do mx = 0, nx
!!!                  adbtc(mx,my,iz,iv,im) =                                              &
!!!                                 ( -( nu_ds(iz,iv,im)*vl(iv)**2                        &
!!!                                       + nu_ps(iz,iv,im)*vp(iz,im)**2 )                &
!!!                                     * 0.25_DP*ksq(mx,my,iz)*Anum(is1)                 &
!!!                                             /Znum(is1)/omg(iz)**2                     &
!!!                                     * (j0(mx,my,iz,im) - j2(mx,my,iz,im))             &
!!!                                   -( nu_hs(iz,iv,im)*vp(iz,im)                        &
!!!                                       - 0.5_DP*nu_ps(iz,iv,im)*vp(iz,im)              &
!!!                                               *(1._DP-vl(iv)**2-vp(iz,im)**2)         &
!!!                                       + 0.5_DP*nu_ds(iz,iv,im)                        &
!!!                                               *(vl(iv)**2/vp(iz,im)-vp(iz,im)) )      &
!!!                                     * dsqrt(ksq(mx,my,iz)*Anum(is1)/tau(is1))/omg(iz) &
!!!                                     * j1(mx,my,iz,im)                                 &
!!!                                   -( nu_ds(iz,iv,im)                                  &
!!!                                               *(2._DP*vl(iv)**2+vp(iz,im)**2)         &
!!!                                       + nu_ps(iz,iv,im)*vp(iz,im)**2 )                &
!!!                                     * 0.25_DP*ksq(mx,my,iz)*Anum(is1)                 &
!!!                                             /Znum(is1)/omg(iz)**2 * j0(mx,my,iz,im)   &
!!!                                 ) * fmx(iz,iv,im)*sgn(is1)*real(iFLR, kind=DP)
!!!                end do 
!!!              end do 
!!!            end do 
!!!          end do 
!!!
!!!        end if
!!!
!!!      end do
!!!
!!!
!!!! --- set v-space functions used in colli_moment
!!!
!!!    vfunc(:,:,:,:,:) = 0._DP
!!!    jfunc(:,:,:,:,:) = 0._DP
!!!
!!!    if ( iFLR == 1 ) then
!!!      !is1 = ranks
!!!      !  do is2 = 0, ns-1
!!!      !    do im = 0, nm
!!!      !      do iv = 1, 2*nv
!!!      !        do iz = -nz, nz-1
!!!      !          do my = ist_y, iend_y
!!!      !            do mx = -nx,nx
!!!      !              vfunc(mx,my,iz,iv,im,is2,1) = & 
!!!      !                              j0(mx,my,iz,im)*c_t0(iz,iv,im,is1,is2,1)/fmx(iz,iv,im) 
!!!      !              vfunc(mx,my,iz,iv,im,is2,2) = & 
!!!      !                              j1(mx,my,iz,im)*vp(iz,im)*c_t0(iz,iv,im,is1,is2,1)/vl(iv)/fmx(iz,iv,im)
!!!      !              vfunc(mx,my,iz,iv,im,is2,3) = & 
!!!      !                              j0(mx,my,iz,im)*c_t0(iz,iv,im,is1,is2,2)/fmx(iz,iv,im)
!!!      !              vfunc(mx,my,iz,iv,im,is2,4) = j0(mx,my,iz,im)*vl(iv)
!!!      !              vfunc(mx,my,iz,iv,im,is2,5) = j1(mx,my,iz,im)*vp(iz,im)
!!!      !              vfunc(mx,my,iz,iv,im,is2,6) = j0(mx,my,iz,im)*(xxa(iz,iv,im)**2-1.5_DP)
!!!      !            end do 
!!!      !          end do 
!!!      !        end do 
!!!      !      end do 
!!!      !    end do 
!!!      !  end do 
!!!      is1 = ranks
!!!        do is2 = 0, ns-1
!!!          do im = 0, nm
!!!            do iv = 1, 2*nv
!!!              do iz = -nz, nz-1
!!!                vfunc(iz,iv,im,is2,1) = c_t0(iz,iv,im,is1,is2,1)/fmx(iz,iv,im) 
!!!                vfunc(iz,iv,im,is2,2) = vp(iz,im)*c_t0(iz,iv,im,is1,is2,1)/vl(iv)/fmx(iz,iv,im)
!!!                vfunc(iz,iv,im,is2,3) = c_t0(iz,iv,im,is1,is2,2)/fmx(iz,iv,im)
!!!                vfunc(iz,iv,im,is2,4) = vl(iv)
!!!                vfunc(iz,iv,im,is2,5) = vp(iz,im)
!!!                vfunc(iz,iv,im,is2,6) = xxa(iz,iv,im)**2-1.5_DP
!!!              end do 
!!!            end do 
!!!          end do 
!!!        end do 
!!!      do im = 0, nm
!!!        do iz = -nz, nz-1
!!!          do my = ist_y, iend_y
!!!            do mx = -nx,nx
!!!              jfunc(mx,my,iz,im,1) = j0(mx,my,iz,im)
!!!              jfunc(mx,my,iz,im,2) = j1(mx,my,iz,im)
!!!              jfunc(mx,my,iz,im,3) = j0(mx,my,iz,im)
!!!              jfunc(mx,my,iz,im,4) = j0(mx,my,iz,im)
!!!              jfunc(mx,my,iz,im,5) = j1(mx,my,iz,im)
!!!              jfunc(mx,my,iz,im,6) = j0(mx,my,iz,im)
!!!            end do 
!!!          end do 
!!!        end do 
!!!      end do 
!!!   
!!!    else 
!!!
!!!      !is1 = ranks
!!!      !  do is2 = 0, ns-1
!!!      !    do im = 0, nm
!!!      !      do iv = 1, 2*nv
!!!      !        do iz = -nz, nz-1
!!!      !          do my = ist_y, iend_y
!!!      !            do mx = -nx,nx
!!!      !              vfunc(mx,my,iz,iv,im,is2,1) = c_t0(iz,iv,im,is1,is2,1)/fmx(iz,iv,im) 
!!!      !              vfunc(mx,my,iz,iv,im,is2,2) = 0._DP
!!!      !              vfunc(mx,my,iz,iv,im,is2,3) = c_t0(iz,iv,im,is1,is2,2)/fmx(iz,iv,im)
!!!      !              vfunc(mx,my,iz,iv,im,is2,4) = vl(iv)
!!!      !              vfunc(mx,my,iz,iv,im,is2,5) = 0._DP
!!!      !              vfunc(mx,my,iz,iv,im,is2,6) = (xxa(iz,iv,im)**2-1.5_DP)
!!!      !            end do 
!!!      !          end do 
!!!      !        end do 
!!!      !      end do 
!!!      !    end do 
!!!      !  end do 
!!!      is1 = ranks
!!!        do is2 = 0, ns-1
!!!          do im = 0, nm
!!!            do iv = 1, 2*nv
!!!              do iz = -nz, nz-1
!!!                vfunc(iz,iv,im,is2,1) = c_t0(iz,iv,im,is1,is2,1)/fmx(iz,iv,im) 
!!!                vfunc(iz,iv,im,is2,2) = 0._DP
!!!                vfunc(iz,iv,im,is2,3) = c_t0(iz,iv,im,is1,is2,2)/fmx(iz,iv,im)
!!!                vfunc(iz,iv,im,is2,4) = vl(iv)
!!!                vfunc(iz,iv,im,is2,5) = 0._DP
!!!                vfunc(iz,iv,im,is2,6) = (xxa(iz,iv,im)**2-1.5_DP)
!!!              end do 
!!!            end do 
!!!          end do 
!!!        end do 
!!!      do im = 0, nm
!!!        do iz = -nz, nz-1
!!!          do my = ist_y, iend_y
!!!            do mx = -nx,nx
!!!              jfunc(mx,my,iz,im,1) = 1._DP
!!!              jfunc(mx,my,iz,im,2) = 0._DP
!!!              jfunc(mx,my,iz,im,3) = 1._DP
!!!              jfunc(mx,my,iz,im,4) = 1._DP
!!!              jfunc(mx,my,iz,im,5) = 0._DP
!!!              jfunc(mx,my,iz,im,6) = 1._DP
!!!            end do 
!!!          end do 
!!!        end do 
!!!      end do 
!!!
!!!    end if
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


! -----------------------------------
! --- Output constants
    if ( rankg == nprocz/2 ) then

      do is1 = 0, ns-1
        do is2 = 0, ns-1 
          write(unit=ocst,fmt="(2I3,SP,256ES24.15e3)") is1, is2, ctheta(is1,is2), calpha(is1,is2), &
                                                                 fcs(is1)/Znum(is1)*ceta(is1,is2), &
                                                               fcs(is1)/Znum(is1)*cgamma(is1,is2), & 
                                                                    cxi(is1,is2), ctauiv(is1,is2), &
                                                                              log_lambda(is1,is2)
! --- Note that, for ns >=3, cgamma(is1,is2) /= cgamma(is2,is1), but dens(is1)*cgamma(is1,is2) = dense(is2)*cgamma(is2,is1)
! ---  due to normalizartion with dens(is). 
        end do
      end do

!! --- for debug
!      iz = -nz
!      do im = 0, nm
!        do iv = 1, 2*nv
!          write(unit=4000,fmt="(2f15.8,SP,256ES24.15e3)") vl(iv), vp(iz,im), & 
!                                    ((nu_h(iz,iv,im,is1,is2), is2 = 0, ns-1), is1 = 0, ns-1)
!          write(unit=4001,fmt="(2f15.8,SP,256ES24.15e3)") vl(iv), vp(iz,im), & 
!                                    ((nu_g(iz,iv,im,is1,is2), is2 = 0, ns-1), is1 = 0, ns-1)
!          write(unit=4002,fmt="(2f15.8,SP,256ES24.15e3)") vl(iv), vp(iz,im), & 
!                                    ((nu_d(iz,iv,im,is1,is2), is2 = 0, ns-1), is1 = 0, ns-1)
!          write(unit=4003,fmt="(2f15.8,SP,256ES24.15e3)") vl(iv), vp(iz,im), & 
!                                    ((nu_p(iz,iv,im,is1,is2), is2 = 0, ns-1), is1 = 0, ns-1)
!          write(unit=4004,fmt="(2f15.8,SP,256ES24.15e3)") vl(iv), vp(iz,im), & 
!                                  ((c_t0(iz,iv,im,is1,is2,1), is2 = 0, ns-1), is1 = 0, ns-1)
!          write(unit=4005,fmt="(2f15.8,SP,256ES24.15e3)") vl(iv), vp(iz,im), & 
!                                  ((c_t0(iz,iv,im,is1,is2,2), is2 = 0, ns-1), is1 = 0, ns-1)
!        end do
!        write (unit=4000,fmt=*)
!        write (unit=4001,fmt=*)
!        write (unit=4002,fmt=*)
!        write (unit=4003,fmt=*)
!        write (unit=4004,fmt=*)
!        write (unit=4005,fmt=*)
!      end do
!
!! --- for debug
!      iz = -nz
!      do im = 0, nm
!        do iv = 1, 2*nv
!          write(unit=5001,fmt="(2f15.8,SP,256ES24.15e3)") vl(iv), vp(iz,im), & 
!                                    ((x_tst(iz,iv,im,is1,is2,1), is2 = 0, ns-1), is1 = 0, ns-1)
!          write(unit=5002,fmt="(2f15.8,SP,256ES24.15e3)") vl(iv), vp(iz,im), & 
!                                    ((x_tst(iz,iv,im,is1,is2,2), is2 = 0, ns-1), is1 = 0, ns-1)
!          write(unit=5003,fmt="(2f15.8,SP,256ES24.15e3)") vl(iv), vp(iz,im), & 
!                                    ((x_tst(iz,iv,im,is1,is2,3), is2 = 0, ns-1), is1 = 0, ns-1)
!          write(unit=5004,fmt="(2f15.8,SP,256ES24.15e3)") vl(iv), vp(iz,im), & 
!                                    ((x_tst(iz,iv,im,is1,is2,4), is2 = 0, ns-1), is1 = 0, ns-1)
!          write(unit=5005,fmt="(2f15.8,SP,256ES24.15e3)") vl(iv), vp(iz,im), & 
!                                    ((x_tst(iz,iv,im,is1,is2,5), is2 = 0, ns-1), is1 = 0, ns-1)
!          write(unit=5006,fmt="(2f15.8,SP,256ES24.15e3)") vl(iv), vp(iz,im), & 
!                                    ((x_tst(iz,iv,im,is1,is2,6), is2 = 0, ns-1), is1 = 0, ns-1)
!          write(unit=6001,fmt="(2f15.8,SP,256ES24.15e3)") vl(iv), vp(iz,im), & 
!                                    ((y_fld(iz,iv,im,is1,is2,1), is2 = 0, ns-1), is1 = 0, ns-1)
!          write(unit=6002,fmt="(2f15.8,SP,256ES24.15e3)") vl(iv), vp(iz,im), & 
!                                    ((y_fld(iz,iv,im,is1,is2,2), is2 = 0, ns-1), is1 = 0, ns-1)
!          write(unit=6003,fmt="(2f15.8,SP,256ES24.15e3)") vl(iv), vp(iz,im), & 
!                                    ((y_fld(iz,iv,im,is1,is2,3), is2 = 0, ns-1), is1 = 0, ns-1)
!          write(unit=6004,fmt="(2f15.8,SP,256ES24.15e3)") vl(iv), vp(iz,im), & 
!                                    ((y_fld(iz,iv,im,is1,is2,4), is2 = 0, ns-1), is1 = 0, ns-1)
!          write(unit=6005,fmt="(2f15.8,SP,256ES24.15e3)") vl(iv), vp(iz,im), & 
!                                    ((y_fld(iz,iv,im,is1,is2,5), is2 = 0, ns-1), is1 = 0, ns-1)
!          write(unit=6006,fmt="(2f15.8,SP,256ES24.15e3)") vl(iv), vp(iz,im), & 
!                                    ((y_fld(iz,iv,im,is1,is2,6), is2 = 0, ns-1), is1 = 0, ns-1)
!        end do
!        write (unit=5001,fmt=*)
!        write (unit=5002,fmt=*)
!        write (unit=5003,fmt=*)
!        write (unit=5004,fmt=*)
!        write (unit=5005,fmt=*)
!        write (unit=5006,fmt=*)
!        write (unit=6001,fmt=*)
!        write (unit=6002,fmt=*)
!        write (unit=6003,fmt=*)
!        write (unit=6004,fmt=*)
!        write (unit=6005,fmt=*)
!        write (unit=6006,fmt=*)
!      end do

    end if

    return

   END SUBROUTINE colli_set_param


!--------------------------------------
  SUBROUTINE colli_LB( ff, phi, cf )
!--------------------------------------
!   Lenard-Bernstein model collsion operator

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: phi
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: cf

    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: vb1e, vb2e, vb1o, vb2o
    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: mb1e, mb2e, mb1o, mb2o
    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: gge, ggo
    integer :: iz
    integer, save :: iflg
    data iflg / 0 /
!$  integer :: nthreads, omp_get_num_threads

      if ( iflg == 0 ) then
        iflg = 1
!$OMP parallel default(shared)
!$OMP master
!$    nthreads = omp_get_num_threads()
!$    if (nthreads > 1) then
!$      nchunk_xy = ((2*nx+1)*(ny+1)-1) / (nthreads-1) + 1
!$      nchunk_yvb = ((ny+1)*(2*nv)*(2*nvb)-1) / (nthreads-1) + 1
!$      nchunk_ymb = ((ny+1)*(nm+1)*(2*nvb)-1) / (nthreads-1) + 1
!$    end if
!$OMP end master
!$OMP end parallel
      end if

      allocate( vb1e(-nx:nx,0:ny,-nz:nz-1,0:nm,1:2*nvb) )
      allocate( vb2e(-nx:nx,0:ny,-nz:nz-1,0:nm,1:2*nvb) )
      allocate( mb1e(-nx:nx,0:ny,-nz:nz-1,1:2*nv,1:2*nvb) )
      allocate( mb2e(-nx:nx,0:ny,-nz:nz-1,1:2*nv,1:2*nvb) )
      !!!allocate( gge(1-nvb:2*nv+nvb,0-nvb:nm+nvb,-nx:nx,0:ny,-nz:nz-1) )
      ! mod by M.Nakata 20200928
      allocate( gge(-nx:nx,0:ny,-nz:nz-1,1-nvb:2*nv+nvb,0-nvb:nm+nvb) )
!     allocate( vb1o(-nx:nx,0:ny,0:nm,1:2*nvb) )
!     allocate( vb2o(-nx:nx,0:ny,0:nm,1:2*nvb) )
!     allocate( mb1o(-nx:nx,0:ny,1:2*nv,1:2*nvb) )
!     allocate( mb2o(-nx:nx,0:ny,1:2*nv,1:2*nvb) )
!     allocate( ggo(1-nvb:2*nv+nvb,0-nvb:nm+nvb,-nx:nx,0:ny) )

!$OMP parallel default(none) &
!$OMP shared(ff,phi,cf,vb1e,vb2e,mb1e,mb2e,gge) &
!$OMP shared(vb1o,vb2o,mb1o,mb2o,ggo) &
!$OMP private(iz)

!$OMP workshare
      vb1e(:,:,:,:,:) = (0._DP, 0._DP)
      vb2e(:,:,:,:,:) = (0._DP, 0._DP)
      mb1e(:,:,:,:,:) = (0._DP, 0._DP)
      mb2e(:,:,:,:,:) = (0._DP, 0._DP)
       gge(:,:,:,:,:) = (0._DP, 0._DP)
!     vb1o(:,:,:,:) = (0._DP, 0._DP)
!     vb2o(:,:,:,:) = (0._DP, 0._DP)
!     mb1o(:,:,:,:) = (0._DP, 0._DP)
!     mb2o(:,:,:,:) = (0._DP, 0._DP)
!      ggo(:,:,:,:) = (0._DP, 0._DP)
      cf(:,:,:,:,:) = (0._DP, 0._DP)
!$OMP end workshare

!!%%% Without overlap %%%
      call colli_LB_buffin_v2(ff, phi, vb1e, mb1e)
      call bndry_vm_sendrecv_v2(vb1e, vb2e, mb1e, mb2e)
      !!!call colli_LB_buffout_v2(ff, phi, vb2e, mb2e, gge)
      !!!call colli_LB_calc_v2(gge, cf)
      ! mod by M.Nakata 20200928
      call colli_LB_buffout_v3_sx(ff, phi, vb2e, mb2e, gge)   
      call colli_LB_calc_v3_sx(gge, cf)
!!%%%%%%%%%%%%%%%%%%%%%%%


!%%% With overlap %%%
!      do iz = -nz, nz-1+3
!        if (mod(iz+nz,2) == 0) then ! even
!!$OMP master
!          if (-nz+1<=iz .and. iz<=nz-1+1) call bndry_vm_sendrecv(vb1o, vb2o, mb1o, mb2o)
!!$OMP end master
!          if (-nz  <=iz .and. iz<=nz-1  ) call colli_LB_buffin(iz, ff, phi, vb1e, mb1e)
!          if (-nz+2<=iz .and. iz<=nz-1+2) call colli_LB_buffout(iz-2, ff, phi, vb2e, mb2e, gge)
!          if (-nz+3<=iz .and. iz<=nz-1+3) call colli_LB_calc(iz-3, ggo, cf)
!        else                        ! odd
!!$OMP master
!          if (-nz+1<=iz .and. iz<=nz-1+1) call bndry_vm_sendrecv(vb1e, vb2e, mb1e, mb2e)
!!$OMP end master
!          if (-nz  <=iz .and. iz<=nz-1  ) call colli_LB_buffin(iz, ff, phi, vb1o, mb1o)
!          if (-nz+2<=iz .and. iz<=nz-1+2) call colli_LB_buffout(iz-2, ff, phi, vb2o, mb2o, ggo)
!          if (-nz+3<=iz .and. iz<=nz-1+3) call colli_LB_calc(iz-3, gge, cf)
!        end if
!!$OMP barrier
!      end do
!%%%%%%%%%%%%%%%%%%%%

!$OMP end parallel

      deallocate( vb1e )
      deallocate( vb2e )
      deallocate( mb1e )
      deallocate( mb2e )
      deallocate( gge )
!     deallocate( vb1o )
!     deallocate( vb2o )
!     deallocate( mb1o )
!     deallocate( mb2o )
!     deallocate( ggo )

  END SUBROUTINE colli_LB


!--------------------------------------
  SUBROUTINE colli_LB_buffin( iz, ff, phi, vb1, mb1 )
!--------------------------------------
!   Lenard-Bernstein model collsion operator

    integer, intent(in) :: iz
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: phi
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,0:nm,1:2*nvb) :: vb1
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,1:2*nv,1:2*nvb) :: mb1

    real(kind=DP) :: cs1
    integer :: mx, my, iv, im

!$OMP master
                                           call clock_sta(1371)
                                         ! call fapp_start("literm_shifts_bufferin",1371,1)
!$OMP end master

      cs1 = sgn(ranks) * Znum(ranks) / tau(ranks) * real(iFLR, kind=DP)

!$OMP do collapse(3) schedule(dynamic,nchunk_ymb)
      do iv = 1, nvb
        do im = 0, nm
          do my = ist_y, iend_y
            do mx = -nx, nx
              vb1(mx,my,im,iv    ) = ff(mx,my,iz,         iv,im) &
                                  + cs1 * fmx(iz,         iv,im) &
                               * phi(mx,my,iz) * j0(mx,my,iz,im)
              vb1(mx,my,im,iv+nvb) = ff(mx,my,iz,2*nv-nvb+iv,im) &
                                  + cs1 * fmx(iz,2*nv-nvb+iv,im) &
                               * phi(mx,my,iz) * j0(mx,my,iz,im)
            end do
          end do
        end do
      end do
!$OMP end do nowait

!$OMP do collapse(3) schedule(dynamic,nchunk_yvb)
      do im = 1, nvb
        do iv = 1, 2*nv
          do my = ist_y, iend_y
            do mx = -nx, nx
              mb1(mx,my,iv,im    ) = ff(mx,my,iz,iv,     im-1) &
                                  + cs1 * fmx(iz,iv,     im-1) &
                      * phi(mx,my,iz) * j0(mx,my,iz,     im-1)
              mb1(mx,my,iv,im+nvb) = ff(mx,my,iz,iv,nm-nvb+im) &
                                  + cs1 * fmx(iz,iv,nm-nvb+im) &
                      * phi(mx,my,iz) * j0(mx,my,iz,nm-nvb+im)
            end do
          end do
        end do
      end do
!$OMP end do nowait

!$OMP master
                                         ! call fapp_stop("literm_shifts_bufferin",1371,1)
                                           call clock_end(1371)
!$OMP end master

  END SUBROUTINE colli_LB_buffin


!--------------------------------------
  SUBROUTINE colli_LB_buffout( iz, ff, phi, vb2, mb2, gg )
!--------------------------------------
!   Lenard-Bernstein model collsion operator

    integer, intent(in) :: iz
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: phi
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,0:nm,1:2*nvb) :: vb2
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,1:2*nv,1:2*nvb) :: mb2
    complex(kind=DP), intent(inout), &
      dimension(1-nvb:2*nv+nvb,0-nvb:nm+nvb,-nx:nx,0:ny) :: gg

    real(kind=DP) :: cs1
    integer  ::  mx, my, iv, im

!$OMP master
                                           call clock_sta(1373)
                                         ! call fapp_start("literm_shifts_bufferout",1373,1)
!$OMP end master

      cs1 = sgn(ranks) * Znum(ranks) / tau(ranks) * real(iFLR, kind=DP)

!$OMP do collapse(2) schedule(dynamic,nchunk_xy)
      do my = ist_y, iend_y
        do mx = -nx, nx
          do im = 0, nm
            do iv = 1, 2*nv
              gg(iv,im,mx,my) = ff(mx,my,iz,iv,im) &
                        + cs1 * fmx(iz,iv,im) * phi(mx,my,iz) * j0(mx,my,iz,im)
            end do
          end do
        end do
      end do
!$OMP end do nowait

      if ( rankv == 0 ) then
!$OMP do collapse(2) schedule(dynamic,nchunk_xy)
        do my = ist_y, iend_y
          do mx = -nx, nx
            do im = 0, nm
              do iv = 1, nvb
                gg(-nvb+iv,im,mx,my) = (0._DP, 0._DP)
                gg(2*nv+iv,im,mx,my) = vb2(mx,my,im,iv+nvb)
              end do
            end do
          end do
        end do
!$OMP end do nowait
      else if ( rankv == nprocv-1 ) then
!$OMP do collapse(2) schedule(dynamic,nchunk_xy)
        do my = ist_y, iend_y
          do mx = -nx, nx
            do im = 0, nm
              do iv = 1, nvb
                gg(-nvb+iv,im,mx,my) = vb2(mx,my,im,iv    )
                gg(2*nv+iv,im,mx,my) = (0._DP, 0._DP)
              end do
            end do
          end do
        end do
!$OMP end do nowait
      else
!$OMP do collapse(2) schedule(dynamic,nchunk_xy)
        do my = ist_y, iend_y
          do mx = -nx, nx
            do im = 0, nm
              do iv = 1, nvb
                gg(-nvb+iv,im,mx,my) = vb2(mx,my,im,iv    )
                gg(2*nv+iv,im,mx,my) = vb2(mx,my,im,iv+nvb)
              end do
            end do
          end do
        end do
!$OMP end do nowait
      end if

      if ( rankm == 0 ) then
!$OMP do collapse(2) schedule(dynamic,nchunk_xy)
        do my = ist_y, iend_y
          do mx = -nx, nx
            do im = 1, nvb
              do iv = 1, 2*nv
                gg(iv,-nvb-1+im,mx,my) = (0._DP, 0._DP)
                gg(iv,nm+im    ,mx,my) = mb2(mx,my,iv,im+nvb)
              end do
            end do
          end do
        end do
!$OMP end do nowait
      else if ( rankm == nprocm-1 ) then
!$OMP do collapse(2) schedule(dynamic,nchunk_xy)
        do my = ist_y, iend_y
          do mx = -nx, nx
            do im = 1, nvb
              do iv = 1, 2*nv
                gg(iv,-nvb-1+im,mx,my) = mb2(mx,my,iv,im    )
                gg(iv,nm+im    ,mx,my) = (0._DP, 0._DP)
              end do
            end do
          end do
        end do
!$OMP end do nowait
      else
!$OMP do collapse(2) schedule(dynamic,nchunk_xy)
        do my = ist_y, iend_y
          do mx = -nx, nx
            do im = 1, nvb
              do iv = 1, 2*nv
                gg(iv,-nvb-1+im,mx,my) = mb2(mx,my,iv,im    )
                gg(iv,nm+im    ,mx,my) = mb2(mx,my,iv,im+nvb)
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

  END SUBROUTINE colli_LB_buffout


!--------------------------------------
  SUBROUTINE colli_LB_calc( iz, gg, cf )
!--------------------------------------
!   Lenard-Bernstein model collsion operator

    integer, intent(in) :: iz
    complex(kind=DP), intent(in), &
      dimension(1-nvb:2*nv+nvb,0-nvb:nm+nvb,-nx:nx,0:ny) :: gg
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: cf

    real(kind=DP) :: nu_s, cef1, cef2, cef3, cef4, cs1, cflr
    integer  ::  mx, my, iv, im


!$OMP master
                                           call clock_sta(1311)
                                         ! call fapp_start("literm_colli",1311,1)
!$OMP end master

! --- Note that nu(ranks) is a bias factor 
      nu_s = nu(ranks)*3._DP*dsqrt(pi)*ctauiv(ranks,ranks)/4._DP

       cs1 = sgn(ranks) * Znum(ranks) / tau(ranks) * real(iFLR, kind=DP)
      cef1 = nu_s / ( 12._DP * dv * dv )
      cef2 = nu_s / ( 12._DP * dv )
      cef3 = nu_s / ( 12._DP * dvp(iz) * dvp(iz) )
      cef4 = nu_s / ( 12._DP * dvp(iz) )
      cflr = nu_s * Anum(ranks) * tau(ranks)  &
                  / ( Znum(ranks) * omg(iz) )**2 * real(iFLR, kind=DP)


      if ( rankm /= 0 ) then

!$OMP do collapse(2) schedule(dynamic,nchunk_xy)
        do my = ist_y, iend_y
          do mx = -nx, nx

            do im = 0, nm
              do iv = 1, 2*nv
                cf(mx,my,iz,iv,im) =                                   &
                          ( -          gg(iv+2,im,mx,my)               &
                            + 16._DP * gg(iv+1,im,mx,my)               &
                            - 30._DP * gg(iv  ,im,mx,my)               &
                            + 16._DP * gg(iv-1,im,mx,my)               &
                            -          gg(iv-2,im,mx,my)               &
                           ) * cef1                                    &
                        + ( -          gg(iv+2,im,mx,my)               &
                            +  8._DP * gg(iv+1,im,mx,my)               &
                            -  8._DP * gg(iv-1,im,mx,my)               &
                            +          gg(iv-2,im,mx,my)               &
                          ) * cef2 * vl(iv)                            &
                        + ( -          gg(iv,im+2,mx,my)               &
                            + 16._DP * gg(iv,im+1,mx,my)               &
                            - 30._DP * gg(iv,im  ,mx,my)               &
                            + 16._DP * gg(iv,im-1,mx,my)               &
                            -          gg(iv,im-2,mx,my)               &
                          ) * cef3                                     &
                        + ( -          gg(iv,im+2,mx,my)               &
                            +  8._DP * gg(iv,im+1,mx,my)               &
                            -  8._DP * gg(iv,im-1,mx,my)               &
                            +          gg(iv,im-2,mx,my)               &
                          ) * cef4 * ( vp(iz,im) + 1._DP / vp(iz,im) ) &
                        + nu_s * 3._DP * gg(iv,im,mx,my)               &      
                        - cflr * ksq(mx,my,iz) * gg(iv,im,mx,my)
              end do
            end do

          end do
        end do
!$OMP end do nowait

      else ! rankm == 0

!$OMP do collapse(2) schedule(dynamic,nchunk_xy)
        do my = ist_y, iend_y
          do mx = -nx, nx

            im = 0
              do iv = 1, 2*nv
                cf(mx,my,iz,iv,im) =                                   &
                          ( -          gg(iv+2,im,mx,my)               &
                            + 16._DP * gg(iv+1,im,mx,my)               &
                            - 30._DP * gg(iv  ,im,mx,my)               &
                            + 16._DP * gg(iv-1,im,mx,my)               &
                            -          gg(iv-2,im,mx,my)               &
                          ) * cef1                                     &
                        + ( -          gg(iv+2,im,mx,my)               &
                            +  8._DP * gg(iv+1,im,mx,my)               &
                            -  8._DP * gg(iv-1,im,mx,my)               &
                            +          gg(iv-2,im,mx,my)               &
                          ) * cef2 * vl(iv)                            &
                        + ( -          gg(iv,im+2,mx,my)               &
                            + 16._DP * gg(iv,im+1,mx,my)               &
                            - 30._DP * gg(iv,im  ,mx,my)               &
                            + 16._DP * gg(iv,im+1,mx,my)               &
                            -          gg(iv,im+2,mx,my)               &
                          ) * cef3 * 2._DP                             &
                        + nu_s * 3._DP * gg(iv,im,mx,my)               &
                        - cflr * ksq(mx,my,iz) * gg(iv,im,mx,my)
              end do

            im = 1
              do iv = 1, 2*nv
                cf(mx,my,iz,iv,im) =                                   &
                          ( -          gg(iv+2,im,mx,my)               &
                            + 16._DP * gg(iv+1,im,mx,my)               &
                            - 30._DP * gg(iv  ,im,mx,my)               &
                            + 16._DP * gg(iv-1,im,mx,my)               &
                            -          gg(iv-2,im,mx,my)               &
                          ) * cef1                                     &
                        + ( -          gg(iv+2,im,mx,my)               &
                            +  8._DP * gg(iv+1,im,mx,my)               &
                            -  8._DP * gg(iv-1,im,mx,my)               &
                            +          gg(iv-2,im,mx,my)               &
                          ) * cef2 * vl(iv)                            &
                        + ( -          gg(iv,im+2,mx,my)               &
                            + 16._DP * gg(iv,im+1,mx,my)               &
                            - 30._DP * gg(iv,im  ,mx,my)               &
                            + 16._DP * gg(iv,im-1,mx,my)               &
                            -          gg(iv,im  ,mx,my)               &
                          ) * cef3                                     &
                        + ( -          gg(iv,im+2,mx,my)               &
                            +  8._DP * gg(iv,im+1,mx,my)               &
                            -  8._DP * gg(iv,im-1,mx,my)               &
                            +          gg(iv,im  ,mx,my)               &
                          ) * cef4 * ( vp(iz,im) + 1._DP / vp(iz,im) ) &
                        + nu_s * 3._DP * gg(iv,im,mx,my)               &   
                        - cflr * ksq(mx,my,iz) * gg(iv,im,mx,my)
              end do

            do im = 2, nm
              do iv = 1, 2*nv
                cf(mx,my,iz,iv,im) =                                   &
                          ( -          gg(iv+2,im,mx,my)               &
                            + 16._DP * gg(iv+1,im,mx,my)               &
                            - 30._DP * gg(iv  ,im,mx,my)               &
                            + 16._DP * gg(iv-1,im,mx,my)               &
                            -          gg(iv-2,im,mx,my)               &
                          ) * cef1                                     &
                        + ( -          gg(iv+2,im,mx,my)               &
                            +  8._DP * gg(iv+1,im,mx,my)               &
                            -  8._DP * gg(iv-1,im,mx,my)               &
                            +          gg(iv-2,im,mx,my)               &
                          ) * cef2 * vl(iv)                            &
                        + ( -          gg(iv,im+2,mx,my)               &
                            + 16._DP * gg(iv,im+1,mx,my)               &
                            - 30._DP * gg(iv,im  ,mx,my)               &
                            + 16._DP * gg(iv,im-1,mx,my)               &
                            -          gg(iv,im-2,mx,my)               &
                          ) * cef3                                     &
                        + ( -          gg(iv,im+2,mx,my)               &
                            +  8._DP * gg(iv,im+1,mx,my)               &
                            -  8._DP * gg(iv,im-1,mx,my)               &
                            +          gg(iv,im-2,mx,my)               &
                          ) * cef4 * ( vp(iz,im) + 1._DP / vp(iz,im) ) &
                        + nu_s * 3._DP * gg(iv,im,mx,my)               &
                        - cflr * ksq(mx,my,iz) * gg(iv,im,mx,my)
              end do
            end do

          end do
        end do
!$OMP end do nowait

      end if

!$OMP master
                                    ! call fapp_stop("literm_colli",1311,1)
                                      call clock_end(1311)
!$OMP end master


  END SUBROUTINE colli_LB_calc



!!--------------------------------------
!  SUBROUTINE colli_LB( ff, phi, cf )
!!--------------------------------------
!!   Lenard-Bernstein model collsion operator
!
!    complex(kind=DP), intent(inout), &
!      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
!    complex(kind=DP), intent(in), &
!      dimension(-nx:nx,0:ny,-nz:nz-1) :: phi
!    complex(kind=DP), intent(out), &
!      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: cf
!
!    complex(kind=DP), dimension(:,:,:,:), allocatable :: vb1e, vb2e, vb1o, vb2o
!    complex(kind=DP), dimension(:,:,:,:), allocatable :: mb1e, mb2e, mb1o, mb2o
!    integer :: iz
!
!      allocate( vb1e(-nx:nx,0:ny,0:nm,1:2*nvb) )
!      allocate( vb2e(-nx:nx,0:ny,0:nm,1:2*nvb) )
!      allocate( mb1e(-nx:nx,0:ny,1:2*nv,1:2*nvb) )
!      allocate( mb2e(-nx:nx,0:ny,1:2*nv,1:2*nvb) )
!      allocate( vb1o(-nx:nx,0:ny,0:nm,1:2*nvb) )
!      allocate( vb2o(-nx:nx,0:ny,0:nm,1:2*nvb) )
!      allocate( mb1o(-nx:nx,0:ny,1:2*nv,1:2*nvb) )
!      allocate( mb2o(-nx:nx,0:ny,1:2*nv,1:2*nvb) )
!
!!$OMP parallel default(none) &
!!$OMP shared(ff,cf,vb1e,vb2e,mb1e,mb2e) &
!!$OMP shared(vb1o,vb2o,mb1o,mb2o) &
!!$OMP private(iz)
!
!!!%%% Without overlap %%%
!!      do iz = -nz, nz-1
!!        call bndry_vm_buffin(iz, ff, vb1e, mb1e)
!!!$OMP barrier
!!!$OMP master
!!        call bndry_vm_sendrecv(vb1e, vb2e, mb1e, mb2e)
!!!$OMP end master
!!!$OMP barrier
!!        call bndry_vm_buffout(iz, vb2e, mb2e, ff)
!!!$OMP barrier
!!        call colli_LB_model_rev(iz, ff, cf)
!!!$OMP barrier
!!      end do
!!!%%%%%%%%%%%%%%%%%%%%%%%
!
!
!!%%% With overlap %%%
!      do iz = -nz, nz-1+3
!        if (mod(iz+nz,2) == 0) then ! even
!!$OMP master
!          if (-nz+1<=iz .and. iz<=nz-1+1) call bndry_vm_sendrecv(vb1o, vb2o, mb1o, mb2o)
!!$OMP end master
!          if (-nz  <=iz .and. iz<=nz-1  ) call bndry_vm_buffin(iz, ff, vb1e, mb1e)
!          if (-nz+2<=iz .and. iz<=nz-1+2) call bndry_vm_buffout(iz-2, vb2e, mb2e, ff)
!          if (-nz+3<=iz .and. iz<=nz-1+3) call colli_LB_model_rev(iz-3, ff, cf)
!        else                        ! odd
!!$OMP master
!          if (-nz+1<=iz .and. iz<=nz-1+1) call bndry_vm_sendrecv(vb1e, vb2e, mb1e, mb2e)
!!$OMP end master
!          if (-nz  <=iz .and. iz<=nz-1  ) call bndry_vm_buffin(iz, ff, vb1o, mb1o)
!          if (-nz+2<=iz .and. iz<=nz-1+2) call bndry_vm_buffout(iz-2, vb2o, mb2o, ff)
!          if (-nz+3<=iz .and. iz<=nz-1+3) call colli_LB_model_rev(iz-3, ff, cf)
!        end if
!!$OMP barrier
!      end do
!!%%%%%%%%%%%%%%%%%%%%
!
!!$OMP end parallel
!
!      deallocate( vb1e )
!      deallocate( vb2e )
!      deallocate( mb1e )
!      deallocate( mb2e )
!      deallocate( vb1o )
!      deallocate( vb2o )
!      deallocate( mb1o )
!      deallocate( mb2o )
!
!  END SUBROUTINE colli_LB
!
!
!!--------------------------------------
!  SUBROUTINE colli_LB_model_rev( iz, ff, cf )
!!--------------------------------------
!!   Lenard-Bernstein model collsion operator
!
!    integer, intent(in) :: iz
!    complex(kind=DP), intent(in), &
!      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
!    complex(kind=DP), intent(inout), &
!      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: cf
!
!    real(kind=DP) :: nu_s, cef1, cef2, cef3, cef4
!    integer  ::  mx, my, iv, im
!
!
!!$OMP master
!                                           call clock_sta(1311)
!                                         ! call fapp_start("literm_colli",1311,1)
!!$OMP end master
!
!! --- Note that nu(ranks) is a bias factor 
!      nu_s = nu(ranks)*3._DP*dsqrt(pi)*ctauiv(ranks,ranks)/4._DP
!!      nu_s = 1.d-3
!
!      cef1   = nu_s / ( 12._DP * dv * dv )
!      cef2   = nu_s / ( 12._DP * dv )
!      cef3   = nu_s / ( 12._DP * dvp(iz) * dvp(iz) )
!      cef4   = nu_s / ( 12._DP * dvp(iz) )
!
!      if( rankm /= 0  ) then
!
!        do im = 0, nm
!!$OMP do schedule(dynamic)
!          do iv = 1, 2*nv
!              do my = ist_y, iend_y
!                do mx = -nx, nx
!                  cf(mx,my,iz,iv,im) =                                        &
!                            ( -          ff(mx,my,iz,iv+2,im)                 &
!                              + 16._DP * ff(mx,my,iz,iv+1,im)                 &
!                              - 30._DP * ff(mx,my,iz,iv  ,im)                 &
!                              + 16._DP * ff(mx,my,iz,iv-1,im)                 &
!                              -          ff(mx,my,iz,iv-2,im)                 &
!                             ) * cef1                                         &
!                           + ( -          ff(mx,my,iz,iv+2,im)                &
!                               +  8._DP * ff(mx,my,iz,iv+1,im)                &
!                               -  8._DP * ff(mx,my,iz,iv-1,im)                &
!                               +          ff(mx,my,iz,iv-2,im)                &
!                             ) * cef2 * vl(iv)                                &
!                           + ( -          ff(mx,my,iz,iv,im+2)                &
!                               + 16._DP * ff(mx,my,iz,iv,im+1)                &
!                               - 30._DP * ff(mx,my,iz,iv,im  )                &
!                               + 16._DP * ff(mx,my,iz,iv,im-1)                &
!                               -          ff(mx,my,iz,iv,im-2)                &
!                             ) * cef3                                         &
!                           + ( -          ff(mx,my,iz,iv,im+2)                &
!                               +  8._DP * ff(mx,my,iz,iv,im+1)                &
!                               -  8._DP * ff(mx,my,iz,iv,im-1)                &
!                               +          ff(mx,my,iz,iv,im-2)                &
!                             ) * cef4 * ( vp(iz,im) + 1._DP / vp(iz,im) )     &
!                           + nu_s * 3._DP * ff(mx,my,iz,iv,im)                &      
!                           - nu_s * ksq(mx,my,iz) * Anum(ranks) * tau(ranks)  &
!                             / ( Znum(ranks) * omg(iz) )**2 * ff(mx,my,iz,iv,im) &
!                               * real(iFLR, kind=DP)
!                end do
!              end do
!          end do
!!$OMP end do nowait
!        end do
!
!      else
!
!        im = 0
!!$OMP do schedule(dynamic)
!          do iv = 1, 2*nv
!              do my = ist_y, iend_y
!                do mx = -nx, nx
!                  cf(mx,my,iz,iv,im) =                                        &
!                             ( -          ff(mx,my,iz,iv+2,im)                &
!                               + 16._DP * ff(mx,my,iz,iv+1,im)                &
!                               - 30._DP * ff(mx,my,iz,iv  ,im)                &
!                               + 16._DP * ff(mx,my,iz,iv-1,im)                &
!                               -          ff(mx,my,iz,iv-2,im)                &
!                             ) * cef1                                         &
!                           + ( -          ff(mx,my,iz,iv+2,im)                &
!                               +  8._DP * ff(mx,my,iz,iv+1,im)                &
!                               -  8._DP * ff(mx,my,iz,iv-1,im)                &
!                               +          ff(mx,my,iz,iv-2,im)                &
!                             ) * cef2 * vl(iv)                                &
!                           + ( -          ff(mx,my,iz,iv,im+2)                &
!                               + 16._DP * ff(mx,my,iz,iv,im+1)                &
!                               - 30._DP * ff(mx,my,iz,iv,im  )                &
!                               + 16._DP * ff(mx,my,iz,iv,im+1)                &
!                               -          ff(mx,my,iz,iv,im+2)                &
!                             ) * cef3 * 2._DP                                 &
!                           + nu_s * 3._DP * ff(mx,my,iz,iv,im)                &
!                           - nu_s * ksq(mx,my,iz) * Anum(ranks) * tau(ranks)  &
!                             / ( Znum(ranks) * omg(iz) )**2 * ff(mx,my,iz,iv,im) &
!                             * real(iFLR, kind=DP)
!                end do
!              end do
!          end do
!!$OMP end do nowait
!
!        im = 1
!
!!$OMP do schedule(dynamic)
!          do iv = 1, 2*nv
!              do my = ist_y, iend_y
!                do mx = -nx, nx
!                  cf(mx,my,iz,iv,im) =                                        &
!                             ( -          ff(mx,my,iz,iv+2,im)                &
!                               + 16._DP * ff(mx,my,iz,iv+1,im)                &
!                               - 30._DP * ff(mx,my,iz,iv  ,im)                &
!                               + 16._DP * ff(mx,my,iz,iv-1,im)                &
!                               -          ff(mx,my,iz,iv-2,im)                &
!                             ) * cef1                                         &
!                           + ( -          ff(mx,my,iz,iv+2,im)                &
!                               +  8._DP * ff(mx,my,iz,iv+1,im)                &
!                               -  8._DP * ff(mx,my,iz,iv-1,im)                &
!                               +          ff(mx,my,iz,iv-2,im)                &
!                             ) * cef2 * vl(iv)                                &
!                           + ( -          ff(mx,my,iz,iv,im+2)                &
!                               + 16._DP * ff(mx,my,iz,iv,im+1)                &
!                               - 30._DP * ff(mx,my,iz,iv,im  )                &
!                               + 16._DP * ff(mx,my,iz,iv,im-1)                &
!                               -          ff(mx,my,iz,iv,im  )                &
!                             ) * cef3                                         &
!                           + ( -          ff(mx,my,iz,iv,im+2)                &
!                               +  8._DP * ff(mx,my,iz,iv,im+1)                &
!                               -  8._DP * ff(mx,my,iz,iv,im-1)                &
!                               +          ff(mx,my,iz,iv,im  )                &
!                             ) * cef4 * ( vp(iz,im) + 1._DP / vp(iz,im) )     &
!                           + nu_s * 3._DP * ff(mx,my,iz,iv,im)                &   
!                           - nu_s * ksq(mx,my,iz) * Anum(ranks) * tau(ranks)  &
!                             / ( Znum(ranks) * omg(iz) )**2 * ff(mx,my,iz,iv,im) &
!                             * real(iFLR, kind=DP)
!                end do
!              end do
!          end do
!!$OMP end do nowait
!
!        do im = 2, nm
!!$OMP do schedule(dynamic)
!          do iv = 1, 2*nv
!              do my = ist_y, iend_y
!                do mx = -nx, nx
!                  cf(mx,my,iz,iv,im) =                                        &
!                             ( -          ff(mx,my,iz,iv+2,im)                &
!                               + 16._DP * ff(mx,my,iz,iv+1,im)                &
!                               - 30._DP * ff(mx,my,iz,iv  ,im)                &
!                               + 16._DP * ff(mx,my,iz,iv-1,im)                &
!                               -          ff(mx,my,iz,iv-2,im)                &
!                             ) * cef1                                         &
!                           + ( -          ff(mx,my,iz,iv+2,im)                &
!                               +  8._DP * ff(mx,my,iz,iv+1,im)                &
!                               -  8._DP * ff(mx,my,iz,iv-1,im)                &
!                               +          ff(mx,my,iz,iv-2,im)                &
!                             ) * cef2 * vl(iv)                                &
!                           + ( -          ff(mx,my,iz,iv,im+2)                &
!                               + 16._DP * ff(mx,my,iz,iv,im+1)                &
!                               - 30._DP * ff(mx,my,iz,iv,im  )                &
!                               + 16._DP * ff(mx,my,iz,iv,im-1)                &
!                               -          ff(mx,my,iz,iv,im-2)                &
!                             ) * cef3                                         &
!                           + ( -          ff(mx,my,iz,iv,im+2)                &
!                               +  8._DP * ff(mx,my,iz,iv,im+1)                &
!                               -  8._DP * ff(mx,my,iz,iv,im-1)                &
!                               +          ff(mx,my,iz,iv,im-2)                &
!                             ) * cef4 * ( vp(iz,im) + 1._DP / vp(iz,im) )     &
!                           + nu_s * 3._DP * ff(mx,my,iz,iv,im)                &
!                           - nu_s * ksq(mx,my,iz) * Anum(ranks) * tau(ranks)  &
!                             / ( Znum(ranks) * omg(iz) )**2 * ff(mx,my,iz,iv,im) & 
!                             * real(iFLR, kind=DP)
!                end do
!              end do
!          end do
!!$OMP end do nowait
!        end do
!
!      end if
!
!!$OMP master
!                                    ! call fapp_stop("literm_colli",1311,1)
!                                      call clock_end(1311)
!!$OMP end master
!
!
!  END SUBROUTINE colli_LB_model_rev




!!!!%%% Subroutines for colli_LB (old version) %%%
!!!!--------------------------------------
!!!  SUBROUTINE colli_LB( ff, phi, cf )
!!!!--------------------------------------
!!!!   Lenard-Bernstein model collsion operator
!!!
!!!    complex(kind=DP), intent(inout), &
!!!      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
!!!    complex(kind=DP), intent(in), &
!!!      dimension(-nx:nx,0:ny,-nz:nz-1) :: phi
!!!    complex(kind=DP), intent(out), &
!!!      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: cf
!!!
!!!    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: vb1, vb2
!!!    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: mb1, mb2
!!!    integer :: im
!!!
!!!      allocate( vb1(-nx:nx,0:ny,-nz:nz-1,1:2*nvb,0:nm) )
!!!      allocate( vb2(-nx:nx,0:ny,-nz:nz-1,1:2*nvb,0:nm) )
!!!      allocate( mb1(-nx:nx,0:ny,-nz:nz-1,1:2*nv,1:2*nvb) )
!!!      allocate( mb2(-nx:nx,0:ny,-nz:nz-1,1:2*nv,1:2*nvb) )
!!!
!!!!$OMP parallel default(none) &
!!!!$OMP shared(ff,cf,vb1,vb2,mb1,mb2) &
!!!!$OMP private(im)
!!!
!!!      call bndry_shifts_m_buffin ( ff, mb1, mb2 )
!!!!$OMP barrier
!!!
!!!!$OMP master
!!!      call bndry_shifts_m_sendrecv ( mb1, mb2 )
!!!!$OMP end master
!!!      do im = 0, nm
!!!        call bndry_shifts_v_buffin ( ff(:,:,:,:,im), vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
!!!      end do
!!!      call colli_zeroset( cf )
!!!!$OMP barrier
!!!
!!!      im = 0
!!!!$OMP master
!!!        call bndry_shifts_v_sendrecv ( vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
!!!!$OMP end master
!!!        call bndry_shifts_m_buffout ( mb2, ff )
!!!!$OMP barrier
!!!
!!!      im = 1
!!!!$OMP master
!!!        call bndry_shifts_v_sendrecv ( vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
!!!!$OMP end master
!!!        call bndry_shifts_v_buffout ( vb2(:,:,:,:,im-1), ff(:,:,:,:,im-1) )
!!!!$OMP barrier
!!!
!!!      im = 2
!!!!$OMP master
!!!        call bndry_shifts_v_sendrecv ( vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
!!!!$OMP end master
!!!        call bndry_shifts_v_buffout ( vb2(:,:,:,:,im-1), ff(:,:,:,:,im-1) )
!!!        call colli_LB_model( ff, im-2, cf(:,:,:,:,im-2) )
!!!!$OMP barrier
!!!
!!!      do im = 3, nm
!!!!$OMP master
!!!        call bndry_shifts_v_sendrecv ( vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
!!!!$OMP end master
!!!        call bndry_shifts_v_buffout ( vb2(:,:,:,:,im-1), ff(:,:,:,:,im-1) )
!!!        call colli_LB_model( ff, im-2, cf(:,:,:,:,im-2) )
!!!!$OMP barrier
!!!      end do
!!!
!!!      call bndry_shifts_v_buffout ( vb2(:,:,:,:,nm), ff(:,:,:,:,nm) )
!!!      call colli_LB_model( ff, nm-1, cf(:,:,:,:,nm-1) )
!!!!$OMP barrier
!!!
!!!      call colli_LB_model( ff, nm, cf(:,:,:,:,nm) )
!!!!$OMP barrier
!!!
!!!!$OMP end parallel
!!!
!!!      deallocate( vb1 )
!!!      deallocate( vb2 )
!!!      deallocate( mb1 )
!!!      deallocate( mb2 )
!!!
!!!  END SUBROUTINE colli_LB
!!!
!!!
!!!!--------------------------------------
!!!  SUBROUTINE colli_LB_model( ff, im, cf )
!!!!--------------------------------------
!!!!   Lenard-Bernstein model collsion operator
!!!
!!!    complex(kind=DP), intent(in), &
!!!      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
!!!    integer, intent(in) :: im
!!!    complex(kind=DP), intent(out), &
!!!      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv)           :: cf
!!!
!!!    real(kind=DP) :: nu_s, cef1, cef2
!!!    real(kind=DP), dimension(-nz:nz-1) :: cef3, cef4
!!!    integer  ::  mx, my, iz, iv
!!!
!!!
!!!!$OMP master
!!!                                           call clock_sta(1311)
!!!                                         ! call fapp_start("literm_colli",1311,1)
!!!!$OMP end master
!!!
!!!! --- Note that nu(ranks) is a bias factor 
!!!      nu_s = nu(ranks)*3._DP*dsqrt(pi)*ctauiv(ranks,ranks)/4._DP
!!!!      nu_s = 1.d-3
!!!
!!!      cef1   = nu_s / ( 12._DP * dv * dv )
!!!      cef2   = nu_s / ( 12._DP * dv )
!!!
!!!      do iz = -nz, nz-1
!!!        cef3(iz)   = nu_s / ( 12._DP * dvp(iz) * dvp(iz) )
!!!        cef4(iz)   = nu_s / ( 12._DP * dvp(iz) )
!!!      end do
!!!
!!!      if( rankm /= 0  ) then
!!!
!!!          do iv = 1, 2*nv
!!!!$OMP do schedule (dynamic)
!!!            do iz = -nz, nz-1
!!!              do my = ist_y, iend_y
!!!                do mx = -nx, nx
!!!                  cf(mx,my,iz,iv) =                                           &
!!!                            ( -          ff(mx,my,iz,iv+2,im)                 &
!!!                              + 16._DP * ff(mx,my,iz,iv+1,im)                 &
!!!                              - 30._DP * ff(mx,my,iz,iv  ,im)                 &
!!!                              + 16._DP * ff(mx,my,iz,iv-1,im)                 &
!!!                              -          ff(mx,my,iz,iv-2,im)                 &
!!!                             ) * cef1                                         &
!!!                           + ( -          ff(mx,my,iz,iv+2,im)                &
!!!                               +  8._DP * ff(mx,my,iz,iv+1,im)                &
!!!                               -  8._DP * ff(mx,my,iz,iv-1,im)                &
!!!                               +          ff(mx,my,iz,iv-2,im)                &
!!!                             ) * cef2 * vl(iv)                                &
!!!                           + ( -          ff(mx,my,iz,iv,im+2)                &
!!!                               + 16._DP * ff(mx,my,iz,iv,im+1)                &
!!!                               - 30._DP * ff(mx,my,iz,iv,im  )                &
!!!                               + 16._DP * ff(mx,my,iz,iv,im-1)                &
!!!                               -          ff(mx,my,iz,iv,im-2)                &
!!!                             ) * cef3(iz)                                     &
!!!                           + ( -          ff(mx,my,iz,iv,im+2)                &
!!!                               +  8._DP * ff(mx,my,iz,iv,im+1)                &
!!!                               -  8._DP * ff(mx,my,iz,iv,im-1)                &
!!!                               +          ff(mx,my,iz,iv,im-2)                &
!!!                             ) * cef4(iz) * ( vp(iz,im) + 1._DP / vp(iz,im) ) &
!!!                           + nu_s * 3._DP * ff(mx,my,iz,iv,im)                &      
!!!                           - nu_s * ksq(mx,my,iz) * Anum(ranks) * tau(ranks)  &
!!!                             / ( Znum(ranks) * omg(iz) )**2 * ff(mx,my,iz,iv,im) &
!!!                               * real(iFLR, kind=DP)
!!!                end do
!!!              end do
!!!            end do
!!!!$OMP end do nowait
!!!          end do
!!!
!!!      else
!!!
!!!          if ( im == 0 ) then
!!!
!!!            do iv = 1, 2*nv
!!!!$OMP do schedule (dynamic)
!!!              do iz = -nz, nz-1
!!!                do my = ist_y, iend_y
!!!                  do mx = -nx, nx
!!!                    cf(mx,my,iz,iv) =                                           &
!!!                               ( -          ff(mx,my,iz,iv+2,im)                &
!!!                                 + 16._DP * ff(mx,my,iz,iv+1,im)                &
!!!                                 - 30._DP * ff(mx,my,iz,iv  ,im)                &
!!!                                 + 16._DP * ff(mx,my,iz,iv-1,im)                &
!!!                                 -          ff(mx,my,iz,iv-2,im)                &
!!!                               ) * cef1                                         &
!!!                             + ( -          ff(mx,my,iz,iv+2,im)                &
!!!                                 +  8._DP * ff(mx,my,iz,iv+1,im)                &
!!!                                 -  8._DP * ff(mx,my,iz,iv-1,im)                &
!!!                                 +          ff(mx,my,iz,iv-2,im)                &
!!!                               ) * cef2 * vl(iv)                                &
!!!                             + ( -          ff(mx,my,iz,iv,im+2)                &
!!!                                 + 16._DP * ff(mx,my,iz,iv,im+1)                &
!!!                                 - 30._DP * ff(mx,my,iz,iv,im  )                &
!!!                                 + 16._DP * ff(mx,my,iz,iv,im+1)                &
!!!                                 -          ff(mx,my,iz,iv,im+2)                &
!!!                               ) * cef3(iz) * 2._DP                             &
!!!                             + nu_s * 3._DP * ff(mx,my,iz,iv,im)                &
!!!                             - nu_s * ksq(mx,my,iz) * Anum(ranks) * tau(ranks)  &
!!!                               / ( Znum(ranks) * omg(iz) )**2 * ff(mx,my,iz,iv,im) &
!!!                               * real(iFLR, kind=DP)
!!!                  end do
!!!                end do
!!!              end do
!!!!$OMP end do nowait
!!!            end do
!!!
!!!          else if ( im == 1 ) then
!!!
!!!            do iv = 1, 2*nv
!!!!$OMP do schedule (dynamic)
!!!              do iz = -nz, nz-1
!!!                do my = ist_y, iend_y
!!!                  do mx = -nx, nx
!!!                    cf(mx,my,iz,iv) =                                           &
!!!                               ( -          ff(mx,my,iz,iv+2,im)                &
!!!                                 + 16._DP * ff(mx,my,iz,iv+1,im)                &
!!!                                 - 30._DP * ff(mx,my,iz,iv  ,im)                &
!!!                                 + 16._DP * ff(mx,my,iz,iv-1,im)                &
!!!                                 -          ff(mx,my,iz,iv-2,im)                &
!!!                               ) * cef1                                         &
!!!                             + ( -          ff(mx,my,iz,iv+2,im)                &
!!!                                 +  8._DP * ff(mx,my,iz,iv+1,im)                &
!!!                                 -  8._DP * ff(mx,my,iz,iv-1,im)                &
!!!                                 +          ff(mx,my,iz,iv-2,im)                &
!!!                               ) * cef2 * vl(iv)                                &
!!!                             + ( -          ff(mx,my,iz,iv,im+2)                &
!!!                                 + 16._DP * ff(mx,my,iz,iv,im+1)                &
!!!                                 - 30._DP * ff(mx,my,iz,iv,im  )                &
!!!                                 + 16._DP * ff(mx,my,iz,iv,im-1)                &
!!!                                 -          ff(mx,my,iz,iv,im  )                &
!!!                               ) * cef3(iz)                                     &
!!!                             + ( -          ff(mx,my,iz,iv,im+2)                &
!!!                                 +  8._DP * ff(mx,my,iz,iv,im+1)                &
!!!                                 -  8._DP * ff(mx,my,iz,iv,im-1)                &
!!!                                 +          ff(mx,my,iz,iv,im  )                &
!!!                               ) * cef4(iz) * ( vp(iz,im) + 1._DP / vp(iz,im) ) &
!!!                             + nu_s * 3._DP * ff(mx,my,iz,iv,im)                &   
!!!                             - nu_s * ksq(mx,my,iz) * Anum(ranks) * tau(ranks)  &
!!!                               / ( Znum(ranks) * omg(iz) )**2 * ff(mx,my,iz,iv,im) &
!!!                               * real(iFLR, kind=DP)
!!!                  end do
!!!                end do
!!!              end do
!!!!$OMP end do nowait
!!!            end do
!!!
!!!          else  ! 2=<im=<nm
!!!
!!!            do iv = 1, 2*nv
!!!!$OMP do schedule (dynamic)
!!!              do iz = -nz, nz-1
!!!                do my = ist_y, iend_y
!!!                  do mx = -nx, nx
!!!                    cf(mx,my,iz,iv) =                                           &
!!!                               ( -          ff(mx,my,iz,iv+2,im)                &
!!!                                 + 16._DP * ff(mx,my,iz,iv+1,im)                &
!!!                                 - 30._DP * ff(mx,my,iz,iv  ,im)                &
!!!                                 + 16._DP * ff(mx,my,iz,iv-1,im)                &
!!!                                 -          ff(mx,my,iz,iv-2,im)                &
!!!                               ) * cef1                                         &
!!!                             + ( -          ff(mx,my,iz,iv+2,im)                &
!!!                                 +  8._DP * ff(mx,my,iz,iv+1,im)                &
!!!                                 -  8._DP * ff(mx,my,iz,iv-1,im)                &
!!!                                 +          ff(mx,my,iz,iv-2,im)                &
!!!                               ) * cef2 * vl(iv)                                &
!!!                             + ( -          ff(mx,my,iz,iv,im+2)                &
!!!                                 + 16._DP * ff(mx,my,iz,iv,im+1)                &
!!!                                 - 30._DP * ff(mx,my,iz,iv,im  )                &
!!!                                 + 16._DP * ff(mx,my,iz,iv,im-1)                &
!!!                                 -          ff(mx,my,iz,iv,im-2)                &
!!!                               ) * cef3(iz)                                     &
!!!                             + ( -          ff(mx,my,iz,iv,im+2)                &
!!!                                 +  8._DP * ff(mx,my,iz,iv,im+1)                &
!!!                                 -  8._DP * ff(mx,my,iz,iv,im-1)                &
!!!                                 +          ff(mx,my,iz,iv,im-2)                &
!!!                               ) * cef4(iz) * ( vp(iz,im) + 1._DP / vp(iz,im) ) &
!!!                             + nu_s * 3._DP * ff(mx,my,iz,iv,im)                &
!!!                             - nu_s * ksq(mx,my,iz) * Anum(ranks) * tau(ranks)  &
!!!                               / ( Znum(ranks) * omg(iz) )**2 * ff(mx,my,iz,iv,im) & 
!!!                               * real(iFLR, kind=DP)
!!!                  end do
!!!                end do
!!!              end do
!!!!$OMP end do nowait
!!!            end do
!!!
!!!          end if
!!!
!!!      end if
!!!
!!!!$OMP master
!!!                                    ! call fapp_stop("literm_colli",1311,1)
!!!                                      call clock_end(1311)
!!!!$OMP end master
!!!
!!!
!!!  END SUBROUTINE colli_LB_model
!!!
!!!
!!!!--------------------------------------
!!!  SUBROUTINE colli_zeroset( cf )
!!!!--------------------------------------
!!!!-------------------------------------------------------------------------------
!!!!
!!!!    zero clear for collision terms 
!!!!
!!!!    by M. Nakata, M. Nunami, April 2014
!!!!
!!!!-------------------------------------------------------------------------------
!!!
!!!    complex(kind=DP), intent(inout), &
!!!      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm)    :: cf
!!!
!!!    integer  ::  mx, my, iz, iv, im
!!!
!!!!$OMP master
!!!                                      call clock_sta(1311)
!!!                                    ! call fapp_start("literm_colli",1311,1)
!!!!$OMP end master
!!!
!!!!$OMP do schedule (dynamic)
!!!        do im = 0, nm
!!!          do iv = 1, 2*nv
!!!            do iz = -nz, nz-1
!!!              do my = ist_y, iend_y
!!!                do mx = -nx, nx
!!!                  cf(mx,my,iz,iv,im) = ( 0._DP, 0._DP )     
!!!                end do
!!!              end do
!!!            end do
!!!          end do
!!!        end do
!!!!$OMP end do nowait 
!!!
!!!
!!!!$OMP master
!!!                                    ! call fapp_stop("literm_colli",1311,1)
!!!                                      call clock_end(1311)
!!!!$OMP end master
!!!
!!!
!!!  END SUBROUTINE colli_zeroset


!!!%%% Subroutines for colli_full %%%
!!!!--------------------------------------
!!!  SUBROUTINE colli_GK_CT( ff, phi, dfdvp, im, cf )
!!!!--------------------------------------
!!!!-------------------------------------------------------------------------------
!!!!
!!!!    Differential and FLR terms of test particle part in gyrokinetic collision
!!!!                                                         with 4th order CFD
!!!!    by M. Nakata, M. Nunami, April 2014
!!!!
!!!!-------------------------------------------------------------------------------
!!!
!!!    complex(kind=DP), intent(in), &
!!!      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
!!!    complex(kind=DP), intent(in), &
!!!      dimension(-nx:nx,0:ny,-nz:nz-1)                       :: phi
!!!    integer, intent(in) :: im
!!!    complex(kind=DP), intent(out), &
!!!      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv)           :: cf
!!!
!!!    complex(kind=DP), intent(in), &
!!!      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb) :: dfdvp
!!!
!!!    real(kind=DP) :: cef1, cef2
!!!    real(kind=DP), dimension(-nz:nz-1) :: cef3
!!!    integer  ::  mx, my, iz, iv, is1
!!!
!!!
!!!!$OMP master
!!!                                      call clock_sta(1311)
!!!                                    ! call fapp_start("literm_colli_ct",1311,1)
!!!!$OMP end master
!!!
!!!
!!!      cef1   = 1._DP / ( 12._DP * dv * dv )
!!!      cef2   = 1._DP / ( 12._DP * dv )
!!!
!!!      do iz = -nz, nz-1
!!!        cef3(iz)   = 1._DP / ( 12._DP * dvp(iz) * dvp(iz) )
!!!      end do
!!!
!!!      if ( rankm == 0 .AND. im == 0 ) then
!!!
!!!        is1 = ranks
!!!!$OMP do schedule (dynamic)
!!!          do iv = 1, 2*nv
!!!            do iz = -nz, nz-1
!!!              do my = ist_y, iend_y
!!!                do mx = -nx, nx
!!!                  cf(mx,my,iz,iv) =                                                &
!!!                             ( -          ff(mx,my,iz,iv-2,im)                     &
!!!                               + 16._DP * ff(mx,my,iz,iv-1,im)                     &
!!!                               - 30._DP * ff(mx,my,iz,iv  ,im)                     &
!!!                               + 16._DP * ff(mx,my,iz,iv+1,im)                     &
!!!                               -          ff(mx,my,iz,iv+2,im)                     &
!!!                             ) * cef1                                              &
!!!                               * (   nu_ps(iz,iv,im)*vl(iv)**2                     & 
!!!                                 ) * 0.5_DP                                        &          
!!!                           + ( - 30._DP * ff(mx,my,iz,iv,im  )                     &
!!!                               + 32._DP * ff(mx,my,iz,iv,im+1)                     &
!!!                               -  2._DP * ff(mx,my,iz,iv,im+2)                     &
!!!                             ) * cef3(iz)                                          &
!!!                               * (   nu_ds(iz,iv,im)*vl(iv)**2                     & 
!!!                                 )                                                 &           
!!!                           + ( +          ff(mx,my,iz,iv-2,im)                     &
!!!                               -  8._DP * ff(mx,my,iz,iv-1,im)                     &
!!!                               +  8._DP * ff(mx,my,iz,iv+1,im)                     &
!!!                               -          ff(mx,my,iz,iv+2,im)                     &
!!!                             ) * cef2 * nu_gs(iz,iv,im)*vl(iv)                     & 
!!!                           + ( nu_hs(iz,iv,im)*xxa(iz,iv,im)**2*2._DP              & 
!!!                                - 0.25_DP*ksq(mx,my,iz)*Anum(is1)*tau(is1)         &
!!!                                 * ( nu_ds(iz,iv,im) * 2._DP*vl(iv)**2 )           &
!!!                                 / (Znum(is1)*omg(iz))**2 * real(iFLR, kind=DP)    &
!!!                             ) * ff(mx,my,iz,iv,im)                                &
!!!                           + adbtc(mx,my,iz,iv,im)*phi(mx,my,iz)
!!!                end do
!!!              end do
!!!            end do
!!!          end do
!!!!$OMP end do nowait
!!!
!!!      else if ( rankm == 0 .AND. im == 1 ) then
!!!
!!!        is1 = ranks
!!!!$OMP do schedule (dynamic)
!!!          do iv = 1, 2*nv
!!!            do iz = -nz, nz-1
!!!              do my = ist_y, iend_y
!!!                do mx = -nx, nx
!!!                  cf(mx,my,iz,iv) =                                                &
!!!                             ( -          ff(mx,my,iz,iv-2,im)                     &
!!!                               + 16._DP * ff(mx,my,iz,iv-1,im)                     &
!!!                               - 30._DP * ff(mx,my,iz,iv  ,im)                     &
!!!                               + 16._DP * ff(mx,my,iz,iv+1,im)                     &
!!!                               -          ff(mx,my,iz,iv+2,im)                     &
!!!                             ) * cef1                                              &
!!!                               * (   nu_ps(iz,iv,im)*vl(iv)**2                     & 
!!!                                   + nu_ds(iz,iv,im)*vp(iz,im)**2                  &
!!!                                 ) * 0.5_DP                                        &          
!!!                           + ( + 16._DP * ff(mx,my,iz,iv,im-1)                     &
!!!                               - 31._DP * ff(mx,my,iz,iv,im  )                     &
!!!                               + 16._DP * ff(mx,my,iz,iv,im+1)                     &
!!!                               -          ff(mx,my,iz,iv,im+2)                     &
!!!                             ) * cef3(iz)                                          &
!!!                               * (   nu_ds(iz,iv,im)*vl(iv)**2                     & 
!!!                                   + nu_ps(iz,iv,im)*vp(iz,im)**2                  &
!!!                                 ) * 0.5_DP                                        &           
!!!                           + ( +          dfdvp(mx,my,iz,iv-2)                     &             
!!!                               -  8._DP * dfdvp(mx,my,iz,iv-1)                     &
!!!                               +  8._DP * dfdvp(mx,my,iz,iv+1)                     &
!!!                               -          dfdvp(mx,my,iz,iv+2)                     &
!!!                             ) * cef2 * vl(iv) * vp(iz,im)                         & 
!!!                               * (   nu_ps(iz,iv,im)                               & 
!!!                                   - nu_ds(iz,iv,im) )                             &
!!!                           + ( +          ff(mx,my,iz,iv-2,im)                     &
!!!                               -  8._DP * ff(mx,my,iz,iv-1,im)                     &
!!!                               +  8._DP * ff(mx,my,iz,iv+1,im)                     &
!!!                               -          ff(mx,my,iz,iv+2,im)                     &
!!!                             ) * cef2 * nu_gs(iz,iv,im)*vl(iv)                     & 
!!!                           + ( dfdvp(mx,my,iz,iv) )                                &
!!!                               * (   nu_gs(iz,iv,im)*vp(iz,im)                     & 
!!!                                   + nu_ds(iz,iv,im)*0.5_DP                        &
!!!                                     * (vl(iv)**2/vp(iz,im)+vp(iz,im))             &
!!!                                 )                                                 &
!!!                           + ( nu_hs(iz,iv,im)*xxa(iz,iv,im)**2*2._DP              & 
!!!                                - 0.25_DP*ksq(mx,my,iz)*Anum(is1)*tau(is1)         &
!!!                                 * ( nu_ds(iz,iv,im)                               &
!!!                                      * (2._DP*vl(iv)**2 + vp(iz,im)**2)           &
!!!                                    + nu_ps(iz,iv,im)*vp(iz,im)**2 )               &
!!!                                 / (Znum(is1)*omg(iz))**2 * real(iFLR, kind=DP)    &
!!!                             ) * ff(mx,my,iz,iv,im)                                &
!!!                           + adbtc(mx,my,iz,iv,im)*phi(mx,my,iz)                    
!!!                end do
!!!              end do
!!!            end do
!!!          end do
!!!!$OMP end do nowait
!!!
!!!      else  ! im=[0,nm] for rankm > 0 and im=[2,nm] nm for rankm = 0  
!!!
!!!        is1 = ranks
!!!!$OMP do schedule (dynamic)
!!!          do iv = 1, 2*nv
!!!            do iz = -nz, nz-1
!!!              do my = ist_y, iend_y
!!!                do mx = -nx, nx
!!!                  cf(mx,my,iz,iv) =                                                &
!!!                             ( -          ff(mx,my,iz,iv-2,im)                     &
!!!                               + 16._DP * ff(mx,my,iz,iv-1,im)                     &
!!!                               - 30._DP * ff(mx,my,iz,iv  ,im)                     &
!!!                               + 16._DP * ff(mx,my,iz,iv+1,im)                     &
!!!                               -          ff(mx,my,iz,iv+2,im)                     &
!!!                             ) * cef1                                              &
!!!                               * (   nu_ps(iz,iv,im)*vl(iv)**2                     & 
!!!                                   + nu_ds(iz,iv,im)*vp(iz,im)**2                  &
!!!                                 ) * 0.5_DP                                        &          
!!!                           + ( -          ff(mx,my,iz,iv,im-2)                     &
!!!                               + 16._DP * ff(mx,my,iz,iv,im-1)                     &
!!!                               - 30._DP * ff(mx,my,iz,iv,im  )                     &
!!!                               + 16._DP * ff(mx,my,iz,iv,im+1)                     &
!!!                               -          ff(mx,my,iz,iv,im+2)                     &
!!!                             ) * cef3(iz)                                          &
!!!                               * (   nu_ds(iz,iv,im)*vl(iv)**2                     & 
!!!                                   + nu_ps(iz,iv,im)*vp(iz,im)**2                  &
!!!                                 ) * 0.5_DP                                        &           
!!!                           + ( +          dfdvp(mx,my,iz,iv-2)                     &             
!!!                               -  8._DP * dfdvp(mx,my,iz,iv-1)                     &
!!!                               +  8._DP * dfdvp(mx,my,iz,iv+1)                     &
!!!                               -          dfdvp(mx,my,iz,iv+2)                     &
!!!                             ) * cef2 * vl(iv) * vp(iz,im)                         & 
!!!                               * (   nu_ps(iz,iv,im)                               & 
!!!                                   - nu_ds(iz,iv,im) )                             &
!!!                           + ( +          ff(mx,my,iz,iv-2,im)                     &
!!!                               -  8._DP * ff(mx,my,iz,iv-1,im)                     &
!!!                               +  8._DP * ff(mx,my,iz,iv+1,im)                     &
!!!                               -          ff(mx,my,iz,iv+2,im)                     &
!!!                             ) * cef2 * nu_gs(iz,iv,im)*vl(iv)                     & 
!!!                           + ( dfdvp(mx,my,iz,iv) )                                &
!!!                               * (   nu_gs(iz,iv,im)*vp(iz,im)                     & 
!!!                                   + nu_ds(iz,iv,im)*0.5_DP                        &
!!!                                     * (vl(iv)**2/vp(iz,im)+vp(iz,im))             &
!!!                                 )                                                 &
!!!                           + ( nu_hs(iz,iv,im)*xxa(iz,iv,im)**2*2._DP              & 
!!!                                - 0.25_DP*ksq(mx,my,iz)*Anum(is1)*tau(is1)         &
!!!                                 * ( nu_ds(iz,iv,im)                               &
!!!                                      * (2._DP*vl(iv)**2 + vp(iz,im)**2)           &
!!!                                    + nu_ps(iz,iv,im)*vp(iz,im)**2 )               &
!!!                                 / (Znum(is1)*omg(iz))**2 * real(iFLR, kind=DP)    &
!!!                             ) * ff(mx,my,iz,iv,im)                                &
!!!                           + adbtc(mx,my,iz,iv,im)*phi(mx,my,iz)                    
!!!                end do
!!!              end do
!!!            end do
!!!          end do
!!!!$OMP end do nowait
!!!
!!!      end if
!!!
!!!!$OMP master
!!!                                    ! call fapp_stop("literm_colli_ct",1311,1)
!!!                                      call clock_end(1311)
!!!!$OMP end master
!!!
!!!
!!!  END SUBROUTINE colli_GK_CT
!!!
!!!
!!!!--------------------------------------
!!!  SUBROUTINE colli_GK_CT6( ff, phi, dfdvp, im, cf )
!!!!--------------------------------------
!!!!-------------------------------------------------------------------------------
!!!!
!!!!    Differential and FLR terms of test particle part in gyrokinetic collision
!!!!                                                         with 6th order CFD
!!!!    by M. Nakata, M. Nunami, April 2014
!!!!
!!!!-------------------------------------------------------------------------------
!!!
!!!    complex(kind=DP), intent(in), &
!!!      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
!!!    complex(kind=DP), intent(in), &
!!!      dimension(-nx:nx,0:ny,-nz:nz-1)                       :: phi
!!!    integer, intent(in) :: im
!!!    complex(kind=DP), intent(inout), &
!!!      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv)           :: cf
!!!
!!!    complex(kind=DP), intent(in), &
!!!      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb) :: dfdvp
!!!
!!!    real(kind=DP) :: cef1, cef2
!!!    real(kind=DP), dimension(-nz:nz-1) :: cef3
!!!    integer  ::  mx, my, iz, iv, is1
!!!
!!!
!!!!$OMP master
!!!                                      call clock_sta(1311)
!!!                                    ! call fapp_start("literm_colli_ct",1311,1)
!!!!$OMP end master
!!!
!!!
!!!      cef1   = 1._DP / ( 90._DP * dv * dv )
!!!      cef2   = 1._DP / ( 60._DP * dv )
!!!
!!!      do iz = -nz, nz-1
!!!        cef3(iz)   = 1._DP / ( 90._DP * dvp(iz) * dvp(iz) )
!!!      end do
!!!
!!!      if ( rankm == 0 .AND. im == 0 ) then
!!!
!!!        is1 = ranks
!!!!$OMP do schedule (dynamic)
!!!          do iv = 1, 2*nv
!!!            do iz = -nz, nz-1
!!!              do my = ist_y, iend_y
!!!                do mx = -nx, nx
!!!                   cf(mx,my,iz,iv) =                                               &
!!!                             ( +           ff(mx,my,iz,iv-3,im)                    &
!!!                               - 13.5_DP * ff(mx,my,iz,iv-2,im)                    &
!!!                               + 135._DP * ff(mx,my,iz,iv-1,im)                    &
!!!                               - 245._DP * ff(mx,my,iz,iv  ,im)                    &
!!!                               + 135._DP * ff(mx,my,iz,iv+1,im)                    &
!!!                               - 13.5_DP * ff(mx,my,iz,iv+2,im)                    &
!!!                               +           ff(mx,my,iz,iv+3,im)                    &
!!!                             ) * cef1                                              &
!!!                               * (   nu_ps(iz,iv,im)*vl(iv)**2                     & 
!!!                                 ) * 0.5_DP                                        &          
!!!                           + ( - 245._DP * ff(mx,my,iz,iv,im  )                    &
!!!                               + 270._DP * ff(mx,my,iz,iv,im+1)                    &
!!!                               -   27_DP * ff(mx,my,iz,iv,im+2)                    &
!!!                               +   2._DP * ff(mx,my,iz,iv,im+3)                    &
!!!                             ) * cef3(iz)                                          &
!!!                               * (   nu_ds(iz,iv,im)*vl(iv)**2                     & 
!!!                                 )                                                 &           
!!!                           + ( -           ff(mx,my,iz,iv-3,im)                    &
!!!                               +   9._DP * ff(mx,my,iz,iv-2,im)                    &
!!!                               -  45._DP * ff(mx,my,iz,iv-1,im)                    &
!!!                               +  45._DP * ff(mx,my,iz,iv+1,im)                    &
!!!                               -   9._DP * ff(mx,my,iz,iv+2,im)                    &
!!!                               +           ff(mx,my,iz,iv+3,im)                    &
!!!                             ) * cef2 * nu_gs(iz,iv,im)*vl(iv)                     & 
!!!                           + ( nu_hs(iz,iv,im)*xxa(iz,iv,im)**2*2._DP              & 
!!!                                - 0.25_DP*ksq(mx,my,iz)*Anum(is1)*tau(is1)         &
!!!                                 * ( nu_ds(iz,iv,im) * 2._DP*vl(iv)**2 )           &
!!!                                 / (Znum(is1)*omg(iz))**2 * real(iFLR, kind=DP)    &
!!!                             ) * ff(mx,my,iz,iv,im)                                &
!!!                           + adbtc(mx,my,iz,iv,im)*phi(mx,my,iz)
!!!                end do
!!!              end do
!!!            end do
!!!          end do
!!!!$OMP end do nowait
!!!
!!!      else if ( rankm == 0 .AND. im == 1 ) then
!!!
!!!        is1 = ranks
!!!!$OMP do schedule (dynamic)
!!!          do iv = 1, 2*nv
!!!            do iz = -nz, nz-1
!!!              do my = ist_y, iend_y
!!!                do mx = -nx, nx
!!!                   cf(mx,my,iz,iv) =                                               &
!!!                             ( +           ff(mx,my,iz,iv-3,im)                    &
!!!                               - 13.5_DP * ff(mx,my,iz,iv-2,im)                    &
!!!                               + 135._DP * ff(mx,my,iz,iv-1,im)                    &
!!!                               - 245._DP * ff(mx,my,iz,iv  ,im)                    &
!!!                               + 135._DP * ff(mx,my,iz,iv+1,im)                    &
!!!                               - 13.5_DP * ff(mx,my,iz,iv+2,im)                    &
!!!                               +           ff(mx,my,iz,iv+3,im)                    &
!!!                             ) * cef1                                              &
!!!                               * (   nu_ps(iz,iv,im)*vl(iv)**2                     & 
!!!                                   + nu_ds(iz,iv,im)*vp(iz,im)**2                  &
!!!                                 ) * 0.5_DP                                        &          
!!!                           + ( + 135._DP  * ff(mx,my,iz,iv,im-1)                   &
!!!                               - 258.5_DP * ff(mx,my,iz,iv,im  )                   &
!!!                               + 136._DP  * ff(mx,my,iz,iv,im+1)                   &
!!!                               - 13.5_DP  * ff(mx,my,iz,iv,im+2)                   &
!!!                               +            ff(mx,my,iz,iv,im+3)                   &
!!!                             ) * cef3(iz)                                          &
!!!                               * (   nu_ds(iz,iv,im)*vl(iv)**2                     & 
!!!                                   + nu_ps(iz,iv,im)*vp(iz,im)**2                  &
!!!                                 ) * 0.5_DP                                        &           
!!!                           + ( -           dfdvp(mx,my,iz,iv-3)                    &             
!!!                               +   9._DP * dfdvp(mx,my,iz,iv-2)                    &
!!!                               -  45._DP * dfdvp(mx,my,iz,iv-1)                    &
!!!                               +  45._DP * dfdvp(mx,my,iz,iv+1)                    &
!!!                               -   9._DP * dfdvp(mx,my,iz,iv+2)                    &
!!!                               +           dfdvp(mx,my,iz,iv+3)                    &
!!!                             ) * cef2 * vl(iv) * vp(iz,im)                         & 
!!!                               * (   nu_ps(iz,iv,im)                               & 
!!!                                   - nu_ds(iz,iv,im) )                             &
!!!                           + ( -           ff(mx,my,iz,iv-3,im)                    &
!!!                               +   9._DP * ff(mx,my,iz,iv-2,im)                    &
!!!                               -  45._DP * ff(mx,my,iz,iv-1,im)                    &
!!!                               +  45._DP * ff(mx,my,iz,iv+1,im)                    &
!!!                               -   9._DP * ff(mx,my,iz,iv+2,im)                    &
!!!                               +           ff(mx,my,iz,iv+3,im)                    &
!!!                             ) * cef2 * nu_gs(iz,iv,im)*vl(iv)                     & 
!!!                           + ( dfdvp(mx,my,iz,iv) )                                &
!!!                               * (   nu_gs(iz,iv,im)*vp(iz,im)                     & 
!!!                                   + nu_ds(iz,iv,im)*0.5_DP                        &
!!!                                     * (vl(iv)**2/vp(iz,im)+vp(iz,im))             &
!!!                                 )                                                 &
!!!                           + ( nu_hs(iz,iv,im)*xxa(iz,iv,im)**2*2._DP              & 
!!!                                - 0.25_DP*ksq(mx,my,iz)*Anum(is1)*tau(is1)         &
!!!                                 * ( nu_ds(iz,iv,im)                               &
!!!                                      * (2._DP*vl(iv)**2 + vp(iz,im)**2)           &
!!!                                    + nu_ps(iz,iv,im)*vp(iz,im)**2 )               &
!!!                                 / (Znum(is1)*omg(iz))**2 * real(iFLR, kind=DP)    &
!!!                             ) * ff(mx,my,iz,iv,im)                                &
!!!                           + adbtc(mx,my,iz,iv,im)*phi(mx,my,iz)                    
!!!                end do
!!!              end do
!!!            end do
!!!          end do
!!!!$OMP end do nowait
!!!
!!!      else if ( rankm == 0 .AND. im == 2 ) then
!!!
!!!        is1 = ranks
!!!!$OMP do schedule (dynamic)
!!!          do iv = 1, 2*nv
!!!            do iz = -nz, nz-1
!!!              do my = ist_y, iend_y
!!!                do mx = -nx, nx
!!!                   cf(mx,my,iz,iv) =                                               &
!!!                             ( +           ff(mx,my,iz,iv-3,im)                    &
!!!                               - 13.5_DP * ff(mx,my,iz,iv-2,im)                    &
!!!                               + 135._DP * ff(mx,my,iz,iv-1,im)                    &
!!!                               - 245._DP * ff(mx,my,iz,iv  ,im)                    &
!!!                               + 135._DP * ff(mx,my,iz,iv+1,im)                    &
!!!                               - 13.5_DP * ff(mx,my,iz,iv+2,im)                    &
!!!                               +           ff(mx,my,iz,iv+3,im)                    &
!!!                             ) * cef1                                              &
!!!                               * (   nu_ps(iz,iv,im)*vl(iv)**2                     & 
!!!                                   + nu_ds(iz,iv,im)*vp(iz,im)**2                  &
!!!                                 ) * 0.5_DP                                        &          
!!!                           + ( - 13.5_DP * ff(mx,my,iz,iv,im-2)                    &
!!!                               + 136._DP * ff(mx,my,iz,iv,im-1)                    &
!!!                               - 245._DP * ff(mx,my,iz,iv,im  )                    &
!!!                               + 135._DP * ff(mx,my,iz,iv,im+1)                    &
!!!                               - 13.5_DP * ff(mx,my,iz,iv,im+2)                    &
!!!                               +           ff(mx,my,iz,iv,im+3)                    &
!!!                             ) * cef3(iz)                                          &
!!!                               * (   nu_ds(iz,iv,im)*vl(iv)**2                     & 
!!!                                   + nu_ps(iz,iv,im)*vp(iz,im)**2                  &
!!!                                 ) * 0.5_DP                                        &           
!!!                           + ( -           dfdvp(mx,my,iz,iv-3)                    &             
!!!                               +   9._DP * dfdvp(mx,my,iz,iv-2)                    &
!!!                               -  45._DP * dfdvp(mx,my,iz,iv-1)                    &
!!!                               +  45._DP * dfdvp(mx,my,iz,iv+1)                    &
!!!                               -   9._DP * dfdvp(mx,my,iz,iv+2)                    &
!!!                               +           dfdvp(mx,my,iz,iv+3)                    &
!!!                             ) * cef2 * vl(iv) * vp(iz,im)                         & 
!!!                               * (   nu_ps(iz,iv,im)                               & 
!!!                                   - nu_ds(iz,iv,im) )                             &
!!!                           + ( -           ff(mx,my,iz,iv-3,im)                    &
!!!                               +   9._DP * ff(mx,my,iz,iv-2,im)                    &
!!!                               -  45._DP * ff(mx,my,iz,iv-1,im)                    &
!!!                               +  45._DP * ff(mx,my,iz,iv+1,im)                    &
!!!                               -   9._DP * ff(mx,my,iz,iv+2,im)                    &
!!!                               +           ff(mx,my,iz,iv+3,im)                    &
!!!                             ) * cef2 * nu_gs(iz,iv,im)*vl(iv)                     & 
!!!                           + ( dfdvp(mx,my,iz,iv) )                                &
!!!                               * (   nu_gs(iz,iv,im)*vp(iz,im)                     & 
!!!                                   + nu_ds(iz,iv,im)*0.5_DP                        &
!!!                                     * (vl(iv)**2/vp(iz,im)+vp(iz,im))             &
!!!                                 )                                                 &
!!!                           + ( nu_hs(iz,iv,im)*xxa(iz,iv,im)**2*2._DP              & 
!!!                                - 0.25_DP*ksq(mx,my,iz)*Anum(is1)*tau(is1)         &
!!!                                 * ( nu_ds(iz,iv,im)                               &
!!!                                      * (2._DP*vl(iv)**2 + vp(iz,im)**2)           &
!!!                                    + nu_ps(iz,iv,im)*vp(iz,im)**2 )               &
!!!                                 / (Znum(is1)*omg(iz))**2 * real(iFLR, kind=DP)    &
!!!                             ) * ff(mx,my,iz,iv,im)                                &
!!!                           + adbtc(mx,my,iz,iv,im)*phi(mx,my,iz)                    
!!!                end do
!!!              end do
!!!            end do
!!!          end do
!!!!$OMP end do nowait
!!!
!!!      else  ! im=[0,nm] for rankm > 0 and im=[3,nm] for rankm = 0  
!!!
!!!        is1 = ranks
!!!!$OMP do schedule (dynamic)
!!!          do iv = 1, 2*nv
!!!            do iz = -nz, nz-1
!!!              do my = ist_y, iend_y
!!!                do mx = -nx, nx
!!!                   cf(mx,my,iz,iv) =                                               &
!!!                             ( +           ff(mx,my,iz,iv-3,im)                    &
!!!                               - 13.5_DP * ff(mx,my,iz,iv-2,im)                    &
!!!                               + 135._DP * ff(mx,my,iz,iv-1,im)                    &
!!!                               - 245._DP * ff(mx,my,iz,iv  ,im)                    &
!!!                               + 135._DP * ff(mx,my,iz,iv+1,im)                    &
!!!                               - 13.5_DP * ff(mx,my,iz,iv+2,im)                    &
!!!                               +           ff(mx,my,iz,iv+3,im)                    &
!!!                             ) * cef1                                              &
!!!                               * (   nu_ps(iz,iv,im)*vl(iv)**2                     & 
!!!                                   + nu_ds(iz,iv,im)*vp(iz,im)**2                  &
!!!                                 ) * 0.5_DP                                        &          
!!!                           + ( +           ff(mx,my,iz,iv,im-3)                    &
!!!                               - 13.5_DP * ff(mx,my,iz,iv,im-2)                    &
!!!                               + 135._DP * ff(mx,my,iz,iv,im-1)                    &
!!!                               - 245._DP * ff(mx,my,iz,iv,im  )                    &
!!!                               + 135._DP * ff(mx,my,iz,iv,im+1)                    &
!!!                               - 13.5_DP * ff(mx,my,iz,iv,im+2)                    &
!!!                               +           ff(mx,my,iz,iv,im+3)                    &
!!!                             ) * cef3(iz)                                          &
!!!                               * (   nu_ds(iz,iv,im)*vl(iv)**2                     & 
!!!                                   + nu_ps(iz,iv,im)*vp(iz,im)**2                  &
!!!                                 ) * 0.5_DP                                        &           
!!!                           + ( -           dfdvp(mx,my,iz,iv-3)                    &             
!!!                               +   9._DP * dfdvp(mx,my,iz,iv-2)                    &
!!!                               -  45._DP * dfdvp(mx,my,iz,iv-1)                    &
!!!                               +  45._DP * dfdvp(mx,my,iz,iv+1)                    &
!!!                               -   9._DP * dfdvp(mx,my,iz,iv+2)                    &
!!!                               +           dfdvp(mx,my,iz,iv+3)                    &
!!!                             ) * cef2 * vl(iv) * vp(iz,im)                         & 
!!!                               * (   nu_ps(iz,iv,im)                               & 
!!!                                   - nu_ds(iz,iv,im) )                             &
!!!                           + ( -           ff(mx,my,iz,iv-3,im)                    &
!!!                               +   9._DP * ff(mx,my,iz,iv-2,im)                    &
!!!                               -  45._DP * ff(mx,my,iz,iv-1,im)                    &
!!!                               +  45._DP * ff(mx,my,iz,iv+1,im)                    &
!!!                               -   9._DP * ff(mx,my,iz,iv+2,im)                    &
!!!                               +           ff(mx,my,iz,iv+3,im)                    &
!!!                             ) * cef2 * nu_gs(iz,iv,im)*vl(iv)                     & 
!!!                           + ( dfdvp(mx,my,iz,iv) )                                &
!!!                               * (   nu_gs(iz,iv,im)*vp(iz,im)                     & 
!!!                                   + nu_ds(iz,iv,im)*0.5_DP                        &
!!!                                     * (vl(iv)**2/vp(iz,im)+vp(iz,im))             &
!!!                                 )                                                 &
!!!                           + ( nu_hs(iz,iv,im)*xxa(iz,iv,im)**2*2._DP              & 
!!!                                - 0.25_DP*ksq(mx,my,iz)*Anum(is1)*tau(is1)         &
!!!                                 * ( nu_ds(iz,iv,im)                               &
!!!                                      * (2._DP*vl(iv)**2 + vp(iz,im)**2)           &
!!!                                    + nu_ps(iz,iv,im)*vp(iz,im)**2 )               &
!!!                                 / (Znum(is1)*omg(iz))**2 * real(iFLR, kind=DP)    &
!!!                             ) * ff(mx,my,iz,iv,im)                                &
!!!                           + adbtc(mx,my,iz,iv,im)*phi(mx,my,iz)                    
!!!                end do
!!!              end do
!!!            end do
!!!          end do
!!!!$OMP end do nowait
!!!
!!!      end if
!!!
!!!!$OMP master
!!!                                    ! call fapp_stop("literm_colli_ct",1311,1)
!!!                                      call clock_end(1311)
!!!!$OMP end master
!!!
!!!
!!!  END SUBROUTINE colli_GK_CT6
!!!
!!!
!!!!--------------------------------------
!!!  SUBROUTINE colli_GK_DT( moment_ab_wk, im, cf )
!!!!--------------------------------------
!!!!-------------------------------------------------------------------------------
!!!!
!!!!    Non-isothermal terms of test particle part and field particle part 
!!!!                                                      in gyrokinetic collision
!!!!
!!!!    by M. Nakata, M. Nunami, April 2014
!!!!
!!!!-------------------------------------------------------------------------------
!!!
!!!    complex(kind=DP), intent(in), &
!!!      dimension(1:6,-nx:nx,0:ny,-nz:nz-1,0:ns-1)       :: moment_ab_wk
!!!    integer, intent(in)                                :: im
!!!    complex(kind=DP), intent(inout), &
!!!      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv)           :: cf
!!!
!!!    integer  ::  mx, my, iz, iv, is1, is2
!!!
!!!
!!!!$OMP master
!!!                                      call clock_sta(1312)
!!!                                    ! call fapp_start("literm_colli_dt",1312,1)
!!!!$OMP end master
!!!
!!!
!!!     if ( iFLR == 1 ) then  ! full-GK
!!!
!!!       is1 = ranks
!!!         do is2 = 0, ns-1
!!!!%%%% do schedule (dynamic)
!!!!$OMP do
!!!           do iv = 1, 2*nv
!!!             do iz = -nz, nz-1
!!!               do my = ist_y, iend_y
!!!                 do mx = -nx, nx
!!!                   cf(mx,my,iz,iv) = cf(mx,my,iz,iv)                     &
!!!                                    + x_tst(1,iz,iv,im,is1,is2)          &
!!!                                     * moment_ab_wk(1,mx,my,iz,is2)      &
!!!                                      * j0(mx,my,iz,im)                  & 
!!!                                    + x_tst(2,iz,iv,im,is1,is2)          &
!!!                                     * moment_ab_wk(2,mx,my,iz,is2)      &
!!!                                      * j1(mx,my,iz,im)                  & 
!!!                                    + x_tst(3,iz,iv,im,is1,is2)          &
!!!                                     * moment_ab_wk(3,mx,my,iz,is2)      &
!!!                                      * j0(mx,my,iz,im)                  & 
!!!                                    + x_tst(4,iz,iv,im,is1,is2)          &
!!!                                     * moment_ab_wk(4,mx,my,iz,is2)      &
!!!                                      * j0(mx,my,iz,im)                  & 
!!!                                    + x_tst(5,iz,iv,im,is1,is2)          &
!!!                                     * moment_ab_wk(5,mx,my,iz,is2)      &
!!!                                      * j1(mx,my,iz,im)                  & 
!!!                                    + x_tst(6,iz,iv,im,is1,is2)          &
!!!                                     * moment_ab_wk(6,mx,my,iz,is2)      &
!!!                                      * j0(mx,my,iz,im)                  
!!!                end do
!!!              end do
!!!            end do
!!!          end do
!!!!%%%% end do nowait
!!!!$OMP end do
!!!        end do
!!!
!!!     else if ( iFLR == 0 ) then ! DK-limit
!!!
!!!       is1 = ranks
!!!         do is2 = 0, ns-1
!!!!%%%% do schedule (dynamic)
!!!!$OMP do
!!!           do iv = 1, 2*nv
!!!             do iz = -nz, nz-1
!!!               do my = ist_y, iend_y
!!!                 do mx = -nx, nx
!!!                   cf(mx,my,iz,iv) = cf(mx,my,iz,iv)                     &
!!!                                    + x_tst(1,iz,iv,im,is1,is2)          &
!!!                                     * moment_ab_wk(1,mx,my,iz,is2)      &
!!!                                    + x_tst(3,iz,iv,im,is1,is2)          &
!!!                                     * moment_ab_wk(3,mx,my,iz,is2)      &
!!!                                    + x_tst(4,iz,iv,im,is1,is2)          &
!!!                                     * moment_ab_wk(4,mx,my,iz,is2)      &
!!!                                    + x_tst(6,iz,iv,im,is1,is2)          &
!!!                                     * moment_ab_wk(6,mx,my,iz,is2)          
!!!                end do
!!!              end do
!!!            end do
!!!          end do
!!!!%%%% end do nowait
!!!!$OMP end do
!!!        end do
!!!
!!!     end if
!!!
!!!
!!!!$OMP master
!!!                                    ! call fapp_stop("literm_colli_dt",1312,1)
!!!                                      call clock_end(1312)
!!!!$OMP end master
!!!
!!!
!!!  END SUBROUTINE colli_GK_DT
!!!
!!!
!!!!--------------------------------------
!!!  SUBROUTINE colli_GK_CF_DT(moment_ba_wk, moment_ab_wk, im, cf)
!!!!--------------------------------------
!!!!-------------------------------------------------------------------------------
!!!!
!!!!    Field particle and non-isothermal parts in gyrokinetic collision
!!!!
!!!!    by M. Nakata, M. Nunami, April 2014
!!!!
!!!!-------------------------------------------------------------------------------
!!!
!!!    complex(kind=DP), intent(in), &
!!!      dimension(1:6,-nx:nx,0:ny,-nz:nz-1,0:ns-1)       :: moment_ba_wk, moment_ab_wk
!!!    integer, intent(in)                                :: im
!!!    complex(kind=DP), intent(inout), &
!!!      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv)           :: cf
!!!
!!!    integer  ::  mx, my, iz, iv, is1, is2
!!!
!!!
!!!!$OMP master
!!!                                      call clock_sta(1313)
!!!                                    ! call fapp_start("literm_colli_cf",1313,1)
!!!!$OMP end master
!!!
!!!     if ( iFLR == 1 ) then  ! full-GK
!!!
!!!       is1 = ranks
!!!         do is2 = 0, ns-1
!!!!%%%% do schedule (dynamic)
!!!!$OMP do
!!!           do iv = 1, 2*nv
!!!             do iz = -nz, nz-1
!!!               do my = ist_y, iend_y
!!!                 do mx = -nx, nx
!!!                   cf(mx,my,iz,iv) = cf(mx,my,iz,iv)                     &
!!!                                    + y_fld(1,iz,iv,im,is1,is2)          &
!!!                                     * moment_ba_wk(1,mx,my,iz,is2)      &
!!!                                      * j0(mx,my,iz,im)                  & 
!!!                                    + y_fld(2,iz,iv,im,is1,is2)          &
!!!                                     * moment_ba_wk(2,mx,my,iz,is2)      &
!!!                                      * j1(mx,my,iz,im)                  & 
!!!                                    + y_fld(3,iz,iv,im,is1,is2)          &
!!!                                     * moment_ba_wk(3,mx,my,iz,is2)      &
!!!                                      * j0(mx,my,iz,im)                  & 
!!!                                    + y_fld(4,iz,iv,im,is1,is2)          &
!!!                                     * moment_ba_wk(4,mx,my,iz,is2)      &
!!!                                      * j0(mx,my,iz,im)                  & 
!!!                                    + y_fld(5,iz,iv,im,is1,is2)          &
!!!                                     * moment_ba_wk(5,mx,my,iz,is2)      &
!!!                                      * j1(mx,my,iz,im)                  & 
!!!                                    + y_fld(6,iz,iv,im,is1,is2)          &
!!!                                     * moment_ba_wk(6,mx,my,iz,is2)      &
!!!                                      * j0(mx,my,iz,im)                  &
!!!
!!!                                    + x_tst(1,iz,iv,im,is1,is2)          &
!!!                                     * moment_ab_wk(1,mx,my,iz,is2)      &
!!!                                      * j0(mx,my,iz,im)                  & 
!!!                                    + x_tst(2,iz,iv,im,is1,is2)          &
!!!                                     * moment_ab_wk(2,mx,my,iz,is2)      &
!!!                                      * j1(mx,my,iz,im)                  & 
!!!                                    + x_tst(3,iz,iv,im,is1,is2)          &
!!!                                     * moment_ab_wk(3,mx,my,iz,is2)      &
!!!                                      * j0(mx,my,iz,im)                  & 
!!!                                    + x_tst(4,iz,iv,im,is1,is2)          &
!!!                                     * moment_ab_wk(4,mx,my,iz,is2)      &
!!!                                      * j0(mx,my,iz,im)                  & 
!!!                                    + x_tst(5,iz,iv,im,is1,is2)          &
!!!                                     * moment_ab_wk(5,mx,my,iz,is2)      &
!!!                                      * j1(mx,my,iz,im)                  & 
!!!                                    + x_tst(6,iz,iv,im,is1,is2)          &
!!!                                     * moment_ab_wk(6,mx,my,iz,is2)      &
!!!                                      * j0(mx,my,iz,im)                  
!!!                end do
!!!              end do
!!!            end do
!!!          end do
!!!!%%%% end do nowait
!!!!$OMP end do
!!!        end do
!!!
!!!     else if ( iFLR == 0 ) then ! DK-limit
!!!
!!!       is1 = ranks
!!!         do is2 = 0, ns-1
!!!!%%%% do schedule (dynamic)
!!!!$OMP do
!!!           do iv = 1, 2*nv
!!!             do iz = -nz, nz-1
!!!               do my = ist_y, iend_y
!!!                 do mx = -nx, nx
!!!                   cf(mx,my,iz,iv) = cf(mx,my,iz,iv)                     &
!!!                                    + y_fld(1,iz,iv,im,is1,is2)          &
!!!                                     * moment_ba_wk(1,mx,my,iz,is2)      &
!!!                                    + y_fld(3,iz,iv,im,is1,is2)          &
!!!                                     * moment_ba_wk(3,mx,my,iz,is2)      &
!!!                                    + y_fld(4,iz,iv,im,is1,is2)          &
!!!                                     * moment_ba_wk(4,mx,my,iz,is2)      &
!!!                                    + y_fld(6,iz,iv,im,is1,is2)          &
!!!                                     * moment_ba_wk(6,mx,my,iz,is2)      &
!!!
!!!                                    + x_tst(1,iz,iv,im,is1,is2)          &
!!!                                     * moment_ab_wk(1,mx,my,iz,is2)      &
!!!                                    + x_tst(3,iz,iv,im,is1,is2)          &
!!!                                     * moment_ab_wk(3,mx,my,iz,is2)      &
!!!                                    + x_tst(4,iz,iv,im,is1,is2)          &
!!!                                     * moment_ab_wk(4,mx,my,iz,is2)      &
!!!                                    + x_tst(6,iz,iv,im,is1,is2)          &
!!!                                     * moment_ab_wk(6,mx,my,iz,is2)          
!!!                end do
!!!              end do
!!!            end do
!!!          end do
!!!!%%%% end do nowait
!!!!$OMP end do
!!!        end do
!!!
!!!     end if
!!!
!!!!$OMP master
!!!                                    ! call fapp_stop("literm_colli_cf",1313,1)
!!!                                      call clock_end(1313)
!!!!$OMP end master
!!!
!!!
!!!  END SUBROUTINE colli_GK_CF_DT
!!!
!!!
!!!!--------------------------------------
!!!  SUBROUTINE colli_GK_CF(moment_ba_wk, im, cf)
!!!!--------------------------------------
!!!!-------------------------------------------------------------------------------
!!!!
!!!!    Field particle part in gyrokinetic collision
!!!!
!!!!    by M. Nakata, M. Nunami, April 2014
!!!!
!!!!-------------------------------------------------------------------------------
!!!
!!!    complex(kind=DP), intent(in), &
!!!      dimension(1:6,-nx:nx,0:ny,-nz:nz-1,0:ns-1)       :: moment_ba_wk
!!!    integer, intent(in)                                :: im
!!!    complex(kind=DP), intent(inout), &
!!!      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv)           :: cf
!!!
!!!    integer  ::  mx, my, iz, iv, is1, is2
!!!
!!!
!!!!$OMP master
!!!                                      call clock_sta(1313)
!!!                                    ! call fapp_start("literm_colli_cf",1313,1)
!!!!$OMP end master
!!!
!!!     if ( iFLR == 1 ) then  ! full-GK
!!!
!!!       is1 = ranks
!!!         do is2 = 0, ns-1
!!!!%%%% do schedule (dynamic)
!!!!$OMP do
!!!           do iv = 1, 2*nv
!!!             do iz = -nz, nz-1
!!!               do my = ist_y, iend_y
!!!                 do mx = -nx, nx
!!!                   cf(mx,my,iz,iv) = cf(mx,my,iz,iv)                     &
!!!                                    + y_fld(1,iz,iv,im,is1,is2)          &
!!!                                     * moment_ba_wk(1,mx,my,iz,is2)      &
!!!                                      * j0(mx,my,iz,im)                  & 
!!!                                    + y_fld(2,iz,iv,im,is1,is2)          &
!!!                                     * moment_ba_wk(2,mx,my,iz,is2)      &
!!!                                      * j1(mx,my,iz,im)                  & 
!!!                                    + y_fld(3,iz,iv,im,is1,is2)          &
!!!                                     * moment_ba_wk(3,mx,my,iz,is2)      &
!!!                                      * j0(mx,my,iz,im)                  & 
!!!                                    + y_fld(4,iz,iv,im,is1,is2)          &
!!!                                     * moment_ba_wk(4,mx,my,iz,is2)      &
!!!                                      * j0(mx,my,iz,im)                  & 
!!!                                    + y_fld(5,iz,iv,im,is1,is2)          &
!!!                                     * moment_ba_wk(5,mx,my,iz,is2)      &
!!!                                      * j1(mx,my,iz,im)                  & 
!!!                                    + y_fld(6,iz,iv,im,is1,is2)          &
!!!                                     * moment_ba_wk(6,mx,my,iz,is2)      &
!!!                                      * j0(mx,my,iz,im)                  
!!!                end do
!!!              end do
!!!            end do
!!!          end do
!!!!%%%% end do nowait
!!!!$OMP end do
!!!        end do
!!!
!!!     else if ( iFLR == 0 ) then ! DK-limit
!!!
!!!       is1 = ranks
!!!         do is2 = 0, ns-1
!!!!%%%% do schedule (dynamic)
!!!!$OMP do
!!!           do iv = 1, 2*nv
!!!             do iz = -nz, nz-1
!!!               do my = ist_y, iend_y
!!!                 do mx = -nx, nx
!!!                   cf(mx,my,iz,iv) = cf(mx,my,iz,iv)                     &
!!!                                    + y_fld(1,iz,iv,im,is1,is2)          &
!!!                                     * moment_ba_wk(1,mx,my,iz,is2)      &
!!!                                    + y_fld(3,iz,iv,im,is1,is2)          &
!!!                                     * moment_ba_wk(3,mx,my,iz,is2)      &
!!!                                    + y_fld(4,iz,iv,im,is1,is2)          &
!!!                                     * moment_ba_wk(4,mx,my,iz,is2)      &
!!!                                    + y_fld(6,iz,iv,im,is1,is2)          &
!!!                                     * moment_ba_wk(6,mx,my,iz,is2)         
!!!                end do
!!!              end do
!!!            end do
!!!          end do
!!!!%%%% end do nowait
!!!!$OMP end do
!!!        end do
!!!
!!!     end if
!!!
!!!!$OMP master
!!!                                    ! call fapp_stop("literm_colli_cf",1313,1)
!!!                                      call clock_end(1313)
!!!!$OMP end master
!!!
!!!
!!!  END SUBROUTINE colli_GK_CF
!!!
!!!
!!!!--------------------------------------
!!!  SUBROUTINE colli_moment_calc( hh, phi, ww )
!!!!--------------------------------------
!!!!-------------------------------------------------------------------------------
!!!!
!!!!    Moment calculations for gyrokinetic collision: local velocity moment part 
!!!!
!!!!    by M. Nakata, M. Nunami, April 2014
!!!!
!!!!-------------------------------------------------------------------------------
!!!
!!!    complex(kind=DP), intent(in), &
!!!      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: hh
!!!    complex(kind=DP), intent(in), &
!!!      dimension(-nx:nx,0:ny,-nz:nz-1)               :: phi
!!!
!!!    !complex(kind=DP), intent(out), &
!!!    complex(kind=DP), intent(inout), &
!!!      dimension(-nx:nx,0:ny,-nz:nz-1,0:ns-1,1:6)        :: ww 
!!!
!!!    real(kind=DP) :: v2a, v2b, dflg
!!!    complex(kind=DP) :: wf1, wf2
!!!    integer :: mx, my, iz, iv, im, is1, is2, ii
!!!
!!!!$OMP master
!!!                                      call clock_sta(1314)
!!!                                    ! call fapp_start("literm_colli_mom",1314,1)
!!!!$OMP end master
!!!
!!!      if ( rankm == 0 ) then
!!!
!!!!$OMP do collapse(2) schedule(dynamic) private(mx,my,iz,iv,im,is2,ii,v2a,v2b,wf1,wf2)
!!!      do ii = 1, 6
!!!        do is2 = 0, ns-1
!!!          do im = 0, nm
!!!            do iv = 1, 2*nv
!!!              do iz = -nz, nz-1
!!!                do my = ist_y, iend_y
!!!                  do mx = -nx, nx
!!!                    !ww(mx,my,iz,is2,ii) = ww(mx,my,iz,is2,ii) + hh(mx,my,iz,iv,im)*vfunc(mx,my,iz,iv,im,is2,ii)
!!!                    ww(mx,my,iz,is2,ii) = ww(mx,my,iz,is2,ii) + hh(mx,my,iz,iv,im)*jfunc(mx,my,iz,im,ii)*vfunc(iz,iv,im,is2,ii)
!!!                  end do
!!!                end do
!!!              end do
!!!            end do
!!!          end do
!!!
!!!! for edge compensation
!!!!         im  = 1
!!!          do im = 1, 1
!!!            do iv = 1, 2*nv
!!!              do iz = -nz, nz-1
!!!                do my = ist_y, iend_y
!!!                  do mx = -nx, nx
!!!                    !v2a = vfunc(mx,my,iz,iv,im,is2,ii)
!!!                    !v2b = vfunc(mx,my,iz,iv,im+1,is2,ii)
!!!                    v2a = jfunc(mx,my,iz,im,ii)  *vfunc(iz,iv,im,is2,ii)
!!!                    v2b = jfunc(mx,my,iz,im+1,ii)*vfunc(iz,iv,im+1,is2,ii)
!!!                    wf1 = hh(mx,my,iz,iv,im)   
!!!                    wf2 = hh(mx,my,iz,iv,im+1)
!!!                    ww(mx,my,iz,is2,ii)  = ww(mx,my,iz,is2,ii)            &
!!!                          - ( - wf1/12._DP*v2a + ( wf2*v2b - wf1*2._DP*v2a )*11._DP/720._DP ) 
!!!                  end do
!!!                end do
!!!              end do
!!!            end do
!!!          end do
!!!        end do
!!!      end do
!!!!$OMP end do nowait
!!!
!!!      else
!!!
!!!!$OMP do collapse(2) schedule(dynamic) private(mx,my,iz,iv,im,is2,ii)
!!!      do ii=1,6
!!!        do is2 = 0, ns-1
!!!          do im = 0, nm
!!!            do iv = 1, 2*nv
!!!              do iz = -nz, nz-1
!!!                do my = ist_y, iend_y
!!!                  do mx = -nx, nx
!!!                    !ww(mx,my,iz,is2,ii) = ww(mx,my,iz,is2,ii) + hh(mx,my,iz,iv,im)*vfunc(mx,my,iz,iv,im,is2,ii)
!!!                    ww(mx,my,iz,is2,ii) = ww(mx,my,iz,is2,ii) + hh(mx,my,iz,iv,im)*jfunc(mx,my,iz,im,ii)*vfunc(iz,iv,im,is2,ii)
!!!                  end do
!!!                end do
!!!              end do
!!!            end do
!!!          end do
!!!        end do
!!!      end do
!!!!$OMP end do nowait
!!!
!!!      end if
!!!
!!!
!!!!$OMP master
!!!                                    ! call fapp_stop("literm_colli_mom",1314,1)
!!!                                      call clock_end(1314)
!!!!$OMP end master
!!!
!!!  END SUBROUTINE colli_moment_calc
!!!
!!!
!!!!--------------------------------------
!!!  SUBROUTINE colli_moment_redc( ww, wn )
!!!!--------------------------------------
!!!!-------------------------------------------------------------------------------
!!!!
!!!!    Moment calculations for gyrokinetic collision: All_reduce_part
!!!!
!!!!    by M. Nakata, M. Nunami, April 2014
!!!!
!!!!-------------------------------------------------------------------------------
!!!
!!!    complex(kind=DP), intent(in), &
!!!      dimension(-nx:nx,0:ny,-nz:nz-1,0:ns-1,1:6)        :: ww
!!! 
!!!    complex(kind=DP), intent(out), &
!!!      dimension(-nx:nx,0:ny,-nz:nz-1,0:ns-1,1:6)        :: wn
!!!
!!!
!!!!$OMP master
!!!                                           call clock_sta(1315)
!!!                                         ! call fapp_start("literm_colli_ar",1315,1)
!!!!$OMP end master
!!!
!!!      call MPI_Allreduce( ww, wn, nxyz*ns*6, MPI_DOUBLE_COMPLEX, &
!!!                          MPI_SUM, vel_comm_world, ierr_mpi )
!!!
!!!!$OMP master
!!!                                         ! call fapp_stop("literm_colli_ar",1315,1)
!!!                                           call clock_end(1315)
!!!!$OMP end master
!!!
!!!  END SUBROUTINE colli_moment_redc
!!!
!!!
!!!!--------------------------------------
!!!  SUBROUTINE colli_comm_alltoall( wm, wn )
!!!!--------------------------------------
!!!!-------------------------------------------------------------------------------
!!!!
!!!!    Inter-species communication of moment quantities for field particle part
!!!!       with MPI_AlltoAll
!!!!
!!!!    by M. Nakata, M. Nunami, April 2014
!!!!
!!!!-------------------------------------------------------------------------------
!!!
!!!    complex(kind=DP), intent(in), &
!!!      dimension(-nx:nx,0:ny,-nz:nz-1,0:ns-1,1:6)  :: wm
!!!    complex(kind=DP), intent(out), &
!!!      dimension(-nx:nx,0:ny,-nz:nz-1,0:ns-1,1:6)  :: wn
!!!
!!!    complex(kind=DP),              & 
!!!      dimension(-nx:nx,0:ny,-nz:nz-1,1:6,0:ns-1)  :: send_buff, recv_buff
!!!
!!!    integer :: mx, my, iz, is, ii
!!!    integer :: datasize, datasize_ns
!!!
!!!
!!!!$OMP master
!!!                                       call clock_sta(1316)
!!!                                     ! call fapp_start("literm_colli_com",1316,1)
!!!!$OMP end master
!!!
!!!
!!!      datasize = (2*nx+1)*(ny+1)*(2*nz)*6
!!!      datasize_ns = (2*nx+1)*(ny+1)*(2*nz)*6*ns
!!!
!!!      if ( vel_rank == 0 ) then
!!!
!!!        do ii = 1, 6
!!!          do is = 0, ns-1
!!!            do iz = -nz, nz-1
!!!              do my = ist_y, iend_y
!!!                do mx = -nx, nx
!!!!                  send_buff(mx,my,iz,ii,is)  = real(ranks,kind=DP) + ii ! for debug
!!!                  send_buff(mx,my,iz,ii,is)  = wm(mx,my,iz,is,ii)
!!!                end do
!!!              end do 
!!!            end do
!!!          end do
!!!        end do
!!!
!!!
!!!          call MPI_Alltoall( send_buff(-nx,ist_y,-nz,1,0), datasize, MPI_DOUBLE_COMPLEX,  &
!!!                             recv_buff(-nx,ist_y,-nz,1,0), datasize, MPI_DOUBLE_COMPLEX,  & 
!!!                             col_comm_world, ierr_mpi  )
!!! 
!!!!! --- for debug 
!!!!        write(unit=8000+ranks,fmt="(I3,SP,6ES24.15e3)") ranks, real(recv_buff(0,1,1,1,0)), real(recv_buff(0,1,1,6,0)), &
!!!!                                                           real(recv_buff(0,1,1,1,1)), real(recv_buff(0,1,1,6,1)), &
!!!!                                                           real(recv_buff(0,1,1,1,2)), real(recv_buff(0,1,1,6,2))
!!!!        write(unit=8000+ranks,fmt=*)
!!!
!!!      end if
!!!
!!!
!!!      call MPI_Bcast( recv_buff(-nx,ist_y,-nz,1,0), datasize_ns, MPI_DOUBLE_COMPLEX, & 
!!!                      0, vel_comm_world, ierr_mpi  ) 
!!!
!!!
!!!      do is = 0, ns-1
!!!        do ii = 1, 6
!!!          do iz = -nz, nz-1
!!!            do my = ist_y, iend_y
!!!              do mx = -nx, nx
!!!                wn(mx,my,iz,is,ii) = recv_buff(mx,my,iz,ii,is)
!!!              end do
!!!            end do 
!!!          end do
!!!        end do         
!!!      end do         
!!!
!!!!$OMP master
!!!                                     ! call fapp_stop("literm_colli_com",1316,1)
!!!                                       call clock_end(1316)
!!!!$OMP end master
!!!
!!!  END SUBROUTINE colli_comm_alltoall
!!!
!!!
!!!!--------------------------------------
!!!  SUBROUTINE colli_dfdvp( ff, dfdvp )
!!!!--------------------------------------
!!!!-------------------------------------------------------------------------------
!!!!
!!!!    calculation of df/dv_perp term with 4th order CFD
!!!!
!!!!    by M. Nakata, M. Nunami, April 2014
!!!!
!!!!-------------------------------------------------------------------------------
!!!
!!!    complex(kind=DP), intent(in), &
!!!      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
!!!
!!!    complex(kind=DP), intent(out), &
!!!      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: dfdvp
!!!
!!!    real(kind=DP), dimension(-nz:nz-1) :: cef4
!!!    integer  ::  mx, my, iz, iv, im
!!!
!!!
!!!!$OMP master
!!!                                      call clock_sta(1317)
!!!                                    ! call fapp_start("literm_colli_dvp",1317,1)
!!!!$OMP end master
!!!
!!!
!!!      do iz = -nz, nz-1
!!!        cef4(iz)   = 1._DP / ( 12._DP * dvp(iz) )
!!!      end do
!!!
!!!
!!!      if ( rankm == 0 ) then 
!!!
!!!        im = 0
!!!!$OMP do schedule (dynamic)
!!!          do iv = 1, 2*nv
!!!            do iz = -nz, nz-1
!!!              do my = ist_y, iend_y
!!!                do mx = -nx, nx
!!!                  dfdvp(mx,my,iz,iv,im) = (0._DP, 0._DP) 
!!!                end do
!!!              end do
!!!            end do
!!!          end do
!!!!$OMP end do nowait
!!!
!!!        im = 1
!!!!$OMP do schedule (dynamic)
!!!          do iv = 1, 2*nv
!!!            do iz = -nz, nz-1
!!!              do my = ist_y, iend_y
!!!                do mx = -nx, nx
!!!                  dfdvp(mx,my,iz,iv,im) =                                     &
!!!                             ( -  8._DP * ff(mx,my,iz,iv,im-1)                &
!!!                               +          ff(mx,my,iz,iv,im  )                &
!!!                               +  8._DP * ff(mx,my,iz,iv,im+1)                &
!!!                               -          ff(mx,my,iz,iv,im+2)                &
!!!                             ) * cef4(iz)
!!!                end do
!!!              end do
!!!            end do
!!!          end do
!!!!$OMP end do nowait
!!!
!!!
!!!        do im = 2, nm
!!!!$OMP do schedule (dynamic)
!!!          do iv = 1, 2*nv
!!!            do iz = -nz, nz-1
!!!              do my = ist_y, iend_y
!!!                do mx = -nx, nx
!!!                  dfdvp(mx,my,iz,iv,im) =                                     &
!!!                             ( +          ff(mx,my,iz,iv,im-2)                &
!!!                               -  8._DP * ff(mx,my,iz,iv,im-1)                &
!!!                               +  8._DP * ff(mx,my,iz,iv,im+1)                &
!!!                               -          ff(mx,my,iz,iv,im+2)                &
!!!                             ) * cef4(iz)
!!!                end do
!!!              end do
!!!            end do
!!!          end do
!!!!$OMP end do nowait
!!!        end do
!!!
!!!      else   
!!!
!!!        do im = 0, nm
!!!!$OMP do schedule (dynamic)
!!!          do iv = 1, 2*nv
!!!            do iz = -nz, nz-1
!!!              do my = ist_y, iend_y
!!!                do mx = -nx, nx
!!!                  dfdvp(mx,my,iz,iv,im) =                                     &
!!!                             ( +          ff(mx,my,iz,iv,im-2)                &
!!!                               -  8._DP * ff(mx,my,iz,iv,im-1)                &
!!!                               +  8._DP * ff(mx,my,iz,iv,im+1)                &
!!!                               -          ff(mx,my,iz,iv,im+2)                &
!!!                             ) * cef4(iz)
!!!                end do
!!!              end do
!!!            end do
!!!          end do
!!!!$OMP end do nowait
!!!        end do
!!!     
!!!      end if 
!!!
!!!
!!!!$OMP master
!!!                                    ! call fapp_stop("literm_colli_dvp",1317,1)
!!!                                      call clock_end(1317)
!!!!$OMP end master
!!!
!!!
!!!  END SUBROUTINE colli_dfdvp
!!!
!!!
!!!!--------------------------------------
!!!  SUBROUTINE colli_dfdvp6( ff, dfdvp )
!!!!--------------------------------------
!!!!-------------------------------------------------------------------------------
!!!!
!!!!    calculation of df/dv_perp term with 6th order CFD
!!!!
!!!!    by M. Nakata, M. Nunami, April 2014
!!!!
!!!!-------------------------------------------------------------------------------
!!!
!!!    complex(kind=DP), intent(in), &
!!!      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
!!!
!!!    complex(kind=DP), intent(out), &
!!!      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: dfdvp
!!!
!!!    real(kind=DP), dimension(-nz:nz-1) :: cef4
!!!    integer  ::  mx, my, iz, iv, im
!!!
!!!
!!!!$OMP master
!!!                                      call clock_sta(1317)
!!!                                    ! call fapp_start("literm_colli_dvp",1317,1)
!!!!$OMP end master
!!!
!!!
!!!      do iz = -nz, nz-1
!!!        cef4(iz)   = 1._DP / ( 60._DP * dvp(iz) )
!!!      end do
!!!
!!!
!!!      if ( rankm == 0 ) then 
!!!
!!!        im = 0
!!!!$OMP do schedule (dynamic)
!!!          do iv = 1, 2*nv
!!!            do iz = -nz, nz-1
!!!              do my = ist_y, iend_y
!!!                do mx = -nx, nx
!!!                  dfdvp(mx,my,iz,iv,im) = (0._DP, 0._DP) 
!!!                end do
!!!              end do
!!!            end do
!!!          end do
!!!!$OMP end do nowait
!!!
!!!        im = 1
!!!!$OMP do schedule (dynamic)
!!!          do iv = 1, 2*nv
!!!            do iz = -nz, nz-1
!!!              do my = ist_y, iend_y
!!!                do mx = -nx, nx
!!!                  dfdvp(mx,my,iz,iv,im) =                                     &
!!!                             ( - 45._DP * ff(mx,my,iz,iv,im-1)                &
!!!                               +  9._DP * ff(mx,my,iz,iv,im  )                &
!!!                               + 44._DP * ff(mx,my,iz,iv,im+1)                &
!!!                               -  9._DP * ff(mx,my,iz,iv,im+2)                &
!!!                               +          ff(mx,my,iz,iv,im+3)                &
!!!                             ) * cef4(iz)
!!!                end do
!!!              end do
!!!            end do
!!!          end do
!!!!$OMP end do nowait
!!!
!!!        im = 2
!!!!$OMP do schedule (dynamic)
!!!          do iv = 1, 2*nv
!!!            do iz = -nz, nz-1
!!!              do my = ist_y, iend_y
!!!                do mx = -nx, nx
!!!                  dfdvp(mx,my,iz,iv,im) =                                     &
!!!                             ( +  9._DP * ff(mx,my,iz,iv,im-2)                &
!!!                               - 46._DP * ff(mx,my,iz,iv,im-1)                &
!!!                               + 45._DP * ff(mx,my,iz,iv,im+1)                &
!!!                               -  9._DP * ff(mx,my,iz,iv,im+2)                &
!!!                               +          ff(mx,my,iz,iv,im+3)                &
!!!                             ) * cef4(iz)
!!!                end do
!!!              end do
!!!            end do
!!!          end do
!!!!$OMP end do nowait
!!!
!!!
!!!!$OMP do schedule (dynamic)
!!!        do im = 3, nm
!!!          do iv = 1, 2*nv
!!!            do iz = -nz, nz-1
!!!              do my = ist_y, iend_y
!!!                do mx = -nx, nx
!!!                  dfdvp(mx,my,iz,iv,im) =                                     &
!!!                             ( -          ff(mx,my,iz,iv,im-3)                &
!!!                               +  9._DP * ff(mx,my,iz,iv,im-2)                &
!!!                               - 45._DP * ff(mx,my,iz,iv,im-1)                &
!!!                               + 45._DP * ff(mx,my,iz,iv,im+1)                &
!!!                               -  9._DP * ff(mx,my,iz,iv,im+2)                &
!!!                               +          ff(mx,my,iz,iv,im+3)                &
!!!                             ) * cef4(iz)
!!!                end do
!!!              end do
!!!            end do
!!!          end do
!!!        end do
!!!!$OMP end do nowait
!!!
!!!      else   
!!!
!!!!$OMP do schedule (dynamic)
!!!        do im = 0, nm
!!!          do iv = 1, 2*nv
!!!            do iz = -nz, nz-1
!!!              do my = ist_y, iend_y
!!!                do mx = -nx, nx
!!!                  dfdvp(mx,my,iz,iv,im) =                                     &
!!!                             ( -          ff(mx,my,iz,iv,im-3)                &
!!!                               +  9._DP * ff(mx,my,iz,iv,im-2)                &
!!!                               - 45._DP * ff(mx,my,iz,iv,im-1)                &
!!!                               + 45._DP * ff(mx,my,iz,iv,im+1)                &
!!!                               -  9._DP * ff(mx,my,iz,iv,im+2)                &
!!!                               +          ff(mx,my,iz,iv,im+3)                &
!!!                             ) * cef4(iz)
!!!                end do
!!!              end do
!!!            end do
!!!          end do
!!!        end do
!!!!$OMP end do nowait
!!!     
!!!      end if 
!!!
!!!
!!!!$OMP master
!!!                                    ! call fapp_stop("literm_colli_dvp",1317,1)
!!!                                      call clock_end(1317)
!!!!$OMP end master
!!!
!!!
!!!  END SUBROUTINE colli_dfdvp6
!!!
!!!
!!!!--------------------------------------
!!!  SUBROUTINE colli_hhset(hh,phi,ff)
!!!!--------------------------------------
!!!    complex(kind=DP), &
!!!      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: hh
!!!    complex(kind=DP), intent(in), &
!!!      dimension(-nx:nx,0:ny,-nz:nz-1)               :: phi
!!!    complex(kind=DP), intent(in), &
!!!      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
!!!    real(kind=DP) :: dflg
!!!    integer :: mx, my, iz, iv, im, is1
!!!
!!!!$OMP master
!!!                                      call clock_sta(1319)
!!!                                    ! call fapp_start("literm_colli_hwset",1319,1)
!!!!$OMP end master
!!!
!!!     !dflg = real(1-icheck,kind=DP)
!!!      dflg = real(1-icheck,kind=DP) * real(iFLR,kind=DP)
!!!
!!!      is1 = ranks
!!!!$OMP do collapse(2) private(mx,my,iz,iv,im)
!!!      do im = 0, nm
!!!        do iv = 1, 2*nv
!!!          do iz = -nz, nz-1
!!!            do my = ist_y, iend_y
!!!              do mx = -nx, nx
!!!                hh(mx,my,iz,iv,im) = ( ff(mx,my,iz,iv,im) + dflg*sgn(is1)*j0(mx,my,iz,im)*phi(mx,my,iz)   &
!!!                                                                * fmx(iz,iv,im)*Znum(is1)/tau(is1) )      &
!!!                                    * vp(iz,im) * dvp(iz) * dv * twopi
!!!              end do
!!!            end do
!!!          end do
!!!        end do
!!!      end do
!!!!$OMP end do nowait
!!!
!!!!$OMP master
!!!                                    ! call fapp_stop("literm_colli_hwset",1319,1)
!!!                                      call clock_end(1319)
!!!!$OMP end master
!!!
!!!
!!!  END SUBROUTINE colli_hhset
!!!
!!!
!!!!--------------------------------------
!!!  SUBROUTINE colli_wwset(ww)
!!!!--------------------------------------
!!!    complex(kind=DP), &
!!!      dimension(-nx:nx,0:ny,-nz:nz-1,0:ns-1,1:6)        :: ww 
!!!    integer :: mx, my, iz, is1, ii
!!!
!!!!$OMP master
!!!                                      call clock_sta(1319)
!!!                                    ! call fapp_start("literm_colli_hwset",1319,1)
!!!!$OMP end master
!!!
!!!!$OMP do collapse(2) private(mx,my,iz,is1,ii)
!!!    do ii = 1, 6
!!!      do is1 = 0, ns-1
!!!        do iz = -nz, nz-1
!!!          do my = ist_y, iend_y
!!!            do mx = -nx, nx
!!!              ww(mx,my,iz,is1,ii) = ( 0._DP, 0._DP )
!!!            end do 
!!!          end do 
!!!        end do 
!!!      end do 
!!!    end do 
!!!!$OMP end do
!!!
!!!!$OMP master
!!!                                    ! call fapp_stop("literm_colli_hwset",1319,1)
!!!                                      call clock_end(1319)
!!!!$OMP end master
!!!
!!!
!!!  END SUBROUTINE colli_wwset
!!!
!!!
!!!!!--------------------------------------
!!!!  SUBROUTINE colli_full( ff, phi, cf ) ! Analytic derivative of J0*phi
!!!!!--------------------------------------
!!!!!   Sugama collision operator
!!!!
!!!!    complex(kind=DP), intent(inout), &
!!!!      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
!!!!    complex(kind=DP), intent(in), &
!!!!      dimension(-nx:nx,0:ny,-nz:nz-1) :: phi
!!!!    complex(kind=DP), intent(out), &
!!!!      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: cf
!!!!
!!!!    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: vb1, vb2
!!!!    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: mb1, mb2
!!!!    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: cft, cff
!!!!    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: wrkm, moment_ab, moment_ba
!!!!    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: moment_ab_wk, moment_ba_wk
!!!!    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: dfdvp
!!!!    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: hh2
!!!!    integer  ::  mx, my, iz, iv, im, ii, nn
!!!!
!!!!      allocate( vb1(-nx:nx,0:ny,-nz:nz-1,1:2*nvb,0:nm) )
!!!!      allocate( vb2(-nx:nx,0:ny,-nz:nz-1,1:2*nvb,0:nm) )
!!!!      allocate( mb1(-nx:nx,0:ny,-nz:nz-1,1:2*nv,1:2*nvb) )
!!!!      allocate( mb2(-nx:nx,0:ny,-nz:nz-1,1:2*nv,1:2*nvb) )
!!!!      allocate( cft(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) )
!!!!      allocate( cff(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) )
!!!!      allocate( wrkm(-nx:nx,0:ny,-nz:nz-1,0:ns-1,1:6) )
!!!!      allocate( moment_ab(-nx:nx,0:ny,-nz:nz-1,0:ns-1,1:6) )
!!!!      allocate( moment_ba(-nx:nx,0:ny,-nz:nz-1,0:ns-1,1:6) )
!!!!      allocate( moment_ab_wk(1:6,-nx:nx,0:ny,-nz:nz-1,0:ns-1) )
!!!!      allocate( moment_ba_wk(1:6,-nx:nx,0:ny,-nz:nz-1,0:ns-1) )
!!!!      allocate( dfdvp(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) )
!!!!      allocate( hh2(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) )
!!!!
!!!!!$OMP parallel default (none) &
!!!!!$OMP shared(ff,phi,dfdvp) &
!!!!!$OMP shared(mb1,mb2,vb1,vb2) &
!!!!!$OMP shared(wrkm,moment_ab,moment_ba,hh2) &
!!!!!$OMP shared(cft,cff) &
!!!!!$OMP private(im,ii) 
!!!!
!!!!      call bndry_shifts_m_buffin ( ff, mb1, mb2 )
!!!!      call colli_hhset(hh2,phi,ff)
!!!!      call colli_wwset(wrkm)
!!!!!$OMP barrier
!!!!
!!!!!----------------------------------------------------- ovlp1
!!!!!$OMP master
!!!!      call bndry_shifts_m_sendrecv ( mb1, mb2 )
!!!!!$OMP end master
!!!!      call colli_moment_calc( hh2, phi, wrkm )
!!!!      call colli_zeroset( cff )
!!!!
!!!!!$OMP barrier
!!!!!-----------------------------------------------------
!!!!
!!!!!----------------------------------------------------- non-ovlp part
!!!!      call bndry_shifts_m_buffout ( mb2, ff )
!!!!!$OMP barrier
!!!!!-----------------------------------------------------
!!!!
!!!!!----------------------------------------------------- ovlp2
!!!!!$OMP master
!!!!      call colli_moment_redc( wrkm, moment_ab )
!!!!!$OMP end master
!!!!
!!!!!!!      call colli_dfdvp6( ff, dfdvp )  ! 6th-order CFD
!!!!      call colli_dfdvp( ff, dfdvp )  ! 4th-order CFD
!!!!!$OMP barrier
!!!!!-----------------------------------------------------
!!!!
!!!!!----------------------------------------------------- non-ovlp part
!!!!      do im = 0, nm
!!!!        call bndry_shifts_v_buffin ( dfdvp(:,:,:,:,im), vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
!!!!      end do
!!!!!$OMP barrier
!!!!!-----------------------------------------------------
!!!!
!!!!!----------------------------------------------------- ovlp3
!!!!      im = 0
!!!!!$OMP master
!!!!        call bndry_shifts_v_sendrecv ( vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
!!!!!$OMP end master
!!!!
!!!! ! --- No calculations appear here in f0.52 (Nakata July2015)
!!!! !!!      call colli_zeroset( cdt )
!!!!
!!!!!$OMP barrier
!!!!
!!!!      do im = 1, nm
!!!!!$OMP master
!!!!        call bndry_shifts_v_sendrecv ( vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
!!!!!$OMP end master
!!!!
!!!!        call bndry_shifts_v_buffout ( vb2(:,:,:,:,im-1), dfdvp(:,:,:,:,im-1) )
!!!!!$OMP barrier
!!!!      end do
!!!!
!!!!        call bndry_shifts_v_buffout ( vb2(:,:,:,:,nm), dfdvp(:,:,:,:,nm) )
!!!!
!!!!!-----------------------------------------------------
!!!!!$OMP end parallel
!!!!
!!!!
!!!!!$OMP parallel default (none) &
!!!!!$OMP shared(ff) &
!!!!!$OMP shared(cft,cff) &
!!!!!$OMP shared(phi,dfdvp,moment_ab,moment_ba,wrkm,moment_ab_wk,moment_ba_wk) &
!!!!!$OMP shared(vb1,vb2,mb1,mb2) &
!!!!!$OMP shared(rankw,ist_y,iend_y) &
!!!!!$OMP private(mx,my,iz,iv,im,nn) 
!!!!
!!!!!----------------------------------------------------- ovlp4
!!!!!$OMP master
!!!!      call colli_comm_alltoall( moment_ab, moment_ba )
!!!!!$OMP end master
!!!!
!!!!      do im = 0, nm
!!!!        call bndry_shifts_v_buffin ( ff(:,:,:,:,im), vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
!!!!      end do
!!!!!$OMP barrier
!!!!
!!!!!$OMP do collapse(2) private(mx,my,iz,nn)
!!!!do nn = 0, ns-1
!!!!  do iz = -nz, nz-1
!!!!    do my = ist_y, iend_y
!!!!      do mx = -nx, nx
!!!!        moment_ab_wk(1,mx,my,iz,nn) = moment_ab(mx,my,iz,nn,1)
!!!!        moment_ab_wk(2,mx,my,iz,nn) = moment_ab(mx,my,iz,nn,2)
!!!!        moment_ab_wk(3,mx,my,iz,nn) = moment_ab(mx,my,iz,nn,3)
!!!!        moment_ab_wk(4,mx,my,iz,nn) = moment_ab(mx,my,iz,nn,4)
!!!!        moment_ab_wk(5,mx,my,iz,nn) = moment_ab(mx,my,iz,nn,5)
!!!!        moment_ab_wk(6,mx,my,iz,nn) = moment_ab(mx,my,iz,nn,6)
!!!!
!!!!        moment_ba_wk(1,mx,my,iz,nn) = moment_ba(mx,my,iz,nn,1)
!!!!        moment_ba_wk(2,mx,my,iz,nn) = moment_ba(mx,my,iz,nn,2)
!!!!        moment_ba_wk(3,mx,my,iz,nn) = moment_ba(mx,my,iz,nn,3)
!!!!        moment_ba_wk(4,mx,my,iz,nn) = moment_ba(mx,my,iz,nn,4)
!!!!        moment_ba_wk(5,mx,my,iz,nn) = moment_ba(mx,my,iz,nn,5)
!!!!        moment_ba_wk(6,mx,my,iz,nn) = moment_ba(mx,my,iz,nn,6)
!!!!      enddo
!!!!    enddo
!!!!  enddo
!!!!enddo
!!!!!$OMP enddo
!!!!!-----------------------------------------------------
!!!!
!!!!!----------------------------------------------------- ovlp5
!!!!      im = 0
!!!!!$OMP master
!!!!        call bndry_shifts_v_sendrecv ( vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
!!!!!$OMP end master
!!!!
!!!!!$OMP barrier
!!!!!----------------------------------------------------- 
!!!!
!!!!!----------------------------------------------------- ovlp6
!!!!      im = 1
!!!!!$OMP master
!!!!        call bndry_shifts_v_sendrecv ( vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
!!!!!$OMP end master
!!!!
!!!!        call bndry_shifts_v_buffout ( vb2(:,:,:,:,im-1), ff(:,:,:,:,im-1) )
!!!!!$OMP barrier
!!!!!----------------------------------------------------- 
!!!!
!!!!!----------------------------------------------------- ovlp7
!!!!      im = 2
!!!!!$OMP master
!!!!        call bndry_shifts_v_sendrecv ( vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
!!!!!$OMP end master
!!!!
!!!!        call bndry_shifts_v_buffout ( vb2(:,:,:,:,im-1), ff(:,:,:,:,im-1) )
!!!!!!!        call colli_GK_CT6( ff, phi, dfdvp(:,:,:,:,im-2), im-2, cft(:,:,:,:,im-2) )  ! 6th-order CFD
!!!!        call colli_GK_CT( ff, phi, dfdvp(:,:,:,:,im-2), im-2, cft(:,:,:,:,im-2) )  ! 4th-order CFD 
!!!!        call colli_GK_CF_DT( moment_ba_wk, moment_ab_wk, im-2, cff(:,:,:,:,im-2) )
!!!!!$OMP barrier
!!!!!----------------------------------------------------- 
!!!!
!!!!!----------------------------------------------------- ovlp8
!!!!      do im = 3, nm
!!!!!$OMP master
!!!!        call bndry_shifts_v_sendrecv ( vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
!!!!!$OMP end master
!!!!
!!!!        call bndry_shifts_v_buffout ( vb2(:,:,:,:,im-1), ff(:,:,:,:,im-1) )
!!!!!!!        call colli_GK_CT6( ff, phi, dfdvp(:,:,:,:,im-2), im-2, cft(:,:,:,:,im-2) )  ! 6th-order CFD
!!!!        call colli_GK_CT( ff, phi, dfdvp(:,:,:,:,im-2), im-2, cft(:,:,:,:,im-2) )  ! 4th-order CFD
!!!!        call colli_GK_CF_DT( moment_ba_wk, moment_ab_wk, im-2, cff(:,:,:,:,im-2) )
!!!!!$OMP barrier
!!!!      end do
!!!!!----------------------------------------------------- 
!!!!
!!!!      call bndry_shifts_v_buffout ( vb2(:,:,:,:,nm), ff(:,:,:,:,nm) )
!!!!!!!      call colli_GK_CT6( ff, phi, dfdvp(:,:,:,:,nm-1), nm-1, cft(:,:,:,:,nm-1) )  ! 6th-order CFD
!!!!      call colli_GK_CT( ff, phi, dfdvp(:,:,:,:,nm-1), nm-1, cft(:,:,:,:,nm-1) )  ! 4th-order CFD
!!!!      call colli_GK_CF_DT( moment_ba_wk, moment_ab_wk, nm-1, cff(:,:,:,:,nm-1) )
!!!!!$OMP barrier
!!!!
!!!!!!!      call colli_GK_CT6( ff, phi, dfdvp(:,:,:,:,nm), nm, cft(:,:,:,:,nm) )  ! 6th-order CFD
!!!!      call colli_GK_CT( ff, phi, dfdvp(:,:,:,:,nm), nm, cft(:,:,:,:,nm) )  ! 4th-order CFD
!!!!      call colli_GK_CF_DT( moment_ba_wk, moment_ab_wk, nm, cff(:,:,:,:,nm) )
!!!!!$OMP barrier
!!!!
!!!!!$OMP end parallel
!!!!
!!!!!$OMP parallel workshare
!!!!      cf(:,:,:,:,:) = cft(:,:,:,:,:) + cff(:,:,:,:,:)
!!!!!$OMP end parallel workshare
!!!!
!!!!      deallocate( vb1 )
!!!!      deallocate( vb2 )
!!!!      deallocate( mb1 )
!!!!      deallocate( mb2 )
!!!!      deallocate( cft )
!!!!      deallocate( cff )
!!!!      deallocate( wrkm )
!!!!      deallocate( moment_ab )
!!!!      deallocate( moment_ba )
!!!!      deallocate( moment_ab_wk )
!!!!      deallocate( moment_ba_wk )
!!!!      deallocate( dfdvp )
!!!!      deallocate( hh2 )
!!!!
!!!!  END SUBROUTINE colli_full
!!!!--------------------------------------
!!!  SUBROUTINE colli_full( ff, phi, cf ) ! Finite difference of J0*phi
!!!!--------------------------------------
!!!!   Sugama collision operator
!!!
!!!    complex(kind=DP), intent(inout), &
!!!      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
!!!    complex(kind=DP), intent(in), &
!!!      dimension(-nx:nx,0:ny,-nz:nz-1) :: phi
!!!    complex(kind=DP), intent(out), &
!!!      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: cf
!!!
!!!    complex(kind=DP),  &
!!!      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: wff
!!!    complex(kind=DP),  &
!!!      dimension(-nx:nx,0:ny,-nz:nz-1) :: wphi
!!!    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: vb1, vb2
!!!    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: mb1, mb2
!!!    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: cft, cff
!!!    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: wrkm, moment_ab, moment_ba
!!!    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: moment_ab_wk, moment_ba_wk
!!!    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: dfdvp
!!!    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: hh2
!!!    integer  ::  mx, my, iz, iv, im, ii, nn
!!!
!!!      allocate( vb1(-nx:nx,0:ny,-nz:nz-1,1:2*nvb,0:nm) )
!!!      allocate( vb2(-nx:nx,0:ny,-nz:nz-1,1:2*nvb,0:nm) )
!!!      allocate( mb1(-nx:nx,0:ny,-nz:nz-1,1:2*nv,1:2*nvb) )
!!!      allocate( mb2(-nx:nx,0:ny,-nz:nz-1,1:2*nv,1:2*nvb) )
!!!      allocate( cft(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) )
!!!      allocate( cff(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) )
!!!      allocate( wrkm(-nx:nx,0:ny,-nz:nz-1,0:ns-1,1:6) )
!!!      allocate( moment_ab(-nx:nx,0:ny,-nz:nz-1,0:ns-1,1:6) )
!!!      allocate( moment_ba(-nx:nx,0:ny,-nz:nz-1,0:ns-1,1:6) )
!!!      allocate( moment_ab_wk(1:6,-nx:nx,0:ny,-nz:nz-1,0:ns-1) )
!!!      allocate( moment_ba_wk(1:6,-nx:nx,0:ny,-nz:nz-1,0:ns-1) )
!!!      allocate( dfdvp(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) )
!!!      allocate( hh2(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) )
!!!
!!!!$OMP parallel do collapse(2)
!!!      do im = 0, nm
!!!        do iv = 1, 2*nv
!!!          do iz = -nz, nz-1
!!!            do my = ist_y, iend_y
!!!              do mx = -nx, nx
!!!                wff(mx,my,iz,iv,im) = ff(mx,my,iz,iv,im) + sgn(ranks) * Znum(ranks) &
!!!                          * fmx(iz,iv,im) / tau(ranks) * j0(mx,my,iz,im) * phi(mx,my,iz) &
!!!                                                                     * real(iFLR, kind=DP)
!!!              end do
!!!            end do
!!!          end do
!!!        end do
!!!      end do
!!!      wphi(:,:,:) = (0._DP, 0._DP)
!!!
!!!!$OMP parallel default (none) &
!!!!$OMP shared(ff,wff,phi,dfdvp) &
!!!!$OMP shared(mb1,mb2,vb1,vb2) &
!!!!$OMP shared(wrkm,moment_ab,moment_ba,hh2) &
!!!!$OMP shared(cft,cff) &
!!!!$OMP private(im,ii) 
!!!
!!!      call bndry_shifts_m_buffin ( wff, mb1, mb2 )
!!!      call colli_hhset(hh2,phi,ff)
!!!      call colli_wwset(wrkm)
!!!!$OMP barrier
!!!
!!!!----------------------------------------------------- ovlp1
!!!!$OMP master
!!!      call bndry_shifts_m_sendrecv ( mb1, mb2 )
!!!!$OMP end master
!!!      call colli_moment_calc( hh2, phi, wrkm )
!!!      call colli_zeroset( cff )
!!!
!!!!$OMP barrier
!!!!-----------------------------------------------------
!!!
!!!!----------------------------------------------------- non-ovlp part
!!!      call bndry_shifts_m_buffout ( mb2, wff )
!!!!$OMP barrier
!!!!-----------------------------------------------------
!!!
!!!!----------------------------------------------------- ovlp2
!!!!$OMP master
!!!      call colli_moment_redc( wrkm, moment_ab )
!!!!$OMP end master
!!!
!!!!!!      call colli_dfdvp6( ff, dfdvp )  ! 6th-order CFD
!!!      call colli_dfdvp( wff, dfdvp )  ! 4th-order CFD
!!!!$OMP barrier
!!!!-----------------------------------------------------
!!!
!!!!----------------------------------------------------- non-ovlp part
!!!      do im = 0, nm
!!!        call bndry_shifts_v_buffin ( dfdvp(:,:,:,:,im), vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
!!!      end do
!!!!$OMP barrier
!!!!-----------------------------------------------------
!!!
!!!!----------------------------------------------------- ovlp3
!!!      im = 0
!!!!$OMP master
!!!        call bndry_shifts_v_sendrecv ( vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
!!!!$OMP end master
!!!
!!! ! --- No calculations appear here in f0.52 (Nakata July2015)
!!! !!!      call colli_zeroset( cdt )
!!!
!!!!$OMP barrier
!!!
!!!      do im = 1, nm
!!!!$OMP master
!!!        call bndry_shifts_v_sendrecv ( vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
!!!!$OMP end master
!!!
!!!        call bndry_shifts_v_buffout ( vb2(:,:,:,:,im-1), dfdvp(:,:,:,:,im-1) )
!!!!$OMP barrier
!!!      end do
!!!
!!!        call bndry_shifts_v_buffout ( vb2(:,:,:,:,nm), dfdvp(:,:,:,:,nm) )
!!!
!!!!-----------------------------------------------------
!!!!$OMP end parallel
!!!
!!!
!!!!$OMP parallel default (none) &
!!!!$OMP shared(wff) &
!!!!$OMP shared(cft,cff) &
!!!!$OMP shared(wphi,dfdvp,moment_ab,moment_ba,wrkm,moment_ab_wk,moment_ba_wk) &
!!!!$OMP shared(vb1,vb2,mb1,mb2) &
!!!!$OMP shared(rankw,ist_y,iend_y) &
!!!!$OMP private(mx,my,iz,iv,im,nn) 
!!!
!!!!----------------------------------------------------- ovlp4
!!!!$OMP master
!!!      call colli_comm_alltoall( moment_ab, moment_ba )
!!!!$OMP end master
!!!
!!!      do im = 0, nm
!!!        call bndry_shifts_v_buffin ( wff(:,:,:,:,im), vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
!!!      end do
!!!!$OMP barrier
!!!
!!!!$OMP do collapse(2) private(mx,my,iz,nn)
!!!do nn = 0, ns-1
!!!  do iz = -nz, nz-1
!!!    do my = ist_y, iend_y
!!!      do mx = -nx, nx
!!!        moment_ab_wk(1,mx,my,iz,nn) = moment_ab(mx,my,iz,nn,1)
!!!        moment_ab_wk(2,mx,my,iz,nn) = moment_ab(mx,my,iz,nn,2)
!!!        moment_ab_wk(3,mx,my,iz,nn) = moment_ab(mx,my,iz,nn,3)
!!!        moment_ab_wk(4,mx,my,iz,nn) = moment_ab(mx,my,iz,nn,4)
!!!        moment_ab_wk(5,mx,my,iz,nn) = moment_ab(mx,my,iz,nn,5)
!!!        moment_ab_wk(6,mx,my,iz,nn) = moment_ab(mx,my,iz,nn,6)
!!!
!!!        moment_ba_wk(1,mx,my,iz,nn) = moment_ba(mx,my,iz,nn,1)
!!!        moment_ba_wk(2,mx,my,iz,nn) = moment_ba(mx,my,iz,nn,2)
!!!        moment_ba_wk(3,mx,my,iz,nn) = moment_ba(mx,my,iz,nn,3)
!!!        moment_ba_wk(4,mx,my,iz,nn) = moment_ba(mx,my,iz,nn,4)
!!!        moment_ba_wk(5,mx,my,iz,nn) = moment_ba(mx,my,iz,nn,5)
!!!        moment_ba_wk(6,mx,my,iz,nn) = moment_ba(mx,my,iz,nn,6)
!!!      enddo
!!!    enddo
!!!  enddo
!!!enddo
!!!!$OMP enddo
!!!!-----------------------------------------------------
!!!
!!!!----------------------------------------------------- ovlp5
!!!      im = 0
!!!!$OMP master
!!!        call bndry_shifts_v_sendrecv ( vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
!!!!$OMP end master
!!!
!!!!$OMP barrier
!!!!----------------------------------------------------- 
!!!
!!!!----------------------------------------------------- ovlp6
!!!      im = 1
!!!!$OMP master
!!!        call bndry_shifts_v_sendrecv ( vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
!!!!$OMP end master
!!!
!!!        call bndry_shifts_v_buffout ( vb2(:,:,:,:,im-1), wff(:,:,:,:,im-1) )
!!!!$OMP barrier
!!!!----------------------------------------------------- 
!!!
!!!!----------------------------------------------------- ovlp7
!!!      im = 2
!!!!$OMP master
!!!        call bndry_shifts_v_sendrecv ( vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
!!!!$OMP end master
!!!
!!!        call bndry_shifts_v_buffout ( vb2(:,:,:,:,im-1), wff(:,:,:,:,im-1) )
!!!!!!        call colli_GK_CT6( ff, phi, dfdvp(:,:,:,:,im-2), im-2, cft(:,:,:,:,im-2) )  ! 6th-order CFD
!!!        call colli_GK_CT( wff, wphi, dfdvp(:,:,:,:,im-2), im-2, cft(:,:,:,:,im-2) )  ! 4th-order CFD 
!!!        call colli_GK_CF_DT( moment_ba_wk, moment_ab_wk, im-2, cff(:,:,:,:,im-2) )
!!!!$OMP barrier
!!!!----------------------------------------------------- 
!!!
!!!!----------------------------------------------------- ovlp8
!!!      do im = 3, nm
!!!!$OMP master
!!!        call bndry_shifts_v_sendrecv ( vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
!!!!$OMP end master
!!!
!!!        call bndry_shifts_v_buffout ( vb2(:,:,:,:,im-1), wff(:,:,:,:,im-1) )
!!!!!!        call colli_GK_CT6( ff, phi, dfdvp(:,:,:,:,im-2), im-2, cft(:,:,:,:,im-2) )  ! 6th-order CFD
!!!        call colli_GK_CT( wff, wphi, dfdvp(:,:,:,:,im-2), im-2, cft(:,:,:,:,im-2) )  ! 4th-order CFD
!!!        call colli_GK_CF_DT( moment_ba_wk, moment_ab_wk, im-2, cff(:,:,:,:,im-2) )
!!!!$OMP barrier
!!!      end do
!!!!----------------------------------------------------- 
!!!
!!!      call bndry_shifts_v_buffout ( vb2(:,:,:,:,nm), wff(:,:,:,:,nm) )
!!!!!!      call colli_GK_CT6( ff, phi, dfdvp(:,:,:,:,nm-1), nm-1, cft(:,:,:,:,nm-1) )  ! 6th-order CFD
!!!      call colli_GK_CT( wff, wphi, dfdvp(:,:,:,:,nm-1), nm-1, cft(:,:,:,:,nm-1) )  ! 4th-order CFD
!!!      call colli_GK_CF_DT( moment_ba_wk, moment_ab_wk, nm-1, cff(:,:,:,:,nm-1) )
!!!!$OMP barrier
!!!
!!!!!!      call colli_GK_CT6( ff, phi, dfdvp(:,:,:,:,nm), nm, cft(:,:,:,:,nm) )  ! 6th-order CFD
!!!      call colli_GK_CT( wff, wphi, dfdvp(:,:,:,:,nm), nm, cft(:,:,:,:,nm) )  ! 4th-order CFD
!!!      call colli_GK_CF_DT( moment_ba_wk, moment_ab_wk, nm, cff(:,:,:,:,nm) )
!!!!$OMP barrier
!!!
!!!!$OMP end parallel
!!!
!!!!$OMP parallel workshare
!!!      cf(:,:,:,:,:) = cft(:,:,:,:,:) + cff(:,:,:,:,:)
!!!!$OMP end parallel workshare
!!!
!!!      deallocate( vb1 )
!!!      deallocate( vb2 )
!!!      deallocate( mb1 )
!!!      deallocate( mb2 )
!!!      deallocate( cft )
!!!      deallocate( cff )
!!!      deallocate( wrkm )
!!!      deallocate( moment_ab )
!!!      deallocate( moment_ba )
!!!      deallocate( moment_ab_wk )
!!!      deallocate( moment_ba_wk )
!!!      deallocate( dfdvp )
!!!      deallocate( hh2 )
!!!
!!!  END SUBROUTINE colli_full
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!--------------------------------------
  SUBROUTINE colli_LB_buffin_v2( ff, phi, vb1, mb1 )
!--------------------------------------
!   Lenard-Bernstein model collsion operator

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: phi
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,0:nm,1:2*nvb) :: vb1
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,1:2*nvb) :: mb1

    real(kind=DP) :: cs1
    integer :: mx, my, iz, iv, im

!$OMP master
                                           call clock_sta(1371)
                                         ! call fapp_start("literm_shifts_bufferin",1371,1)
!$OMP end master

      cs1 = sgn(ranks) * Znum(ranks) / tau(ranks) * real(iFLR, kind=DP)

!!TBI!! !$OMP do collapse(3) schedule(dynamic,nchunk_ymb)
      do iv = 1, nvb
        do im = 0, nm
         do iz = -nz, nz-1
          do my = ist_y, iend_y
            do mx = -nx, nx
              vb1(mx,my,iz,im,iv    ) = ff(mx,my,iz,         iv,im) &
                                  + cs1 * fmx(iz,         iv,im) &
                               * phi(mx,my,iz) * j0(mx,my,iz,im)
              vb1(mx,my,iz,im,iv+nvb) = ff(mx,my,iz,2*nv-nvb+iv,im) &
                                  + cs1 * fmx(iz,2*nv-nvb+iv,im) &
                               * phi(mx,my,iz) * j0(mx,my,iz,im)
            end do
          end do
         end do
        end do
      end do
!!TBI!! !$OMP end do nowait

!!TBI!! !$OMP do collapse(3) schedule(dynamic,nchunk_yvb)
      do im = 1, nvb
        do iv = 1, 2*nv
         do iz = -nz, nz-1
          do my = ist_y, iend_y
            do mx = -nx, nx
              mb1(mx,my,iz,iv,im    ) = ff(mx,my,iz,iv,     im-1) &
                                  + cs1 * fmx(iz,iv,     im-1) &
                      * phi(mx,my,iz) * j0(mx,my,iz,     im-1)
              mb1(mx,my,iz,iv,im+nvb) = ff(mx,my,iz,iv,nm-nvb+im) &
                                  + cs1 * fmx(iz,iv,nm-nvb+im) &
                      * phi(mx,my,iz) * j0(mx,my,iz,nm-nvb+im)
            end do
          end do
         end do
        end do
      end do
!!TBI!! !$OMP end do nowait

!$OMP master
                                         ! call fapp_stop("literm_shifts_bufferin",1371,1)
                                           call clock_end(1371)
!$OMP end master

  END SUBROUTINE colli_LB_buffin_v2


!--------------------------------------
  SUBROUTINE colli_LB_buffout_v2( ff, phi, vb2, mb2, gg )
!--------------------------------------
!   Lenard-Bernstein model collsion operator

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: phi
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,0:nm,1:2*nvb) :: vb2
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,1:2*nvb) :: mb2
    complex(kind=DP), intent(inout), &
      dimension(1-nvb:2*nv+nvb,0-nvb:nm+nvb,-nx:nx,0:ny,-nz:nz-1) :: gg

    real(kind=DP) :: cs1
    integer  ::  mx, my, iz, iv, im

!$OMP master
                                           call clock_sta(1373)
                                         ! call fapp_start("literm_shifts_bufferout",1373,1)
!$OMP end master

      cs1 = sgn(ranks) * Znum(ranks) / tau(ranks) * real(iFLR, kind=DP)

!!TBI!! !$OMP do collapse(2) schedule(dynamic,nchunk_xy)
      do iz = -nz, nz-1
      do my = ist_y, iend_y
        do mx = -nx, nx
          do im = 0, nm
            do iv = 1, 2*nv
              gg(iv,im,mx,my,iz) = ff(mx,my,iz,iv,im) &
                        + cs1 * fmx(iz,iv,im) * phi(mx,my,iz) * j0(mx,my,iz,im)
            end do
          end do
        end do
      end do
      end do
!!TBI!! !$OMP end do nowait

      if ( rankv == 0 ) then
!!TBI!! !$OMP do collapse(2) schedule(dynamic,nchunk_xy)
        do iz = -nz, nz-1
        do my = ist_y, iend_y
          do mx = -nx, nx
            do im = 0, nm
              do iv = 1, nvb
                gg(-nvb+iv,im,mx,my,iz) = (0._DP, 0._DP)
                gg(2*nv+iv,im,mx,my,iz) = vb2(mx,my,iz,im,iv+nvb)
              end do
            end do
          end do
        end do
        end do
!!TBI!! !$OMP end do nowait
      else if ( rankv == nprocv-1 ) then
!!TBI!! !$OMP do collapse(2) schedule(dynamic,nchunk_xy)
        do iz = -nz, nz-1
        do my = ist_y, iend_y
          do mx = -nx, nx
            do im = 0, nm
              do iv = 1, nvb
                gg(-nvb+iv,im,mx,my,iz) = vb2(mx,my,iz,im,iv    )
                gg(2*nv+iv,im,mx,my,iz) = (0._DP, 0._DP)
              end do
            end do
          end do
        end do
        end do
!!TBI!! !$OMP end do nowait
      else
!!TBI!! !$OMP do collapse(2) schedule(dynamic,nchunk_xy)
        do iz = -nz, nz-1
        do my = ist_y, iend_y
          do mx = -nx, nx
            do im = 0, nm
              do iv = 1, nvb
                gg(-nvb+iv,im,mx,my,iz) = vb2(mx,my,iz,im,iv    )
                gg(2*nv+iv,im,mx,my,iz) = vb2(mx,my,iz,im,iv+nvb)
              end do
            end do
          end do
        end do
        end do
!!TBI!! !$OMP end do nowait
      end if

      if ( rankm == 0 ) then
!!TBI!! !$OMP do collapse(2) schedule(dynamic,nchunk_xy)
        do iz = -nz, nz-1
        do my = ist_y, iend_y
          do mx = -nx, nx
            do im = 1, nvb
              do iv = 1, 2*nv
                gg(iv,-nvb-1+im,mx,my,iz) = (0._DP, 0._DP)
                gg(iv,nm+im    ,mx,my,iz) = mb2(mx,my,iz,iv,im+nvb)
              end do
            end do
          end do
        end do
        end do
!!TBI!! !$OMP end do nowait
      else if ( rankm == nprocm-1 ) then
!!TBI!! !$OMP do collapse(2) schedule(dynamic,nchunk_xy)
        do iz = -nz, nz-1
        do my = ist_y, iend_y
          do mx = -nx, nx
            do im = 1, nvb
              do iv = 1, 2*nv
                gg(iv,-nvb-1+im,mx,my,iz) = mb2(mx,my,iz,iv,im    )
                gg(iv,nm+im    ,mx,my,iz) = (0._DP, 0._DP)
              end do
            end do
          end do
        end do
        end do
!!TBI!! !$OMP end do nowait
      else
!!TBI!! !$OMP do collapse(2) schedule(dynamic,nchunk_xy)
        do iz = -nz, nz-1
        do my = ist_y, iend_y
          do mx = -nx, nx
            do im = 1, nvb
              do iv = 1, 2*nv
                gg(iv,-nvb-1+im,mx,my,iz) = mb2(mx,my,iz,iv,im    )
                gg(iv,nm+im    ,mx,my,iz) = mb2(mx,my,iz,iv,im+nvb)
              end do
            end do
          end do
        end do
        end do
!!TBI!! !$OMP end do nowait
      end if

!$OMP master
                                         ! call fapp_stop("literm_shifts_bufferout",1373,1)
                                           call clock_end(1373)
!$OMP end master

  END SUBROUTINE colli_LB_buffout_v2


!--------------------------------------
  SUBROUTINE colli_LB_buffout_v3_sx( ff, phi, vb2, mb2, gg )
!--------------------------------------
!   Lenard-Bernstein model collsion operator
!  This is for SX-Aurora TSUBASA (mod by M. Nakata, 20200928)
!  DO NOT use overlap with OMP 

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: phi
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,0:nm,1:2*nvb) :: vb2
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,1:2*nvb) :: mb2
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: gg

    real(kind=DP) :: cs1
    integer  ::  mx, my, iz, iv, im

!$OMP master
                                           call clock_sta(1373)
                                         ! call fapp_start("literm_shifts_bufferout",1373,1)
!$OMP end master

      cs1 = sgn(ranks) * Znum(ranks) / tau(ranks) * real(iFLR, kind=DP)

!!TBI!! !$OMP do collapse(2) schedule(dynamic,nchunk_xy)
      do im = 0, nm
      do iv = 1, 2*nv
        do iz = -nz, nz-1
          do my = ist_y, iend_y
            do mx = -nx, nx
              gg(mx,my,iz,iv,im) = ff(mx,my,iz,iv,im) &
                        + cs1 * fmx(iz,iv,im) * phi(mx,my,iz) * j0(mx,my,iz,im)
            end do
          end do
        end do
      end do
      end do
!!TBI!! !$OMP end do nowait

      if ( rankv == 0 ) then
!!TBI!! !$OMP do collapse(2) schedule(dynamic,nchunk_xy)
        do im = 0, nm
        do iv = 1, nvb
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                gg(mx,my,iz,-nvb+iv,im) = (0._DP, 0._DP)
                gg(mx,my,iz,2*nv+iv,im) = vb2(mx,my,iz,im,iv+nvb)
              end do
            end do
          end do
        end do
        end do
!!TBI!! !$OMP end do nowait
      else if ( rankv == nprocv-1 ) then
!!TBI!! !$OMP do collapse(2) schedule(dynamic,nchunk_xy)
        do im = 0, nm
        do iv = 1, nvb
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                gg(mx,my,iz,-nvb+iv,im) = vb2(mx,my,iz,im,iv    )
                gg(mx,my,iz,2*nv+iv,im) = (0._DP, 0._DP)
              end do
            end do
          end do
        end do
        end do
!!TBI!! !$OMP end do nowait
      else
!!TBI!! !$OMP do collapse(2) schedule(dynamic,nchunk_xy)
        do im = 0, nm
        do iv = 1, nvb
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                gg(mx,my,iz,-nvb+iv,im) = vb2(mx,my,iz,im,iv    )
                gg(mx,my,iz,2*nv+iv,im) = vb2(mx,my,iz,im,iv+nvb)
              end do
            end do
          end do
        end do
        end do
!!TBI!! !$OMP end do nowait
      end if

      if ( rankm == 0 ) then
!!TBI!! !$OMP do collapse(2) schedule(dynamic,nchunk_xy)
        do im = 1, nvb
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                gg(mx,my,iz,iv,-nvb-1+im) = (0._DP, 0._DP)
                gg(mx,my,iz,iv,nm+im) = mb2(mx,my,iz,iv,im+nvb)
              end do
            end do
          end do
        end do
        end do
!!TBI!! !$OMP end do nowait
      else if ( rankm == nprocm-1 ) then
!!TBI!! !$OMP do collapse(2) schedule(dynamic,nchunk_xy)
        do im = 1, nvb
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                gg(mx,my,iz,iv,-nvb-1+im) = mb2(mx,my,iz,iv,im    )
                gg(mx,my,iz,iv,nm+im) = (0._DP, 0._DP)
              end do
            end do
          end do
        end do
        end do
!!TBI!! !$OMP end do nowait
      else
!!TBI!! !$OMP do collapse(2) schedule(dynamic,nchunk_xy)
        do im = 1, nvb
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                gg(mx,my,iz,iv,-nvb-1+im) = mb2(mx,my,iz,iv,im    )
                gg(mx,my,iz,iv,nm+im ) = mb2(mx,my,iz,iv,im+nvb)
              end do
            end do
          end do
        end do
        end do
!!TBI!! !$OMP end do nowait
      end if

!$OMP master
                                         ! call fapp_stop("literm_shifts_bufferout",1373,1)
                                           call clock_end(1373)
!$OMP end master

  END SUBROUTINE colli_LB_buffout_v3_sx


!--------------------------------------
  SUBROUTINE colli_LB_calc_v3_sx( gg, cf )
!--------------------------------------
!   Lenard-Bernstein model collsion operator

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: gg
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: cf

    real(kind=DP) :: nu_s, cef1, cef2, cef3, cef4, cs1, cflr
    integer  ::  mx, my, iz, iv, im


!$OMP master
                                           call clock_sta(1311)
                                         ! call fapp_start("literm_colli",1311,1)
!$OMP end master

! --- Note that nu(ranks) is a bias factor 
      nu_s = nu(ranks)*3._DP*dsqrt(pi)*ctauiv(ranks,ranks)/4._DP

       cs1 = sgn(ranks) * Znum(ranks) / tau(ranks) * real(iFLR, kind=DP)
      cef1 = nu_s / ( 12._DP * dv * dv )
      cef2 = nu_s / ( 12._DP * dv )


      if ( rankm /= 0 ) then

!!TBI!! !$OMP do collapse(2) schedule(dynamic,nchunk_xy)
        do iz = -nz, nz-1
          cef3 = nu_s / ( 12._DP * dvp(iz) * dvp(iz) )
          cef4 = nu_s / ( 12._DP * dvp(iz) )
          cflr = nu_s * Anum(ranks) * tau(ranks)  &
                      / ( Znum(ranks) * omg(iz) )**2 * real(iFLR, kind=DP)


        do im = 0, nm
          do iv = 1, 2*nv
            do my = ist_y, iend_y
              do mx = -nx, nx
                cf(mx,my,iz,iv,im) =                                   &
                          ( -          gg(mx,my,iz,iv+2,im)               &
                            + 16._DP * gg(mx,my,iz,iv+1,im)               &
                            - 30._DP * gg(mx,my,iz,iv  ,im)               &
                            + 16._DP * gg(mx,my,iz,iv-1,im)               &
                            -          gg(mx,my,iz,iv-2,im)               &
                           ) * cef1                                    &
                        + ( -          gg(mx,my,iz,iv+2,im)               &
                            +  8._DP * gg(mx,my,iz,iv+1,im)               &
                            -  8._DP * gg(mx,my,iz,iv-1,im)               &
                            +          gg(mx,my,iz,iv-2,im)               &
                          ) * cef2 * vl(iv)                            &
                        + ( -          gg(mx,my,iz,iv,im+2)               &
                            + 16._DP * gg(mx,my,iz,iv,im+1)               &
                            - 30._DP * gg(mx,my,iz,iv,im  )               &
                            + 16._DP * gg(mx,my,iz,iv,im-1)               &
                            -          gg(mx,my,iz,iv,im-2)               &
                          ) * cef3                                     &
                        + ( -          gg(mx,my,iz,iv,im+2)               &
                            +  8._DP * gg(mx,my,iz,iv,im+1)               &
                            -  8._DP * gg(mx,my,iz,iv,im-1)               &
                            +          gg(mx,my,iz,iv,im-2)               &
                          ) * cef4 * ( vp(iz,im) + 1._DP / vp(iz,im) ) &
                        + nu_s * 3._DP * gg(mx,my,iz,iv,im)               &      
                        - cflr * ksq(mx,my,iz) * gg(mx,my,iz,iv,im)
              end do
            end do

          end do
        end do

        end do
!!TBI!! !$OMP end do nowait

      else ! rankm == 0

!!TBI!! !$OMP do collapse(2) schedule(dynamic,nchunk_xy)
        do iz = -nz, nz-1
          cef3 = nu_s / ( 12._DP * dvp(iz) * dvp(iz) )
          cef4 = nu_s / ( 12._DP * dvp(iz) )
          cflr = nu_s * Anum(ranks) * tau(ranks)  &
                      / ( Znum(ranks) * omg(iz) )**2 * real(iFLR, kind=DP)

          do iv = 1, 2*nv

            im = 0
            do my = ist_y, iend_y
              do mx = -nx, nx
                cf(mx,my,iz,iv,im) =                                   &
                          ( -          gg(mx,my,iz,iv+2,im)               &
                            + 16._DP * gg(mx,my,iz,iv+1,im)               &
                            - 30._DP * gg(mx,my,iz,iv  ,im)               &
                            + 16._DP * gg(mx,my,iz,iv-1,im)               &
                            -          gg(mx,my,iz,iv-2,im)               &
                          ) * cef1                                     &
                        + ( -          gg(mx,my,iz,iv+2,im)               &
                            +  8._DP * gg(mx,my,iz,iv+1,im)               &
                            -  8._DP * gg(mx,my,iz,iv-1,im)               &
                            +          gg(mx,my,iz,iv-2,im)               &
                          ) * cef2 * vl(iv)                            &
                        + ( -          gg(mx,my,iz,iv,im+2)               &
                            + 16._DP * gg(mx,my,iz,iv,im+1)               &
                            - 30._DP * gg(mx,my,iz,iv,im  )               &
                            + 16._DP * gg(mx,my,iz,iv,im+1)               &
                            -          gg(mx,my,iz,iv,im+2)               &
                          ) * cef3 * 2._DP                             &
                        + nu_s * 3._DP * gg(mx,my,iz,iv,im)               &
                        - cflr * ksq(mx,my,iz) * gg(mx,my,iz,iv,im)
              end do
            end do

            im = 1
            do my = ist_y, iend_y
              do mx = -nx, nx
                cf(mx,my,iz,iv,im) =                                   &
                          ( -          gg(mx,my,iz,iv+2,im)               &
                            + 16._DP * gg(mx,my,iz,iv+1,im)               &
                            - 30._DP * gg(mx,my,iz,iv  ,im)               &
                            + 16._DP * gg(mx,my,iz,iv-1,im)               &
                            -          gg(mx,my,iz,iv-2,im)               &
                          ) * cef1                                     &
                        + ( -          gg(mx,my,iz,iv+2,im)               &
                            +  8._DP * gg(mx,my,iz,iv+1,im)               &
                            -  8._DP * gg(mx,my,iz,iv-1,im)               &
                            +          gg(mx,my,iz,iv-2,im)               &
                          ) * cef2 * vl(iv)                            &
                        + ( -          gg(mx,my,iz,iv,im+2)               &
                            + 16._DP * gg(mx,my,iz,iv,im+1)               &
                            - 30._DP * gg(mx,my,iz,iv,im  )               &
                            + 16._DP * gg(mx,my,iz,iv,im-1)               &
                            -          gg(mx,my,iz,iv,im  )               &
                          ) * cef3                                     &
                        + ( -          gg(mx,my,iz,iv,im+2)               &
                            +  8._DP * gg(mx,my,iz,iv,im+1)               &
                            -  8._DP * gg(mx,my,iz,iv,im-1)               &
                            +          gg(mx,my,iz,iv,im  )               &
                          ) * cef4 * ( vp(iz,im) + 1._DP / vp(iz,im) ) &
                        + nu_s * 3._DP * gg(mx,my,iz,iv,im)               &   
                        - cflr * ksq(mx,my,iz) * gg(mx,my,iz,iv,im)
              end do
            end do

            do im = 2, nm
              do my = ist_y, iend_y
                do mx = -nx, nx
                cf(mx,my,iz,iv,im) =                                   &
                          ( -          gg(mx,my,iz,iv+2,im)               &
                            + 16._DP * gg(mx,my,iz,iv+1,im)               &
                            - 30._DP * gg(mx,my,iz,iv  ,im)               &
                            + 16._DP * gg(mx,my,iz,iv-1,im)               &
                            -          gg(mx,my,iz,iv-2,im)               &
                          ) * cef1                                     &
                        + ( -          gg(mx,my,iz,iv+2,im)               &
                            +  8._DP * gg(mx,my,iz,iv+1,im)               &
                            -  8._DP * gg(mx,my,iz,iv-1,im)               &
                            +          gg(mx,my,iz,iv-2,im)               &
                          ) * cef2 * vl(iv)                            &
                        + ( -          gg(mx,my,iz,iv,im+2)               &
                            + 16._DP * gg(mx,my,iz,iv,im+1)               &
                            - 30._DP * gg(mx,my,iz,iv,im  )               &
                            + 16._DP * gg(mx,my,iz,iv,im-1)               &
                            -          gg(mx,my,iz,iv,im-2)               &
                          ) * cef3                                     &
                        + ( -          gg(mx,my,iz,iv,im+2)               &
                            +  8._DP * gg(mx,my,iz,iv,im+1)               &
                            -  8._DP * gg(mx,my,iz,iv,im-1)               &
                            +          gg(mx,my,iz,iv,im-2)               &
                          ) * cef4 * ( vp(iz,im) + 1._DP / vp(iz,im) ) &
                        + nu_s * 3._DP * gg(mx,my,iz,iv,im)               &
                        - cflr * ksq(mx,my,iz) * gg(mx,my,iz,iv,im)
                end do
              end do
            end do

          end do

        end do
!!TBI!! !$OMP end do nowait

      end if

!$OMP master
                                    ! call fapp_stop("literm_colli",1311,1)
                                      call clock_end(1311)
!$OMP end master


  END SUBROUTINE colli_LB_calc_v3_sx


!--------------------------------------
  SUBROUTINE colli_LB_calc_v2( gg, cf )
!--------------------------------------
!   Lenard-Bernstein model collsion operator

    complex(kind=DP), intent(in), &
      dimension(1-nvb:2*nv+nvb,0-nvb:nm+nvb,-nx:nx,0:ny,-nz:nz-1) :: gg
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: cf

    real(kind=DP) :: nu_s, cef1, cef2, cef3, cef4, cs1, cflr
    integer  ::  mx, my, iz, iv, im


!$OMP master
                                           call clock_sta(1311)
                                         ! call fapp_start("literm_colli",1311,1)
!$OMP end master

! --- Note that nu(ranks) is a bias factor 
      nu_s = nu(ranks)*3._DP*dsqrt(pi)*ctauiv(ranks,ranks)/4._DP

       cs1 = sgn(ranks) * Znum(ranks) / tau(ranks) * real(iFLR, kind=DP)
      cef1 = nu_s / ( 12._DP * dv * dv )
      cef2 = nu_s / ( 12._DP * dv )


      if ( rankm /= 0 ) then

!!TBI!! !$OMP do collapse(2) schedule(dynamic,nchunk_xy)
        do iz = -nz, nz-1
          cef3 = nu_s / ( 12._DP * dvp(iz) * dvp(iz) )
          cef4 = nu_s / ( 12._DP * dvp(iz) )
          cflr = nu_s * Anum(ranks) * tau(ranks)  &
                      / ( Znum(ranks) * omg(iz) )**2 * real(iFLR, kind=DP)

        do my = ist_y, iend_y
          do mx = -nx, nx

            do im = 0, nm
              do iv = 1, 2*nv
                cf(mx,my,iz,iv,im) =                                   &
                          ( -          gg(iv+2,im,mx,my,iz)               &
                            + 16._DP * gg(iv+1,im,mx,my,iz)               &
                            - 30._DP * gg(iv  ,im,mx,my,iz)               &
                            + 16._DP * gg(iv-1,im,mx,my,iz)               &
                            -          gg(iv-2,im,mx,my,iz)               &
                           ) * cef1                                    &
                        + ( -          gg(iv+2,im,mx,my,iz)               &
                            +  8._DP * gg(iv+1,im,mx,my,iz)               &
                            -  8._DP * gg(iv-1,im,mx,my,iz)               &
                            +          gg(iv-2,im,mx,my,iz)               &
                          ) * cef2 * vl(iv)                            &
                        + ( -          gg(iv,im+2,mx,my,iz)               &
                            + 16._DP * gg(iv,im+1,mx,my,iz)               &
                            - 30._DP * gg(iv,im  ,mx,my,iz)               &
                            + 16._DP * gg(iv,im-1,mx,my,iz)               &
                            -          gg(iv,im-2,mx,my,iz)               &
                          ) * cef3                                     &
                        + ( -          gg(iv,im+2,mx,my,iz)               &
                            +  8._DP * gg(iv,im+1,mx,my,iz)               &
                            -  8._DP * gg(iv,im-1,mx,my,iz)               &
                            +          gg(iv,im-2,mx,my,iz)               &
                          ) * cef4 * ( vp(iz,im) + 1._DP / vp(iz,im) ) &
                        + nu_s * 3._DP * gg(iv,im,mx,my,iz)               &      
                        - cflr * ksq(mx,my,iz) * gg(iv,im,mx,my,iz)
              end do
            end do

          end do
        end do

        end do
!!TBI!! !$OMP end do nowait

      else ! rankm == 0

!!TBI!! !$OMP do collapse(2) schedule(dynamic,nchunk_xy)
        do iz = -nz, nz-1
          cef3 = nu_s / ( 12._DP * dvp(iz) * dvp(iz) )
          cef4 = nu_s / ( 12._DP * dvp(iz) )
          cflr = nu_s * Anum(ranks) * tau(ranks)  &
                      / ( Znum(ranks) * omg(iz) )**2 * real(iFLR, kind=DP)

        do my = ist_y, iend_y
          do mx = -nx, nx

            im = 0
              do iv = 1, 2*nv
                cf(mx,my,iz,iv,im) =                                   &
                          ( -          gg(iv+2,im,mx,my,iz)               &
                            + 16._DP * gg(iv+1,im,mx,my,iz)               &
                            - 30._DP * gg(iv  ,im,mx,my,iz)               &
                            + 16._DP * gg(iv-1,im,mx,my,iz)               &
                            -          gg(iv-2,im,mx,my,iz)               &
                          ) * cef1                                     &
                        + ( -          gg(iv+2,im,mx,my,iz)               &
                            +  8._DP * gg(iv+1,im,mx,my,iz)               &
                            -  8._DP * gg(iv-1,im,mx,my,iz)               &
                            +          gg(iv-2,im,mx,my,iz)               &
                          ) * cef2 * vl(iv)                            &
                        + ( -          gg(iv,im+2,mx,my,iz)               &
                            + 16._DP * gg(iv,im+1,mx,my,iz)               &
                            - 30._DP * gg(iv,im  ,mx,my,iz)               &
                            + 16._DP * gg(iv,im+1,mx,my,iz)               &
                            -          gg(iv,im+2,mx,my,iz)               &
                          ) * cef3 * 2._DP                             &
                        + nu_s * 3._DP * gg(iv,im,mx,my,iz)               &
                        - cflr * ksq(mx,my,iz) * gg(iv,im,mx,my,iz)
              end do

            im = 1
              do iv = 1, 2*nv
                cf(mx,my,iz,iv,im) =                                   &
                          ( -          gg(iv+2,im,mx,my,iz)               &
                            + 16._DP * gg(iv+1,im,mx,my,iz)               &
                            - 30._DP * gg(iv  ,im,mx,my,iz)               &
                            + 16._DP * gg(iv-1,im,mx,my,iz)               &
                            -          gg(iv-2,im,mx,my,iz)               &
                          ) * cef1                                     &
                        + ( -          gg(iv+2,im,mx,my,iz)               &
                            +  8._DP * gg(iv+1,im,mx,my,iz)               &
                            -  8._DP * gg(iv-1,im,mx,my,iz)               &
                            +          gg(iv-2,im,mx,my,iz)               &
                          ) * cef2 * vl(iv)                            &
                        + ( -          gg(iv,im+2,mx,my,iz)               &
                            + 16._DP * gg(iv,im+1,mx,my,iz)               &
                            - 30._DP * gg(iv,im  ,mx,my,iz)               &
                            + 16._DP * gg(iv,im-1,mx,my,iz)               &
                            -          gg(iv,im  ,mx,my,iz)               &
                          ) * cef3                                     &
                        + ( -          gg(iv,im+2,mx,my,iz)               &
                            +  8._DP * gg(iv,im+1,mx,my,iz)               &
                            -  8._DP * gg(iv,im-1,mx,my,iz)               &
                            +          gg(iv,im  ,mx,my,iz)               &
                          ) * cef4 * ( vp(iz,im) + 1._DP / vp(iz,im) ) &
                        + nu_s * 3._DP * gg(iv,im,mx,my,iz)               &   
                        - cflr * ksq(mx,my,iz) * gg(iv,im,mx,my,iz)
              end do

            do im = 2, nm
              do iv = 1, 2*nv
                cf(mx,my,iz,iv,im) =                                   &
                          ( -          gg(iv+2,im,mx,my,iz)               &
                            + 16._DP * gg(iv+1,im,mx,my,iz)               &
                            - 30._DP * gg(iv  ,im,mx,my,iz)               &
                            + 16._DP * gg(iv-1,im,mx,my,iz)               &
                            -          gg(iv-2,im,mx,my,iz)               &
                          ) * cef1                                     &
                        + ( -          gg(iv+2,im,mx,my,iz)               &
                            +  8._DP * gg(iv+1,im,mx,my,iz)               &
                            -  8._DP * gg(iv-1,im,mx,my,iz)               &
                            +          gg(iv-2,im,mx,my,iz)               &
                          ) * cef2 * vl(iv)                            &
                        + ( -          gg(iv,im+2,mx,my,iz)               &
                            + 16._DP * gg(iv,im+1,mx,my,iz)               &
                            - 30._DP * gg(iv,im  ,mx,my,iz)               &
                            + 16._DP * gg(iv,im-1,mx,my,iz)               &
                            -          gg(iv,im-2,mx,my,iz)               &
                          ) * cef3                                     &
                        + ( -          gg(iv,im+2,mx,my,iz)               &
                            +  8._DP * gg(iv,im+1,mx,my,iz)               &
                            -  8._DP * gg(iv,im-1,mx,my,iz)               &
                            +          gg(iv,im-2,mx,my,iz)               &
                          ) * cef4 * ( vp(iz,im) + 1._DP / vp(iz,im) ) &
                        + nu_s * 3._DP * gg(iv,im,mx,my,iz)               &
                        - cflr * ksq(mx,my,iz) * gg(iv,im,mx,my,iz)
              end do
            end do

          end do
        end do

        end do
!!TBI!! !$OMP end do nowait

      end if

!$OMP master
                                    ! call fapp_stop("literm_colli",1311,1)
                                      call clock_end(1311)
!$OMP end master


  END SUBROUTINE colli_LB_calc_v2


END MODULE GKV_colli
