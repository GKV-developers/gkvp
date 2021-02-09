MODULE GKV_fld
!-------------------------------------------------------------------------------
!
!    Field solver
!
!    Update history of gkvp_fld.f90
!    --------------
!      gkvp_f0.57 (S. Maeyama, Oct 2020)
!        - Version number f0.57 is removed from filename.
!        - Intent is changed from intent(out) :: phi to intent(inout) :: phi,
!          to keep initialized values at padding iend_y<my.
!
!-------------------------------------------------------------------------------

  use GKV_header
  use GKV_mpienv
  use GKV_intgrl, only: intgrl_fsrf, intgrl_v0_moment, intgrl_v0_moment_ms
  use GKV_clock, only: clock_sta, clock_end

  implicit none

  private

  public   fld_esfield, fld_emfield_ff, fld_emfield_hh, fld_ff2hh, fld_hh2ff


CONTAINS


!--------------------------------------
  SUBROUTINE fld_esfield ( ff, phi )
!--------------------------------------
!     electrostatic field calculation

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: phi

    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: wf
    complex(kind=DP), dimension(:,:,:), allocatable :: nw, ww
    complex(kind=DP), dimension(-nx:nx) :: zf
    integer  ::  mx, my, iz, iv, im


                                           call clock_sta(1210)
                                         ! call fapp_start("esfield_other",1210,1)
      allocate( wf(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) )
      allocate( nw(-nx:nx,0:ny,-nz:nz-1) )
      allocate( ww(-nx:nx,0:ny,-nz:nz-1) )

!$OMP parallel
      do im = 0, nm
!$OMP do
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                wf(mx,my,iz,iv,im) = ff(mx,my,iz,iv,im) * j0(mx,my,iz,im) * sgn(ranks) * fcs(ranks)
              end do
            end do
          end do
        end do
!$OMP end do nowait
      end do
!$OMP end parallel
                                         ! call fapp_stop("esfield_other",1210,1)
                                           call clock_end(1210)

                                           call clock_sta(1220)
                                         ! call fapp_start("esfield_v0moment",1220,1)
      call intgrl_v0_moment_ms ( wf, nw )
                                         ! call fapp_stop("esfield_v0moment",1220,1)
                                           call clock_end(1220)

      if ( rankw == 0 ) then
        do iz = -nz, nz-1
          nw(0,0,iz) = ( 0._DP, 0._DP )       !  zero-zero
        end do
      end if

      if ( ns == 1 ) then
! --- adiabatic model for ITG-ae or ETG-ai
        
        zf = ( 0._DP, 0._DP )

        if ( sgn(0) > 0._DP ) then  ! --- ZF-calculation for ITG-ae 

          ! --- calculation of zonal flow potential
          if ( rankw == 0 ) then

                                           call clock_sta(1210)
                                         ! call fapp_start("esfield_other",1210,1)
            my   = 0
!$OMP parallel do
            do iz = -nz, nz-1
                do mx = -nx, -1
                  ww(mx,my,iz)   = nw(mx,my,iz)                                     &
                               / ( (  1._DP - g0(mx,my,iz)  + tau(0)*tau_ad ) * fctgt(mx) )
                end do
                mx = 0
                  ww(mx,my,iz) = (0._DP, 0._DP)
                do mx = 1, nx
                  ww(mx,my,iz)   = nw(mx,my,iz)                                     &
                               / ( (  1._DP - g0(mx,my,iz)  + tau(0)*tau_ad ) * fctgt(mx) )
                end do
            end do
                                         ! call fapp_stop("esfield_other",1210,1)
                                           call clock_end(1210)
    
                                           call clock_sta(1230)
                                         ! call fapp_start("esfield_fsrf",1230,1)
            call intgrl_fsrf ( ww, zf )
                                         ! call fapp_stop("esfield_fsrf",1230,1)
                                           call clock_end(1230)
    
            zf(0)   = ( 0._DP, 0._DP )

          end if

        end if

                                           call clock_sta(1210)
                                         ! call fapp_start("esfield_other",1210,1)
!$OMP parallel do
        do iz = -nz, nz-1
          do my = ist_y, iend_y
            do mx = -nx, nx
              phi(mx, my,iz) = nw(mx, my,iz) / ( 1._DP - g0(mx, my,iz) + tau(0)*tau_ad )
            end do
          end do
        end do
  
        if ( rankw == 0 ) then
!$OMP parallel do private(my)
          do iz = -nz, nz-1
            my = 0
              do mx = -nx, nx
              phi(mx,my,iz) = ( nw(mx,my,iz) + zf(mx)*tau(0)*tau_ad ) &
                             / ( 1._DP - g0(mx,my,iz) + tau(0)*tau_ad )
              end do
          end do
        end if
                                         ! call fapp_stop("esfield_other",1210,1)
                                           call clock_end(1210)

      else
! --- kinetic model for multi-species

                                           call clock_sta(1210)
                                         ! call fapp_start("esfield_other",1210,1)
!$OMP parallel do
        do iz = -nz, nz-1
          do my = ist_y, iend_y
            do mx = -nx, nx
              phi(mx,my,iz) = nw(mx,my,iz) * fct_poisson(mx,my,iz)
            end do
          end do
        end do

      deallocate( wf )
      deallocate( nw )
      deallocate( ww )
                                         ! call fapp_stop("esfield_other",1210,1)
                                           call clock_end(1210)

      end if

      if ( rankw == 0 ) then
        do iz = -nz, nz-1
          phi(0,0,iz) = ( 0._DP, 0._DP )
        end do
      end if


  END SUBROUTINE fld_esfield


!--------------------------------------
  SUBROUTINE fld_emfield_ff ( ff, Al )
!--------------------------------------
!     magnetic field calculation

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: Al

    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: wf
    complex(kind=DP), dimension(:,:,:), allocatable :: nw
    integer  ::  mx, my, iz, iv, im


                                           call clock_sta(1210)
                                         ! call fapp_start("esfield_other",1210,1)
      allocate( wf(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) )
      allocate( nw(-nx:nx,0:ny,-nz:nz-1) )

!$OMP parallel
      do im = 0, nm
!$OMP do
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                wf(mx,my,iz,iv,im) = ff(mx,my,iz,iv,im) * j0(mx,my,iz,im) * sgn(ranks) * fcs(ranks) &
                                   * sqrt( tau(ranks) / Anum(ranks) ) * vl(iv)
              end do
            end do
          end do
        end do
!$OMP end do nowait
      end do
!$OMP end parallel
                                         ! call fapp_stop("esfield_other",1210,1)
                                           call clock_end(1210)

                                           call clock_sta(1220)
                                         ! call fapp_start("esfield_v0moment",1220,1)
      call intgrl_v0_moment_ms ( wf, nw )
                                         ! call fapp_stop("esfield_v0moment",1220,1)
                                           call clock_end(1220)

      if ( rankw == 0 ) then
        do iz = -nz, nz-1
          nw(0,0,iz) = ( 0._DP, 0._DP )       !  zero-zero
        end do
      end if

                                           call clock_sta(1210)
                                         ! call fapp_start("esfield_other",1210,1)
!$OMP parallel do
        do iz = -nz, nz-1
          do my = ist_y, iend_y
            do mx = -nx, nx

              if ( rankw == 0 .and. mx == 0 .and. my == 0 ) then
                Al(mx,my,iz) = ( 0._DP, 0._DP )
              else
                Al(mx,my,iz) = nw(mx,my,iz) * beta / ksq(mx,my,iz)
              end if

            end do
          end do
        end do

      deallocate( wf )
      deallocate( nw )
                                         ! call fapp_stop("esfield_other",1210,1)
                                           call clock_end(1210)

  END SUBROUTINE fld_emfield_ff


!--------------------------------------
  SUBROUTINE fld_emfield_hh ( hh, Al )
!--------------------------------------
!     magnetic field calculation

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: hh
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: Al

    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: wf
    complex(kind=DP), dimension(:,:,:), allocatable :: nw
    integer  ::  mx, my, iz, iv, im


                                           call clock_sta(1210)
                                         ! call fapp_start("esfield_other",1210,1)
      allocate( wf(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) )
      allocate( nw(-nx:nx,0:ny,-nz:nz-1) )

!$OMP parallel
      do im = 0, nm
!$OMP do
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                wf(mx,my,iz,iv,im) = hh(mx,my,iz,iv,im) * j0(mx,my,iz,im) * sgn(ranks) * fcs(ranks) &
                                   * sqrt( tau(ranks) / Anum(ranks) ) * vl(iv)
              end do
            end do
          end do
        end do
!$OMP end do nowait
      end do
!$OMP end parallel
                                         ! call fapp_stop("esfield_other",1210,1)
                                           call clock_end(1210)

                                           call clock_sta(1220)
                                         ! call fapp_start("esfield_v0moment",1220,1)
      call intgrl_v0_moment_ms ( wf, nw )
                                         ! call fapp_stop("esfield_v0moment",1220,1)
                                           call clock_end(1220)

      if ( rankw == 0 ) then
        do iz = -nz, nz-1
          nw(0,0,iz) = ( 0._DP, 0._DP )       !  zero-zero
        end do
      end if

                                           call clock_sta(1210)
                                         ! call fapp_start("esfield_other",1210,1)
!$OMP parallel do
        do iz = -nz, nz-1
          do my = ist_y, iend_y
            do mx = -nx, nx
              Al(mx,my,iz) = nw(mx,my,iz) * beta * fct_ampere(mx,my,iz)
            end do
          end do
        end do

      deallocate( wf )
      deallocate( nw )
                                         ! call fapp_stop("esfield_other",1210,1)
                                           call clock_end(1210)

      if ( rankw == 0 ) then
        do iz = -nz, nz-1
          Al(0,0,iz) = ( 0._DP, 0._DP )
        end do
      end if


  END SUBROUTINE fld_emfield_hh


!--------------------------------------
  SUBROUTINE fld_ff2hh ( ff, Al, hh )
!--------------------------------------
!     ff -> hh

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1)                       :: Al
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm)           :: hh

    integer :: mx, my, iz, iv, im

!$OMP parallel
      do im = 0, nm
!$OMP do
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                hh(mx,my,iz,iv,im) = ff(mx,my,iz,iv,im)  &
                    + sgn(ranks) * Znum(ranks)  / sqrt( Anum(ranks) * tau(ranks) )  &
                    * fmx(iz,iv,im) * vl(iv) * j0(mx,my,iz,im) * Al(mx,my,iz)
              end do
            end do
          end do
        end do
!$OMP end do nowait
      end do
!$OMP end parallel


  END SUBROUTINE fld_ff2hh


!--------------------------------------
  SUBROUTINE fld_hh2ff ( hh, Al, ff )
!--------------------------------------
!     hh -> ff

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm)           :: hh
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1)                       :: Al
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff

    integer :: mx, my, iz, iv, im

!$OMP parallel
      do im = 0, nm
!$OMP do
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                ff(mx,my,iz,iv,im) = hh(mx,my,iz,iv,im)  &
                      - sgn(ranks) * Znum(ranks)  / sqrt( Anum(ranks) * tau(ranks) )  &
                      * fmx(iz,iv,im) * vl(iv) * j0(mx,my,iz,im) * Al(mx,my,iz)
              end do
            end do
          end do
        end do
!$OMP end do nowait
      end do
!$OMP end parallel


  END SUBROUTINE fld_hh2ff


END MODULE GKV_fld
