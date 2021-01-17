MODULE GKV_shearflow
!-------------------------------------------------------------------------------
!
!    Shearflow convection term  
!
!    Update history of gkvp_shearflow.f90
!    --------------
!      gkvp_f0.57 (S. Maeyama, Oct 2020)
!        - Version number f0.57 is removed from filename.
!      gkvp_f0.55 (M. Nakata, Dec 2018)   
!        - First implementation
!
!-------------------------------------------------------------------------------

  use GKV_header
  use GKV_mpienv
  use GKV_fld,    only: fld_esfield, fld_emfield_ff, fld_ff2hh
  use GKV_tips,  only: tips_reality

  implicit none

  private

  public   shearflow_kxmap


CONTAINS


!--------------------------------------
  SUBROUTINE shearflow_kxmap( time, ff, phi, Al, hh )
!--------------------------------------
!     discrete advection in kx direction due to the mean radial flow shear

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1)                       :: phi, Al
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm)           :: hh
    real(kind=DP), intent(in)           :: time

    complex(kind=DP), dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff_tmp

    integer :: mx, my, iz, iv, im, mx_new, gmy, tloop
    integer, dimension(0:ny) :: my_map, loop_mapping

      tloop = nint(time/dt)
      my_map(:) = -1
      if (rankw == 0) loop_mapping(0) = 1

      do my = ist1_y, iend_y
 
        gmy = my + (ny+1)*rankw
        !!!loop_mapping(my) = nint(kxmin_g/(kymin_g*abs(gamma_e)*dt))/gmy
        loop_mapping(my) = nint(kxmin_g/(kymin_g*gmy*abs(gamma_e)*dt))

        if (mod(tloop+loop_mapping(my),loop_mapping(my)) == 0 ) then 
          my_map(my) = my
        else 
          my_map(my) = -1
        end if 

      end do

      if (maxval(my_map) < 0 ) then 
        return
      else 

        if ( gamma_e > 0._DP ) then

!$OMP parallel do collapse(2) default(none) &
!$OMP shared(ist_y,iend_y,my_map,ff,ff_tmp) &
!$OMP private(mx,my,iz,iv,im,mx_new)
          do im = 0-nvb, nm+nvb
            do iv = 1-nvb, 2*nv+nvb
              do iz = -nz-nzb, nz-1+nzb
                do my = ist_y, iend_y

                  if (my == my_map(my)) then       
                    do mx = -nx+1, nx
                      mx_new = mx - 1  
                      ff_tmp(mx_new,my,iz,iv,im) = ff(mx,my,iz,iv,im)
                    end do 
                    ff_tmp(nx,my,iz,iv,im) = (0._DP, 0._DP)
                  end if

                end do 
              end do 
            end do 
          end do 

        else if (gamma_e < 0._DP) then  

!$OMP parallel do collapse(2) default(none) &
!$OMP shared(ist_y,iend_y,my_map,ff,ff_tmp) &
!$OMP private(mx,my,iz,iv,im,mx_new)
          do im = 0-nvb, nm+nvb
            do iv = 1-nvb, 2*nv+nvb
              do iz = -nz-nzb, nz-1+nzb
                do my = ist_y, iend_y

                  if (my == my_map(my)) then       
                    do mx = -nx, nx-1
                      mx_new = mx + 1  
                      ff_tmp(mx_new,my,iz,iv,im) = ff(mx,my,iz,iv,im)
                    end do 
                    ff_tmp(-nx,my,iz,iv,im) = (0._DP, 0._DP)
                  end if

                end do 
              end do 
            end do 
          end do 
      
        end if 

!$OMP parallel do collapse(2) default(none) &
!$OMP shared(ist_y,iend_y,my_map,ff,ff_tmp) &
!$OMP private(mx,my,iz,iv,im,mx_new)
        do im = 0-nvb, nm+nvb
          do iv = 1-nvb, 2*nv+nvb
            do iz = -nz-nzb, nz-1+nzb
              do my = ist_y, iend_y

                if (my == my_map(my)) then       
                  do mx = -nx, nx
                    ff(mx,my,iz,iv,im) = ff_tmp(mx,my,iz,iv,im)
                  end do
                end if

              end do
            end do
          end do
        end do

        call fld_esfield( ff, phi )
        if ( beta .ne. 0._DP ) then
          call fld_emfield_ff( ff, Al )
        end if
        call fld_ff2hh( ff, Al, hh )

        call tips_reality( hh )

      end if

  END SUBROUTINE shearflow_kxmap

END MODULE GKV_shearflow
