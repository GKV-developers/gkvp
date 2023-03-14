MODULE ring_header

  implicit none

  integer, parameter :: DP = selected_real_kind(14)
!
!  integer :: olog = 10
!
!  real(kind=DP) :: ring_a = 0.5_DP

  real(kind=DP) :: ring_a

  public  ring_a


END MODULE ring_header


MODULE ring_func

  use ring_header
  use GKV_math,   only: math_eli1, math_eli2

  implicit none

  private

  public func_k, func_g, func_psi, func_x, func_eli1, func_eli2


CONTAINS


  FUNCTION func_k( r, z )

    real(kind=DP)             :: func_k
    real(kind=DP), intent(in) :: r, z

      func_k = sqrt( ( 4._DP * ring_a * abs(r) ) / ( ( ring_a + r )**2 + z**2 ) )

  END FUNCTION func_k


  FUNCTION func_g( r, z )

    real(kind=DP)             :: func_g
    real(kind=DP), intent(in) :: r, z

    real(kind=DP) :: eli1, eli2, kk

      kk = func_k( r, z )**2

      call math_eli1( kk, eli1 ) 
!!!        if( icon /= 0 ) print *, "# icon in celi1 = ", icon
      call math_eli2( kk, eli2 ) 
!!!        if( icon /= 0 ) print *, "# icon in celi2 = ", icon

      func_g = ( ( 1._DP - 0.5_DP * kk ) * eli1 - eli2 ) &
             * sqrt( ( ring_a + r )**2 + z**2 ) * 0.5_DP

  END FUNCTION func_g


  FUNCTION func_psi( r, z )

    real(kind=DP)             :: func_psi
    real(kind=DP), intent(in) :: r, z

      func_psi = func_g( r, z ) / func_g( 1._DP, 0._DP )

  END FUNCTION func_psi


  FUNCTION func_x( r, z, psi0 )

    real(kind=DP)             :: func_x
    real(kind=DP), intent(in) :: r, z, psi0

      func_x =  (func_psi( 1._DP, 0._DP ) - func_psi( r, z ) ) / psi0

  END FUNCTION func_x


  FUNCTION func_eli1( r, z )

    real(kind=DP)             :: func_eli1
    real(kind=DP), intent(in) :: r, z

    real(kind=DP) :: eli1, kk

      kk = func_k( r, z )**2

      call math_eli1( kk, eli1 ) 

      func_eli1 = eli1

  END FUNCTION func_eli1


  FUNCTION func_eli2( r, z )

    real(kind=DP)             :: func_eli2
    real(kind=DP), intent(in) :: r, z

    real(kind=DP) :: eli2, kk

      kk = func_k( r, z )**2

      call math_eli2( kk, eli2 ) 

      func_eli2 = eli2

  END FUNCTION func_eli2


END MODULE ring_func


MODULE ring_diff

  use ring_header

  implicit none

  private

  public diff_r, diff_z, diff_rho


CONTAINS


  FUNCTION diff_r( fun, rin, zin )

    real(kind=DP), external :: fun
    real(kind=DP), intent(in)  :: rin, zin

    real(kind=DP) :: diff_r

    real(kind=DP) :: dr1, dr2


      dr1   = abs( rin ) * 1.d-4
      dr2   = dr1 * 2._DP

      if( rin == 0._DP ) then
        diff_r = 0._DP

      else
        diff_r = ( - fun(rin+dr2,zin) + 8._DP*fun(rin+dr1,zin)   &
                   - 8._DP*fun(rin-dr1,zin) + fun(rin-dr2,zin) ) &
               / ( 12._DP * dr1 )

      end if


  END FUNCTION diff_r


  FUNCTION diff_z( fun, rin, zin )

    real(kind=DP), external :: fun
    real(kind=DP), intent(in)  :: rin, zin

    real(kind=DP) :: diff_z

    real(kind=DP) :: dz1, dz2


      dz1   = abs( zin ) * 1.d-4
      dz2   = dz1 * 2._DP

!      if( zin == 0._DP ) then
       if( abs(zin) < 1.d-12) then
        diff_z = 0._DP

      else
        diff_z = ( - fun(rin,zin+dz2) + 8._DP*fun(rin,zin+dz1)   &
                   - 8._DP*fun(rin,zin-dz1) + fun(rin,zin-dz2) ) &
               / ( 12._DP * dz1 )

      end if


  END FUNCTION diff_z


  FUNCTION diff_rho( fun, rin, zin )

    real(kind=DP), external :: fun
    real(kind=DP), intent(in)  :: rin, zin

    real(kind=DP) :: diff_rho

!    real(kind=DP) :: rho, tht
    real(kind=DP) :: tht


!      rho = sqrt( ( rin - ring_a )**2 + zin**2 )
      tht = atan2( zin, rin - ring_a ) 

      diff_rho = diff_r( fun, rin, zin ) * cos( tht ) &
               + diff_z( fun, rin, zin ) * sin( tht )


  END FUNCTION diff_rho


END MODULE ring_diff


MODULE ring_bfld

  use ring_header
  use ring_func
  use ring_diff

  implicit none

  private

  public bfld_br, bfld_bz, bfld_magb, bfld_gradbr, bfld_gradbz


CONTAINS


  FUNCTION bfld_br( r, z )

    real(kind=DP)             :: bfld_br
    real(kind=DP), intent(in) :: r, z

!      bfld_br = diff_z( func_psi, r, z) / r

      bfld_br = diff_z( func_psi, r, z ) / r

  END FUNCTION bfld_br


  FUNCTION bfld_bz( r, z )

    real(kind=DP)             :: bfld_bz
    real(kind=DP), intent(in) :: r, z

!      bfld_bz = diff_r( func_psi, r, z ) / r

      bfld_bz = - diff_r( func_psi, r, z ) / r

  END FUNCTION bfld_bz


  FUNCTION bfld_magb( r, z )

    real(kind=DP)             :: bfld_magb
    real(kind=DP), intent(in) :: r, z

    real(kind=DP) :: br, bz

      br = bfld_br( r, z )
      bz = bfld_bz( r, z )

      bfld_magb = sqrt( br**2 + bz**2 )

!!!      print *, br, bz

  END FUNCTION bfld_magb


  FUNCTION bfld_gradbr( r, z )

    real(kind=DP)             :: bfld_gradbr
    real(kind=DP), intent(in) :: r, z

      bfld_gradbr = diff_r( bfld_magb,  r, z )

  END FUNCTION bfld_gradbr


  FUNCTION bfld_gradbz( r, z )

    real(kind=DP)             :: bfld_gradbz
    real(kind=DP), intent(in) :: r, z

      bfld_gradbz = diff_z( bfld_magb,  r, z )

  END FUNCTION bfld_gradbz


END MODULE ring_bfld


MODULE GKV_ring
!-----------------------------------------------------------------
!
!   Flux tube coordinates in the ring dipole geometry
!
!   Definition of the flux tube coordinates
!   in the ring dipole geometry
!
!     x = (Psi_0 - Psi) / (R_0*B'_0)
!     y = R_0 * phi
!     z = Theta (= arctan(Z/(R-a)) 
!     where (R, phi, Z ) are the cylindorical coordinates
!
!     B'_0 = 1/R * (dPsi/dr - dPsi/dz) | r=1,z=0
!     B_0 = B / B'_0                   | r=1,z=0
!   We also use the normalization of Psi_0*(gradx crs grady) = B_0*R_0^2 
!     with R_0 = 1 and B_0 = 1 as the units
!
!-----------------------------------------------------------------

  use ring_header
!  use gkv_header
  use ring_func
  use ring_diff
  use ring_bfld

  implicit none

  private

  public ring_coordinates


CONTAINS

  SUBROUTINE ring_coordinates( a, tht, bb, ub_dot_grdb, ub_crs_grdb, &
                               gxx, gxy, gxz, gyy, gyz, gzz, rootg, dbdx, dbdz)

  real(kind=DP), intent(in)  :: a, tht
  real(kind=DP), intent(out) :: bb, ub_dot_grdb, ub_crs_grdb
  real(kind=DP), intent(out) :: gxx, gxy, gxz, gyy, gyz, gzz, rootg
!>>
  real(kind=DP), intent(out) :: dbdx, dbdz
  real(kind=DP) :: R0, psi0
  real(kind=DP) :: eps_x, rho_p, rho_m, r_p, r_m, z_p, z_m
!<<
  real(kind=DP) :: r, z, rho, rho1, hh, dh

  real(kind=DP) :: b0
  real(kind=DP) :: gbr, gbz, ubr, ubz!, psi_n, b_rootg_i, ub_dot_grdh

!  real(kind=DP) :: rootg2

  real(kind=DP) :: eps, eps0 = 0.00000001_DP
  integer       :: ic, nc = 20


    ring_a = a

    rho  = 1._DP - ring_a

      hh  = 0._DP
      dh  = acos(-1._DP) * 0.01_DP



! compute for hh
!!!      do while ( hh < abs(tht)-dh*0.5 )
      do while ( hh < abs(tht)-dh )
        hh = hh + dh

        eps = 1._DP
        ic  = 0

        r = rho * cos(hh) + ring_a
        z = rho * sin(hh)

      do while ( eps > eps0  .AND.  ic < nc )

        ic = ic + 1

        rho1 = rho - ( func_psi( r, z ) - 1._DP ) / diff_rho( func_psi, r, z )
        eps  = abs( rho1 - rho ) / abs( rho )
        rho  = rho1

        r = rho * cos(hh) + ring_a
        z = rho * sin(hh)
        if( abs(z) < 1.d-14 ) z = 0._DP

      end do
!!!           print *, "# ic, tht, rho1, psi = ", ic, tht, rho, func_psi( r, z )

      end do

! compute for tht
        eps = 1._DP
        ic  = 0

        r = rho * cos(tht) + ring_a
        z = rho * sin(tht)

      do while ( eps > eps0  .AND.  ic < nc )

        ic = ic + 1

        rho1 = rho - ( func_psi( r, z ) - 1._DP ) / diff_rho( func_psi, r, z )
        eps  = abs( rho1 - rho ) / abs( rho )
        rho  = rho1

        r = rho * cos(tht) + ring_a
        z = rho * sin(tht)
        if( abs(z) < 1.d-14 ) z = 0._DP

      end do

        if( ic == nc ) then
           print *, "# ic, tht, rho1, psi = ", ic, tht, rho, func_psi( r, z )
        end if

!>>
        R0 = 1._DP
        eps_x = 0.0000001_DP
        rho_p = rho + eps_x
        r_p = rho_p * cos(tht) + ring_a
        z_p = rho_p * sin(tht)
        rho_m = rho - eps_x
        r_m = rho_m * cos(tht) + ring_a
        z_m = rho_m * sin(tht)
!<<

        b0 = bfld_magb( 1._DP, 0._DP )

        bb  = bfld_magb  ( r, z )/b0

!>>
        psi0 = b0*R0**2
        dbdx = ( bfld_magb( r_p, z_p ) - bfld_magb( r_m, z_m ) )/b0/( func_x(r_p, z_p, psi0) - func_x(r_m, z_m, psi0) )
!<<

        ubr = bfld_br( r, z )/b0 / bb
        ubz = bfld_bz( r, z )/b0 / bb

        gbr = bfld_gradbr( r, z )/b0
        gbz = bfld_gradbz( r, z )/b0

        ub_dot_grdb = ubr * gbr + ubz * gbz
        ub_crs_grdb = ubz * gbr - ubr * gbz

!        ub_dot_grdh = - ubr * sin( tht ) + ubz * cos( tht )

!        psi_n = func_psi( r, z ) / func_psi( 1._DP, 0._DP )

        gxx   = ( r * bb )**2
        gxy   = 0._DP
        gxz   = - r * bb * ( ubz*sin(tht) + ubr*cos(tht) ) / rho
        gyy   = 1._DP / r**2
        gyz   = 0._DP
        gzz   = 1._DP / rho**2
        rootg = 1._DP / sqrt( gyy*( gxx*gzz - gxz**2 ) )

!>>
        dbdz = ub_dot_grdb * bb * rootg
!<<

!        b_rootg_i = 1._DP / ( bb * rootg )

!        rootg2= r * rho * b0 / ( diff_r( func_psi, r, z )*cos(tht)   &
!                                 + diff_z( func_psi, r, z )*sin(tht) )
!        rootg2= 1._DP / ( ubz * bb * cos(tht)   &
!                        - ubr * bb * sin(tht) )

!! debug
!        write(unit=6, fmt="(1p, 32e15.7)" ) tht, rho, r, z, func_psi(r,z), &
!              bb, ubr, ubz, gbr, gbz, ub_dot_grdb, ub_crs_grdb, ub_dot_grdh,  &
!              psi_n, gxx, gxy, gxz, gyy, gyz, gzz, rootg, b_rootg_i, rootg2
!! debug



  END SUBROUTINE ring_coordinates


END MODULE GKV_ring
