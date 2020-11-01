MODULE GKV_math
!
!  Mathematical functions using SSLII library
!
!  This routine uses the unit number of 99 which 
!    should not be overlapped with others
!
!  T.-H. Watanabe (Feb 2011)
!                 (May 2011) with new functions, 
!                            J1 and zero points of J0
!

  implicit none

  private

  public   math_j0, math_j0zero, math_j1, math_j2, &
           math_i0, math_g0,     math_random

  integer, parameter :: DP = selected_real_kind(14)
  integer, parameter :: ifnc = 99 ! unit number preserved 
                                  ! for this module


CONTAINS


  SUBROUTINE math_j0( x, j0 )
!
!     0th-order Bessel function
!
    real(kind=DP), intent(in)  :: x
    real(kind=DP), intent(out) :: j0
    integer :: ierr

      call dbj0( x, j0, ierr )

    return

  END SUBROUTINE math_j0


  SUBROUTINE math_j0zero( i, j0zero )
!
!     0th-order Bessel function
!
    integer, intent(in)        :: i
    real(kind=DP), intent(out) :: j0zero

    integer, parameter  :: nzero = 4096
    real(kind=DP), save :: j0zeros(nzero)

    integer, save :: isw = 0

      if( i < 1  .or.  i > nzero ) then
        print *, "# range of J0zero is invalid"
        stop
      end if
     
      if( isw == 0 ) then
        include 'Bessel0_Zeros.f90'
      end if

      j0zero = j0zeros(i)
          
    return


  END SUBROUTINE math_j0zero


  SUBROUTINE math_j1( x, j1 )
!
!     1st-order Bessel function
!
    real(kind=DP), intent(in)  :: x
    real(kind=DP), intent(out) :: j1
    integer :: ierr

      call dbj1( x, j1, ierr )

    return

  END SUBROUTINE math_j1


  SUBROUTINE math_j2( x, j2 )
!
!     0th-order Bessel function
!
    real(kind=DP), intent(in)  :: x
    real(kind=DP), intent(out) :: j2
    real(kind=DP)              :: j0, j1
    integer :: ierr

      call dbj0( x, j0, ierr )
      call dbj1( x, j1, ierr )

      if (x /= 0._DP) then 
        j2 = 2._DP*j1/x - j0
      else 
        j2 = 0._DP
      end if

    return

  END SUBROUTINE math_j2


  SUBROUTINE math_i0( x, i0 )
!
!     0th-order modified Bessel function
!
    real(kind=DP), intent(in)  :: x
    real(kind=DP), intent(out) :: i0
    integer :: ierr

      if( 0._DP <= x  .and.  x < 150._DP ) then
        call dbi0( x, i0, ierr )
      else 
        print *, "### math_i0:  x is out of range!"
      end if

    return

  END SUBROUTINE math_i0


  SUBROUTINE math_g0( x, g0 )
!
!     The Gamma_0 function in gyrokinetics
!     defined by G0(x) = I0(x)*exp(-x)
!
    real(kind=DP), intent(in)  :: x
    real(kind=DP), intent(out) :: g0

    real(kind=DP) :: i0

    real(kind=DP)                 :: twopi
    real(kind=DP), dimension(0:5) :: c

      twopi  = 2._DP * 3.141592653589793_DP
     
      c(0) = 1._DP
      c(1) = 0.25_DP
      c(2) = 9._DP / 32._DP
      c(3) = 75._DP / 128._DP
      c(4) = 3675._DP / 2048._DP
      c(5) = 59535._DP / 8192._DP

      if( x < 150._DP ) then
        call math_i0( x, i0 )
        g0 = i0 * exp( -x )
      else 
        g0 = ( c(0)                      &
             + c(1) / ( 2._DP * x )      &
             + c(2) / ( 2._DP * x )**2   &
             + c(3) / ( 2._DP * x )**3   &
             + c(4) / ( 2._DP * x )**4   &
             + c(5) / ( 2._DP * x )**5 ) &
             / sqrt( twopi * x )
      end if

    return

  END SUBROUTINE math_g0


  SUBROUTINE math_random( rr )
!
!     Random number in [0,1]
!

    real(kind=DP), intent(inout), dimension(:) :: rr

    real(kind=4), dimension(:), allocatable :: rr1
    integer(kind=4), save :: iseed 
    integer :: nr, ierr

    data iseed / 211501 /

      nr = size(rr)

      !call ranu2( iseed, rr, nr, ierr )
      !%%% ranu2 is valid for single precision real number %%%
      allocate( rr1(nr) )
      call ranu2( iseed, rr1, nr, ierr )
      rr = real(rr1, kind=DP)
      deallocate( rr1 )

    return

  END SUBROUTINE math_random


END MODULE GKV_math
