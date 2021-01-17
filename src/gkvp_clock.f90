MODULE GKV_clock
!-------------------------------------------------------------------------------
!
!    Elapsed time measurements
!
!    Update history of gkvp_clock.f90
!    --------------
!      gkvp_f0.57 (S. Maeyama, Oct 2020)
!        - Version number f0.57 is removed from filename.
!
!-------------------------------------------------------------------------------

  use GKV_header
  use GKV_mpienv

  implicit none

  private

  public   clock_timer, clock_etime, clock_sta, clock_end, clock_reset

  real(kind=DP), dimension(1:2000), save :: sss, eee, elt
  integer, dimension(1:2000), save :: ccount


CONTAINS


!--------------------------------------
  SUBROUTINE clock_timer( isw, iflg )
!--------------------------------------

    integer, intent(in)  :: isw
    integer, intent(out) :: iflg

    real(kind=DP), save  :: sss0, eee0
    real(kind=DP)        :: ttotl


      if( isw == 0 ) then

        if( rank == 0 ) then
          call clock_etime ( sss0 )
        end if

        call MPI_Bcast( sss0, 1, MPI_DOUBLE_PRECISION, 0, &
                      MPI_COMM_WORLD, ierr_mpi )

        iflg = 0

        sss(1:2000) = 0._DP
        eee(1:2000) = 0._DP
        elt(1:2000) = 0._DP
        ccount(1:2000) = 0

      else if( isw == 1 ) then

        if( rank == 0 ) then
          call clock_etime ( eee0 )
        end if

        call MPI_Bcast( eee0, 1, MPI_DOUBLE_PRECISION, 0, &
                      MPI_COMM_WORLD, ierr_mpi )

        if ( eee0-sss0 > e_limit ) then
          write( olog, * ) " # elapsed time is closing to the limit", eee0-sss0
          iflg = 1
        else
          iflg = 0
        end if

      else if( isw == 2 ) then

        if( rank == 0 ) then
          call clock_etime ( eee0 )
        end if

        call MPI_Bcast( eee0, 1, MPI_DOUBLE_PRECISION, 0, &
                      MPI_COMM_WORLD, ierr_mpi )

        ttotl   = eee0 - sss0
        elt(1310) = sum(elt(1311:1319));               ccount(1310) = ccount(1311)
        elt(1350) = elt(1351) + elt(1352) + elt(1353); ccount(1350) = ccount(1351)
        elt(1360) = elt(1361) + elt(1362) + elt(1363); ccount(1360) = ccount(1361)
        elt(1370) = elt(1371) + elt(1372) + elt(1373); ccount(1370) = ccount(1371)
        elt(1380) = elt(1381) + elt(1382) + elt(1383); ccount(1380) = ccount(1381)
       !elt(1420) = elt(1421) + elt(1422) + elt(1423); ccount(1420) = ccount(1421)
       !elt(1440) = elt(1441) + elt(1442) + elt(1443); ccount(1440) = ccount(1441)
        elt(1520) = elt(1521) + elt(1522) + elt(1523); ccount(1520) = ccount(1521)

        write( olog, * ) ""
        write( olog, * ) " ######### elapsed time (sec) and call count #########"
        write( olog, '(a22,f15.5,i15)' ) "   total            = ", ttotl
        write( olog, '(a22,f15.5,i15)' ) "   pre              = ", elt(1), ccount(1)
        write( olog, '(a22,f15.5,i15)' ) "   timesteploop     = ", elt(2), ccount(2)
        write( olog, '(a22,f15.5,i15)' ) "   post             = ", elt(3), ccount(3)
        write( olog, '(a22,f15.5,i15)' ) "   output           = ", elt(10), ccount(10)
        write( olog, '(a22,f15.5,i15)' ) "   rkg              = ", elt(11), ccount(11)
        write( olog, '(a22,f15.5,i15)' ) "   field            = ", elt(12), ccount(12)
        write( olog, '(a22,f15.5,i15)' ) "   literm           = ", elt(13), ccount(13)
        write( olog, '(a22,f15.5,i15)' ) "   nlterm           = ", elt(14), ccount(14)
        write( olog, '(a22,f15.5,i15)' ) "   zfilter          = ", elt(15), ccount(15)
        write( olog, '(a22,f15.5,i15)' ) "   checkp           = ", elt(16), ccount(16)
        write( olog, '(a22,f15.5,i15)' ) "   colliimp         = ", elt(17), ccount(17)
        write( olog, '(a22,f15.5,i15)' ) "   other            = ", elt(2)-elt(10)-elt(11) &
                                         -elt(12)-elt(13)-elt(14)-elt(15)-elt(16)-elt(17), ccount(2)
        write( olog, * ) " #####################################################"

        write( olog, * ) ""
        write( olog, * ) " ######### elapsed time detail (sec) and call count #########"
        write( olog, '(a22,f15.5,i15)' ) "   field:v0moment   = ", elt(1220), ccount(1220)
        write( olog, '(a22,f15.5,i15)' ) "   field:fsrfave    = ", elt(1230), ccount(1230)
        write( olog, '(a22,f15.5,i15)' ) "   field:other      = ", elt(1210), ccount(1210)
        write( olog, '(a22,f15.5,i15)' ) "   literm:colli     = ", elt(1310), ccount(1310)
        write( olog, '(a22,f15.5,i15)' ) "   literm:perp      = ", elt(1320), ccount(1320)
        write( olog, '(a22,f15.5,i15)' ) "   literm:para      = ", elt(1330), ccount(1330)
        write( olog, '(a22,f15.5,i15)' ) "   literm:other     = ", elt(1340), ccount(1340)
        write( olog, '(a22,f15.5,i15)' ) "   literm:boundf    = ", elt(1350), ccount(1350)
        write( olog, '(a22,f15.5,i15)' ) "   literm:shiftv    = ", elt(1360), ccount(1360)
        write( olog, '(a22,f15.5,i15)' ) "   literm:shiftm    = ", elt(1370), ccount(1370)
        write( olog, '(a22,f15.5,i15)' ) "   literm:bounde    = ", elt(1380), ccount(1380)
        write( olog, '(a22,f15.5,i15)' ) "   nlterm:pack      = ", elt(1410), ccount(1410)
        write( olog, '(a22,f15.5,i15)' ) "   nlterm:transpose = ", elt(1420), ccount(1420)
        write( olog, '(a22,f15.5,i15)' ) "   nlterm:realspcal = ", elt(1430), ccount(1430)
        write( olog, '(a22,f15.5,i15)' ) "   nlterm:transpose = ", elt(1440), ccount(1440)
        write( olog, '(a22,f15.5,i15)' ) "   nlterm:unpack    = ", elt(1450), ccount(1450)
        write( olog, '(a22,f15.5,i15)' ) "   zfilter:calc     = ", elt(1510), ccount(1510)
        write( olog, '(a22,f15.5,i15)' ) "   zfilter:comm     = ", elt(1520), ccount(1520)
        write( olog, * ) " ############################################################"

        write( olog, * ) ""
        write( olog, * ) " ######### elapsed time more detail (sec) and call count #########"
        write( olog, '(a29,f15.5,i15)' ) "   literm:colli:ct         = ", elt(1311), ccount(1311)
        write( olog, '(a29,f15.5,i15)' ) "   literm:colli:dt         = ", elt(1312), ccount(1312)
        write( olog, '(a29,f15.5,i15)' ) "   literm:colli:cf         = ", elt(1313), ccount(1313)
        write( olog, '(a29,f15.5,i15)' ) "   literm:colli:mom        = ", elt(1314), ccount(1314)
        write( olog, '(a29,f15.5,i15)' ) "   literm:colli:ar         = ", elt(1315), ccount(1315)
        write( olog, '(a29,f15.5,i15)' ) "   literm:colli:com        = ", elt(1316), ccount(1316)
        write( olog, '(a29,f15.5,i15)' ) "   literm:colli:dvp        = ", elt(1317), ccount(1317)
        write( olog, '(a29,f15.5,i15)' ) "   literm:colli:0          = ", elt(1318), ccount(1318)
        write( olog, '(a29,f15.5,i15)' ) "   literm:colli:hwset      = ", elt(1319), ccount(1319)
        write( olog, '(a29,f15.5,i15)' ) "   literm:boundf:bufferin  = ", elt(1351), ccount(1351)
        write( olog, '(a29,f15.5,i15)' ) "   literm:boundf:sendrecv  = ", elt(1352), ccount(1352)
        write( olog, '(a29,f15.5,i15)' ) "   literm:boundf:bufferout = ", elt(1353), ccount(1353)
        write( olog, '(a29,f15.5,i15)' ) "   literm:shiftv:bufferin  = ", elt(1361), ccount(1361)
        write( olog, '(a29,f15.5,i15)' ) "   literm:shiftv:sendrecv  = ", elt(1362), ccount(1362)
        write( olog, '(a29,f15.5,i15)' ) "   literm:shiftv:bufferout = ", elt(1363), ccount(1363)
        write( olog, '(a29,f15.5,i15)' ) "   literm:shiftm:bufferin  = ", elt(1371), ccount(1371)
        write( olog, '(a29,f15.5,i15)' ) "   literm:shiftm:sendrecv  = ", elt(1372), ccount(1372)
        write( olog, '(a29,f15.5,i15)' ) "   literm:shiftm:bufferout = ", elt(1373), ccount(1373)
        write( olog, '(a29,f15.5,i15)' ) "   literm:bounde:bufferin  = ", elt(1381), ccount(1381)
        write( olog, '(a29,f15.5,i15)' ) "   literm:bounde:sendrecv  = ", elt(1382), ccount(1382)
        write( olog, '(a29,f15.5,i15)' ) "   literm:bounde:bufferout = ", elt(1383), ccount(1383)
        write( olog, '(a29,f15.5,i15)' ) "   nlterm:backward:Xfft    = ", elt(1421), ccount(1421)
        write( olog, '(a29,f15.5,i15)' ) "   nlterm:backward:shiftXY = ", elt(1422), ccount(1422)
        write( olog, '(a29,f15.5,i15)' ) "   nlterm:backward:Yfft    = ", elt(1423), ccount(1423)
        write( olog, '(a29,f15.5,i15)' ) "   nlterm:forward:Yfft     = ", elt(1441), ccount(1441)
        write( olog, '(a29,f15.5,i15)' ) "   nlterm:forward:shiftYX  = ", elt(1442), ccount(1442)
        write( olog, '(a29,f15.5,i15)' ) "   nlterm:forward:Xfft     = ", elt(1443), ccount(1443)
        write( olog, '(a29,f15.5,i15)' ) "   zfilter:comm:bufferin   = ", elt(1521), ccount(1521)
        write( olog, '(a29,f15.5,i15)' ) "   zfilter:comm:sendrecv   = ", elt(1522), ccount(1522)
        write( olog, '(a29,f15.5,i15)' ) "   zfilter:comm:bufferout  = ", elt(1523), ccount(1523)
        write( olog, * ) " #################################################################"

        write( olog, * ) ""
        write( olog, * ) " ######### For implicit collison solver #########"
        write( olog, '(a23,f15.5,i15)' ) "   colliimp:init     = ", elt(1700), ccount(1700)
        write( olog, '(a23,f15.5,i15)' ) "   colliimp:pack     = ", elt(1710), ccount(1710)
        write( olog, '(a23,f15.5,i15)' ) "   colliimp:alltoall = ", elt(1711), ccount(1711)
        write( olog, '(a23,f15.5,i15)' ) "   colliimp:unpack   = ", elt(1712), ccount(1712)
        write( olog, '(a23,f15.5,i15)' ) "   colliimp:bicgstab = ", elt(1720), ccount(1720)
        write( olog, '(a23,f15.5,i15)' ) "   colliimp:pack     = ", elt(1730), ccount(1730)
        write( olog, '(a23,f15.5,i15)' ) "   colliimp:alltoall = ", elt(1731), ccount(1731)
        write( olog, '(a23,f15.5,i15)' ) "   colliimp:unpack   = ", elt(1732), ccount(1732)
        write( olog, * ) " ################################################"

        iflg = 0

      end if


  END SUBROUTINE clock_timer


!--------------------------------------
  SUBROUTINE clock_etime( ttt )
!--------------------------------------

    real(kind=DP) :: ttt


      ttt = MPI_Wtime()


  END SUBROUTINE clock_etime


!--------------------------------------
  SUBROUTINE clock_sta( id )
!--------------------------------------

    integer, intent(in) :: id


      call clock_etime( sss(id) )


  END SUBROUTINE clock_sta


!--------------------------------------
  SUBROUTINE clock_end( id )
!--------------------------------------

    integer, intent(in) :: id


      call clock_etime( eee(id) )

      elt(id) = elt(id) + eee(id) - sss(id)

      ccount(id)=ccount(id)+1


  END SUBROUTINE clock_end


!--------------------------------------
  SUBROUTINE clock_reset
!--------------------------------------


      elt(10:2000) = 0._DP
      ccount(10:2000) = 0


  END SUBROUTINE clock_reset


END MODULE GKV_clock
