MODULE parameters

  implicit none

  integer, parameter :: nprocw = 8, &
                        nprocz = 4, &
                        nprocv = 3, &
                        nprocm = 2, &
                        nprocs = 4

  integer, parameter :: ll  = 1, &   ! process number per node
                        ls1 = 2, &   ! segment length in z
                        ls2 = 2, &   ! segment length in v
                        ls3 = 2, &   ! segment length in m
                        lg1 = 4, &   ! number of segments in z
                        lg2 = 3, &   ! number of segments in v
                        lg3 = 2, &   ! number of segments in m
                        lp2 = 2, &   ! number of species in v
                        lp3 = 2      ! number of species in m
                        
END MODULE parameters


PROGRAM segmented_mapping
!------------------------------------------------------------------------------
!
!   Create rank-map file for 3D torus network
!
!
!                   -----      
!                 /     /      
!                 -----  | ls3   "a segment"  -  Parallel FFT is performed    
!                |    |  |                       in a segment.
!                |    | / ls2                    ( ll*ls1*ls2*ls3 == nprocw )
!                 -----          
!                  ls1         
!
!
!
!           -----    -----    ------       
!         /     /  /     /  /      /       
!         -----  | -----  | ------  |      
!        |     | ||     | ||     |  |      
!      -----    -----    -----   | /       "a segment group for each species"
!    /     /  /     /  /     /  --    lg3  
!    -----  | -----  | -----  |    /          -  Segments are set up so that
!   |    |  ||    |  ||    |  | --  |            communications in zz, vl, mu
!   |    | / |    | / |    | /   |  |            are performed between
!    -----    -----    -----     | /             adjacent segments.
!    /     /  /     /  /     /  --               ( lg1 == nprocz,
!    -----  | -----  | -----  |                    lg2 == nprocv,
!   |    |  ||    |  ||    |  |   lg2              lg3 == nprocm  )
!   |    | / |    | / |    | /             
!    -----    -----    -----               
!              lg1                         
!
!
!
!          ------------------------- 
!        /                         / 
!       /                         / |    
!      /                         /  |     "species set up"                      
!     /                         /   |                                           
!    /                         /   /         -  Segment groups are arranged
!    -------------------------    /             in vl (or mu) direction
!   |                         |  / /            so that spc_comm_world is
!   |                         | / / |           performed in a cross section. 
!   |                         |/ /  |           ( lp2*lp3 == nprocs )
!    -------------------------  /   |
!    /                         /   /
!    -------------------------    /
!   |                         |  /
!   |                         | /
!   |                         |/
!    ------------------------- 
!
!
!
!------------------------------------------------------------------------------

  use parameters

  implicit none

  character(len=9), dimension(0:ll*ls1*ls2*ls3*lg1*lg2*lg3*lp2*lp3-1) :: c1, c2, c3
  integer :: p1, p2, p3
  integer :: h, i, j, k, l, m, n, o, p, wk


    if ( ll*ls1*ls2*ls3*lg1*lg2*lg3*lp2*lp3 .ne. nprocw*nprocz*nprocv*nprocm*nprocs ) then
      write(*,*) "Total process number is illegal."
      stop
    else if ( ll*ls1*ls2*ls3 .ne. nprocw ) then
      write(*,*) "Segment size is illegal."
      stop
    else if ( lg1 .ne. nprocz ) then
      write(*,*) "Number of segments is illegal."
      stop
    else if ( lg2 .ne. nprocv ) then
      write(*,*) "Number of segments is illegal."
      stop
    else if ( lg3 .ne. nprocm ) then
      write(*,*) "Number of segments is illegal."
      stop
    else if ( lp2*lp3 .ne. nprocs ) then
      write(*,*) "Number of species is illegal."
      stop
    end if


    open(10,file="for_gnuplot.dat")

      wk = 0
      do p = 0, lp3-1    ! mapping for each species
        do o = 0, lp2-1
  
          do n = 0, lg3-1    ! mapping of segments
            do m = 0, lg2-1
              do l = 0, lg1-1
  
                do k = 0, ls3-1    ! mapping in a segment
                  do j = 0, ls2-1
                    do i = 0, ls1-1
                      do h = 1, ll
                        p1 = i+ls1*l
                        p2 = j+ls2*m+ls2*lg2*o
                        p3 = k+ls3*n+ls3*lg3*p
                        write(c1(wk),'(i0)') p1
                        write(c2(wk),'(i0)') p2
                        write(c3(wk),'(i0)') p3
                        write(10,*) p1, p2, p3, "# rank =", wk
                        wk = wk + 1
                      end do
                    end do
                  !  write(10,*)
                  end do
                  write(10,*)
                end do
                write(10,*)
  
              end do
              write(10,*)
            end do
            write(10,*)
          end do
          write(10,*)
  
        end do
      end do

    close(10)


    open(11,file="rankmapfile.dat")

      do wk = 0, ll*ls1*ls2*ls3*lg1*lg2*lg3*lp2*lp3-1
        write(11,'(a)') "("//trim(c1(wk))//","//trim(c2(wk))//","//trim(c3(wk))//")"
      end do

    close(11)


END PROGRAM segmented_mapping

