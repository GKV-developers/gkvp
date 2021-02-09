MODULE GKV_trans
!-------------------------------------------------------------------------------
!
!    Entropy transfer diagnostics
!
!    Update history of gkvp_trans.f90
!    --------------
!      gkvp_f0.60 (S. Maeyama, Jan 2021)
!        - Use fileio module to switch Fortran/NetCDF binary output.
!      gkvp_f0.57 (S. Maeyama, Oct 2020)
!        - Version number f0.57 is removed from filename.
!
!-------------------------------------------------------------------------------

  use GKV_header
  use GKV_mpienv
  use GKV_intgrl, only: intgrl_thet, intgrl_v0_moment
  use GKV_clock, only: clock_sta, clock_end
  use GKV_exb, only: exb_NL_term
  !fj start 202010
  use GKV_fileio
  !fj end 202010

  implicit none

  private

  integer, save, &
    dimension(:), allocatable :: triad_diag_mxt, triad_diag_myt
  integer, save :: nbuff

  public   trans_sum, trans_triad


CONTAINS


!--------------------------------------
  SUBROUTINE trans_sum ( ff, phi, Al, neint, nmint )
!--------------------------------------
!     Check the entropy balance equation

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: phi, Al
    real(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny) :: neint, nmint

    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: wf, wef, wbf
    complex(kind=DP), dimension(:,:,:,:), allocatable :: psi, chi
    complex(kind=DP), dimension(:,:,:), allocatable :: wc3
    complex(kind=DP), dimension(-nx:nx,0:ny) :: wc2
    integer  ::  mx, my, iz, iv, im


      allocate( wf(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) )
      allocate( wef(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) )
      allocate( wbf(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) )
      allocate( psi(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,0:nm) )
      allocate( chi(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,0:nm) )
      allocate( wc3(-nx:nx,0:ny,-nz:nz-1) )

!$OMP parallel workshare
      wf(:,:,:,:,:) = (0._DP, 0._DP)
      psi(:,:,:,:) = (0._DP, 0._DP)
      chi(:,:,:,:) = (0._DP, 0._DP)
      neint(:,:) = 0._DP
      nmint(:,:) = 0._DP
!$OMP end parallel workshare

!$OMP parallel
      do im = 0, nm
!$OMP do
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                psi(mx,my,iz,im) = j0(mx,my,iz,im) * phi(mx,my,iz)
              end do
            end do
          end do
!$OMP end do nowait
      end do

      do im = 0, nm
!$OMP do
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                wf(mx,my,iz,iv,im) = ff(mx,my,iz,iv,im) + sgn(ranks) * Znum(ranks) &
                          * fmx(iz,iv,im) / tau(ranks) * psi(mx,my,iz,im)
              end do
            end do
          end do
        end do
!$OMP end do nowait
      end do
!$OMP end parallel

      call exb_NL_term( wf, psi, chi, wef )  ! Electrostatic part

!$OMP parallel
      do im = 0, nm
!$OMP do
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                chi(mx,my,iz,im) = j0(mx,my,iz,im) * Al(mx,my,iz)
              end do
            end do
          end do
!$OMP end do nowait
      end do
!$OMP end parallel

      call exb_NL_term( wf, psi, chi, wbf )

!$OMP parallel workshare
      wbf(:,:,:,:,:) = wbf(:,:,:,:,:) - wef(:,:,:,:,:) ! Magnetic part
!$OMP end parallel workshare

!$OMP parallel
      do im = 0, nm
!$OMP do
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                wef(mx,my,iz,iv,im) = - fcs(ranks) / Znum(ranks) * wef(mx,my,iz,iv,im)  &
                      * tau(ranks) * conjg( wf(mx,my,iz,iv,im) ) / fmx(iz,iv,im)
              end do
            end do
          end do
        end do
!$OMP end do nowait
      end do
!$OMP end parallel

      call intgrl_v0_moment ( wef, wc3 )

      call intgrl_thet ( wc3, wc2 )

!$OMP parallel do
      do my = ist_y, iend_y
        do mx = -nx, nx
          neint(mx,my) = real( wc2(mx,my), kind=DP )
        end do
      end do

!$OMP parallel
      do im = 0, nm
!$OMP do
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                wbf(mx,my,iz,iv,im) = - fcs(ranks) / Znum(ranks) * wbf(mx,my,iz,iv,im)  &
                      * tau(ranks) * conjg( wf(mx,my,iz,iv,im) ) / fmx(iz,iv,im)
              end do
            end do
          end do
        end do
!$OMP end do nowait
      end do
!$OMP end parallel

      call intgrl_v0_moment ( wbf, wc3 )

      call intgrl_thet ( wc3, wc2 )

!$OMP parallel do
      do my = ist_y, iend_y
        do mx = -nx, nx
          nmint(mx,my) = real( wc2(mx,my), kind=DP )
        end do
      end do


      deallocate( wf )
      deallocate( wef )
      deallocate( wbf )
      deallocate( psi )
      deallocate( chi )
      deallocate( wc3 )


  END SUBROUTINE trans_sum


!--------------------------------------
  SUBROUTINE trans_triad ( time, ff, phi, Al )
!--------------------------------------
!   Triad transfer    T [(kx,ky)|(px,py),(qx,qy)] * delta(kx+px+qx=0,ky+py+qy=0)
!   Setting (kx,ky)=(diagx,diagy) and (qx=-px-kx,qy=-py-ky), 
!   T=T(px,py) represents transfer from (px,py) and (-px-diagx,-py-diagy) to (diagx,diagy).

    real(kind=DP), intent(in)                                     :: time
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1)                             :: phi, Al

    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: gg
    complex(kind=DP), dimension(:,:,:,:), allocatable :: psi, chi, wg
    complex(kind=DP), dimension(:,:,:), allocatable :: wps, wch
    real(kind=DP), dimension(:,:), allocatable :: jkpq_es, jpqk_es, jqkp_es, &
                                                  jkpq_em, jpqk_em, jqkp_em
    real(kind=DP) :: ceff
    integer  ::  mx, my, iz, iv, im, it, mxt, myt
    character(1)   :: srank
    character(3)   :: cnew
    character(4)   :: cmx, cmy
    character(256)   :: env_string
    integer, save ::  iflg

    data iflg / 0 /
    namelist /triad/ mxt, myt


      !%%% Initialize triad_diag_mxt, triad_diag_myt, nbuff %%%
      if( iflg == 0 ) then

        iflg = 1

        close(inml)
        call getenv ( 'fu05',env_string )
        open(inml, file=env_string )

        allocate(triad_diag_mxt(0:num_triad_diag-1))
        allocate(triad_diag_myt(0:num_triad_diag-1))

        do it = 0, num_triad_diag-1
          read(inml, nml=triad)
          triad_diag_mxt(it) = mxt
          triad_diag_myt(it) = myt

          !fj start 202011
          !if ( rank == 0 ) then
          !  write( srank, fmt="(i1.1)" ) ranks
          !  write( cnew,  fmt="(i3.3)" ) inum
          !  write( cmx,   fmt="(i4.4)" ) triad_diag_mxt(it)
          !  write( cmy,   fmt="(i4.4)" ) triad_diag_myt(it)
          !  open( otri, file=trim(f_phi)//"s"//srank//"mx"//cmx//"my"//cmy//".tri."//cnew, & 
          !              form="unformatted", status="replace" )
          !  close( otri )
          !end if
          write( srank, fmt="(i1.1)" ) ranks
          write( cnew,  fmt="(i3.3)" ) inum
          write( cmx,   fmt="(i4.4)" ) triad_diag_mxt(it)
          write( cmy,   fmt="(i4.4)" ) triad_diag_myt(it)
          call fileio_open_tri( trim(f_phi), cmx, cmy, .true. )
          call fileio_close_tri
          !fj end 202011

        end do

        if ( mod(2*nz*(nm+1),nprocw) == 0 ) then
          nbuff = 2*nz*(nm+1)/nprocw
        else
          nbuff = 2*nz*(nm+1)/nprocw + 1
        end if

      end if

      !%%% Set gg, psi, chi %%%
      allocate( gg(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) )
      allocate( psi(-nx:nx,0:ny,-nz:nz-1,0:nm) )
      allocate( chi(-nx:nx,0:ny,-nz:nz-1,0:nm) )
      allocate( wg(-nx:nx,0:global_ny,1:2*nv,0:nbuff-1) )
      allocate( wps(-nx:nx,0:global_ny,0:nbuff-1) )
      allocate( wch(-nx:nx,0:global_ny,0:nbuff-1) )
      allocate( jkpq_es(-nx:nx,-global_ny:global_ny) )
      allocate( jpqk_es(-nx:nx,-global_ny:global_ny) )
      allocate( jqkp_es(-nx:nx,-global_ny:global_ny) )
      allocate( jkpq_em(-nx:nx,-global_ny:global_ny) )
      allocate( jpqk_em(-nx:nx,-global_ny:global_ny) )
      allocate( jqkp_em(-nx:nx,-global_ny:global_ny) )
!$OMP parallel workshare
      gg(:,:,:,:,:) = (0._DP, 0._DP)
      psi(:,:,:,:) = (0._DP, 0._DP)
      chi(:,:,:,:) = (0._DP, 0._DP)
!$OMP end parallel workshare

!$OMP parallel do collapse(2)
      do im = 0, nm
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                psi(mx,my,iz,im) = j0(mx,my,iz,im) * phi(mx,my,iz)
                chi(mx,my,iz,im) = j0(mx,my,iz,im) * Al(mx,my,iz)
              end do
            end do
          end do
      end do
!$OMP parallel do collapse(2)
      do im = 0, nm
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                gg(mx,my,iz,iv,im) = ff(mx,my,iz,iv,im)                       &
                                   + sgn(ranks) * Znum(ranks) * fmx(iz,iv,im) &
                                              / tau(ranks) * psi(mx,my,iz,im)
              end do
            end do
          end do
        end do
      end do

      !%%% Multiply the jacobian on gg for zz,vl,mu integration %%%
      if ( rankm == 0 ) then
        im = 0
!$OMP parallel do private(ceff)
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              ceff = 0._DP
              do my = ist_y, iend_y
                do mx = -nx, nx
                  gg(mx,my,iz,iv,im) = sqrt(ceff) * gg(mx,my,iz,iv,im)
                end do
              end do
            end do
          end do
        im = 1
!$OMP parallel do private(ceff)
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              ceff = vp(iz,im) * dvp(iz) * dv * twopi * (     &
                     1._DP + 1._DP/12._DP + 22._DP / 720._DP )&
                   * rootg(iz) / cfsrf                        &
                   * tau(ranks) / (2._DP * fmx(iz,iv,im))
              do my = ist_y, iend_y
                do mx = -nx, nx
                  gg(mx,my,iz,iv,im) = sqrt(ceff) * gg(mx,my,iz,iv,im)
                end do
              end do
            end do
          end do
        im = 2
!$OMP parallel do private(ceff)
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              ceff = vp(iz,im) * dvp(iz) * dv * twopi * ( &
                     1._DP - 11._DP / 720._DP )           &
                   * rootg(iz) / cfsrf                    &
                   * tau(ranks) / (2._DP * fmx(iz,iv,im))
              do my = ist_y, iend_y
                do mx = -nx, nx
                  gg(mx,my,iz,iv,im) = sqrt(ceff) * gg(mx,my,iz,iv,im)
                end do
              end do
            end do
          end do
!$OMP parallel do private(ceff)
        do im = 3, nm
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              ceff = vp(iz,im) * dvp(iz) * dv * twopi     &
                   * rootg(iz) / cfsrf                    &
                   * tau(ranks) / (2._DP * fmx(iz,iv,im))
              do my = ist_y, iend_y
                do mx = -nx, nx
                  gg(mx,my,iz,iv,im) = sqrt(ceff) * gg(mx,my,iz,iv,im)
                end do
              end do
            end do
          end do
        end do
      else
!$OMP parallel do private(ceff)
        do im = 0, nm
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              ceff = vp(iz,im) * dvp(iz) * dv * twopi     &
                   * rootg(iz) / cfsrf                    &
                   * tau(ranks) / (2._DP * fmx(iz,iv,im))
              do my = ist_y, iend_y
                do mx = -nx, nx
                  gg(mx,my,iz,iv,im) = sqrt(ceff) * gg(mx,my,iz,iv,im)
                end do
              end do
            end do
          end do
        end do
      end if

      !%%% Transpose from (y,z,v,m) to (zm,z,v,m) decomposition %%%
      call trans_triad_transpose(gg, psi, chi, wg, wps, wch)

      do it = 0, num_triad_diag-1

        !%%% Calculate traid coupling %%%
        call trans_triad_coupling(it, wg, wps, wch, jkpq_es, jpqk_es, jqkp_es,&
                                                    jkpq_em, jpqk_em, jqkp_em)

        !%%% Output %%%
        !fj start 202011
        !if ( rank == 0 ) then
        !  write( srank, fmt="(i1.1)" ) ranks
        !  write( cnew,  fmt="(i3.3)" ) inum
        !  write( cmx,   fmt="(i4.4)" ) triad_diag_mxt(it)
        !  write( cmy,   fmt="(i4.4)" ) triad_diag_myt(it)
        !  open( otri, file=trim(f_phi)//"s"//srank//"mx"//cmx//"my"//cmy//".tri."//cnew, & 
        !              form="unformatted", status="unknown", position="append" )
        !  write( unit=otri ) time, jkpq_es, jpqk_es, jqkp_es, jkpq_em, jpqk_em, jqkp_em
        !  close( otri )
        !end if
        write( srank, fmt="(i1.1)" ) ranks
        write( cnew,  fmt="(i3.3)" ) inum
        write( cmx,   fmt="(i4.4)" ) triad_diag_mxt(it)
        write( cmy,   fmt="(i4.4)" ) triad_diag_myt(it)
        call fileio_open_tri( trim(f_phi), cmx, cmy, .false. )
        call fileio_write_tri( jkpq_es, jpqk_es, jqkp_es, jkpq_em, jpqk_em, jqkp_em, time )
        call fileio_close_tri
        !fj end 202011

      end do

      deallocate( gg )
      deallocate( psi )
      deallocate( chi )
      deallocate( wg )
      deallocate( wps )
      deallocate( wch )
      deallocate( jkpq_es )
      deallocate( jpqk_es )
      deallocate( jqkp_es )
      deallocate( jkpq_em )
      deallocate( jpqk_em )
      deallocate( jqkp_em )


  END SUBROUTINE trans_triad


!--------------------------------------
  SUBROUTINE trans_triad_transpose ( gg, psi, chi, wg, wps, wch )
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: gg
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,0:nm)        :: psi, chi
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:global_ny,1:2*nv,0:nbuff-1) :: wg
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:global_ny,0:nbuff-1)        :: wps, wch

    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: send_gg, recv_gg
    complex(kind=DP), dimension(:,:,:,:), allocatable :: send_psi, recv_psi, &
                                                         send_chi, recv_chi
    integer  ::  mx, my, iz, iv, im, izm, ibuff, iprocw, global_my


      allocate( send_gg(-nx:nx,0:ny,1:2*nv,0:nbuff-1,0:nprocw-1) )
      allocate( recv_gg(-nx:nx,0:ny,1:2*nv,0:nbuff-1,0:nprocw-1) )
      allocate( send_psi(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1) )
      allocate( recv_psi(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1) )
      allocate( send_chi(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1) )
      allocate( recv_chi(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1) )

    !%%% PACK: gg -> send_gg %%%
!$OMP parallel workshare
      send_gg(:,:,:,:,:) = (0._DP, 0._DP)
      send_psi(:,:,:,:) = (0._DP, 0._DP)
      send_chi(:,:,:,:) = (0._DP, 0._DP)
!$OMP end parallel workshare
!$OMP parallel do private(izm,ibuff,iprocw)
      do im = 0, nm
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                izm = (2*nz)*im + (iz + nz)
                ibuff = mod(izm, nbuff)
                iprocw = izm / nbuff
                send_gg(mx,my,iv,ibuff,iprocw) = gg(mx,my,iz,iv,im)
              end do
            end do
          end do
        end do
      end do
!$OMP parallel do private(izm,ibuff,iprocw)
      do im = 0, nm
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                izm = (2*nz)*im + (iz + nz)
                ibuff = mod(izm, nbuff)
                iprocw = izm / nbuff
                send_psi(mx,my,ibuff,iprocw) = psi(mx,my,iz,im)
                send_chi(mx,my,ibuff,iprocw) = chi(mx,my,iz,im)
              end do
            end do
          end do
      end do

    !%%% TRANSPOSE: send_gg -> recv_gg %%%
      call MPI_Alltoall( send_gg,                      &
                         (2*nx+1)*(ny+1)*(2*nv)*nbuff, &
                         MPI_DOUBLE_COMPLEX,           &
                         recv_gg,                      &
                         (2*nx+1)*(ny+1)*(2*nv)*nbuff, &
                         MPI_DOUBLE_COMPLEX,           &
                         fft_comm_world,               &
                         ierr_mpi )
      call MPI_Alltoall( send_psi,              &
                         (2*nx+1)*(ny+1)*nbuff, &
                         MPI_DOUBLE_COMPLEX,    &
                         recv_psi,              &
                         (2*nx+1)*(ny+1)*nbuff, &
                         MPI_DOUBLE_COMPLEX,    &
                         fft_comm_world,        &
                         ierr_mpi )
      call MPI_Alltoall( send_chi,              &
                         (2*nx+1)*(ny+1)*nbuff, &
                         MPI_DOUBLE_COMPLEX,    &
                         recv_chi,              &
                         (2*nx+1)*(ny+1)*nbuff, &
                         MPI_DOUBLE_COMPLEX,    &
                         fft_comm_world,        &
                         ierr_mpi )

    !%%% UNPACK: recv_gg -> wg %%%
!$OMP parallel do private(iprocw,my)
      do ibuff = 0, nbuff-1
        do iv = 1, 2*nv
          do global_my = 0, global_ny
            do mx = -nx, nx
              iprocw = global_my / (ny+1)
              my = mod(global_my, ny+1)
              wg(mx,global_my,iv,ibuff) = recv_gg(mx,my,iv,ibuff,iprocw)
            end do
          end do
        end do
      end do
!$OMP parallel do private(iprocw,my)
      do ibuff = 0, nbuff-1
          do global_my = 0, global_ny
            do mx = -nx, nx
              iprocw = global_my / (ny+1)
              my = mod(global_my, ny+1)
              wps(mx,global_my,ibuff) = recv_psi(mx,my,ibuff,iprocw)
              wch(mx,global_my,ibuff) = recv_chi(mx,my,ibuff,iprocw)
            end do
          end do
      end do

      deallocate( send_gg )
      deallocate( recv_gg )
      deallocate( send_psi )
      deallocate( recv_psi )
      deallocate( send_chi )
      deallocate( recv_chi )
                                         !%%% For debug %%%
                                         !iz = nz-1
                                         !im = nm
                                         !call MPI_Allgather(     &
                                         !     psi(-nx,0,iz,im),  &
                                         !     (2*nx+1)*(ny+1),   &
                                         !     MPI_DOUBLE_COMPLEX,&
                                         !     wch(-nx,0,0),      &
                                         !     (2*nx+1)*(ny+1),   &
                                         !     MPI_DOUBLE_COMPLEX,&
                                         !     fft_comm_world,    &
                                         !     ierr_mpi)
                                         !izm = (2*nz)*im+(iz+nz)
                                         !ibuff = mod(izm,nbuff)
                                         !iprocw = izm/nbuff
                                         !if (rankg==iprocw) then
                                         !write(888,*)"#",ibuff,iprocw
                                         !do my = 0, global_ny
                                         !  do mx = -nx, nx
                                         !    write(888,*) mx, my,    &
                                         !     dble(wps(mx,my,ibuff)),&
                                         !    aimag(wps(mx,my,ibuff)),&
                                         !     dble(wch(mx,my,0)),    &
                                         !    aimag(wch(mx,my,0))
                                         !  end do
                                         !  write(888,*)
                                         !end do
                                         !end if
                                         !call MPI_Finalize(ierr_mpi)
                                         !stop
                                         !%%%%%%%%%%%%%%%%%

  END SUBROUTINE trans_triad_transpose


!--------------------------------------
  SUBROUTINE trans_triad_coupling ( it, wg, wps, wch, jkpq_es, jpqk_es, jqkp_es,&
                                                      jkpq_em, jpqk_em, jqkp_em )
!--------------------------------------

    integer, intent(in)                              :: it
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:global_ny,1:2*nv,0:nbuff-1) :: wg
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:global_ny,0:nbuff-1)        :: wps, wch
    real(kind=DP), intent(out), &
      dimension(-nx:nx,-global_ny:global_ny) :: jkpq_es, jpqk_es, jqkp_es, &
                                                jkpq_em, jpqk_em, jqkp_em

    complex(kind=DP), dimension(:,:), allocatable :: gg, psi, chi
    real(kind=DP), dimension(:,:), allocatable :: &
                        wkpq_es, wpqk_es, wqkp_es, wkpq_em, wpqk_em, wqkp_em
    real(kind=DP) :: dky, wky(-global_ny:global_ny), cs1
    integer :: mx, my, px, py, qx, qy ! (mx,my) + (px,py) + (qx,qy) = (0,0)
    integer :: iv, ibuff

      allocate(  gg(-nx:nx,-global_ny:global_ny) )
      allocate( psi(-nx:nx,-global_ny:global_ny) )
      allocate( chi(-nx:nx,-global_ny:global_ny) )
      allocate( wkpq_es(-nx:nx,-global_ny:global_ny) )
      allocate( wpqk_es(-nx:nx,-global_ny:global_ny) )
      allocate( wqkp_es(-nx:nx,-global_ny:global_ny) )
      allocate( wkpq_em(-nx:nx,-global_ny:global_ny) )
      allocate( wpqk_em(-nx:nx,-global_ny:global_ny) )
      allocate( wqkp_em(-nx:nx,-global_ny:global_ny) )

      !-set (mx,my)-
      mx = triad_diag_mxt(it)
      my = triad_diag_myt(it)

      !-set wky-
      dky = ky(1) - ky(0)
      do py = -global_ny, global_ny
        wky(py) = dky * real( py, kind=DP )
      end do

!$OMP parallel workshare
      wkpq_es(:,:) = 0._DP
      wpqk_es(:,:) = 0._DP
      wqkp_es(:,:) = 0._DP
      wkpq_em(:,:) = 0._DP
      wpqk_em(:,:) = 0._DP
      wqkp_em(:,:) = 0._DP
!$OMP end parallel workshare
      cs1 = sqrt( tau(ranks) / Anum(ranks) )
!$OMP parallel default(none) &
!$OMP shared(mx,my,kx,wky,vl,cs1,nbuff,wps,wch,wg) &
!$OMP shared(wkpq_es,wpqk_es,wqkp_es,wkpq_em,wpqk_em,wqkp_em) &
!$OMP private(px,py,qx,qy,iv,ibuff,psi,chi,gg)
!$OMP do reduction(+:wkpq_es) reduction(+:wpqk_es) reduction(+:wqkp_es) &
!$OMP    reduction(+:wkpq_em) reduction(+:wpqk_em) reduction(+:wqkp_em)
      do ibuff = 0, nbuff-1

        !-copy-
        do py = 0, global_ny
          do px = -nx, nx
            psi(px,py) = wps(px,py,ibuff)
            chi(px,py) = wch(px,py,ibuff)
          end do
        end do
        do py = 1, global_ny
          do px = -nx, nx
            psi(-px,-py) = conjg( wps(px,py,ibuff) )
            chi(-px,-py) = conjg( wch(px,py,ibuff) )
          end do
        end do

        do iv = 1, 2*nv

          !-copy-
          do py = 0, global_ny
            do px = -nx, nx
              gg(px,py) = wg(px,py,iv,ibuff)
            end do
          end do
          do py = 1, global_ny
            do px = -nx, nx
              gg(-px,-py) = conjg( wg(px,py,iv,ibuff) )
            end do
          end do

          !-triad coupling among (mx,my)+(px,py)+(qx,qy)=0-
          do py = max(-global_ny-my,-global_ny), min(global_ny,global_ny-my)
            qy = - my - py
            do px = max(-nx-mx,-nx), min(nx,nx-mx)
              qx = - mx - px
             !wkpq(px,py) = wkpq(px,py)                                         &
             !  - (- kx(px) * wky(qy) + wky(py) * kx(qx))                       &
             !  * real((  (psi(px,py) - cs1 * vl(iv) * chi(px,py)) * gg(qx,qy)  &
             !          - (psi(qx,qy) - cs1 * vl(iv) * chi(qx,qy)) * gg(px,py)) &
             !                                              * gg(mx,my), kind=DP)
             !wpqk(px,py) = wpqk(px,py)                                         &
             !  - (- kx(qx) * wky(my) + wky(qy) * kx(mx))                       &
             !  * real((  (psi(qx,qy) - cs1 * vl(iv) * chi(qx,qy)) * gg(mx,my)  &
             !          - (psi(mx,my) - cs1 * vl(iv) * chi(mx,my)) * gg(qx,qy)) &
             !                                              * gg(px,py), kind=DP)
             !wqkp(px,py) = wqkp(px,py)                                         &
             !  - (- kx(mx) * wky(py) + wky(my) * kx(px))                       &
             !  * real((  (psi(mx,my) - cs1 * vl(iv) * chi(mx,my)) * gg(px,py)  &
             !          - (psi(px,py) - cs1 * vl(iv) * chi(px,py)) * gg(mx,my)) &
             !                                              * gg(qx,qy), kind=DP)
              wkpq_es(px,py) = wkpq_es(px,py)             &
                - (- kx(px) * wky(qy) + wky(py) * kx(qx)) &
                * real((  (psi(px,py)) * gg(qx,qy)        &
                        - (psi(qx,qy)) * gg(px,py))       &
                                      * gg(mx,my), kind=DP)
              wpqk_es(px,py) = wpqk_es(px,py)             &
                - (- kx(qx) * wky(my) + wky(qy) * kx(mx)) &
                * real((  (psi(qx,qy)) * gg(mx,my)        &
                        - (psi(mx,my)) * gg(qx,qy))       &
                                      * gg(px,py), kind=DP)
              wqkp_es(px,py) = wqkp_es(px,py)             &
                - (- kx(mx) * wky(py) + wky(my) * kx(px)) &
                * real((  (psi(mx,my)) * gg(px,py)        &
                        - (psi(px,py)) * gg(mx,my))       &
                                      * gg(qx,qy), kind=DP)
              wkpq_em(px,py) = wkpq_em(px,py)                        &
                - (- kx(px) * wky(qy) + wky(py) * kx(qx))            &
                * real((  (- cs1 * vl(iv) * chi(px,py)) * gg(qx,qy)  &
                        - (- cs1 * vl(iv) * chi(qx,qy)) * gg(px,py)) &
                                                 * gg(mx,my), kind=DP)
              wpqk_em(px,py) = wpqk_em(px,py)                        &
                - (- kx(qx) * wky(my) + wky(qy) * kx(mx))            &
                * real((  (- cs1 * vl(iv) * chi(qx,qy)) * gg(mx,my)  &
                        - (- cs1 * vl(iv) * chi(mx,my)) * gg(qx,qy)) &
                                                 * gg(px,py), kind=DP)
              wqkp_em(px,py) = wqkp_em(px,py)                        &
                - (- kx(mx) * wky(py) + wky(my) * kx(px))            &
                * real((  (- cs1 * vl(iv) * chi(mx,my)) * gg(px,py)  &
                        - (- cs1 * vl(iv) * chi(px,py)) * gg(mx,my)) &
                                                 * gg(qx,qy), kind=DP)
            end do
          end do

        end do
      end do
!$OMP end do
!$OMP end parallel

      !-zz,vl,mu integration-
!$OMP parallel workshare
      jkpq_es(:,:) = 0._DP
      jpqk_es(:,:) = 0._DP
      jqkp_es(:,:) = 0._DP
      jkpq_em(:,:) = 0._DP
      jpqk_em(:,:) = 0._DP
      jqkp_em(:,:) = 0._DP
!$OMP end parallel workshare
      call MPI_Allreduce( wkpq_es, jkpq_es, (2*nx+1)*(2*global_ny+1), &
                          MPI_DOUBLE_PRECISION,                       &
                          MPI_SUM, sub_comm_world, ierr_mpi )
      call MPI_Allreduce( wpqk_es, jpqk_es, (2*nx+1)*(2*global_ny+1), &
                          MPI_DOUBLE_PRECISION,                       &
                          MPI_SUM, sub_comm_world, ierr_mpi )
      call MPI_Allreduce( wqkp_es, jqkp_es, (2*nx+1)*(2*global_ny+1), &
                          MPI_DOUBLE_PRECISION,                       &
                          MPI_SUM, sub_comm_world, ierr_mpi )
      call MPI_Allreduce( wkpq_em, jkpq_em, (2*nx+1)*(2*global_ny+1), &
                          MPI_DOUBLE_PRECISION,                       &
                          MPI_SUM, sub_comm_world, ierr_mpi )
      call MPI_Allreduce( wpqk_em, jpqk_em, (2*nx+1)*(2*global_ny+1), &
                          MPI_DOUBLE_PRECISION,                       &
                          MPI_SUM, sub_comm_world, ierr_mpi )
      call MPI_Allreduce( wqkp_em, jqkp_em, (2*nx+1)*(2*global_ny+1), &
                          MPI_DOUBLE_PRECISION,                       &
                          MPI_SUM, sub_comm_world, ierr_mpi )

      deallocate(  gg )
      deallocate( psi )
      deallocate( chi )
      deallocate( wkpq_es )
      deallocate( wpqk_es )
      deallocate( wqkp_es )
      deallocate( wkpq_em )
      deallocate( wpqk_em )
      deallocate( wqkp_em )

  END SUBROUTINE trans_triad_coupling


END MODULE GKV_trans
