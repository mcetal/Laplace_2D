!
! Program that solves 2D Laplace's equation using FIEM
!
!     Domain can be multiply connected, bounded or unbounded in extent
!     Dirichlet or Neumann BCs
!

program LAPLACE_2D

   implicit none
!
! Maximum dimensions
   integer, parameter :: kmax = 10, npmax = 2048, nmax = kmax*npmax
!
! System size
   integer k0, k, nd, nbk
!
! Geometry of Domain
   real(kind=8) :: ak(kmax), bk(kmax)
   integer :: ncyc(kmax)
   complex(kind=16) :: zk(kmax)
!
! Grid Geometry
   integer :: nx, ny
!
! Discretization
   real(kind=8) :: x(nmax), y(nmax), ds_dth(nmax), kappa(nmax)
   complex(kind=16) :: z(nmax), dz(nmax)
!
! Problem type
   logical :: debug, bounded, dirichlet
   
      call INIT_DOMAIN(kmax, k0, k, nd, nbk, nx, ny, ak, bk, zk, ncyc)  
      call BUILD_DOMAIN(k0, k, nd, nbk, ak, bk, zk, ncyc, z, dz, &
                        ds_dth, kappa)
      print *, zk(1), zk(2)
   
end program LAPLACE_2D

!----------------------------------------------------------------------

subroutine INIT_DOMAIN(kmax, k0, k, nd, nbk, nx, ny, ak, bk, zk, ncyc)
!
! Initialize size of system and parameters for domain
! Input:
!    kmax:  Max number of holes
! Returns:
!    k0:  0 if bounded domain, 1 if unbounded domain
!    k:   number of interior contours if bounded, number of contours if 
!         unbounded
!    nd:  number of points per hole
!    nbk: total number of discretization points
!    nx, ny: number of grid points in x and y direction
!    ak, bk: major and minor axes of elliptical holes
!    ncyc: number of starfish arms
!
   implicit none
   integer, intent(in) :: kmax
   integer, intent(out) :: k0, k, nd, nbk, nx, ny
   real(kind=8), intent(out) :: ak(kmax), bk(kmax)
   integer, intent(out) :: ncyc(kmax)
   complex(kind=16), intent(out) :: zk(kmax)
   
      k0 = 0
      k = 0
      if (k + 1 - k0 > kmax) then
         print *, 'Too Many Holes!'
         print *, 'k = ', k
         print *, 'kmax = ', kmax
         stop
      end if
      nd = 64
      nbk = k*nd + (1-k0)*nd
      ak = (/ 1.d0 /)
      bk = (/ 1.d0 /)
      ncyc = (/ 0 /)
      zk = (/ (0.d0,0.d0), (1.d0,1.d0) /)
      
end subroutine INIT_DOMAIN

!----------------------------------------------------------------------

subroutine BUILD_DOMAIN(k0, k, nd, nbk, ak, bk, zk, ncyc, z, dz, &
                        ds_dth, kappa)
!
! Builds domain
