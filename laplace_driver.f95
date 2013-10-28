!
! Program that solves 2D Laplace's equation using FIEM
!
!     Domain can be multiply connected, bounded or unbounded in extent
!     Dirichlet or Neumann BCs
!

program LAPLACE_2D
!
! testing module
   use geometry_mod

   implicit none
!
! Integral equation density and log sources
   real(kind=8) :: mu(nmax), A_log(kmax)
!
! Right hand side to integral equation 
   real(kind=8) :: rhs(nmax+kmax)
!
!  Matrix equation variables for GMRES
!  MAXL is the maximum nubmer of GMRES iterations performed
!       before restarting.
!  LRWORK is the dimension of a real workspace needed by DGMRES.
!  LIWORK is the dimension of an integer workspace needed by DGMRES.
!  GMWORK and IWORK are work arrays used by DGMRES
!
   integer, parameter :: maxl = 50, liwork=30,  & 
                         lrwork=10+(nmax+kmax)*(maxl+6)+maxl*(maxl+3)
   real(kind=8) :: gmwork(lrwork), soln(nmax+kmax)
   integer :: igwork(liwork)
!
! Problem type
   logical :: debug, bounded, dirichlet
   
!
! Initialize print output units
   call PRINI(6, 13)

!
! Initialize Geometry, problem type, and system size
   call INITIALIZE(debug, bounded, dirichlet)
   call INIT_HOLE_GEO() 
   call BUILD_DOMAIN()
   
!
! Set boundary conditions
   call GET_BCS(debug, bounded, rhs)   
   call PRIN2(' rhs = *', rhs, nbk)
!
! Solve integral equation
   call SOLVE (maxl, rhs, lrwork, liwork, soln, mu, A_log, &
               gmwork, igwork)

end program LAPLACE_2D

!----------------------------------------------------------------------

subroutine INITIALIZE(debug, bounded, dirichlet)
!
! reads in k0, k, nd from module and calculates nbk
! describes problem type
! Returns:
!   debug :: true if in debugging mode
!   bounded :: true if bounded domain, false if unbounded
!   dirichlet :: true if Dirichlet BVP, false if Neumann
   use geometry_mod, only: pi, eye, kmax, npmax, nbk, k0, k, nd, h
   implicit none
   logical, intent(out) :: debug, bounded, dirichlet

!
! initialize constants
      pi = 4.d0*datan(1.d0)
      eye = dcmplx(0.d0,1.d0)

!
! initialize number of holes and points per hole
      k0 = 0
      k = 0
      nd = 32
      bounded = k0==0
      print *, 'bounded = ', bounded
!
! Check dimensions are ok
      if (k + 1 - k0 > kmax) then
         print *, 'Too Many Holes!'
         print *, 'k = ', k
         print *, 'kmax = ', kmax
         stop
      elseif (nd > npmax) then
         print *, 'Too Many Points Per Hole!'
         print *, 'nd = ', nd
         print *, 'npmax = ', npmax
         stop
      end if

      nbk = k*nd + (1-k0)*nd
      print *, 'nbk = ', nbk

!
! initialize mesh spacing

      h = 2.d0*pi/nd
      
!
! initialize problem type
      debug = .true.
      dirichlet = .true.

end subroutine INITIALIZE

!----------------------------------------------------------------------

subroutine INIT_HOLE_GEO()
!
! Initialize size of system and parameters for domain
! Returns:
!    ak, bk: major and minor axes of elliptical holes
!    zk:  geometric centre of holes
!    ncyc: number of starfish arms
!
   use geometry_mod, only: k0, k, ak, bk, ncyc, zk
   implicit none
   
      ak(1) = 1.d0
      bk(1) = 1.d0
      ncyc(1) = 0
      zk(1) = dcmplx(0.d0,0.d0)
      
      ak(2) = 0.5d0
      bk(2) = 0.5d0
      ncyc(2) = 0
      zk(2) = dcmplx(0.d0,0.d0)
      
end subroutine INIT_HOLE_GEO

!----------------------------------------------------------------------

subroutine GET_BCS(debug, bounded, rhs)

! Constructs right hand side of integral equation
! Inputs:
!   debug: logical for debugging mode
!   bounded: logical for bounded domain 
! Returns:
!   rhs: right hand size of integral equation

   use geometry_mod, only: k0, k, nd, nbk, zk, z
   implicit none
   logical, intent(in) :: debug, bounded
   real(kind=8), intent(out) :: rhs(nbk)
   real(kind=8) :: U_EXACT
   integer :: i
   
      do i = 1, nbk
	 rhs(i) = U_EXACT(bounded, z(i))
      end do

end subroutine GET_BCS


!----------------------------------------------------------------------

subroutine POLAR_COORD(z, theta, r)      

!
! Given a point z = x + i y, calculate polar coords r and theta
! Inputs: 
!   z: point in complex plane 
! Returns
!   r, theta: polar coordinates of z

   use geometry_mod, only : pi
   implicit none
   complex(kind=8), intent(in) :: z
   real(kind=8), intent(out) :: theta, r

     r = cdabs(z)
     theta = dacos(dreal(z)/r)
     if (dimag(z).lt.0.d0) theta = 2.d0*pi - theta
     
end subroutine POLAR_COORD

!----------------------------------------------------------------------

real(kind=8) function U_EXACT(bounded, z)     

! This is an exact test function. It should be harmonic and well 
! defined in the domain
! It's used to generate boundary conditions and to test accuracy of 
! solutions
! Input:
!   bounded: logical indicated whether domain is bounded or not
!   z:  point in space

   use geometry_mod, only: k0, k, zk
   implicit none
   logical, intent(in) :: bounded
   complex(kind=8), intent(in) :: z 
   complex(kind=8) :: zsrc, zdis
   integer :: n
   real(kind=8) :: theta, r, A, B
   
!
! harmonic mode
      n = 2
      A = 1.d0
      B = 0.d0
!
! if bounded,
!   u = A r^n cos(n theta) + B r^n sin(n theta)
! if unbounded
!   u = A cos(n theta)/r^n + B sin(n theta)/r^n

      zdis = z-zk(1)
      call POLAR_COORD(zdis, theta, r)

      if (bounded) then 
         U_EXACT = (A*dcos(n*theta) + B*dsin(n*theta)) * r**n
       else
         U_EXACT = (A*dcos(n*theta) + B*dsin(n*theta)) / r**n
      end if

end function U_EXACT
