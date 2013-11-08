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
   logical :: debug, dirichlet
!
! Target points
   integer, parameter :: ntar = 20
   real(kind=8) :: u_tar(ntar)
   complex(kind=8) :: z_tar(ntar)
   
!
! Initialize print output units
   call PRINI(6, 13)

!
! Initialize Geometry, problem type, and system size
   call INITIALIZE(debug, dirichlet)
   call INIT_HOLE_GEO() 
   call BUILD_DOMAIN()
   
!
! Get target points
   if (debug) call GET_TARGETS(ntar, z_tar)
   
!
! Set boundary conditions
   call GET_BCS(debug, dirichlet, rhs)   
   call PRIN2(' rhs = *', rhs, nbk)
   
!
! Solve integral equation
   call SOLVE (maxl, rhs, lrwork, liwork, dirichlet, soln, mu, A_log, &
               gmwork, igwork)
               
!
! Check solution at target points
   if (debug) call GET_SOL_TAR(ntar, z_tar, mu, A_log, u_tar)

end program LAPLACE_2D

!----------------------------------------------------------------------

subroutine INITIALIZE(debug, dirichlet)
!
! reads in k0, k, nd from module and calculates nbk
! describes problem type
! Returns:
!   debug :: true if in debugging mode
!   bounded :: true if bounded domain, false if unbounded
!   dirichlet :: true if Dirichlet BVP, false if Neumann
   use geometry_mod, only: pi, eye, kmax, npmax, nbk, k0, k, nd, h, &
                           bounded
   implicit none
   logical, intent(out) :: debug, dirichlet

!
! initialize constants
      pi = 4.d0*datan(1.d0)
      eye = dcmplx(0.d0,1.d0)

!
! initialize number of holes and points per hole
      k0 = 1
      k = 1
      nd = 64
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
      bk(1) = 0.2d0
      ncyc(1) = 3
      zk(1) = dcmplx(0.2d0,-2.d0)
      
      ak(2) = 0.5d0
      bk(2) = 0.5d0
      ncyc(2) = 0
      zk(2) = dcmplx(0.d0,0.d0)
      
end subroutine INIT_HOLE_GEO

   
!----------------------------------------------------------------------

subroutine GET_TARGETS(ntar, z_tar)
!
! Gets a small number of target points in the domain for checking accuracy
! Inputs:
!   ntar: number of target points
! Returns:
!   z_tar: location of target points in complex plain
!
   use geometry_mod, only: pi, k0, k, zk, ak, bk, Z_PLOT
   implicit none
   integer, intent(in) :: ntar
   complex(kind=8), intent(out) :: z_tar(ntar)
! local
   integer :: i
   real(kind=8) :: dth, theta, a_tar, b_tar
   character(32) :: options
   
      dth = 2.d0*pi/ntar
      
      do i = 1, ntar
         theta = dth*(i-1.d0)
         if ( (k0==0) .and. (k==0) ) then
            a_tar = 0.5d0*ak(1)
            b_tar = 0.5d0*bk(1)
         elseif ( (k0==0) .and. (k==1) ) then
            a_tar = 0.5d0*(ak(1) + ak(2))
            b_tar = 0.5d0*(bk(1) + bk(2))
         elseif ( (k0==1) .and. (k==1) ) then
            a_tar = 2.d0*max(ak(1), bk(1))
            b_tar = a_tar
         else
            print *, 'Cannot calculate target points for this geometry'
            print *, 'Errors in solution check may result'
         end if
         z_tar(i) = zk(1) + dcmplx(a_tar*dcos(theta), b_tar*dsin(theta))
      end do
      
      call PRIN2(' z_tar = *', z_tar, 2*ntar)
      
! dump out for plotting
      open (unit = 21, file = 'mat_plots/target_points.m', &
            status = 'unknown')
      options = '''r*'''
      call Z_PLOT( z_tar, ntar, options, 21)
      close(21)
      
end subroutine GET_TARGETS
!----------------------------------------------------------------------

subroutine GET_BCS(debug, dirichlet, rhs)

! Constructs right hand side of integral equation
! Inputs:
!   debug: logical for debugging mode
!   dirichlet: logical, true for Dirichlet BVP
! Returns:
!   rhs: right hand size of integral equation

   use geometry_mod, only: k0, k, nd, nbk, zk, z, bounded
   implicit none
   logical, intent(in) :: debug, dirichlet
   real(kind=8), intent(out) :: rhs(nbk+k)
   real(kind=8) :: U_EXACT
   integer :: i, kbod
   
      do i = 1, nbk
	 rhs(i) = U_EXACT(bounded, z(i))
      end do
      
!   constraints for multiply connected domain
      do kbod = 1, k
         rhs(nbk + kbod) = 0.d0
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
   complex(kind=8) :: zdis
   integer :: n
   real(kind=8) :: theta, r, A, B
   
!
! harmonic mode
      n = 2
      A = 1.d0
      B = 0.5d0
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

!----------------------------------------------------------------------

subroutine GET_SOL_TAR(ntar, z_tar, mu, A_log, u_tar)

! Calculate solution at target points and check accuracy
! Inputs:
!   ntar: number of target points
!   z_tar: location of target points in complex plane
!   mu: density of integral operator
!   A_log: strength of log sources
!   bounded: logical whether domain is bounded or not
! Returns:
!   u_tar: solution

   use geometry_mod, only: k0, k, nd, nbk, pi, h, eye, z, dz, bounded
   implicit none
   integer, intent(in) :: ntar
   real(kind=8), intent(in) :: mu(nbk), A_log(k)
   complex(kind=8), intent(in) :: z_tar(nbk)
   real(kind=8), intent(out) :: u_tar(nbk)
!
! local
   real(kind=8) :: err, u_ex, U_EXACT
   integer :: i, itar
   complex(kind=8) :: zcauchy, z2pii 

      z2pii = 1.d0/(2.d0*pi*eye)

      err = 0.d0
      
      do itar = 1, ntar
      
         u_tar(itar) = 0.d0
         
         do i = 1, nbk
            zcauchy = mu(i)*dz(i)/(z(i) - z_tar(itar))
            zcauchy = h*zcauchy*z2pii
            u_tar(itar) = u_tar(itar) + dreal(zcauchy)
         end do
         
         u_ex = U_EXACT(bounded, z_tar(itar))
         err = max(err,dabs(u_ex-u_tar(itar)))
!!!         call PRINF('itar = *', itar, 1)
!!!         call PRIN2('   u_exact = *', u_ex, 1)
!!!         call PRIN2('   u_tar = *', u_tar(itar), 1)
!!!         call PRIN2('   diff = *', u_ex-u_tar(itar), 1)
         
      end do
      call PRIN2 (' MAX ERROR AT TARGET POINTS = *', err, 1)

end subroutine GET_SOL_TAR
