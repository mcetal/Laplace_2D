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
! Solution on grid
   real(kind=8) :: u_grd(ngrd_max), umin, umax
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
   call BUILD_GRID(i_grd, x_grd, y_grd)
   
!
! Get target points
   if (debug) call GET_TARGETS(ntar, z_tar)
   
!
! Set boundary conditions
   call GET_BCS(debug, dirichlet, rhs)   
!   call PRIN2(' rhs = *', rhs, nbk)
   
!
! Solve integral equation
   igwork(1) = maxl
   call SOLVE (maxl, rhs, lrwork, liwork, dirichlet, soln, mu, A_log, &
               gmwork, igwork)
   
!
! Get solution on the grid
   call GET_SOL_GRID(mu, A_log, i_grd, x_grd, y_grd, u_grd, umin, umax)
!
! Check solution at target points
   if (debug) then
      call GET_SOL_TAR(ntar, z_tar, mu, A_log, u_tar)
      call CHECK_ERROR_GRID(i_grd, x_grd, y_grd, u_grd, umin, umax)
   end if

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
                           bounded, nx, ny, ngrd_max
   implicit none
   logical, intent(out) :: debug, dirichlet

!
! initialize constants
      pi = 4.d0*datan(1.d0)
      eye = dcmplx(0.d0,1.d0)

!
! initialize number of holes and points per hole
      k0 = 1
      k = 2
      nd = 1024
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
      
!
! initialize grid size
      nx = 10
      ny = 10
      if (nx*ny > ngrd_max) then
         print *, 'Too many grid points!'
         print *, 'nx*ny = ', nx*ny
         print *, 'ngrd_max = ', ngrd_max
         stop
      end if

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
      zk(1) = dcmplx(-1.d0, -0.5d0)
      
      ak(2) = 0.5d0
      bk(2) = 0.1d0
      ncyc(2) = 4
      zk(2) = dcmplx(1.d0,0.5d0)
      
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
   use geometry_mod, only: pi, k0, k, zk, ak, bk, Z_PLOT, xmin, xmax, &
                           ymin, ymax
   implicit none
   integer, intent(in) :: ntar
   complex(kind=8), intent(out) :: z_tar(ntar)
! local
   integer :: i
   real(kind=8) :: dth, theta, a_tar, b_tar
   character(32) :: options
   complex(kind=8) :: z_centre, z_corner
   
      dth = 2.d0*pi/ntar

      if ( (k0==0) .and. (k==0) ) then
         a_tar = 0.5d0*ak(1)
         b_tar = 0.5d0*bk(1)
         z_centre = zk(1)
      elseif ( (k0==0) .and. (k==1) ) then
         a_tar = 0.5d0*(ak(1) + ak(2))
         b_tar = 0.5d0*(bk(1) + bk(2))
         z_centre = zk(1)
      elseif ( (k0==1) .and. (k==1) ) then
         a_tar = 2.d0*max(ak(1), bk(1))
         b_tar = a_tar
         z_centre = zk(1)
      elseif (k0 == 1) then
         z_centre = 0.5d0*dcmplx(xmin + xmax, ymin + ymax)
         z_corner = dcmplx(xmax, ymax)
         a_tar = cdabs(z_corner - z_centre)
         b_tar = a_tar   
      else
         print *, 'Cannot calculate target points for this geometry'
         print *, 'Errors in solution check may result'
      end if
      
      do i = 1, ntar
         theta = dth*(i-1.d0)
         z_tar(i) = z_centre + dcmplx(a_tar*dcos(theta), &
                                      b_tar*dsin(theta))
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
   integer :: n, kbod
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

      if ((bounded) .and. (k == 0)) then 
         U_EXACT = (A*dcos(n*theta) + B*dsin(n*theta)) * r**n
      elseif (bounded) then
         U_EXACT = 0.d0
         do kbod = 1, k
            U_EXACT = U_EXACT + 1.d0/(z-zk(kbod+1))
         end do
      elseif (k == 1) then
         U_EXACT = (A*dcos(n*theta) + B*dsin(n*theta)) / r**n
      else
         U_EXACT = 0.d0
         do kbod = 1, k
            U_EXACT = U_EXACT + 1.d0/(z-zk(kbod))
         end do         
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
      call PRIN2 ('Max error at target points = *', err, 1)

end subroutine GET_SOL_TAR

!----------------------------------------------------------------------

subroutine GET_SOL_GRID(mu, A_log, i_grd, x_grd, y_grd, u_grd, umin, umax)

! Calculate solution at grid points 
! Inputs:
!   mu: density of integral operator
!   A_log: strength of log sources
!   i_grd(i,j): flags whether grid point in domain or not
!   x_grd, y_grd: grid points
! Returns:
!   u_grd: solution
!   umin: minimum solution value
!   umax: maximum solution value

   use geometry_mod, only: k0, k, nd, nbk, pi, h, eye, z, dz, bounded, &
                           nx, ny, zk, ds_dth, REAL_GRID_DUMP
   implicit none
   integer, intent(in) :: i_grd(nx,ny)
   real(kind=8), intent(in) :: mu(nbk), A_log(k), x_grd(nx, ny), &
                               y_grd(nx, ny)
   real(kind=8), intent(out) :: u_grd(nx, ny)
!
! local variables
   integer :: i, j, istart, kbod
   real(kind=8) :: umin, umax
   complex(kind=8) :: z_grid
!
! FMM work arrays
   integer :: iprec, ifcharge, ifdipole, ifpot, ifgrad, ifhess, ntarget, &
              ifpottarg, ifgradtarg, ifhesstarg, ier
   real(kind=8) :: source(2,nbk), dipvec(2,nbk), target(2, nx*ny)
   complex(kind=8) :: charge(nbk), dipstr(nbk), pot(nbk), grad(2,nbk), &
                      hess(3,nbk), pottarg(nx*ny), gradtarg(2, nx*ny), &
                      hesstarg(3, nx*ny)

!
!   Define grid target points
      istart = 1
      do i = 1, nx
         do j = 1, ny
            target(1, istart) = x_grd(i,j)
            target(2, istart) = y_grd(i,j)
            istart = istart + 1 
         end do
      end do
      ntarget = istart - 1
      
!
!   Assemble arrays for FMM call
      do i = 1, nbk
         source(1, i) = dreal(z(i))
         source(2, i) = dimag(z(i))
         dipvec(1,i) = dreal(-eye*dz(i))/ds_dth(i)
         dipvec(2,i) = dimag(-eye*dz(i))/ds_dth(i)
         charge(i) = 0.d0
         dipstr(i) = h*mu(i)*ds_dth(i)/(2.d0*pi)
      end do

! set parameters for FMM routine DAPIF2
	
      iprec = 5   ! err < 10^-14
      ifcharge = 0 ! no charges, only dipoles
      ifdipole = 1
      ifpot = 1
      ifgrad = 0
      ifhess = 0
      ifpottarg = 1
      ifgradtarg = 0
      ifhesstarg = 0
      
! call FMM

      call PRINI(0, 13)
      call lfmm2dparttarg(ier, iprec, nbk, source, ifcharge, charge, &
                          ifdipole, dipstr, dipvec, ifpot, pot, ifgrad,  &
                          grad, ifhess, hess, ntarget, target, ifpottarg, &
                          pottarg, ifgradtarg, gradtarg, ifhesstarg, hesstarg)
      call PRINI(6, 13)
      
      if (ier.eq.4) then
         print *, 'ERROR IN FMM: Cannot allocate tree workspace'
         stop
      else if(ier.eq.8) then
         print *, 'ERROR IN FMM: Cannot allocate bulk FMM workspace'
         stop
      else if(ier.eq.16) then
         print *, 'ERROR IN FMM: Cannot allocate multipole expansion workspace' 
         stop
      end if
	  
! unpack into grid
      umax = -1.d10
      umin = 1.d10
      istart = 1
      do i = 1, nx
         do j = 1, ny
            z_grid = dcmplx(x_grd(i, j), y_grd(i, j))
            if (i_grd(i, j) .eq. 1) then  
               u_grd(i, j) = dreal(pottarg(istart))
               do kbod = 1, k
                  u_grd(i, j) = u_grd(i, j) & 
                    + A_log(kbod)*dlog(cdabs(z_grid - zk(kbod + 1 - k0)))
               end do
               umax = max(umax, u_grd(i, j))
               umin = min(umin, u_grd(i, j))
             else
               u_grd(i, j) = -100.d0
            end if
            istart = istart + 1 
         end do
      end do
      
      call PRIN2('Min solution on grid = *', umin, 1)
      call PRIN2('Max solution on grid = *', umax, 1)

      open(unit = 31, file = 'mat_plots/ugrid.m')

      write(31, *) 'ulim = ['
      write(31, '(2(D15.6))') umin, umax
      write(31, *) '];'

      call REAL_GRID_DUMP(u_grd, 31)

      close(31)
         
end subroutine GET_SOL_GRID

!----------------------------------------------------------------------

subroutine CHECK_ERROR_GRID(i_grd, x_grd, y_grd, u_grd, umin, umax)

! Checks error in solution at grid points 
! Inputs:
!   i_grd(i,j): flags whether grid point in domain or not
!   x_grd, y_grd: grid points
!   u_grd: solution
!   umin, umax: min and max of solution

   use geometry_mod, only: nx, ny, bounded
   implicit none
   integer, intent(in) :: i_grd(nx,ny)
   real(kind=8), intent(in) :: x_grd(nx, ny), y_grd(nx, ny), &
                               u_grd(nx, ny), umin, umax
!
! local variables
   integer :: i, j
   real(kind=8) :: err, u_ex, u_inf, U_EXACT
   complex(kind=8) :: z_grid

      err = 0.d0
      u_inf = max(dabs(umin), dabs(umax))
      do i = 1, nx
         do j = 1, ny
            z_grid = dcmplx(x_grd(i, j), y_grd(i, j))
            if (i_grd(i,j) .eq. 1) then  
               u_ex = U_EXACT(bounded, z_grid)
               err = max(err, dabs(u_ex - u_grd(i, j)))
            end if
         end do
      end do
      
      call PRIN2 ('Max error on grid = *', err, 1)
      call PRIN2 ('Max relative error on grid = *', err / u_inf, 1)
         
end subroutine CHECK_ERROR_GRID
