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
   use laplace_system_mod

   implicit none
!
! Integral equation density and log sources
   real(kind=8) :: mu(nmax), mu_res(ibeta*nmax), A_log(kmax)
!
! Right hand side to integral equation 
   real(kind=8) :: rhs(nmax+kmax)
!
! Solution on grid
   real(kind=8) :: u_grd(ngrd_max), umin, umax, &
				ugrd_bad(ngrd_max), umin_bad, umax_bad
!
! Matrix equation solution from GMRES
   real(kind=8) :: soln(nmax+kmax)

!
! Problem type
   logical :: debug
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
   call INITIALIZE(debug)
   call INIT_HOLE_GEO() 
   call BUILD_DOMAIN()
   call BAD_DOMAIN_BNDRY()
   call BUILD_CLOSEEVAL_GRID()
   call BUILD_GRID(i_grd, x_grd, y_grd)
   
!
! Get target points
  ! if (debug) call GET_TARGETS(ntar, z_tar)
   
!
! Set boundary conditions
   call GET_BCS(debug, rhs)   
!   call PRIN2(' rhs = *', rhs, nbk)
   
!
! Solve integral equation
   call SOLVE (rhs, soln, mu, A_log)
   
!
! Get solution on the grid
   call RESAMPLE_DOMAIN()
   call BUILD_BARNETT(mu,mu_res)
   call GET_SOL_GRID(mu, A_log, i_grd, x_grd, y_grd, u_grd, umin, umax)
   call GET_CLOSEEVAL_SOL_GRID(mu_res,A_log,ugrd_bad,umin_bad,umax_bad)
! Check solution at target points
   if (debug) then
    !  call GET_SOL_TAR(ntar, z_tar, mu, A_log, u_tar)
      call CHECK_ERROR_GRID(i_grd, x_grd, y_grd, u_grd, umin, umax)
	  call CHECK_ERROR_CLOSEEVAL_GRID(ugrd_bad,umin_bad, umax_bad)
   end if

end program LAPLACE_2D

!----------------------------------------------------------------------

subroutine INITIALIZE(debug)
!
! reads in k0, k, nd from module and calculates nbk
! describes problem type
! Returns:
!   debug :: true if in debugging mode
!   bounded :: true if bounded domain, false if unbounded
!   dirichlet :: true if Dirichlet BVP, false if Neumann
   use geometry_mod, only: pi, eye, kmax, npmax, nbk, k0, k, nd, h, &
                           bounded, nx, ny, ngrd_max, nr, ntheta, nb, &
						   ndres, nbkres, hres, ibeta, ig, PARAM_DUMP
   use laplace_system_mod, only: dirichlet
   implicit none
   logical, intent(out) :: debug

!
! initialize constants
      pi = 4.d0*datan(1.d0)
      eye = dcmplx(0.d0,1.d0)

!
! initialize number of holes and points per hole
      k0 = 0
      k = 3
      nd = 256
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
      nx = 200
      ny = 200
      if (nx*ny > ngrd_max) then
         print *, 'Too many grid points!'
         print *, 'nx*ny = ', nx*ny
         print *, 'ngrd_max = ', ngrd_max
         stop
      end if

! initialize resampled boundary points
	  ndres = ibeta*nd
	  nbkres = (k + 1 -k0)*ndres
	  hres = 2.d0*pi/ndres
	  ig = 10

! initialize number of boxes 
	  nb = nd/4
! initialize close evaluation grid
	  nr = 5
	  ntheta = 250
		
	  open(unit = 31, file="mat_plots/bad_param.m")
	  call PARAM_DUMP((k-k0+1),nr, ntheta, 31)
	  close(31)

end subroutine INITIALIZE

!----------------------------------------------------------------------
subroutine INIT_HOLE_GEO()
!
! Initialize size of system and parameters for domain
! Returns:
! ak, bk: major and minor axes of elliptical holes
! zk: geometric centre of holes
! ncyc: number of starfish arms
!
use geometry_mod, only: k0, k, ak, bk, ncyc, zk
implicit none
ak(1) = 2.d0
bk(1) = 1.8d0
ncyc(1) = 0
zk(1) = dcmplx(0.d0, 0.0d0)
ak(2) = 0.5d0
bk(2) = 0.3d0
ncyc(2) = 0
zk(2) = dcmplx(-0.75d0, 0.0d0)
ak(3) = 0.1d0
bk(3) = 0.2d0
ncyc(3) = 0
zk(3) = dcmplx(0.75d0, 0.0d0)
ak(4) = 0.1d0
bk(4) = 0.05d0
ncyc(4) = 0
zk(4) = dcmplx(0.0d0, 0.75d0)
end subroutine INIT_HOLE_GEO
!----------------------------------------------------------------------!----------------------------------------------------------------------
  
subroutine GET_BCS(debug, rhs)

! Constructs right hand side of integral equation
! Inputs:
!   debug: logical for debugging mode
!   dirichlet: logical, true for Dirichlet BVP
! Returns:
!   rhs: right hand size of integral equation

   use geometry_mod, only: k0, k, nd, nbk, zk, z, bounded
   use laplace_system_mod, only: dirichlet
   implicit none
   logical, intent(in) :: debug
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
            U_EXACT = U_EXACT + dreal(1.d0/(z-zk(kbod+1)))
         end do
      elseif (k == 1) then
         U_EXACT = (A*dcos(n*theta) + B*dsin(n*theta)) / r**n
      else
         U_EXACT = 0.d0
         do kbod = 1, k
            U_EXACT = U_EXACT + dreal(1.d0/(z-zk(kbod)))
         end do         
      end if

end function U_EXACT
!-----------------------------------------------------------------------------

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
                           nx, ny, zk, ds_dth, REAL_GRID_DUMP, &
                           z_res, dz_res, ibeta, RESAMPLE_DOMAIN, &
                           n_close, i_close, j_close, ic_pnt, X_DUMP 

   implicit none
   integer, intent(in) :: i_grd(nx,ny)
   real(kind=8), intent(in) :: mu(nbk), A_log(k), x_grd(nx, ny), &
                               y_grd(nx, ny)
   real(kind=8), intent(out) :: u_grd(nx, ny)
!
! local variables
   integer :: i, j, istart, kbod
   real(kind=8) :: umin, umax, mu_res(ibeta*nbk)
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
            if (i_grd(i,j) .eq. 2) then   
               target(1, istart) = x_grd(i,j)
               target(2, istart) = y_grd(i,j)
               istart = istart + 1 
            end if
         end do
      end do
      ntarget = istart - 1
      call PRINF (' Number of active points in grid = *', ntarget, 1)
      
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
	  
! For points far enough from boundary, unpack into grid
      umax = -1.d10
      umin = 1.d10
      istart = 1
      do i = 1, nx
         do j = 1, ny
            z_grid = dcmplx(x_grd(i, j), y_grd(i, j))
            if (i_grd(i, j) .eq. 2) then 
               u_grd(i, j) = dreal(pottarg(istart))
               do kbod = 1, k
                  u_grd(i, j) = u_grd(i, j) & 
                    + A_log(kbod)*dlog(cdabs(z_grid - zk(kbod + 1 - k0)))
               end do
			   
               umax = max(umax, u_grd(i, j))
               umin = min(umin, u_grd(i, j))
               istart = istart + 1 
             else
               u_grd(i, j) = -100.d0
            end if
         end do
      end do
      
!
! For points close to a boundary curve, use Barnett's corrections

!!!      call RESAMPLE_DOMAIN ()
!!!      do kbod = k0, k
!!!      call PRINF('kbod = *', kbod, 1)
!!!      do indx = ic_pnt(kbod-k0+1), ic_pnt(kbod-k0+2) - 1
!!!         i = i_close(indx)
!!!         j = j_close(indx)
!!!         call PRINF('indx = *', indx, 1)
!!!         call PRINF(' i = *', i, 1)
!!!         call PRINF(' j = *', j, 1)
!!!         call PRINF('i_grd = *', i_grd(i,j), 1)
!!!      end do
!!!   end do
      
      call PRIN2('Min solution on grid = *', umin, 1)
      call PRIN2('Max solution on grid = *', umax, 1)

      open(unit = 31, file = 'mat_plots/ugrid.m')
      open(unit = 32, file = 'mat_plots/ugrd_plot.m')
      write(31, *) 'ulim = ['
      write(31, '(2(D15.6))') umin, umax
      write(31, *) '];'

	  write(32, *) 'ulim = ['
      write(32, '(2(D15.6))') umin, umax
      write(32, *) '];'


      call REAL_GRID_DUMP(u_grd, 31)
	  call X_DUMP(u_grd, nx*ny, 32)

      close(31)
	  close(32)
         
end subroutine GET_SOL_GRID

!-----------------------------------------------------------------------

subroutine GET_CLOSEEVAL_SOL_GRID(mu_res, A_log,ugrd_bad, & 
							umin_bad, umax_bad)

! Calculate solution at grid points in the bad region. 

   use geometry_mod, only: k0, k, pi, h, nd, nbk, eye, z, dz, bounded, &
                           nb,zgrd_bad, xgrd_bad, ygrd_bad, z0_box, &
						   nx, ny, nr, ntheta, ndres, nbkres, & 
						   z_res, dz_res, hres, ibeta,zk, n_neigh, &
						   neigh_boxes,X_DUMP, GET_NEAR_POINTS 


   use laplace_system_mod, only: cm, p

   implicit none
   real(kind=8), intent(in) :: mu_res(nbkres), A_log(k-k0)
   real(kind=8), intent(out) :: ugrd_bad((k-k0+1)*nr*ntheta), &
							umin_bad, umax_bad


! local variables
   integer :: i, j, ipoint, kbod, im, ibox(nd), &
			  iibox, istart, llimit, rlimit, icl, jcl, ntarget
   complex(kind=8):: zpoint, z0, zcauchy, z2pii
   real(kind=8):: magdz

! FMM work arrays
   integer :: iprec, ifpot, ifgrad, ifhess,  &
              ifpottarg, ifgradtarg, ifhesstarg, ier
   real(kind=8) :: source(2,nbkres), & 
					target(2,(k-k0+1)*nr*ntheta), targ(2,2)
   complex(kind=8) :: dipstr(nbkres), pot(nbkres), grad(2,nbkres), &
                      hess(3,nbkres), pottarg((k-k0+1)*nr*ntheta), & 
				      gradtarg(2, (k-k0+1)*nr*ntheta), &
                      hesstarg(3, (k-k0+1)*nr*ntheta)

	
	
	umin_bad = 1.d10
	umax_bad = -1.d10
    z2pii = 1.d0/(2.d0*pi*eye)

	 do j = 1, ntheta	
	  	ibox(j) = -1
	  	do iibox = 1, nb
		  if((j.ge.(iibox-0.5d0)*ntheta/nb) .and. &
				j.lt.(iibox + 0.5d0)*ntheta/nb) then
				ibox(j) = iibox
			end if
	  	end do			
		if(ibox(j).eq.-1) then
			ibox(j) = nb
		end if
	 end do

!
!   Define grid target points

	ntarget = (k - k0 + 1)*nr*ntheta
	do i = 1, ntarget  
               target(1, i) = xgrd_bad(i)
               target(2, i) = ygrd_bad(i)
              
    end do

    call PRINF (' Number of active points in the bad grid = *', ntarget, 1)
      
!
!   Assemble arrays for FMM call
	
    do i = 1, nbkres
         source(1, i) = dreal(z_res(i))
         source(2, i) = dimag(z_res(i))
         dipstr(i) = -1.d0*hres*mu_res(i)*z2pii*dz_res(i)
    end do

! set parameters for FMM routine lfmm2dparttarg
	
      iprec = 5   ! err < 10^-14
      ifpot = 1
      ifgrad = 0
      ifhess = 0
      ifpottarg = 1
      ifgradtarg = 0
      ifhesstarg = 0
      
! call FMM

      call PRINI(0, 13)
      call zfmm2dparttarg(ier, iprec, nbkres, source, &
                          dipstr, ifpot, pot, ifgrad,  &
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
	  
! For points far enough from boundary, unpack into grid
   
         
	do kbod = k0, k
		do i = 1, nr
			do j = 1,ntheta
				ipoint = (kbod - k0)*nr*ntheta + (i-1)*ntheta + j	
				zpoint = zgrd_bad(ipoint)
				z0 = z0_box(kbod - k0 + 1,ibox(j))	
				ugrd_bad(ipoint) = dreal(pottarg(ipoint))
				do im = 1, p
					ugrd_bad(ipoint) = ugrd_bad(ipoint) + &
						dreal(cm(kbod - k0 + 1, ibox(j), im)*((zpoint - z0)**(im-1)))			
				end do
			 	 
				do icl = 1, n_neigh(kbod - k0 + 1, ibox(j))
						jcl = neigh_boxes(kbod - k0 + 1, ibox(j), icl)
       					zcauchy = mu_res(jcl)*dz_res(jcl)/ &
						(z_res(jcl) - zpoint)
						zcauchy = hres*zcauchy*z2pii
						ugrd_bad(ipoint) = ugrd_bad(ipoint) &
								- dreal(zcauchy)
				end do
				umin_bad = min(umin_bad, ugrd_bad(ipoint))
				umax_bad = max(umax_bad, ugrd_bad(ipoint))
			end do			
		end do
	end do	

	open(unit = 31, file = 'mat_plots/ugrd_bad_plot.m')

      write(31, *) 'ulim = ['
      write(31, '(2(D15.6))') umin_bad, umax_bad
      write(31, *) '];'

      call X_DUMP(ugrd_bad,(k-k0+1)*nr*ntheta, 31)

      close(31)


end subroutine GET_CLOSEEVAL_SOL_GRID
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
            if (i_grd(i,j) .eq. 2) then  
               u_ex = U_EXACT(bounded, z_grid)
               err = max(err, dabs(u_ex - u_grd(i, j)))			  
         !      call PRIN2 ('u_ex = *', u_ex, 1)
         !      call PRIN2 ('  u_grd = *', u_grd(i,j), 1)
            end if
         end do
      end do
      
      call PRIN2 ('Max error on grid = *', err, 1)
      call PRIN2 ('Max relative error on grid = *', err / u_inf, 1)
         
end subroutine CHECK_ERROR_GRID
!------------------------------------------------------------------

subroutine CHECK_ERROR_CLOSEEVAL_GRID(ugrd_bad,umin_bad,umax_bad)

! Checks error in solution at grid points 
! Inputs:
!   i_grd(i,j): flags whether grid point in domain or not
!   x_grd, y_grd: grid points
!   u_grd: solution
!   umin, umax: min and max of solution

   use geometry_mod, only: k, k0, nr, nd,ntheta, zgrd_bad, bounded, &
							Y_PLOT
   implicit none
   real(kind=8), intent(in) :: ugrd_bad((k-k0+1)*nr*ntheta), &
		 umin_bad, umax_bad
!
! local variables
   integer :: i, j, kbod, ipoint, nb
   real(kind=8) :: err, u_ex_bad, u_inf, U_EXACT
   complex(kind=8) :: z_grid


      err = 0.d0
      u_inf = max(dabs(umin_bad), dabs(umax_bad))
	 
	
	  do kbod = k0, k
      	do i = 1, nr
        	 do j = 1, ntheta
				ipoint = kbod*nr*ntheta + (i-1)*ntheta + j
				z_grid = zgrd_bad(ipoint)
               	u_ex_bad = U_EXACT(bounded, z_grid)	
               	err = max(err, dabs(u_ex_bad - ugrd_bad(ipoint)))
				!print 1000, kbod, i, j,u_ex_bad, ugrd_bad(ipoint)
				!1000 format(I3,I3,I5,2(D15.6))
                end do
         end do
      end do
      
      call PRIN2 ('Max error on the bad grid = *', err, 1)
      call PRIN2 ('Max relative error on the bad grid = *', err / u_inf, 1)
	
         
	  close(51)
end subroutine CHECK_ERROR_CLOSEEVAL_GRID
!----------------------------------------------------------------------------
