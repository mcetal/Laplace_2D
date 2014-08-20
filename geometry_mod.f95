!
! comman variables etc for 2-D integral equation solvers
! routines to describe the boundary of the domain and the grid
! make sure a directory "mat_plots/" exists
!
module geometry_mod

   implicit none
   
   save

   integer, parameter :: kmax = 10, npmax = 2050, nmax = kmax*npmax
   integer, parameter :: nx_max = 500, ny_max = 500, ngrd_max = nx_max*ny_max
!
! Number of holes, system size
   integer :: k0, k, nd, nbk
!
! Geometry of Domain
   real(kind=8) :: ak(kmax), bk(kmax)
   integer :: ncyc(kmax)
   complex(kind=8) :: zk(kmax)
!
! Type of domain
   logical :: bounded

!
! size of domain
   real(kind=8) :: xmin, xmax, ymin, ymax
!
! Constants and mesh spacing
   real(kind=8) :: pi, h
   complex(kind=8) :: eye
!
! Discretization of domain
   real(kind=8) :: x(nmax), y(nmax), ds_dth(nmax), kappa(nmax), &
                   theta(nmax)
   complex(kind=8) :: z(nmax), dz(nmax)
   
!
! Boundary of "bad" part of domain - i.e. region close enough to boundary that
! trapezoid rule breaks down
   real(kind=8) :: x_bad(nmax), y_bad(nmax), ds_bad(nmax)
   integer :: nb
   complex(kind=8) :: z_bad(nmax), dz_bad(nmax), z0_box(kmax,nmax/5)

!
! Pointer arrays to points in grid that are in "close" region, and which 
! contour they are close to
   integer :: n_close(kmax), i_close(ngrd_max), j_close(ngrd_max), &
              ic_pnt(kmax+1), neigh_boxes(kmax,nmax/5,nmax/5), &
			  n_neigh(kmax,nmax/5) 
!
! resampled domain variables
   integer, parameter :: ibeta = 4, p = 20, ialpha = 10
   integer ::  ndres, nbkres, ig
   real(kind=8) :: g, hres, x_res(ibeta*nmax), y_res(ibeta*nmax), & 
				  ds_res(ibeta*nmax)
   complex(kind=8) :: z_res(ibeta*nmax), dz_res(ibeta*nmax)
   integer :: inear_box(nmax/5, ibeta*nmax)	


! Grid variables
   integer :: nx, ny, i_grd(ngrd_max), nr, ntheta
   real(kind=8) :: x_grd(ngrd_max), y_grd(ngrd_max)
   real(kind=8) :: xgrd_bad(ngrd_max), ygrd_bad(ngrd_max)
   complex(kind=8):: zgrd_bad(ngrd_max)!, dzgrd_bad(ngrd_max)
      
contains
!----------------------------------------------------------------------
subroutine CURVE_PARAM(th, kbod, xp, yp, xdot, ydot, xddot, yddot)
!
! The initial curve parameterization is given in term of th
! Inputs:
! th: parameter value
! kbod: the particular body index point is on
! Returns:
! xp, yp: coordinates of point
! xdot, ydot: derivative of position wrt parameter
! xddot, yddot: 2nd derivative of position wrt parameter
!
	implicit none
	real(kind=8), intent(in) :: th
	integer, intent(in) :: kbod
	real(kind=8), intent(out) :: xp, yp, xdot, ydot, xddot, yddot
	real(kind=8) :: ai, bi, N, a2, b2, rnorm, radius, rdot, rddot, R, &
	anu, den, ddot, dddot, cs, sn, eps, snn, csn, thetax
	N = ncyc(kbod+1-k0)
	ai = ak(kbod+1-k0)
	bi = bk(kbod+1-k0)
	thetax = 0.d0
	if (N.eq.0) then
	cs = dcos(th-thetax)
	sn = dsin(th-thetax)
	a2 = (ai)**2
	b2 = (bi)**2
	rnorm = dsqrt(b2*cs**2 + a2*sn**2)
	radius = ai*bi/rnorm
	rdot = -ai*bi*cs*sn*(-b2+a2)/rnorm**3
	rddot = -ai*bi*(2.d0*a2*b2*cs**2*sn**2 &
		+a2**2*cs**4 + a2*b2-a2**2+b2**2*cs**4 &
		-2.d0*b2**2*cs**2)/rnorm**5
	xp = radius*dcos(th)
	yp = radius*dsin(th)
	cs = dcos(th)
	sn = dsin(th)
	xdot = rdot*cs - radius*sn
	ydot = rdot*sn + radius*cs
	xddot = rddot*cs - 2.d0*rdot*sn - radius*cs
	yddot = rddot*sn + 2.d0*rdot*cs - radius*sn
	elseif (N.lt.0) then
	R = ai
	anu = bi
	den = (1.d0-2.d0*anu*dcos(2.d0*th)+anu**2) &
	*dsqrt(1.d0+anu**2)
	ddot = 4.d0*anu*dsin(2.d0*th)*dsqrt(1.d0+anu**2)
	dddot = 8.d0*anu*dcos(2.d0*th)*dsqrt(1.d0+anu**2)
	xp = (1.d0-anu**2)*(1.d0-anu)*R*dsqrt(2.d0)*dcos(th)/den
	yp = (1.d0-anu**2)*(1.d0+anu)*R*dsqrt(2.d0)*dsin(th)/den
	xdot = -(1.d0-anu**2)*(1.d0-anu)*R*dsqrt(2.d0)*dsin(th)/den &
		-xp*ddot/den
	xddot = -xp &
	+(1.d0-anu**2)*(1.d0-anu)*R*dsqrt(2.d0) &
	*dsin(th)*ddot/den**2 &
	- xdot*ddot/den - xp*dddot/den &
	+ xp*ddot**2/den**2
	ydot = (1.d0-anu**2)*(1.d0+anu)*R*dsqrt(2.d0)*dcos(th)/den &
		-yp*ddot/den
	yddot = -yp &
		-(1.d0-anu**2)*(1.d0+anu)*R*dsqrt(2.d0) &
		*dcos(th)*ddot/den**2 &
		- ydot*ddot/den - yp*dddot/den &
		+ yp*ddot**2/den**2
	else
		eps = bi
		cs = dcos(th)
		sn = dsin(th)
		snn = dsin(N*th)
		csn = dcos(N*th)
		xp = ai*cs + bi*csn
		yp = ai*sn - bi*snn
		xdot = -ai*sn - N*bi*snn
		xddot = -ai*cs - N**2*bi*csn
		ydot = ai*cs - N*bi*csn
		yddot = -ai*sn + N**2*bi*snn
	end if
end subroutine CURVE_PARAM
!----------------------------------------------------------------------
   
!----------------------------------------------------------------------

subroutine BUILD_DOMAIN()
!
! Input
!    ak, bk:  major and minor axes of elliptical holes
!    zk:  geometric centre of holes
!    ncyc:  number of starfish arms
! Output
!    theta:  array of discrete theta values
!    z:  points on holes
!    dz:  derivative of z wrt theta
!    ds_dth:  |dz|
!    kappa:  curvature
!
   implicit none
   integer :: i, kbod, istart
   real(kind=8) xp, yp, xdot, ydot, xddot, yddot, th, delx, dely
   character(32) :: options
   
!
! open unit for matlab plotting
      open(unit = 31, file = 'mat_plots/geo_contours.m')
      options = '''k'',''LineWidth'',2'
   
!
! construct theta
      do i = 1, nd
         theta(i) = h*(i-1.d0)
      end do
      
      call PRIN2('in build domain, h = *', h, 1)
!
! construct holes
      istart = 0
      xmin = 1.d10
      xmax = -1.d10
      ymin = 1.d10
      ymax = -1.d10
      do kbod = k0, k
         do i = (kbod-k0)*nd+1, (kbod-k0)*nd+nd
            th = theta(i-(kbod-k0)*nd)
            call CURVE_PARAM (th, kbod, xp, yp, xdot, ydot, xddot, yddot)
            z(i) = zk(kbod+1-k0) + dcmplx(xp, yp)
            ds_dth(i) = cdabs(dcmplx(xdot,ydot))    
            if (kbod.eq.0) then
               dz(i) = dcmplx(xdot, ydot)
               kappa(i) = (xdot*yddot-ydot*xddot)/ds_dth(i)**3
             else
               dz(i) = -dcmplx(xdot, ydot)
               kappa(i) = -(xdot*yddot-ydot*xddot)/ds_dth(i)**3
            endif 
            xmin = min(xmin, dreal(z(i)))
            xmax = max(xmax, dreal(z(i)))
            ymin = min(ymin, dimag(z(i)))
            ymax = max(ymax, dimag(z(i)))
         end do   
         call Z_PLOT(z((kbod-k0)*nd+1), nd, options, 31)
!         call PRIN2('z = *', z((kbod-k0)*nd+1), 2*nd)
!         call PRIN2('dz = *', dz((kbod-k0)*nd+1), 2*nd)
!         call PRIN2('kappa = *', kappa((kbod-k0)*nd+1), nd)
         istart = istart + nd  
      end do
      
      close(31)

!
!  Print out dimensions for domain; increase box size if problem is unbounded
      call PRIN2('In BUILD_DOMAIN, xmin = *', xmin, 1)
      call PRIN2('In BUILD_DOMAIN, xmax = *', xmax, 1)
      call PRIN2('In BUILD_DOMAIN, ymin = *', ymin, 1)
      call PRIN2('In BUILD_DOMAIN, ymax = *', ymax, 1)
       
      if (.not. bounded) then
         print *, 'Domain is unbounded, increasing size of domain for grid:'
         delx = xmax - xmin
         dely = ymax - ymin
         xmax = xmax + delx
         xmin = xmin - delx
         ymax = ymax + dely
         ymin = ymin - dely
         call PRIN2('   New xmin = *', xmin, 1)
         call PRIN2('   New xmax = *', xmax, 1)
         call PRIN2('   New ymin = *', ymin, 1)
         call PRIN2('   New ymax = *', ymax, 1)
      end if

end subroutine BUILD_DOMAIN

   
!----------------------------------------------------------------------


subroutine BAD_DOMAIN_BNDRY()
!
! For each contour, define a curve which dilineates where it is safe to use
! the trapezoid rule to evaluate the layer potentials within the domain. 
! Reference: 
! Alex Barnett, EVALUATION OF LAYER POTENTIALS CLOSE TO THE BOUNDARY FOR LAPLACE AND 
! HELMHOLTZ PROBLEMS ON ANALYTIC PLANAR DOMAINS
! SIAM J. Sci. Stat. Comput. 2012
! This subroutine defines the curve Gamma_alpha in this paper
!
! Input
!    z:  boundary of domain
! Output
!    z_bad:  Gamma_alpha for each contour
!    dz_bad:  d Gamma_alpha
!    z0_box: centres for boxes 
!
! local variables
   implicit none
   integer :: i, ibox, kbod, istart, istartb, kmode
   real(kind=8) alpha_bad
   complex(kind=8) theta
   character(32) :: options, optionsb
   
!
! FFT work arrays
   real(kind=8) wsave(4*nd + 15)
   complex(kind=8) :: zf(nd)
   
!
! open unit for matlab plotting
      open(unit = 31, file = 'mat_plots/geo_bad_contours.m')
      options = '''r'',''LineWidth'',1'
      open(unit = 32, file = 'mat_plots/geo_bad_box.m')
      optionsb = '''r*'',''LineWidth'',1'

      call DCFFTI (nd, wsave)
      
      alpha_bad = 2.d0*pi/nb
      call PRIN2 (' alpha_bad = *', alpha_bad, 1)
      istart = 0
      istartb = 0
      do kbod = k0, k
         do i = 1, nd
            zf(i) = z(istart+i)
         end do
         call DCFFTF (nd, zf, wsave)
!
! first use these Fourier coefficients to calculate the box centres
         do ibox = 1, nb
            if (kbod .eq. 0) then 
               theta = ibox*alpha_bad + 0.5d0*eye*alpha_bad
             else
               theta = ibox*alpha_bad - 0.5d0*eye*alpha_bad
            end if
            z0_box(kbod+1,ibox) = zf(1)/nd
            do kmode = 1, nd/2 - 1
               z0_box(kbod+1,ibox) = z0_box(kbod+1,ibox) &
                                + zf(kmode+1)*cdexp(eye*kmode*theta)/nd
               z0_box(kbod+1,ibox) = z0_box(kbod+1,ibox) &
                                + zf(nd-kmode+1)*cdexp(-eye*kmode*theta)/nd
            end do
            call Z_PLOT(z0_box(kbod+1,ibox), 1, optionsb, 32)
         end do
!
! now use the Fourier coefficients to calculate the contours of the bad region
         zf(1) = zf(1)/nd
         do kmode = 1, nd/2 - 1
            if (kbod.eq.0) then 
               zf(kmode+1) = zf(kmode+1)*dexp(-kmode*alpha_bad)/nd
               zf(nd-kmode+1) = zf(nd-kmode+1)*dexp(kmode*alpha_bad)/nd
             else
               zf(kmode+1) = zf(kmode+1)*dexp(kmode*alpha_bad)/nd
               zf(nd-kmode+1) = zf(nd-kmode+1)*dexp(-kmode*alpha_bad)/nd
            end if
            zf(nd/2+1) = 0.d0
         end do
         call DCFFTB (nd, zf, wsave)
         do i = 1, nd
            z_bad(istart+i) = zf(i)
         end do
         call FDIFFF(z_bad(istart+1), dz_bad(istart+1), nd, wsave)
         do i = 1, nd
            ds_bad(istart+i) = cdabs(dz_bad(istart+i))
            if (kbod.ne.0) then
               dz_bad(istart+i) = -dz_bad(istart+i)
            end if
         end do
         call Z_PLOT(z_bad(istart+1), nd, options, 31)
         istart = istart + nd
         istartb = istartb + nd/5
      end do
      close(31)
      close(32)

end subroutine BAD_DOMAIN_BNDRY

!----------------------------------------------------------------------


subroutine BUILD_CLOSEEVAL_GRID()
!
! Define a few points within each box in the "bad domain". 
! Reference: 
! Alex Barnett, EVALUATION OF LAYER POTENTIALS CLOSE TO THE BOUNDARY FOR LAPLACE AND 
! HELMHOLTZ PROBLEMS ON ANALYTIC PLANAR DOMAINS
! SIAM J. Sci. Stat. Comput. 2012
! Input
!    z0_box: Array of box centers.
! Output
!    zgrd_bad:  grid point in each box 
!    dzgrd_bad: d grid point in each box
! local variables
   implicit none
   integer :: i, ibox, kbod, istart, istartb, kmode
   integer :: ipoint, jpoint,inum
   real(kind=8) alpha_bad, hbad
   complex(kind=8) theta, th
   character(32) :: options
   
!
! FFT work arrays
   real(kind=8) wsave(4*nd + 15)
   complex(kind=8) :: zf(nd)

! open unit for matlab plotting
      open(unit = 31, file = 'mat_plots/bad_grid_points.m')
      options = '''r'',''LineWidth'',1'

      call DCFFTI (nd, wsave)
      
      alpha_bad = 2.d0*pi/nb
      istart = 0
      istartb = 0
	  
	  hbad = 2.d0*pi/ntheta
      do kbod = k0, k
         do i = 1, nd
            zf(i) = z(istart+i)
         end do
         call DCFFTF (nd, zf, wsave)
!
! first use these Fourier coefficients to calculate the grid points.
       
          
         do ipoint = 1, nr
       		if(kbod.eq.0) then
            	th = eye*ipoint/(nr + 1.d0)*alpha_bad
            else
                th = - eye*ipoint/(nr + 1.d0)*alpha_bad
            end if
			
            do jpoint = 1, ntheta
                    inum = kbod*nr*ntheta + (ipoint-1)*ntheta + jpoint 
                    th = th + hbad  
                    zgrd_bad(inum) = zf(1)/nd
                    do kmode = 1, nd/2 - 1
                        zgrd_bad(inum) = zgrd_bad(inum) &
                                + zf(kmode+1)*cdexp(eye*kmode*th)/nd
                        zgrd_bad(inum) = zgrd_bad(inum) &
                                + zf(nd-kmode+1)*cdexp(-eye*kmode*th)/nd
                    end do
                    call Z_PLOT(zgrd_bad(inum), 1, options, 31)
               		xgrd_bad(inum) = dreal(zgrd_bad(inum))
					ygrd_bad(inum) = dimag(zgrd_bad(inum)) 
            end do
		!	call FDIFFF(zgrd_bad(inum - ntheta + 1), dzgrd_bad(inum - ntheta +1), ntheta, wsave)
         
		!	do i = 1, ntheta
		!		print *, dzgrd_bad(inum - ntheta +i)
        !    	if (kbod.ne.0) then
        !       		dzgrd_bad(inum - ntheta + i) = -dzgrd_bad(inum - ntheta + i)
        !    	end if
        ! 	end do

       end do
         istart = istart + nd
      end do

	
      open(unit = 32, file = 'mat_plots/xgrid_bad_plot.m')
      open(unit = 33, file = 'mat_plots/ygrid_bad_plot.m')

   
      call X_DUMP(xgrd_bad,(k-k0+1)*nr*ntheta, 32)
      call X_DUMP(ygrd_bad, (k-k0+1)*nr*ntheta,33)

  
      close(32)
      close(33)
 



end subroutine BUILD_CLOSEEVAL_GRID

!----------------------------------------------------------------------

subroutine GET_NEAR_POINTS()

!temporarily assuming equal sized boxes.
! marks points within ig*boxradius of the center of ibox.	

	implicit none

! local variables
	integer:: fac, ibox, jpoint, &
				i, fpoint, istart, kbod, icount
	real(kind=8):: box_rad
	complex(kind=8):: z1, z2 



	fac = ibeta*nd/nb
	icount = 0
	do kbod = k0, k
		do ibox = 1, nb
			z1 = z0_box(kbod- k0 +1, ibox)
			fpoint = (kbod - k0)*ndres +  (ibox - 1)*fac + 1
			box_rad = cdabs(z1 - z_res(fpoint)) 
			istart = 1
			do jpoint = 1,(k - k0 + 1)*ndres			
				z2 = z_res(jpoint)
				if(cdabs(z1 - z2) .le. ig*box_rad) then
					neigh_boxes(kbod + 1, ibox, istart) = jpoint
					istart = istart + 1
					!			print *, kbod, ibox, neigh_boxes(kbod + 1, ibox, istart - 1)
				end if
			end do
			!print *, "istart = ", istart
			if(istart .ne. 2049 ) then
						icount = icount + 1
			end if

			n_neigh(kbod - k0 + 1, ibox) = istart - 1
		end do
	end do
	print *, "Number of boxes which are not neighbours with everyone ", icount
end subroutine GET_NEAR_POINTS

!------------------------------------------------------------------

subroutine RESAMPLE_DOMAIN ()
!
! Resample the boundary points to ibeta*nd points per curve
!
! Input
!    z:  boundary of domain
!   dz:
! Output
!    z_res:  
!
! local variables
   implicit none
   integer :: i, kbod, istart, istartr
   complex(kind=8) :: work(3*nd+3*ibeta*nd+20)
   character(32) :: options
!
! open unit for matlab plotting
      open(unit = 31, file = 'mat_plots/geo_resample.m')
      options = '''g'',''LineWidth'',1'
   
      istart = 0
      istartr = 0
      do kbod = k0, k
         call FINTERC (z(istart+1), z_res(istartr+1), nd, ndres, work)
         call FINTERC (dz(istart+1), dz_res(istartr+1), nd, ndres, work)
         call Z_PLOT(z_res(istartr+1), nd*ibeta, options, 31)
         istart = istart + nd
         istartr = istartr + ndres
      end do
      close(31)
         
end subroutine RESAMPLE_DOMAIN

!----------------------------------------------------------------------

subroutine BUILD_GRID(i_grd, x_grd, y_grd)
!
! This subroutine embeds the domain into a rectangular grid. The coordinates
! of grid points are (x_grd(i,j), y_grd(i,j))
! i_grd(i,j) = 1 if the point is inside the domain
!            = 0 if the point is outside the domain
!
   implicit none
   integer :: i_grd(nx,ny)
   real(kind=8) :: x_grd(nx,ny), y_grd(nx,ny)
   real(kind=8) :: xgrd_good(nx*ny), ygrd_good(nx*ny)
!
! local variables
   integer :: i, j, istart, jstart, kbod, ix, iy, ipot, i_tmp(nx,ny)
   real(kind=8) :: hx, hy, eps, ds_max, xleft, xright, ybot, ytop, x, y, tol

!
! FMM work arrays
   integer :: iprec, ifcharge, ifdipole, ifpot, ifgrad, ifhess, ntarget, &
              ifpottarg, ifgradtarg, ifhesstarg, ier
   real(kind=8) :: source(2,2*nbk), dipvec(2,2*nbk), target(2,nx*ny)
   complex(kind=8) :: charge(2*nbk), dipstr(2*nbk), pot(2*nbk), grad(2,2*nbk), &
                      hess(3,2*nbk), pottarg(nx*ny), gradtarg(2,nx*ny), &
                      hesstarg(3, nx*ny)
   
      hx = (xmax - xmin)/dfloat(nx - 1)
      hy = (ymax - ymin)/dfloat(ny - 1)
      call PRIN2(' hx in BUILD_GRID = *', hx, 1)
      call PRIN2(' hy in BUILD_GRID = *', hy, 1)
   
!
!   Define x and y coordinates of grid
      istart = 1
      do i = 1, nx
         x = xmin + (i - 1.d0) * hx
         do j = 1, ny
            y = ymin + (j - 1.d0) * hy
            x_grd(i, j) = x
            y_grd(i, j) = y
            target(1, istart) = x
            target(2, istart) = y
            istart = istart + 1 
         end do
      end do
      ntarget = istart - 1
      
      i_grd = 0
      i_tmp = 0
            
! set parameters for FMM routine 
    
      iprec = 5   ! err < 10^-14
      ifcharge = 0 ! no charges, only dipoles
      ifdipole = 1
      ifpot = 1
      ifgrad = 0
      ifhess = 0
      ifpottarg = 1
      ifgradtarg = 0
      ifhesstarg = 0

!
!   Determine first, the regions which are inside the domain but 
!   "close" to a contour
!
      istart = 0
      tol = 0.4d0
      n_close = 0.   ! used to count how many points close to each contour
      do kbod = k0, k
!   Assemble arrays for FMM call
         do i = 1, nd
            source(1, i) = dreal(z(istart+i))
            source(2, i) = dimag(z(istart+i))
            dipvec(1,i) = dreal(-eye*dz(istart+i))/ds_dth(istart+i)
            dipvec(2,i) = dimag(-eye*dz(istart+i))/ds_dth(istart+i)
            charge(i) = 0.d0
            dipstr(i) = h*1.d0*ds_dth(istart+i)/(2.d0*pi)

!
!   bad geometry curve
            source(1, i+nd) = dreal(z_bad(istart+i))
            source(2, i+nd) = dimag(z_bad(istart+i))
            dipvec(1,i+nd) = dreal(-eye*dz_bad(istart+i))/ds_bad(istart+i)
            dipvec(2,i+nd) = dimag(-eye*dz_bad(istart+i))/ds_bad(istart+i)
            charge(i+nd) = 0.d0
            dipstr(i+nd) = h*1.d0*ds_bad(istart+i)/(2.d0*pi)
         end do

      
! call FMM

         call PRINI(0, 13)
         call lfmm2dparttarg(ier, iprec, 2*nd, source, ifcharge, charge, &
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
         jstart = 1
         do i = 1, nx
            do j = 1, ny
               ipot = DNINT(dreal(pottarg(jstart)))
               if ((dabs(dreal(pottarg(jstart))-1.d0) .lt. tol) .and. &
                   (kbod.eq.0)) then
                  i_grd(i,j) = 1
                  n_close(1) = n_close(1) + 1
               elseif ((dabs(dreal(pottarg(jstart))+1.d0).lt.tol) .and. &
                       (kbod .ge. 1)) then
                  i_grd(i,j) = -kbod
                  n_close(kbod-k0+1) = n_close(kbod-k0+1) + 1
               end if
               jstart = jstart + 1 
            end do
         end do
         istart = istart + nd
      end do
      open(unit = 31, file = 'mat_plots/igrid_c.m')
      call PRINF ('n_close = *', n_close, k-k0+1)
      call INT_GRID_DUMP(i_grd, 31)

      close(31)
      
! now see which points are well inside the domain, where it is safe to use 
! the trapezoid rule

      istart = 0
      do kbod = k0, k
!   Assemble arrays for FMM call
         do i = 1, nd
            source(1,istart+i) = dreal(z(istart+i))
            source(2,istart+i) = dimag(z(istart+i))
            dipvec(1,istart+i) = dreal(-eye*dz(istart+i))/ds_dth(istart+i)
            dipvec(2,istart+i) = dimag(-eye*dz(istart+i))/ds_dth(istart+i)
            charge(istart+i) = 0.d0
            dipstr(istart+i) = h*1.d0*ds_dth(istart+i)/(2.d0*pi)
         end do
         istart = istart + nd
      end do
      
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
      jstart = 1
	  istart = 1
      tol = 0.4d0
      do i = 1, nx
         do j = 1, ny
            if ((dabs(dreal(pottarg(jstart))-1.d0) .lt. tol).and. &
                (i_grd(i,j) .eq. 0) )then
               i_grd(i,j) = 2
			   xgrd_good(istart) = x_grd(i, j)
			   ygrd_good(istart) = y_grd(i, j)
			   istart = istart + 1
            end if
            jstart = jstart + 1 
         end do
      end do
      
      call ORGANIZE_CLOSE_POINTS (i_grd)
      
      open(unit = 31, file = 'mat_plots/igrid.m')
      open(unit = 32, file = 'mat_plots/xgrid.m')
      open(unit = 33, file = 'mat_plots/ygrid.m')

      open(unit = 34, file = 'mat_plots/xgrd_plot.m')
      open(unit = 35, file = 'mat_plots/ygrd_plot.m')


      call INT_GRID_DUMP(i_grd, 31)
      call REAL_GRID_DUMP(x_grd, 32)
      call REAL_GRID_DUMP(y_grd, 33)
	  call X_DUMP(xgrd_good, istart -1, 34)
	  call X_DUMP(ygrd_good, istart - 1, 35)

      close(31)
      close(32)
      close(33)
	  close(34)
	  close(35)
 
      print *, 'SUCCESSFULLY BUILT GRID'
      
end subroutine BUILD_GRID

!----------------------------------------------------------------------

subroutine POPULATE_INEAR_BOX()

	!local variables
	implicit none
	integer :: nb,ibox, kbod, ipoint, istart, inum
		
	nb = nd/5

	do ibox =  1, nb
		istart = 1
		do kbod = k0, k
			do ipoint = 1, ndres
				inum = kbod*ndres + ipoint
				if(cdabs(z_res(inum) - z0_box(kbod, ibox)).le.g) then
					inear_box(ibox, istart) = inum
					istart = istart + 1
				end if
			end do   
		end do
	end do

end subroutine POPULATE_INEAR_BOX


!------------------------------------------------------------------------

subroutine ORGANIZE_CLOSE_POINTS (i_grd)
!
! This subroutine constructs the arrays i_close and j_close
   implicit none
   integer :: i_grd(nx,ny)
!
! local variables
   integer :: i, j, k_count(kmax), kbod, indx
   
   k_count = 0
   
   ic_pnt(1) = 1
   do kbod = k0+1, k+1
      ic_pnt(kbod-k0+1) = ic_pnt(kbod-k0) + n_close(kbod-k0)
   end do
!!!   call PRINF ('pointer for close points = *', ic_pnt, k-k0+1)
   
   do i = 1, nx
      do j = 1, ny
         if ((i_grd(i,j) .ne. 0) .and. (i_grd(i,j) .ne. 2)) then
            if (i_grd(i,j) .eq. 1) then
               kbod = 0
               indx = ic_pnt(1) + k_count(1)
               k_count(1) = k_count(1) + 1
             else
               kbod = -i_grd(i,j)
               indx = ic_pnt(kbod-k0+1) + k_count(kbod-k0+1)
               k_count(kbod-k0+1) = k_count(kbod-k0+1) + 1
            end if
            i_close(indx) = i
            j_close(indx) = j
         end if 
      end do
   end do
   
!!!   do kbod = k0, k
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
   
end subroutine ORGANIZE_CLOSE_POINTS


!----------------------------------------------------------------------

subroutine BUILD_GRID_OLD(i_grd, x_grd, y_grd)
!
! This subroutine embeds the domain into a rectangular grid. The coordinates
! of grid points are (x_grd(i,j), y_grd(i,j))
! i_grd(i,j) = 1 if the point is inside the domain
!            = 0 if the point is outside the domain
!
   implicit none
   integer :: i_grd(nx,ny)
   real(kind=8) :: x_grd(nx,ny), y_grd(nx,ny)
!
! local variables
   integer :: i, j, istart, kbod, ix, iy, ileft, iright, jbot, jtop
   real(kind=8) :: hx, hy, eps, ds_max, xleft, xright, ybot, ytop, x, y

!
! FMM work arrays
   integer :: iprec, ifcharge, ifdipole, ifpot, ifgrad, ifhess, ntarget, &
              ifpottarg, ifgradtarg, ifhesstarg, ier
   real(kind=8) :: source(2,nbk), dipvec(2,nbk), target(2, nx*ny)
   complex(kind=8) :: charge(nbk), dipstr(nbk), pot(nbk), grad(2,nbk), &
                      hess(3,nbk), pottarg(nx*ny), gradtarg(2, nx*ny), &
                      hesstarg(3, nx*ny)
   
      hx = (xmax - xmin)/dfloat(nx - 1)
      hy = (ymax - ymin)/dfloat(ny - 1)
      call PRIN2(' hx in BUILD_GRID = *', hx, 1)
      call PRIN2(' hy in BUILD_GRID = *', hy, 1)
   
      i_grd = 0
      
!
!   Define x and y coordinates of grid
      istart = 1
      do i = 1, nx
         x = xmin + (i - 1.d0) * hx
         do j = 1, ny
            y = ymin + (j - 1.d0) * hy
            x_grd(i, j) = x
            y_grd(i, j) = y
            target(1, istart) = x
            target(2, istart) = y
            istart = istart + 1 
         end do
      end do
      ntarget = istart - 1
!
!   Determine which grid points are too close to the boundary
!   i_grd(i,j) = -1 if grid point is within eps of boundary point
      istart = 0
      ds_max = 0.d0
      do kbod = k0, k
         ds_max = MAXVAL(ds_dth((kbod-k0)*nd+1:(kbod-k0+1)*nd+1))*h
         call PRINF(' kbod = *', kbod, 1)
         call PRIN2('   ds_max = *', ds_max, 1)
         do i = 1,nd
            eps = 5.d0*ds_max
            xleft = dreal(z(istart+i)) - eps
            xright = dreal(z(istart+i)) + eps
            ybot = dimag(z(istart+i)) - eps
            ytop = dimag(z(istart+i)) + eps
            ileft = IDINT((xleft - xmin)/hx + 1.d0)
               if (ileft.lt.1) ileft = 1
            iright = DNINT((xright - xmin)/hx + 1.d0)
               if (iright.gt.nx) iright = nx
            jbot = IDINT((ybot - ymin)/hy + 1.d0)
               if (jbot.lt.1) jbot = 1
            jtop = DNINT((ytop - ymin)/hy + 1.d0)
               if (jtop.gt.ny) jtop = ny
            do ix = ileft,iright
               do iy = jbot,jtop
                  i_grd(ix,iy) = -1
               end do
            end do
         end do
         istart = istart + nd
      end do
      
!
!   Assemble arrays for FMM call
      do i = 1, nbk
         source(1, i) = dreal(z(i))
         source(2, i) = dimag(z(i))
         dipvec(1,i) = dreal(-eye*dz(i))/ds_dth(i)
         dipvec(2,i) = dimag(-eye*dz(i))/ds_dth(i)
         charge(i) = 0.d0
         dipstr(i) = h*1.d0*ds_dth(i)/(2.d0*pi)
      end do

! set parameters for FMM routine 
    
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
      istart = 1
      do i = 1, nx
         do j = 1, ny
            if (i_grd(i,j) .ge. 0) then  
               i_grd(i, j) = k0 + DNINT(dreal(pottarg(istart)))
             else
               i_grd(i, j) = 0
            end if
            istart = istart + 1 
         end do
      end do
      
      open(unit = 31, file = 'mat_plots/igrid.m')
      open(unit = 32, file = 'mat_plots/xgrid.m')
      open(unit = 33, file = 'mat_plots/ygrid.m')

      call INT_GRID_DUMP(i_grd, 31)
      call REAL_GRID_DUMP(x_grd, 32)
      call REAL_GRID_DUMP(y_grd, 33)

      close(31)
      close(32)
      close(33)
 
   
end subroutine BUILD_GRID_OLD

!------------------
! PLOTTING ROUTINES
!------------------
!----------------------------------------------------------------------

subroutine X_DUMP(x, n, iw)
! This subroutine writes out the vector x 
! 
   implicit none
   integer, intent(in) :: n, iw
   real(kind=8), intent(in) :: x(n)

! local variables
   integer :: i
   
      

      write(iw, *) 'sol = ['
      write(iw, '(D15.6)') (x(i), i = 1, n)
      write(iw, *) '];'
      
end subroutine X_DUMP



!----------------------------------------------------------------------


subroutine PARAM_DUMP(x, y, z, iw)
! This subroutine writes out the vector x 
! 
   implicit none
   integer, intent(in) ::  iw
   integer, intent(in) :: x,y,z

! local variables
   integer :: i
   
      

      write(iw, *) 'K = ', x
	  write(iw, *) 'NR = ', y
      write(iw, *) 'NT = ', z
	  
      
end subroutine PARAM_DUMP



!----------------------------------------------------------------------


subroutine Z_PLOT(z, n, options, iw)
!
! This subroutine writes out the curve specified by its nodes z in 
! complex plane for matlab plotting. It assumes the curve is not closed,  
! i.e. z(1) <> z(n)
! Inputs:
!    z:  points on curve
!    n: number of points on curve
!    options: a string for plotting options for matlab command "plot"
!      this string should look like '''string'''
!    iw:  the fortran unit number on which the output data set is written
! Outputs:
!    none
! 
   implicit none
   integer, intent(in) :: n, iw
   complex(kind=8), intent(in) :: z(n)
   character(32), intent(in) :: options
! local variables
   character(32) :: plot_str, end_bracket
   integer :: i
   
      plot_str = 'plot(x(:,1),x(:,2),'
      end_bracket = ')'

      write(iw, *) 'x = ['
      write(iw, '(2(D15.6))') (z(i), i = 1, n)
      write(iw, '(2(D15.6))') z(1)
      write(iw, *) '];'
      write(iw, *) trim(plot_str) // trim(options) // &
                   trim(end_bracket)
      write(iw, *) 'hold on'
      
end subroutine Z_PLOT

!----------------------------------------------------------------------

subroutine XY_PLOT(x, y, n, options, iw)
!
! This subroutine writes out the curve specified by its nodes x, y
! for matlab plotting. It assumes the curve is not closed, i.e. 
! (x(1), y(1)) <> (x(n), y(n))
! Inputs:
!    (x, y):  points on curve
!    n: number of points on curve
!    options: a string for plotting options for matlab command "plot"
!    iw:  the fortran unit number on which the output data set is written
! Outputs:
!    none
! 
   implicit none
   integer, intent(in) :: n, iw
   real(kind=8), intent(in) :: x(n), y(n)
   character(32), intent(in) :: options
! local variables
   character(32) :: plot_str, end_bracket
   integer :: i
   
      plot_str = 'plot(x(:,1),x(:,2),'
      end_bracket = ')'

      write(iw, *) 'x = ['
      write(iw, '(2(D15.6))') (x(i), y(i), i = 1, n)
      write(iw, *) '];'
      write(iw, *) trim(plot_str) // trim(options) // &
                   trim(end_bracket)
      write(iw, *) 'hold on'
      
end subroutine XY_PLOT


!----------------------------------------------------------------------

subroutine Y_PLOT(y, n, options, iw)
!
! This subroutine writes out the curve specified by its nodes x, y
! for matlab plotting. It assumes the curve is not closed, i.e. 
! (x(1), y(1)) <> (x(n), y(n))
! Inputs:
!    (x, y):  points on curve
!    n: number of points on curve
!    options: a string for plotting options for matlab command "plot"
!    iw:  the fortran unit number on which the output data set is written
! Outputs:
!    none
! 
   implicit none
   integer, intent(in) :: n, iw
   real(kind=8), intent(in) :: y(n)
   character(32), intent(in) :: options
! local variables
   character(32) :: plot_str, end_bracket
   integer :: i
   
      plot_str = 'plot(x,'
      end_bracket = ')'

      write(iw, *) 'x = ['
      write(iw, '(1(D15.6))') (y(i), i = 1, n)
      write(iw, *) '];'
      write(iw, *) trim(plot_str) // trim(options) // &
                   trim(end_bracket)
      write(iw, *) 'hold on'
      
end subroutine Y_PLOT

!----------------------------------------------------------------------

subroutine INT_GRID_DUMP(igrid, iunit)

!
   implicit none
   integer :: iunit, igrid(nx,ny)
!
! local variables
   integer i, j
   
      write(iunit, *) 'NX = ', ny
      write(iunit, *) 'NY = ', nx
      write(iunit, *) 'a = zeros(NX,NY);'
      write(iunit, *) 'sol = ['
      
      do i = 1, nx
         do j = 1, ny
            write(iunit, '(i4,$)') (igrid(i,j))
            write(iunit, '(a)') ''
         end do
      end do
      
      write (iunit, *) '];'
      write (iunit, *) 'a(:) = sol(:);'

end subroutine INT_GRID_DUMP

!----------------------------------------------------------------------

subroutine REAL_GRID_DUMP(rgrid, iunit)

!
   implicit none
   integer :: iunit
   real(kind=8) :: rgrid(nx,ny)
!
! local variables
   integer i, j
   
      write(iunit, *) 'NX = ', ny
      write(iunit, *) 'NY = ', nx
      write(iunit, *) 'a = zeros(NX,NY);'
      write(iunit, *) 'sol = ['
      
      do i = 1, nx
         do j = 1, ny
            write(iunit, '(e20.13,$)') (rgrid(i,j))
            write(iunit, '(a)') ''
         end do
      end do
      
      write (iunit, *) '];'
      write (iunit, *) 'a(:) = sol(:);'

end subroutine REAL_GRID_DUMP

end module geometry_mod
