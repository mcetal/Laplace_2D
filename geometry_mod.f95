!
! comman variables etc for 2-D integral equation solvers
! routines to describe the boundary of the domain and the grid
! make sure a directory "mat_plots/" exists
!
module geometry_mod

   implicit none
   
   save

   integer, parameter :: kmax = 10, npmax = 2048, nmax = kmax*npmax
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
   complex(kind=8) :: z_bad(nmax), dz_bad(nmax)   
!
! Grid variables
   integer :: nx, ny, i_grd(ngrd_max)
   real(kind=8) :: x_grd(ngrd_max), y_grd(ngrd_max)
      
contains
   
!----------------------------------------------------------------------

subroutine CURVE_PARAM(th, kbod, xp, yp, xdot, ydot, xddot, yddot)
!
!  The initial curve parameterization is given in term of th
!  Inputs:
!    th:  parameter value
!    kbod: the particular body index point is on
!  Returns:
!    xp, yp:  coordinates of point
!    xdot, ydot:  derivative of position wrt parameter
!    xddot, yddot:  2nd derivative of position wrt parameter
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
            rddot =  -ai*bi*(2.d0*a2*b2*cs**2*sn**2                &
                           +a2**2*cs**4 + a2*b2-a2**2+b2**2*cs**4  &
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
            den = (1.d0-2.d0*anu*dcos(2.d0*th)+anu**2)  &
                 *dsqrt(1.d0+anu**2)
            ddot = 4.d0*anu*dsin(2.d0*th)*dsqrt(1.d0+anu**2)
            dddot = 8.d0*anu*dcos(2.d0*th)*dsqrt(1.d0+anu**2)
            xp = (1.d0-anu**2)*(1.d0-anu)*R*dsqrt(2.d0)*dcos(th)/den
            yp = (1.d0-anu**2)*(1.d0+anu)*R*dsqrt(2.d0)*dsin(th)/den
            xdot = -(1.d0-anu**2)*(1.d0-anu)*R*dsqrt(2.d0)*dsin(th)/den  &
                   -xp*ddot/den            
            xddot = -xp   &
                    +(1.d0-anu**2)*(1.d0-anu)*R*dsqrt(2.d0)  &
                                  *dsin(th)*ddot/den**2   &
                    - xdot*ddot/den - xp*dddot/den    &
                    + xp*ddot**2/den**2
            ydot = (1.d0-anu**2)*(1.d0+anu)*R*dsqrt(2.d0)*dcos(th)/den  &
                   -yp*ddot/den
            yddot = -yp  &
                    -(1.d0-anu**2)*(1.d0+anu)*R*dsqrt(2.d0)  &
                                  *dcos(th)*ddot/den**2  &
                    - ydot*ddot/den - yp*dddot/den   &
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
!
! local variables
   implicit none
   integer :: i, kbod, istart, kmode
   real(kind=8) alpha_bad
   character(32) :: options
   
!
! FFT work arrays
   real(kind=8) wsave(4*nd + 15)
   complex(kind=8) :: zf(nd)
   
!
! open unit for matlab plotting
      open(unit = 31, file = 'mat_plots/geo_bad_contours.m')
      options = '''r'',''LineWidth'',1'

      call DCFFTI (nd, wsave)
      
      alpha_bad = 5.d0*h
      call PRIN2 (' alpha_bad = *', alpha_bad, 1)
      istart = 0
      do kbod = k0, k
         do i = 1, nd
            zf(i) = z(istart+i)
         end do
         call DCFFTF (nd, zf, wsave)
         zf(1) = zf(1)/nd
         do kmode = 1, nd/2 - 1
            if (kbod.eq.0) then 
               zf(kmode+1) = zf(kmode+1)*dexp(-kmode*alpha_bad)/nd
               zf(nd-kmode+1) = zf(kmode+1)*dexp(kmode*alpha_bad)/nd
             else
               zf(kmode+1) = zf(kmode+1)*dexp(kmode*alpha_bad)/nd
               zf(nd-kmode+1) = zf(kmode+1)*dexp(-kmode*alpha_bad)/nd
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
      end do
      close(31)

end subroutine BAD_DOMAIN_BNDRY

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
               elseif ((dabs(dreal(pottarg(jstart))+1.d0).lt.tol) .and. &
                       (kbod .ge. 1)) then
                  i_grd(i,j) = -kbod
               end if
               jstart = jstart + 1 
            end do
         end do
         istart = istart + nd
      end do
      open(unit = 31, file = 'mat_plots/igrid_c.m')

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
      tol = 0.4d0
      do i = 1, nx
         do j = 1, ny
            if ((dabs(dreal(pottarg(jstart))-1.d0) .lt. tol).and. &
                (i_grd(i,j) .eq. 0) )then
               i_grd(i,j) = 2
            end if
            jstart = jstart + 1 
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
 
      print *, 'SUCCESSFULLY BUILT GRID'
      
end subroutine BUILD_GRID

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
      write(iw, '(2(D15.6))') x(1), y(1)
      write(iw, *) '];'
      write(iw, *) trim(plot_str) // trim(options) // &
                   trim(end_bracket)
      write(iw, *) 'hold on'
      
end subroutine XY_PLOT


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
