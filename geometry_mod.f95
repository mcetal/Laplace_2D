!
! comman variables etc for laplace_2d
!
module geometry_mod

   implicit none
   
   save

   integer, parameter :: kmax = 10, npmax = 2048, nmax = kmax*npmax
!
! Number of holes, system size
   integer :: k0, k, nd, nbk
!
! Geometry of Domain
   real(kind=8) :: ak(kmax), bk(kmax)
   integer :: ncyc(kmax)
   complex(kind=8) :: zk(kmax)
!
! Constants and mesh spacing
   real(kind=8) :: pi, h
   complex(kind=8) :: eye
!
! Discretization of domain
   real(kind=8) :: x(nmax), y(nmax), ds_dth(nmax), kappa(nmax), &
                   theta(nmax)
   complex(kind=16) :: z(nmax), dz(nmax)
      
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
   real(kind=8) xp, yp, xdot, ydot, xddot, yddot, th
   
!
! construct theta
      do i = 1, nd
         theta(i) = h*(i-1.d0)
      end do
      
      call PRIN2('in build domain, h = *', h, 1)
!
! construct holes
      istart = 0
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
         end do   
         call PRIN2('z = *', z((kbod-k0)*nd+1), nd)
         call PRIN2('kappa = *', kappa((kbod-k0)*nd+1), nd)
         istart = istart + nd  
      end do

end subroutine BUILD_DOMAIN

end module geometry_mod
