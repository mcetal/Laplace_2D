
!----------------------------------------------------------------------

subroutine SOLVE (maxl, rhs, lrwork, liwork, soln, mu, A_log, &
                  gmwork, igwork )

!
! Solves integral equation using GMRES
! Inputs:
!   maxl: max number of iterations before restart
!   rhs: rhs of system
!   lrwork: size of real workspace for GMRES
!   liwork: size of integer workspace for GMRES
! Returns:
!   soln: solution
!   mu: density of integral equation
!   A_log: log sources
!   gmwork: real workspace for GMRES
!   igwork: integer workspace for GMRES

   use geometry_mod, only: k0, k, nd, nbk
   implicit none
   integer, intent(in) :: maxl, lrwork, liwork
   real(kind=8), intent(in) :: rhs(nbk)
   integer, intent(out) :: igwork(liwork)
   real(kind=8), intent(out) :: soln(*), mu(nbk), A_log(k), gmwork(lrwork)
   external MATVEC, MSOLVE
   real(kind=4) t0, t1, tsec, timep(2), etime
   integer :: i, itol, isym, norder, ierr, iw, nelt, ia, ja, iter, itmax
   real(kind=8) :: tol, a, err, sb, sx, rw

!
!  parameters for DGMRES - see dgmres.f to see what they mean
      itol = 0
      tol = 1.0d-12
      isym = 0
      do i = 2, liwork
         igwork(i) = 0
      end do
      norder = nbk

!  Preconditioner flag 
      igwork(4) = 0

!  Restart flag, will restart after iwork(5) iterations. This value 
!  should be less than or equal to maxl
      igwork(5) = maxl      

!  provide initial guess soln
      do i = 1, norder
         soln(i) = rhs(i)
      end do

      t0 = etime(timep)
      call PRINI(0, 13)
      call DGMRES(norder, rhs, soln, nelt, ia, ja, a, isym, MATVEC, &
                  MSOLVE, itol, tol, itmax, iter, err,ierr, 6, sb, sx, &
                  gmwork, lrwork, igwork, liwork, rw, iw)
      call PRINI(6, 13)
      call PRINF('  # GMRES ITERATIONS = *',iter,1)
      if (ierr.gt.2) then
         call PRINF('  SOMETHING WRONG IN GMRES, IERR = *',ierr,1)
         call PRINF('  iwork = *',igwork,10)
         stop
        elseif (ierr.ge.0) then
         t1 = etime(timep)
         tsec = t1 - t0
         call PRIN2 (' TIME TAKEN IN GMRES = *', tsec, 1)
      end if

!  unpack RHS into U
      do i = 1,nbk
         mu(i) = soln(i)
      end do 
      call PRIN2(' mu = *', mu, nbk) 

end subroutine SOLVE

!----------------------------------------------------------------------

subroutine MATVEC (N, XX, YY, NELT, IA, JA, A, ISYM)

!
! Required by DGMRES with this precise calling sequence.
! We ignore most of the parameters except N, XX and YY
! and call the fast multipole routine FASMVP to do the actual
! work.
! Inputs:
!   N: number of unknowns
!   XX: current guess to solution
! Returns: 
!   YY: matrix applied to vector XX
! Dummy stuff:
!   IA, JA, A, ISYM
!
   use geometry_mod, only: nmax
   implicit none
   integer, intent(in) :: N, NELT, IA, JA, ISYM
   real(kind=8), intent(in) :: XX(N), A
   real(kind=8), intent(out) :: YY(N)
   real(kind=4) :: timep(2), etime, t0, t1
!
! fmm arrays
   real(kind=8) :: source(2*nmax), dipvec(2*nmax)
   complex(kind=8) :: charge(nmax), dipstr(nmax), pot(nmax), &
                      grad(2*nmax), hess(3*nmax)

      t0 = etime(timep)

      call FASMVP (N, xx, yy, source, dipvec, charge, dipstr, pot, &
                   grad, hess)
      
      t1 = etime(timep)

!      WRITE(13,*) 'TIME IN SECONDS FOR MATVEC = ', t1 - t0
!      WRITE(6,*) 'TIME IN SECONDS FOR MATVEC = ', t1 - t0

end subroutine MATVEC

!----------------------------------------------------------------------

subroutine MSOLVE(N, R, U, NELT, IA, JA, A, ISYM, RWORK, IWORK)

!
!  Another routine required by DGMRES. It allows the use of a
!  preconditioner.
!
   implicit double precision (A-H,O-Z)


end subroutine MSOLVE

!----------------------------------------------------------------------

subroutine FASMVP(n, u, w, source, dipvec, charge, dipstr, pot, grad, &
                  hess)

!
! Calculates the matrix vector product
!     u is the current guess for the density
!     w is (0.5I + K) u
! 
!
   use geometry_mod
   implicit none
   integer, intent(in) :: n
   real(kind=8), intent(in) :: u(n)
   real(kind=8), intent(out) :: w(n), source(2,n), dipvec(2,n)
   complex(kind=8), intent(out) :: charge(n), dipstr(n), pot(n), &
                                   grad(2,n), hess(3,n)
!
! fmm 
   integer :: ier, iprec, nsource, ifcharge, ifdipole, ifpot, ifgrad, ifhess
!
! local work variables
   integer :: i
   real(kind=8) self, cauchy

! set density and source points for fmm call
      
      do i = 1, nbk
         charge(i) = 0.d0
         dipstr(i) = h*u(i)*ds_dth(i)/(2.d0*pi)
         source(1,i) = dreal(z(i))
         source(2,i) = dimag(z(i))
         dipvec(1,i) = dreal(-eye*dz(i))
         dipvec(2,i) = dimag(-eye*dz(i))
      end do

! set parameters for FMM routine DAPIF2
	
      iprec = 5   ! err < 10^-14
      nsource = nbk
      ifcharge = 0 ! no charges, only dipoles
      ifdipole = 1
      ifpot = 1
      ifgrad = 0
      ifhess = 0

! call FMM

      call lfmm2dpartself(ier, iprec, nsource, source, ifcharge, charge, &
                          ifdipole, dipstr, dipvec, ifpot, pot, ifgrad,  &
                          grad, ifhess, hess)
	
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
	  

! Discrete integral operator
      do i = 1, nbk
         self = 0.25d0*h*kappa(i)*ds_dth(i)/pi
         cauchy = self*u(i) + dreal(pot(i))
         w(i) = 0.5d0*u(i) + cauchy
      end do

end subroutine FASMVP
