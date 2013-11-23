
!----------------------------------------------------------------------

subroutine SOLVE (maxl, rhs, lrwork, liwork, dirichlet, soln, mu, A_log, &
                  gmwork, igwork, schur, ipvtbf, bnew)

!
! Solves integral equation using GMRES
! Inputs:
!   maxl: max number of iterations before restart
!   rhs: rhs of system
!   lrwork: size of real workspace for GMRES
!   liwork: size of integer workspace for GMRES
!   dirichlet: logical, true if dirichlet bvp
! Returns:
!   soln: solution
!   mu: density of integral equation
!   A_log: log sources
!   gmwork: real workspace for GMRES
!   igwork: integer workspace for GMRES
!   schur: preconditioner Schur complement
!   ipvtbf: pivots for precondition
!   bnew: rhs for preconditioner

   use geometry_mod, only: k0, k, kmax, nd, nbk, bounded
   implicit none
   integer, intent(in) :: maxl, lrwork, liwork
   real(kind=8), intent(in) :: rhs(nbk)
   logical, intent(in) :: dirichlet
   integer, intent(out) :: igwork(liwork), ipvtbf(kmax)
   real(kind=8), intent(out) :: soln(*), mu(nbk), A_log(k), gmwork(lrwork), &
                                schur(kmax * kmax), bnew(kmax)

! Local variables
   external MATVEC_DIR, MSOLVE
   real(kind=4) t0, t1, tsec, timep(2), etime
   integer :: i, itol, isym, norder, ierr, iw, nelt, ia, ja, iter, itmax, &
              kbod
   real(kind=8) :: tol, a, err, sb, sx, rw
   
!
!  parameters for DGMRES - see dgmres.f to see what they mean
      itol = 0
      tol = 1.0d-12
      isym = 0

!  Restart flag, will restart after iwork(1) iterations. This value 
!  should be less than or equal to maxl
      igwork(1) = maxl      

      do i = 2, liwork
         igwork(i) = 0
      end do
      norder = nbk + k

!  Preconditioner flag 
!  igwork(4) = 0 - no preconditioner
!  igwork(4) < 0 - preconditioner on the left
!  igwork(4) > 0 - preconditioner on the right
      igwork(4) = -1

!  provide initial guess soln
      do i = 1, norder
         soln(i) = rhs(i)
      end do

!  factor preconditioner
      if (igwork(4) < 0) then  
         t0 = etime(timep)
         if ((bounded) .and. (dirichlet)) then  
            call SCHUR_FACTOR_DIR_BNDED(schur, bnew, ipvtbf)
         end if
         stop
         t1 = etime(timep)
         call PRIN2('Time in factoring preconditioner = *', t1 - t0, 1)
      end if
      
      t0 = etime(timep)
      call PRINI(0, 13)
      
      if (dirichlet) then   
         call DGMRES(norder, rhs, soln, nelt, ia, ja, a, isym, MATVEC_DIR, &
                     MSOLVE, itol, tol, itmax, iter, err,ierr, 6, sb, sx, &
                     gmwork, lrwork, igwork, liwork, rw, iw)
      end if
      
      call PRINI(6, 13)
      call PRINF('# GMRES ITERATIONS = *',iter,1)
      if (ierr.gt.2) then
         call PRINF('SOMETHING WRONG IN GMRES, IERR = *',ierr,1)
         call PRINF('iwork = *',igwork,10)
         stop
        elseif (ierr.ge.0) then
         t1 = etime(timep)
         tsec = t1 - t0
         call PRIN2 ('TIME TAKEN IN GMRES = *', tsec, 1)
      end if

!  unpack RHS into U and A_log
      do i = 1,nbk
         mu(i) = soln(i)
      end do 
      do kbod = 1, k
         A_log(kbod) = soln(nbk + kbod)
      end do
!      call PRIN2(' mu = *', mu, nbk) 
      call PRIN2('A_log = *', A_log, k)

end subroutine SOLVE

!----------------------------------------------------------------------

subroutine MATVEC_DIR(N, XX, YY, NELT, IA, JA, A, ISYM)

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
   use geometry_mod, only: nmax, kmax
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
!
! local arrays
   real(kind=8) :: A_log(kmax)

      t0 = etime(timep)
      call FASMVP_DIR(N, xx, yy, A_log, source, dipvec, charge, &
                      dipstr, pot, grad, hess)
      
      t1 = etime(timep)

!      WRITE(13,*) 'TIME IN SECONDS FOR MATVEC = ', t1 - t0
!      WRITE(6,*) 'TIME IN SECONDS FOR MATVEC = ', t1 - t0

end subroutine MATVEC_DIR

!----------------------------------------------------------------------

subroutine MSOLVE(N, R, U, NELT, IA, JA, A, ISYM, RWORK, IWORK)

!
!  Another routine required by DGMRES. It allows the use of a
!  preconditioner.
!
   implicit double precision (A-H,O-Z)


end subroutine MSOLVE

!----------------------------------------------------------------------

subroutine SCHUR_FACTOR_DIR_BNDED(schur, wb, ipvtbf)

! Constructs and factors the Schur complement of the preconditioner 
! matrix for the Dirichlet BVP in bounded domains
! The discrete system takes the form Mx = b where
!     
!           M = |Mp | C |   Mp is N x N (N = total # bdry points)
!               |---|---|   C is N x k, R is k x N and L is k x k.
!               | R | L |
!     Mp is of the form I + Q where Q is an integral operator 
!     with continuous kernel. A good preconditioner, therefore,
!     is
!     
!           B = | I | C |   I is N x N (N = total # bdry points)
!               |---|---|   C is N x k, R is k x N and L is k x k.
!               | R | L |
!
!     The vector U is ordered as 1) dipole densities 2) coeffs Ak.
!
! Returns: 
!   Schur: contains the LU factors of the Schur complement
!   ipvtbf: contains pivoting information from LINPACK
!   wb: this is a work array for DGECO

   use geometry_mod, only: k0, k, kmax, nd, nbk, z, zk, ds_dth
   implicit none
   integer, intent(out) :: ipvtbf(k)
   real(kind=8) :: schur(k, k), wb(k)

! local variables
   integer :: ibod, kbod, i, istart, info, iwork(k)
   real(kind=8) :: sum1, rcond, anorm, work(4*k)
   character*1 :: norm
   
      print *, '** In Preconditioner factoring routine  **'

      istart = nd
      do ibod = 1, k
         do kbod = 1, k
            sum1 = 0.d0
            do i = 1, nd
               sum1 = sum1 &
                      + dlog(cdabs(z(istart+i) - zk(kbod + 1 - k0)))
            end do
            schur(ibod, kbod) = sum1
         end do
         istart = istart + nd
      end do

      call DGETRF(k, k, schur, k, ipvtbf, info) 
      if (info < 0) then
         call PRINF(' DGETRF ERROR: illegal value in argument *', info, 1)
         stop
      else if (info > 0) then
         call PRINF(' DGETRF ERROR: matrix is singular, *', info, 1)
         stop
      end if
      
! Compute the column sum for 1 norm
      anorm = 0.d0
      do kbod = 1, k
         sum1 = 0.d0
         do ibod = 1, k
	    sum1 = sum1 + dabs(schur(ibod, kbod))
	 end do
	 anorm = max(anorm, sum1)
      end do
      call PRIN2(' 1 Norm of Schur complement = *', anorm, 1)
      
! estimate condition number of schur complement
      norm = '1'
      call DGECON(norm, k, schur, k, anorm, rcond, work, iwork, info) 
      if (info < 0) then
         call PRINF(' DGECON ERROR: illegal value in argument *', info, 1)
         stop
      end if
      call prin2(' Condition number of Schur Complement = *', 1.d0/rcond, 1)

!!!      print *, '   After factorization'
!!!      do kbod = 1, k
!!!	 call prinf(' column *', kbod, 1)
!!!	 call prin2(' SCHUR = *',schur(1,kbod), k)
!!!      end do
!!!      call PRINf(' ipvtbf = *', ipvtbf, k)
!!!      call PRIN2(' wb = *', wb, k)

end subroutine SCHUR_FACTOR_DIR_BNDED

!----------------------------------------------------------------------

subroutine FASMVP_DIR(n, u, w, A_log, source, dipvec, charge, dipstr,  &
                      pot, grad, hess)

!
! Calculates the matrix vector product
!     u is the current guess for the density
!     w is (0.5I + K) u + sum A_log 
! 
!
   use geometry_mod
   implicit none
   integer, intent(in) :: n
   real(kind=8), intent(in) :: u(n)
   real(kind=8), intent(out) :: A_log(k), w(n), source(2,n), dipvec(2,n)
   complex(kind=8), intent(out) :: charge(n), dipstr(n), pot(n), &
                                   grad(2,n), hess(3,n)
!
! fmm 
   integer :: ier, iprec, nsource, ifcharge, ifdipole, ifpot, ifgrad, ifhess
!
! local work variables
   integer :: i, kbod, istart
   real(kind=8) self, cauchy, far_field

! set density and source points for fmm call
      
      do i = 1, nbk
         charge(i) = 0.d0
         dipstr(i) = h*u(i)*ds_dth(i)/(2.d0*pi)
         source(1,i) = dreal(z(i))
         source(2,i) = dimag(z(i))
         dipvec(1,i) = dreal(-eye*dz(i))/ds_dth(i)
         dipvec(2,i) = dimag(-eye*dz(i))/ds_dth(i)
      end do
      
! get current values of log sources
      do kbod = 1, k
         A_log(kbod) = u(nbk + kbod)
      end do

! calculate far field if unbounded problem
      far_field = 0.d0
      if (.not. bounded) then
         do i = 1, nbk
            far_field = far_field + u(i) * ds_dth(i)
         end do
         far_field = far_field/(2.d0*pi)
      end if 
         
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
         call prin2(' cauchy = *', cauchy, 1)
         w(i) = 0.5d0*u(i) + cauchy + far_field
         
!      add on log sources
         do kbod = 1, k
            w(i) = w(i)  &
                    + A_log(kbod)*dlog(cdabs(z(i) - zk(kbod + 1 - k0)))
         end do
      end do
      
! Constraints for multiply-connected domain
      if (bounded) then
         istart = nd
         do kbod = 1, k
            do i = 1, nd
               w(nbk + kbod) = w(nbk + kbod) + u(istart+i)
            end do
            istart = istart + nd
         end do
       else
         istart = nd
         w(nbk + 1) = 0.d0
         do kbod = 2, k
            w(nbk + 1) = w(nbk + 1) + A_log(kbod)
         end do
         do kbod = 2, k
            do i = 1, nd
               w(nbk + kbod) = w(nbk + kbod) + u(istart+i)
            end do
            istart = istart + nd
         end do
      end if 

end subroutine FASMVP_DIR
