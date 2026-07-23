!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!

! #define LAPACK

module analysis_stripping_tools
!
! Simulates the inspiral of two stars in a circular orbit caused by gravitational wave
!   radiation.
!   Author: Bernard Field (supervisor: James Wurster & Paul Lasky)
!   Changes for stripping were done by Marat Potashov
!
! :References: e.g. Tong (2015) classical dynamics lecture notes
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - stop_ratio : *ratio of particles crossing CoM to indicate a merger*
!
! :Dependencies: centreofmass, dump_utils, infile_utils, io, physcon, units
!
   implicit none

   !
   ! Runtime parameters
   !

   integer,    public :: nstar(2) = 0 ! give default value in case dump header not read

   real,       public :: com(3)
   real,       public :: vcom(3)

   ! save from previous call (of get_gw_force)
   real, save, public :: evector_old(3)
   real, save, public :: time_old(2)
   real, save, public :: omega_old(3)
   !
   ! subroutines
   !
   public :: get_momentofinertia, calculate_omega, correct_evector, &
      correct_sign_evector, get_internal_quadrupole
   private

contains
!-----------------------------------------------------------------------
!+
! Calculates the moment of inertia
! This is done about the coordinate axes whose origin is at the
!   centre of mass
! Mechanics, Third Edition: Volume 1 (Course of Theoretical Physics)
!   L. Landau, and E. Lifshitz. eq. 32.6
!+
!-----------------------------------------------------------------------
   subroutine get_momentofinertia(xyzh,vxyzu,center_of_mass,vcenter_of_mass,npart,density_cutoff,particlemass,&
      npartused,inertia,principle,evectors,rmax,omega)

      use part, only: rhoh
      use vectorutils, only: cross_product3D

      real,             intent(in)  :: xyzh(:,:)
      real,             intent(in)  :: vxyzu(:,:)
      real,             intent(in)  :: center_of_mass(3)
      real,             intent(in)  :: vcenter_of_mass(3)
      integer,          intent(in)  :: npart
      real,             intent(in)  :: density_cutoff
      real,             intent(in)  :: particlemass
      integer,          intent(out) :: npartused
      real,             intent(out) :: inertia(3,3)
      real,             intent(out) :: principle(3), evectors(3,3)
      real,             intent(out) :: rmax
      real,   optional, intent(out) :: omega(3)

      integer                       :: i
      real                          :: inertia_copy(3,3)
      real                          :: dot_inertia(3,3)
      integer                       :: smallIIndex
      real                          :: smallI
      real                          :: smallIEvector(3)
      real                          :: c(3)
      real                          :: c1(3)
      real                          :: dRdt(3)
! #ifdef LAPACK
      ! real                          :: inertia2(3,3)
! #endif
      real                          :: x,y,z,vx,vy,vz,r2,rmax2

      inertia     = 0.0
      dot_inertia = 0.0
      npartused   = 0
      rmax2       = 0.0
      if (present(omega))&
         omega       = 0.0

!$omp parallel default(none) &
!$omp shared(npart,xyzh,vxyzu,center_of_mass,vcenter_of_mass,particlemass,density_cutoff) &
!$omp private(i,x,y,z,vx,vy,vz,r2) &
!$omp reduction(+:inertia,dot_inertia,npartused) &
!$omp reduction(max:rmax2)
!$omp do
      do i = 1, npart
         if(rhoh(xyzh(4,i),particlemass) > density_cutoff) then
            x = xyzh(1,i) - center_of_mass(1)
            y = xyzh(2,i) - center_of_mass(2)
            z = xyzh(3,i) - center_of_mass(3)
            vx = vxyzu(1,i) - vcenter_of_mass(1)
            vy = vxyzu(2,i) - vcenter_of_mass(2)
            vz = vxyzu(3,i) - vcenter_of_mass(3)
            inertia(1,1) = inertia(1,1) + y**2 + z**2
            inertia(2,2) = inertia(2,2) + x**2 + z**2
            inertia(3,3) = inertia(3,3) + x**2 + y**2
            inertia(1,2) = inertia(1,2) - x*y
            inertia(1,3) = inertia(1,3) - x*z
            inertia(2,3) = inertia(2,3) - y*z
            dot_inertia(1,1) = dot_inertia(1,1) + 2.0*y*vy + 2.0*z*vz
            dot_inertia(2,2) = dot_inertia(2,2) + 2.0*x*vx + 2.0*z*vz
            dot_inertia(3,3) = dot_inertia(3,3) + 2.0*x*vx + 2.0*y*vy
            dot_inertia(1,2) = dot_inertia(1,2) - vx*y - x*vy
            dot_inertia(1,3) = dot_inertia(1,3) - vx*z - x*vz
            dot_inertia(2,3) = dot_inertia(2,3) - vy*z - y*vz
            ! Additional useful values
            npartused    = npartused + 1
            r2           = x*x + y*y + z*z
            rmax2        = max(rmax2, r2)
         endif
      enddo
!$omp enddo
!$omp end parallel
      rmax = sqrt(rmax2)
      !--The symmetric components
      inertia(2,1) = inertia(1,2)
      inertia(3,1) = inertia(1,3)
      inertia(3,2) = inertia(2,3)
      dot_inertia(2,1) = dot_inertia(1,2)
      dot_inertia(3,1) = dot_inertia(1,3)
      dot_inertia(3,2) = dot_inertia(2,3)
      !--Multiply in constant
      inertia      = inertia*particlemass
      dot_inertia  = dot_inertia*particlemass
      inertia_copy = inertia
      !
! #ifdef LAPACK
      ! inertia2 = inertia
! #endif
      !
      !--Find the eigenvectors
      !
#ifndef LAPACK
      !  note: i is a dummy out-integer that we don't care about
      call jacobi(inertia_copy,3,3,principle,evectors,i)
      ! write(*,*) 'Eigenvalues JACOBI:'
      ! do i = 1, 3
      !   write(*,*) i, principle(i)
      ! enddo
      ! write(*,*)
      ! write(*,*) 'Eigenvectors JACOBI:'
      ! do i = 1, 3
      !   write(*,*) i, evectors(:,i)
      ! enddo
      ! write(*,*)
#else
      call eigensystem(inertia_copy,3,principle)
      evectors = inertia_copy
      ! call eigensystem(inertia2,3,principle)
      ! evectors = inertia2

      ! write(*,*) 'Eigenvalues LAPACK:'
      ! do i = 1, 3
      !   write(*,*) i, principle(i)
      ! enddo
      ! write(*,*)
      ! write(*,*) 'Eigenvectors LAPACK:'
      ! do i = 1, 3
      !   write(*,*) i, evectors(:,i)
      ! enddo
      ! write(*,*)
#endif
      !
      if (present(omega)) then
         ! \[
         ! \mathbf{\Omega}_i^\mathrm{orb} =
         !   \sum_{j \neq i}\frac{1}{\lambda_j - \lambda_i}
         !     \left[\Big((\mathbf{e}_j^T \mathbf{\dot{I}} \mathbf{e}_i)\mathbf{e}_j\Big)
         !       \times \mathbf{e}_i \right]
         ! \]
         smallIIndex = minloc(principle, dim=1)
         smallIEvector = evectors(:, smallIIndex)
         smallI = principle(smallIIndex)
         c = matmul(dot_inertia, smallIEvector)
         dRdt = 0.0
         do i = 1, 3
            if(i == smallIIndex) cycle
            c1 = evectors(:, i)
            dRdt = dRdt + (dot_product(c,c1)*c1)/(principle(i) - smallI)
         enddo
         call cross_product3D(dRdt, smallIEvector, omega)
      endif
      !
   end subroutine get_momentofinertia
!-----------------------------------------------------------------------
!+
! LAPACK: DSYEV computes the eigenvalues and, optionally,
!   the left and/or right eigenvectors for SY matrices
! Calls the LAPACK diagonalization subroutine DSYEV
! input:  a(n,n) = real symmetric matrix to be diag
!         n  = size of a
! output: a(n,n) = orthonormal eigenvectors of a
!         v(n) = eigenvalues of a in ascending order
!+
!-----------------------------------------------------------------------
#ifdef LAPACK
   subroutine eigensystem(a,n,v)

      integer, intent(in)    :: n
      real,    intent(inout) :: a(n,n)
      real,    intent(out) :: v(n)

      integer :: lda
      real(kind=8) :: work(3*n-1)
      integer :: lwork
      integer :: info
      integer :: i

      info = 0
      lda = n
      lwork = 3*n-1
      call dsyev('V','U',n,a,lda,v,work,lwork,info)
      if(info < 0) then
         write(*,'(a, i3, a)') "INFO = ", info,&
            " the i-th argument had an illegal value"
      else if(info < 0) then
         write(*,'(a, i3, a)') "INFO = ", info,&
            " the algorithm failed to converge;&
         & i off-diagonal elements of an intermediate tridiagonal&
         & form did not converge to zero."
      endif

   end subroutine eigensystem
#endif
!-----------------------------------------------------------------------
!+
! Calculates the Jacobian
! Source: http://www.fing.edu.uy/if/cursos/fiscomp/extras/numrec/book/f11.pdf
!+
!-----------------------------------------------------------------------
   subroutine jacobi(a,n,np,d,v,nrot)

      integer, intent(in)    :: n,np
      integer, intent(out)   :: nrot
      real,    intent(inout) :: a(np,np)
      real,    intent(out)   :: d(np),v(np,np)
      integer, parameter     :: nmax = 500
!
! Computes all eigenvalues and eigenvectors of a real symmetric matrix, a,
!   which is of size n by n, stored in a physical np by np array.
! On output, elements of a above the diagonal are destroyed.
! d returns the eigenvalues of a in its first n elements.
! v is a matrix with the same logical and physical dimensions as a,
!   whose columns contain, on output, the normalized eigenvectors of a.
! nrot returns the number of Jacobi rotations that were required.
!
      integer :: i,ip,iq,j
      real :: c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)

      do 12, ip=1,n  !Initialize  to  the  identity  matrix.
         do 11, iq=1,n
            v(ip,iq)=0.
11       enddo
         v(ip,ip)=1.
12    enddo
      do 13, ip=1,n
         b(ip)=a(ip,ip)
!Initialize b and d to the diagonal of a.
         d(ip)=b(ip)
         z(ip)=0.
!This  vector  will  accumulate  terms  of  the  form tapq as  in equation  (11.1.14).
13    enddo

      nrot=0
      do 24,i=1,50
         sm=0.
         do 15,ip=1,n-1
!Sum  off-diagonal elements.
            do 14,iq=ip+1,n
               sm=sm+abs(a(ip,iq))
14          enddo
15       enddo
         if(sm==0.)&
            return
!The normal return, which relies on quadratic convergence to machine  underflow.
         if(i < 4) then
            tresh=0.2*sm/n**2
!...on the first  three sweeps.
         else
            tresh=0.
!...thereafter.
         endif
         do 22,ip=1,n-1
            do 21,iq=ip+1,n
               g=100.*abs(a(ip,iq))
!After four sweeps, skip the rotation if the off-diagonal element is small.
               if((i > 4).and.(abs(d(ip))+g==abs(d(ip))).and.(abs(d(iq))+g==abs(d(iq)))) then
                  a(ip,iq)=0.
               elseif (abs(a(ip,iq)) > tresh) then
                  h=d(iq)-d(ip)
                  if(abs(h)+g==abs(h)) then
                     t=a(ip,iq)/h
!t=1/(2(theta))
                  else
                     theta=0.5*h/a(ip,iq)
!Equation  (11.1.10).
                     t=1./(abs(theta)+sqrt(1.+theta**2))
                     if(theta < 0.)t=-t
                  endif
                  c=1./sqrt(1+t**2)
                  s=t*c
                  tau=s/(1.+c)
                  h=t*a(ip,iq)
                  z(ip)=z(ip)-h
                  z(iq)=z(iq)+h
                  d(ip)=d(ip)-h
                  d(iq)=d(iq)+h
                  a(ip,iq)=0.
                  do 16,j=1,ip-1
!Case of rotations 1<=j<p.
                     g=a(j,ip)
                     h=a(j,iq)
                     a(j,ip)=g-s*(h+g*tau)
                     a(j,iq)=h+s*(g-h*tau)
16                enddo
                  do 17,j=ip+1,iq-1
!Case of rotations p<j<q.
                     g=a(ip,j)
                     h=a(j,iq)
                     a(ip,j)=g-s*(h+g*tau)
                     a(j,iq)=h+s*(g-h*tau)
17                enddo
                  do 18,j=iq+1,n
!Case of rotations q<j<=n.
                     g=a(ip,j)
                     h=a(iq,j)
                     a(ip,j)=g-s*(h+g*tau)
                     a(iq,j)=h+s*(g-h*tau)
18                enddo
                  do 19,j=1,n
                     g=v(j,ip)
                     h=v(j,iq)
                     v(j,ip)=g-s*(h+g*tau)
                     v(j,iq)=h+s*(g-h*tau)
19                enddo
                  nrot=nrot+1
               endif
21          enddo
22       enddo
         do 23,ip=1,n
            b(ip)=b(ip)+z(ip)
            d(ip)=b(ip)
!Update d with the  sum of tapq,
            z(ip)=0.
!and  reinitialize z.
23       enddo
24    enddo
      return
   end subroutine jacobi
!-----------------------------------------------------------------------
! Function for finding omega vector from changes of the eigenvector
!-----------------------------------------------------------------------
   function calculate_omega(evector,evector_prev,time,time_prev,omega_prev) result(omega)

      use vectorutils, only: cross_product3D

      real, intent(inout) :: evector(3)
      real, intent(in)    :: evector_prev(3)
      real, intent(in)    :: time
      real, intent(in)    :: time_prev
      real, intent(in)    :: omega_prev(3)

      real                :: omega(3)
      real                :: dRdt(3)

      ! Calculate omega
      dRdt = 0.
      omega = 0.

      if(time > time_prev) then
         ! NB: change sign of evector if need it
         call correct_evector(evector,evector_prev,time,time_prev,omega_prev)
         dRdt = (evector - evector_prev)/(time - time_prev)
         call cross_product3D(evector,dRdt,omega)
         omega = omega/(norm2(evector)**2)
      endif

   end function calculate_omega
!-----------------------------------------------------------------------
   subroutine correct_evector(evector,evector_prev,time,time_prev,omega_prev)

      use vectorutils, only: cross_product3D

      real, intent(inout) :: evector(3)
      real, intent(in)    :: evector_prev(3)
      real, intent(in)    :: time
      real, intent(in)    :: time_prev
      real, intent(in)    :: omega_prev(3)

      real                :: omega(3)
      real                :: dRdt(3)

      ! Calculate omega
      dRdt = 0.
      omega = 0.

      if(time > time_prev) then
         dRdt = (evector - evector_prev)/(time - time_prev)
         call cross_product3D(evector,dRdt,omega)
         ! NB: change sign of evector if need it
         if(omega_prev(3)*omega(3) < 0.0)&
            evector = -1.*evector
      endif

   end subroutine correct_evector
!-----------------------------------------------------------------------
! Correct sign of evectors
!-----------------------------------------------------------------------
   subroutine correct_sign_evector(evector, evector_prev)

      real, intent(inout) :: evector(3)
      real, intent(in)    :: evector_prev(3)

      if (dot_product(evector, evector_prev) < 0.0) then
         evector = -evector
      endif

   end subroutine correct_sign_evector
!-----------------------------------------------------------------------
! Calculate the traceless quadrupole
!-----------------------------------------------------------------------
   subroutine get_internal_quadrupole(xyzh, com, e1, e2, e3, npart, density_cutoff, particlemass, q)

      use part, only: rhoh

      implicit none

      integer, intent(in) :: npart
      real, intent(in)    :: xyzh(:,:)
      real, intent(in)    :: com(3)
      real, intent(in)    :: e1(3), e2(3), e3(3)
      real, intent(in)    :: particlemass
      real, intent(in)    :: density_cutoff
      real, intent(out)   :: q(3,3)

      integer :: i
      real    :: xi(3)
      real    :: coord(3)
      real    :: r2

      q = 0.
      do i = 1, npart
         if(rhoh(xyzh(4,i),particlemass) > density_cutoff) then
            xi = xyzh(1:3,i) - com

            coord(1) = dot_product(xi,e1)
            coord(2) = dot_product(xi,e2)
            coord(3) = dot_product(xi,e3)

            r2 = dot_product(coord,coord)

            q(1,1) = q(1,1) + coord(1)*coord(1) - r2/3.
            q(2,2) = q(2,2) + coord(2)*coord(2) - r2/3.
            q(3,3) = q(3,3) + coord(3)*coord(3) - r2/3.

            q(1,2) = q(1,2) + coord(1)*coord(2)
            q(1,3) = q(1,3) + coord(1)*coord(3)
            q(2,3) = q(2,3) + coord(2)*coord(3)
         endif
      enddo

      q(2,1) = q(1,2)
      q(3,1) = q(1,3)
      q(3,2) = q(2,3)
      q = particlemass*q

      write(*,*) 'Trace = ', q(1,1)+q(2,2)+q(3,3)

   end subroutine get_internal_quadrupole
!-----------------------------------------------------------------------
! Calculate the fifth time derivative of the quadrupole moment
!-----------------------------------------------------------------------
   subroutine dquadrupole5(npart,xyzh,omega,particlemass,d5q)

      integer, intent(in)  :: npart
      real,    intent(in)  :: xyzh(:,:)
      real,    intent(in)  :: omega(3)
      real,    intent(in)  :: particlemass
      real,    intent(out) :: d5q(3,3)

      real                 :: omegasq
      real                 :: omegari
      integer              :: i, ia, ib, ii, ik
      real                 :: coeff(3)

      d5q = 0.d0

      omegasq = dot_product(omega, omega)

!$omp parallel default(none) &
!$omp shared(npart,xyzh,omega,omegasq) &
!$omp private(i,omegari,ii,ia,ib,coeff,ik) &
!$omp reduction(+:d5q)
!$omp do
      do i = 1, npart
         omegari = dot_product(omega, xyzh(1:3,i))

         coeff = 0.d0
         do ii = 1, 3
            do ia = 1, 3
               do ib = 1, 3
                  coeff(ii) = coeff(ii) + levicivita(ii,ia,ib)*omega(ia)*xyzh(ib,i)
               enddo
            enddo
         enddo

         do ii = 1, 3
            do ik = 1, 3
               d5q(ii,ik) = d5q(ii,ik) +&
                  (16.d0*xyzh(ik,i)*omegasq -&
                  15.d0*omega(ik)*omegari)*coeff(ii) +&
                  (16.d0*xyzh(ii,i)*omegasq -&
                  15.d0*omega(ii)*omegari)*coeff(ik)
            enddo
         enddo
      enddo
!$omp enddo
!$omp end parallel

      d5q = d5q*omegasq*particlemass

   end subroutine dquadrupole5
!-----------------------------------------------------------------------
! Returns the Levi-Civita-Symbol (permutation symbol)
!-----------------------------------------------------------------------
   pure function levicivita(i,j,k) result(lc)

      integer, intent(in)           :: i, j, k
      real(8)                       :: lc

      lc = 0.5d0 * (i - j) * (j - k) * (k - i)

   end function levicivita

end module analysis_stripping_tools
