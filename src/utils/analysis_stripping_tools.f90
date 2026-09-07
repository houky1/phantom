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
! :Dependencies: binary_tools, centreofmass, dump_utils, infile_utils, io, physcon, units
!
   use binary_tools, only: get_momentofinertia, correct_sign_evector

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
   ! subroutines (get_momentofinertia, correct_sign_evector re-exported
   ! from binary_tools -- see there for the implementation)
   !
   public :: get_momentofinertia, calculate_omega, correct_evector, &
      correct_sign_evector, get_internal_quadrupole
   private

contains
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
