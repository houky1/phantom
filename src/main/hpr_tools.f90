!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!

! #define LAPACK

module halted_pendulum_tools
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
   use binary_tools, only: get_momentofinertia, correct_sign_evector, L1_point

   implicit none

   !
   ! Runtime parameters
   !

   integer,    public :: nstar(2) = 0 ! give default value in case dump header not read

   real,       public :: com(3)
   real,       public :: vcom(3)


   real, save, public :: evector_old(3)
   real, save, public :: time_old(2)
   real, save, public :: omega_old(3)

   !
   ! subroutines (re-exported from binary_tools)
   !
   public :: get_momentofinertia, correct_sign_evector, L1_point
   private

end module halted_pendulum_tools