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
   use binary_tools, only: get_momentofinertia, correct_sign_evector, L1_point, &
      get_spin_angular_velocities, get_quadrupole_moments

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
   public :: get_momentofinertia, correct_sign_evector, L1_point, get_spin_angular_velocities, get_quadrupole_moments
   
   ! Interface for HPR with sink support
   interface
      subroutine hpr_check_and_apply_with_sinks(npart,xyzh,vxyzu,massoftype,t,applied)
         integer, intent(in)    :: npart
         real,    intent(inout) :: xyzh(:,:)
         real,    intent(inout) :: vxyzu(:,:)
         real,    intent(in)    :: massoftype(:)
         real,    intent(in)    :: t
         logical, intent(out)   :: applied
      end subroutine hpr_check_and_apply_with_sinks
      subroutine hpr_init_csv()
      end subroutine hpr_init_csv
   end interface
   public :: hpr_check_and_apply_with_sinks, hpr_init_csv
   private

end module halted_pendulum_tools
