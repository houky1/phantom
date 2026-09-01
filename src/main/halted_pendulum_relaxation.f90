!
! Halted-Pendulum Relaxation for tidally-locked binary stars
!
! Tracks kinetic energy in the instantaneous corotating frame over time,
! fits a polynomial to detect local maximum (pendulum top),
! and at that moment zeroes residual velocities relative to
! the corotating frame, leaving only the orbital component.
!
! The angular velocity omega of the corotating frame is NOT prescribed —
! it is estimated at each step from the rotation of the vector connecting
! the centres of mass of the two stars (tracked over previous steps).
!
! References: Kaltenborn et al. (2023), ApJ 952, DOI 10.3847/1538-4357/acd75a
!
module halted_pendulum_relaxation


   implicit none
   public :: hpr_check_and_apply
   public :: hpr_init


   ! ring buffer of (time, ekin_corot) pairs for polynomial fit
   integer, parameter :: hpr_maxbuf = 20
   real,    private :: hpr_tbuf(hpr_maxbuf) = 0.
   real,    private :: hpr_ebuf(hpr_maxbuf) = 0.
   integer, private :: hpr_nbuf             = 0
   integer, private :: hpr_napplied         = 0

   ! ring buffer for omega estimation from separation vector rotation
   integer, parameter :: hpr_omegabuf = 10
   real,    private :: hpr_sep_x(hpr_omegabuf) = 0.  ! x-component of separation vector
   real,    private :: hpr_sep_y(hpr_omegabuf) = 0.  ! y-component of separation vector
   real,    private :: hpr_sep_t(hpr_omegabuf) = 0.  ! time of each measurement
   integer, private :: hpr_nsep              = 0
   real,    private :: hpr_omega_current     = 0.     ! last estimated omega
   real,    private :: hpr_omega_previous    = 0.     ! previous omega for checking stability

   real, save :: evector_old(3) = (/1.,0.,0./)
   logical, save :: hpr_first_call = .true.
   logical, save :: hpr_initialized = .false.
   real,    save, private :: hpr_time_last_halt = 0.
   logical, save, private :: hpr_finished_latched = .false.
   private

contains

!
! Initializes the internal state of HPR before starting work.
!
   subroutine hpr_init()
      use io,      only:id,master,iprint
      use options, only:hpr_nfit

      if (hpr_initialized) return

      if (hpr_nfit < 3) then
         if (id==master) write(iprint,"(a,i0,a)") &
            ' HPR WARNING: hpr_nfit = ',hpr_nfit,' is too small (<3), using 3'
         hpr_nfit = 3
      else if (hpr_nfit > hpr_maxbuf) then
         if (id==master) write(iprint,"(a,i0,a,i0)") &
            ' HPR WARNING: hpr_nfit = ',hpr_nfit,' exceeds buffer capacity, clamped to ',hpr_maxbuf
         hpr_nfit = hpr_maxbuf
      endif

      hpr_tbuf  = 0.
      hpr_ebuf  = 0.
      hpr_nbuf  = 0
      hpr_sep_x = 0.
      hpr_sep_y = 0.
      hpr_sep_t = 0.
      hpr_nsep = 0
      hpr_omega_current = 0.
      hpr_napplied       = 0
      evector_old        = (/1.,0.,0./)
      hpr_first_call     = .true.
      hpr_time_last_halt = 0.
      hpr_finished_latched = .false.
      hpr_initialized = .true.

   end subroutine hpr_init

   subroutine hpr_check_and_apply(npart,xyzh,vxyzu,massoftype,t,applied,finished)
      use io,                        only:id,master,iprint
      use options,                   only:use_hpr,hpr_nfit,hpr_ekin_tol,hpr_settle_time
      use part,                      only:igas
      use centreofmass,              only:get_centreofmass
      use halted_pendulum_tools,     only:get_momentofinertia,correct_sign_evector,L1_point

      integer, intent(in)    :: npart
      real,    intent(inout) :: xyzh(:,:)
      real,    intent(inout) :: vxyzu(:,:)
      real,    intent(in)    :: massoftype(:)
      real,    intent(in)    :: t
      logical, intent(out)   :: applied
      logical, intent(out)   :: finished

      real :: com(3),vcom(3)
      real :: inertia(3,3),principle(3),evectors(3,3),rmax
      integer :: npartused,smallIIndex
      real :: density_cutoff,particlemass
      real :: sep(3),omega_vec(3)
      real :: L1(3),L1_projection
      real :: com1(3),m1,com2(3),m2
      real :: ekin_corot,ekin_total
      real :: aa,bb,cc,tmax
      logical :: has_maximum
      integer :: i0,i1


      applied  = .false.
      finished = .false.
      if (.not.use_hpr) return

      if (hpr_first_call) then
         if (id==master) write(iprint,"(a,i0,a)") &
            ' HPR: monitoring corotating-frame kinetic energy, window = ',hpr_nfit,' samples'
         hpr_time_last_halt = t
         hpr_first_call = .false.
      endif

      particlemass   = massoftype(igas)
      density_cutoff = 0.

      ! system centre of mass -- required INPUT to get_momentofinertia
      ! (and to the L1 search below), not something it computes
      call get_centreofmass(com, vcom, npart, xyzh, vxyzu)

      ! Calculate the tensor inertia and evectors
      call get_momentofinertia(xyzh, vxyzu, com, vcom, npart, density_cutoff, particlemass,&
         npartused, inertia, principle, evectors, rmax)
      smallIIndex = minloc(principle, dim=1)
      !Correct sign of evector
      call correct_sign_evector(evectors(:, smallIIndex), evector_old)
      evector_old = evectors(:, smallIIndex)

      ! Exact L1 point from the full SPH potential.
      call L1_point(2, xyzh, particlemass, npart, evector_old, com, hpr_omega_current, rmax, L1_projection, L1)

      ! Split into the two stars using the L1 point just found.
      call split_by_axis(npart, xyzh, massoftype, evector_old, L1_projection, com1, m1, com2, m2)
      sep = com1 - com2

      call update_omega_estimate(sep(1),sep(2),t)
      omega_vec = (/0.,0.,hpr_omega_current/)

      call get_kinetic_energies(npart,xyzh,vxyzu,massoftype,omega_vec,com,ekin_corot,ekin_total)

      call push_buffer(t,ekin_corot)


      if (hpr_nbuf >= hpr_nfit) then
         i0 = hpr_nbuf - hpr_nfit + 1
         i1 = hpr_nbuf
         call fit_quadratic(hpr_tbuf(i0:i1),hpr_ebuf(i0:i1),aa,bb,cc,has_maximum)
         if (has_maximum) then
            tmax = -bb/(2.*aa) + sum(hpr_tbuf(i0:i1))/real(hpr_nfit)
            if (tmax >= hpr_tbuf(i0) .and. tmax <= hpr_tbuf(i1)) then
               if (ekin_corot > hpr_ekin_tol*ekin_total) then
                  call zero_residual_velocity(npart,xyzh,vxyzu,omega_vec,com)
                  hpr_nbuf     = 0
                  hpr_napplied = hpr_napplied + 1
                  applied      = .true.
                  if (id==master) then
                     write(iprint,"(a,i4,a,es14.6,a,es14.6,a,es14.6,a,es14.6)") &
                        ' HPR halt #',hpr_napplied,'  t = ',t,'  Ekin_corot = ',ekin_corot,&
                        '  tmax = ',tmax,'  omega = ',hpr_omega_current
                  endif
               endif
            endif
         endif
      endif

      ! Relaxation is considered complete once no halt has been needed
      ! for a sustained period (hpr_settle_time)
      if (applied) then
         hpr_time_last_halt = t
      else if (.not.hpr_finished_latched .and. hpr_settle_time > 0. .and. &
         t - hpr_time_last_halt > hpr_settle_time) then
         finished = .true.
         hpr_finished_latched = .true.
         use_hpr = .false.   ! stop paying for L1_point etc. from the next call on
         if (id==master) then
            ! Recalculate Ekin_corot for the final output
            call get_kinetic_energies(npart,xyzh,vxyzu,massoftype,omega_vec,com,ekin_corot,ekin_total)
            write(iprint,"(a,i0,a)") &
               ' HPR: relaxation complete after ',hpr_napplied,' halt(s) -- no halt needed for hpr_settle_time'
            write(iprint,"(a,2(1x,es14.6),a,es14.6)") &
               ' HPR: final separation |a|, Omega = ',norm2(sep),hpr_omega_current,&
               '  Ekin_corot = ',ekin_corot
         endif
      endif

   end subroutine hpr_check_and_apply

!----------------------------------------------------------------
!+
!  estimate omega (z-component) from the rotation of the separation
!  vector (sepx,sepy) over the last hpr_omegabuf samples, as the mean
!  of consecutive finite-difference angle estimates. Angle differences
!  are wrapped into (-pi,pi] to avoid a spurious +-2*pi spike when the
!  orbital phase crosses the atan2 branch cut.
!+
!----------------------------------------------------------------
   subroutine update_omega_estimate(sepx,sepy,t)
      real, intent(in) :: sepx,sepy,t
      real :: phi_new,phi_old,dphi,dt_k,omega_sum,omega_raw
      real :: omegas(hpr_omegabuf)  ! array for storing individual omega estimates
      integer :: k, nvalid
      real, parameter :: max_omega_change_factor = 2.0  ! maximum allowed change factor
      real, parameter :: min_valid_dt = 1e-10            ! minimum valid time difference

      if (hpr_nsep < hpr_omegabuf) then
         hpr_nsep = hpr_nsep + 1
         hpr_sep_x(hpr_nsep) = sepx
         hpr_sep_y(hpr_nsep) = sepy
         hpr_sep_t(hpr_nsep) = t
      else
         hpr_sep_x(1:hpr_omegabuf-1) = hpr_sep_x(2:hpr_omegabuf)
         hpr_sep_y(1:hpr_omegabuf-1) = hpr_sep_y(2:hpr_omegabuf)
         hpr_sep_t(1:hpr_omegabuf-1) = hpr_sep_t(2:hpr_omegabuf)
         hpr_sep_x(hpr_omegabuf) = sepx
         hpr_sep_y(hpr_omegabuf) = sepy
         hpr_sep_t(hpr_omegabuf) = t
      endif

      if (hpr_nsep < 2) return   ! not enough points yet; keep previous estimate

      ! Store all individual omega estimates for filtering
      omegas = 0.
      nvalid = 0
      omega_sum = 0.

      do k=2,hpr_nsep
         phi_new = atan2(hpr_sep_y(k),  hpr_sep_x(k))
         phi_old = atan2(hpr_sep_y(k-1),hpr_sep_x(k-1))
         dphi    = wrap_angle(phi_new - phi_old)
         dt_k    = hpr_sep_t(k) - hpr_sep_t(k-1)
         if (dt_k > min_valid_dt) then
            omega_raw = dphi/dt_k
            ! Filter out extreme values (outliers)
            if (abs(omega_raw) < 1.0) then  ! reasonable physical limit
               nvalid = nvalid + 1
               omegas(nvalid) = omega_raw
               omega_sum = omega_sum + omega_raw
            endif
         endif
      enddo

      if (nvalid > 0) then
         omega_raw = omega_sum/real(nvalid)

         ! Smooth the omega estimate using exponential moving average
         if (hpr_omega_previous == 0.) then
            hpr_omega_current = omega_raw
         else
            ! Apply smoothing: 0.7 weight to previous, 0.3 to new estimate
            hpr_omega_current = 0.7*hpr_omega_previous + 0.3*omega_raw

            ! Limit the rate of change
            if (abs(hpr_omega_current - hpr_omega_previous) > &
               max_omega_change_factor * abs(hpr_omega_previous)) then
               ! If change is too large, use weighted average closer to previous
               hpr_omega_current = 0.9*hpr_omega_previous + 0.1*omega_raw
            endif
         endif

         hpr_omega_previous = hpr_omega_current
      endif

   end subroutine update_omega_estimate

!----------------------------------------------------------------
!+
!  wrap an angle difference into (-pi,pi], to avoid a spurious
!  +-2*pi jump in the omega estimate when the orbital phase
!  crosses the atan2 branch cut
!+
!----------------------------------------------------------------
   real function wrap_angle(dphi) result(w)
      real, intent(in) :: dphi
      real, parameter :: pi = 4.*atan(1.)
      real, parameter :: two_pi = 2.*pi

      w = dphi
      ! Use modulo for better numerical stability
      ! Bring angle to [0, 2*pi) first
      w = w - two_pi * floor((w + pi) / two_pi)
      ! Convert to (-pi, pi]
      if (w > pi) w = w - two_pi

   end function wrap_angle

!----------------------------------------------------------------
!+
!  kinetic energy of residual motion (relative to the corotating
!  frame defined by omega_vec) and total (lab-frame) kinetic energy,
!  summed over all live particles, reduced across MPI ranks.
!+
!----------------------------------------------------------------
   subroutine get_kinetic_energies(npart,xyzh,vxyzu,massoftype,omega_vec,com_in,ekin_corot,ekin_total)
      use part,        only:iamtype,iphase,isdead_or_accreted
      use mpiutils,    only:reduceall_mpi
      use vectorutils, only:cross_product3D
      integer, intent(in)  :: npart
      real,    intent(in)  :: xyzh(:,:),vxyzu(:,:),massoftype(:)
      real,    intent(in)  :: omega_vec(3)
      real,    intent(in)  :: com_in(3)
      real,    intent(out) :: ekin_corot,ekin_total
      integer :: i
      real :: mi,r(3),v_lab(3),v_orb(3),v_res(3)

      ekin_corot = 0.
      ekin_total = 0.
      do i=1,npart
         if (isdead_or_accreted(xyzh(4,i))) cycle
         mi    = massoftype(iamtype(iphase(i)))
         r     = xyzh(1:3,i) - com_in
         v_lab = vxyzu(1:3,i)
         call cross_product3D(omega_vec,r,v_orb)
         v_res = v_lab - v_orb
         ekin_corot = ekin_corot + 0.5*mi*sum(v_res**2)
         ekin_total = ekin_total + 0.5*mi*sum(v_lab**2)
      enddo
      ekin_corot = reduceall_mpi('+',ekin_corot)
      ekin_total = reduceall_mpi('+',ekin_total)

   end subroutine get_kinetic_energies

!----------------------------------------------------------------
!+
!  push a new (t,ekin_corot) sample into the ring buffer
!+
!----------------------------------------------------------------
   subroutine push_buffer(t,ekin)
      real, intent(in) :: t,ekin

      if (hpr_nbuf < hpr_maxbuf) then
         hpr_nbuf = hpr_nbuf + 1
         hpr_tbuf(hpr_nbuf) = t
         hpr_ebuf(hpr_nbuf) = ekin
      else
         hpr_tbuf(1:hpr_maxbuf-1) = hpr_tbuf(2:hpr_maxbuf)
         hpr_ebuf(1:hpr_maxbuf-1) = hpr_ebuf(2:hpr_maxbuf)
         hpr_tbuf(hpr_maxbuf) = t
         hpr_ebuf(hpr_maxbuf) = ekin
      endif

   end subroutine push_buffer

!----------------------------------------------------------------
!+
!  zero the residual velocity of every particle, keeping only the
!  orbital component v_orb = Omega x r. fxyzu is NOT recomputed here
!  -- see the note at the call site in hpr_check_and_apply.
!+
!----------------------------------------------------------------
   subroutine zero_residual_velocity(npart,xyzh,vxyzu,omega_vec,com_in)
      use part,        only:isdead_or_accreted
      use vectorutils, only:cross_product3D
      integer, intent(in)    :: npart
      real,    intent(in)    :: xyzh(:,:)
      real,    intent(inout) :: vxyzu(:,:)
      real,    intent(in)    :: omega_vec(3)
      real,    intent(in)    :: com_in(3)
      integer :: i

      do i=1,npart
         if (isdead_or_accreted(xyzh(4,i))) cycle
         call cross_product3D(omega_vec,xyzh(1:3,i)-com_in,vxyzu(1:3,i))
      enddo

   end subroutine zero_residual_velocity

!----------------------------------------------------------------
!+
!  split particles into two groups by which side of `threshold` their
!  projection onto `axis` falls, accumulating each group's mass and
!  centre of mass directly (no per-particle data is copied). MPI-safe,
!  skips dead/accreted particles.
!+
!----------------------------------------------------------------
   subroutine split_by_axis(npart,xyzh,massoftype,axis,threshold,com1,m1,com2,m2)
      use part,     only:iamtype,iphase,isdead_or_accreted
      use mpiutils, only:reduceall_mpi
      integer, intent(in)  :: npart
      real,    intent(in)  :: xyzh(:,:),massoftype(:)
      real,    intent(in)  :: axis(3),threshold
      real,    intent(out) :: com1(3),m1,com2(3),m2
      integer :: i,k
      real :: mi

      com1 = 0.; m1 = 0.
      com2 = 0.; m2 = 0.
      do i=1,npart
         if (isdead_or_accreted(xyzh(4,i))) cycle
         mi = massoftype(iamtype(iphase(i)))
         if (dot_product(xyzh(1:3,i),axis) >= threshold) then
            com1 = com1 + mi*xyzh(1:3,i)
            m1   = m1 + mi
         else
            com2 = com2 + mi*xyzh(1:3,i)
            m2   = m2 + mi
         endif
      enddo
      ! reduceall_mpi is only used here as a scalar reduction elsewhere
      ! in this file (see get_kinetic_energies) -- reduce component-wise
      ! rather than assume the interface is overloaded for rank-1 args
      do k=1,3
         com1(k) = reduceall_mpi('+',com1(k))
         com2(k) = reduceall_mpi('+',com2(k))
      enddo
      m1 = reduceall_mpi('+',m1)
      m2 = reduceall_mpi('+',m2)
      if (m1 > 0.) com1 = com1/m1
      if (m2 > 0.) com2 = com2/m2

   end subroutine split_by_axis

!----------------------------------------------------------------
!+
!  least-squares quadratic fit y = aa*dx**2 + bb*dx + cc, where
!  dx = x - mean(x), solved via Gaussian elimination on the 3x3
!  normal-equations matrix. has_maximum is true iff aa < 0.
!+
!----------------------------------------------------------------
   subroutine fit_quadratic(x,y,aa,bb,cc,has_maximum)
      real,    intent(in)  :: x(:),y(:)
      real,    intent(out) :: aa,bb,cc
      logical, intent(out) :: has_maximum
      real :: xm,dx(size(x)),matrix(3,3),rhs(3),solution(3),factor
      integer :: i,j,k

      xm = sum(x)/real(size(x))
      dx = x - xm
      matrix = 0.
      rhs = 0.
      do i=1,size(x)
         matrix(1,1) = matrix(1,1) + dx(i)**4
         matrix(1,2) = matrix(1,2) + dx(i)**3
         matrix(1,3) = matrix(1,3) + dx(i)**2
         matrix(2,1) = matrix(1,2)
         matrix(2,2) = matrix(2,2) + dx(i)**2
         matrix(2,3) = matrix(2,3) + dx(i)
         matrix(3,1) = matrix(1,3)
         matrix(3,2) = matrix(2,3)
         matrix(3,3) = matrix(3,3) + 1.
         rhs(1) = rhs(1) + y(i)*dx(i)**2
         rhs(2) = rhs(2) + y(i)*dx(i)
         rhs(3) = rhs(3) + y(i)
      enddo

      do k=1,2
         if (abs(matrix(k,k)) <= tiny(matrix(k,k))) then
            has_maximum = .false.
            aa = 0.; bb = 0.; cc = 0.
            return
         endif
         do i=k+1,3
            factor = matrix(i,k)/matrix(k,k)
            do j=k,3
               matrix(i,j) = matrix(i,j) - factor*matrix(k,j)
            enddo
            rhs(i) = rhs(i) - factor*rhs(k)
         enddo
      enddo

      if (abs(matrix(3,3)) <= tiny(matrix(3,3))) then
         has_maximum = .false.
         aa = 0.; bb = 0.; cc = 0.
         return
      endif

      solution(3) = rhs(3)/matrix(3,3)
      solution(2) = (rhs(2)-matrix(2,3)*solution(3))/matrix(2,2)
      solution(1) = (rhs(1)-matrix(1,2)*solution(2)-matrix(1,3)*solution(3))/matrix(1,1)
      aa = solution(1)
      bb = solution(2)
      cc = solution(3)
      has_maximum = aa < 0.

   end subroutine fit_quadratic

end module halted_pendulum_relaxation
