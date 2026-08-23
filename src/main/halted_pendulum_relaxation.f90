module halted_pendulum_relaxation
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
 integer, parameter :: hpr_omegabuf = 5
 real,    private :: hpr_sep_x(hpr_omegabuf) = 0.  ! x-component of separation vector
 real,    private :: hpr_sep_y(hpr_omegabuf) = 0.  ! y-component of separation vector
 real,    private :: hpr_sep_t(hpr_omegabuf) = 0.  ! time of each measurement
 integer, private :: hpr_nsep              = 0
 real,    private :: hpr_omega_current     = 0.     ! last estimated omega

 real, save :: evector_old(3) = (/1.,0.,0./)
 logical, save :: hpr_first_call = .true.
 logical, save :: hpr_initialized = .false.
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
   hpr_initialized = .true.

end subroutine hpr_init

subroutine hpr_check_and_apply(npart,xyzh,vxyzu,massoftype,t,applied)
 use io,                        only:id,master,iprint
 use options,                   only:use_hpr,hpr_nfit,hpr_ekin_tol
 use part,                      only:igas
 use centreofmass,               only:get_centreofmass
 use analysis_stripping_tools,  only:get_momentofinertia,correct_sign_evector
 integer, intent(in)    :: npart
 real,    intent(inout) :: xyzh(:,:)
 real,    intent(inout) :: vxyzu(:,:)
 real,    intent(in)    :: massoftype(:)
 real,    intent(in)    :: t
 logical, intent(out)   :: applied
 
 real :: com(3),vcom(3)
 real :: inertia(3,3),principle(3),evectors(3,3),rmax
 integer :: npartused,smallIIndex
 real :: density_cutoff,particlemass
 real :: com1(3),m1,com2(3),m2,sep(3),omega_vec(3)
 real :: threshold,L1_pos(3),L1_proj
 real :: ekin_corot,ekin_total
 real :: aa,bb,cc,tmax
 logical :: has_maximum
 integer :: i0,i1
 
 
 applied = .false.
 if (.not.use_hpr) return
 
 if (hpr_first_call .and. id==master) then
    write(iprint,"(a,i0,a)") ' HPR: monitoring corotating-frame kinetic energy, window = ',hpr_nfit,' samples'
    hpr_first_call = .false.
 endif
 
 particlemass   = massoftype(igas)
 density_cutoff = 0.
 
 call get_centreofmass(com,vcom,npart,xyzh,vxyzu)
 call get_momentofinertia(xyzh,vxyzu,com,vcom,npart,density_cutoff,particlemass, &
                           npartused,inertia,principle,evectors,rmax)
 smallIIndex = minloc(principle,dim=1)
 call correct_sign_evector(evectors(:,smallIIndex),evector_old)
 evector_old = evectors(:,smallIIndex)
 
 ! It primarily separates two stars
 threshold = dot_product(com,evector_old)
 call split_by_axis(npart,xyzh,massoftype,evector_old,threshold,com1,m1,com2,m2)
 
 if (m1 <= 0. .or. m2 <= 0.) then
    if (id==master) write(iprint,"(a)") ' HPR WARNING: one side of the split is empty -- skipping this step'
    return
 endif
 
 call estimate_L1(m1,com1,m2,com2,evector_old,L1_pos,L1_proj)
 call split_by_axis(npart,xyzh,massoftype,evector_old,L1_proj,com1,m1,com2,m2)
 
 if (m1 <= 0. .or. m2 <= 0.) then
    if (id==master) write(iprint,"(a)") ' HPR WARNING: one side of the refined split is empty -- skipping this step'
    return
 endif
 
 sep = com2 - com1
 
 call update_omega_estimate(sep(1),sep(2),t)
 omega_vec = (/0.,0.,hpr_omega_current/)

 call get_kinetic_energies(npart,xyzh,vxyzu,massoftype,omega_vec,ekin_corot,ekin_total)

 call push_buffer(t,ekin_corot)
 
 
 if (hpr_nbuf >= hpr_nfit) then
    i0 = hpr_nbuf - hpr_nfit + 1
    i1 = hpr_nbuf
    call fit_quadratic(hpr_tbuf(i0:i1),hpr_ebuf(i0:i1),aa,bb,cc,has_maximum)
    if (has_maximum) then
       tmax = -bb/(2.*aa) + sum(hpr_tbuf(i0:i1))/real(hpr_nfit)
       if (tmax >= hpr_tbuf(i0) .and. tmax <= hpr_tbuf(i1)) then
          if (ekin_corot > hpr_ekin_tol*ekin_total) then
             call zero_residual_velocity(npart,xyzh,vxyzu,omega_vec)
             hpr_nbuf     = 0
             hpr_napplied = hpr_napplied + 1
             applied      = .true.
             if (id==master) then
                write(iprint,"(a,i4,4(1x,es14.6))") &
                   ' HPR halt #',hpr_napplied,t,ekin_corot,tmax,hpr_omega_current
             endif
          endif
       endif
    endif
 endif
 
end subroutine hpr_check_and_apply
 
!----------------------------------------------------------------
!+
!  split particles into two groups by the sign of
!  dot_product(position,axis) - threshold, accumulating each group's
!  mass-weighted centre of mass directly (no array duplication).
!+
!----------------------------------------------------------------
subroutine split_by_axis(npart,xyzh,massoftype,axis,threshold,xcm1,m1,xcm2,m2)
 use part,     only:iamtype,iphase,isdead_or_accreted
 use mpiutils, only:reduceall_mpi
 integer, intent(in)  :: npart
 real,    intent(in)  :: xyzh(:,:)
 real,    intent(in)  :: massoftype(:)
 real,    intent(in)  :: axis(3),threshold
 real,    intent(out) :: xcm1(3),m1,xcm2(3),m2
 integer :: i
 real :: mi
 
 xcm1 = 0.; m1 = 0.
 xcm2 = 0.; m2 = 0.
 do i=1,npart
    if (isdead_or_accreted(xyzh(4,i))) cycle
    mi = massoftype(iamtype(iphase(i)))
    if (dot_product(xyzh(1:3,i),axis) >= threshold) then
       xcm1 = xcm1 + mi*xyzh(1:3,i)
       m1   = m1 + mi
    else
       xcm2 = xcm2 + mi*xyzh(1:3,i)
       m2   = m2 + mi
    endif
 enddo
 xcm1(1) = reduceall_mpi('+',xcm1(1)); xcm1(2) = reduceall_mpi('+',xcm1(2)); xcm1(3) = reduceall_mpi('+',xcm1(3))
 xcm2(1) = reduceall_mpi('+',xcm2(1)); xcm2(2) = reduceall_mpi('+',xcm2(2)); xcm2(3) = reduceall_mpi('+',xcm2(3))
 m1 = reduceall_mpi('+',m1)
 m2 = reduceall_mpi('+',m2)
 if (m1 > 0.) xcm1 = xcm1/m1
 if (m2 > 0.) xcm2 = xcm2/m2
 
end subroutine split_by_axis
 
!----------------------------------------------------------------
!+
!  closed-form L1 estimate for two POINT masses m1 (at xcm1) and m2
!  (at xcm2) on a circular orbit, solved along the axis joining them.
!  Standard two-body Roche force-balance equation in the corotating
!  frame, non-dimensionalised with separation a=1, G(m1+m2)=1, q=m2/m1:
!
!     g(f) = 1/f**2 - q/(1-f)**2 - (1+q)*f + q = 0
!
!  where f in (0,1) is the fractional distance of L1 from m1 toward m2.
!  Solved by Newton's method from f0=0.5 (exact for q=1 by symmetry).
!  This is a point-mass approximation: it does not account for the
!  extended density profile of either star, unlike the SPH-potential
!  L1_point() in analysis_stripping.f90 (see module header for why
!  that routine is not used here).
!+
!----------------------------------------------------------------
subroutine estimate_L1(m1,xcm1,m2,xcm2,axis,L1_pos,L1_proj)
 real, intent(in)  :: m1,xcm1(3),m2,xcm2(3),axis(3)
 real, intent(out) :: L1_pos(3),L1_proj
 real :: q,f,g,dgdf,a_sep
 integer :: iter
 integer, parameter :: maxiter = 30
 real,    parameter :: ftol = 1.e-10
 
 a_sep = sqrt(sum((xcm2-xcm1)**2))
 if (a_sep <= 0. .or. m1 <= 0. .or. m2 <= 0.) then
    L1_pos  = 0.5*(xcm1+xcm2)
    L1_proj = dot_product(L1_pos,axis)
    return
 endif
 
 q = m2/m1
 f = 0.5
 do iter=1,maxiter
    g    = 1./f**2 - q/(1.-f)**2 - (1.+q)*f + q
    dgdf = -2./f**3 - 2.*q/(1.-f)**3 - (1.+q)
    if (abs(dgdf) <= tiny(dgdf)) exit
    f = f - g/dgdf
    f = min(max(f,1.e-3),1.-1.e-3)
    if (abs(g) < ftol) exit
 enddo
 
 L1_pos  = xcm1 + f*(xcm2-xcm1)
 L1_proj = dot_product(L1_pos,axis)
 
end subroutine estimate_L1
 
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
 real :: phi_new,phi_old,dphi,dt_k,omega_sum
 integer :: k
 
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
 
 omega_sum = 0.
 do k=2,hpr_nsep
    phi_new = atan2(hpr_sep_y(k),  hpr_sep_x(k))
    phi_old = atan2(hpr_sep_y(k-1),hpr_sep_x(k-1))
    dphi    = wrap_angle(phi_new - phi_old)
    dt_k    = hpr_sep_t(k) - hpr_sep_t(k-1)
    if (dt_k > 0.) omega_sum = omega_sum + dphi/dt_k
 enddo
 hpr_omega_current = omega_sum/real(hpr_nsep-1)
 
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
 
 w = dphi
 do while (w > pi)
    w = w - 2.*pi
 enddo
 do while (w <= -pi)
    w = w + 2.*pi
 enddo
 
end function wrap_angle
 
!----------------------------------------------------------------
!+
!  kinetic energy of residual motion (relative to the corotating
!  frame defined by omega_vec) and total (lab-frame) kinetic energy,
!  summed over all live particles, reduced across MPI ranks.
!+
!----------------------------------------------------------------
subroutine get_kinetic_energies(npart,xyzh,vxyzu,massoftype,omega_vec,ekin_corot,ekin_total)
 use part,        only:iamtype,iphase,isdead_or_accreted
 use mpiutils,    only:reduceall_mpi
 use vectorutils, only:cross_product3D
 integer, intent(in)  :: npart
 real,    intent(in)  :: xyzh(:,:),vxyzu(:,:),massoftype(:)
 real,    intent(in)  :: omega_vec(3)
 real,    intent(out) :: ekin_corot,ekin_total
 integer :: i
 real :: mi,r(3),v_lab(3),v_orb(3),v_res(3)
 
 ekin_corot = 0.
 ekin_total = 0.
 do i=1,npart
    if (isdead_or_accreted(xyzh(4,i))) cycle
    mi    = massoftype(iamtype(iphase(i)))
    r     = xyzh(1:3,i)
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
subroutine zero_residual_velocity(npart,xyzh,vxyzu,omega_vec)
 use part,        only:isdead_or_accreted
 use vectorutils, only:cross_product3D
 integer, intent(in)    :: npart
 real,    intent(in)    :: xyzh(:,:)
 real,    intent(inout) :: vxyzu(:,:)
 real,    intent(in)    :: omega_vec(3)
 integer :: i
 
 do i=1,npart
    if (isdead_or_accreted(xyzh(4,i))) cycle
    call cross_product3D(omega_vec,xyzh(1:3,i),vxyzu(1:3,i))
 enddo
 
end subroutine zero_residual_velocity
 
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