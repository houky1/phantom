!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
!--------------------------------------------------------------------------!
module halted_pendulum_relaxation
!
! Hook for halted-pendulum relaxation during the physical evolution.
!
 use io,       only:id,master,iprint
 use mpiutils, only:reduceall_mpi
 use deriv,     only:derivs
 use part,     only:vxyzu,massoftype,igas,isdead_or_accreted,xyzh,fxyzu,fext,divcurlv,divcurlB, &
                   Bevol,dBevol,rad,drad,radprop,dustprop,ddustprop,dustevol,ddustevol,filfac, &
                   dustfrac,eos_vars,pxyzu,dens,metrics,apr_level
 implicit none
 private
 integer, parameter :: nwindow = 7
 integer, parameter :: max_halts = 10
 integer, save :: nhalts = 0
 real,    save :: time_window(nwindow) = 0.
 real,    save :: ekin_window(nwindow) = 0.
 integer, save :: nstored = 0
 logical, save :: first_call = .true.
 public :: apply_halted_pendulum_relax

contains

subroutine apply_halted_pendulum_relax(npart,time,dt,halted)
 integer, intent(inout)  :: npart
 real,    intent(in)  :: time,dt
 logical, intent(out) :: halted
 real :: ekin,t_sample,tmax,aa,bb,cc
 logical :: has_maximum
 integer :: i

 halted = .false.
 if (nhalts >= max_halts) return

 ekin = 0.
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       ekin = ekin + 0.5*massoftype(igas)*sum(vxyzu(1:3,i)**2)
    endif
 enddo
 ekin = reduceall_mpi('+',ekin)
 t_sample = time + dt

 if (first_call .and. id==master) then
    write(iprint,"(a,i0,a)") ' HPR kinetic-energy window: ',nwindow,' samples'
    first_call = .false.
 endif

 if (nstored < nwindow) then
    nstored = nstored + 1
    time_window(nstored) = t_sample
    ekin_window(nstored) = ekin
 else
    time_window(1:nwindow-1) = time_window(2:nwindow)
    ekin_window(1:nwindow-1) = ekin_window(2:nwindow)
    time_window(nwindow) = t_sample
    ekin_window(nwindow) = ekin
 endif

 if (nstored == nwindow) then
    call fit_quadratic(time_window,ekin_window,aa,bb,cc,has_maximum)
    if (has_maximum) then
       tmax = -bb/(2.*aa) + sum(time_window)/real(nwindow)
       if (tmax >= time_window(1) .and. tmax <= time_window(nwindow)) then
          if (id==master) then
             write(iprint,"(a,i3,3(es14.6))") &
             ' HPR halt, number, t, Ekin, tmax = ',nhalts+1,t_sample,ekin,tmax
          endif
          call reset_halted_pendulum_state(npart,time,dt)
          nhalts = nhalts + 1
          nstored = 0
          halted = (nhalts >= max_halts)
       endif
    endif
 endif

end subroutine apply_halted_pendulum_relax

subroutine reset_halted_pendulum_state(npart,time,dt)
 integer, intent(inout) :: npart
 real,    intent(in) :: time,dt
 real :: dtnew

 vxyzu(1:3,1:npart) = 0.
 call derivs(2,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,dBevol, &
             rad,drad,radprop,dustprop,ddustprop,dustevol,ddustevol,filfac,dustfrac, &
             eos_vars,time,dt,dtnew,pxyzu,dens,metrics,apr_level)

 if (id==master) write(iprint,"(a,es14.6)") ' HPR reset: fxyzu recomputed at t = ',time

end subroutine reset_halted_pendulum_state

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

 solution(3) = rhs(3)/matrix(3,3)
 solution(2) = (rhs(2)-matrix(2,3)*solution(3))/matrix(2,2)
 solution(1) = (rhs(1)-matrix(1,2)*solution(2)-matrix(1,3)*solution(3))/matrix(1,1)
 aa = solution(1)
 bb = solution(2)
 cc = solution(3)
 has_maximum = aa < 0.

end subroutine fit_quadratic

end module halted_pendulum_relaxation
