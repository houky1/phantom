!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!

! #define LAPACK

module binary_tools
!
! Common tools for binary-system analysis, shared between the
!   halted-pendulum/GW-inspiral runtime module (halted_pendulum_tools,
!   used during the actual calculation) and the stripping analysis
!   module (analysis_stripping / analysis_stripping_tools, used in
!   post-processing). Purely computational: this module holds no
!   state of its own (no com/vcom/evector_old/...) -- each of the two
!   consumers keeps its own copy of that, as before.
!   - Moment of inertia tensor
!   - L1 point (exact, from the full SPH gravitational/Roche potential)
!   - Gravitational/Roche potentials and forces used to locate L1
!
! :References: e.g. Tong (2015) classical dynamics lecture notes
!
! :Owner: Daniel Price
!
! :Dependencies: part, vectorutils
!
   implicit none

   type :: l1_context
      real,    allocatable :: xyzh(:,:)
      real                  :: particlemass
      integer               :: npart
      real                  :: axis(3)
      real                  :: com(3)
      real                  :: omega
   end type l1_context
   !
   ! subroutines
   !
   public :: get_momentofinertia, correct_sign_evector, L1_point
   public :: roche_potential, gravitational_potential
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
      else if(info > 0) then
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
! Exact L1 point from the full SPH gravitational/Roche potential.
! axis is the long axis of the binary (e.g. the eigenvector of the
! smallest principal moment of inertia), com is the system centre of
! mass, and omega is the angular speed of the corotating frame -- all
! supplied explicitly by the caller, not read from module state.
!-----------------------------------------------------------------------
   subroutine L1_point(method, xyzh, particlemass, npart, axis, com_in, omega, rmax, L1_proj, L1)

      integer, intent(in)  :: method ! 1 - Newton, 2 - Golden Section Search
      integer, intent(in)  :: npart
      real,    intent(in)  :: xyzh(4,npart)
      real,    intent(in)  :: particlemass
      real,    intent(in)  :: axis(3)
      real,    intent(in)  :: com_in(3)
      real,    intent(in)  :: omega
      real,    intent(in)  :: rmax    ! search bracket is +-rmax around com_in, along axis
      real,    intent(out) :: L1_proj
      real,    intent(out) :: L1(3)

      integer              :: nIter
      integer, parameter   :: maxIter = 50
      real,    parameter   :: threshold = 1e-7
      type(l1_context) :: ctx

      ! for Newton's method
      real                 :: p, pNew, residual, f

      ! for Golden section search method
      real                 :: p1, p2
      integer              :: nest, extr
      real                 :: eps, err

      L1_proj = 0.
      L1 = 0.

      ! package everything the potential/force evaluation needs, to be
      ! threaded explicitly through func(point, ctx, ...) below
      ctx%xyzh         = xyzh(1:4,1:npart)
      ctx%particlemass = particlemass
      ctx%npart        = npart
      ctx%axis         = axis
      ctx%com          = com_in
      ctx%omega        = omega

      ! Initial values for number of iteration
      nIter = 1

      p = 0.

      if(method == 1) then
         ! Initial values for residual
         residual = 1.e10

         ! Keep search iteration until
         ! (a) residual is bigger then a user-defined threshold value, and
         ! (b) iteration number is less than a user-defined maximum iteration number.

         do while ((residual > threshold) .and. (nIter < maxIter))

            ! Search using conventional Newton's method
            call newton_method(p, ctx, pNew, f, residual, gravitational_force_wrapper)

            ! Save for the next search iteration
            p = pNew

            ! Update iteration number
            nIter = nIter + 1

         end do

         L1 = com_in + p*axis

      else

         call bracket_between_stars(ctx, rmax, p1, p2)
         eps = threshold

         p = golden_section_search_method(p1, p2,&
            eps, err, extr,&
            maxIter, Nest, nIter,&
            ctx, roche_potential_wrapper)

         L1 = com_in + p*axis

      endif

      L1_proj = dot_product(L1, axis)

   end subroutine L1_point
!-----------------------------------------------------------------------
! Coarsely scans the axis on each side of the COM (p>0 and p<0) for the
! point of most negative Roche potential -- an approximate location of
! each star -- and returns the interval strictly between them for the
! golden-section search. Falls back to +-rmax if the two are
! degenerate (e.g. too few particles, or both stars on the same side).
!-----------------------------------------------------------------------
   subroutine bracket_between_stars(ctx, rmax, p1, p2)
      type(l1_context), intent(in)  :: ctx
      real,                  intent(in)  :: rmax
      real,                  intent(out) :: p1, p2

      integer, parameter :: nscan = 40
      integer :: iscan
      real :: pscan, potscan
      real :: pot_min_pos, p_min_pos, pot_min_neg, p_min_neg
      real :: margin

      pot_min_pos = huge(1.); p_min_pos = 0.5*rmax
      pot_min_neg = huge(1.); p_min_neg = -0.5*rmax

      do iscan = 1, nscan/2
         pscan = rmax*real(iscan)/real(nscan/2)
         call roche_potential_wrapper(pscan, ctx, potscan)
         if (potscan < pot_min_pos) then
            pot_min_pos = potscan
            p_min_pos = pscan
         endif
      enddo

      do iscan = 1, nscan/2
         pscan = -rmax*real(iscan)/real(nscan/2)
         call roche_potential_wrapper(pscan, ctx, potscan)
         if (potscan < pot_min_neg) then
            pot_min_neg = potscan
            p_min_neg = pscan
         endif
      enddo

      p1 = min(p_min_pos, p_min_neg)
      p2 = max(p_min_pos, p_min_neg)

      if (p2 - p1 < 1.e-6*rmax) then
         ! degenerate scan (e.g. both stars landed on the same side) --
         ! fall back to the full extent rather than a zero-width bracket
         p1 = -rmax
         p2 =  rmax
      else
         margin = 0.15*(p2 - p1)
         p1 = p1 + margin
         p2 = p2 - margin
      endif

   end subroutine bracket_between_stars
   !-----------------------------------------------------------------
   subroutine gravitational_force(p, xyzh, particlemass, npart, axis, com_in, force, dforce)

      use part, only: isdead_or_accreted

      real,           intent (in)  :: p
      integer,        intent (in)  :: npart
      real,           intent (in)  :: xyzh(4,npart)
      real,           intent (in)  :: particlemass
      real,           intent (in)  :: axis(3)
      real,           intent (in)  :: com_in(3)
      real,           intent (out) :: force
      real,           intent (out) :: dforce

      integer                      :: i
      real(kind=8)                 :: point(3)
      real(kind=8)                 :: dpoint(3)
      real(kind=8)                 :: f(3)
      real(kind=8)                 :: df(6)
      real(kind=8)                 :: dr  ! 1/sqrt(r^2)
      real(kind=8)                 :: dr3 ! 1/sqrt(r^2)^3
      real(kind=8)                 :: dr5 ! 1/sqrt(r^2)^5

      force = 0.
      dforce = 0.

      f = 0.
      df = 0.

      point = com_in + p*axis

      do i = 1, npart

         if (isdead_or_accreted(xyzh(4,i))) cycle

         dpoint = point - xyzh(1:3,i)
         dr = 1./norm2(dpoint)
         dr3 = dr*dr*dr
         dr5 = dr3*dr*dr

         f = f - dpoint*dr3

         ! NB: Check for correctness
         df(1) = df(1) + dr5*(3.*dpoint(1)*dpoint(1) - 1.) ! dfx/dx
         df(2) = df(2) + dr5*(3.*dpoint(1)*dpoint(2))      ! dfx/dy = dfy/dx
         df(3) = df(3) + dr5*(3.*dpoint(1)*dpoint(3))      ! dfx/dz = dfz/dx
         df(4) = df(4) + dr5*(3.*dpoint(2)*dpoint(2) - 1.) ! dfy/dy
         df(5) = df(5) + dr5*(3.*dpoint(2)*dpoint(3))      ! dfy/dz = dfz/dy
         df(6) = df(6) + dr5*(3.*dpoint(3)*dpoint(3) - 1.) ! dfz/dz

      enddo

      force = norm2(f)*particlemass
      dforce = norm2(df)*particlemass

   end subroutine gravitational_force
!-----------------------------------------------------------------------
   subroutine gravitational_force_wrapper(p, ctx, force, dforce)

      real,                  intent (in)  :: p
      type(l1_context),  intent (in)  :: ctx
      real,                  intent (out) :: force
      real,                  intent (out) :: dforce

      call gravitational_force(p, ctx%xyzh, ctx%particlemass, ctx%npart, ctx%axis, ctx%com, force, dforce)

   end subroutine gravitational_force_wrapper
!-----------------------------------------------------------------------
   subroutine roche_potential(p, xyzh, particlemass, npart, axis, com_in, omega, potential)

      real,    intent (in)  :: p
      integer, intent (in)  :: npart
      real,    intent (in)  :: xyzh(4,npart)
      real,    intent (in)  :: particlemass
      real,    intent (in)  :: axis(3)
      real,    intent (in)  :: com_in(3)
      real,    intent (in)  :: omega
      real,    intent (out) :: potential

      potential = 0.

      call gravitational_potential(p, xyzh, particlemass, npart, axis, com_in, potential)

      potential = potential - 0.5*(omega*omega)*(p*p)

   end subroutine roche_potential
!-----------------------------------------------------------------------
   subroutine gravitational_potential(p, xyzh, particlemass, npart, axis, com_in, potential)

      use part, only: isdead_or_accreted

      real,           intent (in)  :: p
      integer,        intent (in)  :: npart
      real,           intent (in)  :: xyzh(4,npart)
      real,           intent (in)  :: particlemass
      real,           intent (in)  :: axis(3)
      real,           intent (in)  :: com_in(3)
      real,           intent (out) :: potential

      integer                      :: i
      real(kind=8)                 :: point(3)
      real(kind=8)                 :: dpoint(3)
      real(kind=8)                 :: dr ! 1/sqrt(r^2)

      potential = 0.

      point = com_in + p*axis

!$omp parallel default(none) &
!$omp shared(npart,xyzh,point) &
!$omp private(i,dpoint,dr) &
!$omp reduction(-:potential)
!$omp do
      do i = 1, npart
         if (isdead_or_accreted(xyzh(4,i))) cycle
         dpoint = point - xyzh(1:3,i)
         dr = 1./norm2(dpoint)
         potential = potential - dr
      enddo
!$omp enddo
!$omp end parallel

      potential = potential*particlemass

   end subroutine gravitational_potential
!---------------------------------------------------------------------------------
   subroutine roche_potential_wrapper(p, ctx, potential)

      real,                  intent (in)  :: p
      type(l1_context),  intent (in)  :: ctx
      real,                  intent (out) :: potential

      call roche_potential(p, ctx%xyzh, ctx%particlemass, ctx%npart, ctx%axis, ctx%com, ctx%omega, potential)

   end subroutine roche_potential_wrapper
! Finding the zero of gravitational force by Newton method
! Due to the noise in the function this method is unsuccessful
!-----------------------------------------------------------------------
   subroutine newton_method(p, ctx, pNew, fNew, residual, func)

      real,                  intent(in)  :: p
      type(l1_context),  intent(in)  :: ctx
      real,                  intent(out) :: fNew, residual

      real                 :: pNew, f, df

      interface
         subroutine func(point, ctx, force, dforce)
            import :: l1_context
            real,                  intent (in)  :: point
            type(l1_context),  intent (in)  :: ctx
            real,                  intent (out) :: force
            real,                  intent (out) :: dforce
         end subroutine func
      end interface

      ! compute function value evaluated at x
      call func(p, ctx, f, df)

      ! numerical second derivative
      ! write(*,*) p, f, df
      ! pNew = p + 1.e-4
      ! call func(pNew, fNew, df)
      ! df = (fNew - f)/(pNew - p)
      ! write(*,*) pNew, f, df
      ! stop

      ! Exit if f' is near or become zero
      if(abs(df) < 1.e-12) then
         print *, '[Error: newton_method] Function derivative becomes very close to zero or zero.'
         print *, 'f=',f, 'df/dp =',df
         print *, 'Aborting now in order to avoid division by zero.'
         stop
      end if

      ! Algorithm
      pNew = p - f/df
      fNew = f

      ! Search fails if a newly updated value x is out of the search domain
      ! if((pNew < pBeg) .or. (pNew > pEnd)) then
      !   print *, '[Error: newton_method] pNew',pNew, 'is out of domain.'
      !   print *, 'Failed in search. Aborting now.'
      !   stop
      ! end if

      ! Calculate a new residual
      residual = abs(pNew - p)

   end subroutine newton_method
!-----------------------------------------------------------------------
! The golden-section search is a technique for finding an extremum
!   (minimum or maximum) of a function inside a specified interval.
! The implementation is based on the more robust approach described in
!   V.G. Karmanov Mathematical programming, Moscow: FML, 2008, pp. 134-142.
!-----------------------------------------------------------------------
   real function golden_section_search_method(a, b,&
      eps, err, extr, maxIter, Nest, iter, ctx, func)

      real,    intent(in)  :: a, b        ! left and right boundaries
      ! of the extremum search interval
      real,    intent(in)  :: eps        ! specified accuracy

      real,    intent(out) :: err        ! achieved accuracy
      integer, intent(out) :: extr       ! extremum type: 1 - minimum; -1 - maximum

      integer, intent(in)  :: maxIter
      integer, intent(out) :: iter, Nest ! number of iterations actually made
      ! and their lower bound
      type(l1_context), intent(in) :: ctx

      real                 :: q = 0.5d0*(sqrt(5.0d0)-1.0d0),&
         alpha, fy0, fz0,&
         a0, a1, b0, b1, y0,&
         y1, z0, z1, d0, d1, d2, d3, d10
      integer              :: nfail

      interface
         subroutine func(point, ctx, potential)
            import :: l1_context
            real,                  intent (in)  :: point
            type(l1_context),  intent (in)  :: ctx
            real,                  intent (out) :: potential
         end subroutine func
      end interface

      golden_section_search_method = 0.

      iter = 1
      err = 1.0d0
      nfail = 0
      Nest = int(log(eps/(b-a))/log(q))
      alpha = 0.8d0

      fy0 = 0.; fz0 = 0.; a0 = 0.; a1 = 0.; b0 = 0.; b1 = 0.
      y0 = 0.; y1 = 0.; z0 = 0.; z1 = 0.; d0 = 0.
      d1 = 0.; d2 = 0.; d3 = 0.; d10 = 0.

      call func(a, ctx, fy0)
      call func(a+0.5*(b-a), ctx, fz0)
      if(fy0 > fz0) then
         extr = 1         ! looking for a minimum
      else
         extr = -1        ! looking for a maximum
      end if
! step 1
      a0 = a
      b0 = b
! step 2
      do
         d10 = d1
         d0 = b0-a0
         d1 = q*d0
         d2 = d0-d1
! checking if precision is achieved (algorithm loops)
         if(abs(d1-d10) <= epsilon(1.0d0) .and. iter > maxIter) then
            nfail = nfail+1
! exit after two consecutive non-decreasing precision
            if(nfail >= 2) then
               err = d1
               golden_section_search_method = 0.5d0*(a1+b1)
               return
            end if
         else
            nfail = 0
         end if

         y0 = a0+d2
         z0 = b0-d2

         call func(y0, ctx, fy0)
         call func(z0, ctx, fz0)
! step 3
         do
            iter = iter+1
            d3 = d1-d2

            if(extr*fy0 <= extr*fz0) then
               a1 = a0
               b1 = z0
               z1 = y0
               y1 = a1+d3
               fz0 = fy0
               call func(y1, ctx, fy0)

               if(y1 >= z1) then
                  a0 = a1
                  b0 = b1
                  exit   ! to step 2
               end if
            else
               a1 = y0
               b1 = b0
               y1 = z0
               z1 = b1-d3
               fy0 = fz0
               call func(z1, ctx, fz0)

               if(z1 <= y1) then
                  a0 = a1
                  b0 = b1
                  exit   ! to step 2
               end if

            end if
! step 4
            if(d1 <= eps) then
               golden_section_search_method = 0.5d0*(a1+b1)
               err = d1
               return
            else

               if(d1 <= alpha*d0) then
                  a0 = a1
                  b0 = b1
                  y0 = y1
                  z0 = z1
                  d1 = d2
                  d2 = d3
                  cycle ! to step 3
               else  ! d1 > eps .and. d1 > alpha**d0)
                  a0 = a1
                  b0 = b1
                  exit ! to step 2
               end if

            end if

         end do

      end do

   end function golden_section_search_method
!-----------------------------------------------------------------------
end module binary_tools
