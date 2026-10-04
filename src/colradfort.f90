module colradfort
!module to pretty much call everything.
   use types
   use atomicdata_module      
   use readadf04_module
   use crm_module
   use interpolation_module
   use plasma_module
   use omp_lib
   use input
   use sorting
   use constants_module
   implicit none
   integer                :: thrid
   integer                :: ierr, i, numTempsReq
   real(f64)              :: temp
   real(f64)              :: plt, pltnosob
   real(f64), allocatable :: tempsReq(:)
   real(f64), allocatable :: upsInterp(:)
   real(f64), allocatable :: cascade(:)
   real(f64), allocatable :: pec(:), crm(:, :), col1(:), pops_old(:)
   real(f64), allocatable :: popcoronal(:)
   real(f64), allocatable :: popsnosob(:), pecnosob(:)
   real(f64), allocatable :: sob_old(:)
   !
   real(f64), allocatable :: wavelengthforspectrum(:)
   real(f64), allocatable :: broadspec(:)
   !
   real(f64)              :: atomicDensity, numions
   ! shell boundaries (in units of c) used by getAtomicDensityLocal
   real(f64)              :: shellVelocityOuterC = 0.0_f64, shellVelocityInnerC = 0.0_f64

   integer                :: sob_iter
   real(f64)              :: sob_change, beta_change, beta_change_old = 1.d6
   logical                :: converged
   integer                :: k, j, p, l, ll
   real(f64)              :: t1, t2
   character(len=300)     :: broadmodedefault = 'gaussian'
   integer*8              :: shellnumtemp = 0
   
   real(f64), allocatable :: d_sob(:), d_sob_old(:)
   integer                :: densindex
   contains

   subroutine getadf04
      implicit none
      call cpu_time(t1)
      if (floersHack) then
         call readhack(trim(adf04Path))
      else
         call readadf04(trim(adf04Path))
      end if
      call cpu_time(t2)
      write (*, '(A,ES10.4,A)') '  [timing] adf04 read time                : ', t2 - t1, ' s'
   end subroutine

   subroutine getAtomicDensityLocal()
      !
      ! Atomic number density and number of ions in a (shell of a) homologously
      ! expanding ejecta. All inputs are module variables:
      !
      !   from module input:
      !     massElementSolar         mass of the element (Msun)
      !     fractionOverride         if > 0, atomicDensity = fractionOverride * density
      !     timeSinceExplosionDays   time since explosion (days)
      !     density                  electron density (cm-3); only used if fractionOverride > 0
      !   from this module:
      !     atomicNumber             from the adf04 file
      !     shellVelocityOuterC      outer shell velocity / c
      !     shellVelocityInnerC      inner shell velocity / c (0 for a full sphere)
      !
      ! Outputs (module variables): atomicDensity (cm-3) and numions.
      implicit none
      real(f64)            :: expansion_volume, time_exp_sec
      real(f64), parameter :: c_cgs = 3e10_f64

      time_exp_sec = timeSinceExplosionDays*86400.0_f64

      !total volume.
      expansion_volume = piFourOnThree*(shellVelocityOuterC*c_cgs*time_exp_sec)**3
      expansion_volume = expansion_volume - piFourOnThree*(shellVelocityInnerC*c_cgs*time_exp_sec)**3
      numions = massElementSolar*m_solar_grams/get_mass_grams(atomicNumber)

      atomicDensity = numions/(expansion_volume)

      if (fractionOverride > 0.0_f64) atomicDensity = fractionOverride*density
      print *, 'atomic number density', atomicDensity, 'cm-3 from new routine. edense=', density
   end subroutine

   subroutine colrad()
      !
      ! Runs the collisional-radiative calculation for one set of plasma
      ! conditions and writes popData/pecData/spectrum files.
      !
      ! All inputs are taken from module variables, which the caller must set
      ! beforehand (and call alloc(numwl) once so the output arrays exist):
      !
      !   from module input:
      !     temperature              plasma temperature (K)
      !     density                  electron density (cm-3)
      !     sobolev                  apply Sobolev escape probabilities
      !     timeSinceExplosionDays   time since explosion (days)
      !     wlmin_nm, wlmax_nm       spectrum wavelength range (nm)
      !     numwl                    number of spectrum wavelength points
      !     careful_la               use the careful linear-algebra solver
      !     writeoutrates            write out rates while building the CRM
      !     velocityExpansionC       expansion velocity / c
      !   from this module:
      !     atomicDensity            atomic number density (cm-3)
      !     numions                  number of ions
      !     broadmodedefault         line-broadening mode ('gaussian', 'box')
      !
      ! Outputs are left in the module arrays wavelengthforspectrum (cm) and
      ! broadspec.
      implicit none

      real(f64)            :: dwl
      character(len=20)    :: filesuffix
      integer, allocatable :: pecPointer(:)
      shellVelocityOuterC = velocityExpansionC
      shellVelocityInnerC = 0.0_f64
      call getAtomicDensityLocal
      call prepare_sobfactors(atomicDensity, timeSinceExplosionDays)
      if (.not. allocated(wavelengthforspectrum) .or. .not. allocated(broadspec)) then
         stop 'colrad: spectrum arrays not allocated - call alloc(numwl) first'
      end if
      if (size(wavelengthforspectrum) /= numwl .or. size(broadspec) /= numwl) then
         stop 'colrad: numwl does not match size of spectrum arrays allocated by alloc'
      end if

      tempsReq(1) = temperature
      shellnumtemp = shellnumtemp + 1
      i = 1

      call interpolate_upsilons_calc_rates(temperature)

      sob = 1.0_f64

      call cpu_time(t1)

      call solve_cr_with_continuity(numLevels,density, crm, col1, ierr)

      call BoltzmanPopulation(numlevels, statweight, energies, tempsReq(i), popcoronal)

      call cpu_time(t2)

      write (*, '(A,ES10.4,A)') '  [timing] initial populations : ', t2 - t1, ' s'
      
      converged = .false.
      sob_old   = 1.0_f64
      popsnosob = col1
      call calculate_pec_plt(numLevels, col1, ntran, aval, sob, pec, plt, density, energies)

      if (sobolev) call convergeSobolev(density)

      call cpu_time(t1)
      call calculate_pec_plt(numLevels, col1, ntran, aval, sob, pec, plt, density, energies)
      call cpu_time(t2)
      write (*, '(A,ES10.4,A)') '  [timing] PEC/PLT calculation : ', t2 - t1, ' s'

      call calculate_total_radiative_cascade(numlevels, ntran, aval, cascade)

      if (shellnumtemp < 10) then
         write (filesuffix, '(2I1)') 0, shellnumtemp
      else
         write (filesuffix, '(I2)') shellnumtemp
      end if

      if (mode .eq. 'astro') filesuffix(:) = ''

      open (100, file='popData'//trim(filesuffix))
      do j = 1, numlevels
         write (100, '(I4,3ES11.4)') j, col1(j), popcoronal(j), cascade(j) !, popcoronal(j)
      end do
      close (100)

      open (100, file='pecData'//trim(filesuffix))
      write (100, *) ' Low, Upp,     Sob,    aval,     pec,         wlcm,    popL,    popU,'

      if (sortpec) then
         allocate (pecPointer(size(pec)))
         do j = 1, size(pec)
            pecPointer(j) = j
         end do

         call qsort(pec, size(pec), pecPointer)

         do j = size(pec), 1, -1
            call inverseupperTriangleIndexing(pecPointer(j), numLevels, k, p)
                                write(100,'(2I5,3ES9.2,ES14.7,3ES9.2)')k,p,sob(pecPointer(j)), aval(pecPointer(j)), pec(j), wl_cm(pecPointer(j)), col1(k), col1(p),upsInterp(pecPointer(j))
         end do

         deallocate (pecPointer)

      else

         do j = 1, numLevels - 1
            do k = j + 1, numLevels
               p = upperTriangleIndexing(j, k, numLevels)
               write (100, '(2I5,3ES9.2,ES14.7,3ES9.2)') j, k, sob(p), aval(p), pec(p), wl_cm(p), col1(j), col1(k), upsInterp(p)
            end do
         end do
      end if

      close (100)

      call cpu_time(t1)
      broadspec(:) = 0.0d0
      wavelengthforspectrum(1) = wlmin_nm*1e-7
      wavelengthforspectrum(numwl) = wlmax_nm*1e-7
      dwl = 1e-7*(wlmax_nm - wlmin_nm)/(numwl - 1)
      do j = 2, numwl - 1
         wavelengthforspectrum(j) = wavelengthforspectrum(j - 1) + dwl
      end do
      call broadenedSpectrum(size(wavelengthforspectrum),wavelengthforspectrum,velocityExpansionC,broadspec,ntran,pec,wl_cm,density,numions,broadmodedefault)
      call cpu_time(t2)
      write (*, '(A,ES10.4,A)') '  [timing] spectrum broadening : ', t2 - t1, ' s'

      call cpu_time(t1)
      open (101, file='spectrum'//trim(filesuffix))
      do j = 1, size(wavelengthforspectrum)
         write (101, *) wavelengthforspectrum(j), broadspec(j)
      end do
      close (101)
      call cpu_time(t2)
      write (*, '(A,ES10.4,A)') '  [timing] spectrum write        : ', t2 - t1, ' s'

      close (1)

   end subroutine

   subroutine tempDensScan 
      implicit none 
      real(f64) :: tempGrid( 3)
      real(f64) :: densGrid(100)
      real(f64),allocatable :: popswithsob(:)
      real(f64) :: t1,t2
      !call cpu_time(t1)
      t1 = omp_get_wtime()
      shellVelocityOuterC = velocityExpansionC
      shellVelocityInnerC = 0.0_f64
      call getAtomicDensityLocal
      call prepare_sobfactors(atomicDensity, timeSinceExplosionDays)
      tempGrid(:) = (/(i * 1000       , i=1, 3, 1)/)
      densGrid(:) = (/(10.0_f64 ** (real(i)/real(10))  , i=1,100, 1)/)
      sob_old   = 1.0_f64
      sob       = 1.0_f64

      allocate(popswithsob(numLevels))

      do i = 1, size(tempGrid)
         call interpolate_upsilons_calc_rates(tempGrid(i))
         do densindex = 1, size(densGrid)
            !
            call solve_cr_with_continuity(numLevels,densGrid(densindex), crm, col1, ierr,useSob=.false.)
            if ( (i == 1) .and. (densindex==1)) popswithsob = col1 
            call calculate_pec_plt(numLevels, col1, ntran, aval, sob, pec, plt, densGrid(densindex), energies,useSob=.false.)
            col1 = popswithsob
            if (sobolev) then 
!               call BoltzmanPopulation(numLevels,statweight,energies, tempGrid(i),col1)
               call convergeSobolev(densGrid(densindex))
               !call newtonSobolev(densGrid(densindex))
               call calculate_pec_plt(numLevels, col1, ntran, aval, sob, pec, plt, densGrid(densindex), energies,useSob=.true.)
               popswithsob = col1 
            end if 
            write(50,'(2ES10.3, I4)')  pltnosob,plt,sob_iter
         end do       
      end do
      t2 = omp_get_wtime()
      write(50,*) '#',t2-t1
      write(50,*) '#',latime
   end subroutine

   subroutine newtonSobolev(electron_density)
      implicit none
      real(f64) :: electron_density

      integer,   parameter :: max_newton = 40
      real(f64), parameter :: newton_tol = 1.0e-8_f64   ! on the full Newton step; cheap, convergence is quadratic
      real(f64), parameter :: pop_floor  = 1.0e-10_f64
      real(f64), parameter :: max_drop   = 0.9_f64      ! no level may lose more than 90% per step
      logical,   parameter :: use_continuation = .false. ! tau-scaling 0.01 -> 0.1 -> 1 (try if it stalls)
      logical,   parameter :: fallback   = .true.        ! on failure, run the damped convergeSobolev
      logical,   parameter :: check_jac  = .false.       ! one-off finite-difference test of J
      logical,   parameter :: verify     = .false.        ! one Picard solve at the end as a consistency check
      logical,   parameter :: debug_sob  = .false.

      real(f64), allocatable :: n(:), nt(:), dn(:), F(:), Ft(:), fscale(:), Q(:,:), J(:,:), chk(:)
      integer,   allocatable :: ipiv(:)
      logical,   allocatable :: lmask(:)
      real(f64) :: stages(3), tsec, ts, tol, alpha, Fnorm, Ftnorm, step_rel
      real(f64) :: tt0, tt1, t_asm, t_sol
      integer   :: nstage, istage, it, ils, info, ilev, ntot, nasm
      logical   :: stage_ok

      allocate (n(numLevels), nt(numLevels), dn(numLevels), F(numLevels), Ft(numLevels))
      allocate (fscale(numLevels), Q(numLevels,numLevels), J(numLevels,numLevels))
      allocate (chk(numLevels), ipiv(numLevels), lmask(numLevels))

      tsec      = timeSinceExplosionDays*86400.0_f64
      converged = .false.
      ntot      = 0
      nasm      = 0
      t_asm     = 0.0_f64
      t_sol     = 0.0_f64

      call cpu_time(t1)
      popsnosob = col1
      pecnosob  = pec
      pltnosob  = plt

      n = max(col1, tiny(1.0_f64))
      n = n/sum(n)

      if (use_continuation) then
         nstage = 3
         stages = [0.01_f64, 0.1_f64, 1.0_f64]
      else
         nstage = 1
         stages = 1.0_f64
      end if

      if (check_jac) call jac_check(n, stages(nstage))

      stage_loop: do istage = 1, nstage
         ts  = stages(istage)
         tol = merge(newton_tol, 1.0e-3_f64, istage == nstage)
         stage_ok = .false.

         ! F, fscale and J at the current n; afterwards they are carried over from the line search
         call assemble(n, ts, F, J, .true.)

         newton_loop: do it = 1, max_newton
            ntot  = ntot + 1
            Fnorm = maxval(abs(F)/fscale)

            dn = -F
            call cpu_time(tt0)
            call dgesv(numLevels, 1, J, numLevels, ipiv, dn, numLevels, info)
            call cpu_time(tt1)
            t_sol = t_sol + (tt1 - tt0)
            if (info /= 0) then
               write (*, '(A,I6)') ' [newton] singular Jacobian, info =', info
               exit stage_loop
            end if

            lmask    = n > pop_floor*maxval(n)
            step_rel = maxval(abs(dn)/max(n, tiny(1.0_f64)), lmask)

            ! Backtracking line search on the scaled residual; componentwise limit on drops.
            ! Each trial also builds the Jacobian, so the accepted point needs no further assembly.
            alpha = 1.0_f64
            do ils = 1, 20
               nt = max(n + alpha*dn, (1.0_f64 - max_drop)*n)
               call assemble(nt, ts, Ft, J, .true.)
               Ftnorm = maxval(abs(Ft)/fscale)
               if (Ftnorm < (1.0_f64 - 1.0e-4_f64*alpha)*Fnorm) exit
               alpha = 0.5_f64*alpha
            end do
            n = nt
            F = Ft                           ! fscale and J already correspond to nt

            if (debug_sob) then
               ilev = maxloc(abs(dn)/max(n, tiny(1.0_f64)), 1, lmask)
               write (0, '(A,F5.2,A,I3,A,ES10.3,A,ES10.3,A,F7.4,A,I5)') ' newton ts=', ts, ' it=', it, &
                  ' |F|=', Ftnorm, ' dn/n=', step_rel, ' alpha=', alpha, ' lev=', ilev
            end if

            if (step_rel < tol) then
               stage_ok = .true.
               exit newton_loop
            end if
         end do newton_loop

         if (.not. stage_ok) exit stage_loop
         if (istage == nstage) converged = .true.
      end do stage_loop

      sob_iter = ntot

      if (.not. converged .and. fallback) then
         write (0, '(A,I4,A,I4)') ' [newton] failed; falling back to damped iteration, temp index ', i,' dens index',densindex
         col1 = popsnosob
         call convergeSobolev(electron_density)
         return
      end if

      ! make col1 and sob mutually consistent: sob = G(n)
      col1 = n
      call sobolev_escape(numLevels, ntran, aval, timeSinceExplosionDays, col1, &
                          statweight, wl_cm_cubed, atomicDensity)

      if (verify) then
         chk = col1
         call solve_cr_with_continuity(numLevels, electron_density, crm, chk, ierr, useSob=.true.)
         step_rel = maxval(abs(chk - n)/max(n, tiny(1.0_f64)), n > pop_floor*maxval(n))
         write (*, '(A,I4,A,ES10.3)') ' [newton] iterations:', ntot, '   verification dPop/Pop:', step_rel
      end if

      call cpu_time(t2)
      write (*, '(A,ES10.4,A)') '  [timing] Sobolev (Newton)         : ', t2 - t1, ' s'
      write (*, '(A,ES10.3,A,I4,A,ES10.3,A)') '  [timing]   assemble: ', t_asm, ' s (', nasm, ' calls),  dgesv: ', t_sol, ' s'
      if (.not. converged) write (*, '(A,I4)') 'WARNING: Newton did not converge for temp index ', i

   contains

      subroutine assemble(x, tscale, Fx, Jx, wantJ)
         ! Builds Q exactly as build_crm does (with beta evaluated from x), then
         !   Fx = Q x with row 1 replaced by sum(x) - 1
         !   Jx = dFx/dx (only when wantJ; otherwise Jx is left untouched)
         real(f64), intent(in)    :: x(numLevels), tscale
         real(f64), intent(inout) :: Fx(numLevels), Jx(numLevels,numLevels)
         logical,   intent(in)    :: wantJ
         real(f64) :: c, tau, beta, dbeta, em, w, ta0, ta1
         integer   :: ii, jj, kk

         call cpu_time(ta0)

         Q = electron_density*qmatrix
         if (wantJ) Jx = 0.0_f64

         kk = 1
         do ii = 1, numLevels - 1
            do jj = ii + 1, numLevels
               c   = tscale*sobconst*aval(kk)*wl_cm_cubed(kk)*statweight(jj)*atomicDensity*tsec
               tau = c*(x(ii)/statweight(ii) - x(jj)/statweight(jj))
               if (tau > 1.0e-5_f64) then          ! same branch as sobolev_escape (beta = 1 otherwise)
                  em   = exp(-tau)
                  beta = (1.0_f64-exp(-tau))/tau
                  if (tau < 1.0e-3_f64) then
                     dbeta = -0.5_f64 + tau/3.0_f64 - tau*tau/8.0_f64
                  else
                     dbeta = (em*(1.0_f64 + tau) - 1.0_f64)/(tau*tau)
                  end if
               else
                  beta  = 1.0_f64
                  dbeta = 0.0_f64
               end if
               Q(ii,jj) = Q(ii,jj) + aval(kk)*beta
               if (wantJ .and. abs(dbeta) > 0.0_f64) then
                  w = aval(kk)*x(jj)*dbeta*c       ! A * n_u * dbeta/dtau * c
                  Jx(ii,ii) = Jx(ii,ii) + w/statweight(ii)
                  Jx(ii,jj) = Jx(ii,jj) - w/statweight(jj)
                  Jx(jj,ii) = Jx(jj,ii) - w/statweight(ii)
                  Jx(jj,jj) = Jx(jj,jj) + w/statweight(jj)
               end if
               kk = kk + 1
            end do
         end do
         do jj = 1, numLevels                      ! loss conservation, as in build_crm
            Q(jj,jj) = -sum(Q(:,jj))
         end do

         call dgemv('N', numLevels, numLevels, 1.0_f64, Q, numLevels, x, 1, 0.0_f64, Fx, 1)

         fscale = tiny(1.0_f64)                    ! residual scale: total flux through each level,
         do jj = 1, numLevels                      ! accumulated column by column (unit stride)
            fscale = fscale + abs(Q(:,jj))*abs(x(jj))
         end do

         if (wantJ) Jx = Jx + Q

         Fx(1)     = sum(x) - 1.0_f64              ! continuity row
         fscale(1) = 1.0_f64
         if (wantJ) Jx(1,:) = 1.0_f64

         call cpu_time(ta1)
         t_asm = t_asm + (ta1 - ta0)
         nasm  = nasm + 1
      end subroutine assemble

      subroutine jac_check(x, tscale)
         real(f64), intent(in) :: x(numLevels), tscale
         real(f64), allocatable :: F0(:), F1(:), xp(:), Jfd(:)
         real(f64) :: h
         integer   :: k
         allocate (F0(numLevels), F1(numLevels), xp(numLevels), Jfd(numLevels))
         call assemble(x, tscale, F0, J, .true.)
         do k = 1, numLevels, max(1, numLevels/8)
            h  = 1.0e-6_f64*max(x(k), 1.0e-12_f64)
            xp = x;  xp(k) = xp(k) + h
            call assemble(xp, tscale, F1, J, .false.)       ! leaves J (analytic) untouched
            Jfd = (F1 - F0)/h
            write (0, '(A,I5,A,ES10.3)') ' jac check col', k, '  max rel err:', &
               maxval(abs(Jfd - J(:,k)))/(maxval(abs(J(:,k))) + tiny(1.0_f64))
         end do
      end subroutine jac_check

   end subroutine newtonSobolev

   subroutine convergeSobolev(electron_density)
      implicit none
      real(f64) :: electron_density

      real(f64), parameter :: d_min      = 2.0e-3_f64
      real(f64), parameter :: shrink     = 0.5_f64     ! d -> shrink*d on a non-decaying sign flip
      real(f64), parameter :: grow       = 1.5_f64     ! recovery toward the per-line cap
      real(f64), parameter :: cap_relax  = 1.02_f64    ! cap creeps back up 2% per iteration
      real(f64), parameter :: flip_ratio = 0.6_f64     ! a flip with |f| < ratio*|fprev| is decaying: ignore
      real(f64), parameter :: pop_floor  = 1.0e-10_f64
      real(f64), parameter :: flux_floor = 1.0e-6_f64
      real(f64), parameter :: res_factor = 10.0_f64
      logical,   parameter :: warm_start  = .true.     ! start from the previous grid point's beta
      logical,   parameter :: require_res = .false.    ! also demand the masked undamped residual be small
      logical,   parameter :: debug_sob   = .false.

      real(f64), allocatable, save :: x_keep(:), cap_keep(:)
      real(f64), allocatable :: x(:), f(:), fprev(:), d(:), dcap(:), eff(:), col_old(:), dpop(:)
      logical,   allocatable :: lev_mask(:), tr_mask(:), flip(:)
      real(f64) :: pop_change, res, f_floor
      integer   :: ilev
      logical   :: warm, conv
      logical   :: aggressiveStart 

      allocate (x(ntran), f(ntran), fprev(ntran), d(ntran), dcap(ntran), eff(ntran), flip(ntran), tr_mask(ntran))
      allocate (col_old(size(col1)), dpop(size(col1)), lev_mask(size(col1)))

      f_floor    = 0.1_f64*sob_tol
      fprev      = 0.0_f64
      converged  = .false.
      pop_change = huge(1.0_f64)
      res        = huge(1.0_f64)

      warm = warm_start .and. (i > 1)
      if (warm) warm = allocated(x_keep)
      if (warm) warm = (size(x_keep) == ntran)

      call cpu_time(t1)
      popsnosob = col1
      pecnosob  = pec
      pltnosob  = plt

      if (warm) then
         x    = x_keep
         dcap = min(1.0_f64, 4.0_f64*cap_keep)
      else
         call sobolev_escape(numLevels, ntran, aval, timeSinceExplosionDays, col1, &
                             statweight, wl_cm_cubed, atomicDensity)
         x    = log(sob)
         dcap = 1.0_f64
      end if

      !aggressiveStart = .true.
      !if (aggressiveStart) then 
      !   !tr_mask = aval > 1e3 
      !   sob(:) = 1.0_f64
      !   do sob_iter = 1,ntran 
      !      if (aval(sob_iter) > 1e3 ) sob(sob_iter) = 1e-3
      !   end do
      !   x    = log(sob)
      !   !write(0,*) 'aggressive start'
      !   call solve_cr_with_continuity(numLevels, electron_density, crm, col1, ierr, useSob=.true.)
      !   d = 1.0_f64
      !end if 


      d = dcap

      sob_iter_loop: do sob_iter = 1, max_sob_iter

         sob = exp(x)                            ! beta used in this solve
         col_old = col1
         call solve_cr_with_continuity(numLevels, electron_density, crm, col1, ierr, useSob=.true.)

         lev_mask   = col1 > pop_floor*maxval(col1)
         dpop       = abs(col1 - col_old)/max(col1, tiny(1.0_f64))
         pop_change = maxval(dpop, lev_mask)

         call sobolev_escape(numLevels, ntran, aval, timeSinceExplosionDays, col1, &
                             statweight, wl_cm_cubed, atomicDensity)      ! sob <- G(x)
         f = log(sob) - x                                                 ! undamped residual in ln(beta)

         eff     = sob_weight*sob
         tr_mask = eff > flux_floor*maxval(eff)
         if (.not. any(tr_mask)) tr_mask = .true.
         res = maxval(abs(1.0_f64 - exp(-f)), tr_mask)

         if (debug_sob) then
            ilev = maxloc(dpop, 1, lev_mask)
            write (0, '(A,I4,A,ES10.3,A,I6,A,ES10.3,A,I7,A,ES9.2)') ' it=', sob_iter, &
               ' dpop=', pop_change, ' (lev ', ilev, ') res=', res, &
               ' n(d<0.5)=', count(d < 0.5_f64), ' dmin=', minval(d)
         end if

         conv = (sob_iter > 1) .and. (pop_change < sob_tol)
         if (require_res) conv = conv .and. (res < res_factor*sob_tol)
         if (conv) then
            converged = .true.
            sob = exp(x)                         ! consistent with the col1 just solved
            write (*, '(A,I4)')     ' [sobolev] converged at iter   :', sob_iter
            write (*, '(A,ES10.4)') '        with maximum dPop/Pop   : ', pop_change
            exit sob_iter_loop
         end if

         ! Per-line damping: shrink only on a sign flip that is not decaying
         if (sob_iter > 1) then
            flip = (f*fprev < 0.0_f64) .and. (abs(f) > f_floor) .and. (abs(f) > flip_ratio*abs(fprev))
            where (flip)
               dcap = max(d_min, shrink*d)
               d    = dcap
            elsewhere
               dcap = min(1.0_f64, cap_relax*dcap)
               d    = min(dcap, grow*d)
            end where
         end if

         fprev = f
         x = min(x + d*f, 0.0_f64)               ! beta <= 1

      end do sob_iter_loop

      if (allocated(x_keep)) then
         if (size(x_keep) /= ntran) deallocate (x_keep, cap_keep)
      end if
      if (.not. allocated(x_keep)) allocate (x_keep(ntran), cap_keep(ntran))
      x_keep   = x
      cap_keep = dcap

      call cpu_time(t2)
      write (*, '(A,ES10.4,A)') '  [timing] Sobolev iteration        : ', t2 - t1, ' s'
      if (.not. converged) then
         write (*, '(A,I4,A,I3,A,2ES10.2)') &
            'WARNING: Sobolev did not converge for temp index ', i, &
            ' after ', max_sob_iter, ' iterations; dPop, res:', pop_change, res
      end if

      deallocate (x, f, fprev, d, dcap, eff, flip, tr_mask, col_old, dpop, lev_mask)
   end subroutine convergeSobolev

   subroutine levelscan()
      implicit none 
      ! Uses module variables temperature, density, careful_la, writeoutrates (from input).
      real(f64), allocatable :: crm_copy(:, :), col1_copy(:)

      !if ( allocated(crm) ) deallocate(crm)
      !allocate(crm(numlevels - 1, numlevels - 1))


      tempsReq(1) = temperature
      allocate (crm_copy(numlevels, numlevels))
      sob = 1.0_f64
      call interpolate_upsilons_calc_rates(temperature)

      call build_crm(numLevels, density, crm)

      !call build_cr_matrix(numLevels, ntran, statweight, energies, &
      !                     upsInterp, aval, sob, tempsReq(1), density, crm, col1, ierr, writeoutrates)

      crm_copy(:, :) = crm(:, :)
      col1_copy = col1
      open (32, file='plt_level_convergence.dat')
      do i = 2, numlevels

         !call solve_cr_populations_axb(numLevels, crm, i, col1, ierr, careful_la)
         call solve_cr_with_continuity(numLevels, density, crm, col1, ierr, skipbuild=1,ninclude=i)
         call calculate_pec_plt(numLevels, col1, ntran, aval, sob, pec, plt, density, energies)
         write (32, *) i, plt, col1(2), col1(1)

         crm(:, :) = crm_copy(:, :)
         col1(:) = col1_copy(:)
      end do
      !open(100,file='popData')
      ! do j = 1,numlevels
      !         write(100,'(I4,2ES11.4)') j ,col1(j) !, popcoronal(j)
      ! end do
      !close(100)

      deallocate (crm_copy)
   end subroutine

   subroutine masscontour()
      ! Uses module variables temperature, density, requiredLumo, careful_la,
      ! writeoutrates and verbose (from input).
      use input, only: contourLower, contourUpper
      implicit none
      real(f64) :: electronDensityLocalvary(1000)
      real(f64) :: temperaturevary(1000)
      real(f64) :: thislumo_per_ion
      real(f64) :: num_req, mass_req
      real(f64) :: num_in_one_solar_mass
      integer   :: ii, jj, counterii = 0, counterjj = 0
      real(f64) :: xx

      num_in_one_solar_mass = 1.0*m_solar_grams/get_mass_grams(atomicnumber)

      sob = 1

      open (90, file='contour.out')

      !get central estimate
      call interpolate_upsilons_calc_rates(temperature)
      call getmassestimate(density, mass_req, num_req, &
                           thislumo_per_ion, &
                           num_in_one_solar_mass)

      write (90, '(A, ES14.6,A)') '# Central temp         = ', temperature, ' Kelvin'
      write (90, '(A, ES14.6,A)') '# Central dens         = ', density, ' /cm3'
      write (90, '(A, ES14.6,A)') '# Central estimate = ', mass_req, ' Msun'
      write (90, '(A, I3)') '# Atomic        number         = ', atomicnumber
      write (90, '(A, I3,A)') '# Atomic        charge         = ', ioncharge_plus, ' +'
      write (90, '(A, ES14.6,A)') '# Central estimate = ', energies(contourUpper) - energies(contourLower), ' eV'

      write (90, '(A, I10)') '# ntemp = ', size(temperaturevary)
      write (90, '(A, I10)') '# ndens = ', size(electronDensityLocalvary)

      xx = log10(temps(size(temps))/temps(1))
      xx = 10**(xx/size(temperaturevary))

      temperaturevary(1) = temps(1)
      do ii = 2, size(temperaturevary)
         temperaturevary(ii) = temperaturevary(ii - 1)*xx
      end do
      temperaturevary(size(temperaturevary)) = temps(numtemps)
      !density grid
      electronDensityLocalvary(1) = 3.0
      electronDensityLocalvary(size(electronDensityLocalvary)) = 13.0
      xx = (electronDensityLocalvary(size(electronDensityLocalvary)) - electronDensityLocalvary(1))/size(electronDensityLocalvary)
      do ii = 2, size(electronDensityLocalvary)
         electronDensityLocalvary(ii) = electronDensityLocalvary(ii - 1) + xx
      end do

      electronDensityLocalvary(:) = 10**electronDensityLocalvary(:)

      if (.not. verbose) then

         write (90, '(A)') '# temp vary'

         !vary density
         do ii = 1, size(temperaturevary)
            call interpolate_upsilons_calc_rates(temperaturevary(ii))
            call getmassestimate(density, mass_req,  num_req, &
                                 thislumo_per_ion, &
                                 num_in_one_solar_mass)
            write (90, '(2ES14.6)') mass_req, temperaturevary(ii)
         end do

         write (90, '(A)') '# dens vary'

         call interpolate_upsilons_calc_rates(temperature)
         do ii = 1, size(electronDensityLocalvary)
            call getmassestimate(electronDensityLocalvary(ii), mass_req, num_req, &
                                 thislumo_per_ion, &
                                 num_in_one_solar_mass)
            write (90, '(2ES14.6)') mass_req, electronDensityLocalvary(ii)
         end do

      else

         open (25, file='tempgrid')
         open (26, file='densgrid')
         do jj = 1, size(temperaturevary), 10
            write (25, '(1ES14.6)') temperaturevary(jj)
         end do

         do ii = 1, size(electronDensityLocalvary), 10
            write (26, '(1ES14.6)') electronDensityLocalvary(ii)
         end do

         close (26)
         close (25)

         do jj = 1, size(temperaturevary), 10
            counterjj = counterjj + 1
            counterii = 0
      call interpolate_upsilons_calc_rates(temperaturevary(jj))
            do ii = 1, size(electronDensityLocalvary), 10
               counterii = counterii + 1
             call getmassestimate(electronDensityLocalvary(ii), mass_req, num_req, &
                                    thislumo_per_ion, &
                                    num_in_one_solar_mass)
               write (90, '(2I10,1ES14.6)') counterii, counterjj, mass_req
            end do
         end do

      end if

      close (90)

   end subroutine masscontour

   subroutine lineplot()
      ! Uses module variable requiredLumo (from input).
      use input, only: contourLower, contourUpper
      implicit none
      real(f64) :: electronDensityLocalvary(200)
      real(f64) :: temperaturevary(3)
      real(f64) :: thislumo_per_ion
      real(f64) :: num_req, mass_req
      real(f64) :: num_in_one_solar_mass
      real(f64) :: massdump(200)
      integer :: ii, jj
      real(f64) :: xx
      ! deliberately local (shadow the input values): lineplot always runs with these off
!      logical :: careful_la = .false., writeoutrates = .false.

      num_in_one_solar_mass = 1.0*m_solar_grams/get_mass_grams(atomicnumber)
      sob = 1
      xx = log10(temps(size(temps))/temps(1))
      xx = 10**(xx/size(temperaturevary))
      !temperaturevary(1) = temps(1)
      !do ii =         2, size(temperaturevary)
      !         temperaturevary(ii) = temperaturevary(ii-1) * xx
      !end do

      temperaturevary(1) = 1000
      temperaturevary(2) = 3000
      temperaturevary(3) = 10000

      electronDensityLocalvary(1) = 4.0
      electronDensityLocalvary(size(electronDensityLocalvary)) = 7.0
      xx = (electronDensityLocalvary(size(electronDensityLocalvary)) - electronDensityLocalvary(1))/size(electronDensityLocalvary)
      do ii = 2, size(electronDensityLocalvary)
         electronDensityLocalvary(ii) = electronDensityLocalvary(ii - 1) + xx
      end do

      electronDensityLocalvary(:) = 10**electronDensityLocalvary(:)

      open (88, file='lineplot.dat')

      write (88, '(F10.8)') wl_cm(upperTriangleIndexing(contourLower, contourUpper, numlevels))
      write (88, '(1000ES10.3)') temperaturevary(:)
      write (88, '(1000ES10.3)') electronDensityLocalvary(:)

      do jj = 1, size(temperaturevary)
         call interpolate_upsilons_calc_rates(temperaturevary(jj))
         do ii = 1, size(electronDensityLocalvary)
            call getmassestimate(electronDensityLocalvary(ii), mass_req, num_req, &
                                 thislumo_per_ion, &
                                 num_in_one_solar_mass)
            massdump(ii) = mass_req
         end do
         write (88, '(1000ES10.3)') massdump(:)

      end do

   end subroutine

   subroutine getmassestimate( &
      reqdens, &
      reqmass, &
      num_req, &
      thislumo_per_ion, &
      num_in_one_solar_mass)
      use input, only: contourLower, contourUpper
      implicit none
      real(f64) :: num_req, thislumo_per_ion, num_in_one_solar_mass
      real(f64) :: reqdens, reqmass
 !     logical :: careful_la, writeoutrates
      integer :: pp
      call solve_cr_with_continuity(numLevels,reqdens, crm, col1, ierr)
      pp = upperTriangleIndexing(contourLower, contourUpper, numlevels)
      write (1414, *) '   The line in contour is: ', wl_cm(pp), aval(pp)

      thislumo_per_ion = col1(contourUpper)*aval(pp)*hc_ergcm/wl_cm(pp)

      num_req = requiredlumo/thislumo_per_ion

      reqmass = num_req/num_in_one_solar_mass

   end subroutine
   subroutine convergeSobolevold(electron_density)
      implicit none
      real(f64) :: electron_density      
      real(f64) :: sob_damp 
      real(f64) :: avg 
      converged = .false.
      sob_damp = sob_damp_initial
      call cpu_time(t1)
      popsnosob = col1
      pecnosob  = pec
      pltnosob  = plt

      !write(0,*) 'initial approximation using  ',col1(1)
      !write(0,*) '----------------------------------------------------'
      call sobolev_escape(numLevels, ntran, aval, timeSinceExplosionDays, col1, &
                           statweight, wl_cm_cubed, atomicDensity)

      sob_iter_loop: do sob_iter = 1, max_sob_iter

         call solve_cr_with_continuity(numLevels,electron_density, crm, col1, ierr,useSob=.true.)

         !write(0,*) maxval(col1 - popsnosob)

         sob_old = sob
         avg = sum(sob) / ntran
         !write(0,'(ES12.5)') avg
         call sobolev_escape(numLevels, ntran, aval, timeSinceExplosionDays, col1, &
                              statweight, wl_cm_cubed, atomicDensity)

         sob = sob_damp*sob + (1.0_f64 - sob_damp) * sob_old

         !this is a fairly conservative convergence criterion - basically it asserts that
         !none of the beta's change by more than 0.1%, for sob_tol = 1e-3.
         beta_change = maxval(abs(sob - sob_old)/sob)
         !write(0,*) beta_change
         if (sob_iter > 1 .and. beta_change < sob_tol) then
            converged = .true.
            write (*, '(A,I4,A,ES10.3)') ' [sobolev] converged at iter   :', sob_iter
            write (*, '(A,ES10.4)') '        with maximum dBeta/Beta : ', beta_change
            exit sob_iter_loop
         end if

         !if (sob_iter > 1) then
         !   if (beta_change > beta_change_old) then
         !   ! Oscillating or diverging: reduce step size
         !      sob_damp = max(0.05_f64, sob_damp * 0.5_f64)
         !   else if (beta_change < 0.8_f64 * beta_change_old) then
         !! Monotonically converging well: gradually restore step size
         !      sob_damp = min(sob_damp_initial, sob_damp * 1.05_f64)
         !   end if
         !end if


         beta_change_old = beta_change

         if (mod(sob_iter,10) == 0) then 
            sob_damp = sob_damp * 0.5_f64
            !sob = sob_damp*sob + (1.0_f64 - sob_damp) * sob_old
         end if

      end do sob_iter_loop

      !write(0,*) 'converged approximation using',col1(1)
      !write(0,*) '-----------------------------------------------------'

      call cpu_time(t2)
      write (*, '(A,ES10.4,A)') '  [timing] Sobolev iteration        : ', t2 - t1, ' s'
      if (.not. converged) then
         write (*, '(A,I4,A,I3,A,2ES10.2)') &
            'WARNING: Sobolev did not converge for temp index ', i, &
            ' after ', max_sob_iter, ' iterations', beta_change, beta_change_old
         write(0,*) 'didnt converge'
      end if

   end subroutine
   subroutine alloc
      implicit none
      numTempsReq = 1
      call allocipiv
      allocate (crm(numLevels, numLevels))
      allocate (col1(numLevels), pops_old(numLevels))
      allocate (tempsReq(numTempsReq))
      allocate (upsInterp(ntran))
      allocate (pec(ntran))
      allocate (sob(ntran))
      allocate (sobcoefficient(ntran))
      allocate (sob_tau(ntran))
      allocate (sob_weight(ntran))
      allocate (sob_old(ntran))
      allocate (pecnosob(ntran))
      allocate (popcoronal(numlevels))
      allocate (cascade(numlevels))
      allocate (wavelengthforspectrum(numwl), broadspec(numwl))

   end subroutine

   subroutine dealloc
      call dealloc_atomicdata
      call deallocipiv
      if (allocated(tempsReq)) deallocate (tempsReq)
      if (allocated(upsInterp)) deallocate (upsInterp)
      if (allocated(crm)) deallocate (crm)
      if (allocated(col1)) deallocate (col1)
      if (allocated(pops_old)) deallocate (pops_old)
      if (allocated(popsnosob)) deallocate (popsnosob)
      if (allocated(pec)) deallocate (pec)
      if (allocated(pecnosob)) deallocate (pecnosob)
      if (allocated(sob)) deallocate (sob)
      if (allocated(sob_old)) deallocate (sob_old)
      if (allocated(wavelengthforspectrum)) deallocate (wavelengthforspectrum)
      if (allocated(broadspec)) deallocate (broadspec)
      if (allocated(cascade)) deallocate (cascade)

   end subroutine

end module colradfort