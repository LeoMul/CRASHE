module crm_module
   ! In this module, the reduced collisional radiative matrix is contructed.
   ! For an nlev system, we put everything relative to level 1, and contruct
   ! three arrays. The nlev-1 x nlev-1 reduced CRM, and the rows and columns
   ! corresponding to level 1. While mathematically they are both of dimension
   ! nlev - we only need from index 2 onwards - and only those elements are calculated.
   ! However - the array Qcol1 being the first column of the full matrix, is still
   ! allocated as nlev. This is because we will store the populations including the
   ! ground here later. The populations are normalized s.t sum(pops) = 1.

   !
   !
   use types
   use constants_module
   use Periodic_Table
   use readadf04_module, only: upperTriangleIndexing
   use atomicdata_module,only :qmatrix, aval,sob, numLevels
   implicit none
   integer, allocatable   :: ipiv(:)

contains
   subroutine deallocipiv
      implicit none
      if (allocated(ipiv)) deallocate(ipiv)
   end subroutine

   subroutine allocipiv 
      implicit none 
      call deallocipiv 
      allocate(ipiv(numlevels))
   end subroutine 



   subroutine coronalPopulation(nlev, ntran, g, E, Ups, Aval, Te, Ne, coronalPop)
      implicit none
      integer, intent(in)  :: nlev, ntran
      real(f64), intent(in)  :: g(nlev)
      real(f64), intent(in)  :: E(nlev)
      real(f64), intent(in)  :: Ups(ntran)
      real(f64), intent(in)  :: Aval(ntran)

      real(f64), intent(in)  :: Ne
      real(f64), intent(in)  :: Te
      real(f64)              :: coronalPop(nlev)

      integer  :: i, j, pp
      real(f64) :: kT, sqrt_Te, dE, q_exc, q_deexc, avsum

      coronalPop(1) = 1.0d0
      kT = kB_eV*Te
      sqrt_Te = sqrt(Te)

      do i = 2, nlev
         dE = E(i) - E(1)
         pp = upperTriangleIndexing(1, i, nlev)
         q_deexc = (coll_fac/(g(i)*sqrt_Te))*Ups(pp)

         if (dE/kT < 700.0_f64) then
            q_exc = (coll_fac/(g(1)*sqrt_Te))*Ups(pp)*exp(-dE/kT)
         else
            q_exc = 0.0_f64
         end if

         avsum = 0.0d0

         !print*,'1-->',i,q_exc,q_deexc,q_exc/q_deexc

         do j = 1, i - 1
            pp = upperTriangleIndexing(j, i, nlev)
            avsum = avsum + aval(pp)
         end do
         avsum = avsum + Ne*q_deexc
         coronalPop(i) = Ne*q_exc/avsum
      end do

      avsum = sum(coronalPop(:))
      coronalPop(:) = coronalPop(:)/avsum

   end subroutine

   subroutine BoltzmanPopulation(nlev, g, E, Te, boltzpop)
      implicit none
      integer, intent(in)  :: nlev
      real(f64), intent(in)  :: g(nlev)
      real(f64), intent(in)  :: E(nlev)

      real(f64), intent(in)  :: Te
      real(f64)              :: boltzpop(nlev)

      integer  :: i
      real(f64) :: kT, sqrt_Te, dE, avsum, w1, wi

      boltzpop(1) = 1.0d0
      kT = kB_eV*Te
      sqrt_Te = sqrt(Te)
      w1 = g(1)

      do i = 2, nlev
         wi = g(i)
         dE = E(i) - E(1)
         boltzpop(i) = (wi/w1)*exp(-de/kt)
         !print*,boltzpop(i)
      end do

      avsum = sum(boltzpop(:))
      !print*,avsum
      boltzpop(:) = boltzpop(:)/avsum

   end subroutine

   subroutine     build_crm(nlev, density, Q)
      implicit none
      integer,   intent(in)    :: nlev
      real(f64)                :: density
      real(f64), intent(inout) :: Q(nlev, nlev)
      integer :: ii,jj,kk

      Q(:,:)  = density * qmatrix(:,:)
      kk=1
      do ii = 1, nlev-1
         do jj = ii+1, nlev
            !jj > ii 
            !Q(jj,ii) =  Q(jj,ii) + density * qmatrix(jj,ii)
            !Q(ii,jj) =  Q(ii,jj) + density * qmatrix(ii,jj) + aval(kk) * sob(kk)
            Q(ii,jj) =  Q(ii,jj) + aval(kk) * sob(kk)
            kk = kk + 1
         end do 
      end do
!
      !enforce loss conservation...
      do jj = 1, nlev
         Q(jj, jj) = -sum(Q(:, jj))
      end do

   end subroutine build_crm

   subroutine solve_cr_with_continuity(nlev, density, Q, pops, ierr,skipbuild,ninclude)
      implicit none
      integer,   intent(in)    :: nlev
      real(f64)                :: density
      real(f64), intent(inout) :: Q(nlev, nlev)
      real(f64), intent(out)   :: pops(nlev)
      
      integer, intent(in), optional :: skipbuild,ninclude
      integer :: skipbuildinternal = 0, nincludeInternal 
      integer, intent(out) :: ierr
      

      !if the CRM for this case has already been built, for some reason.
      if (present(skipbuild)) skipbuildinternal = skipbuild

      if (skipbuildinternal == 0) then 
         call build_crm(nlev, density, Q)
      end if 
      ! replace first row with continuity
      !this is easier to maintain than the old version.
      Q(1, :) = 1.0_f64
      pops(:) = 0.0_f64
      pops(1) = 1.0_f64
!     
!     If for some reason the user wants to override the number of levels
!     actually included in the CRM.
      nincludeInternal = nlev
      if (present(ninclude)) nincludeInternal   = ninclude
!
      call dgesv(nincludeInternal, 1, Q, nlev, ipiv, pops, nlev, ierr)
!      write (69, *) info, pops
   end subroutine

   subroutine calculate_total_radiative_cascade(nlev, ntran, avals, cascade)
      integer   :: nlev, ntran
      real(f64) :: avals(ntran), cascade(nlev)

      integer :: ii, jj, pp

      cascade(:) = 0.0_f64

      do ii = 2, nlev
         do jj = 1, ii - 1
            !print*,jj,ii
            pp = upperTriangleIndexing(jj, ii, nlev)
            cascade(ii) = cascade(ii) + avals(pp)
         end do
      end do

   end subroutine

   subroutine sobolev_escape(nlev, ntran, baseAvals, time_exp_days, pops, weights, wl_cm_cubed, atomicDensityLocal)
      !calculates Sobolev escape probability.
      implicit none
      integer :: nlev, ntran
      !
      real(f64) ::  pops(nlev), weights(nlev)
      real(f64) :: baseAvals(ntran), wl_cm_cubed(ntran)
      !
      real(f64) :: atomicDensityLocal
      real(f64) :: tau, time_exp_days, time_exp_sec
      real(f64), parameter :: c_cgs = 3e10_f64
      integer :: pp
      integer :: ii, jj
      !
      time_exp_sec = time_exp_days*86400.0_f64

      print *, 'number density', atomicDensityLocal, 'cm-3'

      sob(:) = 1.0_f64

      !write(0,*) 'sobconst',sobconst,atomicDensityLocal

      do ii = 1, nlev - 1
         do jj = ii + 1, nlev
            pp = upperTriangleIndexing(ii, jj, nlev)

            tau = sobconst * baseAvals(pp) * wl_cm_cubed(pp) * weights(jj) * atomicDensityLocal * time_exp_sec * (pops(ii)/weights(ii) -  pops(jj)/weights(jj))
            if (tau > 1.0e-5_f64) sob(pp) = (1.0_f64 - exp(-tau))/tau
            if ((tau < 0.0_f64)) write (999, *) ii, pops(ii), jj, pops(jj), baseAvals(pp), tau
         end do
      end do
      !write(0,*) minval(sob)
      !
   end subroutine

   subroutine broadenedSpectrum(numWavelengths, &
                                wavelength, &
                                velocityShell, &
                                spectra, &
                                ntran, &
                                pec, &
                                spectralLinesCM, &
                                electron_density, &
                                numIonsLocal, &
                                calcMode &
                                )
      !
      ! Calculates \sum_i nf * PEC * ProfileShape
      ! For Gaussian: nf = 1 / (sqrt(2pi) σ)
      ! For Box:      nf = 1 / (2 * Δλ_max)
      ! Note: An extra factor of (1 / wavelengthCentral) is included in normfactor
      ! to prepare for the E = hc/λ conversion applied at the end.
      !
      implicit none

      character*10, intent(in) :: calcmode
      integer, intent(in)      :: numWavelengths, ntran
      real(f64), intent(in)    :: wavelength(numWavelengths)
      real(f64), intent(in)    :: velocityShell, electron_density, numIonsLocal
      real(f64), intent(in)    :: pec(ntran), spectralLinesCM(ntran)
      real(f64), intent(inout) :: spectra(numWavelengths)

      real(f64) :: ww, wavelengthCentral
      real(f64) :: sig_cm, sigOneOver
      integer   :: ii, jj

      real(f64) :: normfactor
      real(f64) :: thispec, wl_lo, wl_hi, dwl
      real(f64), parameter :: nsigmacut = 4.0_f64
      integer   :: jlo, jhi
      real(f64) :: totalpec
      real(f64), parameter :: pecthreshold = 1e-4_f64
      real(f64) :: peccutoff
      real(f64) :: thisphotonenergy

      dwl = wavelength(2) - wavelength(1)
      !write(0,*) dwl
      totalpec = sum(pec)
      peccutoff = pecthreshold*totalpec

      ! Calculate spectrum based on the selected mode
      do ii = 1, ntran
         wavelengthCentral = spectralLinesCM(ii)

         if (wavelengthCentral < 1e-30_f64) cycle

         thispec = pec(ii)
         if (thispec < peccutoff) cycle

         thisphotonenergy = thispec
         !write(0,*) 'broad mode ', calcmode,trim(adjustl(calcMode)),trim(adjustl(calcMode))=='box'

         ! --- BRANCH: BOX PROFILE (Expanding Shell) ---
         if ((calcMode(1:3) == 'box') .or. trim(adjustl(calcMode)) == 'onion') then
            write (0, *) 'i am doing a box'
            ! For an expanding shell, max Doppler shift is defined by the shell velocity.
            ! Assuming velocityShell here represents v/c (or v_expansion / c).
            wl_lo = wavelengthCentral*(1.0_f64 - velocityShell)
            wl_hi = wavelengthCentral*(1.0_f64 + velocityShell)

            ! Find grid bins
            jlo = max(1, nint((wl_lo - wavelength(1))/dwl) + 1)
            jhi = min(numWavelengths, nint((wl_hi - wavelength(1))/dwl) + 1)

            ! The box height is 1 / Total Width.
            ! Extra (1 / wavelengthCentral) applied for E = hc/λ step.
            normfactor = (1.0_f64/(wl_hi - wl_lo))/wavelengthCentral
            ww = thisphotonenergy*normfactor

            ! Flat profile: Add uniform intensity to all bins within the box
            do jj = jlo, jhi
               spectra(jj) = spectra(jj) + ww
            end do

            ! --- BRANCH: GAUSSIAN PROFILE (Thermal/Microturbulence) ---
         else

            sig_cm = fwhmSigma*wavelengthCentral*velocityShell
            sigOneOver = 1.0_f64/sig_cm

            ! Gaussian normalization.
            ! Extra (1 / wavelengthCentral) applied for E = hc/λ step.
            normfactor = sigOneOver*oneOverSQRTTWOPI/wavelengthCentral

            wl_lo = wavelengthCentral - nSigmaCut*sig_cm
            wl_hi = wavelengthCentral + nSigmaCut*sig_cm

            jlo = max(1, nint((wl_lo - wavelength(1))/dwl) + 1)
            jhi = min(numWavelengths, nint((wl_hi - wavelength(1))/dwl) + 1)

            do jj = jlo, jhi
               ww = wavelength(jj)
               ww = (ww - wavelengthCentral)*sigOneOver
               ww = thisphotonenergy*normfactor*exp(minusHalf*ww*ww)
               spectra(jj) = spectra(jj) + ww
            end do

         end if

      end do

      ! Final scaling: Units conversion (Photons -> Ergs)
      spectra(:) = spectra(:)*(numIonsLocal*electron_density*hc_ergcm*1e-8_f64)

   end subroutine

end module crm_module