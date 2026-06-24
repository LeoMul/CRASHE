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
  use Periodic_Table
  use readadf04_module,only: upperTriangleIndexing
  implicit none
  real(f64), parameter :: kB_eV    = 8.617333262e-5_f64
  real(f64), parameter :: coll_fac = 8.629e-6_f64
  real(f64), parameter :: pi = 4.0_f64 * atan(1.0_f64)
  real(f64), parameter :: piFourOnThree = 4.0_f64 * pi / 3.0_f64
  real(f64), parameter :: sobconst = 1.0_f64 / (8.0 * pi)
  real(f64), parameter :: m_solar_grams = 1.989e+33_f64
  real(f64), parameter :: fwhmSigma = 1._f64/2.355_f64
  real(f64), parameter :: minusHalf = -0.5_f64
  real(f64), parameter :: oneOverSQRTTWOPI = 1._f64/(sqrt( 2.0_f64 * pi))
  real(f64), parameter :: hc_ergcm = 1.98644586e-16_f64 ! in erg cm
contains

  subroutine coronalPopulation(nlev, ntran, g, E, Ups, Aval, Te, Ne, coronalPop)
    implicit none
    integer,  intent(in)  :: nlev, ntran
    real(f64), intent(in)  :: g(nlev)
    real(f64), intent(in)  :: E(nlev)
    real(f64), intent(in)  :: Ups(ntran)
    real(f64), intent(in)  :: Aval(ntran)

    real(f64), intent(in)  :: Ne
    real(f64), intent(in)  :: Te
    real(f64)              :: coronalPop(nlev)

    integer  :: i, j, pp
    real(f64) :: kT, sqrt_Te, dE, q_exc, q_deexc,avsum


    coronalPop(1) = 1.0d0 
    kT      = kB_eV * Te
    sqrt_Te = sqrt(Te)

    do i = 2, nlev 
      dE  = E(i) - E(1)
      pp = upperTriangleIndexing(1,i,nlev)
      q_deexc = (coll_fac / (g(i) * sqrt_Te)) * Ups(pp)

      if (dE / kT < 700.0_f64) then
        q_exc = (coll_fac / (g(1) * sqrt_Te)) * Ups(pp) * exp(-dE / kT)
      else
        q_exc = 0.0_f64
      end if

      avsum = 0.0d0  

      !print*,'1-->',i,q_exc,q_deexc,q_exc/q_deexc

      do j = 1,i-1 
        pp = upperTriangleIndexing(j,i,nlev)
        avsum = avsum + aval(pp)
      end do 
      avsum = avsum + Ne * q_deexc
      coronalPop(i) = Ne * q_exc / avsum
    end do 

    avsum = sum(coronalPop(:))
    coronalPop(:) = coronalPop(:) / avsum 

  end subroutine

  subroutine BoltzmanPopulation(nlev, g, E, Te, boltzpop)
    implicit none
    integer,  intent(in)  :: nlev
    real(f64), intent(in)  :: g(nlev)
    real(f64), intent(in)  :: E(nlev)

    real(f64), intent(in)  :: Te
    real(f64)              :: boltzpop(nlev)

    integer  :: i
    real(f64) :: kT, sqrt_Te, dE,avsum,w1,wi


    boltzpop(1) = 1.0d0 
    kT      = kB_eV * Te
    sqrt_Te = sqrt(Te)
    w1 = g(1)

    do i = 2, nlev 
      wi = g(i)
      dE  = E(i) - E(1)

      boltzpop(i) = (wi  / w1) * exp( -de/kt) 
      !print*,boltzpop(i)
    end do 

    avsum = sum(boltzpop(:))
    !print*,avsum
    boltzpop(:) = boltzpop(:) / avsum 

  end subroutine


  subroutine build_cr_matrix(nlev, ntran, g, E, Ups, Aval, sob, Te, Ne, Q, Qcol1, ierr,writeoutrates)
    !This routine was AI assisted. 
    !Has been tested against ColRadPy. 
    implicit none
    integer,  intent(in)  :: nlev, ntran
    real(f64), intent(in)  :: g(nlev)
    real(f64), intent(in)  :: E(nlev)
    real(f64), intent(in)  :: Ups(ntran)
    real(f64), intent(in)  :: Aval(ntran)
    real(f64), intent(in)  :: sob(ntran)

    real(f64), intent(in)  :: Ne
    real(f64), intent(in)  :: Te
    real(f64), intent(out) :: Q(nlev-1, nlev-1)
    real(f64)              :: Qcol1(nlev)
    real(f64)              :: Qrow1(nlev-1)
    integer,  intent(out) :: ierr

    integer  :: i, j, ii, jj
    real(f64) :: kT, sqrt_Te, dE, q_exc, q_deexc
    logical :: writeoutrates

    ierr    = 0
    kT      = kB_eV * Te
    sqrt_Te = sqrt(Te)
    Q       = 0.0_f64
    Qcol1   = 0.0_f64 
    Qrow1   = 0.0_f64

    !print*,'build_cr_matrix: '
    !print*,'  nlev=',nlev
    !print*,'  ntran=',ntran
    !print*,'  Te=',Te
    !print*,'  Ne=',Ne

    if (writeoutrates) then 
      open(900,file='atomicRates.dat')
      write(900,*)'# Low   Upp   Upsilon      qexc    qdeexc      Aval      Besc Besc*Aval'

    end if 


    do j = 2, nlev
      jj = j - 1
      call cr_matrix_element(1, j, Qrow1(jj), Qcol1(jj)) 
      do i = 2, j - 1
        ii = i - 1
        call cr_matrix_element(i, j, Q(ii,jj), Q(jj,ii))
      end do
    end do

    !do j = 2, nlev
    !  jj = j - 1
    !  call cr_matrix_element(1, j, Qrow1(jj), Qcol1(jj))
    !end do

    !This part takes the off diagonal elements and constructs the diagonal elements.
    do j = 2, nlev
      jj = j - 1
      Q(jj, jj) = 0.0_f64
      do i = 2, nlev
        ii = i - 1
        if (i /= j) Q(jj, jj) = Q(jj, jj) - Q(ii, jj)
      end do
      Q(jj, jj) = Q(jj, jj) - Qrow1(jj)
    end do


    if (writeoutrates) then 
      close(900)
    end if 

    !open(88,file='crmmatrix')
    !do j = 2, nlev
    !  jj = j - 1
    !  write(88,'(10000(ES10.3,1X))') Q(:,jj)
!
    !end do
    !close(88)

    !write(*,'(A)') 'Column sum check (should be ~ 0):'
    !do j = 2, nlev
    !  jj = j - 1
    !  write(*,'(A,I4,A,ES12.4)') '  col ', j, ': ', sum(Q(:,jj)) + Qrow1(jj)
    !end do

  contains

    subroutine cr_matrix_element(iii, jjj, cij, cji)
      implicit none
      integer,  intent(in)    :: iii, jjj
      real(f64), intent(inout) :: cij, cji
      integer :: pp

      pp  = upperTriangleIndexing(iii, jjj, nlev)
      dE  = E(jjj) - E(iii)

      q_deexc = (coll_fac / (g(jjj) * sqrt_Te)) * Ups(pp)

      if (dE / kT < 700.0_f64) then
        q_exc = (coll_fac / (g(iii) * sqrt_Te)) * Ups(pp) * exp(-dE / kT)
      else
        q_exc = 0.0_f64
      end if

      cji = cji + Ne * q_exc
      cij = cij + Ne * q_deexc + Aval(pp) * sob(pp)

      if (writeoutrates) then 
        write(900,'(2I6,6ES10.3)') iii,jjj, Ups(pp), q_exc, q_deexc, Aval(pp), sob(pp), Aval(pp) * sob(pp)              
      end if 


    end subroutine cr_matrix_element

  end subroutine build_cr_matrix

subroutine solve_cr_populations_axb(nlev, Q, numlevelsincluded, Qcol1, ierr, use_expert)
    !This routine was AI assisted - particularly for the use of dgesvx in cases where the 
    !user is worried about accuracy.

    implicit none
    integer,   intent(in)    :: nlev, numlevelsincluded
    real(f64), intent(inout) :: Q(nlev-1, nlev-1)
    real(f64), intent(inout) :: Qcol1(nlev)
    integer,   intent(out)   :: ierr
    logical,   intent(in)    :: use_expert

    integer, allocatable :: ipiv(:)
    integer :: i, info, n, ninc

    ! dgesvx-only variables
    real(f64), allocatable :: AF(:,:), R(:), C(:), B(:), X(:), ferr(:), berr(:)
    real(f64) :: rcond
    character(1) :: equed
    integer, allocatable :: ipiv_ex(:)
    real(f64), allocatable :: work(:)
    integer,   allocatable :: iwork(:)

    n    = nlev - 1
    ninc = numlevelsincluded - 1
    ierr = 0

    do i = 1, numlevelsincluded - 1
        Qcol1(i) = -1.0_f64 * Qcol1(i)
    end do
    do i = numlevelsincluded, nlev
        Qcol1(i) = 0.0d0
    end do

    if (use_expert) then
        !begin gemini contrbution:
        !write(0,*) 'using expensive la'

        allocate(AF(n, n), R(ninc), C(ninc), B(ninc), X(ninc), ferr(1), berr(1), ipiv_ex(ninc))
        allocate(work(4*n), iwork(n))

        !write(0,*) 'ninc=', ninc, ' n=', n, ' shape(Q)=', shape(Q), ' shape(AF)=', shape(AF)
        !write(0,*) 'shape(R)=', shape(R), ' shape(C)=', shape(C)
        !write(0,*) 'shape(B)=', shape(B), ' shape(X)=', shape(X)

        B(1:ninc) = Qcol1(1:ninc)


        ! dgesvx: equilibrates, factors, solves, and gives forward/backward error bounds
        ! and a reciprocal condition number estimate (rcond).
        ! 'E' = equilibrate, factor, solve.  'N' = no transpose.
        call dgesvx('E', 'N',           &
                    ninc, 1,            &   ! matrix size, nrhs
                    Q, n,               &   ! A, lda
                    AF, n,              &   ! factored A (output), ldaf
                    ipiv_ex,            &   ! pivot indices
                    equed,              &   ! equilibration actually applied (output)
                    R, C,               &   ! row/col scale factors (output)
                    B, ninc,            &   ! RHS, ldb
                    X, ninc,            &   ! solution (output), ldx
                    rcond,              &   ! reciprocal condition number (output)
                    ferr, berr,         &   ! forward/backward error bounds (output)
                          work, iwork, &

                    info)

        if (rcond < epsilon(1.0_f64)) then
            write(6,'(A,ES10.3)') 'WARNING: matrix is near-singular, rcond = ', rcond
        end if

        write(6,'(A,ES10.3,A,ES10.3)') &
            'dgesvx: rcond = ', rcond, '  ferr = ', ferr(1)

        Qcol1(1:ninc) = X(1:ninc)

        deallocate(AF, R, C, B, X, ferr, berr, ipiv_ex,work,iwork)
        !end gemini contribution.

    else

        allocate(ipiv(n))
        call dgesv(ninc, 1, Q, n, ipiv, Qcol1, nlev, info)
        deallocate(ipiv)

    end if

    if (info /= 0) then
        write(*,'(A,I4)') 'ERROR: solver failed, info = ', info
        ierr = info
        return
    end if

    do i = ninc, 1, -1
        Qcol1(i+1) = Qcol1(i)
    end do
    Qcol1(1) = 1.0_f64
    Qcol1    = Qcol1 / sum(Qcol1)

    if (any(Qcol1 < 0.0_f64)) then
        stop ' negative pops - numerical stability '
    end if

end subroutine solve_cr_populations_axb

!  subroutine solve_cr_populations_axb(nlev, Q,numlevelsincluded, Qcol1, ierr)
!    implicit none
!    integer,  intent(in)    :: nlev,numlevelsincluded
!    real(f64), intent(inout)    :: Q(nlev-1, nlev-1)
!    real(f64), intent(inout) :: Qcol1(nlev)
!    integer,  intent(out)   :: ierr
!
!    integer, allocatable :: ipiv(:)
!    integer :: i, info, n,ninc
!
!    n    = nlev - 1
!    ninc = numlevelsincluded - 1 
!    ierr = 0
!
!    allocate(ipiv(n))
!
!    !Indexing is off by one, as Qcol1 does not yet include the ground.
!    do i = 1, numlevelsincluded-1
!      Qcol1(i) = -1.0_f64 * Qcol1(i)
!    end do
!
!    do i = numlevelsincluded, nlev 
!      Qcol1(i) = 0.0d0 
!    end do 
!    
!
!
!    !If I wanted to do a study on including a different number of levels,
!    !I would need to change this to something like:
!    ! call dgesv(n_being_used, 1, Q, n, ipiv, Qcol1, nlev, info)
!    !write(*,*) 'Calling dgesv with ', numlevelsincluded, ' levels included.'
!
!    call dgesv(ninc, 1, Q, n, ipiv, Qcol1, nlev, info)
!    
!    !some high lying levels will have unphysical populations. The correct
!    !answer is so small that it doesnt matter. 
!
!    deallocate(ipiv)  !One must consider if we should be pre-allocating ipiv and reusing it.
!                      ! It gets resued at every sobolev iteration.
!
!    if (info /= 0) then
!      write(*,'(A,I4)') 'ERROR: DGESV failed, info = ', info
!      ierr = info
!      return
!    end if
!
!    do i = ninc, 1, -1
!      Qcol1(i+1) = Qcol1(i)
!    end do
!    Qcol1(1) = 1.0_f64
!    !write(*,*) Qcol1(:)
!    Qcol1 = Qcol1 / sum(Qcol1)
!
!    !do i = 1, nlev
!    !  write(*,*) i, Qcol1(i)
!    !end do
!
!    !write(*,'(A,I6)') 'Negative # = ', count(Qcol1 < 0.0_f64)
!!
!    if (any(Qcol1 < 0.0_f64)) then
!      stop ' negative pops - numerical stability '
!    end if
!
!  end subroutine solve_cr_populations_axb

  subroutine calculate_total_radiative_cascade(nlev,ntran,avals,cascade)
    integer   :: nlev, ntran 
    real(f64) :: avals(ntran), cascade(nlev) 

    integer :: ii ,jj, pp  

    cascade(:) = 0.0_f64

    do ii = 2, nlev 
      do jj = 1,ii-1 
        !print*,jj,ii
        pp  = upperTriangleIndexing(jj, ii, nlev)
        cascade(ii) = cascade(ii) + avals(pp)
      end do 
    end do 

  end subroutine

    subroutine  getAtomicDensityLocal(denslocal,&
                                    numIonsLocal,&
                                   ionMassSolar,& 
                                   atomicnumber,& 
                                   velocity_outer,& 
                                   velocity_inner,& 
                                   fractionOverride,& 
                                   time_exp_days,& 
                                   electron_density_local) 
    implicit none 
    real(f64) :: denslocal, ionMassSolar,velocity_outer,velocity_inner,fractionOverride,time_exp_days
    real(f64)  :: expansion_volume,time_exp_sec,electron_density_local,numIonsLocal
    integer :: atomicnumber
    real(f64),parameter :: c_cgs = 3e10_f64

    time_exp_sec = time_exp_days * 86400.0_f64

    !total volume.
    expansion_volume = piFourOnThree * (velocity_outer*c_cgs * time_exp_sec) ** 3
    expansion_volume = expansion_volume - piFourOnThree * (velocity_inner*c_cgs * time_exp_sec) ** 3
    numIonsLocal = ionMassSolar * m_solar_grams/get_mass_grams(atomicnumber)

    denslocal =   numIonsLocal / (expansion_volume)

    if (fractionOverride > 0.0_f64) denslocal = fractionOverride * electron_density_local
    print*,'atomic number density', denslocal,'cm-3 from new routine. edense=', electron_density_local
    write(0,*) denslocal
  end subroutine 
    
  subroutine sobolev_escape(nlev,ntran,baseAvals,sobesc,time_exp_days,pops,weights,wl_cm_cubed,atomicDensityLocal)
    !calculates Sobolev escape probability. 
    implicit none 
    integer :: nlev , ntran
    !
    real(f64) ::  pops(nlev), weights(nlev)
    real(f64) :: baseAvals(ntran),wl_cm_cubed(ntran)
    real(f64) :: sobesc(ntran)
    ! 
    real(f64) :: atomicDensityLocal
    real(f64) :: tau ,time_exp_days,time_exp_sec
    real(f64),parameter :: c_cgs = 3e10_f64
    integer :: pp 
    integer :: ii , jj 
    !
    time_exp_sec = time_exp_days * 86400.0_f64


    
    print*,'number density', atomicDensityLocal,'cm-3'


    sobesc(:) = 1.0_f64

    !write(0,*) 'sobconst',sobconst,atomicDensityLocal

    do ii = 1, nlev-1 
       do jj = ii+1,nlev
          pp  = upperTriangleIndexing(ii,jj,nlev)

          tau = sobconst * baseAvals(pp) * wl_cm_cubed(pp) * weights(jj) * atomicDensityLocal * time_exp_sec * (pops(ii)/weights(ii) -  pops(jj)/weights(jj))
          !print*,tau, sobconst,baseAvals(pp),atomicDensityLocal,time_exp_sec

          if (tau > 1.0e-5_f64) sobesc(pp) = (1.0_f64 - exp(-tau)) / tau  
          if ( (tau < 0.0_f64) ) write(999,*) ii, pops(ii), jj,pops(jj),baseAvals(pp), tau
        end do 
    end do 
    !
  end subroutine

  subroutine broadenedSpectrum(numWavelengths,& 
                                       wavelength,&
                                       velocityShell,&
                                       spectra, &
                                       ntran, &
                                       pec, &
                                       spectralLinesCM,&
                                       electron_density,&
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
    totalpec = sum(pec)
    peccutoff = pecthreshold * totalpec

    ! Calculate spectrum based on the selected mode
    do ii = 1, ntran 
      wavelengthCentral = spectralLinesCM(ii)

      if (wavelengthCentral < 1e-30_f64) cycle

      thispec = pec(ii)
      if (thispec < peccutoff) cycle

      thisphotonenergy = thispec 
      !write(0,*) 'broad mode ', calcmode,trim(adjustl(calcMode)),trim(adjustl(calcMode))=='box'

      ! --- BRANCH: BOX PROFILE (Expanding Shell) ---
      if ( (calcMode(1:3) == 'box') .or. trim(adjustl(calcMode)) == 'onion') then
        write(0,*) 'i am doing a box'
        ! For an expanding shell, max Doppler shift is defined by the shell velocity.
        ! Assuming velocityShell here represents v/c (or v_expansion / c).
        wl_lo = wavelengthCentral * (1.0_f64 - velocityShell)
        wl_hi = wavelengthCentral * (1.0_f64 + velocityShell)

        ! Find grid bins
        jlo = max(1,              nint((wl_lo - wavelength(1)) / dwl) + 1)
        jhi = min(numWavelengths, nint((wl_hi - wavelength(1)) / dwl) + 1)

        ! The box height is 1 / Total Width.
        ! Extra (1 / wavelengthCentral) applied for E = hc/λ step.
        normfactor = (1.0_f64 / (wl_hi - wl_lo)) / wavelengthCentral
        ww = thisphotonenergy * normfactor 

        ! Flat profile: Add uniform intensity to all bins within the box
        do jj = jlo, jhi
          spectra(jj) = spectra(jj) + ww 
        end do

      ! --- BRANCH: GAUSSIAN PROFILE (Thermal/Microturbulence) ---
      else 

        sig_cm = fwhmSigma * wavelengthCentral * velocityShell
        sigOneOver = 1.0_f64 / sig_cm

        ! Gaussian normalization.
        ! Extra (1 / wavelengthCentral) applied for E = hc/λ step.
        normfactor = sigOneOver * oneOverSQRTTWOPI / wavelengthCentral

        wl_lo = wavelengthCentral - nSigmaCut * sig_cm
        wl_hi = wavelengthCentral + nSigmaCut * sig_cm
        
        jlo = max(1,              nint((wl_lo - wavelength(1)) / dwl) + 1)
        jhi = min(numWavelengths, nint((wl_hi - wavelength(1)) / dwl) + 1)

        do jj = jlo, jhi
          ww = wavelength(jj) 
          ww = (ww - wavelengthCentral) * sigOneOver
          ww = thisphotonenergy * normfactor * exp ( minusHalf * ww * ww )
          spectra(jj) = spectra(jj) + ww 
        end do 

      end if

    end do 
    
    ! Final scaling: Units conversion (Photons -> Ergs)
    spectra(:) = spectra(:) * ( numIonsLocal * electron_density * hc_ergcm * 1e-8_f64 ) 
  
  end subroutine

!  subroutine broadenedSpectrum(numWavelengths,& 
!                                       wavelength,&
!                                       velocityShell,&
!                                       spectra, &
!                                       ntran, &
!                                       pec, &
!                                       spectralLinesCM,&
!                                       electron_density,&
!                                       numIonsLocal, & 
!                                       calcMode &
!                                       )
!    !
!    ! Calculates \sum_i nf * PEC * exp ( -0.5 *  [ (λ - λ_0 ) / σ ]^2)
!    ! where all wavelengths are in nm. the norm factor nf - is 1/nf = 10 sqrt(2pi) σ
!    ! where the extra factor of 10 puts it in units of per angstrom.
!    ! units are annoying. Should probably just keep everything in cgs and ship 
!    ! a post processor with astropy, for my own sanity. 
!    !
!    character* 10 :: calcmode
!    integer :: numWavelengths
!    integer :: ntran 
!    real(f64) :: wavelength(numWavelengths),electron_density
!    real(f64) :: spectra(numWavelengths)
!    real(f64) :: pec(ntran)
!    real(f64) :: spectralLinesCM(ntran)
!    real(f64) :: ww, wavelengthCentral
!    real(f64) :: velocityShell, sig_cm, sigOneOver,numIonsLocal
!    integer :: ii, jj 
!
!
!    real(f64) :: normfactor
!    real(f64) :: thispec,wl_lo,wl_hi, nsigmacut=4,dwl
!    integer :: jlo, jhi 
!    real(f64) :: totalpec 
!    real(f64) :: pecthreshold = 1e-4
!    real(f64) :: peccutoff
!    real(f64) :: thisphotonenergy 
!    dwl = wavelength(2) - wavelength(1)
!    totalpec = sum(pec)
!
!    peccutoff = pecthreshold * totalpec
!
!    !velocityShell = velocityFWHM in the case of a Gaussian model, e.g for a single shell.
!
!    ! Calculate Gaussian spectrum
!    do ii = 1, ntran 
!      !cm
!      wavelengthCentral = spectralLinesCM(ii)
!
!      if (wavelengthCentral < 1e-30) cycle
!
!      thispec = pec(ii)
!
!      if (thispec < peccutoff) cycle
!
!      sig_cm = fwhmSigma * wavelengthCentral  * velocityShell
!
!      sigOneOver = 1._f64/sig_cm
!
!      normfactor = sigOneOver * oneOverSQRTTWOPI / wavelengthCentral
!
!      wl_lo = wavelengthCentral - nSigmaCut * sig_cm
!      wl_hi = wavelengthCentral + nSigmaCut * sig_cm
!      jlo = max(1,              nint((wl_lo - wavelength(1)) / dwl) + 1)
!      jhi = min(numWavelengths, nint((wl_hi - wavelength(1)) / dwl) + 1)
!      thisphotonenergy = thispec 
!
!      do jj = jlo, jhi
!
!        ww = wavelength(jj) 
!        ww = (ww - wavelengthCentral) 
!        ww = ww * sigOneOver
!        ww = thisphotonenergy * normfactor * exp ( minusHalf * ww * ww  )
!        spectra(jj) = spectra(jj) + ww 
!
!      end do 
!    end do 
!    
!    spectra(:) = spectra(:) * ( numIonsLocal * electron_density * hc_ergcm  * 1e-8) !1e-8 to get in ergs per ang per s
!    
!    !print*,num_ions, maxval(spectra)
!
!  
!  end subroutine


end module crm_module