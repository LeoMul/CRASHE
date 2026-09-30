module atomicdata_module
   !
   ! Atomic data from an adf04 file. Filled by readadf04 / readhack
   ! (readadf04_module); read by everything else that needs the atomic data.
   !
   ! Transition-indexed arrays (aval, wl_cm, wl_cm_cubed, ups(:, pp)) are indexed
   ! with upperTriangleIndexing(lower, upper, numLevels); the readers allocate
   ! numLevels*(numLevels+1)/2 entries, of which the first ntran are used.
   !
   use types
   implicit none

   integer :: numLevels = 0          ! number of levels
   integer :: numTemps = 0           ! number of temperatures upsilons are tabulated at
   integer :: ntran = 0              ! number of transitions, numLevels*(numLevels-1)/2
   integer :: atomicNumber = 0       ! nuclear charge
   integer :: ioncharge_plus = 0     ! charge of the ion

   real(f64), allocatable :: energies(:)      ! level energies (eV)
   real(f64), allocatable :: statweight(:)    ! statistical weights 2J+1
   real(f64), allocatable :: temps(:)         ! temperatures of the tabulated upsilons (K)
   real(f64), allocatable :: ups(:, :)        ! effective collision strengths, (numTemps, transition)
   real(f64), allocatable :: aval(:)          ! Einstein A coefficients (s-1), per transition
   real(f64), allocatable :: wl_cm(:)         ! transition wavelengths (cm), per transition
   real(f64), allocatable :: wl_cm_cubed(:)   ! wl_cm**3, per transition

   real(f64), allocatable :: qup(:), qdown(:)


contains

   subroutine dealloc_atomicdata
      ! Deallocates the atomic data and resets the counters.
      implicit none
      if (allocated(energies)) deallocate (energies)
      if (allocated(statweight)) deallocate (statweight)
      if (allocated(temps)) deallocate (temps)
      if (allocated(ups)) deallocate (ups)
      if (allocated(aval)) deallocate (aval)
      if (allocated(wl_cm)) deallocate (wl_cm)
      if (allocated(wl_cm_cubed)) deallocate (wl_cm_cubed)
      numLevels = 0
      numTemps = 0
      ntran = 0
   end subroutine dealloc_atomicdata

end module atomicdata_module