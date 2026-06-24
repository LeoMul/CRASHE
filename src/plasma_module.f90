module plasma_module 
    !This module just calculates the PLT. I should probably just merge it with crm_module.
    use, intrinsic :: iso_fortran_env, only: f64 => real64
    use readadf04_module, only: upperTriangleIndexing
    implicit none
    real(f64) :: eV_to_erg = 1.602176634e-12_f64
    contains 
    subroutine calculate_pec_plt(nlev,pops,ntran,avals,sob,pecs,plt,density,energy)
        integer,  intent(in)  :: nlev , ntran 
        real(f64),intent(in)  :: pops(nlev) 
        real(f64),intent(in)  :: energy(nlev) 
        real(f64),intent(in)  :: avals(ntran)
        real(f64),intent(in)  :: sob(ntran)

        real(f64),intent(out) :: pecs(ntran)
        real(f64),intent(out) :: plt 
        real(f64)             :: density
        real(f64)             :: ei, ej 
        integer               :: ii , jj, pp
        !
        plt = 0.0_f64
        do ii = 1, nlev-1 
            ei = energy(ii)
            do jj = ii+1,nlev
                ej = energy(jj)
                pp = upperTriangleIndexing(ii,jj,nlev)
                pecs(pp) = pops(jj) * avals(pp) * sob(pp) / density
                plt  = plt + pecs(pp) * (ej - ei)
                !print*,pecs(pp),pops(jj),density
            end do 
        end do 
        !
        plt = plt * eV_to_erg 
        !
    end subroutine calculate_pec_plt 
!
end module plasma_module 
