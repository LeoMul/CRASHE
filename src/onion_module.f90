module onion_module
    use colradfort 
    use types 
    use input
    implicit none 
    contains 
    subroutine onion 
        implicit none 
        ! dummy variables for testing. 
        integer,parameter   :: numshells = 10 
        real(f64) :: averageElectronDensity  
        real(f64) :: totalMassMsun = 1e-3 
        real(f64) :: v_min, v_max, fraction
        real(f64) :: v_bounds_c(numshells+1)
        integer   :: s 
        real(f64)   :: shell_v_centers(numShells) ! Shell center velocities (for CRM)
        real(f64)   :: shell_electron_density(numShells) 
        real(f64)   :: shell_mass(numShells) 
        real(f64)   :: total_rel_vol, shell_rel_vol, total_rel_electrons
        real(f64)   :: current_avg_ne, scale_ne, total_rel_mass, scale_mass
        
        reaL(f64) :: dwl 
        real(f64) :: wl1  
        real(f64) :: wl2 
        integer   :: nwl  
        real(f64) :: atomicDensityLocal,numions
        real(f64),allocatable :: currentspec(:)
        real(f64),allocatable :: totalspec(:)
        real(f64),allocatable :: wlarray(:)
        character(len=300) :: broadmode
        !velocity law:
        v_min = 0.04_f64
        v_max = 0.30_f64

        !use user imposed electron density
        averageElectronDensity = density
        wl1 = wlmin_nm
        wl2 = wlmax_nm
        nwl = numwl 
        allocate(currentspec(nwl))
        allocate(totalspec(nwl))
        allocate(wlarray(nwl))

        ! Logarithmic spacing ensures high resolution where density is high
        do s = 1, numShells + 1
            fraction = real(s - 1, f64) / real(numShells, f64)
            v_bounds_c(s) = v_min * (v_max / v_min)**fraction
        end do

        do s = 1, numShells
            shell_v_centers(s) = 0.5_f64 * (v_bounds_c(s) + v_bounds_c(s+1))
        end do

        ! 2. Compute Relative Profiles using the v^-3 Power Law at shell centers
        do s = 1, numShells
           shell_electron_density(s) = shell_v_centers(s)**(-3.0_f64)
           shell_mass(s)    = shell_v_centers(s)**(-3.0_f64)
        end do

        ! 3. Normalize Electron Density Profile to match Target Bulk Density
        total_rel_electrons = 0.0_f64
        do s = 1, numShells
           shell_rel_vol = v_bounds_c(s+1)**3 - v_bounds_c(s)**3
           total_rel_electrons = total_rel_electrons + shell_electron_density(s) * shell_rel_vol
        end do

        total_rel_vol = v_bounds_c(numShells+1)**3 - v_bounds_c(1)**3
        current_avg_ne = total_rel_electrons / total_rel_vol
        scale_ne = averageElectronDensity / current_avg_ne

        shell_electron_density(:) = shell_electron_density(:) * scale_ne

        ! 4. Normalize Mass Profile cleanly
        ! Find the sum of the raw v^-3 bin values
        total_rel_mass = sum(shell_mass)
        
        ! Find the scale factor needed to make the array sum up to totalMassMsun
        scale_mass = totalMassMsun / total_rel_mass
        
        ! Apply the scale factor directly to preserve the v^-3 shape per bin
        shell_mass(:) = shell_mass(:) * scale_mass

        total_rel_mass = 0.0_f64
        do s = 1, numShells
            shell_rel_vol = v_bounds_c(s+1)**3 - v_bounds_c(s)**3   ! proportional to physical vol
            total_rel_mass = total_rel_mass + shell_mass(s) * shell_rel_vol
        end do
        scale_mass = totalMassMsun / total_rel_mass
        
        ! Now convert from density shape to actual mass per shell
        do s = 1, numShells
            shell_rel_vol = v_bounds_c(s+1)**3 - v_bounds_c(s)**3
            shell_mass(s) = shell_mass(s) * scale_mass * shell_rel_vol
        end do


        ! 5. Output to file
        open(30, file = 'velocityProfile')
        do s  = 1, numshells
            write(30,'(4ES14.6)') shell_v_centers(s), shell_mass(s), shell_electron_density(s), v_bounds_c(s)
        end do
            write(30,'(4ES14.6)') 0.0,0.0,0.0, v_bounds_c(numshells+1)
        close(30)

        !Assume a uniform temperature now, for easiness.
        !Then, something like , 

        dwl = (wl2 - wl1) / nwl 
        wlarray(1) = wl1 
        do s = 2,nwl 
            wlarray(s) = wlarray(s-1) + dwl 
        end do 

        broadmode(1:3) = 'box'

        do s = 1, numshells
            
            call getAtomicDensityLocal(atomicDensityLocal, numions, shell_mass(s), atomicNumber,v_bounds_c(s+1),v_bounds_c(s),0.0_f64,29.0_f64,shell_electron_density(s) )


            call colrad(3000.0_f64, & 
            shell_electron_density(s),& 
            sobolev, & 
            29.0_f64,&
            atomicDensityLocal,&
            wl1,&
            wl2,&
            nwl,&
            .false.,&
            .false.,&
            shell_v_centers(s),& 
            wlarray,&
            currentspec,&
            numions, &
            broadmode & 
            )
            write(0,*) broadmode(1:3)
            totalspec(:) = totalspec(:) + currentspec
            !call colrad  !with a flag to do the box instead. 
            !add spectra to  full spectra?
            !output?

        end do 
       open(101,file='spectrum')
       do j = 1, nwl
           if (totalspec(j) >0) then 
           if (totalspec(j) < 1e40_f64) write(101,*) wlarray(j)*1e7, totalspec(j)
           end if 
       end do
       close(101)
       deallocate(currentspec)
       deallocate(totalspec)
       deallocate(wlarray)
    end subroutine


end module 