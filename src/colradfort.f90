module colradfort 
!module to pretty much call everything.
    use types
    use readadf04_module
    use crm_module
    use interpolation_module
    use plasma_module
    use omp_lib
    use input, only: mode,contourLower,contourUpper,sortpec
    use sorting
    implicit none
    integer                :: thrid
    integer                :: numlevels, numtemps, ntran, ierr, i, numTempsReq
    real(f64)              :: temp
    real(f64)              :: plt, pltnosob
    real(f64), allocatable :: tempsReq(:)
    real(f64), allocatable :: temps(:)
    real(f64), allocatable :: ups(:,:)
    real(f64), allocatable :: upsInterp(:)
    real(f64), allocatable :: aval(:), cascade(:)
    real(f64), allocatable :: statweight(:)
    real(f64), allocatable :: energies(:)
    real(f64), allocatable :: pec(:) ,crm(:,:), col1(:), pops_old(:)
    real(f64), allocatable :: popcoronal(:)
    real(f64), allocatable :: popsnosob(:), pecnosob(:)
    real(f64), allocatable :: wl_cm(:), wl_cm_cubed(:)
    real(f64), allocatable :: sob(:), sob_old(:)
    !
    real(f64), allocatable :: wavelengthforspectrum(:)
    real(f64), allocatable :: broadspec(:)
    !
    real(f64)              :: atomicDensity,numions
    real(f64)              :: sob_damp    = 0.5_f64
    real(f64), parameter   :: sob_tol     = 1.0e-2_f64
    integer,   parameter   :: max_sob_iter = 9999
    integer                :: sob_iter
    real(f64)              :: sob_change,beta_change,beta_change_old=1.d6
    logical                :: converged
    integer                :: k, j, p
    integer                :: atomicNumber
    integer                :: ioncharge_plus
    real(8) :: t1, t2 
    character(len=300)     ::broadmodedefault = 'gaussian'
    integer*8              :: shellnumtemp=0

    contains 

    subroutine getadf04(adf04Path,floersHack)
        implicit none 
        logical :: floershack
       character(len=*) :: adf04Path

       call cpu_time(t1)

       if (floersHack) then 

       call readhack(trim(adf04Path), numLevels, numTemps, ups, aval, &
                      statweight, energies, temps, wl_cm, wl_cm_cubed,atomicNumber,ioncharge_plus)
        
       else
       call readadf04(trim(adf04Path), numLevels, numTemps, ups, aval, &
                      statweight, energies, temps, wl_cm, wl_cm_cubed,atomicNumber,ioncharge_plus)
       end if 

       call cpu_time(t2)

       write(*,'(A,ES10.4,A)') '  [timing] adf04 read time     : ', t2-t1, ' s'

       ntran = (numLevels * (numLevels - 1)) / 2
    end subroutine

    subroutine colrad(temperature,            & 
                      electronDensityLocal,   & 
                      sobolev,                &
                      timeSinceExplosionDays, &
                      atomicDensityLocal,     &
                      wlmin_nm,               &
                      wlmax_nm,               &
                      numwl,                  &
                      careful_la ,            &
                      writeoutrates,          &
                      velocityExpansionC,     &
                      wlspec,                 &
                      bspec,                  & 
                      numIonsLocal,           &
                      broadmode               &
                      )
       
       real(f64) :: temperature 
       real(f64) :: electronDensityLocal
       real(f64) :: velocityExpansionC     
       real(f64) :: timeSinceExplosionDays 
       real(f64) :: massElementSolar       
       real(f64) :: wlmin_nm, wlmax_nm, dwl
       integer   :: numwl
       real(f64) :: wlspec(numwl)
       real(f64) :: bspec(numwl)
       character(len=20) :: filesuffix
       real(f64) :: atomicDensityLocal,numIonsLocal
    
       logical :: sobolev,careful_la,writeoutrates
       character(len=300) :: broadmode
       integer, allocatable :: pecPointer(:)

       tempsReq(1) = temperature 
        shellnumtemp = shellnumtemp + 1
        i =1
        !write(0,*) velocityExpansionC,massElementSolar
        call interpolate_upsilons(ntran, numTemps, temps, &
                                  temperature, ups, upsInterp)
        sob      = 1.0_f64
        
        call cpu_time(t1)
        call build_cr_matrix(numLevels, ntran, statweight, energies, &
               upsInterp, aval, sob, tempsReq(i), electronDensityLocal, crm, col1, ierr,writeoutrates)
        call solve_cr_populations_axb(numLevels, crm,numLevels, col1, ierr,careful_la)

        !write(0,*) 'col1' , col1(:)

        call BoltzmanPopulation(numlevels,statweight,energies,tempsReq(i),popcoronal)
        
        call cpu_time(t2)
        write(*,'(A,ES10.4,A)') '  [timing] initial populations : ', t2-t1, ' s'
        converged = .false.
        sob      = 1.0_f64
        sob_old  = 1.0_f64
        pops_old = 0.0_f64
        
        popsnosob = col1
        call calculate_pec_plt(numLevels, col1, ntran, aval, sob, pec, plt, electronDensityLocal, energies)
        
        bspec(:) = 0.0d0 
        wlspec(1)     = wlmin_nm * 1e-7 
        wlspec(numwl) = wlmax_nm * 1e-7

        dwl = 1e-7 * (wlmax_nm - wlmin_nm) / (numwl-1)
        do j = 2, numwl-1
            wlspec(j) = wlspec(j-1) + dwl
        end do 
        

        if (sobolev) then 
            popsnosob = col1 
            pecnosob  = pec 
            pltnosob  = plt 
            call sobolev_escape(numLevels, ntran, aval, sob, timeSinceExplosionDays, col1, &
                                statweight, wl_cm_cubed,atomicDensityLocal)
            !write(0,*) 'sob' , sob(:)
            !write(0,*) 'col1' , col1(:)

            call cpu_time(t1)
            sob_iter_loop: do sob_iter = 1, max_sob_iter
                call build_cr_matrix(numLevels, ntran, statweight, energies, &
                                     upsInterp, aval, sob, tempsReq(i), electronDensityLocal, crm, col1, ierr,writeoutrates)
                call solve_cr_populations_axb(numLevels, crm,numLevels, col1, ierr,careful_la)
                !write(0,*) 'crm' , crm(:,:)
                !write(0,*) 'sob' , sob(:)
                sob_old = sob

                call sobolev_escape(numLevels, ntran, aval, sob, timeSinceExplosionDays, col1, &
                                    statweight, wl_cm_cubed,atomicDensityLocal)
                sob     = sob_damp * sob + (1.0_f64 - sob_damp) * sob_old

                !this is a fairly conservative convergence criterion - basically it asserts that 
                !none of the beta's change by more than 0.1%, for sob_tol = 1e-3.
                beta_change = maxval(abs(sob - sob_old)/sob)
                
                !call calculate_pec_plt(numLevels, col1, ntran, aval, sob, pec, plt, dens, energies)
                !print*,plt

                if (sob_iter > 1 .and. beta_change < sob_tol) then
                    converged = .true.
                    write(*,'(A,I4,A,ES10.3)') ' [sobolev] converged at iter   :', sob_iter
                    write(*,'(A,ES10.4)')      '       with maximum dBeta/Beta : ', beta_change                    
                    exit sob_iter_loop
                end if

                !if (beta_change < beta_change_old) then
                !    sob_damp = min(sob_damp * 1.1, 1.0_f64) ! Can go to 1.0
                !else
                !    sob_damp = max(sob_damp * 0.8, 0.1_f64) ! Don't crash to 0
                !end if
                beta_change_old = beta_change


                !print*,beta_change,beta_change_old,sob_damp

            end do sob_iter_loop

            call cpu_time(t2)
            write(*,'(A,ES10.4,A)') '  [timing] Sobolev iteration   : ', t2-t1, ' s'

            if (.not. converged) then
                write(*,'(A,I4,A,I3,A,2ES10.2)') &
                    'WARNING: Sobolev did not converge for temp index ', i, &
                    ' after ', max_sob_iter, ' iterations',beta_change,beta_change_old
            end if
        
        end if 


       call cpu_time(t1)
       call calculate_pec_plt(numLevels, col1, ntran, aval, sob, pec, plt, electronDensityLocal, energies)
       call cpu_time(t2)
       write(*,'(A,ES10.4,A)') '  [timing] PEC/PLT calculation : ', t2-t1, ' s'

       call calculate_total_radiative_cascade(numlevels,ntran,aval,cascade)

       if (shellnumtemp < 10) then 
        write(filesuffix,'(2I1)') 0,shellnumtemp
       else 
        write(filesuffix,'(I2)') shellnumtemp 
       end if 

       if (mode .eq. 'astro') filesuffix(:) =''

       open(100,file='popData'//trim(filesuffix))
        do j = 1,numlevels 
            write(100,'(I4,3ES11.4)') j ,col1(j),popcoronal(j),cascade(j) !, popcoronal(j)
        end do 
       close(100)

       open(100,file='pecData'//trim(filesuffix))
       write(100,*) 'Low, Upp,     Sob,    aval,     pec,         wlcm,    popL,    popU'


       if (sortpec) then 
        allocate(pecPointer(size(pec))) 
        do j=1, size(pec)
            pecPointer(j) = j 
        end do 

        call qsort(pec, size(pec), pecPointer)

        do j = size(pec),1, -1
            write(100,'(3ES14.7)') aval(pecPointer(j)), wl_cm(pecPointer(j)), pec(j)
        end do 

        deallocate(pecPointer)

       else

        do j = 1, numLevels-1
          do k = j+1, numLevels
             p = upperTriangleIndexing(j, k, numLevels)
             write(100,'(2I5,3ES9.2,ES14.7,2ES9.2)') j, k, sob(p), aval(p),pec(p), wl_cm(p), col1(j), col1(k)
         end do
        end do
       end if 

       close(100)

       call cpu_time(t1)
       call broadenedSpectrum(size(wlspec),wlspec,velocityExpansionC,bspec,ntran,pec,wl_cm,electronDensityLocal,numIonsLocal,broadmode)
       call cpu_time(t2)
       write(*,'(A,ES10.4,A)') '  [timing] spectrum broadening : ', t2-t1, ' s'

       call cpu_time(t1)
       open(101,file='spectrum'//trim(filesuffix))
       do j = 1, size(wlspec)
           write(101,*) wlspec(j), bspec(j)
       end do
       close(101)
       call cpu_time(t2)
       write(*,'(A,ES10.4,A)') '  [timing] spectrum write      : ', t2-t1, ' s'

       close(1)

    end subroutine
        
    subroutine levelscan(temperature, electronDensityLocal, careful_la,writeoutrates)
       real(f64) :: temperature 
       real(f64) :: electronDensityLocal,dens
       real(f64),allocatable :: crm_copy(:,:), col1_copy(:)
       logical :: careful_la,writeoutrates
       tempsReq(1) = temperature 
       dens        = electronDensityLocal
       allocate(crm_copy(numlevels-1,numlevels-1))
        sob      = 1.0_f64
       call interpolate_upsilons(ntran, numTemps, temps, &
                                temperature, ups, upsInterp)

       call build_cr_matrix(numLevels, ntran, statweight, energies, &
              upsInterp, aval, sob, tempsReq(1), dens, crm, col1, ierr,writeoutrates)
        
       crm_copy(:,:) = crm(:,:) 
       col1_copy  = col1 
       open(32,file='plt_level_convergence.dat')
       do i = 2 ,numlevels 

         call solve_cr_populations_axb(numLevels, crm, i, col1, ierr,careful_la)

         call calculate_pec_plt(numLevels, col1, ntran, aval, sob, pec, plt, dens, energies)
         write(32,*) i , plt , col1(2),col1(1)

         crm(:,:) = crm_copy(:,:)
         col1(:)  = col1_copy(:) 
       end do 
       !open(100,file='popData')
       ! do j = 1,numlevels 
       !     write(100,'(I4,2ES11.4)') j ,col1(j) !, popcoronal(j)
       ! end do 
       !close(100)

       deallocate(crm_copy)
    end subroutine


    subroutine masscontour (temperature, electronDensityLocal, requiredlumo,careful_la,writeoutrates,verbose)
        implicit none
       real(f64) :: temperature 
       real(f64) :: electronDensityLocal
       real(f64) :: requiredlumo  
       real(f64) :: electronDensityLocalvary(1000)
       real(f64) :: temperaturevary(1000)
       real(f64) :: thislumo_per_ion
       real(f64) :: num_req, mass_req 
       real(f64) :: num_in_one_solar_mass
       integer :: ii,jj ,counterii=0,counterjj=0
       real(f64) :: xx 
       logical :: careful_la,writeoutrates
       logical :: verbose

       num_in_one_solar_mass = 1.0 * m_solar_grams/ get_mass_grams(atomicnumber)

       sob=1
       

     

       open(90,file = 'contour.out')



       !get central estimate 
       call interpolate_upsilons(ntran, numTemps, temps,temperature, ups, upsInterp)
       call getmassestimate(temperature, electronDensityLocal, mass_req,careful_la,writeoutrates)

       write(90,'(A, ES14.6,A)') '# Central temp     = ', temperature,' Kelvin'
       write(90,'(A, ES14.6,A)') '# Central dens     = ', electronDensityLocal,' /cm3'
       write(90,'(A, ES14.6,A)') '# Central estimate = ', mass_req,' Msun'
       write(90,'(A, I3)'    )   '# Atomic  number   = ', atomicnumber
       write(90,'(A, I3,A)'    ) '# Atomic  charge   = ', ioncharge_plus,' +'


       write(90,'(A, I10)'    )   '# ntemp = ',size(temperaturevary)
       write(90,'(A, I10)'    )   '# ndens = ',size(electronDensityLocalvary)

       xx = log10 (temps(size(temps)) / temps(1)) 
       xx = 10**(xx/size(temperaturevary))
       temperaturevary(1) = temps(1)
       do ii =  2, size(temperaturevary)
            temperaturevary(ii) = temperaturevary(ii-1) * xx
       end do 
       temperaturevary(size(temperaturevary)) = temps(numtemps)
       !electronDensityLocal grid
       electronDensityLocalvary(1) = 3.0 
       electronDensityLocalvary(size(electronDensityLocalvary)) = 13.0 
       xx = (electronDensityLocalvary(size(electronDensityLocalvary)) - electronDensityLocalvary(1)) / size(electronDensityLocalvary)
       do ii = 2,size(electronDensityLocalvary)
        electronDensityLocalvary(ii) = electronDensityLocalvary(ii-1) + xx  
       end do 

       electronDensityLocalvary(:) = 10 ** electronDensityLocalvary(:)

       if (.not. verbose) then 
    
        write(90,'(A)') '# temp vary'
        
        !vary electronDensityLocal
        do ii = 1, size(temperaturevary) 
         call interpolate_upsilons(ntran, numTemps, temps,temperaturevary(ii), ups, upsInterp)
         call getmassestimate(temperaturevary(ii), electronDensityLocal, mass_req,careful_la,writeoutrates)
         write(90,'(2ES14.6)') mass_req , temperaturevary(ii)       
        end do 


        write(90,'(A)') '# dens vary'

        call interpolate_upsilons(ntran, numTemps, temps,temperature, ups, upsInterp)
        do ii = 1, size(electronDensityLocalvary) 
         call getmassestimate(temperature, electronDensityLocalvary(ii), mass_req,careful_la,writeoutrates)
         write(90,'(2ES14.6)') mass_req , electronDensityLocalvary(ii)
        end do 

       
       else 
        
        open(25,file='tempgrid')
        open(26,file='densgrid')
        do jj = 1, size(temperaturevary), 10 
            write(25,'(1ES14.6)') temperaturevary(jj)
        end do 

        do ii = 1, size(electronDensityLocalvary), 10 
            write(26,'(1ES14.6)') electronDensityLocalvary(ii)
        end do 

        close(26)
        close(25)


         do jj = 1, size(temperaturevary), 10 
            counterjj = counterjj + 1
            counterii = 0
           call interpolate_upsilons(ntran, numTemps, temps,temperaturevary(jj), ups, upsInterp)
           do ii = 1, size(electronDensityLocalvary), 10 
            counterii = counterii + 1 
            call getmassestimate(temperaturevary(jj), electronDensityLocalvary(ii), mass_req,careful_la,writeoutrates)
            write(90,'(2I10,1ES14.6)') counterii,counterjj, mass_req
            end do 
         end do 


        
       end if 


       close(90)

       contains 

       subroutine getmassestimate(reqtemp, reqdens, reqmass,careful_la,writeoutrates)
        use input, only: contourLower, contourUpper
        implicit none
        
        real(f64) :: reqtemp, reqdens, reqmass
        logical :: careful_la,writeoutrates
        integer :: pp

        call build_cr_matrix(numLevels, ntran, statweight, energies, &
                upsInterp, aval, sob, reqtemp, reqdens, crm, col1, ierr,writeoutrates)
        call solve_cr_populations_axb(numLevels, crm,numLevels, col1, ierr,careful_la)
        pp = upperTriangleIndexing(contourLower,contourUpper,numlevels)
        write(1414,*) '   The line in contour is: ',wl_cm(pp), aval(pp)
        thislumo_per_ion = col1(contourUpper) * aval(pp) * hc_ergcm / wl_cm(pp)
        num_req = requiredlumo / thislumo_per_ion
        reqmass = num_req/num_in_one_solar_mass
       end subroutine

    end subroutine masscontour

    subroutine alloc(numwl)
       implicit none 
       integer :: numwl
       numTempsReq = 1

       allocate(crm(numLevels-1, numLevels-1))
       allocate(col1(numLevels), pops_old(numLevels))
       allocate(tempsReq(numTempsReq))
       allocate(upsInterp(ntran))
       allocate(pec(ntran))
       allocate(sob(ntran))
       allocate(sob_old(ntran))
       allocate(pecnosob(ntran))
       allocate(popcoronal(numlevels))
       allocate(cascade(numlevels))
       allocate(wavelengthforspectrum(numwl),broadspec(numwl))

    end subroutine

    subroutine dealloc
       if (allocated(tempsReq))            deallocate(tempsReq)
       if (allocated(temps))               deallocate(temps)
       if (allocated(ups))                 deallocate(ups)
       if (allocated(upsInterp))           deallocate(upsInterp)
       if (allocated(aval))                deallocate(aval)
       if (allocated(statweight))          deallocate(statweight)
       if (allocated(energies))            deallocate(energies)
       if (allocated(crm))                 deallocate(crm)
       if (allocated(col1))                deallocate(col1)
       if (allocated(pops_old))            deallocate(pops_old)
       if (allocated(popsnosob))           deallocate(popsnosob)
       if (allocated(pec))                 deallocate(pec)
       if (allocated(pecnosob))            deallocate(pecnosob)
       if (allocated(wl_cm))               deallocate(wl_cm)
       if (allocated(wl_cm_cubed))         deallocate(wl_cm_cubed)
       if (allocated(sob))                 deallocate(sob)
       if (allocated(sob_old))             deallocate(sob_old)
       if (allocated(wavelengthforspectrum)) deallocate(wavelengthforspectrum)
       if (allocated(broadspec))            deallocate(broadspec)
       if (allocated(cascade))            deallocate(cascade)

    end subroutine


end module colradfort