program mycrm
    use input
    use colradfort
    use onion_module
    implicit none 

    open(6,file = 'crm.out')

    call getinput 

    call getadf04(adf04path,floersHack)

    call alloc(numwl)
    print*,mode
    if (mode .eq. 'astro') then

        call getAtomicDensityLocal(atomicDensity,numions, massElementSolar ,&
                                    atomicNumber,&
                                    velocityExpansionC,&
                                    0.0_f64,&
                                    fractionOverride,&
                                    timeSinceExplosionDays,&
                                    density)

        call colrad(temperature,& 
                    density,&
                    sobolev,&
                    timeSinceExplosionDays,&
                    atomicDensity,&
                    wlmin_nm, &
                    wlmax_nm, &
                    numwl , &
                    careful_la, &
                    writeoutrates,& 
                    velocityExpansionC,&
                    wavelengthforspectrum,&
                    broadspec,&
                    numions, broadmodedefault)

    else if (mode .eq. 'levelscan') then 
        call levelscan(temperature,density,careful_la,writeoutrates) 
    else if (mode .eq. 'masscontour') then 
        call masscontour (temperature,density, requiredlumo,careful_la,writeoutrates)
    else if (mode .eq. 'onion') then 
        call onion
    else 
        print*, ' Bad calculation mode requested. Check input. '
    end if

    close(6)
    call dealloc

end program

