program mycrm
    use input
    use colradfort
    use onion_module

    open(6,file = 'crm.out')

    call getinput 

    call getadf04(adf04path,floersHack)

    call alloc(numwl)
    print*,mode
    if (mode .eq. 'astro') then
        
        call colrad(temperature,& 
                    density,&
                    sobolev,&
                    velocityExpansionC,&
                    timeSinceExplosionDays,&
                    massElementSolar,&
                    fractionOverride,& 
                    wlmin_nm, &
                    wlmax_nm, &
                    numwl , &
                    careful_la, &
                    writeoutrates,& 
                    wavelengthforspectrum,&
                    broadspec)

    else if (mode .eq. 'levelscan') then 
        call levelscan(temperature,density,careful_la,writeoutrates) 
    else if (mode .eq. 'masscontour') then 
        call masscontour (temperature,density, requiredlumo,careful_la,writeoutrates)
    else if (mode .eq. 'onion') then 
        call onion
    else 
        print*, ' Bad calculation mode requested. Check input. '
    end if 

    call dealloc

end program

