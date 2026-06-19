module input 
    !input 
    use types 
    implicit none
    !default input variables
    real(f64) :: temperature = 3000
    real(f64) :: density     = 3e5 
    !
    !   Velocity - gets used in both opacity and Gaussian spectra
    real(f64) :: velocityExpansionC     = 0.100_f64 
    real(f64) :: wlmin_nm               = 500.0_f64 
    real(f64) :: wlmax_nm               = 7000.0_f64
    integer   :: numwl                  = 10000
    !   Opacity
    logical   :: sobolev                = .false.
    logical   :: careful_la             = .false.
    real(f64) :: timeSinceExplosionDays = 29.00_f64
    real(f64) :: massElementSolar       = 0.005_f64
    real(f64) :: fractionOverride       = 0.0_f64
    real(f64) :: requiredLumo           = 1.0e37_f64
    character(len=300) ::  mode              = 'astro' !calculation modes 
    logical   :: floersHack = .false.
    logical   :: writeoutrates =.false.
    
    !
    character(len=300) :: adf04path =''
    character(len=256) :: error_msg

    integer   :: ioerror
    logical :: inputexists
    namelist /crm_input/ temperature, density, adf04path,sobolev, &
                         massElementSolar,velocityExpansionC,     &
                         fractionOverride,timeSinceExplosionDays, & 
                         wlmin_nm,wlmax_nm,numwl,mode,requiredLumo,& 
                         careful_la,floersHack,writeoutrates
                         
    contains
    subroutine getinput 

        inquire(file='crm_input',exist=inputexists)

        if (.not. inputexists) then 
            call defaultInput
            stop ' no input detected.'
        end if 
        !print*, 'debug - ',numwl
        open(1,file = 'crm_input')
        read(1,nml=crm_input,iostat=ioerror,iomsg=error_msg)

        if (ioerror.ne.0) then 
            write(0,*) error_msg
            !write(0,*) temperature
            !write(0,*) density
            !write(0,*) requiredLumo
            !write(0,*) mode
            !write(0,*) adf04path

            stop 'input io error'
        end if 
        close(1)
        !print*, 'debug - ',numwl

    end subroutine getinput

    subroutine defaultInput 
        print*," &crm_input temperature=3000, density = 3e5, adf04path='' &end \n"
    end subroutine defaultInput

end module input