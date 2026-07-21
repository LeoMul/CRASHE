module readadf04_module
    !this module reads an adf04file.
    !it also has a hacked routine to read data from Floers - for compairson.
    !in doing that - there is an AI generated routine for getting the van-regemorter collision strengths.
    use types
    use input,only: adfreadflag
    contains
        subroutine readadf04(filepath, numlevels, numtemps, ups, aval,statweight,energies, temps, wl_cm, wl_cm_cubed,atomicNumber,iq)
        
        !this routine is horribly slow. The scientific notation format
        !within these files is very hard to parse. 
        
        implicit none 
        integer,parameter :: maxNumLevels = 1000 
        integer,parameter   :: maxIter = 10000
        integer,parameter   :: maxNumTemps = 320

        character(len=*)  :: filepath
        logical           :: ex 
        integer           :: numLevels, numtemps
        !local variables
        character*320       :: templine
        integer             :: yy,zz,ii,jj,offset,lower,upper
        real(f64) :: myj, myenergy, ei, ej
        character*19        :: dummy
        character*1         :: dummy2

        !
        real(f64) , allocatable :: temps (:), logTemp(:)
        real(f64) , allocatable :: ups (:,:)
        real(f64) , allocatable :: aval(:)
        real(f64) , allocatable :: statweight(:)
        real(f64) , allocatable :: energies(:)
        real(f64), allocatable :: energiesTemp(:)
        real(f64), allocatable :: statweightTemp(:)
        real(f64), allocatable :: tempsTemp(:), wl_cm(:), wl_cm_cubed(:)
        integer             :: pp
        integer :: maxNumTransitions,iostat
        character*2 :: iel,IONTRM
        integer :: iq,atomicNumber,iq1 
        real(f64) :: fipot
        character*10 :: transitionFormat = '(2I4,A300)'
        character*22 :: levelFormat      = '(I5,24X,f4.1,1x,f21.4)'

        
        if (adfreadflag /=0) then 
            transitionFormat = '(2I5,A300)'
            levelFormat      = '(I5,40X,f4.1,1x,f21.4)'
        end if 


        !   
        inquire(file=filepath,exist=ex)
        if (.not.ex) then 
            write(*,*) 'file 1 = ',filepath,' not found'
            open(99,file='ifail')
            close(99)
            stop !'file 1 not found.'
        end if 
        open(1,file=filepath)
        !Borrowed from Martin O'Mullane
        571 FORMAT(A2,1X,I2,2I10,F15.0,1X,A2,1X)
        READ(1,571)IEL,IQ,atomicNumber,IQ1,FIPOT,IONTRM          !GET HEADER
        open (20,file='chargequick')
        write(20,*) IEL,iq 
        close(20)
        !read (1,'(A320)') templine
        !print*,templine
        !first find number of levels.

        allocate(energiesTemp(maxNumLevels),statweightTemp(maxNumLevels))

        do ii=1,maxIter
            read (1,levelFormat,iostat=iostat) yy,myj,myenergy
            print*, myj,myenergy
            !read (1,'(I5,24X,f4.1,1x,f21.4)',iostat=iostat) yy,myj,myenergy 
            energiesTemp(ii) = myenergy 
            statweightTemp(ii) = 2.0d0 * myj + 1.0d0
            !print*,yy,myj,myenergy
            if (yy.eq.-1) exit
        end do 
        numLevels = ii-1
        maxNumTransitions = (numLevels * (numLevels + 1)) / 2
        !print*,'Found ',ii-1,' atomic levels in file 1.' 

        allocate(energies(numLevels), statweight(numLevels)) 
        allocate(aval(maxNumTransitions))
        allocate(wl_cm(maxNumTransitions),wl_cm_cubed(maxNumTransitions))
        energies(1:numLevels) = energiesTemp(1:numlevels)
        statweight(1:numLevels) = statweightTemp(1:numlevels)

        do ii = 1, numLevels-1 
           ei = energies(ii)
           do jj = ii+1,numLevels
              ej = energies(jj)
              pp = upperTriangleIndexing(ii,jj,numLevels)
              wl_cm(pp) = 1.0_f64 / (ej - ei)
              wl_cm_cubed(pp) = wl_cm(pp)**3
           end do 
        end do 

        deallocate(statweightTemp,energiesTemp)

        energies(:) = energies(:) / 8065.54429d0 !convert to eV 

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !Read temps
        read(1,'(A17,A320)') dummy,templine
        !print*,templine
        !print*,temperatures2
        allocate(tempsTemp(maxNumTemps))
        !allocate(logTemp(maxNumTemps))

        call StripSpaces(templine)
        !allocate(temps(100))
        !print*,templine
        offset = 0
        do ii = 1,maxNumTemps
            dummy2 = templine(ii+offset:ii+offset)
            tempsTemp(ii) = a7toFloat(templine(ii+offset:ii+offset+6))
            offset = offset +6
            if (dummy2.eq.' ') then 
                exit 
            end if 
        end do
        !print*,'I have found ',ii-1, 'temperatures.'



        numTemps = ii-1
        allocate(temps(numTemps))
        temps(1:numtemps) = tempsTemp(1:numtemps)
        deallocate(tempsTemp)
        logTemp = log10(temps)
        !do ii = 1, numTemps
        !    write(2,*) temps(ii),logTemp(ii)
        !end do 
        allocate(ups(numTemps,maxNumTransitions))
        do ii = 1,maxNumTransitions
            !not sure if this will work with all compilre==
            read(1,transitionFormat) yy,zz,templine 
            if (yy.eq.-1) exit
            lower = min(zz,yy)
            upper = max(zz,yy)
            !pp = upper-1 + numLevels * (lower-1)
            pp = upperTriangleIndexing(lower,upper,numLevels)
            offset = 0
            aval(pp) = a8toFloat(templine(1+offset:1+offset+7))
            offset = offset + 7
            do jj = 2,numTemps+1
                ups(jj-1,pp) = a8toFloat(templine(jj+offset:jj+offset+7))
                offset = offset +7
            end do
            !upsinf(pp) = a8toFloat(templine(jj+offset:jj+offset+7))
        end do 
        !print*,'Found ',ii-1,' atomic transitions in file 1.' 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        close(1)

        contains 


            function a7toFloat(a7string) result(float)
            !converts A.BC+DE 
            ! or      A.BC-DE
            !to A.BC \times 10 ^{DE}
            real(f64) :: mantessa 
            integer :: exponent 
            real(f64) :: float
            character*7 :: a7string
            character*4 :: a4string
            character*3 :: a3string 
            integer :: ierr

            a4string = a7string(1:4)
            a3string = a7string(5:7)

            read(a4string,1) mantessa 
            read(a3string,2,iostat=ierr) exponent
            if (ierr.ne.0) then 
                write(0,*) a3string
                stop 'error in a7tofloat'
            end if 

            float = mantessa * (10.0d0**exponent)

            1 format(F4.1)
            2 format(I3)
        end function

        function a8toFloat(a8string) result(float)
            !converts A.BC+DE 
            ! or      A.BC-DE
            !to A.BC \times 10 ^{DE}
            !or:
            !converts
            !converts -A.BC+DE 
            ! or      -A.BC-DE
            !to -A.BC \times 10 ^{DE}
            real(f64) :: mantessa 
            integer :: exponent 
            real(f64) :: float
            character*8 :: a8string
            character*4 :: a4string
            character*3 :: a3string 
            character*1 :: a1string
            integer :: ierr

            a1string = a8string(1:1)
            a4string = a8string(2:5)
            a3string = a8string(6:8)



            read(a4string,1) mantessa 
            read(a3string,2,iostat=ierr) exponent
            if (ierr.ne.0) then 
                write(0,*) a3string 
                write(0,*) a8string

                stop 'error in a8tofloat'
            end if 

            float = mantessa * (10.0d0**exponent)

            if (a1string .eq. '-') then 
                float = -1.0d0 * float 
            end if 


            1 format(F4.1)
            2 format(I3)
        end function

        subroutine StripSpaces(string)
        !    https://stackoverflow.com/questions/27179549/removing-whitespace-in-string
            character(len=*) :: string
            integer :: stringLen 
            integer :: last, actual

            stringLen = len (string)
            last = 1
            actual = 1

            do while (actual < stringLen)
                if (string(last:last) == ' ') then
                    actual = actual + 1
                    string(last:last) = string(actual:actual)
                    string(actual:actual) = ' '
                else
                    last = last + 1
                    if (actual < last) &
                        actual = last
                endif
            end do

        end subroutine

        end subroutine
        
        function upperTriangleIndexing(lower,upper,rowsize)
        integer :: lower, upper,rowsize 
        integer :: upperTriangleIndexing 
        if(lower.ge.upper) stop 'error'
        upperTriangleIndexing = (lower - 1) * rowsize  + &
            upper - lower - (lower*(lower-1))/2
        end function 

        function triangleCount(L, rowsize)
        integer :: L, rowsize
        integer :: triangleCount
        ! T(L) = number of (lower,upper) pairs with lower <= L
        triangleCount = L*rowsize - (L*(L+1))/2
        end function

        subroutine inverseUpperTriangleIndexing(idx, rowsize, lower, upper)
        integer, intent(in)  :: idx, rowsize
        integer, intent(out) :: lower, upper
        integer :: L
        !integer :: triangleCount
        real(8) :: n, disc, Lreal

        n    = real(rowsize, 8)
        disc = (2.0d0*n - 1.0d0)**2 - 8.0d0*real(idx, 8)
        if (disc .lt. 0.0d0) disc = 0.0d0   ! guard against tiny negative from roundoff

        Lreal = ((2.0d0*n - 1.0d0) - sqrt(disc)) / 2.0d0

        ! initial guess, then correct for floating point error
        L = ceiling(Lreal - 1.0d-9)
        if (L .lt. 1) L = 1

        do while (triangleCount(L, rowsize) .lt. idx)
            L = L + 1
        end do
        do while (L .gt. 1)
            if (triangleCount(L-1, rowsize) .lt. idx) exit
            L = L - 1
        end do

        lower = L
        upper = L + idx - triangleCount(L-1, rowsize)

        end subroutine

        subroutine readhack(filepath, numlevels, numtemps, ups, aval,statweight,energies, temps, wl_cm, wl_cm_cubed,atomicNumber,iq)
            implicit none
            integer           :: numLevels, numtemps, maxNumTransitions
            character(len=*)  :: filepath

            integer :: atomicNumber ,iq ,pp ,iostat,iallowed
            
            real(f64) , allocatable :: temps (:)
            real(f64) , allocatable :: ups (:,:),ei,ej
            real(f64) , allocatable :: aval(:)
            real(f64) , allocatable :: statweight(:)
            real(f64) , allocatable :: energies(:)

            real(f64), allocatable :: wl_cm(:), wl_cm_cubed(:)
            real(f64) :: av
            integer :: ii , jj ,kk ,ll 
            iq = 0 
            atomicNumber = 98 
            numTemps = 2 
            allocate(temps(numtemps))

            temps(1) = 1000.0d0 
            temps(2) = 10000.0d0 

            open(98,file='energies.dat')
            read(98,*) numLevels
            
            maxNumTransitions = (numLevels * (numLevels + 1)) / 2
            !write(0,*) maxNumTransitions
            allocate(aval(maxNumTransitions))
            allocate(wl_cm(maxNumTransitions),wl_cm_cubed(maxNumTransitions))

            allocate(energies(numLevels), statweight(numLevels))

            allocate(ups(numTemps,maxNumTransitions))

            do ii = 1 , numLevels
                read(98, * ) jj, energies(ii),jj 
                statweight(ii) = 2 * jj + 1
                !write(0,*) jj, energies(ii),jj 
            end do 

            open(546456,file='av.dat',action='read')
            
            do ii = 1, numLevels-1 
            ei = energies(ii)
            do jj = ii+1,numLevels
                ej = energies(jj)
                !write(0,*)ei,ej
                pp = upperTriangleIndexing(ii,jj,numLevels)
                !write(0,*) pp 
                wl_cm(pp) = 1.0_f64 / (ej - ei)
                wl_cm_cubed(pp) = wl_cm(pp)**3

                ups(:,pp) = 100 * statweight(ii) * statweight(jj)
                !write(0,*) ups(:,pp)

                read(546456,'(2I5,1X,ES10.2,I2)') kk,ll,av,iallowed
                if (iallowed==1) call vanregemorter(ii, jj, energies, aval, wl_cm, statweight, temps, ups, numTemps, numLevels, pp)
                aval(pp) = av 
                !write(0,*) aval(pp)

            end do 
            end do 
        
            energies(:) = energies(:) / 8065.54429d0 !convert to eV 


        end subroutine

        subroutine vanregemorter(ii, jj, energies, aval, wl_cm, statweight, temps, ups, numTemps, numLevels, pp)
            !AI generated for quick tests and comparson with Floers. Needs to be reviewed.
            implicit none
            integer,  intent(in)    :: ii, jj, numTemps, numLevels, pp
            real(f64), intent(in)   :: energies(:), aval(:), wl_cm(:), statweight(:), temps(:)
            real(f64), intent(inout):: ups(:,:)     

            ! Constants (CGS)
            real(f64), parameter :: h        = 6.62607015d-27  ! erg·s
            real(f64), parameter :: c        = 2.99792458d10   ! cm/s
            real(f64), parameter :: kB       = 1.380649d-16    ! erg/K
            real(f64), parameter :: e_charge = 4.80326d-10     ! esu (statcoulombs)
            real(f64), parameter :: me       = 9.10938d-28     ! g
            real(f64), parameter :: Ry_erg   = 2.17987d-11     ! 1 Ry in erg
            real(f64), parameter :: pi       = 3.14159265358979d0       

            real(f64) :: dE_erg, f_osc, u, gbar, q_vr, sigma_vr
            real(f64) :: prefactor
            integer   :: kk     

            ! Energy gap in erg (energies stored in eV after conversion, but called before that
            ! conversion in the main loop, so still in cm^-1 here)
            dE_erg = (energies(jj) - energies(ii)) * h * c   ! ΔE = hν, from cm^-1      

            ! Absorption oscillator strength from Einstein A coefficient:
            !   f_ij = (m_e c / 8 pi^2 e^2) * (g_j / g_i) * lambda^2 * A_ji
            !        = (3 m_e c^3) / (8 pi^2 e^2 nu^2) * (g_j/g_i) * A_ji
            f_osc = (me * c**3) / (8.0d0 * pi**2 * e_charge**2) &
                  * (statweight(jj) / statweight(ii)) &
                  * aval(pp) * wl_cm(pp)**2     

            ! Van Regemorter prefactor: Q = (8 pi / sqrt(3)) * (e^2/(m_e c)) * (Ry/dE) * f * p(u)
            !   simplifies to: Q [cm^2] = C0 * (Ry/dE)^2 * f * gbar * (Ry/kT)^... 
            ! Standard form for the thermally-averaged rate coefficient → upsilon:
            !   Upsilon_ij = (g_i * 8 pi / sqrt(3)) * (Ry/dE)^2 * f_ij * Int_0^inf p(u) e^-u du
            ! where the Gaunt factor integral is tabulated.  For a simple implementation:
            !   gbar(u) ~ max(0.2, 0.276 * exp(u) * E1(u))    (Seaton 1962 fit)
            ! and Upsilon = (8 pi / sqrt(3)) * g_i * (Ry/dE)^2 * f * gbar_thermal       

            prefactor = (8.0d0 * pi / sqrt(3.0d0)) &
                      * statweight(ii) &
                      * (Ry_erg / dE_erg)**2 &
                      * f_osc       

            do kk = 1, numTemps
                ! Dimensionless energy ratio u = dE / kT
                u = dE_erg / (kB * temps(kk))       

                ! Thermally-averaged Gaunt factor gbar (Seaton 1962 approximation)
                gbar = max(0.2d0, 0.276d0 * exp(u) * expint1(u))        

                ups(kk, pp) = prefactor * gbar
            end do      

        end subroutine vanregemorter        


        ! Exponential integral E1(x) via rational approximation (Abramowitz & Stegun 5.1.56)
        real(f64) function expint1(x)
        !AI generated for quick tests and comparson with Floers. Needs to be reviewed.

            implicit none
            real(f64), intent(in) :: x
            real(f64) :: t      

            if (x <= 1.0d0) then
        ! Logarithmic series (A&S 5.1.53)
        expint1 = -log(x) - 0.57721566d0 &
                  + x*(1.0d0 + x*(-0.25d0 + x*(0.055555d0 + x*(-0.010416d0 + x*8.33d-4))))
            else
        ! Continued-fraction / rational approximation (A&S 5.1.56)
        t = 1.0d0 / x
        expint1 = exp(-x) * t * (1.0d0 + t*(-1.0d0 + t*(2.0d0 + t*(-6.0d0 + t*24.0d0)))) 
            end if      

        end function expint1

end module