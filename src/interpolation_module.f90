module interpolation_module 
    use types
    implicit none
    contains 
    subroutine interpolate_upsilons(ntran, ntemps_adf04, temps_adf04, temp_req, ups_adf04, ups_interp)
        implicit none
        integer :: ntran 
        integer :: ntemps_adf04 
        real(f64):: ups_adf04 (ntemps_adf04,ntran)
        real(f64):: ups_interp(ntran)
        real(f64):: temps_adf04(ntemps_adf04)
        real(f64):: temp_req

        !local Variables
        real(f64) :: log_temps_adf04(ntemps_adf04)
        real(f64) :: log_temp_req
        real(f64) :: yy(ntemps_adf04+1)
        integer  :: ii  

        ! 
        ! Safety check
        if (temp_req < minval(temps_adf04)) then 
            stop ' Requested temperature below minimum adf04 temperature.'
        end if 
        if (temp_req > maxval(temps_adf04)) then 
            stop ' Requested temperature above maximum adf04 temperature.'
        end if 
        !
        !Take logs for easier interpolation 
        log_temps_adf04 = log10(temps_adf04)
        log_temp_req   = log10(temp_req  )

        do ii = 1, ntran 
            call spline(log_temps_adf04,ups_adf04(:,ii),ntemps_adf04,0.0d0,0.0d0,yy)
            call splint(log_temps_adf04,ups_adf04(:,ii),yy,ntemps_adf04,log_temp_req,ups_interp(ii))
        end do 

    end subroutine
!
SUBROUTINE spline(x,y,n,yp1,ypn,y2)
    !Numerical recipes fortran 77 - originally by William H. Press
    !https://github.com/wangvei/nrf77/blob/master/spline.f - Jon Lighthall
    INTEGER n,NMAX
    DOUBLE PRECISION yp1,ypn,x(n),y(n),y2(n)
    PARAMETER (NMAX=500)
    INTEGER i,k
    DOUBLE PRECISION p,qn,sig,un,u(NMAX)
    !print*,'hello',n,yp1,ypn
    if (yp1.gt..99d30) then
      y2(1)=0.d0
      u(1)=0.d0
    else
      y2(1)=-0.5d0
      u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    endif
    do 11 i=2,n-1
      sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
      p=sig*y2(i-1)+2.d0
      y2(i)=(sig-1.d0)/p
      u(i)=(6.d0*((y(i+1)-y(i))/(x(i+ 1)                        &
     -x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*   &
      u(i-1))/p
11    continue
    if (ypn.gt..99d30) then
      qn=0.d0
      un=0.d0
    else
      qn=0.5d0
      un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    endif
    y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
    do 12 k=n-1,1,-1
      y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
    return
    END SUBROUTINE
!
    SUBROUTINE splint(xa,ya,y2a,n,x,y)
    !Numerical recipes fortran 77 - originally by William H. Press
    !https://github.com/wangvei/nrf77/blob/master/splint.f - Jon Lighthall
    INTEGER n
    DOUBLE PRECISION x,y,xa(n),y2a(n),ya(n)
    INTEGER k,khi,klo
    DOUBLE PRECISION a,b,h
    klo=1
    khi=n
1     if (khi-klo.gt.1) then
      k=(khi+klo)/2
      if(xa(k).gt.x)then
        khi=k
      else
        klo=k
      endif
    goto 1
    endif
    h=xa(khi)-xa(klo)
!    if (h.eq.0.d0) print *, ' bad xa input in splint'
    if (abs(h).lt.1d-30) print *, ' bad xa input in splint' !lpm compiler warning fix.

    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))* &
      (h**2)/6.d0
    return
    END  SUBROUTINE

end module interpolation_module