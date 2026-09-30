module interpolation_module
   use types
   use constants_module
   use atomicdata_module
   use readadf04_module,only: upperTriangleIndexing
   implicit none
contains

   subroutine qrates_from_ups(qrateup,qratedown,upsilon, gl, gu, el, eu, sqrttemp,kT)
      implicit none 
      real(f64),intent(in)   :: upsilon,gl,gu,el,eu,sqrttemp,kt
      real(f64),intent(out)  :: qrateup,qratedown
      
      real(f64) :: xx 
      xx = (eu - el) / kT

      qratedown = coll_fac * upsilon / sqrttemp

      if (xx < 0.20_f64) then 
         qrateup = qratedown * (1.0 - xx) / gl 
      else if (xx < 700.0_f64) then 
         qrateup = qratedown * exp(-xx) /  gl 
      else 
         qrateup = 0.0_f64 
      end if 
      qratedown = qratedown / gu

   end subroutine

   subroutine interpolate_upsilons_calc_rates(temp_req)
      implicit none
      real(f64) :: temp_req
      real(f64) :: upsinterp
      real(f64) :: ei,ej,gi,gj
      real(f64) :: roottemp 
      real(f64) :: KT
      !local Variables
      real(f64) :: log_temps_adf04(numTemps)
      real(f64) :: log_temp_req
      real(f64) :: yy(numTemps + 1)
      integer  :: ii,jj,tt

      !integer :: pp 

      !
      ! Safety check
      if (temp_req < minval(temps)) then
         stop ' Requested temperature below minimum adf04 temperature.'
      end if
      if (temp_req > maxval(temps)) then
         stop ' Requested temperature above maximum adf04 temperature.'
      end if
      !

      if (.not. allocated(qup))   allocate(qup(ntran))
      if (.not. allocated(qdown)) allocate(qdown(ntran))

      !Take logs for easier interpolation
      log_temps_adf04 = log10(temps)
      log_temp_req = log10(temp_req)

      KT = kB_eV * temp_req
      roottemp = sqrt(temp_req)
      tt = 1 
      do ii = 1, numLevels-1
         gi = statweight(ii)
         ei = energies(ii)
         do jj = ii+1, numLevels
            gj = statweight(jj)
            ej = energies(jj)
            call spline(log_temps_adf04, ups(:, tt), numTemps,     0.0d0, 0.0d0, yy)
            call splint(log_temps_adf04, ups(:, tt), yy, numTemps, log_temp_req, upsinterp)
            call qrates_from_ups(qup(tt), qdown(tt), upsinterp,gi,gj,ei,ej,roottemp,KT)
            !pp = upperTriangleIndexing(ii,jj,numLevels)
            !write(0,*) qup(tt), qdown(tt)
            tt = tt+1

         end do 

      end do

   end subroutine
!
   SUBROUTINE spline(x, y, n, yp1, ypn, y2)
      !Numerical recipes fortran 77 - originally by William H. Press
      !https://github.com/wangvei/nrf77/blob/master/spline.f - Jon Lighthall
      INTEGER n, NMAX
      DOUBLE PRECISION yp1, ypn, x(n), y(n), y2(n)
      PARAMETER(NMAX=500)
      INTEGER i, k
      DOUBLE PRECISION p, qn, sig, un, u(NMAX)
      !print*,'hello',n,yp1,ypn
      if (yp1 .gt. .99d30) then
         y2(1) = 0.d0
         u(1) = 0.d0
      else
         y2(1) = -0.5d0
         u(1) = (3.d0/(x(2) - x(1)))*((y(2) - y(1))/(x(2) - x(1)) - yp1)
      end if
      do 11 i = 2, n - 1
         sig = (x(i) - x(i - 1))/(x(i + 1) - x(i - 1))
         p = sig*y2(i - 1) + 2.d0
         y2(i) = (sig - 1.d0)/p
         u(i) = (6.d0*((y(i + 1) - y(i))/(x(i + 1) &
                                          - x(i)) - (y(i) - y(i - 1))/(x(i) - x(i - 1)))/(x(i + 1) - x(i - 1)) - sig* &
                 u(i - 1))/p
11       continue
         if (ypn .gt. .99d30) then
            qn = 0.d0
            un = 0.d0
         else
            qn = 0.5d0
            un = (3.d0/(x(n) - x(n - 1)))*(ypn - (y(n) - y(n - 1))/(x(n) - x(n - 1)))
         end if
         y2(n) = (un - qn*u(n - 1))/(qn*y2(n - 1) + 1.d0)
         do 12 k = n - 1, 1, -1
            y2(k) = y2(k)*y2(k + 1) + u(k)
12          continue
            return
            END SUBROUTINE
!
            SUBROUTINE splint(xa, ya, y2a, n, x, y)
               !Numerical recipes fortran 77 - originally by William H. Press
               !https://github.com/wangvei/nrf77/blob/master/splint.f - Jon Lighthall
               INTEGER n
               DOUBLE PRECISION x, y, xa(n), y2a(n), ya(n)
               INTEGER k, khi, klo
               DOUBLE PRECISION a, b, h
               klo = 1
               khi = n
1              if (khi - klo .gt. 1) then
                  k = (khi + klo)/2
                  if (xa(k) .gt. x) then
                     khi = k
                  else
                     klo = k
                  end if
                  goto 1
               end if
               h = xa(khi) - xa(klo)
!    if (h.eq.0.d0) print *, ' bad xa input in splint'
               if (abs(h) .lt. 1d-30) print *, ' bad xa input in splint' !lpm compiler warning fix.

               a = (xa(khi) - x)/h
               b = (x - xa(klo))/h
               y = a*ya(klo) + b*ya(khi) + ((a**3 - a)*y2a(klo) + (b**3 - b)*y2a(khi))* &
                   (h**2)/6.d0
               return
            END SUBROUTINE

end module interpolation_module
