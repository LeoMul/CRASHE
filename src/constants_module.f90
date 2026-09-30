module constants_module 
    use types 
    implicit none 
   real(f64), parameter :: kB_eV = 8.617333262e-5_f64
   real(f64), parameter :: coll_fac = 8.629e-6_f64
   real(f64), parameter :: pi = 4.0_f64*atan(1.0_f64)
   real(f64), parameter :: piFourOnThree = 4.0_f64*pi/3.0_f64
   real(f64), parameter :: sobconst = 1.0_f64/(8.0*pi)
   real(f64), parameter :: m_solar_grams = 1.989e+33_f64
   real(f64), parameter :: fwhmSigma = 1._f64/2.355_f64
   real(f64), parameter :: minusHalf = -0.5_f64
   real(f64), parameter :: oneOverSQRTTWOPI = 1._f64/(sqrt(2.0_f64*pi))
   real(f64), parameter :: hc_ergcm = 1.98644586e-16_f64 ! in erg cm
end module constants_module