module Periodic_Table
    use types    
    implicit none
    ! Conversion factor: 1 atomic mass unit (u) in grams
    real(f64), parameter :: AMU_TO_GRAMS = 1.6605390666e-24_f64

    ! Standard Atomic Weights (g/mol or u)
    ! Data based on IUPAC standard values. 
    ! For unstable elements (brackets), the mass of the longest-lived isotope is used.
    
    !This array was AI generated. It agrees well with my desk mat.
    real(f64), parameter :: ATOMIC_MASSES(1:118) = [ &
        1.0080, 4.0026, 6.9410, 9.0122, 10.811, 12.011, 14.007, 15.999, 18.998, 20.180, & ! 1-10
        22.990, 24.305, 26.982, 28.085, 30.974, 32.065, 35.453, 39.948, 39.098, 40.078, & ! 11-20
        44.956, 47.867, 50.942, 51.996, 54.938, 55.845, 58.933, 58.693, 63.546, 65.380, & ! 21-30
        69.723, 72.630, 74.922, 78.971, 79.904, 83.798, 85.468, 87.620, 88.906, 91.224, & ! 31-40
        92.906, 95.950, 98.000, 101.07, 102.91, 106.42, 107.87, 112.41, 114.82, 118.71, & ! 41-50
        121.76, 127.60, 126.90, 131.29, 132.91, 137.33, 138.91, 140.12, 140.91, 144.24, & ! 51-60
        145.00, 150.36, 151.96, 157.25, 158.93, 162.50, 164.93, 167.26, 168.93, 173.05, & ! 61-70
        174.97, 178.49, 180.95, 183.84, 186.21, 190.23, 192.22, 195.08, 196.97, 200.59, & ! 71-80
        204.38, 207.20, 208.98, 209.00, 210.00, 222.00, 223.00, 226.00, 227.00, 232.04, & ! 81-90
        231.04, 238.03, 237.00, 244.00, 243.00, 247.00, 247.00, 251.00, 252.00, 257.00, & ! 91-100
        258.00, 259.00, 262.00, 267.00, 270.00, 271.00, 270.00, 277.00, 281.00, 281.00, & ! 101-110
        285.00, 286.00, 289.00, 289.00, 293.00, 293.00, 294.00, 294.00                  & ! 111-118
    ]

contains

    ! Returns mass of a single atom in grams
    function get_mass_grams(Z) result(m_grams)
        integer, intent(in) :: Z
        real(f64) :: m_grams
        
        if (Z >= 1 .and. Z <= 118) then
            m_grams = ATOMIC_MASSES(Z) * AMU_TO_GRAMS
        else
            m_grams = 0.0_f64
            write(*,*) "Warning: Atomic number Z=", Z, " out of range [1-118]."
            stop 
        end if
    end function get_mass_grams

    ! Returns the atomic weight in 'u' (same as g/mol)
    function get_atomic_weight(Z) result(aw)
        integer, intent(in) :: Z
        real(f64) :: aw
        if (Z >= 1 .and. Z <= 118) then
            aw = ATOMIC_MASSES(Z)
        else
            aw = 0.0_f64
        end if
    end function get_atomic_weight

end module Periodic_Table