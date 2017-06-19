!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Parameters module
!
! - All the parameters are defined in this module.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

MODULE parameters_mod

!-----------------------------------------------------------------------------------------
! Variable declaration
!-----------------------------------------------------------------------------------------
IMPLICIT NONE

PUBLIC

! Declare numbers in 64bit floating point
! http://en.wikipedia.org/wiki/Double_precision_floating-point_format
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15,300)

REAL(dp), PARAMETER :: PI   = 2*ASIN(1.0_dp)     ! [-], the constant pi
REAL(dp), PARAMETER :: g    = 9.81_dp            ! [m s-2], gravitational acceleration
REAL(dp), PARAMETER :: R    = 8.3144598_dp       ! [J mol-1 K-1], universal gas constant
REAL(dp), PARAMETER :: NA   = 6.022140857e23_dp  ! [molec mol-1], Avogadro's number 
REAL(dp), PARAMETER :: Mair = 28.96e-3_dp        ! [kg mol-1], mean molar mass of air
REAL(dp), PARAMETER :: kb   = 1.38064852e-23_dp  ! [m2 kg s-2 K-1], Boltzmann constant
REAL(dp), PARAMETER :: Cp = 1012.0_dp            ! [J kg-1 K-1], air specific heat at constant pressure,
                                                 ! in reality, it has slight temperature dependency
REAL(dp), PARAMETER :: Omega = 2*PI/(24.0_dp*60.0_dp*60.0_dp)  ! [rad s-1], Earth angular speed

! Latitude and longitude of Hyytiala:
REAL(dp), PARAMETER :: latitude_deg = 61.8455_dp  ! [degN]
REAL(dp), PARAMETER :: longitude_deg = 24.2833_dp  ! [degE]
REAL(dp), PARAMETER :: latitude = latitude_deg * PI/180.0_dp  ! [rad]
REAL(dp), PARAMETER :: longitude = longitude_deg * PI/180.0_dp  ! [rad]

REAL(dp), PARAMETER :: fcor = 2*Omega*SIN(latitude)  ! Coriolis parameter

! chemistry
real(dp):: O2, N2 ![cm^{-3}] concentrations
real(dp), parameter:: H2O=1.0d16                     ! [cm^{-3}] water vapor concentration
real(dp), parameter:: Dm=0.0538*0.001                ! [kg cm^{-2}] 
real(dp), parameter:: emi_fac=100.0*10.0**(-12)/3.6  ! [kg kg^{-1} s^{-1}]

real(dp), parameter:: a    = 0.0027
real(dp), parameter:: cl1  = 1.006
real(dp), parameter:: ct1  = 95000.0     ! [J mol^{-1}]
real(dp), parameter:: ct2  = 230000.0    ! [J mol^{-1}]
real(dp), parameter:: Ts   = 303.15      ! [K]
real(dp), parameter:: Tm   = 314.0       ! [K]
real(dp), parameter:: beta = 0.09        ! [K^{-1}]

real(dp), parameter:: iso_mass_mol=68.0  ! [g mol^{-1}]
real(dp), parameter:: alp_mass_mol=136.0 ! [g mol^{-1}]

real(dp), parameter:: p00=1.01325d5      ! [Pa]


! aerosol and mixing
REAL(dp), PARAMETER :: vonk = 0.4_dp      ! von Karman constant, dimensionless

END MODULE parameters_mod
