!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Meteorology module
! - Compute eddy diffusivity coefficients
! - Compute wind fields
! - Compute theta field
! - Mix chemicals and aerosols
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

MODULE meteorology_mod

!-----------------------------------------------------------------------------------------
! Load modules
!-----------------------------------------------------------------------------------------
USE parameters_mod
USE time_mod
USE grid_mod

!-----------------------------------------------------------------------------------------
! Variable declaration
!-----------------------------------------------------------------------------------------
IMPLICIT NONE  ! This applies in all subprograms inside the module

PRIVATE
PUBLIC :: ua, va, theta, k_closure  ! basic meteorology variables
PUBLIC :: meteorology_init, surface_values, update_meteo  ! functions

! Some constants
REAL(dp), PARAMETER :: lambda = 300.0_dp  ! maximum mixing length, meters
REAL(dp), PARAMETER :: vonk = 0.4_dp      ! von Karman constant, dimensionless

! Meteorological variables
REAL(dp), DIMENSION(nz) :: ua, va, theta             ! wind veocity components, potential temperature
REAL(dp), DIMENSION(nz) :: u_new, v_new, theta_new   ! temp var for time update 
REAL(dp), DIMENSION(nz-1) :: k_closure, Ri_num       ! K of K-theory, Richardson number
real(dp) :: dudz1, dudz2, dvdz1, dvdz2, dtdz1, dtdz2 ! calculs intermediates

! For convenient
INTEGER :: I, J  ! used for loop


CONTAINS

  ! set the values of K and Ri
  subroutine calculate_k_and_ri()
    implicit none
    !calculate the K(z) of the K-theory for closure problem
    do I = 1,(nz-1)
       k_closure(I)=5.0_dp ![m^2 s^{-1}]
       Ri_num(I) = 0.0     ![] not used yet
    end do
  end subroutine calculate_k_and_ri

  
  ! update the wind speed and the potential temperature
  subroutine update_meteo() ! CHECK: is dt the default time step for meteorological processes
    implicit none
    
    ! boundary conditions
    u_new(1)=0.0
    v_new(1)=0.0
    u_new(nz)=ua(nz)
    v_new(nz)=va(nz)
    theta_new(1)=theta(1) ! the actual surface temperature updated by the previous call of surface_values
    theta_new(nz)=theta(nz)

    !for every other levels
    do I = 2,(nz-1)
       ! u component
       dudz1=(ua(I+1)-ua(I))/(h(I+1)-h(I))
       dudz2=(ua(I)-ua(I-1))/(h(I)-h(I-1))
       u_new(I)=ua(I) + dt*(fcor*(va(I)-va(nz)) + k_closure(I)*(dudz1-dudz2)/(0.5*(h(I+1)-h(I-1))) )

       ! v component
       dvdz1=(va(I+1)-va(I))/(h(I+1)-h(I))
       dvdz2=(va(I)-va(I-1))/(h(I)-h(I-1))
       v_new(I)=va(I) + dt*(-fcor*(ua(I)-ua(nz)) + k_closure(I)*(dvdz1-dvdz2)/(0.5*(h(I+1)-h(I-1))) )

       ! potential temperature
       dtdz1=(theta(I+1)-theta(I))/(h(I+1)-h(I))
       dtdz2=(theta(I)-theta(I-1))/(h(I)-h(I-1))
       theta_new(I)= theta(I)  + dt*k_closure(I)*(dtdz1-dtdz2)/(0.5*(h(I+1)-h(I-1))) 
    end do

    ! update the arrays (wind speed and potential temperature)
    ua=u_new
    va=v_new
    theta=theta_new
  end subroutine update_meteo
  
  




!-----------------------------------------------------------------------------------------
! Set initial values for ua, va, theta.
!-----------------------------------------------------------------------------------------
  SUBROUTINE meteorology_init()
    implicit none
    ! Wind velocity
    ua = 0.0_dp
    ua(nz) = 10.0_dp
    ua(2:nz-1) = (ua(nz)-ua(1)) * h(2:nz-1)/h(nz) + ua(1)

    va = 0.0_dp

    ! Potential temperature
    !theta(1) = 273.15_dp + 25.0_dp ! CHECK: shouldn't it be the actual surface temperature?
    theta = 273.15_dp + 25.0_dp
    theta(nz) = 273.15_dp + 30.0_dp
    !theta(2:nz-1) = (theta(nz)-theta(1)) * h(2:nz-1)/h(nz) + theta(1) ! assume to be a linear function of the height

    !K of the K-theory for closure
    call calculate_k_and_ri()
  END SUBROUTINE meteorology_init


!-----------------------------------------------------------------------------------------
! Get surface temperature at specific time.
!
! Note: can also get water concentrantion, in ppt, if modify this
! subroutine to also use column 8.
!
! Data is taken from:
! http://www.atm.helsinki.fi/~junninen/smartSearch/smartSearch.php
!-----------------------------------------------------------------------------------------
  SUBROUTINE surface_values(temperature, time)
    REAL(dp), INTENT(in)            :: time ! input, in seconds
    REAL(dp), INTENT(out)           :: temperature ! output, in Kelvin
    LOGICAL, SAVE                   :: first_time = .TRUE.
    REAL(dp), DIMENSION(8,50), SAVE :: surface_data
    REAL(dp), DIMENSION(50), SAVE   :: temperature_data
    REAL(dp), PARAMETER             :: seconds_in_day = 24*60*60
    REAL(dp), PARAMETER             :: seconds_in_30min = 30*60
    INTEGER                         :: index
    REAL(dp) :: time24h, time30min, time24plus15, temp1, temp2, x

    ! Only when called for the first time, read in data from file
    ! With this trick, we don't need to open the file in the main program
    IF (first_time) THEN
       OPEN(30, file='input/hyytiala_2011_8_10_t_h2o.dat', status='old')
       READ(30, *) surface_data
       temperature_data(1:50) = surface_data(7,1:50) ! in Celcius
       first_time = .FALSE.
    END IF

    time24h = MODULO(time, seconds_in_day) ! time modulo 24 hours
    time24plus15 = time24h + 15*60 ! time since 23:45 previous day
    time30min = MODULO(time24plus15, seconds_in_30min)
    index = 1 + FLOOR(time24plus15/seconds_in_30min)

    temp1 = temperature_data(index)
    temp2 = temperature_data(index + 1)
    x = time30min/seconds_in_30min

    ! linear interpolation between previous and next temperature data value
    temperature = temp1 + x*(temp2 - temp1) + 273.15_dp ! now in Kelvin

  END SUBROUTINE surface_values


END MODULE meteorology_mod
