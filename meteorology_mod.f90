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

USE chemf !contains the number of concentrations neq and the concentrations
USE aerosol_mod ! 

!-----------------------------------------------------------------------------------------
! Variable declaration
!-----------------------------------------------------------------------------------------
IMPLICIT NONE  ! This applies in all subprograms inside the module

PRIVATE

! TO DO: move the global variable to one module shared by every module (it could include the parameter module)
PUBLIC :: ua, va, theta, k_closure, k_closure_m, k_closure_h, Ri_num  ! basic meteorology variables
PUBLIC :: meteorology_init, surface_values, update_meteo  ! functions
public :: concentrations ! WARNING: could be part of the chemistry model?
public :: particle_conc

! Some constants
REAL(dp), PARAMETER :: lambda = 300.0_dp  ! maximum mixing length, meters
! REAL(dp), PARAMETER :: vonk = 0.4_dp      ! von Karman constant, dimensionless

! Meteorological variables
REAL(dp), DIMENSION(nz) :: ua, va, theta                      ! wind veocity components, potential temperature
REAL(dp), DIMENSION(nz) :: u_new, v_new, theta_new            ! temp var for time update 
REAL(dp), DIMENSION(nz-1) :: k_closure, L, h_half             ! K of K-theory, L for model 2
real(dp) :: dudz1, dudz2, dvdz1, dvdz2, dtdz1, dtdz2          ! calculs intermediates
REAL(dp), DIMENSION(nz-1) :: k_closure_m, k_closure_h, Ri_num ! K of K-theory, Richardson number, L for model 2
real(dp) :: f_m, f_h

! Mixing of the chemistry parameter (the concentrations)
real(dp), dimension(nz,neq):: concentrations, concentrations_new
real(dp), dimension(neq):: dcdz1, dcdz2

! aerosol mixing
real(dp), dimension(nz,nr_bins):: particle_conc, particle_conc_new
real(dp), dimension(nr_bins):: dpdz1, dpdz2

! For convenient
INTEGER :: I, J  ! used for loop
integer :: model


CONTAINS

  
  ! update the wind speed and the potential temperature
  subroutine update_meteo(model) ! CHECK: is dt the default time step for meteorological processes
    implicit none
    integer :: model
    ! boundary conditions
    u_new(1)=0.0
    v_new(1)=0.0
    u_new(nz)=ua(nz)
    v_new(nz)=va(nz)
    theta_new(1)=theta(1) ! the actual surface temperature updated by the previous call of surface_values
    theta_new(nz)=theta(nz)
    concentrations_new(nz,:)=0.0
    particle_conc_new(nz,:)=particle_conc(nz,:)
    ! TO DO: set the new particle concentration in the top layer to the previous state (steady state on nz)
    

    ! for every other levels
    if (model==1) then
       ! model version 1
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
       
    else if (model==2) then
       ! model version 2
       ! update K
       do I = 1,(nz-1)
          dudz1=(ua(I+1)-ua(I))/(h(I+1)-h(I))
          dvdz1=(va(I+1)-va(I))/(h(I+1)-h(I))
          k_closure(I)=(L(I)**2)*sqrt(dudz1**2+dvdz1**2) ![m^2 s^{-1}]
          ! k_closure_t(I,counter+1)
          Ri_num(I) = 0.0     ![] not used yet
       end do

       ! compute evolution of the wind speed and temperature
       do I = 2,(nz-1)
          ! u component
          dudz1=(ua(I+1)-ua(I))/(h(I+1)-h(I))
          dudz2=(ua(I)-ua(I-1))/(h(I)-h(I-1))
          u_new(I)=ua(I) + dt*(fcor*(va(I)-va(nz)) + (k_closure(I)*dudz1-k_closure(I-1)*dudz2)/(0.5*(h(I+1)-h(I-1))) )

          ! v component
          dvdz1=(va(I+1)-va(I))/(h(I+1)-h(I))
          dvdz2=(va(I)-va(I-1))/(h(I)-h(I-1))
          v_new(I)=va(I) + dt*(-fcor*(ua(I)-ua(nz)) + (k_closure(I)*dvdz1-k_closure(I-1)*dvdz2)/(0.5*(h(I+1)-h(I-1))) )

          ! potential temperature
          dtdz1=(theta(I+1)-theta(I))/(h(I+1)-h(I))
          dtdz2=(theta(I)-theta(I-1))/(h(I)-h(I-1))
          theta_new(I)= theta(I)  + dt*(k_closure(I)*dtdz1-k_closure(I-1)*dtdz2)/(0.5*(h(I+1)-h(I-1)))
       end do
    else
       ! model 3

       ! update K and Ri
       do I = 1,(nz-1)
          ! compute derivative
          dudz1=(ua(I+1)-ua(I))/(h(I+1)-h(I))
          dvdz1=(va(I+1)-va(I))/(h(I+1)-h(I))
          dtdz1=(theta(I+1)-theta(I))/(h(I+1)-h(I))
          
          ! Richardson number
          Ri_num(I) = (2.*g/(theta(I+1)+theta(I))) * dtdz1/(dudz1**2+dvdz1**2)
          
          ! K value
          k_closure(I)=(L(I)**2)*sqrt(dudz1**2+dvdz1**2) ![m^2 s^{-1}]

          ! correction terms
          if (Ri_num(I)<0.0) then
             f_m=sqrt((1.0-16.0*Ri_num(I)))
             f_h=f_m*sqrt(f_m)
          else if (Ri_num(I)>=0.2) then
             f_m=0.1
             f_h=f_m
          else
             f_m=max((1.0-5.0*Ri_num(I))**2,0.1)
             f_h=f_m
          end if

          ! correction of the K values
          k_closure_m(I)=f_m*k_closure(I)
          k_closure_h(I)=f_h*k_closure(I)  
       end do
       
       ! compute evolution of the wind speed and temperature
       do I = 2,(nz-1)
          ! u component
          dudz1=(ua(I+1)-ua(I))/(h(I+1)-h(I))
          dudz2=(ua(I)-ua(I-1))/(h(I)-h(I-1))
          u_new(I)=ua(I) + dt*(fcor*(va(I)-va(nz)) + (k_closure_m(I)*dudz1-k_closure_m(I-1)*dudz2)/(0.5*(h(I+1)-h(I-1))) )

          ! v component
          dvdz1=(va(I+1)-va(I))/(h(I+1)-h(I))
          dvdz2=(va(I)-va(I-1))/(h(I)-h(I-1))
          v_new(I)=va(I) + dt*(-fcor*(ua(I)-ua(nz)) + (k_closure_m(I)*dvdz1-k_closure_m(I-1)*dvdz2)/(0.5*(h(I+1)-h(I-1))) )

          ! potential temperature
          dtdz1=(theta(I+1)-theta(I))/(h(I+1)-h(I))
          dtdz2=(theta(I)-theta(I-1))/(h(I)-h(I-1))
          theta_new(I)= theta(I)  + dt*(k_closure_h(I)*dtdz1-k_closure_h(I-1)*dtdz2)/(0.5*(h(I+1)-h(I-1)))

          ! concentrations
          dcdz1=(concentrations(I+1,:)-concentrations(I,:))/(h(I+1)-h(I))
          dcdz2=(concentrations(I,:)-concentrations(I-1,:))/(h(I)-h(I-1))
          concentrations_new(I,:)= concentrations(I,:)  + dt*(k_closure_h(I)*dcdz1-k_closure_h(I-1)*dcdz2)/(0.5*(h(I+1)-h(I-1)))

          ! concentrations
          dpdz1=(particle_conc(I+1,:)-particle_conc(I,:))/(h(I+1)-h(I))
          dpdz2=(particle_conc(I,:)-particle_conc(I-1,:))/(h(I)-h(I-1))
          particle_conc_new(I,:)= particle_conc(I,:)  + dt*(k_closure_h(I)*dpdz1-k_closure_h(I-1)*dpdz2)/(0.5*(h(I+1)-h(I-1)))
       end do
       
       ! set the concentration flux toward the ground to 0.0
       concentrations_new(1,:)=concentrations_new(2,:)
       particle_conc_new(1,:)=particle_conc_new(2,:)
       
    end if

    ! update the arrays (wind speed and potential temperature)
    ua=u_new
    va=v_new
    theta=theta_new
    concentrations=concentrations_new
    particle_conc=particle_conc_new
    
    
  end subroutine update_meteo
  
  
  ! set the values of K and Ri
  subroutine calculate_k_and_ri(model)
    implicit none
    integer:: model
    
    !calculate the K(z) of the K-theory for closure problem
    if(model==1) then
       do I = 1,(nz-1)
          k_closure(I)=5.0_dp ![m^2 s^{-1}]
          Ri_num(I) = 0.0     ![] not used yet
       end do
    else if(model==2) then
       do I = 1,(nz-1)
          h_half(I)=0.5*(h(I)+h(I+1))
          L(I)=vonk*h_half(I)/(1.+(vonk*h_half(I)/lambda))
       end do
    else
       do I = 1,(nz-1)
          h_half(I)=0.5*(h(I)+h(I+1))
          L(I)=vonk*h_half(I)/(1.+(vonk*h_half(I)/lambda))
       end do
    end if
  end subroutine calculate_k_and_ri



!-----------------------------------------------------------------------------------------
! Set initial values for ua, va, theta.
!-----------------------------------------------------------------------------------------
  SUBROUTINE meteorology_init(model)
    implicit none
    integer:: model
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
    call calculate_k_and_ri(model)
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
