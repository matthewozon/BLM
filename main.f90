!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Main program
!
! - Simulate emissions and chemical reactions of gases, aerosol processes as well as 
!   transport of gases and aerosol particles within the planetary boundary layer with a
!   column model.
! - Check Fortran conventions at http://www.fortran90.org/src/best-practices.html
! - Check code conventions at
!   http://www.cesm.ucar.edu/working_groups/Software/dev_guide/dev_guide/node7.html
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

PROGRAM main

!-----------------------------------------------------------------------------------------
! Load modules
!-----------------------------------------------------------------------------------------
USE parameters_mod
USE time_mod
USE grid_mod

USE meteorology_mod ! this module contains all the varaibles used in the mixing model (wind speed, temperature and concentrations)

! contains thetime  evolution of the chemistry: chemistry_step(Conc,time1,time2,O2_in,N2_in,M_in,H2O_in,TEMP_in,exp_coszen)
USE chemf

!-----------------------------------------------------------------------------------------
! Variable declaration
!-----------------------------------------------------------------------------------------
IMPLICIT NONE

! Used for loops
INTEGER :: I, J

! model selection
integer :: model=3 ! REMINDER: chemistry only available for the model 3


! TO DO: move those variables where they belong
real(dp):: exp_coszen, Emi_iso, Emi_alp, M_air ! M_air is really similar to Mair... change the name
real(dp),dimension(nz)::air_pressure, TEMP
real(dp):: tmpEmiIso,tmpEmiAlp
logical, parameter:: chem_on=.true.
tmpEmiIso=0.0
tmpEmiAlp=0.0


!-----------------------------------------------------------------------------------------
! Initialization
!-----------------------------------------------------------------------------------------
CALL time_init()              ! initialize time, time step, date

CALL meteorology_init(model)  ! initialize ua, va and theta

! chemistry init
concentrations=0.0

CALL open_files()             ! open files for future use

CALL write_files(time)        ! write initial values



!-----------------------------------------------------------------------------------------
! Start main loop
!-----------------------------------------------------------------------------------------
DO WHILE (time <= time_end)
   !---------------------------------------------------------------------------------------
   ! Meteorology
   !---------------------------------------------------------------------------------------
   ! Set lower boundary condition
   CALL surface_values(theta(1), time)
   

  !---------------------------------------------------------------------------------------
  ! Chemistry
  !---------------------------------------------------------------------------------------
  ! Start to calculate chemical reactions only after some time to save the computation
  ! time
  IF (time >= time_start_chemistry) THEN
    ! Calculate chemistry every dt_chem
     IF ( MOD( NINT((time - time_start_chemistry)*1000.0_dp), NINT(dt_chem*1000.0_dp) ) == 0 ) THEN
        ! calculate exp_coszen
        exp_coszen = get_exp_coszen(time,daynumber,latitude)
        ! calculate temperature
        TEMP=theta-(g/Cp)*(h-h(1))
        ! calculate the pressure (Mair=P/(kT))
        air_pressure=barometric_law(p00, TEMP, h, nz)

        if (chem_on) then
           ! for each layer, compute the chemistry
           do I = 2,nz-1
              ! set emision values
              if (I==2) then
                 call emission_rate_alpha(TEMP(2))
                 call emission_rate_isoprene(TEMP(2),exp_coszen)
                 ! for the record
                 !write(*,*) Emi_iso, " ", Emi_alp
                 tmpEmiIso = Emi_iso
                 tmpEmiAlp = Emi_alp
              else
                 ! set emission to 0.0
                 Emi_alp = 0.0
                 Emi_iso = 0.0
              end if
              ! set the concentrations of Mair, O2 N2
              M_air=air_pressure(I)/(kb*TEMP(I))
              O2=0.21*M_air
              N2=0.78*M_air
              !reset the constants
              concentrations(I,1)  = 24. * M_air / 1.E9     ! O3 concentration
              concentrations(I,4)  = 0.0                    ! REST
              concentrations(I,5)  = 0.2 * M_air / 1.E9     ! NO2
              concentrations(I,6)  = 0.07 * M_air / 1.E9    ! NO
              concentrations(I,9)  = 100. * M_air / 1.E9    ! CO
              concentrations(I,10) = 0.0                    ! CO2
              concentrations(I,11) = 1759. * M_air / 1.E9   ! CH4
              !concentrations(I,13) = 2. * M_air / 1.E9      ! C5H8  ! remove this one at some point
              concentrations(I,20) = 0.5 * M_air / 1.E9     ! SO2
              !concentrations(I,22) = 0.0                    ! H2SO4 ! remove this one at some point
              !concentrations(I,23) = 2. * M_air / 1.E9      ! alpha-p
              !concentrations(I,25) = 0.0                    ! ELVOC ! remove this one at some point
              CALL chemistry_step(concentrations(I,:),time,time+dt,O2,N2,M_air,H2O,Emi_iso,Emi_alp,TEMP(I),exp_coszen)
           end do
        else
           call emission_rate_alpha(TEMP(2))
           call emission_rate_isoprene(TEMP(2),exp_coszen)
           write(*,*) Emi_iso, " ", Emi_alp
           tmpEmiIso = Emi_iso
           tmpEmiAlp = Emi_alp
        end if
    END IF  ! time - time_start_chemistry, dt_chem
  END IF  ! time >= time_start_chemistry

  !---------------------------------------------------------------------------------------
  ! Aerosol
  !---------------------------------------------------------------------------------------
  ! Start to calculate aerosol processes only after some time to save the computation time
  IF (time >= time_start_aerosol) THEN  ! start to calculate aerosol after time_start_aerosol
    ! Calculate aerosol every dt_aero
    IF ( MOD( NINT((time - time_start_aerosol)*1000.0_dp), NINT(dt_aero*1000.0_dp) ) == 0 ) THEN
    END IF  ! time - time_start_aerosol, dt_aero
  END IF  ! time >= time_start_aerosol
 

  ! it's probably already done it the meteo updateforce an open boundary and null flux toward the ground
  concentrations(nz,:)=0.0
  concentrations(1,:)=concentrations(2,:)
  
  ! update wind velocity and potential temperature
  call update_meteo(model)
 

  !---------------------------------------------------------------------------------------
  ! Ending loop actions
  !---------------------------------------------------------------------------------------
  ! Advance to next time step
  time = time + dt

  ! Write data to files and time infomation to screen
  IF ( MOD( NINT((time - time_start)*1000.0_dp), NINT(dt_output*1000.0_dp) ) == 0 ) THEN
    WRITE(*, '(a8, f8.3, a6)') 'time = ', time/one_hour, '  hours'
    CALL write_files(time)
  END IF

  ! Count loop number
  counter = counter + 1

END DO

!-----------------------------------------------------------------------------------------
! Finalization
!-----------------------------------------------------------------------------------------
! Close all the opened files
CALL close_files()

! Count total time steps
WRITE(*,*) counter, 'time steps'


CONTAINS


!-----------------------------------------------------------------------------------------
! Open files.
!-----------------------------------------------------------------------------------------
SUBROUTINE open_files()
  CHARACTER(255), PARAMETER :: outdir = 'output'

  OPEN(11, FILE = TRIM(ADJUSTL(outdir))//'/time.dat'  , STATUS = 'REPLACE', ACTION = 'WRITE')
  OPEN(12, FILE = TRIM(ADJUSTL(outdir))//'/h.dat'     , STATUS = 'REPLACE', ACTION = 'WRITE')
  OPEN(13, FILE = TRIM(ADJUSTL(outdir))//'/ua.dat'    , STATUS = 'REPLACE', ACTION = 'WRITE')
  OPEN(14, FILE = TRIM(ADJUSTL(outdir))//'/va.dat'    , STATUS = 'REPLACE', ACTION = 'WRITE')
  OPEN(15, FILE = TRIM(ADJUSTL(outdir))//'/theta.dat' , STATUS = 'REPLACE', ACTION = 'WRITE')
  OPEN(16, FILE = TRIM(ADJUSTL(outdir))//'/kt.dat'    , STATUS = 'REPLACE', ACTION = 'WRITE')
  OPEN(17, FILE = TRIM(ADJUSTL(outdir))//'/kmt.dat'   , STATUS = 'REPLACE', ACTION = 'WRITE')
  OPEN(18, FILE = TRIM(ADJUSTL(outdir))//'/kht.dat'   , STATUS = 'REPLACE', ACTION = 'WRITE')
  OPEN(19, FILE = TRIM(ADJUSTL(outdir))//'/Ri.dat'    , STATUS = 'REPLACE', ACTION = 'WRITE')

  ! chemistry
  !OPEN(20, FILE = TRIM(ADJUSTL(outdir))//'/concentration.dat',status='replace',action='write')
  !OPEN(21, FILE = TRIM(ADJUSTL(outdir))//'/radiation.dat'    ,status='replace',action='write')
  OPEN(22, FILE = TRIM(ADJUSTL(outdir))//'/emi_iso.dat'       ,status='replace',action='write')
  OPEN(23, FILE = TRIM(ADJUSTL(outdir))//'/emi_alp.dat'       ,status='replace',action='write')

  OPEN(24, FILE = TRIM(ADJUSTL(outdir))//'/OH.dat'          ,status='replace',action='write')
  OPEN(25, FILE = TRIM(ADJUSTL(outdir))//'/HO2.dat'         ,status='replace',action='write')
  OPEN(26, FILE = TRIM(ADJUSTL(outdir))//'/H2SO4.dat'       ,status='replace',action='write')
  OPEN(27, FILE = TRIM(ADJUSTL(outdir))//'/isoprene.dat'    ,status='replace',action='write')
  OPEN(28, FILE = TRIM(ADJUSTL(outdir))//'/alpha.dat'       ,status='replace',action='write')
  OPEN(29, FILE = TRIM(ADJUSTL(outdir))//'/ELVOC.dat'       ,status='replace',action='write')
END SUBROUTINE open_files


!-----------------------------------------------------------------------------------------
! Save data to files.
!-----------------------------------------------------------------------------------------
SUBROUTINE write_files(time)
  REAL(dp) :: time  ! current time
  CHARACTER(255) :: outfmt0, outfmt1, outfmt2, outfmt3
  ! Get output format for arrays with nz layers
  WRITE(outfmt0, '(a, i3, a)') '(', nz, 'es25.16e3)'
  WRITE(outfmt1, '(a, i3, a)') '(', nz-1, 'es25.16e3)'
  !WRITE(outfmt2, '(a, i3, a)') '(', neq, 'es25.16e3)'
  WRITE(outfmt3, '(a, i3, a)') '(', 1, 'es25.16e3)'

  ! Only save h one time at the beginning
  IF (time == time_start) THEN
    WRITE(12, *) h
  END IF

  WRITE(11, '(f8.4)') time/(24*one_hour)   ! [day], time
  WRITE(13, outfmt0) ua                    ! [m s-1], u wind
  WRITE(14, outfmt0) va                    ! [m s-1], v wind
  WRITE(15, outfmt0) theta                 ! [K], potential temperature
  WRITE(16, outfmt1) k_closure             ! [m^2 s^{-1}], k_closure K
  WRITE(17, outfmt1) k_closure_m           ! [m^2 s^{-1}], k_closure K_m
  WRITE(18, outfmt1) k_closure_h           ! [m^2 s^{-1}], k_closure K_h
  WRITE(19, outfmt1) Ri_num                ! [], Richardson number Ri
  !WRITE(20, outfmt2) Conc                  ! [m^{-3}], concentrations
  !WRITE(21, outfmt3) exp_coszen            ! [], radiation
  WRITE(22, outfmt3) tmpEmiIso             ! [], radiation
  WRITE(23, outfmt3) tmpEmiAlp             ! [], radiation
  WRITE(24, outfmt0) concentrations(:,3)   ! [cm^{-3}...hopfully], OH concentration
  WRITE(25, outfmt0) concentrations(:,8)   ! [cm^{-3}...hopfully], HO2 concentration
  WRITE(26, outfmt0) concentrations(:,21)  ! [cm^{-3}...hopfully], H2SO4 concentration
  WRITE(27, outfmt0) concentrations(:,13)  ! [cm^{-3}...hopfully], isoprene concentration
  WRITE(28, outfmt0) concentrations(:,23)  ! [cm^{-3}...hopfully], alpha-pinene concentration
  WRITE(29, outfmt0) concentrations(:,25)  ! [cm^{-3}...hopfully], ELVOC concentration
  
  
END SUBROUTINE write_files


!-----------------------------------------------------------------------------------------
! Close opened files
!-----------------------------------------------------------------------------------------
SUBROUTINE close_files()
  CLOSE(11)
  CLOSE(12)
  CLOSE(13)
  CLOSE(14)
  CLOSE(15)
  CLOSE(16)
  CLOSE(17)
  CLOSE(18)
  CLOSE(19)
  !CLOSE(20)
  !CLOSE(21)
  CLOSE(22)
  CLOSE(23)
  CLOSE(24)
  CLOSE(25)
  CLOSE(26)
  CLOSE(27)
  CLOSE(28)
  CLOSE(29)
END SUBROUTINE close_files

  !-------------------------------------------------------
  ! emission rates
  !-------------------------------------------------------
  subroutine emission_rate_isoprene(temp,exp_cos)
    !implicit none
    real(dp)::temp,exp_cos
    real(dp):: Cl,Ct,PAR 
    PAR=1000.0*exp_cos
    Cl=a*cl1*PAR/sqrt(1+a**2*PAR**2)
    Ct=exp(ct1*(temp-Ts)/(R*temp*Ts))/(1+exp(ct2*(temp-Tm)/(R*temp*Ts)))
    Emi_iso=Dm*emi_fac*Cl*Ct*NA/iso_mass_mol ! 10.0**(20)*
  end subroutine emission_rate_isoprene

  subroutine emission_rate_alpha(temp)
    !implicit none
    real(dp)::temp
    Emi_alp=Dm*emi_fac*exp(beta*(temp-Ts))*NA/alp_mass_mol ! 10.0**(20)*/2.0
  end subroutine emission_rate_alpha
  

!---------------------------------------------------------------------------------------
! computation of the light exposure
!---------------------------------------------------------------------------------------
REAL(dp) FUNCTION get_hourangle(time)
    REAL(dp), INTENT(in) :: time
    REAL(dp), PARAMETER :: one_day = 24*one_hour
    get_hourangle = MODULO(time,one_day)/one_day * 2 * pi - pi
  END FUNCTION get_hourangle

  REAL(dp) FUNCTION solar_zenith_angle(hourangle,daynumber,latitude)
    ! http://en.wikipedia.org/wiki/Solar_elevation_angle
    ! http://en.wikipedia.org/wiki/Position_of_the_Sun
    INTEGER, INTENT(in) :: daynumber
    REAL(dp), INTENT(in) :: hourangle,latitude
    REAL(dp) :: declination,elevation
    REAL(dp), PARAMETER :: to_rad = pi/180.

    declination = -23.44 * to_rad * COS(2 * pi * (daynumber + 10)/365.)
    elevation = COS(hourangle)*COS(declination)*COS(latitude) &
      + SIN(declination)*SIN(latitude)
    solar_zenith_angle = pi/2. - elevation
    ! Notes:
    ! - Not tested near equador or on the southern hemisphere.
    ! - solar_zenith_angle can be larger than pi/2, it just means
    !   the sun is below horizon.
    ! - solar_zenith_angle assumes time is in local solar time, which
    !   is usually not exactly true
  END FUNCTION solar_zenith_angle

  REAL(dp) FUNCTION get_exp_coszen(time,daynumber,latitude)
    REAL(dp), INTENT(in) :: time,latitude
    INTEGER, INTENT(in) :: daynumber
    REAL(dp) :: hourangle,zenith,coszen
    hourangle = get_hourangle(time)
    zenith = solar_zenith_angle(hourangle,daynumber,latitude)
    coszen = COS(zenith)
    IF (coszen > 0) THEN  ! sun is above horizon
       get_exp_coszen = EXP(-0.575/coszen)
    ELSE
       get_exp_coszen = 0
    ENDIF
  END FUNCTION get_exp_coszen


  !----------------------------------------------------------
  ! compute the pressure for each layer
  !----------------------------------------------------------
  FUNCTION barometric_law(p00, tempK, h, nz) RESULT(p)
    INTEGER, INTENT(IN) :: nz
    REAL(dp), INTENT(IN) :: p00, tempK(nz), h(nz)
    REAL(dp) :: p(nz)
    REAL(dp) :: dh(nz)
    REAL(dp), PARAMETER :: g = 9.81_dp, Md=28.96_dp*10.0**(-3), Rgas=8.3144598_dp
    INTEGER :: I
  
    dh(2:nz) = h(2:nz) - h(1:nz-1)
  
    p(1) = p00
    DO I=2, nz
      p(I) = p(I-1)*EXP(-Md*g/(Rgas*(tempK(I-1)+tempK(I))/2.0d0)*dh(I))
    END DO
  END FUNCTION barometric_law


  
END PROGRAM main
