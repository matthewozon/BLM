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

! include the aerosol module
USE aerosol_mod

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
real(dp),dimension(nz)::air_pressure, TEMP, PN_atm, PM_atm
real(dp):: tmpEmiIso,tmpEmiAlp
logical, parameter:: chem_on=.true.
logical, parameter:: aerosol_on=.true.
integer:: aerosol_init_on
real(dp), dimension(nr_bins) :: layer1_particles
real(dp), dimension(nr_cond,nz) :: cond_sink
real(dp) :: DSWF

tmpEmiIso=0.0
tmpEmiAlp=0.0
aerosol_init_on=1



!-----------------------------------------------------------------------------------------
! Initialization
!-----------------------------------------------------------------------------------------
CALL time_init()              ! initialize time, time step, date

CALL meteorology_init(model)  ! initialize ua, va and theta

! chemistry init
concentrations=0.0

! initialize  aerosol ! only the first layer contains particle at the beginning... or not?
do I = 1,nz
   CALL Aerosol_init(diameter, particle_mass, particle_volume, particle_conc(I,:), &
        particle_density, nucleation_coef, molecular_mass, molar_mass, &
        molecular_volume, molecular_dia, mass_accomm, cond_sink(:,I))
   PN_atm(I)=PN
   PM_atm(I)=PM
end do

layer1_particles=particle_conc(1,:)


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
              M_air=10.0**(-6)*air_pressure(I)/(kb*TEMP(I)) !
              ! write(*,*) air_pressure(I), " " , TEMP(I), " ", air_pressure(I)/(kb*TEMP(I)), " ", M_air
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
              !write(*,*) "O2 ", O2, ", N2 ", N2, ", O3 " , concentrations(I,1), ", H2O ", H2O 
              CALL chemistry_step(concentrations(I,:),time,time+dt_chem,O2,N2,M_air,H2O,Emi_iso,Emi_alp,cond_sink(1,I),cond_sink(2,I),TEMP(I),exp_coszen)
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
     if (aerosol_on) then
        if (aerosol_init_on==1) then
           aerosol_init_on=0
           do I = 1,nz
              call init_part_conc(particle_conc(nz,:), diameter)
              PN_atm(I)=PN
              PM_atm(I)=PM
           end do
        end if
        
        ! Calculate aerosol every dt_aero
        IF ( MOD( NINT((time - time_start_aerosol)*1000.0_dp), NINT(dt_aero*1000.0_dp) ) == 0 ) THEN
           ! compute temperature and pressure for each layer
           TEMP=theta-(g/Cp)*(h-h(1))
           ! calculate the pressure (Mair=P/(kT))
           air_pressure=barometric_law(p00, TEMP, h, nz)
           !do I =1,nz-1
           !   write(*,*) sum(particle_conc(I,:))
           !   write(*,*) sum(particle_conc(I,:)*particle_mass)
           !end do

           ! dry deposition
           DSWF = 6D2 * exp_coszen
           call dry_dep_velocity(diameter,particle_density,TEMP(2),air_pressure(2),DSWF, & 
                0.5_dp*(Ri_num(1)+Ri_num(2)),sqrt(ua(2)**2+va(2)**2), mass_accomm, v_dep, vd_SO2, vd_O3, vd_HNO3)
           ! write(*,*) v_dep
           ! particle_conc(2,:) = particle_conc(2,:)*exp(-v_dep/h(2)*dt_aero)
           ! particle_conc(1,:) = particle_conc(2,:)
           !concentrations(2,20) = concentrations(2,20)*exp(-vd_SO2/h(2)*dt_aero)
           !concentrations(1,20) = concentrations(2,20)
           !concentrations(2,1) = concentrations(2,1)*exp(-vd_O3/h(2)*dt_aero)
           !concentrations(1,1) = concentrations(2,1)
           !concentrations(2,17) = concentrations(2,17)*exp(-vd_HNO3/h(2)*dt_aero)
           !concentrations(1,17) = concentrations(2,17)
           
           do I =1,nz-1 ! the first layer shouldn't be computed because it is equal to the second layer (but for plotting purposes, I will keep the first layer computation)
              ! set the concentrations of Mair
              ! M_air=(air_pressure(I)/(kb*TEMP(I)))/(NA*140000.0)
              ! write(*,*) Mair
              !pause
              ! vapour concentration
              cond_vapour(1)=concentrations(I,21)*1.0D6 ! H2SO4
              cond_vapour(2)=concentrations(I,25)*1.0D6 ! ELVOC

              ! nucleation
              call Nucleation(dt_aero,nucleation_coef,cond_vapour(1),particle_conc(I,:))

              ! coagulation
              call Coagulation(dt_aero, particle_conc(I,:), diameter, TEMP(I),air_pressure(I),particle_mass)

              ! condensation ! why a difference arise in the particle number, particle mass
              call Condensation(dt_aero, TEMP(I), air_pressure(I), mass_accomm, molecular_mass, &
                   molecular_volume, molar_mass, molecular_dia, particle_mass, particle_volume, &
                   particle_conc(I,:), diameter, cond_vapour, cond_sink(:,I))

              if (I==2) then
                 layer1_particles=particle_conc(I,:)
                 !write(*,*) layer1_particles
                 !pause
              end if
              call PNandPM(particle_conc(I,:))
              PN_atm(I)=PN
              PM_atm(I)=PM
           end do
        end IF     
     END IF  ! time - time_start_aerosol, dt_aero
  END IF  ! time >= time_start_aerosol
 

  ! by removing those boundary condition, the mixing knows a non-null flux of concentration in the first layer. Will it be considerable? It looks like it's not really important
  ! it's probably already done it the meteo updateforce an open boundary and null flux toward the ground
  ! concentrations(1,:)=concentrations(2,:) ! impose the null flux in the first layer

  !particle_conc(nz,:)=0.0 ! no outflow
  ! particle_conc(1,:)=particle_conc(2,:) ! impose the null flux in the first layer
   
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


  ! aerosol
  OPEN(20, FILE = TRIM(ADJUSTL(outdir))//'/particle.dat'       ,status='replace',action='write')
  OPEN(21, FILE = TRIM(ADJUSTL(outdir))//'/diameters.dat'      ,status='replace',action='write')

  OPEN(45, FILE = TRIM(ADJUSTL(outdir))//'/PN.dat'      ,status='replace',action='write')
  OPEN(46, FILE = TRIM(ADJUSTL(outdir))//'/PM.dat'      ,status='replace',action='write')
  OPEN(47, FILE = TRIM(ADJUSTL(outdir))//'/CS.dat'      ,status='replace',action='write')
END SUBROUTINE open_files


!-----------------------------------------------------------------------------------------
! Save data to files.
!-----------------------------------------------------------------------------------------
SUBROUTINE write_files(time)
  REAL(dp) :: time  ! current time
  CHARACTER(255) :: outfmt0, outfmt1, outfmt2, outfmt3, outfmt4
  ! Get output format for arrays with nz layers
  WRITE(outfmt0, '(a, i3, a)') '(', nz, 'es25.16e3)'
  WRITE(outfmt1, '(a, i3, a)') '(', nz-1, 'es25.16e3)'
  !WRITE(outfmt2, '(a, i3, a)') '(', neq, 'es25.16e3)'
  WRITE(outfmt3, '(a, i3, a)') '(', 1, 'es25.16e3)'
  WRITE(outfmt4, '(a, i3, a)') '(', nr_bins, 'es25.16e3)'

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
  WRITE(20, outfmt4) layer1_particles      ! [m^{-3}...hopfully], particle concentration

  WRITE(45, outfmt0) PN_atm                ! [m^{-3}...hopfully], particle concentration
  WRITE(46, outfmt0) PM_atm                ! [m^{-3}...hopfully], particle concentration
  if (time>time_start_aerosol) then
     WRITE(21, outfmt4) diameter
  end if
  WRITE(47, outfmt0) cond_sink(1,:)        ! [m^{-3}...hopfully], particle concentration
  
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
  CLOSE(20)
  CLOSE(21)
  CLOSE(22)
  CLOSE(23)
  CLOSE(24)
  CLOSE(25)
  CLOSE(26)
  CLOSE(27)
  CLOSE(28)
  CLOSE(29)
  CLOSE(45)
  CLOSE(46)
  CLOSE(47)
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
