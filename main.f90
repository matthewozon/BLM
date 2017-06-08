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

USE meteorology_mod

!-----------------------------------------------------------------------------------------
! Variable declaration
!-----------------------------------------------------------------------------------------
IMPLICIT NONE

! Used for loops
INTEGER :: I, J

! model selection
integer :: model=3

!-----------------------------------------------------------------------------------------
! Initialization
!-----------------------------------------------------------------------------------------
CALL time_init()         ! initialize time, time step, date

CALL meteorology_init(model)  ! initialize ua, va and theta

CALL open_files()        ! open files for future use

CALL write_files(time)   ! write initial values


!-----------------------------------------------------------------------------------------
! Start main loop
!-----------------------------------------------------------------------------------------
DO WHILE (time <= time_end)
   !---------------------------------------------------------------------------------------
   ! Meteorology
   !---------------------------------------------------------------------------------------
   ! Set lower boundary condition
   CALL surface_values(theta(1), time)
   ! update wind velocity and potential temperature
   call update_meteo(model)
   


  !---------------------------------------------------------------------------------------
  ! Chemistry
  !---------------------------------------------------------------------------------------
  ! Start to calculate chemical reactions only after some time to save the computation
  ! time
  IF (time >= time_start_chemistry) THEN
    ! Calculate chemistry every dt_chem
    IF ( MOD( NINT((time - time_start_chemistry)*1000.0_dp), NINT(dt_chem*1000.0_dp) ) == 0 ) THEN
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
END SUBROUTINE open_files


!-----------------------------------------------------------------------------------------
! Save data to files.
!-----------------------------------------------------------------------------------------
SUBROUTINE write_files(time)
  REAL(dp) :: time  ! current time
  CHARACTER(255) :: outfmt0, outfmt1
  ! Get output format for arrays with nz layers
  WRITE(outfmt0, '(a, i3, a)') '(', nz, 'es25.16e3)'
  WRITE(outfmt1, '(a, i3, a)') '(', nz-1, 'es25.16e3)'

  ! Only save h one time at the beginning
  IF (time == time_start) THEN
    WRITE(12, *) h
  END IF

  WRITE(11, '(f8.4)') time/(24*one_hour)  ! [day], time
  WRITE(13, outfmt0) ua                  ! [m s-1], u wind
  WRITE(14, outfmt0) va                  ! [m s-1], v wind
  WRITE(15, outfmt0) theta               ! [K], potential temperature
  WRITE(16, outfmt1) k_closure           ! [m^2 s^{-1}], k_closure K
  WRITE(17, outfmt1) k_closure_m         ! [m^2 s^{-1}], k_closure K_m
  WRITE(18, outfmt1) k_closure_h         ! [m^2 s^{-1}], k_closure K_h
  WRITE(19, outfmt1) Ri_num              ! [], Richardson number Ri
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
END SUBROUTINE close_files

END PROGRAM main
