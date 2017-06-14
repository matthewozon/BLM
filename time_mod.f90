!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Time module
!
! - Time variables
!
! - Time functions
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

MODULE Time_Mod

!-----------------------------------------------------------------------------------------
! Load modules
!-----------------------------------------------------------------------------------------
USE Parameters_Mod

!-----------------------------------------------------------------------------------------
! Variable declaration
!-----------------------------------------------------------------------------------------
IMPLICIT NONE

PUBLIC

INTEGER, PARAMETER :: one_hour = 60*60  ! [s], one hour in seconds

REAL(dp) :: time  ! [s], current time
REAL(dp) :: time_start, time_end  ! [s], start and end time

REAL(dp) :: dt  ! [s], time step for main loop, usually is equal to meteorology time step
REAL(dp) :: dt_chem  ! [s], time step for chemistry calculation
REAL(dp) :: dt_aero  ! [s], time step for aerosol processes
REAL(dp) :: dt_output  ! [s], time step for output

REAL(dp) :: time_start_chemistry  ! [s], time to start calculating chemistry
REAL(dp) :: time_start_aerosol ! [s], time to start calculating aerosol

INTEGER :: daynumber_start  ! [day], start Julian day
INTEGER :: daynumber  ! [day], current Julian day

INTEGER :: counter  ! [-], counter of time steps


CONTAINS


!-----------------------------------------------------------------------------------------
! Set initial values for time_start, time_end, time steps, day number, counter
!-----------------------------------------------------------------------------------------
SUBROUTINE Time_Init()
  ! Simulation time period
  time_start = 0.0d0
  time_end = 5.0d0 * 24.0d0 * one_hour

  ! Time steps
  dt = 0.5d0
  dt_chem = 10.0d0
  dt_aero = 10.0d0
  dt_output = 3600.0d0

  ! Get the Julian date of Aug. 10, 2011
  daynumber_start = 31+28+31+30+31+30+31+10

  ! Start to run chemistry module after 3 days to save computation time
  time_start_chemistry = 3*24*one_hour

  ! Start to run aerosol evolution after 3 days
  time_start_aerosol = 3*24*one_hour

  ! Current time and date
  time = time_start
  daynumber = daynumber_start

  ! Count how many time steps the code runs
  counter = 0
END SUBROUTINE Time_Init

END MODULE Time_Mod
