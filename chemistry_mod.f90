!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Chemistry module
!
! - Emissions
!
! - Chemistry
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

MODULE chemistry_mod

USE parameters_mod
USE time_mod
USE grid_mod

IMPLICIT NONE

PRIVATE

! Emission variables

! Chemistry variables

! Stuff needed for DLSODE, do not put them into the subroutine, iwork will cause problems
INTEGER, PARAMETER  :: itol  = 1       ! so that atol is a scalar, not array
INTEGER, PARAMETER  :: itask = 1       ! for normal computation from t to tout
INTEGER, PARAMETER  :: iopt  = 1       ! for no optional inputs
INTEGER, PARAMETER  :: lrw   = 22+neq * MAX(16, neq+9) ! real workspace size
INTEGER, PARAMETER  :: liw   = 20+neq  ! integer workspace size
INTEGER, PARAMETER  :: mf    = 22      ! stiff, no user-provided jacobian
REAL(dp), PARAMETER :: rtol = 1d-5     ! relative tolerance
REAL(dp), PARAMETER :: atol = 1d-3     ! absolute tolerance

REAL(dp) :: rwork(lrw)   ! real workspace
INTEGER  :: iwork(liw)   ! integer workspace


CONTAINS


SUBROUTINE chemistry_init()
END SUBROUTINE chemistry_init


END MODULE chemistry_mod
