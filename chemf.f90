MODULE chemf

  ! TO DO: implement the chemistry mixing for model 1 and 2
  use parameters_mod

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: chemistry_step, neq
  ! Whoever uses this module, can "see" only the subroutine
  ! 'chemistry_step' and variable 'neq'. Other names of all things
  ! are internal to this module.

  ! INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(20,300)
  ! DLSODA is in f77 DOUBLE PRECISION, so here we try match that with this

  ! Module global variables
  INTEGER, PARAMETER :: nreact = 36 ! number of chemical reactions
  INTEGER, PARAMETER :: nspec  = 25 ! number of chemical species

  REAL(dp), DIMENSION(nreact) :: k_rate  ! array of rate coefficients
  REAL(dp) :: Mair_chem, TEMP !, O2, N2, H2O, Emi_iso, Emi_alp

  ! H2O is the same as in parameter_mod, it can be used as a global variable, no need to pass it as an argument
  ! O2 and N2 are already defined in param_mod, no need to declare them in this file (no interference)
  ! Emi_iso and Emi_alp can be defined in param_mod, but their declaration must be deleted from the main
  ! TEMP has no reason to be declared in this scope... but it can be tolerated

  ! stuff needed for DLSODA
  INTEGER, PARAMETER  :: neq   = nspec   ! number of differential equations
  INTEGER, PARAMETER  :: itol  = 1       ! so that atol is a scalar, not array
  INTEGER, PARAMETER  :: itask = 1       ! for normal computation from t to tout
  INTEGER, PARAMETER  :: iopt  = 0       ! for no optional inputs
  INTEGER, PARAMETER  :: lrw   = 22+neq * MAX(16, neq+9) ! real workspace size
  INTEGER, PARAMETER  :: liw   = 20+neq  ! integer workspace size
  INTEGER, PARAMETER  :: mf    = 22      ! stiff, no user-provided jacobian
  REAL(dp), PARAMETER :: rtol = 1d-5     ! relative tolerance
  REAL(dp), PARAMETER :: atol = 1d-3     ! absolute tolerance
  REAL(dp) :: rwork(lrw)   ! real workspace
  INTEGER  :: iwork(liw)   ! integer workspace

  CONTAINS

    SUBROUTINE chemistry_step(Conc,time1,time2,O2_in,N2_in,M_in,Emi_iso_in,Emi_alp_in,CS_h2so4,CS_elvoc,TEMP_in,exp_coszen)
      REAL(dp), INTENT(inout) :: Conc(neq)
      REAL(dp), INTENT(in)    :: time1, time2, O2_in, N2_in, M_in, TEMP_in, Emi_iso_in, Emi_alp_in, CS_h2so4, CS_elvoc !, H2O_in
      REAL(dp), INTENT(in)    :: exp_coszen

      ! for DLSODA:
      INTEGER  :: istate  ! a flag
      REAL(dp) :: time1b
      REAL(dp) :: dummy

      ! We cannot give time1 as input to DLSODA, because DLSODA
      ! modifies the start time variable it is given, but we don't
      ! want to modify here the time1 that the main program gave
      ! us. So let's make a new variable here.

      istate = 1
      time1b = time1

      O2 = O2_in
      N2 = N2_in
      Mair_chem = M_in
      ! H2O = H2O_in
      TEMP =TEMP_in

!      Emi_alp=Emi_alp_in
!      Emi_iso=Emi_iso_in

      CALL calculate_k(exp_coszen,Emi_alp_in,Emi_iso_in,CS_h2so4,CS_elvoc)

      CALL DLSODE (f_lsode, neq, Conc, time1b, time2, itol, rtol, atol, itask, &
                   istate, iopt, rwork, lrw, iwork, liw, dummy, mf)

    END SUBROUTINE chemistry_step


    SUBROUTINE calculate_k(exp_coszen,Emi_alp_in,Emi_iso_in,CS_h2so4,CS_elvoc)

       REAL(dp), INTENT(in) :: exp_coszen, Emi_alp_in, Emi_iso_in, CS_h2so4, CS_elvoc

       k_rate(1)  = 3.83D-5*exp_coszen                                                                                ! O3 = O1D + O2 :
       k_rate(2)  = 1.63D-10*EXP(60/TEMP)                                                                             ! O1D + H2O = OH + OH
       k_rate(3)  = 2.15D-11*EXP(110/TEMP)                                                                            ! O1D + N2 = O3 + REST
       k_rate(4)  = 3.30D-11*EXP(55/TEMP)                                                                             ! O1D + O2 = O3
       k_rate(5)  = 1.67D-2*exp_coszen                                                                                ! NO2 = NO + O3 + REST
       k_rate(6)  = 1.47D-4*exp_coszen                                                                                ! CH2O = HO2 + REST
       k_rate(7)  = 2.40D-13                                                                                          ! OH + CO = HO2 + CO2 + REST
       k_rate(8)  = 2.45D-12*EXP(-1775/TEMP)                                                                          ! OH + CH4 = CH3O2 + REST
       k_rate(9)  = 2.7D-11*EXP(390/TEMP)                                                                             ! isoprene + OH = RO2
       k_rate(10) = 2.40D-11                                                                                          ! OH + MVK = HO2 + CH2O + REST
       k_rate(11) = 3.50D-12*EXP(250/TEMP)                                                                            ! HO2 + NO = OH + NO2
       k_rate(12) = 2.80D-12*EXP(300/TEMP)                                                                            ! CH3O2 + NO = HO2 + NO2 + CH2O + REST
       k_rate(13) = 1.00D-11                                                                                          ! RO2 + NO = HO2 + NO2 + CH2O + MVK
       k_rate(14) = 5.50D-12*EXP(125/TEMP)                                                                            ! OH + CH2O = HO2 + REST
       k_rate(15) = ((2.2D-13*EXP(600/TEMP))+(1.9D-33*EXP(980/TEMP)*Mair_chem))*(1+(1+1.4D-21*EXP(2200/TEMP)*H2O))    ! 2HO2 = H2O2 + O2
       k_rate(16) = 4.10D-13*EXP(750/TEMP)                                                                            ! CH3O2 + HO2 = REST
       k_rate(17) = 1.50D-11                                                                                          ! RO2 + HO2 = REST
       k_rate(18) = 3.50D-12*EXP(340/TEMP)                                                                            ! OH + NO2 = HNO3
       k_rate(19) = 3.00D-12*EXP(-1500/TEMP)                                                                          ! NO + O3 = NO2 + O2
       k_rate(20) = 4.80D-11*EXP(250/TEMP)                                                                            ! OH + HO2 = H2O + O2
       k_rate(21) = 2.90D-12*EXP(-160/TEMP)                                                                           ! OH + H2O2 = H2O + HO2
       k_rate(22) = 1.80D-11*EXP(110/TEMP)                                                                            ! NO + NO3 = NO2 + NO2
       k_rate(23) = 1.40D-13*EXP(-2470/TEMP)                                                                          ! NO2 + O3 = NO3 + O2
       k_rate(24) = (0.35*(3.6D-30*(TEMP/300)**(-4.1)*Mair_chem)*(1.9D-12*(TEMP/300)**0.2)) &
                  / ((3.6D-30*(TEMP/300)**(-4.1)*Mair_chem)+(1.9D-12*(TEMP/300)**0.2))                                ! NO2 + NO3 = N2O5
       k_rate(25) = (0.35*(1.3D-3*(TEMP/300)**(-3.5)*EXP(-11000/TEMP)*Mair_chem)*(9.7D14*(TEMP/300)**0.1*EXP(-11080/TEMP))) &
                  /((1.3D-3*(TEMP/300)**(-3.5)*EXP(-11000/TEMP)*Mair_chem)+(9.7D14*(TEMP/300)**0.1*EXP(-11080/TEMP))) ! N2O5 = NO2 + NO3
       k_rate(26) = 2.50D-22                                                                                          ! N2O5 + H2O = HNO3 + HNO3
       k_rate(27) = 1.80D-39                                                                                          ! N2O5 + H2O + H2O = HNO3 + HNO3 + H2O
       k_rate(28) = 2.03D-16*(TEMP/300)**4.57*EXP(693/TEMP)                                                           ! HO2 + O3 = OH + O2 + O2
       k_rate(29) = 1.5D-12                                                                                           ! SO2 + OH = H2SO4
       k_rate(30) = CS_h2so4 !0.001                                                                                   ! CS for H2SO4 and HNO3
       k_rate(31) = Emi_alp_in                                                                                        ! Emission rate of alpha-pinene
       k_rate(32) = Emi_iso_in                                                                                        ! Emission rate of isoprene
       k_rate(33) = 1.2D-11*EXP(440/TEMP)                                                                             ! alp-pinene + OH
       k_rate(34) = 6.3D-16*EXP(-580/TEMP)                                                                            ! alp-pinene + O3
       k_rate(35) = 1.03D-14*EXP(-1995/TEMP)                                                                          ! isoprene + O3
       k_rate(36) = CS_elvoc !0.001                                                                                   ! CS for ELVOC

    END SUBROUTINE calculate_k


   SUBROUTINE f_lsode(neq, time, Conc, Conc_dot)

     ! This computes the right-hand side of the system Conc' = f(t,Conc), where Conc and f are vectors of length neq.
     ! f cannot have other inputs, since DLSODA assumes that f is called exacly with these input arguments.

     INTEGER,  INTENT(in)  :: neq
     REAL(dp), INTENT(in)  :: time
     REAL(dp), INTENT(in)  :: Conc(neq)
     REAL(dp), INTENT(out) :: Conc_dot(neq)

      ! 1 = O3
      Conc_dot(1)  = 0.

      ! 2 = O1D
      Conc_dot(2)  = k_rate(1)*Conc(1) - k_rate(2)*Conc(2)*H2O - k_rate(3)*Conc(2)*N2 - k_rate(4)*Conc(2)*O2 !OK

      ! 3 = OH
      Conc_dot(3)  = 2.0*k_rate(2)*Conc(2)*H2O-k_rate(7)*Conc(3)*Conc(9)-k_rate(8)*Conc(3)*Conc(11)-k_rate(9)*Conc(3)*Conc(13)&
           -k_rate(10)*Conc(3)*Conc(15) + k_rate(11)*Conc(8)*Conc(6) - k_rate(14)*Conc(3)*Conc(7) &
           - k_rate(18)*Conc(3)*Conc(5)  - k_rate(20)*Conc(3)*Conc(8) &
           - k_rate(21)*Conc(3)*Conc(16) + k_rate(28)*Conc(8)*Conc(1) &
           - k_rate(29)*Conc(3)*Conc(20) - k_rate(33)*Conc(3)*Conc(23)

      ! 4 = REST
      Conc_dot(4)  = 0.

      ! 5 = NO2
      Conc_dot(5)  = 0.

      ! 6 = NO
      Conc_dot(6)  = 0.

      ! 7 = CH2O
      Conc_dot(7)  = -k_rate(6)*Conc(7) &
           + k_rate(10)*Conc(3)*Conc(15) + k_rate(12)*Conc(12)*Conc(6) + k_rate(13)*Conc(14)*Conc(6) - k_rate(14)*Conc(3)*Conc(7)

      ! 8 = HO2
      Conc_dot(8)  = k_rate(6)*Conc(7) + k_rate(7)*Conc(3)*Conc(9) &
           + k_rate(10)*Conc(3)*Conc(15) - k_rate(11)*Conc(8)*Conc(6) + k_rate(12)*Conc(12)*Conc(6) + k_rate(13)*Conc(14)*Conc(6)&
           + k_rate(14)*Conc(3)*Conc(7) - 2.0*k_rate(15)*Conc(8)**2   - k_rate(16)*Conc(12)*Conc(8) - k_rate(17)*Conc(14)*Conc(8)&
           - k_rate(20)*Conc(3)*Conc(8) &
           + k_rate(21)*Conc(3)*Conc(16) - k_rate(28)*Conc(8)*Conc(1)

      ! 9 = CO
      Conc_dot(9)  = 0.

      ! 10 = CO2
      Conc_dot(10) = 0.

      ! 11 = CH4
      Conc_dot(11) = 0.

      ! 12 = CH3O2
      Conc_dot(12) = -(k_rate(12)*Conc(6) + k_rate(16)*Conc(8))*Conc(12)&
           + k_rate(8)*Conc(3)*Conc(11)

      ! 13 = Isoprene
      Conc_dot(13) = -k_rate(9)*Conc(3)*Conc(13) - k_rate(35)*Conc(13)*Conc(1)&
           + k_rate(32) ! must be constant in BM

      ! 14 = RO2
      Conc_dot(14) = -(k_rate(13)*Conc(6) + k_rate(17)*Conc(8))*Conc(14)&
           + k_rate(9)*Conc(3)*Conc(13)

      ! 15 = MVK
      Conc_dot(15) = - k_rate(10)*Conc(3)*Conc(15) + k_rate(13)*Conc(14)*Conc(6)

      ! 16 = H2O2
      Conc_dot(16) = - k_rate(21)*Conc(3)*Conc(16) + k_rate(15)*Conc(8)**2

      ! 17 = HNO3
      Conc_dot(17) = - k_rate(30)*Conc(17)&
           + k_rate(18)*Conc(3)*Conc(5) + 2.0*k_rate(26)*Conc(19)*H2O + 2.0*k_rate(27)*Conc(19)*H2O**2

      ! 18 = NO3
      Conc_dot(18) = -(k_rate(22)*Conc(6)+k_rate(24)*Conc(5))*Conc(18)&
           +k_rate(23)*Conc(5)*Conc(1)+k_rate(25)*Conc(19) !OK

      ! 19 = N2O5
      Conc_dot(19) = -(k_rate(25)+k_rate(26)*H2O+k_rate(27)*H2O**2)*Conc(19)&
           +k_rate(24)*Conc(5)*Conc(18)  ! OK

      ! 20 = SO2
      Conc_dot(20) = 0.

      ! 21 = H2SO4
      Conc_dot(21) = -k_rate(30)*Conc(21) + k_rate(29)*Conc(3)*Conc(20) ! CHECK if it's OK

      !22 = H2SO4_P
      Conc_dot(22) = k_rate(30)*Conc(21)

      ! 23 = Alpha-pinene
      Conc_dot(23) = k_rate(31) - k_rate(33)*Conc(3)*Conc(23) - k_rate(34)*Conc(1)*Conc(23) !must be constant in BM

      !24 = HNO3_P
      Conc_dot(24) = k_rate(30)*Conc(17)

      !25 = ELVOC
      Conc_dot(25) = 0.05 * k_rate(33) * Conc(3) * Conc(23) + &
                       0.1 * k_rate(34) * Conc(1) * Conc(23) - &
                             k_rate(36) * Conc(25)

    END SUBROUTINE f_lsode

END MODULE chemf
