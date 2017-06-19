MODULE aerosol_mod

  use parameters_mod
  IMPLICIT NONE

  PUBLIC

  !! ====================== Definition of variables =====================================================================
  
    ! so that numbers will be in 64bit floating point
  ! http://en.wikipedia.org/wiki/Double_precision_floating-point_format

  INTEGER, PARAMETER ::  nr_bins = 100           ! Number of particle size bins
  INTEGER, PARAMETER ::  nr_cond = 2             ! Number of condensable vapours
    
  REAL(DP), DIMENSION(nr_bins) :: diameter  ,&    ! Diameter of each size bin
       particle_mass                        ,&    ! mass of one particle in each size bin
       particle_volume                      ,&    ! volume concentration in each size bin
       coag_loss                            ,&    ! coagulation loss rate of particles in each size bin
       v_dep                                      ! Dry deposition velocity of particles
  real(dp) :: last_dia, last_vol                  ! extra bin for losses due to out of range
  
  
       
  REAL(DP), DIMENSION(nr_cond) :: molecular_mass   ,&   ! molecular mass of the condensing vapours [kg/#]
       molecular_volume                            ,&   ! Molecule volume of condensable vapours [m^3]
       molecular_dia                               ,&   ! Molecule diameter of condensable vapours [m]
       molar_mass                                  ,&   ! Molar mass of condensable vapours [kg/m^3]
       cond_vapour                                    ! Concentration of condensable vapours [molec/m^3]
  

  REAL(DP) ::         PN, PM                  ! Total particle number and mass concentration [cm^-3] 
       
  REAL(DP) :: vd_SO2,vd_O3, vd_HNO3           ! Dry deposition velocity of SO2, O3 & HNO3     

  REAL(DP)    particle_density              ,&    ! [Kg]
              nucleation_coef               ,&    ! Nucleation coefficient 
              mass_accomm                         ! mass accomodation coefficient 

  REAL(dp) :: nucleation_rate                     ! #/(m3*s)

  REAL(DP) :: simu_hours                    ,&    ! total simulation time in hours
              timestep                      ,&    ! Model time step [s]
              temperature                   ,&    ! Temperature [K]
              Richards_nr10m                ,&    ! Richards number at the reference altitude 10 m
              wind_speed10m                 ,&    ! Wind speed at the reference altitude 10 m
              pressure                      ,&    ! Pressure [Pa]
              DSWF                          ,&    ! Downward Shortwave Radiation Flux (W/m^2)
              Mixing_height                       ! Boundary layer mixing height [m]
 
  !! ======================= Programe starts ===========================================================================

CONTAINS

SUBROUTINE Aerosol_init(diameter, particle_mass, particle_volume, particle_conc, &
  particle_density, nucleation_coef, molecular_mass, molar_mass, &
  molecular_volume, molecular_dia, mass_accomm,cond_sink)

    !! ====================== Definition of variables =====================================================================

    REAL(DP), DIMENSION(nr_bins), INTENT(OUT) :: diameter       , &    ! diamter of each size bin
                                                 particle_mass  , &    ! mass of one particle
                                                 particle_volume  , &  ! volume of one particle
                                                 particle_conc         ! number concentration

    REAL(DP), DIMENSION(nr_cond), INTENT(OUT) :: molecular_mass ,& ! molecular mass of the condensing vapours [kg/#]
                                           molecular_volume ,&     ! [m3]
                                           molecular_dia, &        ! [m]
                                           molar_mass, &           ! molar mass of the condensing vapours [kg/mol]
                                           cond_sink               ! condensation sink of all the condensable

    REAL(DP), INTENT(OUT) :: nucleation_coef, mass_accomm

    REAL(DP), DIMENSION(nr_cond) :: density             ! Bulk density of condensing vapours [kg/m^3]
 
    REAL(DP)  ::  particle_density
    
    INTEGER :: i
    mass_accomm = 1D0   ! Mass accommodation coefficient 
    
    nucleation_coef = 1D-20
    
    ! Particle diameters between 2D-9 and 2.5D-6 m:
    diameter(1)=2D-9 
    DO i=2,nr_bins
       diameter(i)=diameter(i-1)*(2.5D-6/diameter(1))**(1D0/(nr_bins-1))
    END DO

    ! compute the diameter and the volume of the last+1 theoretcal bin
    last_dia = diameter(nr_bins)*(2.5D-6/diameter(1))**(1D0/(nr_bins-1))
    last_vol = 1D0/6D0 * pi * last_dia**3
      
    call init_part_conc(particle_conc,diameter)
    
    particle_density = 1.4D3                                        ! Assumed fixed particle density [kg/m^3]
    particle_volume = 1D0/6D0 * pi * diameter**3                      ! Single particle volume (m^3)
    particle_mass=  1D0/6D0 * pi * diameter**3 * particle_density     ! [kg]

    density = (/1.84D3, 1.4D3/)                                     ! density of sulphuric acid and SOA
    molar_mass = (/0.098D0, 0.3D0/)                                 ! H2SO4 and ELVOC
    molecular_mass = molar_mass / NA                                ! molecular mass [kg]
    molecular_volume = molecular_mass / density                     ! molecular volume [m^3]
    molecular_dia = (6D0 * molecular_volume / pi )**(1D0/3D0)       ! molecular diameter [m]

    cond_sink = 0.001

    call PNandPM(particle_conc)

  END SUBROUTINE Aerosol_init

  SUBROUTINE PNandPM(particle_conc_NM)
    real(dp),dimension(nr_bins)::particle_conc_NM
    PN=sum(particle_conc_NM)
    PM=sum(particle_conc_NM*particle_mass)
  end SUBROUTINE PNandPM

  subroutine init_part_conc(particle_conc, diameter)
    REAL(DP), DIMENSION(nr_bins), INTENT(IN) :: diameter               ! diamter of each size bin
    REAL(DP), DIMENSION(nr_bins), INTENT(OUT) :: particle_conc         ! number concentration

    ! Assume an initial particle number concentration of 1 m^-3
    particle_conc = 1D0
    ! add 200 cm^-3 200 nm sized accumulation mode particles
    where((abs(diameter-2D-7)-MINVAL(abs(diameter-2D-7)))<1D-20)  particle_conc=2D8
  end subroutine init_part_conc
  
  
  
SUBROUTINE Nucleation(timestep,nucleation_coef,h2so4_conc,particle_conc_nuc) ! (Add input and output variables here)
  real(dp), intent(in)::timestep, nucleation_coef, h2so4_conc
  real(dp), dimension(nr_bins)::particle_conc_nuc
  ! Consider how kinetic H2SO4 nucleation influence the number concentrations of particles 
  ! in the fist size bin particle_conc(1) within one model time step
  particle_conc_nuc(1)=particle_conc_nuc(1)+timestep*nucleation_coef*h2so4_conc**2
END SUBROUTINE Nucleation


SUBROUTINE Condensation(timestep, temperature, pressure, mass_accomm, molecular_mass, &
  molecular_volume, molar_mass, molecular_dia, particle_mass, particle_volume, &
  particle_conc_cond, diameter, cond_vapour, cond_sink) ! Add more variables if you need it
  
    REAL(DP), DIMENSION(nr_bins), INTENT(IN) :: diameter, particle_mass
    REAL(DP), DIMENSION(nr_cond), INTENT(IN) :: molecular_mass, molecular_dia, &
    molecular_volume, molar_mass
    REAL(DP), INTENT(IN) :: timestep, temperature, pressure, mass_accomm

    REAL(DP), DIMENSION(nr_bins), INTENT(INOUT) :: particle_conc_cond

    REAL(DP), DIMENSION(2), INTENT(IN) :: cond_vapour  ! condensing vapour concentrations, which is H2SO4 and organics (ELVOC) [#/m^3]
    
    REAL(DP), DIMENSION(nr_bins), INTENT(IN)   :: particle_volume
    
    REAL(DP), DIMENSION(nr_bins)   ::  slip_correction, diffusivity, speed_p, &
    particle_conc_new 
    
    REAL(DP), DIMENSION(nr_cond)   ::  diffusivity_gas, speed_gas 
    
    REAL(DP) :: dyn_visc, l_gas, dens_air

    ! real(dp), intent(in):: M_air

    real(dp), DIMENSION(nr_cond), INTENT(out) :: cond_sink
    
    INTEGER :: j

    real(DP), dimension(nr_cond,nr_bins)::kn_number, beta_fs ! Knudsen number and Fuchs-Sutugin correction coefficient
    real(DP), dimension(nr_cond,nr_bins)::CR  ! collision rate
    real(DP), dimension(nr_bins)::particle_volume_new
    real(dp)::xj ! bin fraction
    
    ! Add more variabels as you need it...

    dyn_visc = 1.8D-5*(temperature/298D0)**0.85D0  ! dynamic viscosity of air
    dens_air=Mair*pressure/(R*temperature)        ! Air density
    l_gas=2D0*dyn_visc/(pressure*SQRT(8D0*Mair/(pi*R*temperature))) ! Gas mean free path in air (m)

    slip_correction = 1D0+(2D0*l_gas/(diameter))*&
    (1.257D0+0.4D0*exp(-1.1D0/(2D0*l_gas/diameter))) ! Cunninghams slip correction factor (Seinfeld and Pandis eq 9.34) 
    
    diffusivity = slip_correction*kb*temperature/(3D0*pi*dyn_visc*diameter)   ! Diffusivity for the different particle sizes m^2/s
    speed_p = SQRT(8D0*kb*temperature/(pi*particle_mass))                     ! speed of particles (m/s)
  
    diffusivity_gas=5D0/(16D0*NA*molecular_dia**2D0*dens_air)*&
    SQRT(R*temperature*Mair/(2D0*pi)*(molar_mass+Mair)/molar_mass)            ! Diffusivity of condensable vapours (m^2 s^-1)

    ! Thermal velocity of vapour molecule
    speed_gas=SQRT(8D0*kb*temperature/(pi*molecular_mass)) ! speed of H2SO4 molecule

    ! Calculate the Fuchs-Sutugin correction factor:
    do j= 1,nr_cond
       kn_number(j,:) = 2.0_dp*(3.0*(diffusivity_gas(j)+diffusivity)/sqrt(speed_gas(j)**2+speed_p**2))/(diameter+molecular_dia(j))
       beta_fs(j,:) = 0.75_dp*mass_accomm*(1.0_dp+kn_number(j,:))/(kn_number(j,:)**2+0.283_dp*kn_number(j,:)*mass_accomm+kn_number(j,:)+0.75_dp*mass_accomm)
    end do
    
    ! Calculate the Collision rate (CR [m^3/2]) between gas molecules (H2SO4 and ELVOC) and the particles:
    do j= 1,nr_cond
       CR(j,:) = 2.0_dp*pi*(diameter+molecular_dia(j))*(diffusivity+diffusivity_gas(j))*beta_fs(j,:)
    end do
    
    ! Calculate the new single particle volume after condensation (particle_volume_new):
    particle_volume_new=particle_volume
    do j= 1,nr_cond
       particle_volume_new=particle_volume_new+timestep*CR(j,:)*molecular_volume(j)*cond_vapour(j)
    end do
    
    ! Use the full-stationary method to divide the particles between the existing size bins (fixed diameter grid):
    particle_conc_new=0D0 ! Initialise a new vector with the new particle concentrations

    DO j = 1,nr_bins-1
       ! Add equations that redistributes the particle number concentration 
       !in size bin 1 to nr_bins-1 to the fixed volume (diameter) grid
       xj = (particle_volume(j+1)-particle_volume_new(j))/(particle_volume(j+1)-particle_volume(j))

       ! avoid negative or too rapid growth. Sometimes negative values appear when the particle sizes are small and it becomes unstable.
       if (xj<0.0) then
          xj=0.0
       end if
       if (xj>1.0) then
          xj=1.0
       end if

       ! distribute the new formed particle into the fix bins
       particle_conc_new(j)=particle_conc_new(j)+xj*particle_conc_cond(j)
       particle_conc_new(j+1)=particle_conc_new(j+1)+(1.0-xj)*particle_conc_cond(j)
    END DO

    ! for the last bin, maybe consider an out-of-range loss (should not account for much if the volume range include large particle)
    xj = (last_vol-particle_volume_new(nr_bins))/(last_vol-particle_volume(nr_bins))
    ! avoid negative or too rapid growth. Sometimes negative values appear when the particle sizes are small and it becomes unstable.
    if (xj<0.0) then
       xj=0.0
    end if
    if (xj>1.0) then
       xj=1.0
    end if
    ! compute the loss by condensation (no volume conservation)
    particle_conc_new(nr_bins)=particle_conc_new(nr_bins)+xj*particle_conc_cond(nr_bins)
    
    ! Update the particle concentration in the particle_conc vector:
    particle_conc_cond=particle_conc_new

    ! compute the condensation sink here
    do j = 1,nr_cond
       cond_sink(j) = sum(particle_conc_cond*CR(j,:)) ! particle [nr_bins], CR [nr_cond nr_bins]
    end do
    
END SUBROUTINE Condensation

SUBROUTINE Coagulation(timestep, particle_conc_coag, diameter, &
  temperature,pressure,particle_mass) !,M_air) ! Add more variables if you need it
  
    REAL(DP), DIMENSION(nr_bins), INTENT(IN) :: diameter
    REAL(DP), DIMENSION(nr_bins), INTENT(INOUT) :: particle_conc_coag
    REAL(DP), INTENT(IN) :: timestep
    REAL(DP), DIMENSION(nr_bins), INTENT(IN) :: particle_mass       ! mass of one particle                                 
    REAL(DP), INTENT(IN) :: temperature, pressure
    
    REAL(DP), DIMENSION(nr_bins,nr_bins) :: coagulation_coef        ! coagulation coefficients [m^3/s]

    REAL(DP), DIMENSION(nr_bins) :: slip_correction, diffusivity, dist, speed_p, &
    Beta_Fuchs, free_path_p

    REAL(DP) ::       dyn_visc, &                                   ! dynamic viscosity, kg/(m*s)
                      l_gas                                         ! Gas mean free path in air
    INTEGER  :: i, j

    REAL(DP), DIMENSION(nr_bins) :: particle_conc_new

    ! real(dp), intent(in)::M_air

    ! The Coagulation coefficient is calculated according to formula 13.56 in Seinfield and Pandis (2006), Page 603
    
    dyn_visc = 1.8D-5*(temperature/298.)**0.85                                              ! Dynamic viscosity of air

    l_gas=2D0*dyn_visc/(pressure*SQRT(8D0*Mair/(pi*R*temperature)))                        ! Gas mean free path in air (m)

    slip_correction = 1D0+(2D0*l_gas/(diameter))*&
         (1.257D0+0.4D0*exp(-1.1D0/(2D0*l_gas/diameter)))                                        ! Cunninghams slip correction factor (Seinfeld and Pandis eq 9.34)

    diffusivity = slip_correction*kb*temperature/(3D0*pi*dyn_visc*diameter)                 ! Diffusivity for the different particle sizes m^2/s

    speed_p = SQRT(8D0*kb*temperature/(pi*particle_mass))                                   ! Speed of particles (m/s)

    free_path_p = 8D0*diffusivity/(pi*speed_p)                                              ! Particle mean free path (m)

    dist = (1D0/(3D0*diameter*free_path_p))*((diameter+free_path_p)**3D0 &
         -(diameter**2D0+free_path_p**2D0)**(3D0/2D0))-diameter                    ! mean distance from the center of a sphere reached by particles leaving the sphere's surface (m)

    DO i = 1,nr_bins
       Beta_Fuchs = 1D0/((diameter+diameter(i))/(diameter+diameter(i)+&
            2D0*(dist**2D0+dist(i)**2D0)**0.5D0)+8D0*(diffusivity+diffusivity(i))/&
            (((speed_p**2D0+speed_p(i)**2D0)**0.5D0)*(diameter+diameter(i))))                    ! Fuchs correction factor from Seinfeld and Pandis, 2006, p. 600

       coagulation_coef(i,:) = 2D0*pi*Beta_Fuchs*(diameter*diffusivity(i)+&
            diameter*diffusivity+diameter(i)*diffusivity+diameter(i)*diffusivity(i))             ! coagulation rates between two particles of all size combinations  (m^3/s)    
    END DO

    ! Write equations that considers how the particle number concentration in each size bin 
    !(particle_conc) is influenced by the coagulation sink (loss of smaller particles when
    ! they collide with larger ones)

    ! You can first calculate the loss (loss1) due to self-coagulation between particles in the same size bin
    ! and then calculate the loss (loss2) due to coagulation with larger particles
    ! Then add the two loss terms together loss = loss1 + loss2

    particle_conc_new = particle_conc_coag
    do i = 1,nr_bins
       particle_conc_new(i) = particle_conc_new(i)&
            - timestep*particle_conc_coag(i)*sum(coagulation_coef(i,(i+1):nr_bins)*particle_conc_coag((i+1):nr_bins))
    end do
    particle_conc_coag=particle_conc_new
    
    !particle_conc_new=particle_conc_coag
    !do i = 1,nr_bins
    !   do j = (i+1),nr_bins
    !      particle_conc_new(i) = particle_conc_new(i)&
    !           - timestep*particle_conc_coag(i)*coagulation_coef(i,j)*particle_conc_coag(j)
    !   end do
    !end do
    !particle_conc_coag=particle_conc_new
END SUBROUTINE Coagulation
  
SUBROUTINE dry_dep_velocity(diameter,particle_density,temperature,pressure,DSWF, & 
   Richards_nr10m,wind_speed10m) !,M_air) ! Add more variables if you need it
   
      REAL(dp), DIMENSION(nr_bins), INTENT(IN) :: diameter
      
      REAL(dp), INTENT(IN) :: temperature, pressure, Richards_nr10m, DSWF, &
      wind_speed10m, particle_density
            
      REAL(dp) :: z0m, r_coll, a_landuse, j_landuse, v_kinematic,dyn_visc,l_gas,Pr,beta,&
      gam,zr,u_friction,dens_air, L_Ob, raO3, raSO2, raHNO3, raisoprene, raapinene
      
      ! Specific parameters for the surface resistance of gases:
      REAL(dp) :: rj,rlu,rac, &
       DiffusivityH2O, D_ratio_SO2, D_ratio_O3, D_ratio_HNO3, D_ratio_isoprene, D_ratio_apinene, &
       DiffusivitySO2, DiffusivityO3, DiffusivityHNO3, Diffusivityisoprene, Diffusivityapinene,&
       z_roughSO2, z_roughO3, z_roughHNO3, z_roughisoprene, z_roughapinene, &
       ScSO2, ScO3, ScHNO3, Scisoprene, Scapinene, &
       rbSO2, rbO3, rbHNO3, rbisoprene, rbapinene, &
       H_effSO2, H_effO3, H_effHNO3, H_effisoprene, H_effapinene, &
       f0_SO2, f0_O3, f0_HNO3, f0_isoprene, f0_apinene, &
       rclSO2, rclO3, rgsSO2, rgsO3

      ! real(dp), intent(in)::M_air
            
       dens_air = Mair*pressure/(R*temperature)    ! Air density (kg/m^3)
       dyn_visc = 1.8D-5*(temperature/298.)**0.85   ! dynamic viscosity of air (kg/(m*s))
       v_kinematic = dyn_visc/dens_air              ! kinematic viscosity of air (m^2/s)
       
       zr=10D0                                      ! Reference height [m]
       L_Ob=zr/Richards_nr10m                       ! Monin-Obukhov length scale
       z0m = 0.9D0 ! Surface roughness length for momentum evergreen, needleleaf trees (m)     
       u_friction=vonk*wind_speed10m/(log(zr/z0m))     ! Friction velocity (Eq. 16.67 from Seinfeld and Pandis, 2006)
 
       ! Land use category paramaters from Seinfeld and Pandis, 2006 Table 19.2: 
       r_coll = 2D-3 ! radius of collector evergreen, needleleaf trees
       
       ! coefficients based on land use categories (evergreen, needleleaf trees)
       a_landuse = 1D0
       j_landuse = 0.56D0
    
       Pr = 0.95D0        ! Turbulent Prandtl number (when ka = 0.4 (Hogstrom, 1988))
       beta = 7.8D0       ! When ka = 0.4 (Hogstrom, 1988)
       gam = 11.6D0       ! When ka = 0.4 (Hogstrom, 1988)

       ! Calculate the particle sedimentation velocity:
               
       ! Calculation of aerodynamic resistance for particles for:
       ! stable boundary layer (Ri>1D-6)
       ! neutral boundary layer (abs(Ri)<1D-6
       ! unstable boundary layer Ri<-1D-6
       
        
       ! Calculate the quasi-laminar resistance (rb) for particles:
       
       ! Calculate the dry deposition velocity for particles:
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Calculate the dry deposition velocity for O3, SO2, HNO3, isoprene and a-pinene: 
       
       ! Resistance components used when calculating the surface resistance for gases, 
       ! table 19.3 Seinfeld and Pandis, 2006: 
      
      ! The minimum, bulk canopy stomatal resistance for water vapor:
      rj = 130D0 ! (s/m) Summer, evergreen, needleleaf
      
      ! The resistance of the outer surfaces in the upper canopy
      rlu = 2000D0 ! (s/m) Summer, evergreen, needleleaf
      
      ! transfer resistance on the ground (that depends only on canopy height)
      rac = 2000D0 ! (s/m) Summer, evergreen, needleleaf
      
      ! resistance for uptake by soil, leaf litter, and so on at the ground SO2
      rgsSO2 = 500D0 ! (s/m) Summer, evergreen, needleleaf
      
      ! restistance for uptake by soil, leaf litter, and so on at the ground, O3
      rgsO3 = 200D0 ! (s/m) Summer, evergreen, needleleaf
      
      ! resistance for uptake by leaves,twigs, and other exposed surfaces, SO2
      rclSO2 = 2000D0 ! (s/m) Summer, evergreen, needleleaf  
      
      ! resistance for uptake by leaves, twigs, and other exposed surfaces, O3
      rclO3 = 1000D0 ! (s/m) Summer, evergreen, needleleaf    
        
      ! Diffusion coefficients of selected gases
      DiffusivityH2O = 0.234D-4 ! Diffusion coefficient of water vapor in air (m^2/s), table 16.2 Seinfeld and Pandis
      
      ! ratio between diffusivity of water vapor and SO2, O3 or HNO3 from table 19.4 Seinfeld and Pandis
      D_ratio_SO2 = 1.89D0
      D_ratio_O3 = 1.63D0
      D_ratio_HNO3 = 1.87D0
      D_ratio_isoprene = 2.7D0 ! Estimated 
      D_ratio_apinene = 4D0  ! Estimated 
      
      
      DiffusivitySO2 = DiffusivityH2O/D_ratio_SO2    ! Diffusivity of SO2 (m^2/s)
      DiffusivityO3 = DiffusivityH2O/D_ratio_O3      ! Diffusivity of O3 (m^2/s)
      DiffusivityHNO3 = DiffusivityH2O/D_ratio_HNO3  ! Diffusivity of HNO3 (m^2/s)
      Diffusivityisoprene = DiffusivityH2O/D_ratio_isoprene  ! Diffusivity of isoprene (m^2/s)
      Diffusivityapinene = DiffusivityH2O/D_ratio_apinene  ! Diffusivity of apinene (m^2/s)
      
      ! Calculate the aerodynamic resistance for O3, SO2, HNO3, isoprene & a-pinene (ra) in similar way as
      ! for particles: 
       		
      ! Calculate the quasi-laminar resistance for O3, SO2, HNO3, isoprene & a-pinene (rb):      
      
      ! Calculation of surface resistance for O3, SO2, HNO3, isoprene & a-pinene (rc):                         
      
       ! Effective Henry's lay const:
       H_effSO2 = 1D5   ! M atm^-1
       H_effO3 = 1D-2   ! M atm^-1
       H_effHNO3 = 1D14 ! M atm^-1
       H_effisoprene = 1.2D-2 ! M atm^-1
       H_effapinene = 3D-2 ! M atm^-1
       
       ! Noramlized reactivity, table 19.4 from Seinfeld and Pandis, 2006:
       f0_SO2 = 0D0
       f0_O3 = 1D0
       f0_HNO3 = 0D0
       f0_isoprene = 0D0
       f0_apinene = 0D0
      
      ! Calculate the bulk canopy stomatal resistance (rst)
      
      ! Calculate the combined stomatal and mesophyll resistance (rsm):
      
      ! Calculate the resistance of the outer surfaces in the upper canopy (rlu):
      
      ! Calculate the resistance to transfer by buoyant convection (rdc):
      
      ! Calculate the resistance of the exposed surfaces in the lower portions of 
      ! structures of the canopy (rcl): 
              
      ! Calculate the resistance of the exposed surfaces on the groud 
      !(soil,leaf litter, ground) (rgs):
      
      ! Combine all resistances in order to get the total surface resistance 
      ! for O3, SO2, HNO3, isoprene and a-pinene (rc):
       
      ! Finally calculate the dry deposition velocity of SO2, O3, HNO3, isoprene and a-pinene:                                        

END SUBROUTINE dry_dep_velocity



!-----------------------------------------------------------------------------------------
! Save data to files.
!-----------------------------------------------------------------------------------------
!SUBROUTINE write_files(iter)
!  CHARACTER(255) :: outfmt0, outfmt1
!  integer::iter
!  ! Get output format for arrays
!  WRITE(outfmt0, '(a, i3, a)') '(', nr_bins, 'es25.16e3)'
!  WRITE(outfmt1, '(a, i3, a)') '(', nr_cond, 'es25.16e3)'
!
!  ! Only save h one time at the beginning
!  IF (iter == 1) THEN
!    WRITE(11, outfmt0) particle_volume
!  END IF
!  WRITE(12, outfmt0) particle_conc         ! [# m^{-3}], particle concentration
!  
!  
!END SUBROUTINE write_files



END MODULE aerosol_mod


