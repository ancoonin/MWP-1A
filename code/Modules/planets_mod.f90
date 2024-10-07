module planets_mod
!______________________________________________________________________________________________________________________!
use constants_mod

   real :: radius
   real :: mass
   real :: rhoi
   real :: rhow
   real :: gacc
   real :: omega
   real :: acoef, ccoef
   real :: moiA, moiC
   real :: kf

   contains

   subroutine earth_init

      radius = 6.371E6              ! Radius of the Earth (m)
      mass = 5.976E24               ! Mass of the Earth (kg)
      rhoi = 920.0                  ! Density of ice (kg/m^3)
      rhow = 1000.0                 ! Density of fresh water (kg/m^3)
      gacc = 9.80665                ! Acceleration due to gravity at the Earth's surface (m/s^2)
      omega = 7.292e-5              ! Rotation rate of the Earth (rad/s)
      moiA=0.3296145*mass*radius**2 ! Principal moment of inertia of the Earth
      moiC=0.3307007*mass*radius**2 ! Principal moment of inertia of the Earth
      kf = 0.9342+0.008             ! Fluid (Tidal) Love number

   end subroutine earth_init

   subroutine mars_init

      radius = 3.3899E6             ! Radius of Mars (m)
      mass = 6.4185E23              ! Mass of Mars (kg)
      rhoi = 1220.0                 ! Density of ice mix on Mars (kg/m^3)
      rhow = 1000.0                 ! Density of fresh water (kg/m^3)
      gacc = 3.713                  ! Acceleration due to gravity at Mars's surface (m/s^2)
      omega = 7.08819118E-5         ! Rotation rate of Mars (rad/s)
      moiA=0.363914*mass*radius**2  ! Principal moments of inertia of Mars
      moiC=0.365905*mass*radius**2  ! Principal moments of inertia of Mars
      !----------------------------------------------------------------------------------!
      ! Lithospheric thickness and corresponding kf (based on Zharkov and
      ! Gudkova, 2005)
      !    15    |    48    |    58    |    84    |    110   |    164    |
      !    200    [km]
      ! 1.203959 | 1.127384 | 1.110566 | 1.065899 | 1.023186 | 0.9458358 |
      ! 0.8986673
      !----------------------------------------------------------------------------------!
      kf = 0.899                    ! Fluid (Tidal) Love number

   end subroutine mars_init

end module planets_mod
