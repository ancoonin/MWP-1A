module user_specs_mod

!Planetary model directory

character(*), parameter :: planetFolder = 'INPUT_OTHERS/'

! Common file extensions
!   ! Since the code does not explicitly specify the precision (single or
!   ! double) of the input/output files, an
!   !  extension could be used to designate the type. Alternatively, the
!   !  files could also be delimited text files (to
!   !  be implemented later). Specify the common extension here
!   !  (applicable to the names of ALL input/output files):
   character(4), parameter :: ext = '.dat'    ! '.sgl' | '.dbl' | '.txt'| ''
   ! ... and their file type:
   character(*), parameter :: fType = 'text'  ! 'binary' | 'text'

   ! Planet selection

   character(*), parameter :: whichPlanet      = 'earth'                  !'earth', 'Mars', etc.
   character(*), parameter :: iceModel         = 'iceload'              ! Name of the ice files without incremental suffix
   character(*), parameter :: topoModel        = 'etopo2_512_ice5g'       !Bedrock topography at time = 0;
                                                                          !  without the filename extension
   character(*), parameter :: timeArray        = 'times'                  ! Name of the time array file
   character(*), parameter :: topoModelInit    = 'tgrid1'  ! Name of the initial topography file (if known)


  ! Model parameters
   integer, parameter :: nsteps = 108          ! Number of ice files
   integer, parameter :: numOuter = 1          ! Number of outer loops (enter 99 for convergence check)
   logical, parameter :: initialTopo = .true.  ! .true. to input known initialtopography
                                                !  .false. to use modern
                                                !  topography as first guess for
                                                !  topography

   logical, parameter :: calcRG = .false.       ! .true. to calculate the radial and geoid displacements; 
                                                !    note that
                                                !    the "true" option only
                                                !    works for a fixed number of
                                                !    outer loops
                                                !    (i.e., no convergence
                                                !    checks!).
                                                ! .false. to only calculate RSL.
   logical, parameter :: tGridOut = .false.      ! .true. to output topography grids
   logical, parameter :: tpw = .true.           ! .true. to incorporate rotational feedback
                                                ! .false. for non-rotating
                                                ! planet
   logical, parameter :: iceVolume = .false.     ! .true. to output ice volume at each time step
   logical, parameter :: checkMarine = .true.  ! .true. to check for floating marine-based ice
                                                ! .false. to assume all ice is
                                                ! grounded
   integer, parameter :: norder = 256           ! Max spherical harmonic degree/order
   integer, parameter :: npam = 500             ! Max relaxation modes
   integer, parameter :: nglv = 512             ! Number of GL points in latitude

   real, parameter :: epsilon1 = 5.0E-5         ! Inner loop convergence criterion
   real, parameter :: epsilon2 = 1.0E-5         ! Outer loop convergence criterion
                                                !  (if doing a convergence check
                                                !  for outer loop, see below)

end module user_specs_mod

