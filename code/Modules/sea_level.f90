module sea_level

use spharmt  ! Spherical harmonic transform module
use user_specs_mod ! User specifications
use constants_mod ! Physical and mathematical parameters
use planets_mod ! Earth specific parameters

! ====================================================================== !
!                                                                        !
!                                                                        !
! Coonin (2023): This module contains the "checkpoint" forward sea-level !
!                subroutine used in the coupled NA + MWP-1A inversion    !
!                                                                        !
! Pseudo-spectral ice age sea level model, using algorithm of Kendall    !
! et al, 2005: On post-glacial sea level - II...                         !
!                                                                        !
! Original version of sl_model written by Sam Goldberg, Harvard          !
! University EPS '16, senior thesis advised by Jerry Mitrovica           !
! 12/2015-1/2016                                                         !
!                                                                        !                       
! Since, modified by Erik Chan, Holly Han in 2018 and finally by A.N.    !
! Coonin in 2023 for compatibility with Neighborhood Algorithm for       !
! MWP-1A sea-level fingerprint with transient viscoelastic solid Earth   !
! deformation.                                                           !               
!                                                                        ! 
! See general notes below for more information                           !
!                                                                        !               
!                                                                        !
! ====================================================================== !


! general notes about sl_model:

! Pseudo-spectral ice age sea level model, using algorithm of Kendall et al,
! 2005: On post-glacial sea level - II...
! Original version written by Sam Goldberg, Harvard University EPS '16 for
! senior thesis
! Advised by Jerry Mitrovica
! 12/2015-1/2016

! Modified in 2018 by Erik Chan, McGill University, Post Doctoral Fellow
! Modified in 2018 by Holly Han, McGill University, PhD Student

! CHANGES (prior to tracking with Git on 2018-01-15):
!  2017-11-22:
!  > Replaced "form = 'binary'" with "form = 'unformatted', access = 'stream',
!  status = '<old>/<replace>', where the
!     last option depends on specific read/write needs." This is to bring the
!     code into compliance with GNU compilers.
!  > Added "status = 'old'" to READ statement for planetary model Love number
!  files (line 174 at the time of this
!     writing).
!  > Ensure everything fits within 120 columns.
!  > Cleaned up the indentations for do...enddo, if...endif, etc.
!  by Erik N.-H. Chan
!  2018-01-12:
!  > Added if..else block around lines 519-525 (at the time of writing) to avoid
!  NaN (= 0 / 0 or dividing by 0).
!  by Erik N.-H. Chan

! For all variables except SL, the prefix d- is equivalent to δ- in Kendall,
! i.e. incremental change over one timestep.
! The prefix delta- is equivalent to Δ, i.e. total change since first timestep.
! SL is different - dSL refers to spatially heterogeneous change only (the
! script SL in Kendall), and deltaSL is the
!  total change including eustatic - both since first timestep
! The suffix -xy denotes a quantity in the spatial domain; -lm denotes a
! quantity in the spectral domain.
! The prefix old- generally refers to the previous timestep, used to calculate
! increments, although sometimes
!  to previous iteration to check convergence.
! The suffix -0 generally refers to the value at first timestep, used to
! calculate total changes.
! nouter refers to the iteration of the outer loop (full time series).
! ninner refers to the inner loop within each timestep to converge the ocean
! load.
! n refers to the timestep, between 1 and nsteps.
! Other than the above, I have endeavored to be faithful to the notation of
! Kendall et al.

! The INPUT directory should contain the following files:
!  - Ice: named 'ice1','ice2','ice3',...,'ice(nsteps)', containing ice thickness
!  in metres on a Gauss-Legendre grid.
!  ********************************************************************************************************************
!  *!!! CAUTION: * This code assumes that the first ice file starts with suffix
!  1. If you have an ice file with       *
!  *               suffix 0, modify the variable 'n' to 'n-1' on
!  'write(numstr,'(I3)') n' (LINE 301)!!!]              *
!  ********************************************************************************************************************
!  - Times: named 'times', containing the times, in years, of each ice files.
!  Present is 0, past is negative, etc.
!  - Present-day bedrock topography: named defined in the "topoModel" variable
!  in the user_specs_mod module, with
!     values in metres on a Gauss-Legendre grid,
!  - Love numbers in JXM's maxwell.f output format: name is specified in the
!  "planetModel" variable in user_specs_mod.
! The OUTPUT directory will contain the following files:
!  - 'SL1','SL2',...: Total sea-level change from the start of the simulation,
!  in metres, on a Gauss-Legendre grid,
!  - 'times':  Times in years of each output file,
!  - 'R1','G1','R2','G2',...: Radial and Geoid components of RSL change (only if
!  parameter calcRG is enabled).
! NOTE ON GRID FILES: All of the grid files (i.e., ice, topography, sea-level,
! radial, geoid) could be in ASCII text or
!  flat binary formats. For these files, the order of values in read and write
!  operations are as follows: Along lines
!  of equal latitude, from 0 degrees towards 360 degrees of longitude; each of
!  these lines are then ordered from 0
!  colatitude towards 180 (i.e., from North to South).
!
!
!============ End of general notes=========================================== 


implicit none

contains

subroutine read_config(planetModel,path2loveinputs,path2SLinputs,path2iceinputs)

CHARACTER(LEN=255) :: planetModel,path2loveinputs,path2SLinputs,path2iceinputs

!NOTE THAT THIS FILE MUST BE LOCATED IN THE DIRECTORY YOU ARE EXECUTING FROM!
open(unit=1,file='path_config',form='formatted',access='sequential',status='old')

read(1,*) !do nothing
read(1,*) !do nothing

read(1,*) !planetModel to use in sl code:
read(1,*) planetModel
read(1,*)

read(1,*) !path to planetFolder ('INPUT_OTHERS') which contains love number inputs for SL code
read(1,*) path2loveinputs
read(1,*)

read(1,*) !path to initial topography file and all variables saved for SL checkpoint
read(1,*) path2SLinputs
read(1,*)

read(1,*) !path to driver inputs (pre-MWP-1A ice masks for each ice sheet)
read(1,*) path2iceinputs
read(1,*)

  close(1)

return
end subroutine read_config

subroutine sl_model(icegrid,ntimesduringMWP,times,predicted_data)


!======================================================================================================================!
!                                                      MAIN BLOCK                                                      !
!                                                                                                                      !
!______________________________________________________________________________________________________________________!

!=======================================================================================================MAIN
!subroutine sl_model(icegrid,ntimesduringMWP,times,predicted_data)

! INPUTS: icegrid(:,:,:), ntimesduringMWP, times
! OUTPUTS: predicted_data(:)

!______________________________________________________________________________________________________________________!
implicit none

!==========================================================
!                        VARIABLES
!___________________(Edit with caution)____________________

! Inputs
real, allocatable :: icexy(:,:,:)

real, allocatable,intent(in) :: times(:)

real, allocatable :: rprime(:,:)
real, allocatable :: r(:,:)
real, allocatable :: s(:,:)

real, allocatable :: ke(:)
real, allocatable :: he(:)

real, allocatable :: rprimeT(:,:)
real, allocatable :: rT(:,:)

real, allocatable :: kTE(:)
real, allocatable :: hTE(:)

! Model calculations
real, allocatable :: glw_matrix(:,:)

real, allocatable :: initTopo(:,:)

real, allocatable :: sl(:,:,:)
real, allocatable :: topoxy(:,:,:)
real, allocatable :: cxy(:,:,:)
                                                         !  total spatially
                                                         !  heterogeneous sea
                                                         !  level change
real, allocatable :: deltaslxy(:,:)
real, allocatable :: dslxy(:,:)

                                                        ! Grounded ice
                                                        !thickness
real, allocatable :: icestarxy(:,:)

                                                        ! Grounded ice mask,
                                                        !ice-free ocean function
real, allocatable :: beta(:,:)
real, allocatable :: cstarxy(:,:)
real, allocatable :: cstar0(:,:)

!initial guess to ocean function for iterative convergence
real, allocatable :: cxy_guess(:,:)

! Projections used to
!calculate loads and shoreline migration
real, allocatable :: tOxy(:,:)
real, allocatable :: rOxy(:,:)
real, allocatable :: tTxy(:,:)

complex, allocatable :: cstarlm(:,:)
complex, allocatable :: oldcstarlm(:,:)
complex, allocatable :: tOlm(:,:)
complex, allocatable :: rOlm(:,:)
complex, allocatable :: dSlm(:,:)
complex, allocatable :: olddSlm(:,:)
complex, allocatable :: icestarlm(:,:)
complex, allocatable :: dicestarlm(:,:)
complex, allocatable :: deltaicestarlm(:,:)
complex, allocatable :: oldicestarlm(:,:)
complex, allocatable :: icestar0(:,:)
complex, allocatable :: t0lm(:,:)
complex, allocatable :: oldt0lm(:,:)
complex, allocatable :: tTlm(:,:)
complex, allocatable :: oldtTlm(:,:)
complex, allocatable :: dsllm(:,:)
complex, allocatable :: deltasllm(:,:)

!! Big arrays of changes in loads, used in Love number viscous response 

complex, allocatable :: dicestar(:,:,:)
complex, allocatable :: dS(:,:,:)
complex, allocatable :: deltaicestar(:,:,:)
complex, allocatable :: deltaS(:,:,:)

real :: conserv                                          ! Uniform geoid shift
                                                         !(ΔΦ/g)
real :: ttl, ekhl                                        ! Used in Love number
                                                         ! calculations

real, allocatable :: lovebeta(:,:)

complex, allocatable :: viscous(:,:)

real :: xi, zeta                                         ! Convergence checks

real, allocatable :: ice_volume(:)

! For calculating R and G separately

real, allocatable :: rrxy(:,:)
real, allocatable :: rr(:,:,:)

complex, allocatable :: rrlm(:,:)

real, allocatable :: lovebetarr(:,:)

complex :: viscousrr

! Rotation calculations
! real, parameter :: --------------------------------->     ! Fluid Love number
                                                            ! is defined in the planetary modules
                                                            ! Rotational
                                                            ! driving
complex, allocatable :: lambda(:,:)
complex, allocatable :: oldlambda(:,:)
complex, allocatable :: lambda0(:,:)

                                                            ! Big arrays of
                                                            ! changes in rotational driving
complex, allocatable :: dlambda(:,:,:)
complex, allocatable :: deltalambda(:,:,:)

                                                            ! Rotational
                                                            !perturbation vector
real, allocatable :: mm(:)
real, allocatable :: sum_m(:)
real, allocatable :: oldm(:)

real, allocatable :: dm(:,:)

                                                            ! Load
                                                            ! perturbations to MOI
real, allocatable :: il(:,:)
real, allocatable :: oldil(:,:)
real, allocatable :: sum_il(:,:)

real, allocatable :: dil(:,:,:)
                                                            ! sea-level change
                                                            ! from rotational feedbacks
complex, allocatable :: dsl_rot(:,:)
complex, allocatable :: rr_rot(:,:)

real :: ekhTE                                               ! Used in Love
                                                            ! number calculations

real, allocatable :: lovebetatt(:)
real, allocatable :: lovebetattrr(:)

real :: betatt, betattprime                                 ! Used in Love number calculations
complex :: viscoustt,viscousttrr                            ! Used in Love number calculations

! Miscellaneous variables
integer :: ninner,nouter                                    ! Iteration of inner and outer loops
integer :: i,j,k,l,m,n,nn                                   ! Do-loop indices
integer :: counti, countf,counti1,countf1,countrate         ! Computation timing
type(sphere) :: spheredat                                   ! SH transform data to be passed to subroutines

! For JXM's code to read in Love numbers
integer :: legord(norder),nmod(norder),nmodes(norder),ll,nm
real :: xn

real, allocatable :: elast(:,:)
real, allocatable :: asymv(:,:)
real, allocatable :: telast(:,:)
real, allocatable :: tasymv(:,:)

real, allocatable :: resh(:,:)
real, allocatable :: resl(:,:)
real, allocatable :: resk(:,:)
real, allocatable :: tresh(:,:)
real, allocatable :: tresl(:,:)
real, allocatable :: tresk(:,:)

real :: taurr,taurt,dmx,dmy

! AC additions (2023):

! 512x1024 Gauss-Legendre grid coordinates for each data site
integer,parameter :: tahiti_idx_lat = 307
integer,parameter :: tahiti_idx_lon = 599
integer,parameter :: barbados_idx_lat = 219
integer,parameter :: barbados_idx_lon = 855
integer,parameter :: sunda_idx_lat = 244
integer,parameter :: sunda_idx_lon = 310
integer,parameter :: hyd_idx_lat = 312
integer,parameter :: hyd_idx_lon = 428
integer,parameter :: nog_idx_lat = 305
integer,parameter :: nog_idx_lon = 417
integer,parameter :: nws_idx_lat = 94
integer,parameter :: nws_idx_lon = 1008

CHARACTER(LEN=255) :: planetModel,path2loveinputs,path2SLinputs,path2iceinputs
integer, parameter :: first_tstep = 72 ! first MWP1A timestep

real, allocatable :: icestar0_real(:,:)
real, allocatable :: icestar0_imag(:,:)
real, allocatable :: dS_slice(:,:)
real, allocatable :: deltaS_slice(:,:)
real, allocatable :: dicestar_slice(:,:)

real, allocatable :: oldicestarlm_real(:,:)
real, allocatable :: oldicestarlm_imag(:,:)
real, allocatable :: oldcstarlm_real(:,:)
real, allocatable :: oldcstarlm_imag(:,:)

real, allocatable :: lambda0_real(:,:)
real, allocatable :: lambda0_imag(:,:)
real, allocatable :: oldlambda_imag(:,:)
real, allocatable :: oldlambda_real(:,:)

real, allocatable :: dlambda_real(:,:,:)
real, allocatable :: dlambda_imag(:,:,:)
real, allocatable :: deltalambda_real(:,:,:)
real, allocatable :: deltalambda_imag(:,:,:)

real, allocatable :: dS_real(:,:,:)
real, allocatable :: dS_imag(:,:,:)
real, allocatable :: deltaS_real(:,:,:)
real, allocatable :: deltaS_imag(:,:,:)
real, allocatable :: dicestar_real(:,:,:)
real, allocatable :: dicestar_imag(:,:,:)

real, allocatable :: dil_slice(:,:)
real, allocatable :: dlambda_slice(:,:)
real, allocatable :: deltalambda_slice(:,:)

integer, intent(in) :: ntimesduringMWP
real, dimension(6,ntimesduringMWP), intent(out) :: predicted_data
real, dimension(ntimesduringMWP) :: RSL_scotland, RSL_scotland_mask
logical :: SLO_check

real,allocatable,intent(in) :: icegrid(:,:,:)

integer :: ios       !ios = failure flag for memory deallocation
integer :: convflag  !convflag = sea level model failed to converge in inner loop 

character(LEN=3) :: tindx
! Planetary values
if (whichPlanet == 'earth' .or. whichPlanet == 'Earth' .or. whichPlanet =='EARTH') then
   call earth_init
elseif (whichPlanet == 'mars' .or. whichPlanet == 'Mars' .or. whichPlanet =='MARS') then
   call mars_init
else
   write(*,'(A)') 'The parameters for the planet you entered are not built in.'
   write(*,'(A)') 'Please check your spelling for the variable whichPlanet in the user_specs_mod module.'
   write(*,'(A)') 'If you preferred, you could create a new subroutine for the new planet in the planets_mod module.'
   write(*,'(A)') 'Terminating: program sl_model'
   stop
endif

! Sanity check
if (.not. initialTopo .and. NumOuter == 1) then ! If only 1 outer loop but the
                                                !initial topography option is turned off !--AC
   write(*,'(A)') 'ERROR: You have set the number of outer loops "numOuter" (i.e., topography iteration) to 1.'
   write(*,'(A)') '       However, this is only designed to work with a known initial topography. '
   write(*,'(A)') '       You should set "initialTopo" to .true., and'
   write(*,'(A)') '       provide an initial topography file in "topoModelInit".'
   write(*,'(A)') '!!!---- Program sl_model will now be terminated. ----!!!'
   call abort
elseif (initialTopo .and. NumOuter > 1) then    ! If initial topography is
                                                !turned on, but there are multiple outer loops !--AC
   write(*,'(A)') 'NOTE: You have set the number of outer loops "numOuter" (i.e., topography iteration) to more than 1.'
   write(*,'(A)') '      Since you have set "initialTopo" to .true. and provided an initial topography file,'
   write(*,'(A)') '      this is not needed and may take extra computational time.'
endif

call system_clock(count = counti, count_rate = countrate)  ! Total computation time
call spharmt_init(spheredat, 2*nglv, nglv, norder, radius) ! Initialize spheredat (for SH transform subroutines)

call read_config(planetModel,path2loveinputs,path2SLinputs,path2iceinputs) !--AC

!--AC: allocating memory for large local arrays
allocate(icexy(nglv,2*nglv,nsteps),r(npam,norder),s(npam,norder),rprime(npam,norder))
allocate(ke(norder),he(norder),rprimeT(npam,norder),rT(npam,norder),kTE(norder))
allocate(hTE(norder),initTopo(nglv,2*nglv),glw_matrix(nglv,2*nglv),icestarxy(nglv,2*nglv))
allocate(sl(nglv,2*nglv,nsteps),topoxy(nglv,2*nglv,nsteps),cxy(nglv,2*nglv,nsteps))
allocate(deltaslxy(nglv,2*nglv),dslxy(nglv,2*nglv),beta(nglv,2*nglv),cstarxy(nglv,2*nglv))
allocate(cstar0(nglv,2*nglv),cxy_guess(nglv,2*nglv),tOxy(nglv,2*nglv),rOxy(nglv,2*nglv))
allocate(tTxy(nglv,2*nglv),cstarlm(0:norder,0:norder),oldcstarlm(0:norder,0:norder))
allocate(tOlm(0:norder,0:norder),rOlm(0:norder,0:norder),dSlm(0:norder,0:norder))
allocate(olddSlm(0:norder,0:norder),icestarlm(0:norder,0:norder),dicestarlm(0:norder,0:norder))
allocate(deltaicestarlm(0:norder,0:norder),oldicestarlm(0:norder,0:norder))
allocate(icestar0(0:norder,0:norder),t0lm(0:norder,0:norder),oldt0lm(0:norder,0:norder))
allocate(tTlm(0:norder,0:norder),oldtTlm(0:norder,0:norder),dsllm(0:norder,0:norder))
allocate(deltasllm(0:norder,0:norder),dicestar(0:norder,0:norder,nsteps))
allocate(dS(0:norder,0:norder,nsteps),deltaicestar(0:norder,0:norder,nsteps))
allocate(deltaS(0:norder,0:norder,nsteps),lovebeta(nsteps,norder))
allocate(viscous(0:norder,0:norder),ice_volume(nsteps),rrxy(nglv,2*nglv))
allocate(rr(nglv,2*nglv,nsteps),rrlm(0:norder,0:norder),lovebetarr(nsteps,norder))
allocate(lambda(0:2,0:2),oldlambda(0:2,0:2),lambda0(0:2,0:2),dlambda(0:2,0:2,nsteps))
allocate(deltalambda(0:2,0:2,nsteps),mm(3),sum_m(3),oldm(3),dm(3,nsteps))
allocate(il(3,3),oldil(3,3),sum_il(3,3),dil(3,3,nsteps))
allocate(dsl_rot(0:2,0:2),rr_rot(0:2,0:2),lovebetatt(nsteps),lovebetattrr(nsteps))
allocate(elast(3,norder),asymv(3,norder),telast(3,norder),tasymv(3,norder))
allocate(resk(npam,norder),resl(npam,norder),resh(npam,norder))
allocate(tresk(npam,norder),tresl(npam,norder),tresh(npam,norder))
allocate(icestar0_real(0:norder,0:norder),icestar0_imag(0:norder,0:norder))
allocate(dS_slice(0:norder,0:norder),deltaS_slice(0:norder,0:norder))
allocate(dicestar_slice(0:norder,0:norder),oldicestarlm_real(0:norder,0:norder))
allocate(oldicestarlm_imag(0:norder,0:norder),oldcstarlm_real(0:norder,0:norder))
allocate(oldcstarlm_imag(0:norder,0:norder),lambda0_real(0:2,0:2),lambda0_imag(0:2,0:2))
allocate(oldlambda_real(0:2,0:2),oldlambda_imag(0:2,0:2))
allocate(dlambda_real(3,3,first_tstep-1),dlambda_imag(3,3,first_tstep-1))
allocate(deltalambda_real(3,3,first_tstep-1),deltalambda_imag(3,3,first_tstep-1))
allocate(dS_real(0:norder,0:norder,first_tstep-1),dS_imag(0:norder,0:norder,first_tstep-1))
allocate(deltaS_real(0:norder,0:norder,first_tstep-1),deltaS_imag(0:norder,0:norder,first_tstep-1))
allocate(dicestar_real(0:norder,0:norder,first_tstep-1),dicestar_imag(0:norder,0:norder,first_tstep-1))
allocate(dil_slice(3,3),dlambda_slice(3,3), deltalambda_slice(3,3))

!initializing allocated arrays 
icexy(:,:,:)=0.0
r(:,:)=0.0
s(:,:)=0.0
rprime(:,:)=0.0
ke(:)=0.0
he(:)=0.0
rprimeT(:,:)=0.0
rT(:,:)=0.0
kTE(:)=0.0
hTE(:)=0.0
initTopo(:,:)=0.0
glw_matrix(:,:)=0.0
icestarxy(:,:)=0.0
sl(:,:,:)=0.0
topoxy(:,:,:)=0.0
cxy(:,:,:)=0.0
deltaslxy(:,:)=0.0
dslxy(:,:)=0.0
beta(:,:)=0.0
cstarxy(:,:)=0.0
cstar0(:,:)=0.0
cxy_guess(:,:)=0.0
tOxy(:,:)=0.0
rOxy(:,:)=0.0
tTxy(:,:)=0.0
cstarlm(:,:)=cmplx(0,0)
oldcstarlm(:,:)=cmplx(0,0)
tOlm(:,:)=cmplx(0,0)
rOlm(:,:)=cmplx(0,0)
dSlm(:,:)=cmplx(0,0)
olddSlm(:,:)=cmplx(0,0)
icestarlm(:,:)=cmplx(0,0)
dicestarlm(:,:)=cmplx(0,0)
deltaicestarlm(:,:)=cmplx(0,0)
oldicestarlm(:,:)=cmplx(0,0)
icestar0(:,:)=cmplx(0,0)
t0lm(:,:)=cmplx(0,0)
oldt0lm(:,:)=cmplx(0,0)
tTlm(:,:)=cmplx(0,0)
oldtTlm(:,:)=cmplx(0,0)
dsllm(:,:)=cmplx(0,0)
deltasllm(:,:)=cmplx(0,0)
dicestar(:,:,:)=cmplx(0,0)
dS(:,:,:)=cmplx(0,0)
deltaicestar(:,:,:)=cmplx(0,0)
deltaS(:,:,:)=cmplx(0,0)
lovebeta(:,:)=0.0
viscous(:,:)=cmplx(0,0)
ice_volume(:)=0.0
rrxy(:,:)=0.0
rr(:,:,:)=0.0
rrlm(:,:)=cmplx(0,0)
lovebetarr(:,:)=0.0
lambda(:,:)=cmplx(0,0)
oldlambda(:,:)=cmplx(0,0)
lambda0(:,:)=cmplx(0,0)
dlambda(:,:,:)=cmplx(0,0)
deltalambda(:,:,:)=cmplx(0,0)
mm(:)=0.0
sum_m(:)=0.0
oldm(:)=0.0
dm(:,:)=0.0
il(:,:)=0.0
oldil(:,:)=0.0
sum_il(:,:)=0.0
dil(:,:,:)=0.0
dsl_rot(:,:)=cmplx(0,0)
rr_rot(:,:)=cmplx(0,0)
lovebetatt(:)=0.0
lovebetattrr(:)=0.0
elast(:,:)=0.0
asymv(:,:)=0.0
telast(:,:)=0.0
tasymv(:,:)=0.0
resk(:,:)=0.0
resl(:,:)=0.0
resh(:,:)=0.0
tresk(:,:)=0.0
tresl(:,:)=0.0
tresh(:,:)=0.0
icestar0_real(:,:)=0.0
icestar0_imag(:,:)=0.0
dS_slice(:,:)=0.0
deltaS_slice(:,:)=0.0
dicestar_slice(:,:)=0.0
oldicestarlm_real(:,:)=0.0
oldicestarlm_imag(:,:)=0.0
oldcstarlm_real(:,:)=0.0
oldcstarlm_imag(:,:)=0.0
lambda0_real(:,:)=0.0
lambda0_imag(:,:)=0.0
oldlambda_real(:,:)=0.0
oldlambda_imag(:,:)=0.0
dlambda_real(:,:,:)=0.0
dlambda_imag(:,:,:)=0.0
deltalambda_real(:,:,:)=0.0
deltalambda_imag(:,:,:)=0.0
dS_real(:,:,:)=0.0
dS_imag(:,:,:)=0.0
deltaS_real(:,:,:)=0.0
deltaS_imag(:,:,:)=0.0
dicestar_real(:,:,:)=0.0
dicestar_imag(:,:,:)=0.0
dil_slice(:,:)=0.0
dlambda_slice(:,:)=0.0
deltalambda_slice(:,:)=0.0

!-----------------------------------------------------------
!                     Read input files
!-----------------------------------------------------------
write(*,'(A)') 'Reading input files'

! Ice loads are directly passed as subroutine arguments

icexy(:,:,first_tstep:nsteps)=icegrid(:,:,1:ntimesduringMWP)

! Initial topography
if (initialTopo) then      ! If the initial topography is known...
   write(*,*) 'Reading in the known initial topo file'
   if (fType == 'text') then        ! If the files are space-delimited text files
      open(unit = 1, file = trim(path2SLinputs)//'tgrid1.dat', form ='formatted', access = 'sequential', &
      & status = 'old') !AC: reading in converged initial topo from full run to present
      
      do i=1,512
      read(1,*) (initTopo(i,j),j=1,1024)
      enddo

   close(1)
   endif
endif

!-----------------------------------------------------------
!  Read in Love numbers 
! (in the form of JXM's output from 'maxwell.f')
!-----------------------------------------------------------

open(unit = 2, file = trim(path2loveinputs)//planetFolder//planetModel, status ='old')
! Following code borrowed from JXM
read(2,*)
do j = 1,norder
   read(2,*) legord(j), nmodes(j)
   nm = nmodes(j)
   xn = real(legord(j))
   ll = legord(j)
   nmod(ll) = nm
   read(2,*) (s(i,ll), i=1, nm)
   read(2,*) elast(1,ll), elast(2,ll), taurr, taurt, elast(3,ll)
   read(2,*) asymv(1,ll), asymv(2,ll), taurr, taurt, asymv(3,ll)
   read(2,*) (resh(i,ll), i=1,nm)
   read(2,*) (resl(i,ll), i=1,nm)
   read(2,*) (resk(i,ll), i=1,nm)
   if (xn .lt. 1.5) cycle
   read(2,*) telast(1,ll), telast(2,ll), dmx, dmy, telast(3,ll)
   read(2,*) tasymv(1,ll), tasymv(2,ll), dmx, dmy, tasymv(3,ll)
   read(2,*) (tresh(i,ll), i=1, nm)
   read(2,*) (tresl(i,ll), i=1, nm)
   read(2,*) (tresk(i,ll), i=1, nm)
enddo
close(2)

! Divide by l (numbers are multiplied by l in JXM's output)
do l = 1,norder
   resl(:,l) = resl(:,l) / real(l)
   resk(:,l) = resk(:,l) / real(l)
   tresl(:,l) = tresl(:,l) / real(l)
   tresk(:,l) = tresk(:,l) / real(l)
enddo

! Assign Love number inputs to the letters used in Kendall et al
do l = 1,norder
   he(l) = elast(1,l)
   ke(l) = elast(3,l) / real(l)
   hTE(l) = telast(1,l)
   kTE(l) = telast(3,l) / real(l)
enddo
rprime(:,:) = resk(:,:)
r(:,:) = resh(:,:)
rprimeT(:,:) = tresk(:,:)
rT(:,:) = tresh(:,:)

!AC: Reading in presaved information related to ice + ocean loading prior to MWP-1A

open(unit = 99, file = trim(path2SLinputs)//'saved_vars.dat', form ='formatted', status='old')

    read(99,*) !'real(icestar0)'
        do i=0,norder
        read(99,"(257ES17.8E2)") (icestar0_real(i,j),j=0,norder)
        enddo

    read(99,*) !'imag(icestar0)'
        do i=0,norder
        read(99,"(257ES17.8E2)") (icestar0_imag(i,j),j=0,norder)
        enddo

    icestar0(:,:)=icestar0_real(:,:)+ii*icestar0_imag

    read(99,*) !'cxy_guess'
        do i=1,512
        read(99,"(1024ES17.8E2)") (cxy_guess(i,j),j=1,1024)
        enddo

    read(99,*) !'cstar0'
        do i=1,512
        read(99,"(1024ES17.8E2)") (cstar0(i,j),j=1,1024)
        enddo

    read(99,*) !'real(lambda0)'
        do i=0,2
        read(99,"(3ES17.8E2)") (lambda0_real(i,j),j=0,2)
        enddo

    read(99,*) !'imag(lambda0)'
        do i=0,2
        read(99,"(3ES17.8E2)") (lambda0_imag(i,j),j=0,2)
        enddo

    lambda0(:,:)=lambda0_real(:,:)+ii*lambda0_imag(:,:)

    read(99,*) !'real(oldicestarlm)'
        do i=0,norder
        read(99,"(257ES17.8E2)") (oldicestarlm_real(i,j),j=0,norder)
        enddo

    read(99,*) !'imag(oldicestarlm)'
        do i=0,norder
        read(99,"(257ES17.8E2)") (oldicestarlm_imag(i,j),j=0,norder)
        enddo

    oldicestarlm(:,:)=oldicestarlm_real(:,:)+ii*oldicestarlm_imag(:,:)


    read(99,*) !'real(oldcstarlm)'
        do i=0,norder
        read(99,"(257ES17.8E2)") (oldcstarlm_real(i,j),j=0,norder)
        enddo

    read(99,*) !'imag(oldcstarlm)'
        do i=0,norder
        read(99,"(257ES17.8E2)") (oldcstarlm_imag(i,j),j=0,norder)
        enddo

    oldcstarlm(:,:)=oldcstarlm_real(:,:)+ii*oldcstarlm_imag(:,:)

    read(99,*) !'oldil'
        do i=1,3
        read(99,"(3ES17.8E2)") (oldil(i,j),j=1,3)
        enddo

    read(99,*) !'dil'
        do k=1,first_tstep-1
        read(99,*) !'timestep k'
            do i=1,3
            read(99,"(3ES17.8E2)") (dil_slice(i,j),j=1,3)
            enddo
        dil(:,:,k)=dil_slice(:,:)
        enddo

    read(99,*) !'oldm'
        read(99,"(3ES17.8E2)") oldm

    read(99,*) !'real(oldlambda)'
        do i=0,2
        read(99,"(3ES17.8E2)") (oldlambda_real(i,j),j=0,2)
        enddo

    read(99,*) !'imag(oldlambda)'
        do i=0,2
        read(99,"(3ES17.8E2)") (oldlambda_imag(i,j),j=0,2)
        enddo

    oldlambda(:,:)=oldlambda_real(:,:)+ii*oldlambda_imag(:,:)


    read(99,*) !'real(dlambda)'
        do k=1,first_tstep-1
        read(99,*) !'timestep k'
            do i=1,3
            read(99,"(3ES17.8E2)") (dlambda_slice(i,j),j=1,3)
            enddo
        dlambda_real(:,:,k)=dlambda_slice(:,:)
        enddo

    read(99,*) !'imag(dlambda)'
        do k=1,first_tstep-1
        read(99,*) !'timestep k'
            do i=1,3
            read(99,"(3ES17.8E2)") (dlambda_slice(i,j),j=1,3)
            enddo
        dlambda_imag(:,:,k)=dlambda_slice(:,:)
        enddo

    dlambda(:,:,1:first_tstep-1)=dlambda_real(:,:,:)+ii*dlambda_imag(:,:,:)


   read(99,*) !'real(deltalambda)'
        do k=1,first_tstep-1
        read(99,*) !'timestep k'
            do i=1,3
            read(99,"(3ES17.8E2)") (deltalambda_slice(i,j),j=1,3)
            enddo
        deltalambda_real(:,:,k)=deltalambda_slice(:,:)
        enddo

    read(99,*) !'imag(deltalambda)'
        do k=1,first_tstep-1
        read(99,*) !'timestep k'
            do i=1,3
            read(99,"(3ES17.8E2)") (deltalambda_slice(i,j),j=1,3)
            enddo
        deltalambda_imag(:,:,k)=deltalambda_slice(:,:)
        enddo

    deltalambda(:,:,1:first_tstep-1)=deltalambda_real(:,:,:)+ii*deltalambda_imag(:,:,:)

  read(99,*) !'dm'
    do i=1,3
    read(99,"(71ES17.8E2)") (dm(i,j),j=1,first_tstep-1) !AC hardcoded number of preCheckpoint timesteps : 71
    enddo

close(99)


open(unit = 98, file = trim(path2SLinputs)//'dS.dat', form = 'formatted',status='old')

    read(98,*) !'real(dS)'
        do k=1,first_tstep-1
        read(98,*) !'timestep k'
            do i=0,norder
            read(98,"(257ES17.8E2)") (dS_slice(i,j),j=0,norder)
            enddo
        dS_real(:,:,k)=dS_slice(:,:)
        enddo

    read(98,*) !'imag(dS)'
        do k=1,first_tstep-1
        read(98,*) !'timestep k'
            do i=0,norder
            read(98,"(257ES17.8E2)") (dS_slice(i,j),j=0,norder)
            enddo
        dS_imag(:,:,k)=dS_slice(:,:)
        enddo

    dS(:,:,1:first_tstep-1)=dS_real(:,:,:)+ii*dS_imag(:,:,:)


close(98)

open(unit = 97, file = trim(path2SLinputs)//'deltaS.dat', form = 'formatted',status='old')

    read(97,*) !'real(deltaS)'
        do k=1,first_tstep-1
        read(97,*) !'timestep k'
            do i=0,norder
            read(97,"(257ES17.8E2)") (deltaS_slice(i,j),j=0,norder)
            enddo
        deltaS_real(:,:,k)=deltaS_slice(:,:)
        enddo

    read(97,*) !'imag(deltaS)'
        do k=1,first_tstep-1
        read(97,*) !'timestep k'
            do i=0,norder
            read(97,"(257ES17.8E2)") (deltaS_slice(i,j),j=0,norder)
            enddo
        deltaS_imag(:,:,k)=deltaS_slice(:,:)
        enddo

    deltaS(:,:,1:first_tstep-1)=deltaS_real(:,:,:)+ii*deltaS_imag(:,:,:)

    !reading in dicestar

    read(97,*) !'real(dicestar)'
        do k=1,first_tstep-1
        read(97,*) !'timestep k'
            do i=0,norder
            read(97,"(257ES17.8E2)") (dicestar_slice(i,j),j=0,norder)
            enddo
        dicestar_real(:,:,k)=dicestar_slice(:,:)
        enddo

    read(97,*) !'imag(dicestar)'
        do k=1,first_tstep-1
        read(97,*) !'timestep k'
            do i=0,norder
            read(97,"(257ES17.8E2)") (dicestar_slice(i,j),j=0,norder)
            enddo
        dicestar_imag(:,:,k)=dicestar_slice(:,:)
        enddo

    dicestar(:,:,1:first_tstep-1)=dicestar_real(:,:,:)+ii*dicestar_imag(:,:,:)

close(97)

!Deallocating memory from the variables used solely to read in the presaved
!info regarding ice + ocean loading prior to MWP-1A

  deallocate(icestar0_real,icestar0_imag,dS_slice,STAT=ios)
  deallocate(dicestar_slice,oldicestarlm_real,oldicestarlm_imag,STAT=ios)
  deallocate(oldcstarlm_real,oldcstarlm_imag,STAT=ios)
  deallocate(lambda0_real,lambda0_imag,oldlambda_imag,oldlambda_real,STAT=ios)
  deallocate(dlambda_real,dlambda_imag,deltalambda_real,deltalambda_imag,STAT=ios)
  deallocate(dS_real,dS_imag,deltaS_real,deltaS_imag,dicestar_real,dicestar_imag,STAT=ios)
  deallocate(dil_slice,dlambda_slice,deltalambda_slice,STAT=ios)

  convflag = 0

!===========================================================
!                       CALCULATIONS
!___________________________________________________________

! Initialize topography (STEP 1)

do n = 1,nsteps
   topoxy(:,:,n) = initTopo(:,:) ! (eq. 48)
enddo


! Calculate ocean function from initial topography as first guess
do j = 1,2*nglv
   do i = 1,nglv
      if (initTopo(i,j) < epsilon(0.0)) then
         cxy(i,j,:) = 1
      else
         cxy(i,j,:) = 0
      endif
   enddo
enddo

! Decompose initial topography (STEP 2) (used to check convergence of outer
! loop)
call spat2spec(topoxy(:,:,1), t0lm, spheredat)

!-----------------------------------------------------------
!>>>>>>>>>>>>>>>>>>> Start of outer loop <<<<<<<<<<<<<<<<<<<
!-----------------------------------------------------------
do nouter = 1,numOuter ! Outer loop

   write(*,'(A,I2)') 'Outer loop iteration ', nouter
   call system_clock(counti1) ! Time for each outer loop


   !-----------------------------------------------------------
   !>>>>>>>>>>>>>>>>>>> Start of time loop <<<<<<<<<<<<<<<<<<<<
   !-----------------------------------------------------------
   do n = first_tstep,nsteps ! Time iteration --AC
      write(*,'(A,I4,A,EN15.4E2,A)') '  Timestep', n, ', ', times(n), ' years'
     
      ! Calculate icestar (STEP 3) (eq.43)
      if (checkMarine) then
         do j = 1,2*nglv
            do i = 1,nglv
               if (topoxy(i,j,n) > epsilon(0.0)) then
               ! If not marine...
                  icestarxy(i,j) = icexy(i,j,n)
               elseif (icexy(i,j,n) > (abs(topoxy(i,j,n)) * rhow / rhoi)) then
               !...else if marine, but thick enough to be grounded
                  icestarxy(i,j) = icexy(i,j,n)
               else
               !...if floating ice
                  icestarxy(i,j) = 0
               endif
            enddo
         enddo
      else ! If not checking for floating ice
         icestarxy(:,:) = icexy(:,:,n)
      endif

      ! Calculate beta (STEP 3) (eq. 44)

      do j = 1,2*nglv
         do i = 1,nglv
            if (icestarxy(i,j) < epsilon(0.0)) then
               beta(i,j) = 1
            else
               beta(i,j) = 0
            endif
         enddo
      enddo

     call spat2spec(icestarxy, icestarlm, spheredat) ! Decompose ice field

      if (n == 1) then
         dicestarlm(:,:) = 0.0          ! No change at first timestep
         icestar0(:,:) = icestarlm(:,:) ! Save initial to calculate total change
      else

         dicestarlm(:,:) = icestarlm(:,:) - oldicestarlm(:,:) ! Incremental change
      endif
      oldicestarlm(:,:) = icestarlm(:,:)                   ! Save to calculate increment on next time step
      deltaicestarlm(:,:) = icestarlm(:,:) - icestar0(:,:) ! Total change since time 0
      dicestar(:,:,n) = dicestarlm(:,:)         ! Save into big matrix (each slice for each time step)
      deltaicestar(:,:,n) = deltaicestarlm(:,:) ! Save into big matrix (each slice for each time step)

      ! Calculate cstar (STEP 4)
      if (numOuter == 1 .and. n > 1) then
      ! If running with only 1 outer loop, use the converged ocean function from
      ! the previous timestep (if available)
      ! as an initial guess to the current ocean function
         cxy(:,:,n) = cxy_guess(:,:)
      endif

      ! Calculate Ocean*Beta function at the current timestep
      cstarxy(:,:) = cxy(:,:,n) * beta(:,:) ! (eq. 65)

      if (n == 1) then
         cstar0 = cstarxy ! Value at time0, used in tO calculation below for shoreline migration
         call spat2spec(cstarxy, cstarlm, spheredat) ! Decompose cstar
         oldcstarlm(:,:) = cstarlm(:,:) ! So that there will be no change at time 0
      elseif (n==first_tstep) then ! so that oldcstarlm value we load in from preMWP1A
                                   ! isn't overwritten with zeros --AC
         call spat2spec(cstarxy, cstarlm, spheredat) ! Decompose cstar
      else
         oldcstarlm(:,:) = cstarlm(:,:) ! Save previous to use in eq. 79
         call spat2spec(cstarxy, cstarlm, spheredat) ! Decompose cstar
      endif

      ! Calculate tO (STEP 4), topography correction term
      tOxy(:,:) = topoxy(:,:,1) * (cstarxy(:,:) - cstar0(:,:)) ! (eq. 70)

      call spat2spec(tOxy, tOlm, spheredat) ! Decompose tO

      ! Calculate initial dS guess for first inner-iteration (STEP 5)
      if (numOuter == 1) then
      ! If running only 1 outer loop...

         if (n == 1) then ! If this is the first time step there has been no previous ocean loading changes
            dSlm(:,:) = (0.0,0.0)
         else             ! guess ocean loading change to be eustatic
            dSlm(:,:) =  (cstarlm(:,:) / cstarlm(0,0)) * ((-rhoi / rhow) * dicestarlm(0,0))
         endif

         elseif (numOuter > 1 .and. nouter == 1) then
         !...elseif it is the first outer loop of multiple outer loops

         if (n == 1) then  ! If this is the first time step

            dSlm(:,:) = (0.0,0.0)                      ! No difference on first timestep
            tTxy(:,:) = initTopo(:,:) * cstarxy(:,:)   ! (eq. 80) - tT is the script T in Kendall
            call spat2spec(tTxy, tTlm, spheredat)      ! Decompose tT
            oldtTlm(:,:) = tTlm(:,:)                   ! Save to calculate difference on next timestep

         else

            tTxy(:,:) = initTopo(:,:) * cstarxy(:,:)   ! (eq. 80) - tT is the script T in Kendall
            call spat2spec(tTxy, tTlm, spheredat)      ! Decompose tT
            oldtTlm(:,:) = tTlm(:,:)                   ! Save to calculate difference on next timestep
            dSlm(:,:) = (oldcstarlm(:,:) / oldcstarlm(0,0)) &
                         * ((-rhoi / rhow) * dicestarlm(0,0) + tTlm(0,0) - oldtTlm(0,0)) &
                         - (tTlm(:,:) - oldtTlm(:,:))  ! (eq. 79)

         endif

      else
      !...elseif it is the second outer loop and onwards
         dSlm(:,:) = dS(:,:,n) ! Use converged load from last outer loop (eq.82)

      endif

      !-----------------------------------------------------------
      !>>>>>>>>>>>>>>>>>>> Start of inner loop <<<<<<<<<<<<<<<<<<<
      !-----------------------------------------------------------
      ninner = 1
      
do ! Inner loop
   !-----\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/    Rotation  \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/----!

     if (tpw) then ! If incorporating rotation
            ! MOI perturbations (Matsuyama et al, 2006) (only need (1,3),(2,3),
            ! and (3,3))
            il(3,3) = (-8 * pi * radius**4 / (3 * sqrt(5.0))) * real(rhow * deltaS(2,0,n) + rhoi * deltaicestar(2,0,n))
            il(1,3) = (8 * pi * radius**4 / (sqrt(30.0))) * real(rhow * deltaS(2,1,n) + rhoi * deltaicestar(2,1,n))
            il(2,3) = (-8 * pi * radius**4 / (sqrt(30.0))) &
                     & * real(-1.0 * ii * (rhow * deltaS(2,1,n) + rhoi * deltaicestar(2,1,n)))

            if (n == 1) then
               dil(:,:,n) = il(:,:)
            else
               dil(:,:,n) = il(:,:) - oldil(:,:)
            endif

            ! Calculate m (rotation vector perturbation)
            !  - From Jerry's written notes, also Mitrovica, Wahr, Matsuyama,
            !  and Paulson (2005) eq. 2
            sum_il(:,:) = 0.0
            sum_m(:) = 0.0
            do nn = 1,n-1 ! Sum over all previous timesteps
               betatt = 0.0
               betattprime = 0.0
               do k = 1,nmod(2) ! Sum over k=1,K
                  betatt = betatt + (rprime(k,2) / s(k,2)) &
                           & * ( 1.0 - exp(-1.0 * s(k,2) * (times(n) - times(nn)) / 1000.0) )
                  betattprime = betattprime + &
                              & (rprimeT(k,2) / s(k,2)) * ( 1.0 - exp(-1.0 * s(k,2) * (times(n) - times(nn)) / 1000.0) )
               enddo
               sum_il(:,:) = sum_il(:,:) + dil(:,:,nn) * betatt
               sum_m(:) = sum_m(:) + dm(:,nn) * betattprime
            enddo

            do i = 1,2
               mm(i) = (1 / (1 - kTE(2) / kf)) * &
                     & ( (1 / (moiC - moiA)) * ((1 + kE(2)) * il(i,3) + sum_il(i,3)) + (1 / kf) * sum_m(i) )
            enddo

            mm(3) = (1.0 / moiC) * ( (1 + kE(2)) * il(3,3) + sum_il(3,3) )

            if (n == 1) then
               dm(:,n) = mm(:)
            else
               dm(:,n) = mm(:) - oldm(:)
            endif

            ! Calculate lambda (rotational driving) from m (Milne and Mitrovica
            ! 1998)
            lambda(0,0) = (radius**2 * omega**2 / 3.0) * ((mm(1)**2 + mm(2)**2 + mm(3)**2) + 2.0 * mm(3))
            lambda(2,0) = (radius**2 * omega**2 / (6.0 * sqrt(5.0))) &
                        & * (mm(1)**2 + mm(2)**2 - 2.0 * mm(3)**2 - 4.0 * mm(3))
            lambda(2,1) = (radius**2 * omega**2 / sqrt(30.0)) * ((mm(1) - ii * mm(2)) * (1.0 + mm(3)))
            lambda(2,2) = (radius**2 * omega**2 / sqrt(5.0 * 24.0)) * (mm(2)**2 - mm(1)**2 + 2.0 * ii * mm(1) * mm(2))


            if (n == 1) then
              lambda0(:,:) = lambda(:,:)
              dlambda(:,:,n) = (0.0,0.0)
              deltalambda(:,:,n) = (0.0,0.0)
            else
              dlambda(:,:,n) = lambda(:,:) - oldlambda(:,:)
              deltalambda(:,:,n) = lambda(:,:) - lambda0(:,:)
            endif

            ! Calculate effect on sea level (Love numbers) (Kendall)
            dsl_rot(:,:) = (0.0,0.0)
            if (n /= 1) then
               ekhTE = 1 + kTE(2) - hTE(2)
               do nn = 1,n-1 ! Sum over all previous timesteps
                  lovebetatt(nn) = 0.0
                  do k = 1,nmod(2) ! Sum over k=1,K
                     lovebetatt(nn) = lovebetatt(nn) + ((rprimeT(k,2) - rT(k,2)) / s(k,2)) &
                                      * ( 1 - exp(-1.0 * s(k,2) * (times(n) - times(nn)) / 1000.0) ) ! (eq. B27)
                  enddo
               enddo
               do m = 0,2
                  viscoustt = (0.0,0.0)
                  do nn = 1,n-1 ! Sum the loads over all previous timesteps to get the sum on the 4th line of eq. B28
                     viscoustt = viscoustt + lovebetatt(nn)*dlambda(2,m,nn)
                  enddo
                  dsl_rot(2,m) = (ekhTE * (deltalambda(2,m,n-1) + dlambda(2,m,n)) + viscoustt) / gacc ! (eq. B28/B25)
               enddo
            endif

         if (nouter == numOuter .and. calcRG) then ! For R calculations
               rr_rot(:,:) = (0.0,0.0)
               if (n /= 1) then
                  do nn = 1,n-1 ! Sum over all previous timesteps
                     lovebetattrr(nn) = 0.0
                     do k = 1,nmod(2) ! Sum over k=1,K
                        lovebetattrr(nn) = lovebetattrr(nn) &
                                          + ((rT(k,2)) / s(k,2)) &
                                          * (1 - exp(-1.0 * s(k,2) * (times(n) - times(nn)) / 1000.0)) ! (eq. B27)
                     enddo
                  enddo
                  do m = 0,2
                     viscousttrr = (0.0,0.0)
                     do nn = 1,n-1 ! Sum the loads over all previous timesteps to get the sum on the 4th line of eq. B28
                        viscousttrr = viscousttrr + lovebetattrr(nn) * dlambda(2,m,nn)
                     enddo
                     rr_rot(2,m) = (hTE(2) * (deltalambda(2,m,n-1) + dlambda(2,m,n)) + viscousttrr) / gacc
                  enddo
               endif
         endif

     else ! If not incorporating rotation
            dsl_rot(:,:) = (0.0,0.0)
            rr_rot(:,:) = (0.0,0.0)
     endif ! TPW


         !-----/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\    End Rotation /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\----!

         ! Calculate beta (eq. B14) - viscous response factor
         do l = 1,norder
            do nn = 1,n-1 ! Sum over all previous timesteps
               lovebeta(nn,l) = 0.0
               do k = 1,nmod(l) ! Sum over k=1,K
                  lovebeta(nn,l) = lovebeta(nn,l) + ((rprime(k,l) - r(k,l)) / s(k,l)) &
                                  * (1 - exp(-1.0 * s(k,l) * (times(n) - times(nn)) / 1000.0)) ! (eq. B14)
               enddo
            enddo
         enddo

         if (nouter == numOuter .and. calcRG) then ! For R calculations
            do l = 1,norder
               do nn = 1,n-1 ! Sum over all previous timesteps
                  lovebetarr(nn,l) = 0.0
                  do k = 1,nmod(l) ! Sum over k=1,K (modes)
                     lovebetarr(nn,l) = lovebetarr(nn,l) &
                                    & + ((r(k,l)) / s(k,l)) * (1 - exp(-1.0 * s(k,l) * (times(n) - times(nn)) / 1000.0))
                                    ! (eq. B14)
                  enddo
               enddo
            enddo
         endif

        ! Compute viscous response outside inner loop, since it doesn't depend on current timestep
         viscous(:,:) = (0.0,0.0)
         do l = 1,norder
            do m = 0,l
               do nn = 1,n-1 ! Sum the loads over all previous timesteps to get the sum on the second line of eq. B18
                  viscous(l,m) = viscous(l,m) + lovebeta(nn,l) * (rhoi * dicestar(l,m,nn) + rhow * dS(l,m,nn))
               enddo
            enddo
         enddo


         ! Calculate change in sea level from change in load (STEP 6)
         if (n == 1) then
            dsllm(:,:) = (0.0,0.0) ! No change on first timestep
         else
            do l = 1,norder
               ttl = (4.0 * pi * (radius**3)) / (mass * (2 * real(l) + 1.0)) !(eq. B16)
               ekhl = 1 + ke(l) - he(l) ! (eq. B13)
               do m = 0,l
                  ! Total  =
                  ! _____________________________elastic_____________________________ + _viscous_
                  dsllm(l,m) = ttl * ekhl * (rhoi * deltaicestar(l,m,n) + rhow * deltaS(l,m,n-1) + rhow * dSlm(l,m)) &
                                + ttl * viscous(l,m) ! (eq. B18)
               enddo
            enddo
         endif

         ! Add rotational effects (calculated above)
         dsllm(0:2,0:2) = dsllm(0:2,0:2) + dsl_rot(0:2,0:2)

         ! Convert dSL (total spatially heterogeneous change since time0) to spatial domain
         call spec2spat(dslxy, dsllm, spheredat)

         ! Compute r0 (STEP 7)
         rOxy(:,:) = dslxy(:,:) * cstarxy(:,:) ! (eq. 68)
         call spat2spec(rOxy, rOlm, spheredat) ! Decompose r0

         ! Compute conservation term (STEP 8)
         conserv = real((1 / cstarlm(0,0)) * (-1.0 * (rhoi / rhow) * deltaicestar(0,0,n) - rOlm(0,0) + tOlm(0,0))) ! (eq. 78)

         ! Compute change in ocean load again (STEP 9)
         olddSlm = dSlm ! Save previous-iterate load to check convergence
         if (n == 1) then
            dSlm(:,:) = (0.0,0.0)     ! No change at first timestep
         else
            dSlm(:,:) = -1.0 * deltaS(:,:,n-1) + rOlm(:,:) + conserv * cstarlm(:,:) - tOlm(:,:) ! (eq. 73)
         endif

         dS(:,:,n) = dSlm(:,:) ! Save converged ocean load into big matrix

         ! Save total load changes, used in Love number calculation
         if (n == 1) then
            deltaS(:,:,n) = dSlm(:,:)
         else
            deltaS(:,:,n) = deltaS(:,:,n-1) + dSlm(:,:)
         endif

        ! Update the total sea level change (STEP 10) (eq. 85)
         deltasllm(:,:) = dsllm(:,:) ! Spatially heterogeneous component
         deltasllm(0,0) = deltasllm(0,0) + conserv ! Add uniform conservation term to (0,0)
         call spec2spat(deltaslxy, deltasllm, spheredat) ! Synthesize deltasl

         ! Save into big array
         sl(:,:,n) = deltaslxy(:,:) !total sea-level change from the beginning to a current timestep
         ! Calculate convergence criterion for inner loop
         if ( abs(sum(abs(dSlm)) - sum(abs(olddSlm))) < epsilon(0.0) .and. abs(sum(abs(olddSlm))) < epsilon(0.0)) then
            xi = 0 ! Otherwise xi = 0 / 0 = NaN at the first loop of the first timestep.
         elseif (abs(sum(abs(olddSlm))) < epsilon(0.0)) then
            xi = abs( (sum(abs(dSlm)) - sum(abs(olddSlm))) / (epsilon(0.0) * 10)) ! Avoid dividing by 0
         else
            xi = abs( (sum(abs(dSlm)) - sum(abs(olddSlm))) / sum(abs(olddSlm)) )! (eq. 83)
         endif


         if (xi <= epsilon1) then
         ! if converged

            if (numOuter == 1) then ! If running only one outer loop...
               cxy_guess(:,:) = cxy(:,:,n) ! Save the ocean function for the next timetesp
               exit
            else                    !...else if running multiple outer loops
               exit
            endif

         endif

        if (xi > epsilon1 .AND. numouter == 1) then
         ! If not converged and numouter==1, update the following variables
         ! within the inner loop.

            ! Refine guess to the topography fields
            topoxy(:,:,n) = initTopo(:,:) - deltaslxy(:,:) ! (eq. 39)

            ! Update ocean function
            do j = 1,2*nglv
               do i = 1,nglv
                  if (topoxy(i,j,n) >= 0.0) then
                     cxy(i,j,n) = 0.0
                  else
                     cxy(i,j,n) = 1.0
                  endif
               enddo
            enddo


         ! Refine guess to ocean function with marine ice check
            cstarxy(:,:) = cxy(:,:,n) * beta(:,:) ! (eq. 65)
            call spat2spec(cstarxy, cstarlm, spheredat) ! Decompose cstar

            ! Refine guess to the topography correction fields
            tOxy(:,:) = initTopo * (cstarxy(:,:) - cstar0(:,:)) ! (eq. 70)
            call spat2spec(tOxy, tOlm, spheredat) ! Decompose tO

        endif


            !previously capped at 9999
         if (ninner == 50) then ! If no convergence after a huge number of iterations
            write(*,*)
            write(*,'(A,I5,A)') 'WARNING: The inner loop failed to converge after the limit of ', ninner, ' iterations.'
            write(*,'(A,ES15.3,A)') '         The variable xi finished with a value of ', xi, ', resulting from '
            write(*,'(A,ES15.3)') '         sum(abs(olddSlm)) = ',sum(abs(olddSlm))
            write(*,'(A,ES15.3)') '            sum(abs(dSlm)) = ',sum(abs(dSlm))
            write(*,*)
            write(*,'(A)') '!!!---- Program sl_model will now be terminated.----!!!'
            !call abort ! Terminate program
            !exit ! DEBUG line: Continue inner loop despite non-convergence.
            !Normal operation: Enable the 2 lines above.
            convflag=1


         else
            convflag=0
         endif

         ninner = ninner + 1


enddo ! End inner loop
      !-----------------------------------------------------------
      !<<<<<<<<<<<<<<<<<<<< End of inner loop >>>>>>>>>>>>>>>>>>>>
      !-----------------------------------------------------------

      write(*,'(A,I4,A)') '  ', ninner, ' inner-loop iterations'

      ! Assign rotation-related quantities for the next time step (for the
      ! calculation of various incremental values)
      if (tpw) then
         oldil(:,:) = il(:,:)
         oldm(:) = mm(:)
         oldlambda(:,:) = lambda(:,:)
      endif

      if (nouter == numOuter .and. calcRG) then ! For R calculations
         if (n == 1) then
            rrlm(:,:) = (0.0,0.0)! No change on first timestep
         else
            do l = 1,norder
               ttl = (4.0 * pi * (radius**3)) / (mass * (2.0 * real(l) + 1.0)) !(eq. B16)
               do m = 0,l
                  ! Viscous response
                  viscousrr = (0.0,0.0)
                  do nn = 1,n-1 ! Sum the loads over all previous timesteps to get the sum on the second line of eq. B18
                     viscousrr = viscousrr + lovebetarr(nn,l) * (rhoi *dicestar(l,m,nn) + rhow * dS(l,m,nn))
                  enddo
                  rrlm(l,m) = ttl * he(l) * (rhoi * deltaicestar(l,m,n) + rhow * deltaS(l,m,n-1) + rhow * dSlm(l,m)) &
                           & + ttl * viscousrr
               enddo
            enddo
         endif
         rrlm(0:2,0:2) = rrlm(0:2,0:2) + rr_rot(0:2,0:2)
         call spec2spat(rrxy, rrlm, spheredat)
         rr(:,:,n) = rrxy(:,:)
      endif

   enddo ! End time loop, n = 1,nsteps

   !-----------------------------------------------------------
   !<<<<<<<<<<<<<<<<<<<< End of time loop >>>>>>>>>>>>>>>>>>>>>
   !-----------------------------------------------------------

   ! Update topography fields, ocean functions (STEP 11)
   if (numOuter == 1) then
   ! If running only one outer loop...

      do n = 1,nsteps
         topoxy(:,:,n) = initTopo(:,:) - sl(:,:,n) ! (eq. 12) initTopo is the known initial topography
      enddo

   else
   !...elseif running multiple outer loops (i.e. for an ice-age simulation)

      oldt0lm = t0lm ! Save old initial topography to check convergence
      do n = 1,nsteps
         topoxy(:,:,n) = initTopo(:,:) + sl(:,:,nsteps) - sl(:,:,n) ! (eq. 39) initTopo is the present topography
         ! Update ocean function based on the new topography
         do j = 1,2*nglv
            do i = 1,nglv
               if (topoxy(i,j,n) >= 0.0) then
                  cxy(i,j,n) = 0.0
               else
                  cxy(i,j,n) = 1.0
               endif
            enddo
         enddo
      enddo

   endif


   ! Decompose initial topography to check convergence (STEP 12)
   call spat2spec(topoxy(:,:,1), t0lm, spheredat)

   ! Ice volume at each time step
   if (iceVolume) then
      glw_matrix = spread(spheredat%weights,2,2*nglv)/(2*nglv) !dividing each element by 2*nglv to keep the sum to be 2
      do n = 1,nsteps
         ice_volume(n) = (2*pi*radius**2)*sum(glw_matrix(:,:)*icexy(:,:,n))
      enddo
   endif

   ! call spharmt_destroy to deallocate memory from spharm_init
   call spharmt_destroy(spheredat)

   call system_clock(countf1) ! Time for each loop
   write(*,'(A,F6.2,A)') '  Loop time ', real(countf1 - counti1) / real(countrate), ' seconds'

   ! Calculate convergence criterion (if using convergence check for outer loop,
   ! defined to be when numOuter == 99.
   !  This means that convergence is assumed to happen in less than 100
   !  iterations. Modify if necessary.)
   if (numOuter == 99) then
      zeta = abs( (sum(abs(t0lm)) - sum(abs(oldt0lm))) / sum(abs(oldt0lm)) ) !(eq. 86)
      if (zeta < epsilon2) then ! If converged
         exit ! Get out of outer loop
      endif
   endif


enddo ! End outer loop, nouter = 1,numOuter
!-----------------------------------------------------------
!<<<<<<<<<<<<<<<<<<<< End of outer loop >>>>>>>>>>>>>>>>>>>>
!-----------------------------------------------------------

! If set to run until topography convergence (assumed to happen in less than 100
! iterations)
if (numOuter == 99) then
   write(*,*)
   write(*,*) '  **************'
   write(*,*) '  * Converged! *'
   write(*,*) '  **************'
   write(*,*)
   write(*,'(I3,A)') nouter, ' topography iterations'
   write(*,*)
endif


!===========================================================
!                          OUTPUT
!___________________________________________________________

write(*,'(A)') 'Writing output files...'

!site- and time- specific MWP1A sea level change:

!order is: 'Tahiti', 'Barbados', 'Sunda Shelf',
!          'HYD' (Great Barrier Reef),'NOG'(Great Barrier Reef),
!          'NWS' (NW Scotland)

if (convflag.eq.0) then !i.e. the sl run was successful

!sea level change timeseries relative to the 
!first checkpoint timestep at each site:

! 1. Tahiti:
predicted_data(1,:)=sl(tahiti_idx_lat,tahiti_idx_lon,first_tstep:nsteps) - sl(tahiti_idx_lat,tahiti_idx_lon,first_tstep)

! 2. Barbados:
predicted_data(2,:)=sl(barbados_idx_lat,barbados_idx_lon,first_tstep:nsteps) - sl(barbados_idx_lat,barbados_idx_lon,first_tstep)

! 3. Sunda:
predicted_data(3,:)=sl(sunda_idx_lat,sunda_idx_lon,first_tstep:nsteps) - sl(sunda_idx_lat,sunda_idx_lon,first_tstep)

! 4. HYD:
predicted_data(4,:)=sl(hyd_idx_lat,hyd_idx_lon,first_tstep:nsteps) - sl(hyd_idx_lat,hyd_idx_lon,first_tstep)

! 5. NOG:
predicted_data(5,:)=sl(nog_idx_lat,nog_idx_lon,first_tstep:nsteps) - sl(nog_idx_lat,nog_idx_lon,first_tstep)

! 6. NW Scotland:
predicted_data(6,:)=sl(nws_idx_lat,nws_idx_lon,first_tstep:nsteps) - sl(nws_idx_lat,nws_idx_lon,first_tstep)

!check for sea level oscillation by ensuring RSL has a negative slope on average over MWP-1A window

    if (predicted_data(6,37)-predicted_data(6,1).ge.0.0) then

    predicted_data(:,:)=20000.0

    endif


else

predicted_data(:,:)=10000.0

endif

  ! deallocate memory associated with large local arrays:
  deallocate(icexy,r,s,rprime,ke,he,rprimeT,rT, STAT=ios)
  deallocate(kTE,hTE,initTopo,sl,topoxy, STAT=ios)
  deallocate(cxy,deltaslxy,dslxy,icestarxy,glw_matrix,STAT=ios)
  deallocate(cstar0,cstarxy,beta,cxy_guess,STAT=ios)
  deallocate(tOxy,rOxy,tTxy,cstarlm,oldcstarlm,STAT=ios)
  deallocate(tOlm,rOlm,dSlm,olddSlm,STAT=ios)
  deallocate(icestarlm,dicestarlm,deltaicestarlm,STAT=ios)
  deallocate(oldicestarlm,icestar0,t0lm,oldt0lm,STAT=ios)
  deallocate(tTlm,oldtTlm,dsllm,deltasllm,STAT=ios)
  deallocate(dicestar,dS,deltaicestar,deltaS,STAT=ios)
  deallocate(lovebeta,viscous,ice_volume,rrxy,rr,rrlm,STAT=ios)
  deallocate(lovebetarr,lambda,oldlambda,lambda0,STAT=ios)
  deallocate(dlambda,deltalambda,mm,sum_m,oldm,dm,STAT=ios)
  deallocate(il,oldil,sum_il,dil,dsl_rot,rr_rot,STAT=ios)
  deallocate(lovebetatt,lovebetattrr,STAT=ios)
  deallocate(elast,asymv,telast,tasymv,STAT=ios)
  deallocate(resk,resl,resh,tresk,tresl,tresh,STAT=ios)

call system_clock(countf) ! Total time
write(*,'(A,F7.2,A)') 'Done! Total time ', real(countf - counti) / real(countrate), ' seconds'

return
end subroutine sl_model

end module sea_level

