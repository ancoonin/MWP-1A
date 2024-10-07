! ====================================================================== !
! mwp_na_driver.f90                                                      !
! Coonin 2023                                                            !
!                                                                        !
! Driver program that calls the Neighbourhood algorithm routines of      !
! Sambridge et al. (1999) to search a parameter space.                   !
! Currently set up for Sea Level Fingerprint of MWP-1A                   !
!                                                                        !
! Uses various modules : sea_level.f90, user_specs_mod.f90,              !
!                        planets_mod.f90, constants_mod.f90              !
! ====================================================================== !

program mwp_na_driver

  use spharmt
  use user_specs_mod
  use constants_mod
  use planets_mod
  use sea_level

  implicit none
  include 'mwp_na_param.inc'


  call na

stop
end program

! ====================================================================== !
! subroutine user_init()                                                 !
! Coonin 2023                                                            !
!                                                                        !
! inputs: nd_max (from mwp_na_param.inc)                                 !
! outputs: ranges, scales, nd                                            !
! ====================================================================== !

subroutine user_init(nd,ranges,scales)

  ! Declare Variables
  implicit none
  include 'mwp_na_param.inc'
  include 'prem.l75.umvm7.lmvm7.inc'

  real, dimension(2,nd_max) :: ranges
  real, dimension(nd_max+1) :: scales
  integer :: nd
  ! HL: common blocks are part of the declaration and go first.

  !! HL: Because of "implicit none" *everything* has to be declared
  !!     So while they aren't in the original NA scripts, let's do it here
  !!     I googled this:

  ! "The rule is that undeclared variables which have the first letter I, J, K,
  ! L, M or N are considered INTEGER-s, and undeclared variables which start in
  ! A through H and O through Z are considered REAL-s."
  ! everytime we see the common block, we'll copy and paste this;

  integer ::  lu_na,lu_out,lu_sum,lu_det,lu_sob,lu_dis,lu_nad
  real ::  verbose,debug,timing,summary
  integer :: iproc,nproc,lroot
  real, dimension(ndata) :: observed_data,predicted_data,SLIP_error
  integer :: lu_mod

! Info and Logical unit common
! blocks used by NA routines:

  common /NA_IO/lu_na,lu_out,lu_sum,lu_det,lu_sob,lu_dis,lu_nad,verbose,&
         debug,timing,summary

  common /NAMPI/iproc,nproc,lroot

!  common /mwp_com/observed_data, predicted_data, SLIP_error,lu_mod

  ! Define scales (non-dimensionalizations) and ranges
  ! for each free parameter --AC

  scales(1)=-1 ! scales(1) = 0 this makes all scales = 1
              ! scales(1) = -1 this makes the scale for each parameter the difference between the min max range (prior variance)
              ! anything else, still need to specify each scale factor

 !-------------------------------------------------------------------------             
!run should converge on the "true model" --AC

 !e.s.l. contributions from each ice sheet:
 !scales(2:3)=20.0 ! scaling by fraction of approx average MWP1A sea level rise
 ranges(1,1:2)=0.0 !min possible contribution of e.s.l. from each ice sheet
 
  ! for background maxwell earth model: prem.l75.umvm7.lmvm7
  ! read in from prem.inc file in include directory:

 !ranges(2,1)=LISmaxesl !upper bound for LIS
 !ranges(2,2)= GISmaxesl  !upper bound for GIS
 ranges(2,1)= EISmaxesl  !upper bound for EIS
 ranges(2,2)= WAISmaxesl !upper bound for WAIS
 !ranges(2,5)= WLmaxesl !upper bound for WL

 ! Note that the max possible contributions of e.s.l. from each ice sheet vary by earth model 

 !-------------------------------------------------------------------------             

 !onset of melting for each ice sheet:

 !scales(4:5)=650.0 ! range of potential start times for individual ice sheets to melt in years (14.675 ka to 14.025 ka)
 ranges(1,3:4)=-14675.0  ! earliest possible onset is just after 14.675 ka
 ranges(2,3:4)=-14025.0  ! latest possible onset

 !-------------------------------------------------------------------------             

 !duration of melting for each ice sheet:
 
 !scales(6:7)=500.0 ! range of potential durations of ice melting < 500 years)
 ranges(1,5:6)=25.0 ! if there is nonzero melting, the shortest duration it can occur over is 1 timestep
 ranges(2,5:6)=500.0 ! max duration of melting for each ice sheet

  ! Set up logical units
  lu_out = 6
  lu_mod = 12

  nd=nd_max
return
end subroutine user_init

! ====================================================================== !
! subroutine forward()                                                   !
! Coonin 2023                                                            !
!                                                                        !
! inputs: iiter,imod,model                                               !
! outputs: lppd, model_overwrite                                         !
! ====================================================================== !

subroutine forward(nd,model,lppd,model_overwrite)

  use spharmt
  use user_specs_mod
  use constants_mod
  use planets_mod
  use sea_level

! declare Variables
  implicit none
  include 'mwp_na_param.inc'
  include 'mwp_param.inc' !contains values for f_ice and SLIP sea level data
  include 'prem.l75.umvm7.lmvm7.inc'


  ! input/output
  integer, intent(in) :: nd

  real, intent(in), dimension(nd) :: model
  real,intent(out), dimension(nd) :: model_overwrite


  ! internal
  !character(len=255) :: fmt1,fmt2
  CHARACTER(LEN=255) :: planetModel,path2loveinputs,path2SLinputs,path2iceinputs
  real, dimension(2) :: per_tstep
  real, dimension(2) :: f_ice
  real, dimension(ndata,ntimesduringMWP) :: observed_data,predicted_data
  real, dimension(ndata,ntimesduringMWP) :: observed4misfit,predicted4misfit
  real, dimension(ndata)::SLIP_error
  real :: misfitval,aval
  real :: lppd
  integer :: i, j, k, ios
  logical :: THERE
  real, dimension(ntimesduringMWP+1) :: tsteps, tdiff    ! Timesteps in checkpoint sl subroutine(years)
  integer, dimension(2) :: start_ix, end_ix, nactivetimesteps
  real, dimension(2) :: start_time, end_time, dur
  !Info and Logical unit common block used by NA routines
  integer ::  lu_na,lu_out,lu_sum,lu_det,lu_sob,lu_dis,lu_nad
  real ::  verbose,debug,timing,summary
 ! integer :: iproc,nproc,lroot
  integer :: lu_mod

  !real, dimension(100) :: times
  real, allocatable :: times (:)                ! Timesteps of ice model (years)
  
  !real, dimension(nlat,nlon,ntimesduringMWP) :: icegrid, ice_load
  real, allocatable :: icegrid(:,:,:)
  real, allocatable :: LIS(:,:,:)
  real, allocatable :: GIS(:,:,:)
  real, allocatable :: EIS(:,:,:)
  real, allocatable :: WAIS(:,:,:)
  real, allocatable :: WL(:,:,:)


  !real, dimension(nlat,nlon) :: contribution
  real, allocatable :: contribution(:,:)
  
  !real, dimension(nlat,nlon) :: preMWP_icegrid, baseload
  real, allocatable :: PREMWP_icegrid(:,:)
  real, allocatable :: baseload(:,:)
        
  real, allocatable :: LIS_preMWP(:,:)
  real, allocatable :: GIS_preMWP(:,:)
  real, allocatable :: EIS_preMWP(:,:)
  real, allocatable :: WAIS_preMWP(:,:)
  real, allocatable :: WL_preMWP(:,:)

  !!AC testing:
  real,dimension(512,1024):: check
  character(LEN=2) :: tindx

common /NA_IO/lu_na,lu_out,lu_sum,lu_det,lu_sob,lu_dis,&
                 lu_nad,verbose,debug,timing,summary

!common /mwp_com/observed_data, predicted_data, SLIP_error,lu_mod

  ! begin code:

call read_config(planetModel,path2loveinputs,path2SLinputs,path2iceinputs) !--AC

   ! write(fmt1,'(I4)') nlon
   ! fmt2 = '('//trim(fmt1)//'E15.7)'

  allocate (icegrid(nlat,nlon,ntimesduringMWP))
  allocate(times(108))
  allocate(preMWP_icegrid(nlat,nlon),contribution(nlat,nlon),baseload(nlat,nlon))
  allocate(LIS_preMWP(nlat,nlon),GIS_preMWP(nlat,nlon),EIS_preMWP(nlat,nlon),WAIS_preMWP(nlat,nlon),WL_preMWP(nlat,nlon))
  allocate(LIS(nlat,nlon,ntimesduringMWP),GIS(nlat,nlon,ntimesduringMWP),EIS(nlat,nlon,ntimesduringMWP))
  allocate(WAIS(nlat,nlon,ntimesduringMWP),WL(nlat,nlon,ntimesduringMWP))

  !initializing arrays to avoid NaNs
  icegrid(:,:,:)=0.0
  LIS(:,:,:)=0.0
  GIS(:,:,:)=0.0
  EIS(:,:,:)=0.0
  WAIS(:,:,:)=0.0
  WL(:,:,:)=0.0


  !times -- see below
  preMWP_icegrid(:,:)=0.0
  contribution(:,:)=0.0
  baseload(:,:)=0.0
  LIS_preMWP(:,:)=0.0
  GIS_preMWP(:,:)=0.0
  EIS_preMWP(:,:)=0.0
  WAIS_preMWP(:,:)=0.0
  WL_preMWP(:,:)=0.0


  ! find nearest timesteps for start and end of melting based on the sampled onset and duration for each ice sheet

  times = (/-122000.0, -120000.0, -118000.0, -116000.0, -114000.0, -112000.0, &
   -110000.0, -108000.0, -106000.0, -104000.0, -102000.0, -100000.0, &
   -98000.0, -96000.0, -94000.0, -92000.0, -90000.0, -88000.0, &
   -86000.0, -84000.0, -82000.0, -80000.0, -78000.0, -76000.0, &
   -74000.0, -72000.0, -70000.0, -68000.0, -66000.0, -64000.0, &
   -62000.0, -60000.0, -58000.0, -56000.0, -54000.0, -52000.0, &
   -50000.0, -48000.0, -46000.0, -44000.0, -42000.0, -40000.0, &
   -38000.0, -36000.0, -34000.0, -32000.0, -31000.0, -30000.0, &
   -29000.0, -28000.0, -27000.0, -26000.0, -25000.0, -24000.0, &
   -23000.0, -22000.0, -21000.0, -20000.0, -19000.0, -18250.0, &
   -17500.0, -16750.0, -16250.0, -15750.0, -15500.0, -15250.0, &
   -15125.0, -15000.0, -14875.0, -14750.0, -14675.0, -14650.0, &
   -14625.0, -14600.0, -14575.0, -14550.0, -14525.0, -14500.0, &
   -14475.0, -14450.0, -14425.0, -14400.0, -14375.0, -14350.0, &
   -14325.0, -14300.0, -14275.0, -14250.0, -14225.0, -14200.0, &
   -14175.0, -14150.0, -14125.0, -14100.0, -14075.0, -14050.0, &
   -14025.0, -14000.0, -13975.0, -13950.0, -13925.0, -13900.0, &
   -13875.0, -13850.0, -13825.0, -13800.0, -13775.0, -13750.0/)

  tsteps(:)=times(first_tstep-1:size(times))

        do i=1,2
        !onset of melting
tdiff=ABS(tsteps-model(i+2))
        start_ix(i)=MINLOC(tdiff,1) !start_ix is an integer. 
  start_time(i)=tsteps(start_ix(i))

  
  !end of melting
  tdiff=ABS(tsteps-start_time(i)-model(i+4))
  end_ix(i)=MINLOC(tdiff,1)
  end_time(i)=tsteps(end_ix(i))

  !since we are only running the SL code to 13.75 ka, cap duration of melting
  
          if (end_time(i)>MAXVAL(times)) then
                end_time(i)=MAXVAL(times)
          endif 
  end_ix(i)=end_ix(i)-1 !we start the checkpoint at the timestep following 14.675 ka so this ensures melting ceases by the correct time
  !within the checkpoint script

  enddo

  do i=1,2
        if (model(i).lt.0.01) then !if the NA chooses e.s.l. less than 1cm
                                   !set to zero 
        model_overwrite(i)=0.0
        else 
        model_overwrite(i)=model(i)
        endif
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  dur(:) = end_time(:)-start_time(:)
  model_overwrite(3:4)=start_time(:)
  model_overwrite(5:6)=dur(:)

  write(*,*) 'adjusted model from NA for SL: ', model_overwrite(:)

  ! set up multiplication factor for 1 meter of eustatic melt:
  ! values are from mwp_param.inc !
  !f_ice(1)=NAm_esl ! North American Ice Sheet
  !f_ice(2)=GL_esl ! Greenland Ice Sheet
  f_ice(1)=EIS_esl ! Scandinavian Ice Sheet
  f_ice(2)=WA_esl ! West Antarctic Ice Sheet
  !f_ice(5)=EA_esl ! East Antarctic Ice Sheet

 ! load in pre-MWP ice grid (14.675 ka extent)
  open(unit=99,file=trim(path2iceinputs)//'ice_7g_14675ka.dat',form='formatted',&
       access='sequential',status='old')
        do i=1,nlat
        read(99,*) (preMWP_icegrid(i,j),j=1,nlon)
        enddo
  close(1)

  !now load individual 14.675 ka ice grids for each ice sheet

!  ! Laurentide Ice Sheet (LIS)
!  open(unit=97,file=trim(path2iceinputs)//'LIS_14.675_ka.dat',form='formatted',&
!     access='sequential',status='old')
!        do i=1,nlat
!        read(97,*) (LIS_preMWP(i,j),j=1,nlon)
!        enddo
!
!  close(97)


!! Greenland Ice Sheet (GIS)
!
!  open(unit=96,file=trim(path2iceinputs)//'GIS_14.675_ka.dat',form='formatted',&
!   access='sequential',status='old')
!        do i=1,nlats
!        read(96,*) (GIS_preMWP(i,j),j=1,nlon)
!        enddo
!  close(96)

 ! Eurasian IS)

  open(unit=95,file=trim(path2iceinputs)//'EIS_14.675_ka.dat',form='formatted',&
   access='sequential',status='old')
        do i=1,nlat
        read(95,*) (EIS_preMWP(i,j),j=1,nlon)
        enddo
  close(95)

! West Antarctic Ice Sheet (WAIS)

  open(unit=94,file=trim(path2iceinputs)//'WAIS_14.675_ka.dat',form='formatted',&
   access='sequential',status='old')
        do i=1,nlat
        read(94,*) (WAIS_preMWP(i,j),j=1,nlon)
        enddo
  close(94)

!! Wilkes Land, East Antarctica (W:IS)
!
!  open(unit=93,file=trim(path2iceinputs)//'WL_14.675_ka.dat',form='formatted',&
!   access='sequential',status='old')
!        do i=1,nlat
!!        read(93,*) (WL_preMWP(i,j),j=1,nlon)
!        enddo
!  close(93)
!
! add any ice that falls outside of the 5 ice masks when constructing each ice load file
! to ensure that this mass is not artificially lost during MWP1a: 

!baseload(:,:)=preMWP_icegrid(:,:) - LIS_preMWP(:,:) - GIS_preMWP(:,:) - EIS_preMWP(:,:) - WAIS_preMWP(:,:) - WL_preMWP(:,:)
baseload(:,:)=preMWP_icegrid(:,:) - EIS_preMWP(:,:) - WAIS_preMWP(:,:)

!EURASIAN ICE SHEET LOSS
i=1
nactivetimesteps(i)=NINT(dur(i)/25.0) !number of 25 year timesteps over which melting will occur
per_tstep(i)=model_overwrite(i)/real(nactivetimesteps(i)) ! m of e.s.l. to remove during each active melting timestep
contribution(:,:)=EIS_preMWP(:,:)*per_tstep(i)*f_ice(i);! ice thickness map to shave off during each active melting timestep 

          do j=1,ntimesduringMWP
                if (j.lt.start_ix(i)) then
                EIS(:,:,j)=EIS_preMWP(:,:) ! the ice sheet should not be melting yet
                else if (j.eq.start_ix(i)) then
                EIS(:,:,j)=EIS_preMWP(:,:) - contribution ! first active melting timestep
                else if ((j.gt.start_ix(i)) .AND.( j.LE.end_ix(i))) then
                EIS(:,:,j)=EIS(:,:,j-1) - contribution ! subsequent active timestep
                else if (j.gt.end_ix(i)) then
                EIS(:,:,j)=EIS(:,:,j-1) !no additional melting
                else
                write(*,*) 'somethings wrong with the melting time discretization for EIS'
                endif
          enddo

! next: WAIS loss
i=2
nactivetimesteps(i)=NINT(dur(i)/25.0) !number of 25 year timesteps over which melting will occur
per_tstep(i)=model_overwrite(i)/real(nactivetimesteps(i)) ! m of e.s.l. to remove during each active melting timestep
contribution(:,:)=WAIS_preMWP(:,:)*per_tstep(i)*f_ice(i);! ice thickness map to shave off during each active melting timestep 

          do j=1,ntimesduringMWP
                if (j.lt.start_ix(i)) then
                WAIS(:,:,j)=WAIS_preMWP(:,:) ! the ice sheet should not be melting yet
                else if (j.eq.start_ix(i)) then
                WAIS(:,:,j)=WAIS_preMWP(:,:) - contribution ! first active melting timestep
                else if ((j.gt.start_ix(i)) .AND.( j.LE.end_ix(i))) then
                WAIS(:,:,j)=WAIS(:,:,j-1) - contribution ! subsequent active timestep
                else if (j.gt.end_ix(i)) then
                WAIS(:,:,j)=WAIS(:,:,j-1) !no additional melting
                else
                write(*,*) 'somethings wrong with the melting time discretization for WAIS'
                endif
          enddo


          !sum contribution from each ice sheet at each time step to get total ice load for each time step
          
          do j=1,ntimesduringMWP
          icegrid(:,:,j)=EIS(:,:,j)+WAIS(:,:,j)+baseload(:,:)

          enddo


      deallocate(preMWP_icegrid,contribution,baseload,STAT=ios)
      deallocate(LIS,GIS,EIS,WAIS,WL,STAT=ios)
      deallocate(LIS_preMWP,GIS_preMWP,EIS_preMWP,WAIS_preMWP,WL_preMWP,STAT=ios)

!call sl subroutine
call sl_model(icegrid,ntimesduringMWP,times,predicted_data)

deallocate (icegrid,STAT=ios) !--AC
deallocate (times,STAT=ios) !--AC

!MWP1A magnitude of sea level change from SLIP data (located in mwp_param.inc)
observed_data(1,:)=dsl_Tahiti !Tahiti, Barbados, Sunda Shelf, HYD, NOG
observed_data(2,:)=dsl_Barbados
observed_data(3,:)=dsl_Sunda
observed_data(4,:)=dsl_HYD
observed_data(5,:)=dsl_NOG
observed_data(6,:)=dsl_NWS

!std in magnitude of MWP1A sea level change (m)

SLIP_error(1)=err_Tahiti
SLIP_error(2)=err_Barbados
SLIP_error(3)=err_Sunda
SLIP_error(4)=err_HYD
SLIP_error(5)=err_NOG
SLIP_error(6)=err_NWS


! calculating misfit
 misfitval=0.0

do i=1,ndata !considering all sites

        do j=1,ntimesduringMWP
        aval=(observed_data(i,j)-predicted_data(i,j))/SLIP_error(i)
        misfitval=misfitval+aval**2
        enddo
enddo

lppd=misfitval

write(*,*) 'misfit value for current model: ',lppd
!write(*,*) 'predicted_data for current model: ',predicted_data

! writing file with all predicted_data values for each forward model run:

       INQUIRE( FILE="predicted_data", EXIST=THERE ) 
         if ( THERE ) then 
          open(2, file="predicted_data", status="old", position="append", action="write")
         else
          open(2, file="predicted_data", status="new", action="write")
         endif
         do i=1,ndata
         write(2,'(38E15.7)') (predicted_data(i,j),j=1,ntimesduringMWP),lppd
         enddo
         close(2)


return
end subroutine forward

! ====================================================================== !
! subroutine writemodels()                                               !
! Coonin 2022                                                            !
!                                                                        !
! Use to write out models produced by Neighbourhood algorithm            !
!                                                                        !
! inputs: nd (number of dimensions of parameter space),                  !
!       ntot (number of models generated by NA),                         !
!       models (models generated by na),                                 !
!       misfit (array of model misfits (-log(ppd)),                      !
!       ns1 (initial sample size used by NA),                            !
!       ns2 (normal sample size used by NA),                             !
!       itmax (number of iterations)                                     !
! outputs:                                                               !
! ====================================================================== !

subroutine writemodels(nd,ntot,models,misfit,ns1,ns2,itmax, &
                       nh)!, header)

  include 'mwp_param.inc'
  include 'mwp_na_param.inc'

! NA variables and arrays
  real,dimension(nd,ntot) :: models !changed from real*8
  real :: misfit(ntot) 
  real :: mfitmin
  real :: mfitminc
  real :: mfitmean
  real, dimension(nd) :: rmodel,model_overwrite

  real :: lppd

  !character(len=255) :: header
  integer :: nd, ns, np, mopt, ns1, ns2, itmax, nh, ntot !--AC
  integer ::  lu_na,lu_out,lu_sum,lu_det,lu_sob,lu_dis,lu_nad
  real ::  verbose,debug,timing,summary
 ! integer :: iproc,nproc,lroot
  integer :: lu_mod, jj, lu_out2
  real, dimension(ndata) :: observed_data,predicted_data,SLIP_error

!Info and Logical unit common block used by NA routines

  common /NA_IO/lu_na,lu_out,lu_sum,lu_det,lu_sob,lu_dis,&
                lu_nad,verbose,debug,timing,summary

!  common /mwp_com/ observed_data, predicted_data, SLIP_error,lu_mod

!write out models at each iteration
  mfitmin = misfit(1)
  ns = ns1
  np = 0
  mopt = 1

!turn off writing to standard out by setting lu to zero
  !lu_out2 = lu_out
  lu_out2 = 0 ! RFI example had this setting --AC

!writing file to output model info --AC
!open(unit=lu_mod,file='mwp_output',form='formatted',&
!       access='sequential',status='unknown')
!
!  write(lu_mod,*) ns1,' Number of samples in starting pool'
!  write(lu_mod,*) ns2,' Number of new samples per iteration'
!  write(lu_mod,*) itmax,' Number of iterations'
!
!loop over iterations
  do it=1,itmax+1
   mfitminc = misfit(np+1)
   mfitmean = 0.0

    !find minimum and mean misfit
      do i=1,ns
         jj = np + i
         if(misfit(jj).lt.mfitmin)then
           mfitmin = misfit(jj)
           mopt = jj
         end if
         mfitminc = min(mfitminc,misfit(jj))
         mfitmean = mfitmean + misfit(jj)
      end do
      mfitmean = mfitmean/ns

    !write out all models from this iteration to mwp_output file.
      !write(lu_mod,801) it-1, mfitmin, mfitmean, mfitminc
!       write(lu_mod,*) 'it:',it-1, 'mfitmin:', mfitmin, 'mfitmean:', mfitmean, 'mfitminc:', mfitminc
!      do i=1,ns
!           jj = np + i
!           call display_model(lu_mod, i, models(1,jj), nd, misfit(jj) )
!
!      end do

    np = np + ns
    ns = ns2

    ! Write out summary of optimum model to SL model file and NA summary file

    call output_summary(lu_out2, lu_sum, it-1, models(1,mopt), &
                       np, mfitmin, mfitmean, mfitminc, mopt)
  end do

! Write out final model

 ! call display_final(lu_out, models(1,mopt), nd, mfitmin)
 ! call display_final(lu_mod, models(1,mopt), nd)!, mfitmin)
!close(lu_mod)

  call display_final(lu_sum, models(1,mopt), nd)!, mfitmin)

    do j=1,nd
       rmodel(j) = models(j,mopt) !--AC 
    end do

! repeat forward modelling for optimum model

 call forward(nd,rmodel,lppd,model_overwrite)

!write RFI component of header for NAD file (can leave blank and set nh to 0)
  nh = 0
!--AC commented out to test
  !801      format( 'iteration:',i5,',  misfit: min=',ES10.8, &
  !          ', mean=',ES10.8,', minc=',ES10.8 )

return
end subroutine writemodels

! ====================================================================== !
! subroutine display_model()                                             !
! Coonin 2022                                                            !
!                                                                        !
!  writes out the model contained in the array rmodel and misfit value   !
!  to the file assigned to logical unit lu                               !
!                                                                        !
!       Note: Calls no other routines                                    !
! ====================================================================== !

subroutine display_model(lu,imod,rmodel,moddim,misfitval)

  include 'mwp_param.inc'
  integer :: moddim
  real, dimension(moddim) :: rmodel
  real :: misfitval
  integer :: lu,imod
  !Write out model 'rmodel'
  !write(lu,801) imod, misfitval
  write(lu,fmt=*) 'model: ', imod,'misfit value: ', misfitval

  do j=1,moddim
  !write(lu,811) rmodel(j)
  write(lu,fmt=*) rmodel(j)
  end do

  !801   format( '  model:',i5,',   misfit value:',ES10.8 )
  !811   format(ES10.8)

return
end subroutine display_model

! ====================================================================== !
! subroutine output_summary()                                            !
! Coonin 2022                                                            !
!                                                                        !
!  writes out a summary of the current iteration of the Neighborhood     !
!  algorithm together with the model 'rmodel'                            !
!                                                                        !
!       Note: Calls no other routines                                    !
! ====================================================================== !
subroutine output_summary(lu_out, lu_sum, it, rmodel, &
            ntot, mfitmin, mfitmean, mfitminc, mopt)

  include 'mwp_param.inc'
  include 'mwp_na_param.inc'

  real,dimension(nd_max) :: rmodel
  real :: mfitmin, mfitmean, mfitminc
  logical :: lw
  integer :: lu_sum, lu_out, it,  ntot, mopt
!write out headers for files with summary of optimization performance

  lw = .true.
  if(lu_out.eq.0)lw = .false.
  !write(lu_sum,801)
  write(lu_sum,*) 'it: ', it, 'Nsampled: ',ntot
  write(lu_sum,*) 'Mfitmin: ',mfitmin,'Mfitmean: ',mfitmean,'Mfitmeanc: ',mfitminc, 'Mopt: ', mopt

  if(lw)write(lu_out,801)

  !write(lu_sum,811) it,ntot,mfitmin,mfitmean, &
  !                  mfitminc,mopt

  if(lw)write(lu_out,811) it,ntot,mfitmin,mfitmean, &
                    mfitminc,mopt

  do j=1,nd_max
  !write(lu_sum,812) rmodel(j)
  write(lu_sum,*) rmodel(j)
  if(lw)write(lu_out,812) rmodel(j)
  enddo

  write(lu_sum,*)
  if(lw)write(lu_out,*)

  801   format( 3x,'It',1x,'Nsampled',3x,'Mfitmin',&
             2x,'Mfitmean',1x,'Mfitmeanc',1x, &
                3x,'Mopt')
  811   format( i5,i9,3ES10.8,1x,i7)
  812   format(ES10.8)

return
end subroutine output_summary

! ====================================================================== !
! subroutine display_final()                                             !
! Coonin 2022                                                            !
!                                                                        !
! writes out the final model containedrmodel to the file                 !
! assigned logical unit lu                                               !
!                                                                        !
!       Note: Calls no other routines                                    !
! ====================================================================== !

subroutine display_final(lu, rmodel, moddim)
        !, misfitval)

   ! include 'mwp_param.inc'
   ! include 'mwp_na_param.inc'
    integer :: lu, moddim
    real,dimension(moddim) :: rmodel
    ! Write out final model

    write(lu,'(A)') '*** Final model ***'

    do j=1,moddim
      write(lu,*) rmodel(j)
    end do

    !811   format(ES10.8)

return
end subroutine display_final


!====================================================================================!

