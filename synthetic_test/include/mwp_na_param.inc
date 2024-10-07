!	Array sizes for neighbourhood algorithm routines
!
!	PARAMETER		MEANING
!                                       The NA routines use the
!                                       following include file to define
!                                       all parameters controlling
!                                       memory required by the arrays.
!
!                                       The following is a description
!                                       of each parameter in na_param.in
!                                       that you are allowed to change to suit
!                                       your application. They are
!                                       checked against input parameters
!                                       read in from file na.in and an
!                                       error is reported if any values
!                                       are too small.
!
!       PARAMETER               MEANING
!
!       nd_max                  Maximum number of variables, i.e. the
!                               dimension of the parameter space.
!
!       nit_max                 Maximum number of iterations
!
!       nsample_max             Maximum number of models generated
!                               in each sample
!
!       nh_max                  Maximum size of header of NAD output file
!                               (leave unchanged if you choose not to add
!                                any material to the NAD header. See manual)
!
!	maxseq			Maximum number of random sequences generated
!				by quasi random number generator. This
!				value can be set to 1 if the quasi
!				random number generator is not used.
!
!-----------------------------------------------------------------------
!
!				The following parameters are fixed
!				and should not be changed.
!
!	nmod_max		Maximum number of models to be generated
!			        (determined by nsample_max and nit_max)
!
!	nsleep_max		Maximum number of samples skipped over
!				(Currently fixed at 1. Do not change)
!!       !integer, parameter :: maxseq = nd_max*nsample_max
!-----------------------------------------------------------------------
!

       integer, parameter :: nsample_max = 100
       integer, parameter :: nit_max = 5000
       integer, parameter :: nd_max = 6
       integer, parameter :: ntimesduringMWP = 37
       integer, parameter :: nlat = 512
       integer, parameter :: nlon = 1024
       integer, parameter :: nh_max = 1000
       integer, parameter :: nsleep_max = 1
       integer, parameter :: maxseq = 60
       integer, parameter :: nmod_max = nsample_max*(nit_max+1)
       integer, parameter :: ndata= 6
       integer, parameter :: first_tstep = 72 ! first timestep in times.dat after checkpoint
                                        ! i.e. the timestep following 14.675 ka

