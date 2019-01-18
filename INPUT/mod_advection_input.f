C MODULE THAT CONTAINS INPUT FOR TRANSPORT MODEL

      module mod_advection_input
      
      implicit none
      
      integer npoints_sqrt,npoints,isolve,i_pseudo,i_full,i_eddy,regime
     & ,sim,ii,jj,nbins,i_pvbin,jj1,i_start,coord_range
      real*8 H1,H2,Rd,basinscale,visc,U_0,beta
      real*8 release_interval,release_length
      integer release_no
      integer k_save
      real*8 dt
      integer eof_option
      
      
      
c DYNAMICAL MODEL PARAMETERS
      
      parameter(ii=512,jj=512,jj1 = 3*jj)
      parameter(basinscale=520.D5
     & ,H1=1.D5,H2=3.D5
     & ,beta = 2.D-13
     & ,U_0=6.D0
     & ,Rd= 25.D5
     & ,coord_range =3)
      
C TRANSPORT MODEL PARAMETERS 

      parameter(npoints=5000)
      parameter(isolve=1 ! if isolve = 0 :bicubic, if isolve = 1 :2Dcubic, is isolve = 2: old non-divergent preserving cubic method
     & ,i_eddy=0
     & ,i_pseudo=0
     & ,i_full=1
     & ,release_interval = 200.
     & ,release_length = 1000.
     & ,release_no = 9
     & ,k_save = 1 ! save every k_save days
     & ,dt = 1800
     & ,nbins = 10
     & ,i_pvbin = 0 ! set = 1 if you wish to bin domain according to PV, otherwise the domain is binned uniformly 
     & ,i_start = 0) ! set = 1 if you wish to restart
     
C SET EOF_OPTION TO 1 IF YOU WISH TO ADD EOF MODES TOGETHER
c OTHERWISE SET TO 0

      parameter(eof_option = 0)
     
      character*(*),parameter :: home_dir = 
     & '/home/clustor2/ma/j/jp1115/DATA/1/'

     
       character*(*), parameter :: file_name = 
     &   trim(home_dir) // 'QG/QG.nc'
       character*(*), parameter :: ave_name =
     &   trim(home_dir) // 'QG/QG_ave.nc'
       character*(*), parameter :: full_traj_file =
     &   trim(home_dir) // 'TRAJ/UNIFORM_BINS/
     &EOF_trajectories_minus_1-10_full_new.nc'  
       character*(*), parameter :: traj_file =
     &  trim(home_dir) // 'TRAJ/UNIFORM_BINS/
     &EOF_trajectories_minus_1-10_pseudo_new.nc'
        character*(*), parameter :: eddy_traj_file =
     & trim(home_dir) //'TRAJ/UNIFORM_BINS/
     &EOF_trajectories_minus_1-10_eddy_new.nc'
       character*(*), parameter :: eof_file =
     &  trim(home_dir) // 'QG/psi_minus_1-10.nc'
       character*(*), parameter :: pseudo_name =
     &   trim(home_dir) // 'TRAJ/UNIFORM_BINS/NEW
     &pseudo_uniform_bins_trajectories.nc' 
       character*(*), parameter :: eddy_name =
     &   trim(home_dir) // 'TRAJ/PV_BINS/
     &eddy_PV_bins_trajectories.nc'
       character*(*), parameter :: full_name =
     &   trim(home_dir) // 'TRAJ/UNIFORM_BINS/NEW
     &full_uniform_bins_trajectories.nc'
     
      end module mod_advection_input
