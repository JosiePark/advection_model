C CODE THAT CALCULATE THE LAGRANGIAN AUTOCORRELATION FUNCTION

      program autocorrelation
      
      use mod_vel_variance_netcdf
      use mod_qg2_netcdf
      use mod_traj_netcdf
      use mod_trajbin_netcdf
      
      implicit none
      
      integer nbins,nrel,npoints,t_len
      
      real*8, dimension(:,:,:,:), allocatable :: traj
      
      character*(*),parameter :: home_dir = 
     & '/home/clustor2/ma/j/jp1115/DATA/1/'
      
      character*(*), parameter :: traj_file = 
     & trim(home_dir) //'TRAJ/UNIFORM_BINS/
     &full_uniform_bins_trajectories.nc' 
     
      character*(*), parameter :: var_file = 
     & trim(home_dir) //'STATS/SIGMA/
     &velocity_variance.nc' 
     
      integer,parameter :: ii = 512,jj = 512
      integer,parameter :: ntau = 200
      
      real*8 sigma1(2,2,ii,jj),sigma2(2,2,ii,jj)
     
      integer b,k,t
      
      ! READ DIMENSIONS IN TRAJECTORY FILE
      
      print*,'traj_file = ',traj_file
      call read_binned_file_dimensions(traj_file
     & ,nbins,nrel,npoints,t_len)
      print*,'nbins,nrel,npoints,t_len=',nbins,nrel,npoints,t_len
      allocate(traj(2,npoints,2,t_len))
      print*,'trajectory allocated'
      
      ! READ VELOCITY VARIANCE FILE
      
      call read_vel_variance(var_file,sigma1,sigma2,ii,jj)
      
      
      do b = 1,nbins
        do k = 1,nrel
        ! READ TRAJECTORIES ONE BIN AND RELEASE AT A TIME
        call read_trajbin_file(traj_file,b,k,npoints,t_len,ii,traj)
        print*,'Read trajectory for bin bm release k = ', b,k
        
        ! CALCULATE LAGRANGIAN VELOCITY
        
        
        enddo
      enddo
      
      
      
      
      ! CALCULATE AUTOCORRELATION
      
      ! AUTOCORRELATION IS AVERAGED OVER BINS
      
      end program autocorrelation
