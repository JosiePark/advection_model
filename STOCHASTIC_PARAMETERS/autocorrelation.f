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
      real*8,parameter :: dt = 1.
      
      real*8 sigma1(2,2,ii,jj),sigma2(2,2,ii,jj)
      real*8 a1(ii,jj),a2(ii,jj),b1(ii,jj),b2(ii,jj),c1(ii,jj),c2(ii,jj)
     & ,d1(ii,jj),d2(ii,jj) 
      real*8, dimension(:,:,:,:), allocatable :: u
      real*8, dimension(:,:,:,:), allocatable :: R
     
      integer b,k,t,tau,d,m,l,d1,d2
      
      ! READ DIMENSIONS IN TRAJECTORY FILE
      
      print*,'traj_file = ',traj_file
      call read_binned_file_dimensions(traj_file
     & ,nbins,nrel,npoints,t_len)
      print*,'nbins,nrel,npoints,t_len=',nbins,nrel,npoints,t_len
      allocate(traj(2,npoints,2,t_len))
      allocate(u(2,npoints,2,t_len-1))
      allocate(R(2,2,npoints,2,ntau))
      print*,'trajectory allocated'
      
      ! READ VELOCITY VARIANCE FILE
      
      call read_vel_variance(var_file,sigma1,sigma2,ii,jj)
      
      ! FIND INTERPOLATION COEFFICIENTS OF SIGMA1 AND SIGMA2    
      
      call cubic_coeff_x(ii,jj,sigma1,a1,b1,c1,d1)
      call cubic_coeff_x(ii,jj,sigma2,a2,b2,c2,d2)
      
      
      do b = 1,nbins
        do k = 1,nrel
        ! READ TRAJECTORIES ONE BIN AND RELEASE AT A TIME
        call read_trajbin_file(traj_file,b,k,npoints,t_len,ii,traj)
        print*,'Read trajectory for bin bm release k = ', b,k
        
        ! CALCULATE LAGRANGIAN VELOCITY
        
        do d = 1,2
        do l = 1,2
            u(d,1:npoints,l,1:t_len-1) = (traj(d,1:npoints,l,2:t_len)
     &        - traj(d,1:npoints,l,1:t_len-1))/dt
        enddo
        enddo
        
        print*,'lagrangian velocity calculated'
        
        
        ! CALCULATE AUTOCORRELATION
        
        do d1 = 1,2
        do d2 = 1,2
        do l = 1,2
        do tau = 0,ntau-1
            do t = 1,(nt-tau-1)
            
                R(d1,d2,1:npoints,l,tau+1) = R(d1,d2,1:npoints,l,tau+1) 
     &                + u(d1,1:npoints,l,t)*u(d2,1:npoints,l,t+tau)
            
            enddo
        enddo
        enddo
        enddo
        
        ! INTERPOLATE SIGMA1 AND SIGMA2 TO FIND VALUES AT PARTICLE LOCATIONS
        
        do l = 1,2
        do n = 1,npoints
        
            call cubic_poly_x(ii,jj,traj(1,n,
        enddo
        enddo
        
        
        
        enddo
      enddo
      
      
      
      
      ! CALCULATE AUTOCORRELATION
      
      ! AUTOCORRELATION IS AVERAGED OVER BINS
      
      end program autocorrelation
