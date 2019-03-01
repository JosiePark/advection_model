      program markov1_model
      
!      use stochastic_parameters
      use mod_vel_variance_netcdf
      use mod_theta_netcdf
      use mod_2dcubic
      use mod_qg2_netcdf
      
      implicit none
      
      character*(*),parameter :: home_dir = 
     & '/home/clustor2/ma/j/jp1115/DATA/1/'
     
      character*(*), parameter :: ave_file = 
     &   trim(home_dir) // 'QG/QG_ave.nc'
      character*(*), parameter :: theta_file = 
     &   trim(home_dir) // 'STATS/THETA/theta-1.nc'
      character*(*), parameter :: sigma_file = 
     &   trim(home_dir) // 'STATS/SIGMA/velocity_variance.nc'
     
      integer ii,jj,nbins
      parameter(ii = 512,jj = 512,nbins = 10)
      real*8 psi1_av(ii,jj),psi2_av(ii,jj)
      real*8 a1(ii,jj),b1(ii,jj),c1(ii,jj),d1(ii,jj)
      real*8 a2(ii,jj),b2(ii,jj),c2(ii,jj),d2(ii,jj)
      real*8 sigma1(ii,jj),sigma2(ii,jj),theta(nbins,2,2)
      real*8 time_av
      
C READ TIME-AVERAGED STREAM FUNCTION

      call read_ave_file(ave_file,ii,jj,psi1_av,psi2_av,time_av)
      
C CALCULATE COEFFICIENTS FOR 2D-CUBIC SPATIAL INTERPOLATION

      call cubic_coeff_x(ii,jj,psi1_av,a1,b1,c1,d1)
      call cubic_coeff_x(ii,jj,psi2_av,a2,b2,c2,d2)
      
c CALCULATE COEFFICIENTS FOR MARKOV-1 MODEL

c 1. Markov-1 time tensor (defined for each bin, read from file)
      call read_theta_netcdf(theta_file,nbins,theta)
      
c 2. Velocity variance (read from file)
      call read_vel_variance(sigma_file,sigma1,sigma2,ii,jj)
      
c 3. Derivatives of the velocity variance
c 4. Random forcing ampltiude
c      b1 = sqrt(2*sigma1/theta1)
c      b2 = sqrt(2*sigma2/theta2)

c calculate coefficients for 2d-cubic spatial interpolation
c for the velocity variance and derivatives
      
c DETERMINE BINS

c      do b = 1,nbins+1
      
c      bin_corners(b) = (b-1)*(dfloat(jj)/nbins)
      
c      enddo

      
cc CALCULATE BIN CENTRES

c      do b = 1,nbins
      
c        bin_centres(b) = .5*(bin_corners(b)+bin_corners(b+1))
      
c      enddo
      


cC GENERATE RANDOM LAGRANGIAN PARTICLES

c      n = 0
c      iseed = 123456789

c      do n = 1,npoints
c      do b = 1,nbins
      
c        x1(n,b) = ran1(iseed)*dfloat(jj)
c        y1(n,b) = ran1(iseed)*dfloat(jj)/nbins + bin_corners(b)
c        x2(n,b) = x1(n,b)
c        y2(n,b) = y1(n,b)
        
c        x1_traj(n,b) = x1(n,b)
c        x2_traj(n,b) = x2(n,b)
c        y1_traj(n,b) = y1(n,b)
c        y2_traj(n,b) = y2(n,b)
      
c      enddo
c      enddo
      
cc WRITE INITIAL POSITIONS TO FILE

c      call create_diffusion_trajfile(file_name,npoints,nbins)
      
c      nrec = 1 
      
c      call write_diffusion_trajfile(file_name,npoints,nbins
c     & ,x1,y1,x2,y2,dfloat(0),nrec)
     
c      nrec = nrec + 1
      
      
cc DETERMINE TIME ARRAY

c      t_len  = int((max_run*86400.)/dt) ! number of time steps
c      allocate(time(t_len))
      
c      time(1) = 0.
c      k_s = 0
      
cC NON DIMENSIONALISE TIME

c      scale = basinscale/dfloat(ii)
c      uscale = 1
c      tscale = scale/uscale

c      dt_nondim = dt/tscale

cC MAIN CYCLE

c      do t = 2,t_len
      
      



      
c      enddo
      
      end program markov1_model
