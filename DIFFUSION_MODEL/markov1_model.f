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
     
      integer ii,jj,nbins,npoints
      parameter(ii = 512,jj = 512,nbins = 10,npoints=4000)
      real*8 psi1_av(ii,jj),psi2_av(ii,jj)
      real*8 a1(ii,jj),b1(ii,jj),c1(ii,jj),d1(ii,jj)
      real*8 a2(ii,jj),b2(ii,jj),c2(ii,jj),d2(ii,jj)
      real*8 asigma1(2,2,ii,jj),bsigma1(2,2,ii,jj)
     & ,csigma1(2,2,ii,jj),dsigma1(2,2,ii,jj)
      real*8 asigma2(2,2,ii,jj),bsigma2(2,2,ii,jj)
     & ,csigma2(2,2,ii,jj),dsigma2(2,2,ii,jj)
      real*8 sigma1(2,2,ii,jj),sigma2(2,2,ii,jj),theta1(nbins,2,2)
      real*8 time_av
      real*8 dsigma1(2,2,2,ii,jj),dsigma2(2,2,2,ii,jj)
      
      real*8 x1(npoints,nbins),x2(npoints,nbins),y1(npoints,nbins)
      real*8 u2(npoints,nbins)
      real*8 ufluc_old(2,npoints,nbins),ufluc_new(2,npoints,nbins)
      real*8 drift(2)
      
      integer k,l,m,n,i,j
      
c DETERMINE BIN WIDTHS 

      do b = 1,nbins+1
        bin_corners(b) = (b-1)*(dfloat(jj)/nbins)
      enddo
      
      do b = 1,nbins
        bin_centres(b) = .5*(bin_corners(b)+bin_corners(b+1))
      enddo
      
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
      
      do i = 1,2
      do j = 1,2
      
      call cubic_coeff_x(ii,jj,sigma1(i,j,:,:),asigma1(i,j,:,:)
     & ,bsigma1(i,j,:,:),csigma1(i,j,:,:),dsigma1(i,j,:,:))
      call cubic_coeff_x(ii,jj,sigma2(i,j,:,:),asigma2(i,j,:,:)
     & ,bsigma2(i,j,:,:),csigma2(i,j,:,:),dsigma2(i,j,:,:))
      
      enddo
      enddo
      
c 3. Derivatives of the velocity variance stored at each grid point

      call variance_derivative(sigma1,ii,jj,dsigma1)
      call variance_derivative(sigma2,ii,jj,dsigma2)

c 4. Random forcing ampltiude

      call random_forcing_1(sigma1,theta1,ii,jj,nbins,b1)
      call random_forcing_1(sigma2,theta2,ii,jj,nbins,b2)
      
c b1 and b2 are defined across each zonal bin, so need to interpolate accordingly

c calculate coefficients for 2d-cubic spatial interpolation
c for the velocity variance and derivatives



      
C GENERATE RANDOM LAGRANGIAN PARTICLES

      iseed = 123456789

      do n = 1,npoints
      do b = 1,nbins
      
        x1(n,b) = ran1(iseed)*dfloat(jj)
        y1(n,b) = ran1(iseed)*dfloat(jj)/nbins + bin_corners(b)
        x2(n,b) = x1(n,b)
        y2(n,b) = y1(n,b)
        
        x1_traj(n,b) = x1(n,b)
        x2_traj(n,b) = x2(n,b)
        y1_traj(n,b) = y1(n,b)
        y2_traj(n,b) = y2(n,b)
      
      enddo
      enddo
      
c WRITE INITIAL POSITIONS TO FILE

      call create_diffusion_trajfile(file_name,npoints,nbins)
      
      nrec = 1 
      
      call write_diffusion_trajfile(file_name,npoints,nbins
     & ,x1,y1,x2,y2,dfloat(0),nrec)
     
      nrec = nrec + 1
      
      
c DETERMINE TIME ARRAY

      t_len  = int((max_run*86400.)/dt) ! number of time steps
      allocate(time(t_len))
    
      time(1) = 0.
      k_s = 0
      
C NON DIMENSIONALISE TIME

      scale = basinscale/dfloat(ii)
      uscale = 1
      tscale = scale/uscale

      dt_nondim = dt/tscale
      
c DETERMINE INITIAL CONDITIONS, I.E WHAT IS u'? set to 0 to test

      ufluc = 0.

C MAIN CYCLE

      do t = 2,t_len
      
      do b = 1,nbins
      do p = 1,npoints
c INTERPOLATE PARAMETERS TO FIND VALUES AT PARTICLE LOCATIONS
C DETERMINE DRIFT CORRECTION TERM

        do d = 1,2
            drift(d) = .5*(dsigma(d,d,d) - dsigma(d,d,d)*(u_mean
     & ufluc_old(d,n,b))*dinvsigma(d,d,d)*ufluc_old(d,n,b) - dsigma(d,d,d)            
      
      
      



      
      enddo
      
      end program markov1_model
