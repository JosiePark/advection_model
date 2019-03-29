      program markov1_model
      
      use mod_stochastic_parameters
      use mod_vel_variance_netcdf
      use mod_theta_netcdf
      use mod_2dcubic
      use mod_qg2_netcdf
      use mod_random
      use mod_diffusion_netcdf
      
      implicit none
      
      character*(*),parameter :: home_dir = 
     & '/home/clustor2/ma/j/jp1115/DATA/1/'
     
      character*(*), parameter :: ave_file = 
     &   trim(home_dir) // 'QG/QG_ave.nc'
      character*(*), parameter :: theta_file = 
     &   trim(home_dir) // 'STATS/THETA/theta-1.nc'
      character*(*), parameter :: sigma_file = 
     &   trim(home_dir) // 'STATS/SIGMA/velocity_variance.nc'
      character*(*), parameter :: file_name =
     &   trim(home_dir) // 'TRAJ/DIFFUSION/markov1.nc' 
     
      integer ii,jj,nbins,npoints
      parameter(ii = 512,jj = 512,nbins = 10,npoints=4000)
      real*8 psi1_av(ii,jj),psi2_av(ii,jj)
      real*8 a1(ii,jj),b1(ii,jj),c1(ii,jj),d1(ii,jj)
      real*8 a2(ii,jj),b2(ii,jj),c2(ii,jj),d2(ii,jj)
      real*8 tmp1(2,2,ii,jj),tmp2(2,2,ii,jj)
      real*8 asigma1(2,ii,jj),bsigma1(2,ii,jj)
     & ,csigma1(2,ii,jj),dsigma1(2,ii,jj)
      real*8 asigma2(2,ii,jj),bsigma2(2,ii,jj)
     & ,csigma2(2,ii,jj),dsigma2(2,ii,jj)
      real*8 adsigma1(2,ii,jj),bdsigma1(2,ii,jj)
     & ,cdsigma1(2,ii,jj),ddsigma1(2,ii,jj)
     & ,adsigma2(2,ii,jj),bdsigma2(2,ii,jj)
     & ,cdsigma2(2,ii,jj),ddsigma2(2,ii,jj)
      real*8 sigma1(2,ii,jj),sigma2(2,ii,jj),theta(nbins,2,2)
      real*8 time_av
      real*8 d_sigma1(2,ii,jj),d_sigma2(2,ii,jj) ! derivatives of sigma
      
      real*8 x1(npoints,nbins),x2(npoints,nbins),y1(npoints,nbins)
     & ,y2(npoints,nbins)
     
      real*8 x1_traj(npoints,nbins),y1_traj(npoints,nbins)
     & ,x2_traj(npoints,nbins), y2_traj(npoints,nbins) 
      real*8 u2(npoints,nbins)
      real*8 ufluc_old(2,npoints,nbins),ufluc_new(2,npoints,nbins)
      real*8 drift(2)
      
      real*8 bin_corners(nbins+1),bin_centres(nbins)
      real*8 tmp(2,2,ii,jj)
      
      real*8 forcing_b1(2,ii,jj),forcing_b2(2,ii,jj)
      real*8 ab1(2,ii,jj),bb1(2,ii,jj),cb1(2,ii,jj),db1(2,ii,jj)
      real*8 ab2(2,ii,jj),bb2(2,ii,jj),cb2(2,ii,jj),db2(2,ii,jj)
      
      integer k,l,m,n,i,j,b,iseed,t,p,d
      
      real*8 u_meanufluc_old,dinvsigma(2,2,2)
      
      integer nrec,t_len,k_s
      real*8 dt
      real*8 max_run
      parameter(max_run = 1000)
      
      real*8 scale,basinscale,uscale,tscale,dt_nondim,ufluc
      parameter(basinscale = 520.d5)
      
      real*8, allocatable, dimension(:) :: time
    
      
      
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
      call read_vel_variance(sigma_file,tmp1,tmp2,ii,jj)
      sigma1(1,:,:) = tmp(1,1,:,:)
      sigma2(2,:,:) = tmp(2,2,:,:)
      
      print*,'vel variance read'
      
      do i = 1,2
      
      call cubic_coeff_x(ii,jj,sigma1(i,:,:),asigma1(i,:,:)
     & ,bsigma1(i,:,:),csigma1(i,:,:),dsigma1(i,:,:))
      call cubic_coeff_x(ii,jj,sigma2(i,:,:),asigma2(i,:,:)
     & ,bsigma2(i,:,:),csigma2(i,:,:),dsigma2(i,:,:))
      
      enddo
      
c 3. Derivatives of the velocity variance stored at each grid point

      call variance_derivative(sigma1,ii,jj,d_sigma1)
      call variance_derivative(sigma2,ii,jj,d_sigma2)
      
      print*,'variance derivative calculated'
      
      do i = 1,2
      
      call cubic_coeff_x(ii,jj,d_sigma1(i,:,:),adsigma1(i,:,:)
     & ,bdsigma1(i,:,:),cdsigma1(i,:,:),ddsigma1(i,:,:))
      call cubic_coeff_x(ii,jj,d_sigma2(i,:,:),adsigma2(i,:,:)
     & ,bdsigma2(i,:,:),cdsigma2(i,:,:),ddsigma2(i,:,:))
      
      enddo

c 4. Random forcing ampltiude

      call random_forcing_1(sigma1,theta(:,:,1),ii,jj,nbins,bin_centres,
     & forcing_b1) ! layer 1
      call random_forcing_1(sigma2,theta(:,:,2),ii,jj,nbins,bin_centres,
     & forcing_b2) ! layer 2
     
      print*,'random forcing calculated'
      
      do i = 1,2
      
      call cubic_coeff_x(ii,jj,forcing_b1(i,:,:),ab1(i,:,:)
     & ,bb1(i,:,:),cb1(i,:,:),db1(i,:,:))
      call cubic_coeff_x(ii,jj,forcing_b2(i,:,:),ab2(i,:,:)
     & ,bb2(i,:,:),cb2(i,:,:),db2(i,:,:))
      
      enddo 
      
      stop


      
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
c INTERPOLATE PARAMETERS TO FIND VALUES AT PARTICLE LOCATIONS
C DETERMINE DRIFT CORRECTION TERM

      ! sigma
      ! dsigma
      ! theta
      ! b
      ! mean velocity          
      
      
      
      
      end program markov1_model
