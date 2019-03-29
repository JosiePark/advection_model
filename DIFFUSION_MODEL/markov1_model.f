      program markov1_model
      
      use mod_stochastic_parameters
      use mod_vel_variance_netcdf
      use mod_theta_netcdf
      use mod_2dcubic
      use mod_qg2_netcdf
      use mod_random
      use mod_diffusion_netcdf
      
      implicit none
      
      character*(*), parameter :: home_dir = 
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
      
      real*8 a_dinvsigma1(2,ii,jj),b_dinvsigma1(2,ii,jj)
     & ,c_dinvsigma1(2,ii,jj),d_dinvsigma1(2,ii,jj)
     & ,a_dinvsigma2(2,ii,jj),b_dinvsigma2(2,ii,jj)
     & ,c_dinvsigma2(2,ii,jj),d_dinvsigma2(2,ii,jj)
      real*8 d_invsigma1(2,ii,jj),d_invsigma2(2,ii,jj)
      
      real*8 x1(npoints,nbins),x2(npoints,nbins),y1(npoints,nbins)
     & ,y2(npoints,nbins)
     
      real*8 x1_traj(npoints,nbins),y1_traj(npoints,nbins)
     & ,x2_traj(npoints,nbins), y2_traj(npoints,nbins) 
      real*8 drift(2)
      
      real*8 bin_corners(nbins+1),bin_centres(nbins)
      real*8 tmp(2,2,ii,jj)
      
      real*8 forcing_b1(2,ii,jj),forcing_b2(2,ii,jj)
      real*8 ab1(2,ii,jj),bb1(2,ii,jj),cb1(2,ii,jj),db1(2,ii,jj)
      real*8 ab2(2,ii,jj),bb2(2,ii,jj),cb2(2,ii,jj),db2(2,ii,jj)
      
      real*8 sigma1_tmp(2),sigma2_tmp(2),d_sigma1_tmp(2),d_sigma2_tmp(2)
      real*8 theta1_tmp(2),theta2_tmp(2)
      real*8 d_invsigma1_tmp(2),d_invsigma2_tmp(2)
      real*8 b1_tmp(2),b2_tmp(2)
      real*8 u1,u2,v1,v2
      real*8 psi1_x(4),psi2_x(4)
      
      real*8 drift1(2),drift2(2)
      real*8 u1_fluc_old(2,npoints,nbins)
     & ,u2_fluc_old(2,npoints,nbins)
     & ,u1_fluc_new(2,npoints,nbins)
     & ,u2_fluc_new(2,npoints,nbins)
      real*8 u_top(2), u_bottom(2)
      
      integer k,l,m,n,i,j,b,iseed,t,p,d
      
      real*8 k_save
      
      parameter(k_save = 1)
         
      integer nrec,t_len,k_s
      real*8 dt
      real*8 max_run
      parameter(max_run = 1000)
      parameter(dt = 2160.)
      
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
      
c 3.a. Derivatives of the velocity variance stored at each grid point

      call variance_derivative(sigma1,ii,jj,d_sigma1)
      call variance_derivative(sigma2,ii,jj,d_sigma2)
      
      print*,'variance derivative calculated'
      
      do i = 1,2
      
      call cubic_coeff_x(ii,jj,d_sigma1(i,:,:),adsigma1(i,:,:)
     & ,bdsigma1(i,:,:),cdsigma1(i,:,:),ddsigma1(i,:,:))
      call cubic_coeff_x(ii,jj,d_sigma2(i,:,:),adsigma2(i,:,:)
     & ,bdsigma2(i,:,:),cdsigma2(i,:,:),ddsigma2(i,:,:))
      
      enddo
      
c 3.b Derivatives of the inverse of the velocity variance at each grid point

      call variance_derivative(1./sigma1,ii,jj,d_invsigma1)
      call variance_derivative(1./sigma2,ii,jj,d_invsigma2)

      do i = 1,2
      
      call cubic_coeff_x(ii,jj,d_invsigma1(i,:,:),a_dinvsigma1(i,:,:)
     & ,b_dinvsigma1(i,:,:),c_dinvsigma1(i,:,:),d_dinvsigma1(i,:,:))
      call cubic_coeff_x(ii,jj,d_invsigma2(i,:,:),adsigma2(i,:,:)
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

      u1_fluc_old = 0.
      u2_fluc_old = 0.

C MAIN CYCLE

      do t = 2,t_len

            time(t) = time(t-1) + dt/86400. ! time in days
            k_s = k_s + dt/86400.
            
      print*,'time =',time(t)
      do b = 1,nbins
      do n = 1,npoints
c INTERPOLATE PARAMETERS TO FIND VALUES AT PARTICLE LOCATIONS
C DETERMINE DRIFT CORRECTION TERM

      do d = 1,2
      ! sigma
      
      call cubic_interp_full(ii,jj,asigma1(d,:,:),bsigma1(d,:,:)
     & ,csigma1(d,:,:),dsigma1(d,:,:)
     & ,x1(n,b),y1(n,b),sigma1_tmp(d))
      call cubic_interp_full(ii,jj,asigma2(d,:,:),bsigma2(d,:,:)
     & ,csigma2(d,:,:),dsigma2(d,:,:)
     & ,x2(n,b),y2(n,b),sigma2_tmp(d))
     
     
      ! dsigma
      
      call cubic_interp_full(ii,jj,adsigma1(d,:,:),bdsigma1(d,:,:)
     & ,cdsigma1(d,:,:),ddsigma1(d,:,:)
     & ,x1(n,b),y1(n,b),d_sigma1_tmp(d))
      call cubic_interp_full(ii,jj,adsigma2(d,:,:),bdsigma2(d,:,:)
     & ,cdsigma2(d,:,:),ddsigma2(d,:,:)
     & ,x2(n,b),y2(n,b),d_sigma2_tmp(d))
     
      ! theta
      
      
      call diffusivity_interp_1d(bin_centres,nbins,theta(:,d,1)
     &   ,ii,y1(n,b),theta1_tmp(d))
      call diffusivity_interp_1d(bin_centres,nbins,theta(:,d,2)
     &   ,ii,y2(n,b),theta2_tmp(d))
      
      ! b
      
      call cubic_interp_full(ii,jj,ab1(d,:,:),bb1(d,:,:)
     & ,cb1(d,:,:),db1(d,:,:)
     & ,x1(n,b),y1(n,b),b1_tmp(d))
      call cubic_interp_full(ii,jj,ab2(d,:,:),bb2(d,:,:)
     & ,cb2(d,:,:),db2(d,:,:)
     & ,x2(n,b),y2(n,b),b2_tmp(d))
      
      ! d_invsigma
     
      call cubic_interp_full(ii,jj,a_dinvsigma1(d,:,:)
     & ,b_dinvsigma1(d,:,:)
     & ,c_dinvsigma1(d,:,:),d_dinvsigma1(d,:,:)
     & ,x1(n,b),y1(n,b),d_invsigma1_tmp(d))
      call cubic_interp_full(ii,jj,a_dinvsigma2(d,:,:)
     & ,b_dinvsigma2(d,:,:)
     & ,c_dinvsigma2(d,:,:),d_dinvsigma2(d,:,:)
     & ,x2(n,b),y2(n,b),d_invsigma2_tmp(d))
     
      
      enddo    
      
      ! mean velocity 
      
      call cubic_poly_x(ii,jj,x1(n,b),y1(n,b)
     & ,a1,b1,c1,d1,psi1_x)
      call cubic_poly_x(ii,jj,x2(n,b),y2(n,b)
     & ,a2,b2,c2,d2,psi2_x)
      
      call vel(ii,jj,psi1_x,a1,b1,c1,d1,x1(n,b),y1(n,b),u1,v1)    
      call vel(ii,jj,psi2_x,a2,b2,c2,d2,x2(n,b),y2(n,b),u2,v2) 
      
      u_top(1) = u1
      u_top(2) = v1
      u_bottom(1) = u2
      u_bottom(2) = v2
      
      ! calculate drift correction term
      
      do d = 1,2
        drift1(d) = .5*d_sigma1_tmp(d) 
     &   - .5*sigma1_tmp(d)*(u1_fluc_old(d,n,b) + u_top(d))
     & *d_invsigma1_tmp(d)*u1_fluc_old(d,n,b)
     & -.5*sigma1_tmp(d)*d_invsigma1_tmp(d)*u1_fluc_old(d,n,b)**2   
     
       drift2(d) = .5*d_sigma2_tmp(d) 
     &   - .5*sigma2_tmp(d)*(u2_fluc_old(d,n,b) + u_bottom(d))
     & *d_invsigma2_tmp(d)*u2_fluc_old(d,n,b)
     & -.5*sigma2_tmp(d)*d_invsigma2_tmp(d)*u2_fluc_old(d,n,b)**2   
      
      enddo
      
      ! update velocity fluctuation
      
      do d = 1,2
      
       u1_fluc_new(d,n,b) = u1_fluc_old(d,n,b) 
     &  - (u1_fluc_old(d,n,b)/theta1_tmp(d)
     & + drift1(d))*dt_nondim 
     & + b1_tmp(d)*random_normal()*dsqrt(dt_nondim)  
       
       u2_fluc_new(d,n,b) = u2_fluc_old(d,n,b)
     &   - (u2_fluc_old(d,n,b)/theta2_tmp(d)
     & + drift2(d))*dt_nondim 
     & + b2_tmp(d)*random_normal()*dsqrt(dt_nondim) 
     
      enddo
      
      ! update particle location
      
      x1(n,b) = x1(n,b) + dt_nondim*(u1_fluc_new(1,n,b) + u1)
      x2(n,b) = x2(n,b) + dt_nondim*(u2_fluc_new(1,n,b) + u2)
      y1(n,b) = y1(n,b) + dt_nondim*(u1_fluc_new(2,n,b) + v1)
      y2(n,b) = y2(n,b) + dt_nondim*(u2_fluc_new(2,n,b) + v2)
      
      x1_traj(n,b) = x1_traj(n,b) + dt_nondim*(u1_fluc_new(1,n,b) + u1)
      x2_traj(n,b) = x2_traj(n,b) + dt_nondim*(u2_fluc_new(1,n,b) + u2)
      y1_traj(n,b) = y1_traj(n,b) + dt_nondim*(u1_fluc_new(2,n,b) + v1)
      y2_traj(n,b) = y2_traj(n,b) + dt_nondim*(u2_fluc_new(2,n,b) + v2)
      
      ! treat boundaries
      
            if(x1(n,b).le.0) then
                x1(n,b) = x1(n,b) + dfloat(ii)
            endif
            if(x1(n,b).ge.ii) then
                x1(n,b) = x1(n,b) - dfloat(ii)
            endif
            if(x2(n,b).le.0) then
                x2(n,b) = x2(n,b) + dfloat(ii)
            endif
            if(x2(n,b).ge.ii) then
                x2(n,b) = x2(n,b) - dfloat(ii)
            endif
            
            if(y1(n,b).le.0) then
                y1(n,b) = y1(n,b) + dfloat(jj)
            endif
            if(y1(n,b).ge.jj) then
                y1(n,b) = y1(n,b) - dfloat(jj)
            endif
            if(y2(n,b).le.0) then
                y2(n,b) = y2(n,b) + dfloat(jj)
            endif
            if(y2(n,b).ge.jj) then
                y2(n,b) = y2(n,b) - dfloat(jj)
            endif
            
            ! update velocoty fluctuations
            
      do d = 1,2
        u1_fluc_old(d,n,b) = u1_fluc_new(d,n,b)
        u2_fluc_old(d,n,b) = u2_fluc_new(d,n,b)
      enddo
      
      enddo
      enddo
      
      ! save
      
      if(k_s + 0.001 .ge. k_save) then
      
        print*,'Saving at time',time(t)
        
        call write_diffusion_trajfile(file_name,npoints,nbins
     & ,x1_traj,y1_traj,x2_traj,y2_traj,time(t),nrec)
     
      nrec = nrec + 1
      
      k_s = 0.
      
      endif
            
      enddo
      
      
      end program markov1_model
