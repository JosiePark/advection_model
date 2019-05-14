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
     &   trim(home_dir) // 'STATS/THETA/theta.nc'
      character*(*), parameter :: sigma_file = 
     &   trim(home_dir) // 'STATS/SIGMA/velocity_variance.nc'
      character*(*), parameter :: file_name =
     &   trim(home_dir) // 'TRAJ/DIFFUSION/markov1_new.nc' 
      character*(*), parameter :: d_sigma_file =
     &   trim(home_dir) // 'TRAJ/DIFFUSION/d_sigma.nc' 
      character*(*), parameter :: d_invsigma_file =
     &   trim(home_dir) // 'TRAJ/DIFFUSION/d_invsigma.nc'
      character*(*), parameter :: b_file =
     &   trim(home_dir) // 'TRAJ/DIFFUSION/random_forcing.nc' 
     
      integer ii,jj,nbins,npoints
      parameter(ii = 512,jj = 512,nbins = 10,npoints=1000)
      real*8 psi1_av(ii,jj),psi2_av(ii,jj)
      real*8 a1(ii,jj),b1(ii,jj),c1(ii,jj),d1(ii,jj)
      real*8 a2(ii,jj),b2(ii,jj),c2(ii,jj),d2(ii,jj)
      real*8 tmp1(2,2,ii,jj),tmp2(2,2,ii,jj)
      real*8 asigma11(ii,jj),bsigma11(ii,jj)
     & ,csigma11(ii,jj),dsigma11(ii,jj) ! zonal, top layer
      real*8 asigma12(ii,jj),bsigma12(ii,jj)
     & ,csigma12(ii,jj),dsigma12(ii,jj) ! meridional, top layer
      real*8 asigma21(ii,jj),bsigma21(ii,jj)
     & ,csigma21(ii,jj),dsigma21(ii,jj) ! zonal, bottom layer
      real*8 asigma22(ii,jj),bsigma22(ii,jj)
     & ,csigma22(ii,jj),dsigma22(ii,jj) ! meridional, bottom layer
      real*8 adsigma11(ii,jj),bdsigma11(ii,jj)
     & ,cdsigma11(ii,jj),ddsigma11(ii,jj)
      real*8 adsigma12(ii,jj),bdsigma12(ii,jj)
     & ,cdsigma12(ii,jj),ddsigma12(ii,jj)
      real*8 adsigma21(ii,jj),bdsigma21(ii,jj)
     & ,cdsigma21(ii,jj),ddsigma21(ii,jj)
      real*8 adsigma22(ii,jj),bdsigma22(ii,jj)
     & ,cdsigma22(ii,jj),ddsigma22(ii,jj)
      real*8 sigma11(ii,jj),sigma12(ii,jj),sigma21(ii,jj),sigma22(ii,jj)
      real*8 theta(2,2,nbins)
      real*8 theta11(nbins),theta12(nbins),theta21(nbins),theta22(nbins)
      real*8 time_av
      real*8 d_sigma12(ii,jj),d_sigma11(ii,jj)
     & ,d_sigma22(ii,jj),d_sigma21(ii,jj) ! derivatives of sigma
      
      real*8 x1(npoints,nbins),x2(npoints,nbins),y1(npoints,nbins)
     & ,y2(npoints,nbins)
     
      real*8 x1_traj(npoints,nbins),y1_traj(npoints,nbins)
     & ,x2_traj(npoints,nbins), y2_traj(npoints,nbins) 
      real*8 drift(2)
      
      real*8 bin_corners(nbins+1),bin_centres(nbins)
      real*8 tmp(2,2,ii,jj)
      
      real*8 forcing_b1(2,ii,jj),forcing_b2(2,ii,jj)
      
      real*8 sigma11_tmp,sigma12_tmp,sigma21_tmp,sigma22_tmp
      real*8 d_sigma11_tmp,d_sigma12_tmp,d_sigma21_tmp,d_sigma22_tmp
      real*8 theta11_tmp,theta12_tmp,theta21_tmp,theta22_tmp
      real*8 b11_tmp,b12_tmp,b21_tmp,b22_tmp
      real*8 u1,u2,v1,v2
      real*8 psi1_x(4),psi2_x(4)
      
      real*8 drift11,drift21,drift12,drift22
      real*8 u1_fluc_old(2,npoints,nbins)
     & ,u2_fluc_old(2,npoints,nbins)
     & ,u1_fluc_new(2,npoints,nbins)
     & ,u2_fluc_new(2,npoints,nbins)
      real*8 u_top(2), u_bottom(2)
      
      integer k,l,m,n,i,j,b,iseed,t,p,d
      
      real*8 k_save
      real*8 sigma11_mean,sigma12_mean,sigma21_mean,sigma22_mean
      real*8 sum_tmp
      
      parameter(k_save = 1)
         
      integer nrec,t_len
      real*8 k_s
      real*8 dt
      real*8 max_run
      parameter(max_run = 1000.)
      !parameter(max_run = 1.)
      parameter(dt = 2160.)
      
      real*8 scale,basinscale,uscale,tscale,dt_nondim,ufluc,U0
      parameter(basinscale = 520.d5)
      parameter(U0 = 6.d0)
      
      real*8, allocatable, dimension(:) :: time
      
      real*8 start_time,stop_time
    
      
      
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

      call cpu_time(start_time)

      call cubic_coeff_x(ii,jj,psi1_av,a1,b1,c1,d1)
      call cubic_coeff_x(ii,jj,psi2_av,a2,b2,c2,d2)
      
      call cpu_time(stop_time)
      
      print*,'calculating coefficients =',stop_time - start_time
      
      
c CALCULATE COEFFICIENTS FOR MARKOV-1 MODEL

c 1. Markov-1 time tensor (defined for each bin, read from file)
      call read_theta_netcdf(theta_file,nbins,theta)
      
      print*,'theta=',theta
      

      
c 2. Velocity variance (read from file) 
      call read_vel_variance(sigma_file,tmp1,tmp2,ii,jj)
      sigma11 = tmp1(1,1,:,:) ! top layer - zonal
      sigma12 = tmp1(2,2,:,:) ! top layer - meridional
      sigma21 = tmp2(1,1,:,:)
      sigma22 = tmp2(2,2,:,:)
      
      print*,'vel variance read'
      
      call cpu_time(start_time)
      
      call cubic_coeff_x(ii,jj,sigma11,asigma11
     & ,bsigma11,csigma11,dsigma11)
      call cubic_coeff_x(ii,jj,sigma12,asigma12
     & ,bsigma12,csigma12,dsigma12)
      call cubic_coeff_x(ii,jj,sigma21,asigma21
     & ,bsigma21,csigma21,dsigma21)
      call cubic_coeff_x(ii,jj,sigma22,asigma22
     & ,bsigma22,csigma22,dsigma22)
      
      
      call cpu_time(stop_time)
      
      print*,'calculating coefficients =',stop_time - start_time
      
c 3.a. Derivatives of the velocity variance stored at each grid point

c      call variance_derivative(sigma1,ii,jj,d_sigma1)
c      call variance_derivative(sigma2,ii,jj,d_sigma2)

c        ! save the derivative to file
        
c      call create_variable_netcdf(d_sigma_file,ii,jj)
c      call write_variable_netcdf(d_sigma_file,d_sigma1,d_sigma2,ii,jj)
        
c      print*,'variance derivative calculated'
      
c      do i = 1,2
      
c      call cubic_coeff_x(ii,jj,d_sigma1(i,:,:),adsigma1(i,:,:)
c     & ,bdsigma1(i,:,:),cdsigma1(i,:,:),ddsigma1(i,:,:))
c      call cubic_coeff_x(ii,jj,d_sigma2(i,:,:),adsigma2(i,:,:)
c     & ,bdsigma2(i,:,:),cdsigma2(i,:,:),ddsigma2(i,:,:))
      
c      enddo
      
c 3.b Derivatives of the inverse of the velocity variance at each grid point

c      call variance_derivative(1./sigma1,ii,jj,d_invsigma1)
c      call variance_derivative(1./sigma2,ii,jj,d_invsigma2)
      
c      call create_variable_netcdf(d_invsigma_file,ii,jj)
c      call write_variable_netcdf(d_invsigma_file
c     & ,d_invsigma1,d_invsigma2,ii,jj)

c      do i = 1,2
      
c      call cubic_coeff_x(ii,jj,d_invsigma1(i,:,:),a_dinvsigma1(i,:,:)
c     & ,b_dinvsigma1(i,:,:),c_dinvsigma1(i,:,:),d_dinvsigma1(i,:,:))
c      call cubic_coeff_x(ii,jj,d_invsigma2(i,:,:),adsigma2(i,:,:)
c     & ,bdsigma2(i,:,:),cdsigma2(i,:,:),ddsigma2(i,:,:))
      
c      enddo

c 4. Random forcing ampltiude

c      call random_forcing_1(sigma1,theta(:,:,1),ii,jj,nbins,bin_centres,
c     & forcing_b1) ! layer 1
c      call random_forcing_1(sigma2,theta(:,:,2),ii,jj,nbins,bin_centres,
c     & forcing_b2) ! layer 2
     
c      call create_variable_netcdf(b_file,ii,jj)
c      call write_variable_netcdf(b_file,forcing_b1,forcing_b2,ii,jj)
     
c      print*,'random forcing calculated'
      
c      do i = 1,2
      
c      call cubic_coeff_x(ii,jj,forcing_b1(i,:,:),ab1(i,:,:)
c     & ,bb1(i,:,:),cb1(i,:,:),db1(i,:,:))
c      call cubic_coeff_x(ii,jj,forcing_b2(i,:,:),ab2(i,:,:)
c     & ,bb2(i,:,:),cb2(i,:,:),db2(i,:,:))
      
c      enddo 
      


      
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
      k_s = 0.
      
C NON DIMENSIONALISE TIME

      scale = basinscale/dfloat(ii)
      uscale = 1
      tscale = scale/uscale

      dt_nondim = dt/tscale
      
      theta11 = theta(1,1,:)*86400./tscale
      theta12 = theta(1,2,:)*86400./tscale
      theta21 = theta(2,1,:)*86400./tscale
      theta22 = theta(2,2,:)*86400./tscale
      
c DETERMINE INITIAL CONDITIONS, 
C u' IS INITIALISED AS A GAUSSIAN RANDOM VARIABLE WITH ZERO MEAN AND VARIANCE SIGMA MEAN

      sum_tmp = 0.
      do i = 1,ii
      do j = 1,jj
        sum_tmp = sum_tmp + sigma11(i,j)
      enddo
      enddo
      sigma11_mean = sum_tmp/(ii*jj)
      
      sum_tmp = 0.
      do i = 1,ii
      do j = 1,jj
        sum_tmp = sum_tmp + sigma12(i,j)
      enddo
      enddo
      sigma12_mean = sum_tmp/(ii*jj)
      
      sum_tmp = 0.
      do i = 1,ii
      do j = 1,jj
        sum_tmp = sum_tmp + sigma21(i,j)
      enddo
      enddo
      sigma21_mean = sum_tmp/(ii*jj)
      
      
      sum_tmp = 0.
      do i = 1,ii
      do j = 1,jj
        sum_tmp = sum_tmp + sigma22(i,j)
      enddo
      enddo
      sigma22_mean = sum_tmp/(ii*jj)

      do n = 1,npoints
      do b = 1,nbins
        u1_fluc_old(1,n,b) = random_normal()*dsqrt(sigma11_mean)
        u1_fluc_old(2,n,b) = random_normal()*dsqrt(sigma12_mean)
        u2_fluc_old(1,n,b) = random_normal()*dsqrt(sigma21_mean)
        u2_fluc_old(2,n,b) = random_normal()*dsqrt(sigma22_mean)        
      enddo
      enddo
      

C MAIN CYCLE

      do t = 2,t_len

            time(t) = time(t-1) + dt/86400. ! time in days
            k_s = k_s + dt/86400.
            !print*,'k_s =',k_s
            
      !print*,'time =',time(t)
      
      do b = 1,nbins
      do n = 1,npoints
      !print*,'b=',b
      !print*,'x1,y1,x2,y2=',x1(n,b),y1(n,b),x2(n,b),y2(n,b)
c INTERPOLATE PARAMETERS TO FIND VALUES AT PARTICLE LOCATIONS
C DETERMINE DRIFT CORRECTION TERM

      
      ! sigma and dsigma
      
      call sigma_interpolation(ii,jj,asigma11,bsigma11
     & ,csigma11,dsigma11,1
     & ,x1(n,b),y1(n,b),sigma11_tmp,d_sigma11_tmp)
      call sigma_interpolation(ii,jj,asigma12,bsigma12
     & ,csigma12,dsigma12,2
     & ,x1(n,b),y1(n,b),sigma12_tmp,d_sigma12_tmp)
      call sigma_interpolation(ii,jj,asigma21,bsigma21
     & ,csigma21,dsigma21,1
     & ,x2(n,b),y2(n,b),sigma21_tmp,d_sigma21_tmp)
      call sigma_interpolation(ii,jj,asigma22,bsigma22
     & ,csigma22,dsigma22,2
     & ,x2(n,b),y2(n,b),sigma22_tmp,d_sigma22_tmp)
     
      !print*,'sigma,dsigma=',sigma11_tmp,d_sigma11_tmp
      !print*,'sigma,dsigma=',sigma12_tmp,d_sigma12_tmp
      !print*,'sigma,dsigma=',sigma21_tmp,d_sigma21_tmp
      !print*,'sigma,dsigma=',sigma22_tmp,d_sigma22_tmp
     
      ! theta
      
      call diffusivity_interp_1d(bin_centres,nbins,theta11
     &   ,ii,y1(n,b),theta11_tmp)
      call diffusivity_interp_1d(bin_centres,nbins,theta12
     &   ,ii,y1(n,b),theta12_tmp)
      call diffusivity_interp_1d(bin_centres,nbins,theta21
     &   ,ii,y2(n,b),theta21_tmp)
      call diffusivity_interp_1d(bin_centres,nbins,theta22
     &   ,ii,y2(n,b),theta22_tmp)
     
      ! random forcing 
      
      b11_tmp = sqrt(2*sigma11_tmp/theta11_tmp)
      b12_tmp = sqrt(2*sigma12_tmp/theta12_tmp)
      b22_tmp = sqrt(2*sigma22_tmp/theta22_tmp)
      b21_tmp = sqrt(2*sigma21_tmp/theta21_tmp)
      
      
      call cubic_poly_x(ii,jj,x1(n,b),y1(n,b)
     & ,a1,b1,c1,d1,psi1_x)
      call cubic_poly_x(ii,jj,x2(n,b),y2(n,b)
     & ,a2,b2,c2,d2,psi2_x)
      
      call vel(ii,jj,psi1_x,a1,b1,c1,d1,x1(n,b),y1(n,b),u1,v1)    
      call vel(ii,jj,psi2_x,a2,b2,c2,d2,x2(n,b),y2(n,b),u2,v2) 
      
      u_top(1) = u1 + u0
      u_top(2) = v1
      u_bottom(1) = u2
      u_bottom(2) = v2
      
      ! calculate drift correction term
      
c        drift11 = .5*d_sigma11_tmp*(u1_fluc_old(1,n,b)
c     &   *(u_top(1) + 2*u1_fluc_old(1,n,b))/sigma11_tmp + 1)
c        drift12 = .5*d_sigma12_tmp*(u1_fluc_old(2,n,b)
c     &   *(u_top(2) + 2*u1_fluc_old(2,n,b))/sigma12_tmp + 1)
c        drift21 = .5*d_sigma21_tmp*(u2_fluc_old(1,n,b)
c     &   *(u_bottom(1) + 2*u2_fluc_old(1,n,b))/sigma21_tmp + 1)
c        drift22 = .5*d_sigma22_tmp*(u2_fluc_old(2,n,b)
c     &   *(u_bottom(2) + 2*u2_fluc_old(2,n,b))/sigma22_tmp + 1)


       drift11 = .5*(1+u1_fluc_old(1,n,b)**2/sigma11_tmp)*d_sigma11_tmp
       drift12 = .5*(1+u1_fluc_old(2,n,b)**2/sigma12_tmp)*d_sigma12_tmp
       drift21 = .5*(1+u2_fluc_old(1,n,b)**2/sigma21_tmp)*d_sigma21_tmp
       drift22 = .5*(1+u2_fluc_old(2,n,b)**2/sigma22_tmp)*d_sigma22_tmp

c       drift11 = d_sigma11_tmp
c       drift12 = d_sigma12_tmp
c       drift21 = d_sigma21_tmp
c       drift22 = d_sigma22_tmp
   
c       drift11 = 0.
c       drift12 = 0.
c       drift21 = 0.
c       drift22 = 0.
      
      ! update velocity fluctuation
      
       u1_fluc_new(1,n,b) = u1_fluc_old(1,n,b) 
     &  + (-u1_fluc_old(1,n,b)/theta11_tmp
     & + drift11)*dt_nondim 
     & + b11_tmp*random_normal()*dsqrt(dt_nondim)  
     
       u1_fluc_new(2,n,b) = u1_fluc_old(2,n,b) 
     &  + (-u1_fluc_old(2,n,b)/theta12_tmp
     & + drift12)*dt_nondim 
     & + b12_tmp*random_normal()*dsqrt(dt_nondim)  
     
       u2_fluc_new(1,n,b) = u2_fluc_old(1,n,b) 
     & + (-u2_fluc_old(1,n,b)/theta21_tmp
     & + drift21)*dt_nondim 
     & + b21_tmp*random_normal()*dsqrt(dt_nondim)  
     
       u2_fluc_new(2,n,b) = u2_fluc_old(2,n,b) 
     &  + (-u2_fluc_old(2,n,b)/theta22_tmp
     & + drift22)*dt_nondim 
     & + b22_tmp*random_normal()*dsqrt(dt_nondim)  
     
      !print*,'u1_fluc =',u1_fluc_new(1,n,b)
      !print*,'u1_mean=',u_top(1)
      
       !print*,'theta,a,b=',theta11_tmp,drift11,b11_tmp
       !print*,'u1_old,u1_new=',u1_fluc_old(1,n,b),u1_fluc_new(1,n,b)
       !print*,'x1_old,y1_old =', x1(n,b),y1(n,b)
       
     

      ! update particle location
      
      x1(n,b) = x1(n,b) + dt_nondim*(u1_fluc_new(1,n,b) + u1+u0)
      x2(n,b) = x2(n,b) + dt_nondim*(u2_fluc_new(1,n,b) + u2)
      y1(n,b) = y1(n,b) + dt_nondim*(u1_fluc_new(2,n,b) + v1)
      y2(n,b) = y2(n,b) + dt_nondim*(u2_fluc_new(2,n,b) + v2)
      
      print*,'x1_new,y1_new =',x1(n,b),y1(n,b)
      
      x1_traj(n,b) = x1_traj(n,b) + dt_nondim*(u1_fluc_new(1,n,b)+u1+u0)
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
            
            if ((x1(n,b).le.0) .or. (x1(n,b) .ge. ii)) then
                print*,'x1(n,b) = ',x1(n,b)
                print*,'u1 = ',u1_fluc_new(1,n,b)
                print*,theta11_tmp,d_sigma11_tmp
     &                ,sigma11_tmp,drift11,b11_tmp,u_top(1)
                stop
            endif
            if ((x2(n,b).le.0) .or. (x2(n,b) .ge. ii)) then
                print*,'x2(n,b) = ',x2(n,b)
                print*,'u2 = ',u2_fluc_new(1,n,b)
                stop
            endif
            if ((y1(n,b).le.0) .or. (y1(n,b) .ge. jj)) then
                print*,'y1(n,b) = ',y1(n,b)
                print*,'v1 = ',u1_fluc_new(2,n,b)
                stop
            endif
            if ((y2(n,b).le.0) .or. (y2(n,b) .ge. jj)) then
                print*,'y2(n,b) = ',y2(n,b)
                print*,'v2 = ',u2_fluc_new(2,n,b)
                stop
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
