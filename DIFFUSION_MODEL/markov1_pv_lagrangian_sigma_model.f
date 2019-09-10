      program markov1_pv_lagrangian_markov1_model
      
      use mod_stochastic_parameters
      use mod_vel_variance_netcdf
      use mod_theta_netcdf
      use mod_2dcubic
      use mod_qg2_netcdf
      use mod_random
      use mod_diffusion_netcdf
      use mod_traj_netcdf
      
      implicit none
      
      character*(*), parameter :: home_dir = 
     & '/home/clustor2/ma/j/jp1115/DATA/2/'                                      /'
     
      character*(*), parameter :: ave_file = 
     &   trim(home_dir) // 'QG/QG_ave.nc'
      character*(*), parameter :: pv_theta_file = 
     &   trim(home_dir) // 'STATS/THETA/markov1_theta_PV.nc'
      character*(*), parameter :: theta_file = 
     &   trim(home_dir) // 'STATS/THETA/pseudo_theta_osc.nc'
      character*(*), parameter :: pv_sigma_file = 
     &   trim(home_dir) // 'STATS/SIGMA/new_pv_lagrangian_sigma.nc'
      character*(*), parameter :: sigma_file = 
     &   trim(home_dir) // 'STATS/SIGMA/new_lagrangian_sigma.nc'
      character*(*), parameter :: file_name =
     &   trim(home_dir) // 
     &  'TRAJ/FINAL_MARKOV1/pv_bin_pv_diff_spd_T_markov1.nc' 
       character*(*),parameter :: bin_file =
     & trim(home_dir) // 'TRAJ/PV_BINS/test_bin_width.nc'
       character*(*),parameter :: diffusivity_file =
     & trim(home_dir) // 'STATS/DIFFUSIVITY/
     &test_full_PVDISP_MEAN_DIFF.nc' 
       character*(*), parameter :: diff_file =
     &   trim(home_dir) // 'STATS/DIFFUSIVITY/pseudo_new_DIFF.nc' 
     
      integer ii,jj,nbins,npoints
      parameter(ii = 512,jj = 512,nbins = 10,npoints=1000)
      real*8 psi1_av(ii,jj),psi2_av(ii,jj)
      real*8 a1(ii,jj),b1(ii,jj),c1(ii,jj),d1(ii,jj)
      real*8 a2(ii,jj),b2(ii,jj),c2(ii,jj),d2(ii,jj)
      real*8 tmp1(2,2,ii,jj),tmp2(2,2,ii,jj)
      real*8 pv_sigma(nbins)
      real*8, allocatable, dimension(:,:,:,:) :: sigma
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
      real*8 theta_pv(nbins),theta_pv_tmp, d_theta_pv_tmp
      real*8 time_av
      real*8 d_sigma12(ii,jj),d_sigma11(ii,jj)
     & ,d_sigma22(ii,jj),d_sigma21(ii,jj) ! derivatives of sigma
      real*8 d_theta12_tmp, d_K12_tmp
      
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
      real*8 K_tmp,KPV_tmp
      
      real*8,allocatable,dimension(:) :: KPV
      real*8,allocatable,dimension(:,:,:) :: K_diff
      
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
      
      real*8 bin_width(nbins), bin_boundaries(nbins+1)
      real*8 uniform_bins(nbins)
      
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
      
C NON DIMENSIONALISE TIME

      scale = basinscale/dfloat(ii)
      uscale = 1
      tscale = scale/uscale

      dt_nondim = dt/tscale
    
c  FIND BIN BOUNDARIES

      call read_bin_width(bin_file,nbins,bin_width)
      
c      print*,'bin_width = ',bin_width
      
c CALCULATE BIN CENTRES

      bin_boundaries(1) = 0.
      do b = 2,nbins+1
        bin_boundaries(b) = bin_boundaries(b-1) + bin_width(b-1)
      enddo

      do b = 1,nbins
        bin_centres(b) = .5*(bin_boundaries(b)+bin_boundaries(b+1))
      enddo
c      print*,'bin_boundaries =',bin_boundaries
      print*,'bin_centre=',bin_centres
      
c DETERMINE BINS

      do b = 1,nbins+1
      
      bin_corners(b) = (b-1)*(dfloat(jj)/nbins)
      
      enddo
      
c      print*,'bin_corners=',bin_corners

      
c CALCULATE BIN CENTRES

      do b = 1,nbins
      
        uniform_bins(b) = .5*(bin_corners(b)+bin_corners(b+1))
      
      enddo
      
c      print*,'uniform_bins = ',uniform_bins
      
C READ TIME-AVERAGED STREAM FUNCTION

      call read_ave_file(ave_file,ii,jj,psi1_av,psi2_av,time_av)
      
C CALCULATE COEFFICIENTS FOR 2D-CUBIC SPATIAL INTERPOLATION

      call cpu_time(start_time)

      call cubic_coeff_x(ii,jj,psi1_av,a1,b1,c1,d1)
      
      call cpu_time(stop_time)
      
c      print*,'calculating coefficients =',stop_time - start_time
      
      
c CALCULATE COEFFICIENTS FOR MARKOV-1 MODEL

c 1. Markov-1 time tensor (defined for each bin, read from file)
      call read_pv_theta_netcdf(pv_theta_file,nbins,theta12)
      !theta12 = theta12*86400./tscale
      theta_pv = theta12*86400./tscale
      
      call read_theta_netcdf(theta_file,nbins,theta)
      theta11 = theta(1,1,:)*86400./tscale
      theta12 = theta(1,2,:)*86400./tscale
      !print*,'theta12 = ',theta(1,2,:)
      
      !print*, 'theta12 = ',theta12

      
c 2. Velocity variance (read from file) 
      call read_lagrangian_sigma(sigma_file,sigma)
      
      print*,'sigma read'
      
c read PV lagrangian velocity variance

      call read_pv_lagrangian_sigma(pv_sigma_file,nbins,pv_sigma)
      
      print*,'pv vel variance read'
      
c 3. PV diffusivity coefficient

      call read_PV_diffusivity(diffusivity_file,KPV)
      call read_diffusivity(diff_file,K_diff)
      
      print*,' diffusivity read'

      
c SCALE THE DIFFUSIVITY

      KPV = KPV*(tscale/86400.)/((scale/(10**5))**2)
      K_diff = K_diff*(tscale/86400.)/((scale/(10**5))**2)
      
      print*,'K_PV = ',KPV
      
      
C GENERATE RANDOM LAGRANGIAN PARTICLES

      iseed = 123456789

      do n = 1,npoints
      do b = 1,nbins
      
        x1(n,b) = ran1(iseed)*dfloat(jj)
        !y1(n,b) = ran1(iseed)*dfloat(jj)/nbins + bin_corners(b)
        y1(n,b) = ran1(iseed)*bin_width(b) + bin_boundaries(b)
        
        x1_traj(n,b) = x1(n,b)
        y1_traj(n,b) = y1(n,b)

      
      enddo
      enddo
      
      
c WRITE INITIAL POSITIONS TO FILE

      call create_pvdiffusion_trajfile(file_name,npoints,nbins)
      
      nrec = 1 
      
      call write_pvdiffusion_trajfile(file_name,npoints,nbins
     & ,y1,dfloat(0),nrec)
     
      nrec = nrec + 1
      
      
c DETERMINE TIME ARRAY

      t_len  = int((max_run*86400.)/dt) ! number of time steps
      allocate(time(t_len))      
    
      time(1) = 0.
      k_s = 0.
      
c DETERMINE INITIAL CONDITIONS, 
C u' IS INITIALISED AS A GAUSSIAN RANDOM VARIABLE WITH ZERO MEAN AND VARIANCE SIGMA MEAN

c      sum_tmp = 0.
c      do i = 1,ii
c      do j = 1,jj
c        sum_tmp = sum_tmp + sigma11(i,j)
c      enddo
c      enddo
c      sigma11_mean = sum_tmp/(ii*jj)
      
c      sum_tmp = 0.
c      do i = 1,ii
c      do j = 1,jj
c        sum_tmp = sum_tmp + sigma12(i,j)
c      enddo
c      enddo
c      sigma12_mean = sum_tmp/(ii*jj)
      
      do n = 1,npoints
      do b = 1,nbins
        u1_fluc_old(1,n,b) = random_normal()*dsqrt(sigma(1,1,1,b))
        u1_fluc_old(2,n,b) = random_normal()*dsqrt(pv_sigma(b))       
      enddo
      enddo
      

C MAIN CYCLE

      do t = 2,t_len

            time(t) = time(t-1) + dt/86400. ! time in days
            k_s = k_s + dt/86400.
      
      do b = 1,nbins
      do n = 1,npoints
c INTERPOLATE PARAMETERS TO FIND VALUES AT PARTICLE LOCATIONS
C DETERMINE DRIFT CORRECTION TERM

      
      ! sigma and dsigma
       !print*,' bin=',b
       !print*,'y1(n,b)=',y1(n,b)
       
      
c      call diffusivity_interp_1d(uniform_bins,nbins,sigma(1,1,1,:)
c     &   ,ii,y1(n,b),sigma11_tmp)
c      call diffusivity_interp_1d(bin_centres,nbins,pv_sigma
c     &   ,ii,y1(n,b),sigma12_tmp)
     
c      d_sigma11_tmp = 0.
c      call diffusivity_derivative(bin_centres,nbins
c     &                  ,pv_sigma,ii,y1(n,b),d_sigma12_tmp)
     
      
     
      ! theta
      
      call diffusivity_interp_1d(uniform_bins,nbins,theta11
     &   ,ii,y1(n,b),theta11_tmp)
     
      !print*,'interpolating theta12'c
c      call diffusivity_interp_1d(bin_centres,nbins,theta12
c     &   ,ii,y1(n,b),theta12_tmp)
     
      !print*,'theta12_tmp = ',theta12_tmp
      !if (theta12_tmp .le. 1.D-3) then
      !  theta12_tmp =  3.
      !endif
      !print*,'theta12_tmp = ',theta12_tmp
      call diffusivity_interp_1d(uniform_bins,nbins,theta12
     &   ,ii,y1(n,b),theta12_tmp)
     
      call diffusivity_interp_1d(bin_centres,nbins,theta_pv
     &   ,ii,y1(n,b),theta_pv_tmp)
     
      
     
      

      
      ! KPV
      
      call diffusivity_interp_1d(bin_centres,nbins,KPV
     & , ii,y1(n,b),KPV_tmp) 
     
      call diffusivity_interp_1d(uniform_bins,nbins,K_diff(1,:,1)
     & , ii,y1(n,b),K_tmp) 
     
      call diffusivity_derivative(bin_centres,nbins
     &                  ,KPV,ii,y1(n,b),d_K12_tmp)
    
c      call diffusivity_derivative(bin_centres,nbins
c     &                  ,theta12,ii,y1(n,b),d_theta12_tmp)
     
            call diffusivity_derivative(bin_centres,nbins
     &                  ,theta_pv,ii,y1(n,b),d_theta_pv_tmp)

            call diffusivity_derivative(uniform_bins,nbins
     &                  ,theta12,ii,y1(n,b),d_theta12_tmp)
     
     
cc      ! random forcing 
      
      sigma11_tmp = K_tmp/theta11_tmp
      sigma12_tmp = KPV_tmp/theta12_tmp

      
      d_sigma11_tmp = 0.
      d_sigma12_tmp = (theta12_tmp*d_K12_tmp - KPV_tmp*d_theta12_tmp)
     & /(theta12_tmp**2)  
      
    
      b11_tmp = sqrt(2*sigma11_tmp/theta11_tmp)
      b12_tmp = sqrt(2*sigma12_tmp/theta12_tmp) 
      
      
      call cubic_poly_x(ii,jj,x1(n,b),y1(n,b)
     & ,a1,b1,c1,d1,psi1_x)
      
      call vel(ii,jj,psi1_x,a1,b1,c1,d1,x1(n,b),y1(n,b),u1,v1)    
      
      u_top(1) = u1 + u0
      u_top(2) = v1
      
      ! calculate drift correction term
      
c       drift11 = d_sigma11_tmp
c       drift12 = d_sigma12_tmp
      
       drift11 = .5*(1+u1_fluc_old(1,n,b)**2/sigma11_tmp)*d_sigma11_tmp
       drift12 = .5*(1+u1_fluc_old(2,n,b)**2/sigma12_tmp)*d_sigma12_tmp
 
      ! update velocity fluctuation
      
c       theta11_tmp = K_tmp/sigma11_tmp
c       theta12_tmp = KPV_tmp/sigma12_tmp
      
       u1_fluc_new(1,n,b) = u1_fluc_old(1,n,b) 
     &  + (-u1_fluc_old(1,n,b)/theta11_tmp
     & + drift11)*dt_nondim 
     & + b11_tmp*random_normal()*dsqrt(dt_nondim)  
     
       u1_fluc_new(2,n,b) = u1_fluc_old(2,n,b) 
     &  + (-u1_fluc_old(2,n,b)/theta12_tmp
     & + drift12)*dt_nondim 
     & + b12_tmp*random_normal()*dsqrt(dt_nondim)  
     

      ! update particle location
      
      x1(n,b) = x1(n,b) + dt_nondim*(u1_fluc_new(1,n,b) + u1+u0)
      y1(n,b) = y1(n,b) + dt_nondim*(u1_fluc_new(2,n,b) + v1)
      
      x1_traj(n,b) = x1_traj(n,b) + dt_nondim*(u1_fluc_new(1,n,b)+u1+u0)
      y1_traj(n,b) = y1_traj(n,b) + dt_nondim*(u1_fluc_new(2,n,b) + v1)
      
      ! treat boundaries
      
            if(x1(n,b).le.0) then
                x1(n,b) = x1(n,b) + dfloat(ii)
            endif
            if(x1(n,b).ge.ii) then
                x1(n,b) = x1(n,b) - dfloat(ii)
            endif
            
            if(y1(n,b).le.0) then
                y1(n,b) = y1(n,b) + dfloat(jj)
            endif
            if(y1(n,b).ge.jj) then
                y1(n,b) = y1(n,b) - dfloat(jj)
            endif
            
            if ((x1(n,b).le.0) .or. (x1(n,b) .ge. ii)) then
                print*,'x1(n,b) = ',x1(n,b)
                print*,'u1 = ',u1_fluc_new(1,n,b)
                print*,theta11_tmp,d_sigma11_tmp
     &                ,sigma11_tmp,drift11,b11_tmp,u_top(1)
                stop
            endif
            if ((y1(n,b).le.0) .or. (y1(n,b) .ge. jj)) then
                print*,'y1(n,b) = ',y1(n,b)
                print*,'v1 = ',u1_fluc_new(2,n,b)
                print*,'theta=',theta12_tmp
                print*,'K =',KPV_tmp
                print*,'Kx=',K_tmp
                print*,'sigma=',sigma12_tmp
                stop
            endif
            
            ! update velocoty fluctuations
            
      do d = 1,2
        u1_fluc_old(d,n,b) = u1_fluc_new(d,n,b)
      enddo
      
      enddo
      enddo
      
      ! save
      
      if(k_s + 0.001 .ge. k_save) then
      
        print*,'Saving at time',time(t)
        
        call write_pvdiffusion_trajfile(file_name,npoints,nbins
     & ,y1_traj,time(t),nrec)
     
      nrec = nrec + 1
      
      k_s = 0.
      
      endif
            
      enddo
      
      end program 
