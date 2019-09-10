c MODEL THAT RUNS A DIFFUSION MODEL USING EDDY DIFFUSIVITIES 
C CALCULATED FROM LAGRANGIAN TRAJECTORIES

      program theta_diffusion_model
      
      use mod_diff_input
      use mod_diff_constants
      use mod_diffusivity_functions
      use mod_qg2_netcdf
      use mod_random
      use mod_variables
      use mod_2dcubic
      use mod_diffusion_netcdf
      use mod_1dinterp
      use mod_stochastic_parameters
      
c REWRITE SO THAT THE DERIVED SIFFUSIVTY IS READ FROM A NETCDF FILE THAT
C IS PREVIOUSLY WRITTEN
      
      implicit none
      
      character*(*), parameter :: theta_file = 
     &   trim(home_dir) // 'STATS/THETA/theta.nc'
      character*(*), parameter :: sigma_file = 
     &   trim(home_dir) // 'STATS/SIGMA/velocity_variance.nc'
     
      real*8 theta(2,2,nbins),tmp1(2,2,ii,jj),tmp2(2,2,ii,jj)
      real*8 sigma11(ii,jj),sigma12(ii,jj),sigma21(ii,jj),sigma22(ii,jj)
      real*8 asigma11(ii,jj),bsigma11(ii,jj),csigma11(ii,jj)
     & ,dsigma11(ii,jj)
      real*8 asigma12(ii,jj),bsigma12(ii,jj),csigma12(ii,jj)
     & ,dsigma12(ii,jj)
      real*8 asigma21(ii,jj),bsigma21(ii,jj),csigma21(ii,jj)
     & ,dsigma21(ii,jj)
      real*8 asigma22(ii,jj),bsigma22(ii,jj),csigma22(ii,jj)
     & ,dsigma22(ii,jj)
     
      real*8 K11_tmp,K12_tmp,K21_tmp,K22_tmp

C READ TIME-AVERAGED STREAM FUNCTION

      call read_ave_file(ave_file,ii,jj,psi1_av,psi2_av,time_av)
      
C CALCULATE COEFFICIENTS FOR 2D-CUBIC SPATIAL INTERPOLATION

      call cubic_coeff_x(ii,jj,psi1_av,a1,b1,c1,d1)
      call cubic_coeff_x(ii,jj,psi2_av,a2,b2,c2,d2)
      
      print*,'cubic_coeff_done'
      
C READ THE DIFFUSIVITY TENSOR

      call read_diffusivity(diff_file,K)
      
      print*, 'diffusivity read'
      
c READ THE LAGRANGIAN TIME SCALE
      
      call read_theta_netcdf(theta_file,nbins,theta)
      
      print*,'theta=',theta
      
c Velocity variance (read from file) 
      call read_vel_variance(sigma_file,tmp1,tmp2,ii,jj)
      sigma11 = tmp1(1,1,:,:) ! top layer - zonal
      sigma12 = tmp1(2,2,:,:) ! top layer - meridional
      sigma21 = tmp2(1,1,:,:)
      sigma22 = tmp2(2,2,:,:)
      
      call cubic_coeff_x(ii,jj,sigma11,asigma11
     & ,bsigma11,csigma11,dsigma11)
      call cubic_coeff_x(ii,jj,sigma12,asigma12
     & ,bsigma12,csigma12,dsigma12)
      call cubic_coeff_x(ii,jj,sigma21,asigma21
     & ,bsigma21,csigma21,dsigma21)
      call cubic_coeff_x(ii,jj,sigma22,asigma22
     & ,bsigma22,csigma22,dsigma22)
    
      
c DETERMINE BINS

      do b = 1,nbins+1
      
      bin_corners(b) = (b-1)*(dfloat(jj)/nbins)
      
      enddo

      
c CALCULATE BIN CENTRES

      do b = 1,nbins
      
        bin_centres(b) = .5*(bin_corners(b)+bin_corners(b+1))
      
      enddo
      


C GENERATE RANDOM LAGRANGIAN PARTICLES

      n = 0
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

C MAIN CYCLE

      do t = 2,t_len
            
            time(t) = time(t-1) + dt/86400. ! time in days
            k_s = k_s + dt/86400.
 
            do b = 1,nbins
            do n = 1,npoints

C SPATIALLY INTERPOLATE TO FIND TIME-AVERAGED VELOCITY AT PARTICLE LOCATION


                  call cubic_poly_x(ii,jj,x1(n,b),y1(n,b)
     &                  ,a1,b1,c1,d1,psi1_x)
                  call cubic_poly_x(ii,jj,x2(n,b),y2(n,b)
     &                  ,a2,b2,c2,d2,psi2_x)
     
                  call vel(ii,jj,psi1_x,a1,b1,c1,d1,x1(n,b),y1(n,b)
     &             ,u1,v1)
                  call vel(ii,jj,psi2_x,a2,b2,c2,d2,x2(n,b),y2(n,b)
     &             ,u2,v2)
     
                  !print*,'velocity calculated'
     

     
c INTERPOLATE THETA

      call diffusivity_interp_1d(bin_centres,nbins,theta11
     &   ,ii,y1(n,b),theta11_tmp)
      call diffusivity_interp_1d(bin_centres,nbins,theta12
     &   ,ii,y1(n,b),theta12_tmp)
      call diffusivity_interp_1d(bin_centres,nbins,theta21
     &   ,ii,y2(n,b),theta21_tmp)
      call diffusivity_interp_1d(bin_centres,nbins,theta22
     &   ,ii,y2(n,b),theta22_tmp)

c INTERPOLATE SIGMA

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

c CALCULATE DIFFUSIVITY

      K11_tmp = sigma11_tmp*theta11_tmp
      K12_tmp = sigma12_tmp*theta12_tmp
      K22_tmp = sigma21_tmp*sigma21_tmp
      K22_tmp = sigma22_tmp*sigma22_tmp
                  
c APPROXIMATE DERIVATIVE OF THE DIFFUSIVITY - I.E THE DRIFT TERM    
              
                  dKdx1 = 0.
                  dKdx2 = 0.
                  call diffusivity_derivative(bin_centres,nbins
     &                  ,K(2,:,1),ii,y1(n,b),dKdy1)
                  call diffusivity_derivative(bin_centres,nbins
     &                  ,K(2,:,2),ii,y2(n,b),dKdy2)
     

c SCALE DIFFUSIVTY AND ITS DERIVATIVE

                 Kx1 = Kx1*(tscale/86400.)/((scale/(10**5))**2)
                 Ky1 = Ky1*(tscale/86400.)/((scale/(10**5))**2)
                 Kx2 = Kx2*(tscale/86400.)/((scale/(10**5))**2)
                 Ky2 = Ky2*(tscale/86400.)/((scale/(10**5))**2)
                 
                 dKdy1 = dKdy1*(tscale/86400.)/((scale/(10**5)))
                 dKdy2 = dKdy2*(tscale/86400.)/((scale/(10**5)))
     
C ADVECT PARTICLES
      
                  x1_diff = (u1+U_0+dKdx1)*dt_nondim + 
     &                  dsqrt(2*Kx1)*random_normal()*dsqrt(dt_nondim)
                  if(isnan(x1_diff)) then
                    print*,'dKdx1=',dKdx1
                    print*,'Kx1 =',Kx1
                    stop
                  endif
                  y1_diff = (v1+dKdy1)*dt_nondim + 
     &                  dsqrt(2*Ky1)*random_normal()*dsqrt(dt_nondim)
                  if(isnan(y1_diff)) then
                    print*,'dKdy1=',dKdy1
                    print*,'Ky1 =',Ky1
                    stop
                  endif
                  
                  x2_diff = (u2+dKdx2)*dt_nondim + 
     &                  dsqrt(2*Kx2)*random_normal()*dsqrt(dt_nondim)
                  if(isnan(x2_diff)) then
                    print*,'dKdx2=',dKdx2
                    print*,'Kx2 =',Kx2
                    stop
                  endif
                  y2_diff =  (v2+dKdy2)*dt_nondim + 
     &                  dsqrt(2*Ky2)*random_normal()*dsqrt(dt_nondim)
                  if(isnan(y2_diff)) then
                    print*,'dKdy2=',dKdy2
                    print*,'Ky2 =',Ky2
                    print*,'y_diff=',y2_diff
                    print*,'v2=',v2
                    print*,'y2 = ',y2(n,b)
                    stop
                  endif
     
                  x1(n,b) = x1(n,b) + x1_diff
                  x2(n,b) = x2(n,b) + x2_diff
                  y1(n,b) = y1(n,b) + y1_diff
                  y2(n,b) = y2(n,b) + y2_diff
                  
                  x1_traj(n,b) = x1_traj(n,b) + x1_diff
                  x2_traj(n,b) = x2_traj(n,b) + x2_diff
                  y1_traj(n,b) = y1_traj(n,b) + y1_diff
                  y2_traj(n,b) = y2_traj(n,b) + y2_diff
                  
      
      
c APPLY DOUBLY PERIODIC BOUNDARY CONDITIONS


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
            
            

            enddo
            enddo
            
      


C WRITE TO NETCDF FILE
      
      if(k_s + 0.001 .ge. k_save) then
      
        print*,'Saving at time',time(t)
        
        call write_diffusion_trajfile(file_name,npoints,nbins
     & ,x1_traj,y1_traj,x2_traj,y2_traj,time(t),nrec)
     
      nrec = nrec + 1
      
      k_s = 0.
      
      endif
            
      enddo
      
      end program
