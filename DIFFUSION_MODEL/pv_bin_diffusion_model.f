c MODEL THAT RUNS A DIFFUSION MODEL USING EDDY DIFFUSIVITIES 
C CALCULATED FROM LAGRANGIAN TRAJECTORIES

      program diffusion_model
      
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
      use mod_traj_netcdf
      
      implicit none

C READ TIME-AVERAGED STREAM FUNCTION

      call read_ave_file(ave_file,ii,jj,psi1_av,psi2_av,time_av)
      
C CALCULATE COEFFICIENTS FOR 2D-CUBIC SPATIAL INTERPOLATION

      call cubic_coeff_x(ii,jj,psi1_av,a1,b1,c1,d1)
      call cubic_coeff_x(ii,jj,psi2_av,a2,b2,c2,d2)
      
      print*,'cubic_coeff_done'
      
C READ THE DIFFUSIVITY TENSOR

      

      call read_diffusivity(diff_file,K)
      
      print*, 'diffusivity read'
      
      call read_bin_width(bin_file,nbins,bin_width)
      
      print*,'bin_width = ',bin_width
      
c CALCULATE BIN CENTRES

      bin_boundaries(1) = 0.
      do b = 2,nbins+1
        bin_boundaries(b) = bin_boundaries(b-1) + bin_width(b-1)
      enddo

      do b = 1,nbins
        bin_centres(b) = .5*(bin_boundaries(b)+bin_boundaries(b+1))
      enddo
      print*,'bin_boundaries =',bin_boundaries
      print*,'bin_centre=',bin_centres
    
      
c DETERMINE BINS

      do b = 1,nbins+1
      
      bin_corners(b) = (b-1)*(dfloat(jj)/nbins)
      
      enddo

      
c CALCULATE BIN CENTRES

      do b = 1,nbins
      
        uniform_bins(b) = .5*(bin_corners(b)+bin_corners(b+1))
      
      enddo
      


C GENERATE RANDOM LAGRANGIAN PARTICLES

      n = 0
      iseed = 123456789

      do n = 1,npoints
      do b = 1,nbins
      
        x1(n,b) = ran1(iseed)*dfloat(jj)
        y1(n,b) = ran1(iseed)*bin_width(b) + bin_boundaries(b)
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
     

     
c INTERPOLATE THE EDDY-DIFFUSIVITY AT THE PARTICLE LOCATION

                  call diffusivity_interp_1d(uniform_bins,nbins,K(1,:,1)
     &                  ,ii,y1(n,b),Kx1)
                  call diffusivity_interp_1d(uniform_bins,nbins,K(2,:,1)
     &                  ,ii,y1(n,b),Ky1)
                  call diffusivity_interp_1d(uniform_bins,nbins,K(1,:,2)
     &                  ,ii,y2(n,b),Kx2)
                  call diffusivity_interp_1d(uniform_bins,nbins,K(2,:,2)
     &                  ,ii,y2(n,b),Ky2)
     

                  
c APPROXIMATE DERIVATIVE OF THE DIFFUSIVITY - I.E THE DRIFT TERM    
              
                  dKdx1 = 0.
                  dKdx2 = 0.
                  call diffusivity_derivative(uniform_bins,nbins
     &                  ,K(2,:,1),ii,y1(n,b),dKdy1)
                  call diffusivity_derivative(uniform_bins,nbins
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
                    print*,'y1=',y1(n,b)
                    print*,'K = ',K(2,:,1)
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
