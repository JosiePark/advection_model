c program that runs a diffusion model for the pv mapped dispersion

      program pv_diffusion_model
      
      use mod_diff_input
      use mod_random
      use mod_qg2_netcdf
      use mod_2dcubic
      use mod_traj_netcdf
      use mod_diffusivity_functions
      use mod_1dinterp
      use mod_diffusion_netcdf
      use mod_diff_constants
      
      implicit none
      
c DEFINE DIMENSIONAL VARIABLES

      uscale=1.
      scale=basinscale/dfloat(ii) 
      tscale=scale/uscale
      SS=(scale/Rd)**2
      
      dt_nondim = dt/tscale

      S1=SS/(1.+H1/H2)
      S2=SS-S1
      beta_nondim=beta*scale*scale/uscale
      
      BETA_NONDIM_U1=BETA_NONDIM+U_0*S1
      BETA_NONDIM_U2=BETA_NONDIM-U_0*S2
        
        do j = 1,jj1
            y_c(j) = dfloat(j-1)
        enddo
        
c READ TIME-AVERAGED STREAM FUNCTION

      call read_ave_file(ave_file,ii,jj,psi1_av,psi2_av,time_av)
      
C CALCULATE COEFFICIENTS FOR 2D-CUBIC SPATIAL INTERPOLATION

      call cubic_coeff_x(ii,jj,psi1_av,a1,b1,c1,d1)
      call cubic_coeff_x(ii,jj,psi2_av,a2,b2,c2,d2)

      print*,'cubic coefficients read'
c  FIND BIN BOUNDARIES

      call read_bin_width(bin_file,nbins,bin_width)
      
      print*,'bin_width = ',bin_width
      
c CALCULATE BIN CENTRES

      bin_boundaries(1) = 0.
      do b = 2,nbins+1
        bin_boundaries(b) = bin_boundaries(b-1) + bin_width(b)
      enddo

      do b = 1,nbins
        bin_centres(b) = .5*(bin_boundaries(b)+bin_boundaries(b+1))
      enddo
      
c DETERMINE BINS

      do b = 1,nbins+1
      
      bin_corners(b) = (b-1)*(dfloat(jj)/nbins)
      
      enddo

      
c CALCULATE BIN CENTRES

      do b = 1,nbins
      
        uniform_bins(b) = .5*(bin_corners(b)+bin_corners(b+1))
      
      enddo
      
C READ THE DIFFUSIVITY TENSOR

      call read_PV_diffusivity(pvdiff_file,KPV)
      
      print*, 'diffusivity read'
      
      call read_diffusivity(diff_file,K)
      
      print*, 'full diffusivity read'

c  SEED PARTICLES UNIFORMALLY IN EACH BIN

      n = 0
      iseed = 123456789
      
      do n = 1,npoints
      do b = 1,nbins
      
        y1(n,b) = ran1(iseed)*bin_width(b)
     &             + bin_boundaries(b)
        x1(n,b) = ran1(iseed)*dfloat(ii)
      
      enddo
      enddo
      
      print*,'particles seeded'
      
      call create_pvdiffusion_trajfile(file_name,npoints,nbins)
      
      nrec = 1 
      
      call write_pvdiffusion_trajfile(file_name,npoints,nbins
     & ,y1,dfloat(0),nrec)
     
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
            print*,'time = ',time(t)
            print*,'k_s =', k_s
      
            do b = 1,nbins
            do n = 1,npoints
            
            !! need to calculate a zonally-averaged time-averaged stream function
            
            
            
      
C SPATIALLY INTERPOLATE TO FIND TIME-AVERAGED VELOCITY AT PARTICLE LOCATION


                  call cubic_poly_x(ii,jj,x1(n,b),y1(n,b)
     &                  ,a1,b1,c1,d1,psi1_x)
     
                  call vel(ii,jj,psi1_x,a1,b1,c1,d1,x1(n,b),y1(n,b)
     &             ,u1,v1)
     

     
c INTERPOLATE THE EDDY-DIFFUSIVITY AT THE PARTICLE LOCATION

                  call diffusivity_interp_1d(uniform_bins,nbins,K(1,:,1)
     &                  ,ii,y1(n,b),Kx1)
                  call diffusivity_interp_1d(bin_centres,nbins,KPV
     &                  ,ii,y1(n,b),Ky1)
                  

                  
c APPROXIMATE DERIVATIVE OF THE DIFFUSIVITY - I.E THE DRIFT TERM    
              
                  dKdx1 = 0.
                  call diffusivity_derivative(bin_centres,nbins
     &                  ,KPV,ii,y1(n,b),dKdy1)

     
                
c DETERMINE WHICH EDDY DIFFUSIVITY CONSTANT TO DIFFUSE PARTICLE WITH

c                do bb = 1,nbins
c                    if (y1(n,bb) .le. bin_corners(bb+1)) then
c                        Kx1 = K(1,bb,1)
c                        Ky1 = K(2,bb,1)
c                        !print*,'b=',bb
c                        goto 10
c                    endif
c                enddo
c10      continue
c                do bb = 1,nbins
c                    if (y2(n,bb) .le. bin_corners(bb+1)) then
c                        Kx2 = K(1,b,2)
c                        Ky2 = K(2,b,2)
c                        !print*,'b=',bb
c                        goto 11
c                    endif
c                enddo
c11      continue

C ADVECT PARTICLES
      
                  x1_diff = (u1+U_0+dKdx1)*dt_nondim + 
     &                  dsqrt(2*Kx1)*random_normal()*dt_nondim
                  y1_diff = (v1+dKdy1)*dt_nondim + 
     &                  dsqrt(2*Ky1)*random_normal()*dt_nondim
                  
     
                  x1(n,b) = x1(n,b) + x1_diff
                  y1(n,b) = y1(n,b) + y1_diff
                  
                  x1_traj(n,b) = x1_traj(n,b) + x1_diff
                  y1_traj(n,b) = y1_traj(n,b) + y1_diff
      
      
c APPLY DOUBLY PERIODIC BOUNDARY CONDITIONS


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

            
            

            enddo
            enddo
            
      


C WRITE TO NETCDF FILE
      
      if(k_s + 0.001 .ge. k_save) then
      
        print*,'Saving at time',time(t)
        
        call write_pvdiffusion_trajfile(file_name,npoints,nbins
     & ,y1_traj,time(t),nrec)
     
      nrec = nrec + 1
      
      k_s = 0.
      
      endif
            
      enddo

      end program pv_diffusion_model
