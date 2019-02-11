      program markov1_model
      
      use stochastic_parameters
      
      implicit none
      
C READ TIME-AVERAGED STREAM FUNCTION

      call read_ave_file(ave_file,ii,jj,psi1_av,psi2_av,time_av)
      
C CALCULATE COEFFICIENTS FOR 2D-CUBIC SPATIAL INTERPOLATION

      call cubic_coeff_x(ii,jj,psi1_av,a1,b1,c1,d1)
      call cubic_coeff_x(ii,jj,psi2_av,a2,b2,c2,d2)
      
c CALCULATE COEFFICIENTS FOR MARKOV-1 MODEL
      
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
      
c CALCULATE NECESSARY STOCHASTIC PARAMETERS for previous and current :

c. 1. autocorrelation


      
      enddo
      
      end program markov1_model
