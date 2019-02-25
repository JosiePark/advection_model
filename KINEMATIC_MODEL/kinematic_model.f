      program kinematic_model
      
      use mod_random
      
      implicit none
      
      real*8 factor,pi,A,T,c,lambda,max_time,dt,max_time_secs,t_tot
      real*8 centre
      real*8 basinscale,scale,uscale,tscale,dt_nondim
      real*8 time_o,time_dim
      integer kx,ii,jj
      parameter(ii = 512, jj = 512)
      parameter(max_time = 1000.)
      parameter(dt = 2160.)
      parameter(basinscale = 520.d5)
      
      real*8 psi(ii,jj)
      data pi/3.14159265358979323846D0/
      
      scale = basinscale/dfloat(ii)
      uscale = 1.
      tscale = scale/uscale
      
    
C ------- KINEMATIC MODEL PARAMETERS --------------

      factor = 3.*pi ! width scaling
      kx = 2 ! wave number
      A = 526.1 ! amplitude of wave
      T = 48.*24.*60.*60. ! period in secs
      lambda = float(ii)/kx ! wavelength
      c = lambda/T ! propagation speed
      centre = 302 ! centre in grid points
      
c -------- TIME PARAMETERS -------------

      max_time_secs = max_time*86400
      t_tot = int(max_time_secs/dt) ! total numer of time steps
      t_tot_day = int(86400/dt) ! total number of time steps per day
      dt_nondim = dt/tscale
      
      
C ------ RANDOMLY GENERATE PARTICLES -----------

        iseed = 102
        
        nrec = 0
        do j = 1,npoints_sqrt
        do i = 1,npoints_sqrt
            n = n+1
            x0(n) = ran1(iseed)*dfloat(ii)
            y0(n) = ran1(iseed)*dfloat(jj)

        enddo
        enddo
            

C ------- MAIN CYCLE - ITERATE IN TIME ------------

      do t = 1,t_tot
        time_o = (t-1)*dt_nondim
        time_dim = time_o*tscale/86400.
C ------ CONSTRUCT KINEMATIC FIELD --------
        do i = 1,ii
            do j = 1,jj
                psi(i,j) = A/(cosh(factor*(y-y_centre)**2*cos(kx*(xnd-ct)))
            enddo
        enddo
        
      open(11,file = 'kinematic_psi.dat',form = 'unformatted')
      write(11) psi
      close(11)
C ------ CALCULATE VELOCITY FIELDS AT THE GRID POINTS BY ANALYTICALLY DIFFERENTIATING ------

      do n = 1,npoints
      
C ------ ADVECT PARTICLES ANALYTICALLY -----

      enddo
      enddo

      end program kinematic_model
