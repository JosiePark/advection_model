      program kinematic_test
      
      use mod_random
      use mod_kinematic_netcdf
      use mod_kinematic_advection
      
      implicit none
      
      real*8 factor,pi,A,T,c,lambda,max_time,dt,max_time_secs,t_tot,kx
      real*8 y_centre
      real*8 basinscale,scale,uscale,tscale,dt_nondim
      real*8 time_o,time_dim,time_secs
      integer ii,jj,k_rec
      parameter(ii = 512, jj = 512)
      parameter(max_time = 1000.)
      parameter(dt = 2160.)
      parameter(basinscale = 520.d5)
      
      integer i,j,n,t_tot_day,iseed
      
      integer npoints,npoints_sqrt
      parameter(npoints_sqrt = 500, npoints = npoints_sqrt**2)
      
      real*8 save_time,eps
      parameter(save_time = 1.) ! save every save_time days
      
      real*8 x0(npoints),y0(npoints),x_traj(npoints),y_traj(npoints)
      real*8 x_diff, y_diff
      
      character*(*),parameter :: home_dir = 
     & '/home/clustor2/ma/j/jp1115/DATA/1/'
      
      character*(*), parameter :: psi_file = 
     &   trim(home_dir) // 'KINEMATIC/psi_eof_1-2.nc'
      character*(*), parameter :: traj_file = 
     &   trim(home_dir) // 'KINEMATIC/traj_eof_1_2.nc'
      
      real*8 psi(ii,jj),x_c(ii),y_c(jj),xnd(ii),ynd(jj),cnd
      real*8 k_o
      data pi/3.14159265358979323846D0/
      
      scale = basinscale/dfloat(ii)
      uscale = 1.
      tscale = scale/uscale
      
    
C ------- KINEMATIC MODEL PARAMETERS --------------

      factor = 3.*pi ! width scaling
      kx = 2 ! wave number
      A = 526.1 ! amplitude of wave
      T = 48. ! period in days
      lambda = float(ii)/(2*kx) ! wavelength
      c = lambda/T ! propagation speed
      y_centre = 302. ! centre in grid points
      cnd = c*2*pi/ii - pi
    
c ------- DEFINE SPATIAL VARIABLES ----------------

      do i = 1,ii
        x_c(i) = dfloat(i-1)
      enddo
      xnd = x_c*2*pi/ii - pi
      
      do j = 1,jj
        y_c(j) = dfloat(j-1)
      enddo
      ynd = y_c*4*pi/jj - 2*pi
      
      
c -------- TIME PARAMETERS -------------

      max_time_secs = max_time*86400
      t_tot = int(max_time_secs/dt) ! total numer of time steps
      t_tot_day = int(86400/dt) ! total number of time steps per day
      dt_nondim = dt/tscale
      
c -------- CREATE FILES ------------------

      call create_kinematic_field(psi_file,ii,jj)
      call create_kinematic_traj(traj_file,npoints)
      print*,'created files'
      
        k_o = 0.
            
C ------- MAIN CYCLE - ITERATE IN TIME ------------

      do t = 1,t_tot
        time_o = (t-1)*dt_nondim
        time_dim = time_o*tscale/86400.
        time_secs = (t-1)*dt
        print*,'Time = ',time_dim
        print*,'k_o = ',k_o/86400
C ------ CONSTRUCT KINEMATIC FIELD --------
        do i = 1,ii
        do j = 1,jj
            psi(i,j) = A/(cosh(factor*(y_c(j)-y_centre)/ii))**2
     &       *cos(kx*(xnd(i)-kx*cnd*time_secs/86400.))
        enddo
        enddo
        
        if (t .eq. 1) then
        k_rec = 1

       call write_kinematic_field(psi_file,ii,jj,psi,k_rec,time_dim)
       k_rec = k_rec + 1
       else
       k_o = k_o + 1
C ------ CALCULATE VELOCITY FIELDS AT THE GRID POINTS BY ANALYTICALLY DIFFERENTIATING ------

      

c ------- WRITE TO FILE ------------------

      if (k_o .eq. t_tot_day*save_time) then
        k_o = 0.
        call write_kinematic_field(psi_file,ii,jj,psi,k_rec,time_dim)
        k_rec = k_rec + 1  
        print*,'file saved at', time_dim,'days' 
      endif 
      endif
      enddo

      end program kinematic_test
