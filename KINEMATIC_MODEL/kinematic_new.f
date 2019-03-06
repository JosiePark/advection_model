      program kinematic_model
      
      use mod_random
      use mod_kinematic_netcdf
      use mod_kinematic_advection
      use mod_2dcubic
      use mod_rk4
      
      implicit none
      
      real*8 factor,pi,A,T,c,lambda,max_time,dt,max_time_secs,t_tot,kx
      real*8 y_centre,time_day
      real*8 basinscale,scale,uscale,tscale,dt_nondim
      real*8 time_o,time_dim,time_secs,time_half
      integer ii,jj,k_rec
      parameter(ii = 512, jj = 512)
      parameter(max_time = 100.)
      parameter(dt = 2160.)
      parameter(basinscale = 520.d5)
      
      real*8 psi(ii,jj),x_c(ii),y_c(jj),xnd(ii),ynd(jj),cnd
      real*8 psi_old(ii,jj),psi_half(ii,jj), psi_new(ii,jj)
      real*8 a_new(ii,jj),b_new(ii,jj),c_new(ii,jj),d_new(ii,jj)
      real*8 a_half(ii,jj),b_half(ii,jj),c_half(ii,jj),d_half(ii,jj)
      real*8 a_old(ii,jj),b_old(ii,jj),c_old(ii,jj),d_old(ii,jj)
      real*8 k_o
      real*8 u(ii,jj), v(ii,jj)
      
      integer i,j,n,t_tot_day,iseed
      
      integer npoints,npoints_sqrt
      parameter(npoints_sqrt = 50, npoints = npoints_sqrt**2)
      
      real*8 save_time,eps
      parameter(save_time = 1.) ! save every save_time days
      
      real*8 x0(npoints),y0(npoints),x_traj(npoints),y_traj(npoints)
      real*8 x_diff, y_diff
      
      character*(*),parameter :: home_dir = 
     & '/home/clustor2/ma/j/jp1115/DATA/1/'
      
      character*(*), parameter :: psi_file = 
     &   trim(home_dir) // 'KINEMATIC/psi_eof_1-2.nc'
      character*(*), parameter :: u_file = 
     &   trim(home_dir) // 'KINEMATIC/u_eof_1-2.nc'
      character*(*), parameter :: v_file = 
     &   trim(home_dir) // 'KINEMATIC/v_eof_1-2.nc'
      character*(*), parameter :: traj_file = 
     &   trim(home_dir) // 'KINEMATIC/traj_eof_1_2.nc'
      
      
      
      data pi/3.14159265358979323846D0/
      
      print*,'Using 2d cubic interpolation'
      
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
      call create_kinematic_field(u_file,ii,jj)
      call create_kinematic_field(v_file,ii,jj)
      call create_kinematic_traj(traj_file,npoints)
      
      print*,'created files'
      
        k_o = 0.
            
C ------- MAIN CYCLE - ITERATE IN TIME ------------

      do t = 1,t_tot
        time_o = (t-1)*dt_nondim
        time_dim = time_o*tscale/86400.
        time_secs = (t-1)*dt
        time_day = time_secs/86400.
        
        if (t .eq. 1) then
        k_rec = 1

C ------ RANDOMLY GENERATE PARTICLES -----------

        iseed = 102
        n = 0
        
        do j = 1,npoints_sqrt
        do i = 1,npoints_sqrt
            n = n+1
            x0(n) = ran1(iseed)*dfloat(ii)
            y0(n) = ran1(iseed)*dfloat(jj)

        enddo
        enddo
        x_traj = x0
        y_traj = y0
       
       call write_kinematic_traj(traj_file,npoints
     & ,x0,y0,dfloat(0),k_rec)
       
       
c -------  CONSTRUCT INITIAL KINEMATIC FIELD -----------------

      do i = 1,ii
      do j = 1,jj

       psi_new(i,j) = A/(cosh(factor*(y_c(j)-y_centre)/ii))**2
     &       *cos(kx*(xnd(i)-kx*cnd*time_secs/86400.))
     
       call kinematic_velocity(ii,jj,dfloat(i-1),dfloat(j-1)
     &  ,pi,time_secs,a,factor,cnd,kx,y_centre,u(i,j),v(i,j))
     
      enddo
      enddo
      
    
      
c -------- CONSTRUCT CUBIC INTERPOLATION COEFFICIENTS ---------

      call cubic_coeff_x(ii,jj,psi_new,a_new,b_new,c_new,d_new)
      
      call write_kinematic_field(psi_file,ii,jj,psi_new,k_rec,time_dim)
      call write_kinematic_field(u_file,ii,jj,u,k_rec,time_dim)
      call write_kinematic_field(v_file,ii,jj,v,k_rec,time_dim)
      
      k_rec = k_rec + 1
     
       else
       k_o = k_o + 1
       
c ------ CONSTRUCT NEW AND HALF KINEMATIC FIELDS --------------

      psi_old = psi_new
      
      do i = 1,ii
      do j = 1,jj

       psi_new(i,j) = A/(cosh(factor*(y_c(j)-y_centre)/ii))**2
     &       *cos(kx*(xnd(i)-kx*cnd*time_secs/86400.))
       call kinematic_velocity(ii,jj,dfloat(i-1),dfloat(j-1)
     &  ,pi,time_secs,a,factor,cnd,kx,y_centre,u(i,j),v(i,j))
     
      enddo
      enddo
      
      time_half = time_secs - dt/2.
      
      do i = 1,ii
      do j = 1,jj

       psi_half(i,j) = A/(cosh(factor*(y_c(j)-y_centre)/ii))**2
     &       *cos(kx*(xnd(i)-kx*cnd*time_half/86400.))
     
      enddo
      enddo
      
c -----  FIND CUBIC INTERPOLATION COEFFICIENTS ------------------

      a_old = a_new
      b_old = b_new
      c_old = c_new
      d_old = d_new
      
      call cubic_coeff_x(ii,jj,psi_half,a_half,b_half,c_half,d_half)
      call cubic_coeff_x(ii,jj,psi_new,a_new,b_new,c_new,d_new)


C ------ CALCULATE VELOCITY FIELDS AT THE GRID POINTS BY ANALYTICALLY DIFFERENTIATING ------

      do n = 1,npoints
      
C ------ ADVECT PARTICLES ANALYTICALLY -----

        call rk4_2dcubic(ii,jj,x0(n),y0(n),dt_nondim,dfloat(0)
     & ,a_old,b_old,c_old,d_old
     & ,a_half,b_half,c_half,d_half
     & ,a_new,b_new,c_new,d_new
     & ,x_diff,y_diff)   
        
        x0(n) = x0(n) + x_diff
        y0(n) = y0(n) + y_diff
        
        x_traj(n) = x_traj(n) + x_diff
        y_traj(n) = y_traj(n) + y_diff
        !print*,'x_diff,y_diff = ',x_diff,y_diff
        !print*,'x_traj,y_traj = ',x_traj,y_traj
c ------ TREAT BOUNDARIES --------------

        if (x0(n) .gt. ii) then
            x0(n) = x0(n) - ii
        else if (x0(n) .lt. 0) then
            x0(n) = x0(n) + ii
        endif
        
        if (y0(n) .gt. jj) then
            y0(n) = y0(n) - jj
        else if (y0(n) .lt. 0) then
            y0(n) = y0(n) + jj
        endif
        
        if (x0(n) .gt. ii) then
            print*,'x0=',x0(n)
            exit
        else if (x0(n) .lt. 0) then
            print*,'x0=',x0(n)
            exit
        endif
        
        if (y0(n) .gt. jj) then
            print*,'y0=',y0(n)
            exit
        else if (y0(n) .lt. 0) then
            print*,'y0=',y0(n)
            exit
        endif

      enddo
      

c ------- WRITE TO FILE ------------------

      if (k_o .eq. t_tot_day*save_time) then
        k_o = 0.
        call write_kinematic_traj(traj_file,npoints,x_traj,y_traj
     &   ,time_dim,
     & k_rec)
        call write_kinematic_field(psi_file,ii,jj
     &   ,psi_new,k_rec,time_dim)
        call write_kinematic_field(u_file,ii,jj,u,k_rec,time_dim)
        call write_kinematic_field(v_file,ii,jj,v,k_rec,time_dim)
        
        k_rec = k_rec + 1  
        print*,'file saved at', time_dim,'days' 
      endif 
      endif
      enddo

      end program kinematic_model
