      program kinematic_model
      
      use mod_random
      use mod_kinematic_netcdf
      use mod_kinematic_advection
      
      implicit none
      
      real*8 factor,pi,A,T,c,lambda,max_time,dt,max_time_secs,t_tot,kx
      real*8 y_centre,time_day
      real*8 basinscale,scale,uscale,tscale,dt_nondim
      real*8 time_o,time_dim,time_secs
      integer ii,jj,k_rec
      parameter(ii = 512, jj = 512)
      parameter(max_time = 1000.)
      parameter(dt = 2160.)
      parameter(basinscale = 520.d5)
      
      integer i,j,n,t_tot_day,iseed
      
      integer npoints,npoints_sqrt
      parameter(npoints_sqrt = 50, npoints = npoints_sqrt**2)
      
      real*8 save_time,eps
      parameter(save_time = 1.) ! save every save_time days
      
      real*8 x0(npoints),y0(npoints),x_traj(npoints),y_traj(npoints)
      real*8 x_diff(npoints), y_diff(npoints)
      real*8 u(ii,jj), v(ii,jj)
      
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
      character*(*), parameter :: diff_file = 
     &   trim(home_dir) // 'KINEMATIC/diff_eof_1_2.nc' 
      
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
      lambda = float(ii)/(kx) ! wavelength
      c = lambda/T ! propagation speed
      y_centre = 302. ! centre in grid points
      cnd = c*2*pi/ii - pi
      !cnd = 0.
      !c = 0.
    
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
      print*,'dt = ',dt
      print*,'dt_nondim = ',dt_nondim
      
c -------- CREATE FILES ------------------

      call create_kinematic_field(psi_file,ii,jj)
      call create_kinematic_field(u_file,ii,jj)
      call create_kinematic_field(v_file,ii,jj)
      call create_kinematic_traj(traj_file,npoints)
      call create_kinematic_traj(diff_file,npoints)
      print*,'created files'
      
        k_o = 0.
            
C ------- MAIN CYCLE - ITERATE IN TIME ------------

      do t = 1,t_tot
        time_o = (t-1)*dt_nondim
        time_dim = time_o*tscale/86400.
        time_secs = (t-1)*dt
        time_day = time_secs/86400.
        print*,'Time = ',time_dim
        print*,'k_o = ',k_o/86400
C ------ CONSTRUCT KINEMATIC FIELD --------
        do i = 1,ii
        do j = 1,jj
            psi(i,j) = A/(cosh(factor*(y_c(j)-y_centre)/ii))**2
     &       *cos(kx*(xnd(i)-cnd*time_secs/86400.))
     
            call kinematic_velocity(ii,jj,dfloat(i-1),dfloat(j-1)
     &  ,pi,time_secs,a,factor,cnd,kx,y_centre,u(i,j),v(i,j))
        enddo
        enddo
        
        if (t .eq. 1) then
        k_rec = 1

C ------ RANDOMLY GENERATE PARTICLES -----------

        iseed = 102

        do n = 1,npoints
            x0(n) = ran1(iseed)*dfloat(ii)
            y0(n) = ran1(iseed)*dfloat(jj)
            !y0(n) = 310.
            x_traj(n) = x0(n)
            y_traj(n) = y0(n)
            x_diff(n) = 0.
            y_diff(n) = 0.

        enddo
       call write_kinematic_traj(traj_file,npoints
     & ,x0,y0,dfloat(0),k_rec)
       call write_kinematic_traj(diff_file,npoints
     & ,x_diff,y_diff,dfloat(0),k_rec)
       call write_kinematic_field(psi_file,ii,jj,psi,k_rec,time_dim)
       call write_kinematic_field(u_file,ii,jj,u,k_rec,time_dim)
       call write_kinematic_field(v_file,ii,jj,v,k_rec,time_dim)


       k_rec = k_rec + 1
       else
       k_o = k_o + 1
C ------ CALCULATE VELOCITY FIELDS AT THE GRID POINTS BY ANALYTICALLY DIFFERENTIATING ------

       do n = 1,npoints
       !print*,'n=',n
      
C ------ ADVECT PARTICLES ANALYTICALLY -----

        call rk4_kinematic(ii,jj,x0(n),y0(n),pi
     &   ,time_secs,a,factor,cnd,kx,y_centre,dt,dt_nondim
     &        ,x_diff(n),y_diff(n)) ! t is not correct input, need input in days
        
        x0(n) = x0(n) + x_diff(n)
        y0(n) = y0(n) + y_diff(n)
        
        x_traj(n) = x_traj(n) + x_diff(n)
        y_traj(n) = y_traj(n) + y_diff(n)
        !print*,'x_traj(n) = ',n,x_traj(n)
        !print*,'y_traj(n) = ',n,y_traj(n)
        !print*,'x_diff,y_diff = ',x_diff(n),y_diff(n)
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

      enddo
      

c ------- WRITE TO FILE ------------------

      if (k_o .eq. t_tot_day*save_time) then
        k_o = 0.
        call write_kinematic_traj(traj_file,npoints,x_traj,y_traj
     &   ,time_dim,
     & k_rec)
        call write_kinematic_traj(diff_file,npoints,x_diff,y_diff
     &   ,time_dim,
     & k_rec)
      call write_kinematic_field(psi_file,ii,jj,psi,k_rec,time_dim)
       call write_kinematic_field(u_file,ii,jj,u,k_rec,time_dim)
       call write_kinematic_field(v_file,ii,jj,v,k_rec,time_dim)
        k_rec = k_rec + 1  
        print*,'file saved at', time_dim,'days' 
      endif 
      endif
      enddo

      end program kinematic_model
