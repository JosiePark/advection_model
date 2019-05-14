      program kinematic_model
      
      use mod_random
      use mod_kinematic_netcdf
      use mod_kinematic_advection
      
      implicit none
      
      real*8 factor,pi,A,T,c,lambda,max_time,dt,max_time_secs,t_tot,kx
      real*8 y_centre,time_day
      real*8 basinscale,scale,uscale,tscale,dt_nondim
      real*8 time_o,time_dim,time_secs
      real*8 bin_width,u0
      integer layer,regime
      integer ii,jj,k_rec,nbins
      parameter(ii = 512, jj = 512)
      parameter(max_time = 1000)
      !parameter(dt = 2160.)
      parameter(dt = 2160.)
      parameter(basinscale = 520.d5)
      parameter(nbins = 10)
      parameter(u0 = 0.d0)
      parameter(layer = 1)
      parameter(regime = 2)
      
      integer i,j,n,t_tot_day,iseed,b,k
      
      integer npoints,npoints_sqrt
      parameter(npoints_sqrt = 50, npoints = npoints_sqrt**2)
      
      real*8 save_time,eps
      parameter(save_time = 1.) ! save every save_time days
      
      real*8 x0(nbins,npoints),y0(nbins,npoints)
     & ,x_traj(nbins,npoints),y_traj(nbins,npoints)
      real*8 x_diff, y_diff
      real*8 u(ii,jj), v(ii,jj)
      
      character*(*),parameter :: home_dir = 
     & '/home/clustor2/ma/j/jp1115/DATA/2/'
      
      character*(*), parameter :: psi_file = 
     &   trim(home_dir) // 'TRAJ/KINEMATIC/psi_bottom.nc'
      character*(*), parameter :: u_file = 
     &   trim(home_dir) // 'TRAJ/KINEMATIC/u_bottom.nc'
      character*(*), parameter :: v_file = 
     &   trim(home_dir) // 'TRAJ/KINEMATIC/v_bottom.nc'
      character*(*), parameter :: traj_file = 
     &   trim(home_dir) // 'TRAJ/KINEMATIC/traj_top.nc'
      character*(*), parameter :: diff_file = 
     &   trim(home_dir) // 'TRAJ/KINEMATIC/diff_bottom.nc'
      
      real*8 psi(ii,jj),x_c(ii),y_c(jj),xnd(ii),ynd(jj),cnd
      real*8 k_o
      data pi/3.14159265358979323846D0/
      
      print*,'Running kinematic model'
      
      scale = basinscale/dfloat(ii)
      uscale = 1.
      tscale = scale/uscale
      
    
C ------- KINEMATIC MODEL PARAMETERS --------------

      if (regime .eq. 1) then

          if (layer .eq. 1) then

            factor = 3.*pi ! width scaling
            A = 526.1 ! amplitude of wave
            
          else
          
            factor = 6.54
            A = 146.27
        
          endif
          
      elseif (regime .eq. 2) then
      
          if (layer .eq. 1) then

            factor = 7.393939 ! width scaling
            A = 609.8531 ! amplitude of wave
            
          else
          
            factor = 5.47474747
            A = 191.977480738
        
          endif
      
      endif
      
      kx = 2 ! wave number
      if (regime .eq. 1) then
       T = 48. ! period in days
      else
       T = 44.
      endif
      lambda = float(ii)/(kx) ! wavelength
      c = lambda/T ! propagation speed
      !c = 1.
      if (regime .eq. 1) then
       y_centre = 302. ! centre in grid points
      else
       y_centre = 282.
      endif
      cnd = c*2*pi/ii
      bin_width = dfloat(jj)/nbins
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

      
c -------- CREATE FILES ------------------

      call create_kinematic_field(psi_file,ii,jj)
      call create_kinematic_field(u_file,ii,jj)
      call create_kinematic_field(v_file,ii,jj)
      call create_kinematic_traj(traj_file,npoints,nbins)
      print*,'created files'
      
        k_o = 0.
            
C ------- MAIN CYCLE - ITERATE IN TIME ------------

      do k = 1,t_tot
        time_o = (k-1)*dt_nondim
        time_dim = time_o*tscale/86400.
        time_secs = (k-1)*dt
        time_day = time_secs/86400.
        
C ------ CONSTRUCT KINEMATIC FIELD --------
        do i = 1,ii
        do j = 1,jj
            psi(i,j) = A/(cosh(factor*(y_c(j)-y_centre)/ii))**2
     &       *cos(kx*(xnd(i)-cnd*time_secs/86400.))
     
            call kinematic_velocity(ii,jj,dfloat(i-1),dfloat(j-1)
     &  ,pi,time_secs,a,factor,cnd,kx,y_centre,u(i,j),v(i,j))
        enddo
        enddo
        
        if (k .eq. 1) then
        k_rec = 1

C ------ RANDOMLY GENERATE PARTICLES -----------

        iseed = 102

        do b = 1,nbins
        do n = 1,npoints
            x0(b,n) = ran1(iseed)*dfloat(jj)
            y0(b,n) = ran1(iseed)*bin_width + (b-1)*bin_width
            !y0(n) = 310.
            x_traj(b,n) = x0(b,n)
            y_traj(b,n) = y0(b,n)
        enddo
        enddo
       call write_kinematic_traj(traj_file,npoints,nbins
     & ,x0,y0,dfloat(0),k_rec)
       call write_kinematic_field(psi_file,ii,jj,psi,k_rec,time_dim)
       call write_kinematic_field(u_file,ii,jj,u,k_rec,time_dim)
       call write_kinematic_field(v_file,ii,jj,v,k_rec,time_dim)


       k_rec = k_rec + 1
       else
       k_o = k_o + 1
C ------ CALCULATE VELOCITY FIELDS AT THE GRID POINTS BY ANALYTICALLY DIFFERENTIATING ------

       do b = 1,nbins
       do n = 1,npoints
       !print*,'n=',n
      
C ------ ADVECT PARTICLES ANALYTICALLY -----

        call rk4_kinematic(ii,jj,x0(b,n),y0(b,n),pi
     &   ,time_secs,u0,a,factor,cnd,kx,y_centre,dt,dt_nondim
     &        ,x_diff,y_diff) 
        
        x0(b,n) = x0(b,n) + x_diff
        y0(b,n) = y0(b,n) + y_diff
        
        x_traj(b,n) = x_traj(b,n) + x_diff
        y_traj(b,n) = y_traj(b,n) + y_diff

c ------ TREAT BOUNDARIES --------------

        if (x0(b,n) .gt. ii) then
            x0(b,n) = x0(b,n) - ii
        else if (x0(b,n) .lt. 0) then
            x0(b,n) = x0(b,n) + ii
        endif
        
        if (y0(b,n) .gt. jj) then
            y0(b,n) = y0(b,n) - jj
        else if (y0(b,n) .lt. 0) then
            y0(b,n) = y0(b,n) + jj
        endif

      enddo
      enddo
      

c ------- WRITE TO FILE ------------------

      if (k_o .eq. t_tot_day*save_time) then
        k_o = 0.
        call write_kinematic_traj(traj_file,npoints,nbins,x_traj,y_traj
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
