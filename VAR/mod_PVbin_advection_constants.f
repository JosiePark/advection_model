      module mod_PVbin_advection_constants

      use mod_advection_input
      
      implicit none
      
      integer iseed,n,nrec,i,j,k,m,t,k_index,p,bin_index,kk
     & ,idx
      
      real*8 uscale,tscale,scale,PV_x(4),y_index(4),y_c(jj1)
      
c Lagrangian particle positions      
      
      real*8 x0,y0,x1(npoints),y1(npoints),x2(npoints)
     & ,y2(npoints),x1_eddy(npoints),y1_eddy(npoints),x2_eddy(npoints)
     & ,y2_eddy(npoints)
     & ,x1_pseudo(npoints),y1_pseudo(npoints),x2_pseudo(npoints)
     & ,y2_pseudo(npoints),x1_mean(npoints),y1_mean(npoints)
     & ,x2_mean(npoints),y2_mean(npoints)
     
    
c coordinate of Lagrangian particle position
     
      integer x1_coord(npoints),y1_coord(npoints),x2_coord(npoints)
     & ,y2_coord(npoints),x1_eddy_coord(npoints),y1_eddy_coord(npoints)
     & ,x2_eddy_coord(npoints),y2_eddy_coord(npoints)
     & ,x1_pseudo_coord(npoints),y1_pseudo_coord(npoints)
     & ,x2_pseudo_coord(npoints),y2_pseudo_coord(npoints)
     & ,x1_mean_coord(npoints),x2_mean_coord(npoints)
     & ,y1_mean_coord(npoints),y2_mean_coord(npoints)
     
c Define new particle positions for the varying time release, so the particle
c position are defined in a matrix sorted by number of particles and number of 
c releases
c npoints will now be the number of particles per release

      real*8, dimension(:,:), allocatable :: x1r,y1r,x2r
     & ,y2r,x1r_eddy,y1r_eddy,x2r_eddy
     & ,y2r_eddy,x1r_pseudo,y1r_pseudo,x2r_pseudo
     & ,y2r_pseudo
     & ,x1r_mean,x2r_mean,y1r_mean,y2r_mean
     
     
      integer, dimension(:,:), allocatable :: x1r_coord,y1r_coord
     & ,x2r_coord
     & ,y2r_coord,x1r_eddy_coord,y1r_eddy_coord,x2r_eddy_coord
     & ,y2r_eddy_coord,x1r_pseudo_coord,y1r_pseudo_coord
     & ,x2r_pseudo_coord
     & ,y2r_pseudo_coord,x1r_mean_coord,y1r_mean_coord
     & ,x2r_mean_coord
     & ,y2r_mean_coord

     
      real*8, dimension(:), allocatable :: release_time ! contains the release time for each realisation
      integer, dimension(:), allocatable :: nrel ! records the number of times the trajectory has been written 
                                                 ! to the netcdf file for each realisation,
                                                 ! this is so we know where to write the trajectory
                                                 
      real*8,dimension(:), allocatable :: time,time_o ! time_o is offline time array
      real*8,dimension(:,:), allocatable :: grad
      real*8,dimension(:),allocatable :: time_dim
      
c file names

      character*100 eddy_file,pseudo_file
     & ,full_file
      character*100 str
                                                       
      
c Names of time variables

      integer t_len,k_old,time_interp,t_tot,step_tim,k_new,k_half
      integer max_time_lag,t_start,release_start,start_info(2)
      real*8 time_cubic(4), time_half, dt05, dt6, time_step
     & ,new_tim,time_day,dt_nondim,release_length_secs,time_av
      integer,dimension(:),allocatable :: k_s
      integer l_tot_day
      integer k_t
c Names of model variables

      real*8 psi1(ii,jj,4),psi2(ii,jj,4)
     & ,psi1_half(ii,jj),psi2_half(ii,jj)
     & ,psi1_av(ii,jj),psi2_av(ii,jj),psi1_old(ii,jj)
     & ,psi2_old(ii,jj), psi1_new(ii,jj), psi2_new(ii,jj)
     & ,psi1_eddy_old(ii,jj), psi2_eddy_old(ii,jj)
     & ,psi1_eddy_half(ii,jj),psi2_eddy_half(ii,jj)
     & ,psi1_eddy_new(ii,jj), psi2_eddy_new(ii,jj)
     & ,rel1(ii,jj),rel2(ii,jj),zeta1(ii,jj),zeta2(ii,jj)
     & ,S1,S2,beta1_y(jj),beta_nondim_u1,SS,beta_nondim
     & ,beta_nondim_u2
     & ,psi1_eof(ii,jj,4),psi2_eof(ii,jj,4)
     & ,psi1_eof_half(ii,jj), psi2_eof_half(ii,jj)
     & ,psi1_eof_old(ii,jj), psi2_eof_old(ii,jj)
     & ,psi1_eof_new(ii,jj), psi2_eof_new(ii,jj)
      
c Variables needed for spatial interpolation (i.e polynomial coefficients)

      real*8 a1_old(ii,jj),b1_old(ii,jj),c1_old(ii,jj),d1_old(ii,jj)
     & ,a2_old(ii,jj),b2_old(ii,jj),c2_old(ii,jj),d2_old(ii,jj)
     & ,a1_half(ii,jj),b1_half(ii,jj), c1_half(ii,jj),d1_half(ii,jj)
     & ,a2_half(ii,jj),b2_half(ii,jj),c2_half(ii,jj),d2_half(ii,jj)
     & ,a1_new(ii,jj),b1_new(ii,jj),c1_new(ii,jj),d1_new(ii,jj)
     & ,a2_new(ii,jj),b2_new(ii,jj),c2_new(ii,jj),d2_new(ii,jj)
     
     & ,a1_eddy_old(ii,jj),b1_eddy_old(ii,jj)
     & ,c1_eddy_old(ii,jj),d1_eddy_old(ii,jj)
     & ,a2_eddy_old(ii,jj),b2_eddy_old(ii,jj)
     & ,c2_eddy_old(ii,jj),d2_eddy_old(ii,jj)
     & ,a1_eddy_half(ii,jj),b1_eddy_half(ii,jj)
     & ,c1_eddy_half(ii,jj),d1_eddy_half(ii,jj)
     & ,a2_eddy_half(ii,jj),b2_eddy_half(ii,jj)
     & ,c2_eddy_half(ii,jj),d2_eddy_half(ii,jj)
     & ,a1_eddy_new(ii,jj),b1_eddy_new(ii,jj)
     & ,c1_eddy_new(ii,jj),d1_eddy_new(ii,jj)
     & ,a2_eddy_new(ii,jj),b2_eddy_new(ii,jj)
     & ,c2_eddy_new(ii,jj),d2_eddy_new(ii,jj)
     
     & ,a1_av(ii,jj),b1_av(ii,jj),c1_av(ii,jj),d1_av(ii,jj)
     & ,a2_av(ii,jj),b2_av(ii,jj),c2_av(ii,jj),d2_av(ii,jj)
     
     & ,a(ii,-3:jj+3),b(ii,-3:jj+3),c(ii,-3:jj+3),d(ii,-3:jj+3)
     
     & ,psi1_x(4),psi2_x(4)
     
     & ,M1_old(ii,jj,4,4),M2_old(ii,jj,4,4)
     & ,M1_half(ii,jj,4,4),M2_half(ii,jj,4,4)
     & ,M1_new(ii,jj,4,4),M2_new(ii,jj,4,4)
     
     & ,M1_eddy_old(ii,jj,4,4),M2_eddy_old(ii,jj,4,4)
     & ,M1_eddy_half(ii,jj,4,4),M2_eddy_half(ii,jj,4,4)
     & ,M1_eddy_new(ii,jj,4,4),M2_eddy_new(ii,jj,4,4)
     
     & ,M1_av(ii,jj,4,4),M2_av(ii,jj,4,4)
     
c Variables needed for RK4 advection

      real*8 u1_point0, v1_point0, u2_point0, v2_point0
     & ,u1_point1, v1_point1, u2_point1, v2_point1
     & ,u1_point2, v1_point2, u2_point2, v2_point2
     & ,u1_point3, v1_point3, u2_point3, v2_point3
      
     & ,xt1,yt1,xt2,yt2,x_diff1,y_diff1,x_diff2,y_diff2
     & ,x_av_diff1,y_av_diff1,x_av_diff2,y_av_diff2
     & ,x_disp1(npoints),y_disp1(npoints),x_disp2(npoints)
     & ,y_disp2(npoints)
     
c Bin variables
        
      real*8 PV_bin(nbins),d_bin,Y_bin(nbins),j_bin(jj1)
      real*8 PV_bar(jj1),rel1_av(ii,jj),rel2_av(ii,jj),zeta1_av(ii,jj)
     & ,zeta2_av(ii,jj),PV_0,y_map,PV_bar_domain(jj) ,PV_av(ii,jj1)
     & ,pv(ii,jj)
      integer bin_count(nbins),points_count,bin_stop(nbins),index
      
      end module mod_PVbin_advection_constants
