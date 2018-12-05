c module that contains variable declaration for PV_MAPPING

      module mod_PVmap_constants
      
      use mod_PVmap_input
      
      implicit none 
      
c iterative variables and indices
      integer i,j,kk,idx,t,k,n,exact,m,k_index,j_index,cubic_index(4)
     & ,i_index,p
      
c dynamical model parameters
      real*8 uscale,scale,tscale,SS,beta_nondim,beta_nondim_u1
     & ,beta_nondim_u2,S1,S2
     
      real*8, dimension(:), allocatable :: x_c,y_c 
      
     
c time data
      integer qg_t_len,traj_t_len
      real*8,dimension(:),allocatable :: qg_time
      real*8 time_av
      real*8,dimension(:,:),allocatable :: traj_time
     
c dynamical model variables
      real*8 psi1_av(ii,jj), psi2_av(ii,jj)
     & ,rel1_av(ii,jj),rel2_av(ii,jj)
     & ,zeta1_av(ii,jj), zeta2_av(ii,jj)
     & ,PV_av(ii,jj1),beta1_y(jj1),PV_bar(jj1)
     & ,psi1(ii,jj,4),psi2(ii,jj,4), rel1(ii,jj)
     & ,rel2(ii,jj),zeta1(ii,jj),zeta2(ii,jj)
     & ,PV(ii,jj1)
     
c transport model variables
      integer nrel,npoints,nbins
      
      real*8,dimension(:),allocatable :: x1,y1,y_loc
      integer,dimension(:),allocatable :: x1_coord,y1_coord
      real*8,dimension(:,:),allocatable :: bin_width,bin_boundaries
      
c interpolation variables
      real*8 time_cubic(4),psi1_interp(ii,jj),psi2_interp(ii,jj)
      real*8 a(ii,jj1)
     &  ,b(ii,jj1),c(ii,jj1)
     &  ,d(ii,jj1),PV_x(4),y_index(4)
       
      real*8,dimension(:),allocatable :: PV_0,PV_t
c mapped variable
      real*8,dimension(:),allocatable :: y0_map
      real*8,dimension(:),allocatable :: y_map
      
      
      end module mod_PVmap_constants
