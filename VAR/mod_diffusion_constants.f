C MODULE THAT CONTAINS ALL THE VARIABLE DECLARATIONS FOR 
C DIFFUSION MODEL.

      module mod_diff_constants
      
      use mod_diff_input
      
      implicit none
      
C ITERATION VARIABLES

      integer t,i,j,n,iseed,b,p,bb
      
C TIME VARIABLES

      integer t_len,step_time
      integer nrec
      real*8 new_time,dt_nondim,k_s
      real*8,allocatable, dimension(:) :: time
      
C FILE NAMES

      character(len=8) :: fmt
      character(len=8) :: x_name
      
      character*100 file_loc
      
c TRAJECTORY VARIABLES

      real*8 x1(npoints,nbins),y1(npoints,nbins),x2(npoints,nbins),
     & y2(npoints,nbins),x1_traj(npoints,nbins), x2_traj(npoints,nbins),
     & y1_traj(npoints,nbins),y2_traj(npoints,nbins),
     & x1_diff,x2_diff,y1_diff,y2_diff, time_av,y_c(jj)
      real*8,allocatable,dimension(:,:,:) :: K
      real*8,allocatable,dimension(:) :: KPV
      real*8 Kx1,Ky1,Kx2,Ky2,dKdx1,dKdy1,dKdx2,dKdy2
      
      real*8 bin_corners(nbins+1),bin_centres(nbins)
     & ,bin_boundaries(nbins+1),bin_width(nbins),uniform_bins(nbins)
    
      
c CABARET VARIABLES

      real*8 psi1_av(ii,jj), psi2_av(ii,jj)
      real*8 uscale,scale,tscale,SS,S1,S2,beta_nondim
     & ,beta_nondim_u1,beta_nondim_u2
     
C INTERPOLATION VARIABLES 
      
      real*8 a1(ii,jj), b1(ii,jj), c1(ii,jj), d1(ii,jj)
     & ,a2(ii,jj),b2(ii,jj),c2(ii,jj),d2(ii,jj)
     & ,psi1_x(4), psi2_x(4)
     & ,u1,v1,u2,v2
      
      
      end module mod_diff_constants
