      module mod_pvdiffusion_input
      
      implicit none
      
c FILE NAMES
      
      character*(*),parameter :: home_dir = 
     & '/home/clustor2/ma/j/jp1115/DATA/1/'
     
      character*(*),parameter :: traj_file =
     & trim(home_dir) // 'TRAJ/PV_BINS/full_PV_bins_trajectories.nc'
      character*(*),parameter :: ave_file = 
     & trim(home_dir) // 'QG/QG_ave.nc' 
      character*(*),parameter :: qg_file = 
     & trim(home_dir) // 'QG/QG.nc'
      character*(*),parameter :: diffusivity_file =
     & trim(home_dir) // 'STATS/DIFFUSIVITY/full_PV_MEAN_DIFF.nc'
      character*(*),parameter :: diffusion_file = 
     & trim(home_dir) // 'TRAJ/DIFFUSION/
     &full_PV_diffusion_trajectories.nc'
       character*(*),parameter :: bin_file =
     & trim(home_dir) // 'TRAJ/PV_BINS/bin_width.nc'
     
c PARAMETERS

      integer ii,jj,jj1
      
      real*8 basinscale,H1,H2,beta,U_0,Rd,coord_range
      
      parameter(ii=512,jj=512,jj1 = 3*jj)
      parameter(basinscale=520.D5
     & ,H1=1.D5,H2=3.D5
     & ,beta = 2.D-13
     & ,U_0=6.D0
     & ,Rd= 25.D5
     & ,coord_range =3)
     
      integer k_save,npoints_sqrt,npoints,nbins
      real*8 dt,max_run
     
      parameter(npoints_sqrt = 1,npoints=4000)
      parameter(dt = 3600.) ! how often in seconds you wish to advect particle
      parameter(max_run = 2000.) ! length of experiment in days
      parameter(k_save = 1) ! how often you wish to save trajectory data in days
      parameter(nbins = 10)
     
      end module
