C MODULE THAT CONTAINS THE INPUTS FOR THE DIFFUSION MODEL

      module mod_diff_input
      
      implicit none
      
      integer npoints_sqrt,npoints,k_save,sim,nbins,ii,jj,regime,jj1
     & ,coord_range
      real*8 dt,max_run,basinscale,U_0,Rd,H1,H2,beta
      
      parameter(npoints_sqrt = 1,npoints=4000)
      parameter(ii=512,jj=512,jj1= 3*jj)
      parameter(dt = 3600.) ! how often in seconds you wish to advect particle
      parameter(max_run = 2000.) ! length of experiment in days
      parameter(k_save = 1) ! how often you wish to save trajectory data in days
      parameter(sim = 2)
      parameter(nbins = 10)
      parameter(regime = 1)
      parameter(basinscale = 520.D5)
      parameter(U_0 = 6.D0)
    
      parameter(H1=1.D5,H2=3.D5
     & ,beta = 2.D-13
     & ,Rd= 25.D5
     & ,coord_range =3)
      
        character*(*),parameter :: home_dir = 
     & '/home/clustor2/ma/j/jp1115/DATA/2/'

     
       character*(*), parameter :: file_name = 
     &   trim(home_dir) // 'TRAJ/DIFFUSION/
     &full_pv_diffusion_trajectories.nc'
       character*(*), parameter :: ave_file =
     &   trim(home_dir) // 'QG/QG_ave.nc'
       character*(*), parameter :: diff_file =
     &   trim(home_dir) // 'STATS/DIFFUSIVITY/full_MEAN_DIFF.nc'
       character*(*), parameter :: pvdiff_file =
     &   trim(home_dir) // 'STATS/DIFFUSIVITY/full_PVDISP_MEAN_DIFF.nc' 
       character*(*),parameter :: traj_file =
     & trim(home_dir) // 'TRAJ/UNIFORM_BINS/
     &full_uniform_bins_trajectories.nc'
       character*(*),parameter :: bin_file =
     & trim(home_dir) // 'TRAJ/PV_BINS/bin_width.nc'
       character*(*),parameter :: test_file =
     & trim(home_dir) // 'diffusivity_test.nc'
      
      
      end module mod_diff_input
