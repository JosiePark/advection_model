c CONTAINS INPUT PARAMETERS FOR THE DYNAMICAL MODEL 
C TO BE USED WITH PV_MAPPING

      module mod_PVmap_input
          
      implicit none
      
      integer ii,jj,regime,coord_min,coord_max,coord_range,jj1,sim
      real*8 basinscale,H1,H2,Rd,visc,visc_bot,U_0,beta
      
      
      parameter(ii=512,jj=512)                                                !! Also change ii, jj in solv_ell_mike

      parameter(basinscale=520.D5
     & ,H1=1.D5,H2=3.D5,Rd=25.D5
     & ,visc=1.D4
     & ,U_0=6.D0
     & ,regime=1
     & ,sim = 3
     & ,beta = 2.D-13
     & ,coord_min = -3
     & ,coord_max = 3
     & ,coord_range = coord_max - coord_min
     & ,jj1 = int((coord_range+1)*jj))
     
      character*(*),parameter :: home_dir = 
     & '/home/clustor2/ma/j/jp1115/DATA/2/'

     
       character*(*), parameter :: qg_file = 
     &   trim(home_dir) // 'QG/QG.nc'
       character*(*), parameter :: ave_file =
     &   trim(home_dir) // 'QG/QG_ave.nc'
       character*(*), parameter :: traj_file =
     &  trim(home_dir) // 'TRAJ/PV_BINS/
     &test_full_pv_bins_trajectories.nc'
       character*(*), parameter :: disp_file =
     &  trim(home_dir) // 'TRAJ/PV_BINS/
     &test_full_PV_mapped_dispersion.nc'
       character*(*), parameter :: bin_file = 
     & trim(home_dir) // 'TRAJ/PV_BINS/
     &bin_width.nc'  
       

      end module mod_PVmap_input
