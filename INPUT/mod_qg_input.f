c MODULE THAT CONTAINS INPUT FOR DYNAMICAL MODEL

      module mod_qg_input
      
      implicit none
      
      integer ii,jj,ips,max_spin,max_time,istart
     & ,regime,istart_ave,ii1,jj1,ii2,jj2
      real*8 basinscale,H1,H2,h_max,Rd,visc,visc_bot,U_0
     & ,read_tim,time_out,time_save,cfl 
      
      
      parameter(ips=1,max_spin=0,max_time=20000+max_spin  
     & ,ii=512,jj=512)                                                   !! Also change ii, jj in solv_ell_mike

      parameter(basinscale=520.D5
     & ,H1=1.D5,H2=3.D5,h_max=1.D2,Rd=25.D5
     & ,visc=1.D4
     & ,U_0=6.D0
     & ,istart=1
     & ,istart_ave = 1
     & ,regime=2
     & ,read_tim=0.+0.0001
     & ,TIME_OUT=10.D0     ! accumulate data every TIME_OUT days
     & ,TIME_SAVE=1.D0
     & ,CFL=0.4 )

      parameter(ii1=ii-1,jj1=jj-1,ii2=ii-2,jj2=jj-2)
     
      character*(*),parameter :: home_dir = 
     & '/work/jp1115/saves/'
     
      character*(*),parameter :: file_name = 'QG_new.nc'
      character*(*),parameter :: ave_file = 'QG_ave_new.nc'


      end module mod_qg_input
