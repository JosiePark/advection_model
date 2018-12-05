C CODE THAT READS A SNAPSHOT OF THE STREAMFUNCTION AND WRITES TO A NEW NETCDF FILE
C CREATES A RESTART FILE SO THAT CABARET CAN BE RUN AFRESH

      program write_restart_file
      
      use MOD_qg2_netcdf
      use MOD_qg_input
      use MOD_constants
      use MOD_variables
      
      implicit none
      
      integer t_tot,t_len
      real*8,allocatable,dimension(:) :: time
      
      
      character (LEN=*), parameter:: OLD_NAME=
     & '/rds/general/user/jp1115/home/WORK/DATA/1/QG/QG.nc' !! Name the old NetCDF file
      character (LEN=*), parameter:: NEW_NAME= 'QG_condensed.nc'

      t_tot = 100
      
      call create_netcdf_file(new_name,ii,jj,basinscale)
      call read_time(old_name,time,t_len)
      
      do i = 1,t_tot
      
      print*,'time=',time(i)
      
      call read_netcdf(old_name,psi1,psi2,ii,jj,time(i)+0.001
     & ,new_tim,step_tim)
     
      print*,'new_tim = ',new_tim
      print*,'t=',i
      print*,'step_tim = ',step_tim
      
      call energy(H1,H2,S1,S2,ii,jj,psi1,psi2,ekin1,ekin2,epot)
     
      call write_netcdf(new_name,psi1,psi2,
     +    epot,ekin1,ekin2,ii,jj,new_tim,i)
     
      enddo
     
      end program write_restart_file
