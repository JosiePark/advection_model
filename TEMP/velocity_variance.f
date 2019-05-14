C CODE THAT CALCULATES VELOCITY VARIANCE TENSOR AND WRITES TO FILE

      program velocity_variance
      
      use mod_vel_variance_netcdf
      use mod_qg2_netcdf
      use mod_variables
      
      implicit none
      
      integer ii,jj,nlayers
      parameter(ii = 512,jj=512,nlayers = 2)
      
      real*8 psi1_av(ii,jj),psi2_av(ii,jj),sigma1(2,2,ii,jj),
     & sigma2(2,2,ii,jj),psi1(ii,jj),psi2(ii,jj),u1(2,ii,jj),u2(2,ii,jj)
     
      integer t_len
      real*8 ,dimension(:), allocatable :: time
      real*8 time_av,new_tim
      integer step_tim
      
      integer t,i,j,k,m
      
      character*(*), parameter :: home_dir =
     & '/home/clustor2/ma/j/jp1115/DATA/1/'
     
      character*(*), parameter :: ave_file =
     &   trim(home_dir) // 'QG/QG_ave.nc'
      character*(*), parameter :: qg_file = 
     &   trim(home_dir) // 'QG/QG.nc'
      character*(*), parameter :: file_name = 
     &   trim(home_dir) // 'STATS/SIGMA/velocity_variance.nc'
     
      print*,'start_running'
      
      
C READ TIME AVERAGED STREAM FUNCTION

      do i = 1,ii
      do j = 1,jj
      do k = 1,2
      do m = 1,2
      sigma1(k,m,i,j) = 0.
      sigma2(k,m,i,j) = 0 .
      enddo
      enddo
      enddo
      enddo

      call read_ave_file(ave_file,ii,jj,psi1_av,psi2_av,time_av)
      
      print*,'psi_ave read'

C READ TIME

       call read_time(qg_file,time,t_len)
       
       print*,'time_read'

C AT EACH TIME STEP, FIND VELOCITY FLUCTUATION 

        do t = 1,t_len
        
            call read_netcdf(qg_file,psi1,psi2,ii,jj,time(t),
     +  new_tim, step_tim)
     
C CALCULATE VELOCITY AT EACH GRID POINT

            call vel_from_psi(ii,jj,psi1-psi1_av,psi2-psi2_av
     &       ,u1(1,:,:),u1(2,:,:),u2(1,:,:),u2(2,:,:))
            
            
C CALCULATE PRODUCT

            do i = 1,ii
            do j = 1,jj
                do k = 1,2
                do m = 1,2
                 sigma1(k,m,i,j) = sigma1(k,m,i,j)+u1(k,i,j)*u1(m,i,j)
                 sigma2(k,m,i,j) = sigma2(k,m,i,j)+u2(k,i,j)*u2(m,i,j)
                enddo
                enddo
            enddo
            enddo

        enddo

C TAKE TIME AVERAGE

        sigma1 = sigma1/(t_len)
        sigma2 = sigma2/(t_len)

C WRITE TO FILE

      call create_vel_variance(file_name,ii,jj)
      call write_vel_variance(file_name,sigma1,sigma2,ii,jj)
      
      end program velocity_variance
