c Code that creates time averaged qg file from different files

      program create_ave
      
      use mod_qg2_netcdf
      use mod_time_interp
      
      implicit none
      
      integer ii,jj
      parameter(ii=512,jj=512)
      real*8 psi1(ii,jj), psi2(ii,jj)
     & ,psi1_temp(ii,jj,4), psi2_temp(ii,jj,4)
     & ,psi1_av(ii,jj),psi2_av(ii,jj),basinscale 
     & ,tscale,new_tim
     & ,time_nondim  
     
      real*8, allocatable, dimension(:) :: time
      integer t_len,step_time,t,i,j
      parameter(basinscale = 520.D5) 
     
       character (LEN=*), parameter:: file_name = 
     & "/work/jp1115/saves/1/QG/QG.nc"
        character (LEN=*), parameter:: ave_name = 'QG_ave_new.nc'
        
        
        call read_time(file_name,time,t_len)
        
      do i = 1,ii
      do j = 1,jj
        psi1_av(i,j) = 0.
        psi2_av(i,j) = 0.
      enddo
      enddo

      do t = 1,t_len
      
        print*,'time = ',time(t)
      
        if (t .eq. 1) then
            call read_psi(file_name,ii,jj,1,psi1_temp,psi2_temp)
            psi1 = psi1_temp(:,:,1)
            psi2 = psi2_temp(:,:,1)
        else
        call read_netcdf(FILE_NAME,psi1,psi2,ii,jj,time(t),
     +  new_tim, step_time)
        endif
        
        
        if (t.ne.1) then
        
            psi1_av = psi1_av + (time(t)-time(t-1))*psi1
            psi2_av = psi2_av + (time(t)-time(t-1))*psi2
            
        endif
        
      enddo
      
      psi1_av = psi1_av/time(t_len)
      psi2_av = psi2_av/time(t_len)
      
      
      call create_ave_file(ave_name,ii,jj,basinscale)
      
      call write_ave_file(ave_name,ii,jj,psi1_av,psi2_av,time(t_len))
      
      end program create_ave
