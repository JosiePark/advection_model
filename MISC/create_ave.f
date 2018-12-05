c Code that creates time averaged qg file from different files

      program create_ave
      
      use mod_read_psi_av
      use mod_qg2_netcdf
      
      implicit none
      
      integer ii,jj
      parameter(ii=512,jj=512)
      real*8 psi1_av1(ii,jj), psi2_av1(ii,jj)
     & ,psi1_av2(ii,jj),psi2_av2(ii,jj),basinscale 
     & ,psi1_av(ii,jj),psi2_av(ii,jj),time,tscale
     & ,time_nondim  
      parameter(basinscale = 520.D5) 
      
      character (LEN=*), parameter:: file_name1=
     & "/media/josiepark/Seagate Expansion Drive"//
     & "/PhD/DATA/Saves/3/AVE/QG_ave.nc"    

      character (LEN=*), parameter:: file_name2=
     & "/media/josiepark/Seagate Expansion Drive"//
     & "/PhD/DATA/Saves/3/AVE/QG_ave(1).nc"
     
      character (LEN=*), parameter:: file_name=
     & "/media/josiepark/Seagate Expansion Drive"//
     & "/PhD/DATA/Saves/3/AVE/QG_ave(tot).nc"
     
      print*,file_name1
      call psi_av(file_name1,ii,jj,psi1_av1,psi2_av1)
      
      call psi_av(file_name2,ii,jj,psi1_av2,psi2_av2)
      
      time = 10000.
      tscale = basinscale/dfloat(jj)
      time_nondim = time*86400/tscale
      
      psi1_av = (psi1_av1+psi1_av2)/time_nondim
      psi2_av = (psi2_av1+psi2_av2)/time_nondim
      
      
      call create_ave_file(file_name,ii,jj,basinscale)
      
      call write_ave_file(file_name,ii,jj,psi1_av,psi2_av,time)
      
      end program create_ave
