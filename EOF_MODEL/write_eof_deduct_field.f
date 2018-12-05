c CODE THE READS EOF AND PC FILES AND WRITES A COMBINED EOF FIELD TO A NETCDF FILE

      program write_eof_field
      
      use mod_qg2_netcdf
      use mod_eof_netcdf

      implicit none
      
c NAME FILES

      character*(*),parameter :: home_dir = 
     & '/home/clustor2/ma/j/jp1115/DATA/2/'

       character*(*), parameter :: qg_file = 
     &   trim(home_dir) // 'QG/QG.nc'
       character*(*), parameter :: ave_file =
     &   trim(home_dir) // 'QG/QG_ave.nc'
      character*(*), parameter :: eof_file =
     &  trim(home_dir) // 'STATS/EOF/eof.nc'
      character*(*), parameter :: pc_file =
     &  trim(home_dir) // 'STATS/EOF/pc.nc'
      character*(*), parameter :: new_file = 
     &   trim(home_dir) // 'QG/psi_minus_235678.nc'
      
c TIME VARIABLES
      integer t_len
      real*8, allocatable, dimension(:) :: time
      real*8 time_av,new_tim
      integer eof_option,step_tim
      
c EOF VARIABLES

      integer nmodes,ii,jj
      parameter(nmodes = 6,ii=512,jj=512)
      parameter(eof_option = 0)

      real*8 eof(nmodes,ii,jj,2),PC(nmodes),eof_field(ii,jj,2)
      real*8 psi1_av(ii,jj),psi2_av(ii,jj),psi1(ii,jj),psi2(ii,jj)
      integer modes(nmodes)
      
C LOOPING VARIABLES
      integer i,j,l,m,t
      

      
C READ TIME
      call read_time(qg_file,time,t_len)
      print*,'t_len=',t_len
      
      modes = (/2,3,5,6,7,8/)
      
C READ EOF.NC  

c      call read_eof_netcdf(eof_file,ii,jj,nmodes,2,eof)

      do m = 1,nmodes
           call read_one_eof_netcdf(eof_file,ii,jj,modes(m),
     &           2,eof(m,:,:,:))
      enddo
      
c READ TIME-AVERAGED STREAM FUNCTION

      call read_ave_file(ave_file,ii,jj,psi1_av,psi2_av,time_av)
    
c CREATE EOF_FIELD FILE

      call create_eof_netcdf(new_file,ii,jj,2,t_len)
      
c LOOP THROUGH EACH TIME STEP  

      do t = 1,t_len
      
c INITILIASE EOF_FIELD

      do i = 1,ii
      do j = 1,jj
c        eof_field(i,j,1) = psi1_av(i,j)
c        eof_field(i,j,2) = psi2_av(i,j)
        do l = 1,2
            eof_field(i,j,l) = 0.
        enddo
      enddo
      enddo
      
        print*,'time step = ',t

C READ PC

c      call read_pc_netcdf(pc_file,t,nmodes,PC)

      do m = 1,nmodes
           call read_pc_netcdf(pc_file,t,modes(m),PC(m))
      enddo

C ADD TO EOFS

      do m = 1,nmodes
        eof_field = eof_field + PC(m)*eof(m,:,:,:)
      enddo
      
      if (eof_option .eq. 0) then ! deduct eof_field from the instantaneous velocity field
      
      call read_netcdf(qg_file,psi1,psi2,ii,jj,time(t),
     +  new_tim, step_tim)
     
      eof_field(:,:,1) = psi1 - eof_field(:,:,1)
      eof_field(:,:,2) = psi2 - eof_field(:,:,2)
      
      endif
      

      print*,'eof_field created'


C WRITE TO FILE

      call write_eof_netcdf(new_file,ii,jj,2,t,eof_field)

C END LOOP   

      enddo
      
      end program write_eof_field
