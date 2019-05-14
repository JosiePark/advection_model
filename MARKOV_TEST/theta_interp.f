      program theta_interp
      
      use mod_theta_netcdf
      use mod_diffusivity_functions
      
      implicit none
      
      character*(*), parameter :: home_dir = 
     & '/home/clustor2/ma/j/jp1115/DATA/1/'
      
      character*(*), parameter :: theta_file = 
     &   trim(home_dir) // 'STATS/THETA/theta-1.nc'
      character*(*), parameter :: interp_file =
     &   trim(home_dir) // 'STATS/THETA/theta_interp.dat'
     
      integer ii,jj,nbins
      parameter(ii = 512, jj = 512, nbins = 10)
      real*8 theta(2,2,nbins) ! theta as read from the file
      real*8 theta_c(2,2,ii) ! interpolated theta
      real*8 y_c(jj)
      real*8 bin_corners(nbins+1),bin_centres(nbins)
      
      integer i,j,m,n,B
      
c DEFINE Y VALUES AT GRID POINTS

      do j = 1,jj
        y_c(j) = dfloat(j-1)
      enddo
      
c CALCULATE BIN CENTRES

      do b = 1,nbins+1
        bin_corners(b) = (b-1)*(dfloat(jj)/nbins)
      enddo
      
      do b = 1,nbins
        bin_centres(b) = .5*(bin_corners(b)+bin_corners(b+1))
      enddo
      
c READ THETA

      call read_theta_netcdf(theta_file,nbins,theta)
      
c INTERPOLATE THETA
     
      do m = 1,2
      do n = 1,2
      do j = 1,jj
        call diffusivity_interp_1d(bin_centres,nbins
     & ,theta(m,n,:),jj,y_c(j),theta_c(m,n,j))
      enddo
      enddo
      enddo
      
c WRITE INTERPOLATED THETA TO FILE 

      open(1,file = interp_file,status = 'new')
      
      do j = 1,jj
        write(1,*) theta_c(1,1,j)
      enddo
      
      close(1)
      
      end program theta_interp
