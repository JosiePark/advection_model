      module mod_theta_netcdf
      
      use mod_netcdf_error
      
      contains
      
      subroutine read_theta_netcdf(file_name,nbins,theta)
      
      implicit none
      include 'netcdf.inc'
      
      character*(*) file_name
      integer nbins
      real*8 theta(2,2,nbins)
      
      character*(*) theta_name,bin_name,layer_name,dim_name
      parameter(theta_name = 'Theta')
      parameter(bin_name = 'Bin')
      parameter(layer_name = 'Layer')
      parameter(dim_name = 'Dimension')
      
      integer theta_varid,bin_dimid,layer_dimid,dim_dimid
      integer ndims
      parameter(ndims = 3)
      integer dimid(ndims),start(ndims),count(ndims)
      
      integer ncid, retval
      
      retval = nf_open(file_name,nf_nowrite,ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_inq_varid(ncid,theta_name,theta_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      start(1) = 1
      start(2) = 1
      start(3) = 1
      count(1) = 2
      count(2) = 2
      count(3) = nbins
      
      retval = nf_get_vara_double(ncid,theta_varid,start,count,theta)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      
      end subroutine read_theta_netcdf
      
      end module mod_theta_netcdf
