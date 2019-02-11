c subroutine to write velocity variance to file

      module mod_vel_variance_netcdf
      
      use mod_netcdf_error
      
      contains
      
      subroutine create_vel_variance(file_name,ii,jj)
      
      implicit none
      include 'netcdf.inc'
      
      character*(*) file_name
      integer ii,jj
      
      integer ndims,nlayers
      parameter(ndims = 5, nlayers = 2)
      integer start(ndims),count(ndims),dimid(ndims)
      
      ! dimensions
      integer x_dimid,y_dimid,lvl_dimid,dim1_dimid,dim2_dimid
      character*(*) x_name,y_name,lvl_name,dim1_name,dim2_name
      parameter(x_name = 'X', y_name = 'Y', lvl_name = 'Layer',
     & dim1_name = 'Dimension 1',dim2_name = 'Dimension 2')
     
      ! variables
      integer var_varid
      character*(*) var_name
      parameter(var_name = 'Velocity Variance')
      
      integer ncid,retval
      
      retval = nf_create(file_name,nf_clobber,ncid)
      if (retval .ne. nf_noerr) call handle_error(retval)
      
      ! define dimensions
      retval = nf_def_dim(ncid,lvl_name,nlayers,lvl_dimid)
      if (retval .ne. nf_noerr) call handle_error(retval)
      retval = nf_def_dim(ncid,x_name,ii,x_dimid)
      if (retval .ne. nf_noerr) call handle_error(retval)
      retval = nf_def_dim(ncid,y_name,jj,y_dimid)
      if (retval .ne. nf_noerr) call handle_error(retval)
      retval = nf_def_dim(ncid,dim1_name,2,dim1_dimid)
      if (retval .ne. nf_noerr) call handle_error(retval)
      retval = nf_def_dim(ncid,dim2_name,2,dim2_dimid)
      if (retval .ne. nf_noerr) call handle_error(retval)
      
      ! define variables
      
      dimid(1) = dim1_dimid
      dimid(2) = dim2_dimid
      dimid(3) = x_dimid
      dimid(4) = y_dimid
      dimid(5) = lvl_dimid
      retval = nf_def_var(ncid,var_name,nf_double,ndims,dimid,var_varid)
      
      retval = nf_enddef(ncid)                                        !! End define mode    
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_close(ncid)                                         !! Close file
      if (retval .ne. nf_noerr) call handle_err(retval)
   
      print *,'!!! Netcdf file created i.e. ', FILE_NAME, '!!!'
      print *,'!!!                                         !!!'
      
      
      end subroutine create_vel_variance
      
c ------------------------------------------------------------------
      
      subroutine write_vel_variance(file_name,sigma1,sigma2,ii,jj)
      
      implicit none
      include 'netcdf.inc'
    
      character*(*) file_name
      integer ii,jj
      
      integer ndims,nlayers
      parameter(ndims = 5, nlayers = 2)
      integer start(ndims),count(ndims),dimid(ndims)
      
      ! dimensions
      integer x_dimid,y_dimid,lvl_dimid,dim_dimid
      character*(*) x_name,y_name,lvl_name,dim1_name,dim2_name
      parameter(x_name = 'X', y_name = 'Y', lvl_name = 'Layer',
     & dim1_name = 'Dimension 1',dim2_name = 'Dimension 2')
     
      ! variables
      integer var_varid
      character*(*) var_name
      parameter(var_name = 'Velocity Variance')
      
      integer ncid,retval
      
      real*8 sigma1(2,2,ii,jj),sigma2(2,2,ii,jj)
      
      retval = nf_open(FILE_NAME, nf_write, ncid)                     !! Open file in write mode
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_inq_varid(ncid,var_name,var_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
      start(5) = 1
      
      count(1) = 2
      count(2) = 2
      count(3) = ii
      count(4) = jj
      count(5) = 1
      
      retval = nf_put_vara_double(ncid, var_varid, start, count,      !! Update Data
     +        sigma1)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
      start(5) = 2
      
      retval = nf_put_vara_double(ncid, var_varid, start, count,      !! Update Data
     +        sigma2)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = nf_close(ncid)                                         !! Close NetCDF file                
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        print *,'!!! Netcdf file modified i.e. ', FILE_NAME, '!!!'
      
      
      end subroutine write_vel_variance
    
c -------------------------------------------------------------------
      
c      subroutine read_vel_variance
      
c      implicit none
      
c      end subroutine read_vel_variance
      
      end module MOD_vel_variance_netcdf
