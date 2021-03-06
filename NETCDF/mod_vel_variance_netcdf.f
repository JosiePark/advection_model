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
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      ! define dimensions
      retval = nf_def_dim(ncid,lvl_name,nlayers,lvl_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid,x_name,ii,x_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid,y_name,jj,y_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid,dim1_name,2,dim1_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid,dim2_name,2,dim2_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      ! define variables
      
      dimid(1) = dim1_dimid
      dimid(2) = dim2_dimid
      dimid(3) = x_dimid
      dimid(4) = y_dimid
      dimid(5) = lvl_dimid
      retval = nf_def_var(ncid,var_name,nf_double,ndims,dimid,var_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
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
      
      subroutine read_vel_variance(file_name,sigma1,sigma2,ii,jj)
      
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
      
      retval = nf_open(file_name, nf_nowrite, ncid)                   !! Open file in read mode
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
      
      retval = nf_get_vara_double(ncid,var_varid,start,count,sigma1)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      start(5) = 2
      
      retval = nf_get_vara_double(ncid,var_varid,start,count,sigma2)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      
      
      end subroutine read_vel_variance
      
c ---------------------------------------------------------------------------

      subroutine create_variable_netcdf(file_name,ii,jj)
      
      implicit none
      include 'netcdf.inc'
      
      character*(*) file_name
      integer ii,jj
      
      integer ndims,nlayers
      parameter(ndims = 4, nlayers = 2)
      integer start(ndims),count(ndims),dimid(ndims)
      
      ! dimensions
      integer x_dimid,y_dimid,lvl_dimid,dim_dimid
      character*(*) x_name,y_name,lvl_name,dim_name
      parameter(x_name = 'X', y_name = 'Y', lvl_name = 'Layer',
     & dim_name = 'Dimension')
     
      ! variables
      integer var_varid
      character*(*) var_name
      parameter(var_name = 'Variable')
      
      integer ncid,retval
      
      retval = nf_create(file_name,nf_clobber,ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      ! define dimensions
      retval = nf_def_dim(ncid,lvl_name,nlayers,lvl_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid,x_name,ii,x_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid,y_name,jj,y_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid,dim_name,2,dim_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      ! define variables
      
      dimid(1) = dim_dimid
      dimid(2) = x_dimid
      dimid(3) = y_dimid
      dimid(4) = lvl_dimid
      retval = nf_def_var(ncid,var_name,nf_double,ndims,dimid,var_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_enddef(ncid)                                        !! End define mode    
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_close(ncid)                                         !! Close file
      if (retval .ne. nf_noerr) call handle_err(retval)
   
      print *,'!!! Netcdf file created i.e. ', FILE_NAME, '!!!'
      print *,'!!!        !!!'
      
      end subroutine create_variable_netcdf
      
c -----------------------------------------------------------------------------
      
      subroutine write_variable_netcdf(file_name,var_layer1,var_layer2
     & ,ii,jj)
      
      implicit none
      include 'netcdf.inc'
      
      character*(*) file_name
      integer ii,jj
      
      integer ndims,nlayers
      parameter(ndims = 4, nlayers = 2)
      integer start(ndims),count(ndims),dimid(ndims)
      
      ! dimensions
      integer x_dimid,y_dimid,lvl_dimid,dim_dimid
      character*(*) x_name,y_name,lvl_name,dim_name
      parameter(x_name = 'X', y_name = 'Y', lvl_name = 'Layer',
     & dim_name = 'Dimension')
     
      ! variables
      integer var_varid
      character*(*) var_name
      parameter(var_name = 'Variable')
      
      integer ncid,retval
      
      real*8 var_layer1(2,ii,jj),var_layer2(2,ii,jj)
      
      retval = nf_open(FILE_NAME, nf_write, ncid)                     !! Open file in write mode
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_inq_varid(ncid,var_name,var_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
      
      count(1) = 2
      count(2) = ii
      count(3) = jj
      count(4) = 1
      
      retval = nf_put_vara_double(ncid, var_varid, start, count,      !! Update Data
     +        var_layer1)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
      start(4) = 2
      
      retval = nf_put_vara_double(ncid, var_varid, start, count,      !! Update Data
     +        var_layer2)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = nf_close(ncid)                                         !! Close NetCDF file                
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        print *,'!!! Netcdf file modified i.e. ', FILE_NAME, '!!!'
     
      end subroutine write_variable_netcdf 
      
c ----------------------------------------------------------------------
      
      subroutine read_lagrangian_sigma(file_name,sigma)
      
      implicit none
      include 'netcdf.inc'
      
      character*(*) file_name
      integer nbins
      
      integer retval,ncid,dimid(4),dim_len(4),i
      
      integer sigma_varid
      
      real*8,allocatable,dimension(:,:,:,:) :: sigma
      real*8,allocatable,dimension(:,:,:,:) :: new_sigma
      
      integer count(4),start(4)
      
      character*(*) sigma_name
      parameter(sigma_name = 'Lagrangian Velocity Variance')
      
      retval = nf_open(file_name,nf_nowrite,ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      print*,'sigma file found'
      
      print*,ncid
      
      retval = nf_inq_varid(ncid,sigma_name,sigma_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      print*,'sigma found'
      
      retval = nf_inq_vardimid(ncid,sigma_varid,dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      do i = 1,4
        
        retval = nf_inq_dimlen(ncid,dimid(i),dim_len(i))
        if (retval .ne. nf_noerr) call handle_err(retval)
        
      enddo
      
      print*,'dimension lengths = ',dim_len
      
      do i = 1,4
        
        count(i) = dim_len(i)
        start(i) = 1
        
      enddo
      
      allocate(sigma(dim_len(1),dim_len(2),dim_len(3),dim_len(4)))
      allocate(new_sigma(dim_len(4),dim_len(3),dim_len(2),dim_len(1)))
      
      retval = nf_get_vara_double(ncid,sigma_varid,start,count,sigma)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      new_sigma = reshape(sigma, 
     &     (/dim_len(4),dim_len(3),dim_len(2),dim_len(1)/),
     &      ORDER = (/4,3,2,1/))
      
      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
          
      end subroutine read_lagrangian_sigma
      
c ------------------------------------------------------------------------
c ------------------------------------------------------------------------
      
      subroutine read_pv_lagrangian_sigma(file_name,nbins,sigma)
      
      implicit none
      include 'netcdf.inc'
      
      character*(*) file_name
      integer nbins
      real*8 sigma(nbins)
      
      integer retval,ncid,sigma_varid,start,count
      
      character*(*) sigma_name
      parameter(sigma_name = 'Lagrangian Velocity Variance')
      
      retval = nf_open(file_name,nf_nowrite,ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      print*,'sigma file found'
      
      retval = nf_inq_varid(ncid,sigma_name,sigma_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      print*,'sigma found'
      
      start = 1
      count = nbins
      
      retval = nf_get_vara_double(ncid,sigma_varid,start,count,sigma)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      
      
      end subroutine read_pv_lagrangian_sigma
      
      end module mod_vel_variance_netcdf
