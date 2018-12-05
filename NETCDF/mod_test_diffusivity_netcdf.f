C MODULE THAT WRITES INTERPOLATED DIFFUSIVITY AND DIFFERENTIATION TO FILE

      module mod_test_diffusivity_netcdf
      
      use mod_netcdf_error

      implicit none
      
      contains
      
      subroutine create_diffusivity_file(file_name,ii)
      
      implicit none
      include 'netcdf.inc'
      
      character*(*) file_name
      integer ii
      
      character*(*) diffusivity_name, derivative_name
      parameter(diffusivity_name = 'Diffusivity')
      parameter(derivative_name = 'Derivative')
      
      character*(*),parameter :: layer_name = 'Layer'
      character*(*),parameter :: dim_name = 'Dimension'
      character*(*),parameter :: y_name = 'Grid Point'
      
      integer diffusivity_varid,derivative_varid
      integer y_dimid,dim_dimid,dimid(3),count(3),start(3),layer_dimid
      
      integer retval,ncid
      
      retval = nf_create(file_name,nf_clobber,ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_dim(ncid,layer_name,2,layer_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_dim(ncid,dim_name,2,dim_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_dim(ncid,y_name,ii,y_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      dimid(1) = dim_dimid
      dimid(2) = layer_dimid
      dimid(3) = y_dimid
      
      retval = nf_def_var(ncid,diffusivity_name,nf_double,3,dimid,
     & diffusivity_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)      
      
      retval = nf_def_var(ncid,derivative_name,nf_double,3,dimid,
     & derivative_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)    
      
      retval = nf_enddef(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_close(ncid)                                         !! Close file
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      end subroutine create_diffusivity_file
      
      subroutine write_diffusivity_file(file_name,Kx1,Ky1,Kx2,Ky2,
     & dKdx1,dKdy1,dKdx2,dKdy2,i)
      
      implicit none
      include 'netcdf.inc'
      
      character*(*) file_name
      integer ii,i
      
      character*(*) diffusivity_name, derivative_name
      parameter(diffusivity_name = 'Diffusivity')
      parameter(derivative_name = 'Derivative')
      
      character*(*),parameter :: layer_name = 'Layer'
      character*(*),parameter :: dim_name = 'Dimension'
      character*(*),parameter :: y_name = 'Grid Point'
      
      integer diffusivity_varid,derivative_varid
      integer y_dimid,dim_dimid,dimid(3),count(3),start(3),layer_dimid
      
      integer retval,ncid
      
      real*8 Kx1,Ky1,Kx2,Ky2,dKdx1,dKdy1,dKdx2,dKdy2
      
      retval = nf_open(file_name,nf_write,ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_inq_varid(ncid,diffusivity_name,diffusivity_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_inq_varid(ncid,derivative_name,derivative_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      start(1) = 1
      start(2) = 1
      start(3) = i
      
      count(1) = 1
      count(2) = 1
      count(3) = 1
      
      retval = nf_put_vara_double(ncid,diffusivity_varid
     &    ,start,count,Kx1)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid,derivative_varid
     &    ,start,count,dKdx1)
      if (retval .ne. nf_noerr) call handle_err(retval)
     
      start(1) = 2
      
      retval = nf_put_vara_double(ncid,diffusivity_varid
     &    ,start,count,Ky1)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid,derivative_varid
     &    ,start,count,dKdy1)
      if (retval .ne. nf_noerr) call handle_err(retval)
     
      start(1) = 1
      start(2) = 2
      
      retval = nf_put_vara_double(ncid,diffusivity_varid
     &    ,start,count,Kx2)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid,derivative_varid
     &    ,start,count,dKdx2)
      if (retval .ne. nf_noerr) call handle_err(retval)
     
      start(1) = 2
      
      retval = nf_put_vara_double(ncid,diffusivity_varid
     &    ,start,count,Ky2)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid,derivative_varid
     &    ,start,count,dKdy2)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_close(ncid)                                         !! Close NetCDF file                
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      
      
      end subroutine write_diffusivity_file

      end module mod_test_diffusivity_netcdf
