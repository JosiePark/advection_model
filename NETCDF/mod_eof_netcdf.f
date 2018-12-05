C MODULE THAT READS EOFS, PCS AND CREATES AND WRITES AND READS NEW EOF FIELD FILE
C MODULE THAT READS THE EOF MODES
C USED TO INTERPOLATE ONTO A FINER GRID

      module mod_eof_netcdf
      
      use mod_netcdf_error
      
      implicit none
      
      contains
      
c ---------------------------------------------------------------------
c ---------------------------------------------------------------------
c ---------------------- READ EOF --------------------------------------
      
      subroutine read_eof_netcdf(file_name,ii,jj,nmodes,nlayers,eof)
      
      implicit none
      
      include 'netcdf.inc'
      
      integer ii,jj,nmodes,nlayers,start(4),count(4),i
      real*8 eof(ii,jj,nlayers,nmodes)
      
      character*(*) file_name
      
      character*(*) eof_name
      parameter(eof_name = 'EOFs')
      
      integer retval,ncid
      integer eof_varid
      
      retval = nf_open(FILE_NAME, nf_write, ncid)                     !! Open file in write mode
      if (retval .ne. nf_noerr) call handle_err(retval)
        
      retval = nf_inq_varid(ncid, eof_NAME, eof_varid)                !! Request varid
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      do i = 1,4
            start(i) = 1
      enddo
      
      count(1) = nmodes
      count(2) = ii
      count(3) = jj
      count(4) = nlayers
      
      
      retval = nf_get_vara_double(ncid,eof_varid,start,count,eof)
      if (retval .ne. nf_noerr) call handle_err(retval)
        
      retval = nf_close(ncid)                                         !! Close NetCDF file                
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      !print*,'EOF READ'
      
      end subroutine read_eof_netcdf
      
c ---------------------------------------------------------------------
c ---------------------------------------------------------------------
c ---------------------- READ EOF --------------------------------------
      
      subroutine read_one_eof_netcdf(file_name,ii,jj,mode,nlayers,eof)
      
      implicit none
c reads a select mode
      
      include 'netcdf.inc'
      
      integer ii,jj,mode,nlayers,start(4),count(4),i
      real*8 eof(ii,jj,nlayers)
      
      character*(*) file_name
      
      character*(*) eof_name
      parameter(eof_name = 'EOFs')
      
      integer retval,ncid
      integer eof_varid
      
      retval = nf_open(FILE_NAME, nf_write, ncid)                     !! Open file in write mode
      if (retval .ne. nf_noerr) call handle_err(retval)
        
      retval = nf_inq_varid(ncid, eof_NAME, eof_varid)                !! Request varid
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      start(1) = mode
      
      do i = 2,4
            start(i) = 1
      enddo
      
      count(1) = 1
      count(2) = ii
      count(3) = jj
      count(4) = nlayers
      
      
      retval = nf_get_vara_double(ncid,eof_varid,start,count,eof)
      if (retval .ne. nf_noerr) call handle_err(retval)
        
      retval = nf_close(ncid)                                         !! Close NetCDF file                
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      !print*,'EOF READ'
      
      end subroutine read_one_eof_netcdf

c ---------------------------------------------------------------------
c ---------------------------------------------------------------------
C ----------------------- CREATE EOF FILE ------------------------------
      
      subroutine create_eof_netcdf(file_name,ii,jj,nlayers,nt)
      
      implicit none
      
      include 'netcdf.inc'
      
      integer ii,jj,nlayers,nt
      
      character*(*) file_name
      
      character*(*) eof_name,lvl_name,x_name,y_name,t_name
      parameter(eof_name = 'Stream Function')
      parameter(lvl_name = 'Layer')
      parameter(x_name = 'X')
      parameter(y_name = 'Y')
      parameter(t_name = 'Time')
      
      integer lvl_dimid,x_dimid,y_dimid,dimids(4),ndims
     & ,t_dimid
      parameter(ndims=4)
      
      integer ncid,retval,eof_varid
      
      retval = nf_create(FILE_NAME, nf_clobber, ncid)                 !! Create netcdf file
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_dim(ncid, Lvl_NAME, nlayers, lvl_dimid)           !! Define Dimensions in NetCDF file
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid, X_NAME, ii, x_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid, Y_NAME, jj, y_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid, t_NAME, nt, t_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      dimids(1) = x_dimid
      dimids(2) = y_dimid
      dimids(3) = lvl_dimid
      dimids(4) = t_dimid
        
        retval = nf_def_var(ncid, eof_NAME, NF_DOUBLE, NDIMS, dimids,   !! Creating Main varialbes
     +        eof_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
      retval = nf_close(ncid)                                         !! Close file
      if (retval .ne. nf_noerr) call handle_err(retval)
   
      print *,'!!! Netcdf file created i.e. ', FILE_NAME, '!!!'
      print *,'!!!                                         !!!'
 
      end subroutine create_eof_netcdf
      
      
c ---------------------------------------------------------------------
c ---------------------------------------------------------------------
C -----------------------WRITE TO EOF FILE ----------------------------
      
      subroutine write_eof_netcdf(file_name,ii,jj,nlayers,time,eof)
      
      implicit none
      
      include 'netcdf.inc'
      
      character*(*) file_name
      integer ii,jj,time,nlayers,i
      real*8 eof(ii,jj,nlayers)
      
      character*(*) eof_name
      parameter(eof_name = 'Stream Function')
      
      integer ncid,retval
      
      integer eof_varid,start(4),count(4)
      
      retval = nf_open(FILE_NAME, nf_write, ncid)                     !! Open file in write mode
      if (retval .ne. nf_noerr) call handle_err(retval)
        
      retval = nf_inq_varid(ncid, eof_NAME, eof_varid)                !! Request varid
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      !print*,'size eof in entcdf = ',size(eof)
      
      do i = 1,3
        start(i) = 1
      enddo
      
      start(4) = time
      
      count(1) = ii
      count(2) = jj 
      count(3) = nlayers
      count(4) = 1
      
      retval = nf_put_vara_double(ncid, eof_varid, start, count,      !! Update Data
     +        eof)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
        retval = nf_close(ncid)                                         !! Close NetCDF file                
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        !print *,'!!! Netcdf file modified i.e. ', FILE_NAME, '!!!'
        !print *,'!!!                                         !!!'
      
      end subroutine write_eof_netcdf
      
c ---------------------------------------------------------------------
c ---------------------------------------------------------------------
C ----------------------- READ PC ----------------------------

      subroutine read_pc_netcdf(file_name,time,nmodes,PC)
      
      implicit none
      
      include 'netcdf.inc'
      
      integer ii,jj,nmodes,start(2),count(2),i,time
      real*8 PC(nmodes)
      
      character*(*) file_name
      
      character*(*) PC_name
      parameter(PC_name = 'PCs')
      
      integer retval,ncid
      integer pc_varid
      
      retval = nf_open(FILE_NAME, nf_write, ncid)                     !! Open file in write mode
      if (retval .ne. nf_noerr) call handle_err(retval)
        
      retval = nf_inq_varid(ncid, pc_NAME, pc_varid)                !! Request varid
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      start(1) = time
      start(2) = 1
      
      count(1) = 1
      count(2) = nmodes

      
      retval = nf_get_vara_double(ncid,pc_varid,start,count,pc)
      if (retval .ne. nf_noerr) call handle_err(retval)
        
      retval = nf_close(ncid)                                         !! Close NetCDF file                
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      !print*,'PC READ'
      
      end subroutine read_pc_netcdf
      
c ---------------------------------------------------------------------
c ---------------------------------------------------------------------
C ----------------------- READ PC ----------------------------

      subroutine read_one_pc_netcdf(file_name,time,mode,PC)
      
      implicit none
      
      include 'netcdf.inc'
      
      integer ii,jj,mode,start(2),count(2),i,time
      real*8 PC
      
      character*(*) file_name
      
      character*(*) PC_name
      parameter(PC_name = 'PCs')
      
      integer retval,ncid
      integer pc_varid
      
      retval = nf_open(FILE_NAME, nf_write, ncid)                     !! Open file in write mode
      if (retval .ne. nf_noerr) call handle_err(retval)
        
      retval = nf_inq_varid(ncid, pc_NAME, pc_varid)                !! Request varid
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      start(1) = time
      start(2) = mode
      
      count(1) = 1
      count(2) = 1

      
      retval = nf_get_vara_double(ncid,pc_varid,start,count,pc)
      if (retval .ne. nf_noerr) call handle_err(retval)
        
      retval = nf_close(ncid)                                         !! Close NetCDF file                
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      !print*,'PC READ'
      
      end subroutine read_one_pc_netcdf
      
      end module mod_eof_netcdf
      
