      module mod_kinematic_netcdf
      
      use mod_netcdf_error
      
      contains
      
c ---------------------------------------------------------------------
      
      subroutine create_kinematic_field(file_name,ii,jj)
      
      implicit none
      include 'netcdf.inc'
      
      character*(*) file_name
      integer ii,jj
      
      integer ncid,t_dimid,x_dimid,y_dimid,retval
      integer ndims
      parameter (ndims = 3)
      integer dimid(ndims)
      integer t_varid,psi_varid
      character*(*) x_name,y_name,psi_name,t_name
      parameter(x_name = 'X',y_name = 'Y',t_name = 'Time')
      parameter(psi_name = 'Stream Function')
      
      retval = nf_create(file_name,nf_clobber,ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_dim(ncid,x_name,ii,x_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid,y_name,jj,y_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid,t_name,nf_unlimited,t_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      dimid(1) = x_dimid
      dimid(2) = y_dimid
      dimid(3) = t_dimid
      
      retval = nf_def_var(ncid,t_name,nf_double,1,t_dimid,t_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid,psi_name,nf_double,ndims,dimid,psi_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_enddef(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      print*,'NETCDF file created', file_name     
      
      end subroutine create_kinematic_field
      
c ---------------------------------------------------------------------

      subroutine write_kinematic_field(file_name,ii,jj,psi,nrec,time)
      
      implicit none
      include 'netcdf.inc'
      
      character*(*) file_name
      integer ii,jj
      real*8 psi(ii,jj)
      integer nrec
      real*8 time
      
      integer ncid,retval,ndims
      parameter(ndims = 3)
      integer dimid(ndims),start(ndims),count(ndims)
      
      integer x_dimid,y_dimid,t_dimid
      integer t_varid,psi_varid
      
      character*(*) x_name,y_name,t_name,psi_name
      parameter(x_name = 'X',y_name = 'Y', t_name = 'Time')
      parameter(psi_name = 'Stream Function')
      
      retval = nf_open(file_name,nf_write,ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_inq_varid(ncid,psi_name,psi_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      print*,'Stream Function found'
      retval = nf_inq_varid(ncid,t_name,t_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      print*,' Time found'
      
      start(1) = 1
      start(2) = 1
      start(3) = nrec
      count(1) = ii
      count(2) = jj
      count(3) = 1
      
      retval = nf_put_vara_double(ncid,psi_varid,start,count,psi)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_put_vara_double(ncid,t_varid,nrec,1,time)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      print*,'NETCDF file modified', file_name
      
      end subroutine write_kinematic_field
      
c ---------------------------------------------------------------------

      subroutine create_kinematic_traj(file_name,npoints)
      
      implicit none
      include 'netcdf.inc'
      
      character*(*) file_name
      integer npoints
      
      integer ncid,retval
      
c Dimension IDs

      integer t_dimid, p_dimid, loc_dimid, rec_dimid
      integer t_varid, traj_varid
      character*(*) t_name, p_name, loc_name, traj_name, rec_name
      parameter(t_name = 'Time (days)')
      parameter(p_name = 'Particle Number')
      parameter(loc_name = 'Location of Particle')
      parameter(rec_name = 'Time')
      parameter(traj_name = 'Trajectories')
      
      integer ndims
      parameter(ndims = 3)
      integer dimid(ndims)
      
c create NETCDF file

      retval = nf_create(file_name,nf_clobber,ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
c Define Dimensions

      retval = nf_def_dim(ncid,p_name,npoints,p_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_dim(ncid,loc_name,2,loc_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_dim(ncid,rec_name,nf_unlimited,rec_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
c Define Variables

      dimid(1) = p_dimid
      dimid(2) = loc_dimid
      dimid(3) = rec_dimid

      retval = nf_def_var(ncid,traj_name,nf_double
     & ,ndims,dimid,traj_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_var(ncid,t_name,nf_double,1,rec_dimid,t_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_enddef(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)


      end subroutine create_kinematic_traj
      
c ---------------------------------------------------------------------

      subroutine write_kinematic_traj(file_name,npoints,x,y,time
     & ,step_time)
      
      implicit none
      include 'netcdf.inc'
      
      character*(*) file_name
      integer npoints
      real*8 x(npoints),y(npoints),time,traj(npoints,2)
      integer step_time
      
      integer ncid,retval,ndims
      parameter(ndims = 3)
      integer dimid(ndims)
      
      integer traj_varid,t_varid
      
      character*(*) traj_name,t_name
      parameter(traj_name = 'Trajectories')
      parameter(t_name = 'Time (days)')
      
      integer start(ndims),count(ndims)
      
      retval = nf_open(file_name,nf_write,ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
c Inquire variable IDs
      
      retval = nf_inq_varid(ncid,traj_name,traj_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_inq_varid(ncid,t_name,t_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
c Update data

      start(1) = 1
      start(2) = 1
      start(3) = step_time
      count(1) = npoints
      count(2) = 2
      count(3) = 1
      
      traj(:,1) = x
      traj(:,2) = y
      
      retval = nf_put_vara_double(ncid,traj_varid,start,count,traj)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_put_vara_double(ncid,t_varid,step_time,1,time)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      print*,'NETCDF file created', file_name   
      
      end subroutine write_kinematic_traj
      
      end module mod_kinematic_netcdf
