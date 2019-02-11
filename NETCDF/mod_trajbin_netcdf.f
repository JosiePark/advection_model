C MODULE THAT CONTAINS SUBROUTINES THAT READ TRAJECTORIES A BIN AT A TIME

      module mod_trajbin_netcdf
      
      use mod_netcdf_error
      
      contains
      
      subroutine read_trajbin_file(file_name
     & ,bin_no,rel_no,npoints,nt,ii,traj)
     
      implicit none
      include 'netcdf.inc'
     
      integer bin_no,npoints,rel_no,nt,ii
      character*(*) file_name
      real*8 traj(2,npoints,2,nt),tmp(2,npoints,2,nt)
      integer coord(2,npoints,2,nt)
      
      integer ndims
      parameter(ndims = 6)
      integer start(ndims),count(ndims)
      integer ncid,retval
      
c Dimension names
      
      character*(*) t_name,part_name,bin_name,rel_name,layer_name
     & ,loc_name,rec_name,t_units,x_units,y_units,units
      parameter(t_name = 'Time')
      parameter(part_name = 'Particle Number')
      parameter(bin_name = 'Bin Number')
      parameter(rel_name = 'Release Number')
      parameter(layer_name = 'Layer')
      parameter(loc_name = 'Location of Particle')
      parameter(rec_name = 'Record Time')
      parameter(t_units = 'Days')
      parameter(x_units = 'Grid Points')
      parameter(y_units = 'Grid Point')
      parameter(units = 'units')
      
c Variable names
      
      character*(*) bin,traj_name,coord_name,start_name
      parameter(bin = 'Bin Boundaries')
      parameter(traj_name = 'Trajectories')
      parameter(coord_name = 'Domain Coordinate')
      parameter(start_name = 'Restart Point')
      
c variable ids
      
      integer traj_varid,coord_varid,t_varid
      
      print*,'start'
      
c Open netcdf file in read mode
      retval = nf_open(file_name,nf_nowrite,ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      print*,'file_opened'
      
c request the variable ids of the variables in the file
      retval = nf_inq_varid(ncid,traj_name,traj_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_inq_varid(ncid,t_name,t_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)  
      
      retval = nf_inq_varid(ncid,coord_name,coord_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
c Now tell the program where to read the data in the NETCDF file
      count(1) = 1
      count(2) = 2
      count(3) = npoints
      count(4) = 1
      count(5) = 2
      count(6) = nt
      start(1) = bin_no
      start(2) = 1
      start(3) = 1
      start(4) = rel_no ! tells the program to start writing at the release number
      start(5) = 1
      start(6) = 1
      
      retval = nf_get_vara_double(ncid,traj_varid,start,count,tmp)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_get_vara_int(ncid,coord_varid,start,count,coord)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      ! calculate trajectory and return
      
      traj = tmp + dble(ii*coord)
      
      retval = nf_close(ncid)                                         !! Close file
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      end subroutine read_trajbin_file
      
      
      end module mod_trajbin_netcdf
