c module that writes the PV mapped dispersion

      module mod_PVdisp_netcdf
      
      use mod_netcdf_error
      
      contains
      
      subroutine create_disp_netcdf(file_name,npoints,nrel)
      
      implicit none
      include 'netcdf.inc'
      
      character*(*) file_name
      integer npoints,nrel
      
      integer ncid,retval
      
      character*(*) disp_name,part_name,rel_name,t_name,PV_name
      parameter(disp_name = 'PV mapped Y')
      parameter(PV_name = 'PV at particle')
      parameter(part_name = 'Particle Number')
      parameter(rel_name = 'Release Number')
      parameter(t_name = 'Time')
      
      integer disp_dimid,part_dimid,rel_dimid,t_dimid,PV_dimid
      integer disp_varid,part_varid,rel_varid,t_varid,PV_varid
      
      integer ndims
      parameter(ndims = 3)
      
      integer dimid(ndims)
c create NETCDF file
      retval = nf_create(file_name,nf_clobber,ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
c Define dimensions in the NETCDF files
      retval = nf_def_dim(ncid,part_name,npoints,part_dimid) ! defines particle numbers
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_dim(ncid,t_name,nf_unlimited,t_dimid) ! time, allowed to go to infinity
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_dim(ncid,rel_name,nrel,rel_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      dimid(1) = part_dimid
      dimid(2) = rel_dimid
      dimid(3) = t_dimid
      
      retval = nf_def_var(ncid,disp_name,nf_double,ndims,dimid
     &      ,disp_varid)
      if (retval .ne. nf_noerr) then 
      print *, 'Y_PV not created'
      call handle_err(retval)
      endif
      
      retval = nf_def_var(ncid,PV_name,nf_double,ndims,dimid
     &      ,PV_varid)
      if (retval .ne. nf_noerr) then 
      print *, 'PV not created'
      call handle_err(retval)
      endif
      
      retval = nf_enddef(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_close(ncid)                                         !! Close file
      if (retval .ne. nf_noerr) call handle_err(retval)
   
      end subroutine
      
c ------------------------------------------------------------------------
c ------------------------------------------------------------------------
      
      subroutine write_disp_netcdf(file_name,npoints,k
     &      ,Y_map,PV,step_time)
      
      implicit none
      include 'netcdf.inc'
      
      character*(*) file_name
      integer npoints,k,step_time
      
      real*8 Y_map(npoints),PV(npoints)
      
      integer ncid,retval
      
      character*(*) disp_name,part_name,rel_name,t_name,PV_name
      parameter(PV_name = 'PV at particle')
      parameter(disp_name = 'PV mapped Y')
      parameter(part_name = 'Particle Number')
      parameter(rel_name = 'Release Number')
      parameter(t_name = 'Time')
      
      integer disp_dimid,part_dimid,rel_dimid,t_dimid,PV_dimid
      integer disp_varid,part_varid,rel_varid,t_varid,PV_varid
      
      integer ndims
      parameter(ndims = 3)
      integer count(ndims),start(ndims)
      
      integer dimid(ndims)
      
c open file in read mode
      retval = nf_open(file_name, nf_write,ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_inq_varid(ncid,disp_name,disp_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_inq_varid(ncid,PV_name,PV_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      count(1) = npoints
      count(2) = 1
      count(3) = 1
      
      start(1) = 1
      start(2) = k
      start(3) = step_time
      
      retval = nf_put_vara_double(ncid,disp_varid,start,count,Y_map)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_put_vara_double(ncid,PV_varid,start,count,PV)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      print*, 'NETCDF FILE EDITED'
      
      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
    
      
      
      end subroutine
      
      
c ------------------------------------------------------------------------------
c ----------------- SUBROUTINE THAT CREATES BINNED PV DISPERSION FILE ----------
C ------------------------------------------------------------------------------

      subroutine create_binned_PVdisp_file(file_name,npoints,nbins,nrel
     & ,bin_boundaries)
      
      implicit none
      include 'netcdf.inc'
      
      character*(*) file_name
      integer npoints,nbins,nrel
      
      integer ncid,retval
      
      real*8 bin_boundaries(nbins)
      
      character*(*) disp_name,part_name,rel_name,t_name,PV_name,bin_name
     & ,bin
      parameter(disp_name = 'PV mapped Y')
      parameter(PV_name = 'PV at particle')
      parameter(part_name = 'Particle Number')
      parameter(rel_name = 'Release Number')
      parameter(t_name = 'Time')
      parameter(bin_name = 'Bin Number')
      parameter(bin = 'Bin Boundaries')
      
      integer disp_dimid,part_dimid,rel_dimid,t_dimid,PV_dimid,bin_dimid
      integer disp_varid,part_varid,rel_varid,t_varid,PV_varid,bin_varid
      
      integer ndims
      parameter(ndims = 4)
      
      integer dimid(ndims)
c create NETCDF file
      retval = nf_create(file_name,nf_clobber,ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
c Define dimensions in the NETCDF files
      retval = nf_def_dim(ncid,part_name,npoints,part_dimid) ! defines particle numbers
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_dim(ncid,t_name,nf_unlimited,t_dimid) ! time, allowed to go to infinity
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_dim(ncid,rel_name,nrel,rel_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_dim(ncid,bin_name,nbins,bin_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      dimid(1) = bin_dimid
      dimid(2) = part_dimid
      dimid(3) = rel_dimid
      dimid(4) = t_dimid
      
      retval = nf_def_var(ncid,disp_name,nf_double,ndims,dimid
     &      ,disp_varid)
      if (retval .ne. nf_noerr) then 
      print *, 'Y_PV not created'
      call handle_err(retval)
      endif
      
      retval = nf_def_var(ncid,PV_name,nf_double,ndims,dimid
     &      ,PV_varid)
      if (retval .ne. nf_noerr) then 
      print *, 'PV not created'
      call handle_err(retval)
      endif
      
      retval = nf_def_var(ncid,bin,nf_double,1,bin_dimid,bin_varid)
            if (retval .ne. nf_noerr) then 
      print *, 'Bin boundaries not created'
      call handle_err(retval)
      endif

      retval = nf_enddef(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_close(ncid)                                         !! Close file
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_open(file_name, nf_write,ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_put_vara_double(ncid,bin_varid,1,1,bin_boundaries)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      end subroutine create_binned_PVdisp_file
      
c ----------------------------------------------------------------------------
c ---------  SUBROUTINE THAT WRITES THE PV MAPPED DISPERSION ---------------
C ------------------------- IN BINS TO THE NETCDF FILE ------------------
C -----------------------------------------------------------------------------

      subroutine write_binned_PVdisp_file(file_name,npoints,k,b
     &      ,Y_map,PV,step_time)
      
      implicit none
      include 'netcdf.inc'
      
      character*(*) file_name
      integer npoints,k,step_time,b
      
      real*8 Y_map(npoints),PV(npoints)
      
      integer ncid,retval
      
      character*(*) disp_name,part_name,rel_name,t_name,PV_name,bin_name
      parameter(PV_name = 'PV at particle')
      parameter(bin_name = 'Bin Number')
      parameter(disp_name = 'PV mapped Y')
      parameter(part_name = 'Particle Number')
      parameter(rel_name = 'Release Number')
      parameter(t_name = 'Time')
      
      integer disp_dimid,part_dimid,rel_dimid,t_dimid,PV_dimid,bin_dimid
      integer disp_varid,part_varid,rel_varid,t_varid,PV_varid
      
      integer ndims
      parameter(ndims = 4)
      integer count(ndims),start(ndims)
      
      integer dimid(ndims)
      
c open file in read mode
      retval = nf_open(file_name, nf_write,ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_inq_varid(ncid,disp_name,disp_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_inq_varid(ncid,PV_name,PV_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      count(1) = 1
      count(2) = npoints
      count(3) = 1
      count(4) = 1
      
      start(1) = b
      start(2) = 1
      start(3) = k
      start(4) = step_time
      
      retval = nf_put_vara_double(ncid,disp_varid,start,count,Y_map)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_put_vara_double(ncid,PV_varid,start,count,PV)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      print*, 'NETCDF FILE EDITED'
      print*, 'bin = ',b
      print*, 'k=',k
      print*, 't=',step_time
      print*,'Y_map(1) = ',Y_map(1)
      
      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      end subroutine write_binned_PVdisp_file
      
      subroutine read_PV_binned_trajectories(file_name,npoints,k,p
     & ,step_time
     & ,x1,y1
     & ,x1_coord,y1_coord)

c reads the bin_no bin
c and the rel_no release
c and the time index, t_index
c returns the trajectory
c and domain coordinate
      
      implicit none
      include 'netcdf.inc'
      
      
      character*(*) file_name
      
      integer ncid,retval,npoints,nbins,p,step_time,k
      
      real*8 x1(npoints),y1(npoints)
      integer x1_coord(npoints),y1_coord(npoints)
      real*8 x2(npoints),y2(npoints)
      integer x2_coord(npoints),y2_coord(npoints)
      
c Dimension names
      
      character*(*) t_name,part_name,bin_name,rel_name,layer_name
     & ,loc_name,rec_name,coor_name,start_name
      parameter(t_name = 'Time')
      parameter(part_name = 'Particle Number')
      parameter(bin_name = 'Bin Number')
      parameter(rel_name = 'Release Number')
      parameter(layer_name = 'Layer')
      parameter(loc_name = 'Location of Particle')
      parameter(rec_name = 'Record Time')
      parameter(coor_name = 'Domain Coordinate')
      parameter(start_name = 'Restart Point')
      
c Variable names
      
      character*(*) bin,traj_name,coord_name
      parameter(bin = 'Bin Boundaries')
      parameter(traj_name = 'Trajectories')
      parameter(coord_name = 'Domain Coordinate')
      
c Dimension IDs

      integer bin_dimid,rel_dimid,t_dimid,part_dimid,layer_dimid
     & ,loc_dimid,rec_dimid,coor_varid,start_varid
      
c Variable ids

      integer traj_varid,t_varid
      
      integer ndims
      parameter(ndims = 6)
      
      integer start(ndims),count(ndims)
      
      retval = nf_open(FILE_NAME, nf_nowrite, ncid)                   !! Open file in read mode
      if (retval .ne. nf_noerr) call handle_err(retval)
              

        retval = nf_inq_varid(ncid,traj_name,traj_varid)
        if (retval .ne. nf_noerr) then
        print *, 'trajectories not found'
        call handle_err(retval)
        endif

        retval = nf_inq_varid(ncid,t_name,t_varid)
        if (retval .ne. nf_noerr) then 
        print *, 'time not found'
        call handle_err(retval)
        endif

        retval = nf_inq_varid(ncid,coor_name,coor_varid)
        if (retval .ne. nf_noerr) then 
        print *, 'coordinate not found'
        call handle_err(retval)
        endif
        
        retval = nf_inq_varid(ncid,start_name,start_varid)
        if (retval .ne. nf_noerr) then 
        print *, 'start info not found'
        call handle_err(retval)
        endif
        
      count(1) = 1
      count(2) = 1
      count(3) = npoints
      count(4) = 1
      count(5) = 1
      count(6) = 1
      start(1) = p
      start(2) = 1
      start(3) = 1
      start(4) = k
      start(5) = 1
      start(6) = step_time ! tells the program to start writing from step time, in the time variable, but writes in all columns in the location and particle numbers
      
      retval = nf_get_vara_double(ncid,traj_varid,start,count,x1)
      if (retval .ne. nf_noerr) call handle_err(retval)
      !print*,'found traj'
      retval = nf_get_vara_int(ncid,coor_varid,start,count,x1_coord)
      if (retval .ne. nf_noerr) call handle_err(retval)
      !print*,'found coord'
      
      
      start(2) = 2
      retval = nf_get_vara_double(ncid,traj_varid,start,count,y1)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_get_vara_int(ncid,coor_varid,start,count,y1_coord)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      start(5) = 2
      retval = nf_get_vara_double(ncid,traj_varid,start,count,y2)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_get_vara_int(ncid,coor_varid,start,count,y2_coord)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      start(2) = 1
      retval = nf_get_vara_double(ncid,traj_varid,start,count,x2)
      if (retval .ne. nf_noerr) call handle_err(retval)
      !print*,'found traj'
      retval = nf_get_vara_int(ncid,coor_varid,start,count,x2_coord)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      
      
        retval = nf_close(ncid)                                         !! Close file
        if (retval .ne. nf_noerr) call handle_err(retval)
      
      
      end subroutine read_PV_binned_trajectories
      
      
      
      end module
