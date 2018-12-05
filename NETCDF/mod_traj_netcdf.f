      module MOD_traj_NETCDF

       use mod_netcdf_error
       
      contains
      
      
c ---------------------------------------------------------------------
c ----------- CREATE NETCDF FILE THAT STORES TRAJECTORY FOR A SINGLE RELEASE ----------
C ----------------------------------------------------------------------------------

      subroutine create_traj_file(file_name,npoints)
c Creates a file called file_name, and takes the number of lagrangian particles as the input
      implicit none
      include 'netcdf.inc'
      character*(*) file_name
      integer npoints, ncid, retval,nlvls
      parameter(nlvls = 2)

      
c Define names of variables
      character*(*) t_name, loc_name, part_name,traj_name
      character*(*) coor_name,lvl_name
      parameter (t_name = 'Time',loc_name ='Location of Particle')
      parameter(part_name = 'Particle Number')
      parameter(traj_name = 'Trajectories')
      parameter(coor_name = 'Domain Coordinate')
      parameter(lvl_name = 'Layer')

      
c Define units
      character*(*) t_units, x_units, y_units, units
      parameter (t_units = 'Days',x_units = 'cms',y_units = 'cms')
      parameter (units = 'units')
      
c initilisae dimension ids
      integer t_dimid, loc_dimid, part_dimid, traj_dimid, coor_dimid
      integer lvl_dimid

c initialise variable ids
      integer t_varid, loc_varid,part_varid, traj_varid,coor_varid

      
c number of dimensions in space(2: x & y)
      integer spatial_dims
      parameter (spatial_dims = 2)
c number of dimensions in the particle location variables : space, particle number, time, layer
      integer ndims
      parameter (ndims = 4)
      integer dimid(ndims)

      
c create NETCDF file
      retval = nf_create(file_name,nf_clobber,ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
c Define dimensions in the NETCDF files
      retval = nf_def_dim(ncid,part_name,npoints,part_dimid) ! defines particle numbers
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_dim(ncid,t_name,nf_unlimited,t_dimid) ! time, allowed to go to infinity
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_dim(ncid,loc_name,spatial_dims,loc_dimid) ! spatial dimensions for the location
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_dim(ncid,lvl_name,nlvls,lvl_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
c Define variables
      retval = nf_def_var(ncid,t_name,nf_double,1,t_dimid,t_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_put_att_text(ncid,t_varid,units,len(t_units)
     & ,t_units)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      dimid(1) = loc_dimid
      dimid(2) = part_dimid
      dimid(3) = lvl_dimid
      dimid(4) = t_dimid
      
      
c Define Trajectory matrix, which stores locations of all npoints lagrangian particles at time t

      retval = nf_def_var(ncid,traj_name,nf_double,ndims,dimid
     & ,traj_varid)
      if (retval .ne. nf_noerr) then 
      print *, 'traj not created'
      call handle_err(retval)
      endif
      
      retval = nf_def_var(ncid,coor_name,nf_int,ndims,dimid,
     & coor_varid)
      if (retval .ne. nf_noerr) call handle_err(retval) 
      
      retval = nf_put_att_text(ncid,traj_varid,units
     & ,len(x_units),x_units)
      if (retval .ne. nf_noerr) call handle_err(retval)
    
      
      retval = nf_enddef(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      
     
      
      retval = nf_close(ncid)                                         !! Close file
      if (retval .ne. nf_noerr) call handle_err(retval)
   
      print *,'!!! Netcdf file created i.e. ', FILE_NAME, '!!!'
      print *,'!!!                                         !!!'

      end subroutine create_traj_file
      
c -------------------------------------------------------------------------------------- 
c -------------------- SUBROUTINE TO WRITE TRAJECTORIES TO A NETCDF FILE ------------------------------
      subroutine write_traj_file(file_name,npoints,x1,y1,x2,y2
     & ,x1_coord
     & ,y1_coord,x2_coord,y2_coord,save_time
     & ,step_time)
c writes the new particle locations x,y to the NETCDF file at the save_time
      implicit none
      include 'netcdf.inc'
      
      character*(*) file_name
      integer npoints,retval,ncid,i,j, step_time,traj_varid,eddy_varid
     & ,x1_coord(npoints),y1_coord(npoints),coor_varid,coor_eddy_varid
     & ,x2_coord(npoints),y2_coord(npoints)
      real*8 x1(npoints),y1(npoints), save_time,x2(npoints),y2(npoints)
      
c Temporary trajectory variable
      real*8 traj_temp(2,npoints,2)
      
      integer coor(2,npoints,2)
      
      integer nlvls
      parameter(nlvls = 2)
      
c Define names
      character*(*) traj_name, t_name,coor_name,lvl_name
      parameter(traj_name = 'Trajectories')
      parameter(t_name = 'Time')
      parameter(coor_name = 'Domain Coordinate')
      parameter(lvl_name = 'Layer')

      
c number of dimensions in space(2: x & y)
      integer spatial_dims, t_varid, lvl_dimid
      parameter (spatial_dims = 2)
c number of dimensions in the particle location variables : location, coordinate, particle number, time
      integer ndims
      parameter (ndims = 4)
      integer dimid(ndims)
      integer start(ndims), count(ndims)
      
c open file in read mode
      retval = nf_open(file_name, nf_write,ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
c request the variable ids of the variables in the file
      retval = nf_inq_varid(ncid,traj_name,traj_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_inq_varid(ncid,t_name,t_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)  
      
      retval = nf_inq_varid(ncid,coor_name,coor_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      

c Write x, y to the temporary trajectory variable
      do j = 1,npoints
              traj_temp(1,j,1) = x1(j)
              traj_temp(2,j,1) = y1(j)
              traj_temp(1,j,2) = x2(j)
              traj_temp(2,j,2) = y2(j)
              coor(1,j,1) = x1_coord(j)
              coor(2,j,1) = y1_coord(j)
              coor(1,j,2) = x2_coord(j)
              coor(2,j,2) = y2_coord(j)
      end do
      
c Now tell the program where to write the data in the NETCDF file
      count(1) = spatial_dims
      count(2) = npoints
      count(3) = nlvls
      count(4) = 1
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = step_time ! tells the program to start writing from step time, in the time variable, but writes in all columns in the location and particle numbers

c Update data
      retval = nf_put_vara_double(ncid,traj_varid,start,count,traj_temp)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_int(ncid,coor_varid,start,count,coor)
      if (retval .ne. nf_noerr) call handle_err(retval)
      

c Update time
      retval = nf_put_vara_double(ncid, t_varid, step_time, 1,
     +       save_time)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_close(ncid)                                         !! Close NetCDF file                
      if (retval .ne. nf_noerr) call handle_err(retval)
        
      print *,'!!! Netcdf file modified i.e. ', FILE_NAME, '!!!'
      print *,'!!!                                         !!!'
      end subroutine write_traj_file
      
      
c --------------------------------------------------------------
c SUBROUTINE TO READ TRAJECTORY FILE FROM READ TIME ------------
        subroutine read_traj_file(file_name,npoints,x1,y1,x2,y2
     &   ,x1_coord,y1_coord,x2_coord,y2_coord 
     & ,read_tim,new_tim,step_tim)
     
        implicit none
        include 'netcdf.inc'
        
        integer npoints,x1_coord(npoints),y1_coord(npoints),step_tim
        integer x2_coord(npoints),y2_coord(npoints)
        character*(*) file_name
        real*8 x1(npoints),y1(npoints),read_tim,new_tim
        real*8 x2(npoints),y2(npoints)
        integer retval,ncid, traj_varid, t_varid, traj_dimid, t_dimid
     & ,coor_varid,coor_dimid,t_len,i
        real *8, dimension (:), allocatable :: tp_time
        
        character*(*) traj_name, t_name,coor_name,lvl_name
        parameter(traj_name = 'Trajectories')
        parameter(t_name = 'Time')
        parameter(coor_name = 'Domain Coordinate')
        parameter(lvl_name = 'Layer')
        
c number of dimensions in space(2: x & y)
      integer spatial_dims
      parameter (spatial_dims = 2)
      integer nlvls
      parameter (nlvls = 2)
c number of dimensions in the particle location variables : location, coordinate, particle number, time
      integer ndims
      parameter (ndims = 4)
      integer dimid(ndims)
      integer start(ndims), count(ndims)
        
        retval = nf_open(file_name,nf_nowrite,ncid) !open file in read mode
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
        
        retval = nf_inq_dimid(ncid,traj_name,traj_dimid)
        if (retval .ne. nf_noerr) then 
        print *, 'traj dimid not found'
        call handle_err(retval)
        endif
        retval = nf_inq_dimid(ncid,t_name,t_dimid)
        if (retval .ne. nf_noerr)then
        print *, 'time dimid not found'
         call handle_err(retval)
         endif
        retval = nf_inq_dimid(ncid,coor_name,coor_dimid)
        if (retval .ne. nf_noerr) then 
        print *,'coordinate dimid not found'
        call handle_err(retval)
        endif
        
        allocate(tp_time(t_len))
        retval = nf_get_vara_double(ncid,t_varid,1,t_len,tp_time)
        if (retval .ne. nf_noerr) then 
        print *, 'time not read'
        call handle_err(retval)
        endif
        
        if(read_tim.ne.0.) then
          if(read_tim.gt.tp_time(t_len)) then
            print *,'Maximum time saved at t = ',tp_time(t_len),' Days'
            print *,'Starting simulation from this time .......' 
            new_tim = tp_time(t_len)
            step_tim = t_len
          else
            do i = 1, t_len
              if(tp_time(i).ge.read_tim) then
                step_tim = i
                new_tim = tp_time(i)
                print *,'Starting simulation from t = ',new_tim,' Days'
                stop
              end if
            end do
          end if
        else
          new_tim = tp_time(t_len)
          step_tim = t_len
          print *,'Starting simulation from t = ',new_tim,' Days'
        end if  
        
      count(1) = 1
      count(2) = npoints
      count(3) = 1
      count(4) = 1
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = step_tim ! tells the program to start writing from step time, in the time variable, but writes in all columns in the location and particle numbers
      
      retval = nf_get_vara_double(ncid,traj_varid,start,count,x1)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_get_vara_int(ncid,coor_varid,start,count,x1_coord)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      
      start(1) = 2
      retval = nf_get_vara_double(ncid,traj_varid,start,count,y1)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_get_vara_int(ncid,coor_varid,start,count,y1_coord)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      start(1) = 1
      start(3) = 2
      
      retval = nf_get_vara_double(ncid,traj_varid,start,count,x2)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_get_vara_int(ncid,coor_varid,start,count,x2_coord)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      start(1) = 2
      
      retval = nf_get_vara_double(ncid,traj_varid,start,count,y2)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_get_vara_int(ncid,coor_varid,start,count,y2_coord)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      
        retval = nf_close(ncid)                                         !! Close file
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        deallocate (tp_time)
        
       end subroutine read_traj_file
 
c ---------------------------------------------------------------------------
c ---------------------------------------------------------------------------
c SUBROUTINE THAT CREATES NETCDF FILES THAT STORE THE LOCATIONS OF LAGRANGIAN PARTICLES AT DIFFERENT RELEASE TIMES

      subroutine create_release_file(file_name,npoints,release_no)
c npoints is the number of lagrangian particles per release
c release_no is the number of releases

      implicit none
      include 'netcdf.inc'
      
      character*(*) file_name
      
      integer npoints,release_no,retval,ncid
      
c Define names of variables
      character*(*) t_name, loc_name, part_name,traj_name
      character*(*) coor_name,lvl_name,rel_name,rec_name
      parameter (t_name = 'Time',loc_name ='Location of Particle')
      parameter(part_name = 'Particle Number')
      parameter(traj_name = 'Trajectories')
      parameter(coor_name = 'Domain Coordinate')
      parameter(lvl_name = 'Layer')
      parameter(rel_name = 'Release Number')
      parameter(rec_name = 'Record Time')
      
c Define units
      character*(*) t_units, x_units, y_units, units
      parameter (t_units = 'Days',x_units = 'grid point'
     &      ,y_units = 'grid point')
      parameter (units = 'units')
c initilisae dimension ids
      integer t_dimid, loc_dimid, part_dimid, traj_dimid, coor_dimid
      integer lvl_dimid, rel_dimid,rec_dimid

c initialise variable ids
      integer t_varid,loc_varid,part_varid,traj_varid,coor_varid
     & ,rel_varid,rec_varid

      
c number of dimensions in space(2: x & y)
      integer spatial_dims
      parameter (spatial_dims = 2)
c number of dimensions in the particle location variables : space, particle number, realisation number, time, layer
      integer ndims
      parameter (ndims = 5)
      integer dimid(ndims)
c number of layers
      integer nlvls
      parameter (nlvls=2)
c dimension of the time array 
      integer tdimid(2)
      

      
c create NETCDF file
      retval = nf_create(file_name,nf_clobber,ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
c Define dimensions in the NETCDF files
      retval = nf_def_dim(ncid,part_name,npoints,part_dimid) ! defines particle numbers
      if (retval .ne. nf_noerr) call handle_err(retval)
      
c      retval = nf_def_dim(ncid,t_name,nf_unlimited,t_dimid) ! time, allowed to go to infinity
c      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_dim(ncid,loc_name,spatial_dims,loc_dimid) ! spatial dimensions for the location
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_dim(ncid,lvl_name,nlvls,lvl_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)  
      
      retval = nf_def_dim(ncid,rel_name,release_no,rel_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval) 
      
      retval = nf_def_dim(ncid,rec_name,nf_unlimited,rec_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
c Define variables
c Time -->
c Need time to be size of number of releases so that it stores 
c time record for each realisation
      
      dimid(1) = loc_dimid
      dimid(2) = part_dimid
      dimid(3) = rel_dimid
      dimid(4) = lvl_dimid
      dimid(5) = rec_dimid
      
      tdimid(1) = rel_dimid
      tdimid(2) = rec_dimid
      
      
c Define Trajectory matrix, which stores locations of all npoints lagrangian particles at time t

      retval = nf_def_var(ncid,traj_name,nf_double,ndims,dimid
     & ,traj_varid)
      if (retval .ne. nf_noerr) then 
      print *, 'traj not created'
      call handle_err(retval)
      endif
      
      retval = nf_def_var(ncid,coor_name,nf_int,ndims,dimid,
     & coor_varid)
      if (retval .ne. nf_noerr) call handle_err(retval) 
      
      retval = nf_put_att_text(ncid,traj_varid,units
     & ,len(x_units),x_units)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_var(ncid,t_name,nf_double,2,tdimid,t_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_put_att_text(ncid,t_varid,units,len(t_units)
     & ,t_units)
      if (retval .ne. nf_noerr) call handle_err(retval)
    
      
      retval = nf_enddef(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      
      retval = nf_close(ncid)                                         !! Close file
      if (retval .ne. nf_noerr) call handle_err(retval)
   
      print *,'!!! Netcdf file created i.e. ', FILE_NAME, '!!!'
      print *,'!!!                                         !!!'

      end subroutine create_release_file  
      
c ----- WRITE TO THE RELEASE FILE -------------------
c -----------------------------------------------------
c ---------------------------------------
      
      subroutine write_release_file(file_name,npoints,release_no,x1,y1
     & ,x2,y2,x1_coord,y1_coord,x2_coord,y2_coord
     & ,release,save_time,step_time)

c release_no is the number of realisations
c release is the realisation you wish to write to
      
      implicit none
      include 'netcdf.inc'
      
      character*(*) file_name
      integer npoints,retval,ncid,i,j, step_time,traj_varid
     & ,x1_coord(npoints),y1_coord(npoints),coor_varid
     & ,x2_coord(npoints),y2_coord(npoints)
      integer release,release_no
      real*8 x1(npoints),y1(npoints), save_time,x2(npoints),y2(npoints)
      
c Temporary trajectory variable
      real*8 traj_temp(2,npoints,2)
      
      integer coor(2,npoints,2)
      
      integer nlvls
      parameter(nlvls = 2)
      
c Define names
      character*(*) traj_name, t_name,coor_name,lvl_name, rel_name
     & ,rec_name      
      parameter(traj_name = 'Trajectories')
      parameter(t_name = 'Time')
      parameter(coor_name = 'Domain Coordinate')
      parameter(lvl_name = 'Layer')
      parameter(rel_name = 'Release Number') 
      parameter(rec_name = 'Record Time')

      
c number of dimensions in space(2: x & y)
      integer spatial_dims, t_varid, lvl_dimid
      parameter (spatial_dims = 2)
c number of dimensions in the particle location variables : location, coordinate, particle number, time
      integer ndims
      parameter (ndims = 5)
      integer dimid(ndims)
      integer start(ndims), count(ndims)
      
      integer tdimid(2)
      integer start1(2),count1(2)
      
c open file in write mode
      retval = nf_open(file_name, nf_write,ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
c request the variable ids of the variables in the file
      retval = nf_inq_varid(ncid,traj_name,traj_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_inq_varid(ncid,t_name,t_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)  
      
      retval = nf_inq_varid(ncid,coor_name,coor_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      

c Write x, y to the temporary trajectory variable
      do j = 1,npoints
              traj_temp(1,j,1) = x1(j)
              traj_temp(2,j,1) = y1(j)
              traj_temp(1,j,2) = x2(j)
              traj_temp(2,j,2) = y2(j)
              coor(1,j,1) = x1_coord(j)
              coor(2,j,1) = y1_coord(j)
              coor(1,j,2) = x2_coord(j)
              coor(2,j,2) = y2_coord(j)
      end do
      
c Now tell the program where to write the data in the NETCDF file
      count(1) = spatial_dims
      count(2) = npoints
      count(3) = 1
      count(4) = nlvls
      count(5) = 1
      start(1) = 1
      start(2) = 1
      start(3) = release ! tells the program to start writing at the release number
      start(4) = 1
      start(5) = step_time ! tells the program to start writing from step time, in the time variable, but writes in all columns in the location and particle numbers

c Update data
      retval = nf_put_vara_double(ncid,traj_varid,start,count,traj_temp)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_int(ncid,coor_varid,start,count,coor)
      if (retval .ne. nf_noerr) call handle_err(retval)
      

c Update time

      count1(1) = 1
      count1(2) = 1
      start1(1) = release
      start1(2) = step_time
      
      retval = nf_put_vara_double(ncid, t_varid, start1, count1,
     +       save_time)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_close(ncid)                                         !! Close NetCDF file                
      if (retval .ne. nf_noerr) call handle_err(retval)
        
      print *,'!!! Netcdf file modified i.e. ', FILE_NAME, '!!!'
      
      
      end subroutine write_release_file   
      
c -----------------------------------------------------------------------
c ----------------------------------------------------------------------

C SUBROUTINE THAT READS THE TIME ARRAY FROM THE TRAJECTORY RELEASE FILES
C RETURNS A TWO DIMENSIONAL MATRIX TIME
      
      subroutine read_release_time(file_name,time,t_len,nrel)
      
      implicit none
      
            include 'netcdf.inc'
            
      character*(*) t_name, loc_name, part_name,traj_name
      character*(*) coor_name,lvl_name,rel_name,rec_name
      parameter (t_name = 'Time',loc_name ='Location of Particle')
      parameter(part_name = 'Particle Number')
      parameter(traj_name = 'Trajectories')
      parameter(coor_name = 'Domain Coordinate')
      parameter(lvl_name = 'Layer')
      parameter(rel_name = 'Release Number')
      parameter(rec_name = 'Record Time')
      

            character*(*) file_name
            integer t_len,nrel,npoints ! nrel is the number of releases
                                 ! t_len is the length of a release
                                 ! npoints is the number of particles in each experiment
            real*8 time(nrel,t_len)
            
        
            integer retval, ncid
            
            integer t_dimid, t_varid,n_dimid,n_varid
            integer rec_dimid,rel_dimid,rec_varid,rel_varid
            
            integer ndims
            parameter(ndims=3)
            !integer count(ndims), start(ndims)
            
              retval = nf_open(FILE_NAME, nf_nowrite, ncid)                   !! Open file in read mode
              if (retval .ne. nf_noerr) call handle_err(retval)
              

              retval = nf_inq_varid(ncid,t_name,t_varid)
              if (retval .ne. nf_noerr) then
                print *, 'time varid not found'
                call handle_err(retval)
                end if


              retval = nf_get_var_double(ncid,t_varid,time)   ! Read time first
              if (retval .ne. nf_noerr) call handle_err(retval)
              
              retval = nf_close(ncid)
              if (retval .ne. nf_noerr) call handle_err(retval)
              
              !deallocate(tp_time,u,v)
              
              write(*,*) 'binned time read'
      
      end subroutine read_release_time
      
C --------------------------------------------------------
C SUBROUTINE THAT READS THE SIZE OF THE TIME ARRAY FROM THE TRAJECTORY RELEASE FILES
C RETURNS A TWO DIMENSIONAL MATRIX TIME WITH THE NUMBER OF PARTICLES AS WELL
      
      subroutine read_release_time_length(file_name,t_len,nrel,npoints)
      

      
      implicit none
      
            include 'netcdf.inc'
            
      character*(*) t_name, loc_name, part_name,traj_name
      character*(*) coor_name,lvl_name,rel_name,rec_name
      parameter (t_name = 'Time',loc_name ='Location of Particle')
      parameter(part_name = 'Particle Number')
      parameter(traj_name = 'Trajectories')
      parameter(coor_name = 'Domain Coordinate')
      parameter(lvl_name = 'Layer')
      parameter(rel_name = 'Release Number')
      parameter(rec_name = 'Record Time')
      
                  character*(*) file_name
            integer t_len,nrel,npoints ! nrel is the number of releases
                                 ! t_len is the length of a release
                                 ! npoints is the number of particles in each experiment
        
            integer retval, ncid
            
            integer t_dimid, t_varid,n_dimid,n_varid
            integer rec_dimid,rel_dimid,rec_varid,rel_varid
            
            integer ndims
            parameter(ndims=3)
            !integer count(ndims), start(ndims)
            
              retval = nf_open(FILE_NAME, nf_nowrite, ncid)                   !! Open file in read mode
              if (retval .ne. nf_noerr) call handle_err(retval)
              

              retval = nf_inq_varid(ncid,t_name,t_varid)
              if (retval .ne. nf_noerr) then
                print *, 'time varid not found'
                call handle_err(retval)
                end if
              
              retval = nf_inq_dimid(ncid,rel_name,rel_dimid)
              if (retval .ne. nf_noerr) then
                print *, 'rel dimid not found'
               call handle_err(retval)
               end if
              retval = nf_inq_dimid(ncid,rec_name,rec_dimid)
              if (retval .ne. nf_noerr) then
                print *, 'rec dimid not found'
               call handle_err(retval)
               end if
              retval = nf_inq_dimid(ncid,part_name,n_dimid)
              if (retval .ne. nf_noerr) then
                print *, 'particle dimid not found'
               call handle_err(retval)
               end if

               
              retval = nf_inq_dimlen(ncid,rec_dimid,t_len)
              if (retval .ne. nf_noerr) call handle_err(retval)
              print*,t_len
              retval = nf_inq_dimlen(ncid,rel_dimid,nrel)
              if (retval .ne. nf_noerr) call handle_err(retval)
              print*,nrel
              retval = nf_inq_dimlen(ncid,n_dimid,npoints)
              if (retval .ne. nf_noerr) call handle_err(retval)
              print*,npoints
              
              retval = nf_close(ncid)                   !! close file
              if (retval .ne. nf_noerr) call handle_err(retval)
              
      
      end subroutine read_release_time_length
      
c --------------------------------------------------------------
c SUBROUTINE TO READ TRAJECTORY DATA IN THE RELEASE FILE
        subroutine read_release_file(file_name,npoints,x1,y1
     &   ,x1_coord,y1_coord 
     & ,rel_no,t_index)
     
c rel_no indicates the release experiment number
c t_index indicates the location in the time array
     
        implicit none
        include 'netcdf.inc'
        
        
        integer,intent(out) :: x1_coord(npoints),y1_coord(npoints)
        integer,intent(in) :: npoints,rel_no,t_index
        character*(*),intent(in) :: file_name
        real*8,intent(out) :: x1(npoints),y1(npoints)
        real*8 read_tim,new_tim
        integer retval,ncid, traj_varid, t_varid, traj_dimid, t_dimid
     & ,coor_varid,coor_dimid,t_len,i
        real *8, dimension (:), allocatable :: tp_time
        
      character*(*) t_name, loc_name, part_name,traj_name
      character*(*) coor_name,lvl_name,rel_name,rec_name
      parameter (t_name = 'Time',loc_name ='Location of Particle')
      parameter(part_name = 'Particle Number')
      parameter(traj_name = 'Trajectories')
      parameter(coor_name = 'Domain Coordinate')
      parameter(lvl_name = 'Layer')
      parameter(rel_name = 'Release Number')
      parameter(rec_name = 'Record Time')
        
c number of dimensions in space(2: x & y)
      integer spatial_dims
      parameter (spatial_dims = 2)
      integer nlvls
      parameter (nlvls = 2)
c number of dimensions in the particle location variables : location, coordinate, particle number, time
      integer ndims
      parameter (ndims = 5)
      integer dimid(ndims)
      integer start(ndims), count(ndims)
      
      !print*,'start reading release file'
      !print*,'file_name',file_name
      
      ncid = 0

        retval = nf_open(file_name,nf_nowrite,ncid) !open file in read mode
        !print*,'retval obtained'
        if (retval .ne. nf_noerr) call handle_err(retval)
        !print*,'file opened'
        
        retval = nf_inq_varid(ncid,traj_name,traj_varid)
        if (retval .ne. nf_noerr) then
        print *, 'trajectories not found'
        call handle_err(retval)
        endif
        !print*,'trajectories found'
        retval = nf_inq_varid(ncid,t_name,t_varid)
        if (retval .ne. nf_noerr) then 
        print *, 'time not found'
        call handle_err(retval)
        endif
        !print*,'time found'
        retval = nf_inq_varid(ncid,coor_name,coor_varid)
        if (retval .ne. nf_noerr) then 
        print *, 'coordinate not found'
        call handle_err(retval)
        endif
        !print*,'cooridinate found'
        

      count(1) = 1
      count(2) = npoints
      count(3) = 1
      count(4) = 1
      count(5) = 1
      start(1) = 1
      start(2) = 1
      start(3) = rel_no
      start(4) = 1
      start(5) = t_index ! tells the program to start writing from step time, in the time variable, but writes in all columns in the location and particle numbers
      
      retval = nf_get_vara_double(ncid,traj_varid,start,count,x1)
      if (retval .ne. nf_noerr) call handle_err(retval)
      !print*,'found traj'
      retval = nf_get_vara_int(ncid,coor_varid,start,count,x1_coord)
      if (retval .ne. nf_noerr) call handle_err(retval)
      !print*,'found coord'
      
      
      start(1) = 2
      retval = nf_get_vara_double(ncid,traj_varid,start,count,y1)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_get_vara_int(ncid,coor_varid,start,count,y1_coord)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      
      
        retval = nf_close(ncid)                                         !! Close file
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        !print*, 'Particles read'
        
        
       end subroutine
       
c -----------------------------------------------------------------------------
c --------------------  CREATE FILE FOR BINNED TRAJECTORIES -----------------
C -------------------------------------------------------------------------------

      subroutine create_binned_file(file_name,nbins,npoints,release_no)
      
c nbins is the number of bins
c npoints is the number of particles per bin
c nrel is the number of releases
c bin_boundaries contains the upper Y coordinate of the bin boundaries
      
      implicit none
      include 'netcdf.inc'
      
      
      character*(*), intent(in) :: file_name
      integer, intent(in) :: nbins,npoints,release_no
      real*8 bin_boundaries(nbins)
      
      integer ncid,retval,j
c Dimension names
      
      character*(*) t_name,part_name,bin_name,rel_name,layer_name
     & ,loc_name,rec_name,restart_name,t_units,x_units,y_units,units
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
      parameter(restart_name = 'Restart Info')
      
      
c Variable names
      
      character*(*) bin,traj_name,coord_name,width_name,start_name
      parameter(bin = 'Bin Boundaries')
      parameter(traj_name = 'Trajectories')
      parameter(coord_name = 'Domain Coordinate')
      parameter(width_name = 'Number of Particles')
      parameter(start_name = 'Restart Point')

c Dimension sizes

      integer nlvls,spatial_dims,ndims
      parameter(nlvls = 2,spatial_dims = 2,ndims = 6)
      integer dimid(ndims), tdimid(2),bdimid(2)
      
c Dimension ids

      integer bin_dimid,rel_dimid,t_dimid,part_dimid,layer_dimid
     & ,loc_dimid,rec_dimid,start_dimid
      
c Variable ids

      integer bin_varid,traj_varid,coord_varid,t_varid,width_varid
     & ,start_varid
      
c create NETCDF file
      
      retval = nf_create(file_name,nf_clobber,ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
c Define dimensions in the NETCDF files
      retval = nf_def_dim(ncid,part_name,npoints,part_dimid) ! defines particle numbers
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_dim(ncid,loc_name,spatial_dims,loc_dimid) ! spatial dimensions for the location
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_dim(ncid,layer_name,nlvls,layer_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)  
      
      retval = nf_def_dim(ncid,rel_name,release_no,rel_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval) 
      
      retval = nf_def_dim(ncid,rec_name,nf_unlimited,rec_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_dim(ncid,bin_name,nbins,bin_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_dim(ncid,restart_name,2,start_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      
c Define variables
c Time -->
c Need time to be size of number of releases so that it stores 
c time record for each realisation
      
      dimid(1) = bin_dimid
      dimid(2) = loc_dimid
      dimid(3) = part_dimid
      dimid(4) = rel_dimid
      dimid(5) = layer_dimid
      dimid(6) = rec_dimid
      
      tdimid(1) = rel_dimid
      tdimid(2) = rec_dimid
      
      bdimid(1) = bin_dimid
      bdimid(2) = rel_dimid
      

      
      
c Define Trajectory matrix, which stores locations of all npoints lagrangian particles at time t

c Trajectory matrix

      retval = nf_def_var(ncid,traj_name,nf_double,ndims,dimid
     & ,traj_varid)
      if (retval .ne. nf_noerr) then 
      print *, 'traj not created'
      call handle_err(retval)
      endif
      
c Coordinate Matrix
      
      retval = nf_def_var(ncid,coord_name,nf_int,ndims,dimid,
     & coord_varid)
      if (retval .ne. nf_noerr) call handle_err(retval) 
      
      retval = nf_put_att_text(ncid,traj_varid,units
     & ,len(x_units),x_units)
      if (retval .ne. nf_noerr) call handle_err(retval)

c Time Matrix
      
      retval = nf_def_var(ncid,t_name,nf_double,2,tdimid,t_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_put_att_text(ncid,t_varid,units,len(t_units)
     & ,t_units)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
c Record of Bins

      retval = nf_def_var(ncid,bin,nf_double,2,bdimid,bin_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_put_att_text(ncid,bin_varid,units
     & ,len(x_units),x_units)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_var(ncid,width_name,nf_int,2,bdimid,width_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
c Record of stopping time

      retval = nf_def_var(ncid,start_name,nf_int,2,start_dimid
     & ,start_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_enddef(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      
      retval = nf_close(ncid)                                         !! Close file
      if (retval .ne. nf_noerr) call handle_err(retval)
   
      print *,'!!! Netcdf file created i.e. ', FILE_NAME, '!!!'
      print *,'!!!                                         !!!'
  
      end subroutine create_binned_file
      
c --------------------------------------------------------------------------------
c ---------------------- WRITE TO BINNED TRAJECTORY FILE -----------------------
C ---------------------------------------------------------------------------------

      subroutine write_binned_file(file_name,npoints,bin_no,release
     & ,x1,y1
     & ,x2,y2,x1_coord,y1_coord,x2_coord,y2_coord
     & ,save_time,step_time)
     
      implicit none
      include 'netcdf.inc'
      
c write to a single bin for a single release for npoints particles at time save_time

      character*(*), intent(in) :: file_name
      integer,intent(in) :: npoints,bin_no,release,step_time
      real*8,intent(in) :: x1(npoints),y1(npoints)
     & ,x2(npoints),y2(npoints)
      integer, intent(in) :: x1_coord(npoints),y1_coord(npoints)
     & ,x2_coord(npoints),y2_coord(npoints)
      real*8, intent(in) :: save_time
      
      integer ncid,retval,j
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
      
      real*8 traj_temp(2,npoints,2)
      integer coor(2,npoints,2),start_info(2)
      
c Dimension sizes

      integer nlvls,spatial_dims,ndims
      parameter(nlvls = 2,spatial_dims = 2,ndims = 6)
      integer dimid(ndims), tdimid(2)
      
c Dimension ids

      integer bin_dimid,rel_dimid,t_dimid,part_dimid,layer_dimid
     & ,loc_dimid
      
c Variable ids

      integer bin_varid,traj_varid,coord_varid,t_varid,start_varid
      
c Write to file 

      integer count(ndims), start(ndims),count1(2),start1(2)
      
c open file in write mode
      retval = nf_open(file_name, nf_write,ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
c request the variable ids of the variables in the file
      retval = nf_inq_varid(ncid,traj_name,traj_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_inq_varid(ncid,t_name,t_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)  
      
      retval = nf_inq_varid(ncid,coord_name,coord_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_inq_varid(ncid,start_name,start_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      

c Write x, y to the temporary trajectory variable
      do j = 1,npoints
              traj_temp(1,j,1) = x1(j)
              traj_temp(2,j,1) = y1(j)
              traj_temp(1,j,2) = x2(j)
              traj_temp(2,j,2) = y2(j)
              coor(1,j,1) = x1_coord(j)
              coor(2,j,1) = y1_coord(j)
              coor(1,j,2) = x2_coord(j)
              coor(2,j,2) = y2_coord(j)
      end do
      
      start_info(1) = release
      start_info(2) = step_time
      
c Now tell the program where to write the data in the NETCDF file
      count(1) = 1
      count(2) = spatial_dims
      count(3) = npoints
      count(4) = 1
      count(5) = nlvls
      count(6) = 1
      start(1) = bin_no
      start(2) = 1
      start(3) = 1
      start(4) = release ! tells the program to start writing at the release number
      start(5) = 1
      start(6) = step_time ! tells the program to start writing from step time, in the time variable, but writes in all columns in the location and particle numbers

c Update data
      retval = nf_put_vara_double(ncid,traj_varid,start,count,traj_temp)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_int(ncid,coord_varid,start,count,coor)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_put_vara_double(ncid,start_varid,1,1,start_info)
      if (retval .ne. nf_noerr) call handle_err(retval)
      

c Update time

      count1(1) = 1
      count1(2) = 1
      start1(1) = release
      start1(2) = step_time
      
      retval = nf_put_vara_double(ncid, t_varid, start1, count1,
     +       save_time)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_close(ncid)                                         !! Close NetCDF file                
      if (retval .ne. nf_noerr) call handle_err(retval)
        
      !print *,'!!! Netcdf file modified i.e. ', FILE_NAME, '!!!'
      
      end subroutine write_binned_file
      
      
c -------------------------------------------------------------------------
c --------------  SUBROUTINE THAT WRITES BIN WIDTHS -------------------
C -----------------------------------------------------------------------

      subroutine write_width(file_name,bin_count,nbins,k)
      
      implicit none
      include 'netcdf.inc'
      
      character*(*) file_name
      integer nbins
      integer bin_count(nbins)
      
      integer ncid,retval,k,start(2),count(2)
      
      character*(*) bin_name,width_name
      parameter(bin_name = 'Bin Number')
      parameter(width_name = 'Number of Particles')
      
      integer bin_dimid,width_varid
      
      retval = nf_open(file_name, nf_write,ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
c request the variable ids of the variables in the file
      retval = nf_inq_varid(ncid,width_name,width_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      start(1) = 1
      start(2) = k
      
      count(1) = nbins
      count(2) = 1
      
      retval = nf_put_vara_int(ncid,width_varid,start,count,bin_count)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_close(ncid)                                         !! Close NetCDF file                
      if (retval .ne. nf_noerr) call handle_err(retval)
    
      print*,'width written'
      
      end subroutine write_width
      
c ------------------------------------------------------------------------------
c -------------- SUBROUTINE THAT WRITE THE BIN BOUNDARIES TO FILE --------------
C -------------------------------------------------------------------------------

      subroutine write_bin_boundaries(file_name,rel_no,nbins
     & ,bin_boundaries)
     
      implicit none
      include 'netcdf.inc'
      
      character*(*) file_name
      integer rel_no,nbins
      real*8 bin_boundaries(nbins)
      
      character*(*) bin_name,bin,rel_name
      parameter(bin = 'Bin Boundaries')
      parameter(bin_name = 'Bin Number')
      parameter(rel_name = 'Release Number')
      
      integer ncid,retval
      
      integer bin_dimid,rel_dimid,bin_varid
      integer start(2),count(2)
      
c open file in write mode
      retval = nf_open(file_name, nf_write,ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
c request the variable ids of the variables in the file
      retval = nf_inq_varid(ncid,bin,bin_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
c write bin boundaries to file

      start(1) = 1
      start(2) = rel_no
      
      count(1) = nbins
      count(2) = 1

       retval = nf_put_vara_double(ncid,bin_varid
     &  ,start,count,bin_boundaries)
      
      
       retval = nf_close(ncid)                                         !! Close file
        if (retval .ne. nf_noerr) call handle_err(retval)
      
      
     
       end subroutine write_bin_boundaries
      
c ----------------------------------------------------------------------------
c ------------- SUBROUTINE THAT READS THE BINNED TRAJECTORY FILE -------------
c ----------------------------------------------------------------------------

      subroutine read_binned_file(file_name,npoints,nbins
     & ,x1,y1,x2,y2
     & ,x1_coord,y1_coord
     & ,x2_coord,y2_coord,start_info)

c reads the bin_no bin
c and the rel_no release
c and the time index, t_index
c returns the trajectory
c and domain coordinate
      
      implicit none
      include 'netcdf.inc'
      
      
      character*(*) file_name
      integer start_info(2)
      
      integer ncid,retval,npoints,nbins
      
      real*8 x1(nbins,npoints),y1(nbins,npoints)
      integer x1_coord(nbins,npoints),y1_coord(nbins,npoints)
      real*8 x2(nbins,npoints),y2(nbins,npoints)
      integer x2_coord(nbins,npoints),y2_coord(nbins,npoints)
      
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
        
        retval = nf_get_vara_int(ncid,start_varid,1,2,start_info)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
      count(1) = nbins
      count(2) = 1
      count(3) = npoints
      count(4) = 1
      count(5) = 1
      count(6) = 1
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = start_info(1)
      start(5) = 1
      start(6) = start_info(2)! tells the program to start writing from step time, in the time variable, but writes in all columns in the location and particle numbers
      
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
      
      
      
      end subroutine read_binned_file
      
c ----------------------------------------------------------------------------
c ----------- SUBROUTINE THAT READS THE TIME FROM THE BINNED TRAJECTORY FILE -------
c ----------- ALSO RETURNS THE VALUES OF THE BIN BOUNDARIES ---------------------
C --------------------------------------------------------------------------------

      subroutine read_binned_file_time(file_name,nrel,nbins
     & ,t_len,time)
      
      implicit none
      include 'netcdf.inc'
      
      character*(*) file_name
      integer nrel,t_len,nbins
      integer ncid,retval
      
      real*8 time(nrel,t_len)
      
      character*(*) t_name,rel_name
      parameter(t_name = 'Time')
      parameter(rel_name = 'Release Number')

      
      integer t_dimid,rel_dimid,t_varid,bin_varid
      integer start(2),count(2)
      
      
      retval = nf_open(FILE_NAME, nf_nowrite, ncid)                   !! Open file in read mode
      if (retval .ne. nf_noerr) call handle_err(retval)
              

      retval = nf_inq_varid(ncid,t_name,t_varid)
      if (retval .ne. nf_noerr) then
      print *, 'time varid not found'
      call handle_err(retval)
      end if

      

      start(1) = 1
      start(2) = 1
      count(1) = nrel
      count(2) = t_len
      retval = nf_get_vara_double(ncid,t_varid,start,count,time)   ! Read time first
      if (retval .ne. nf_noerr) call handle_err(retval)
      print*,'traj time read'
      

              
      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
  
      end subroutine read_binned_file_time
      
c -----------------------------------------------------------------------------
c ------------- SUBROUTINE THAT READS THE DIMENSIONS OF THE TRAJECTORY DATA ---
C ----------------------------------------------------------------------------
      
      subroutine read_binned_file_dimensions(file_name
     & ,nbins,nrel,npoints,t_len)
      
      implicit none
      include 'netcdf.inc'
      
      character*(*) file_name
      integer nbins,nrel,npoints,t_len
      
      integer ncid,retval
      
c Dimension names
      
      character*(*) t_name,part_name,bin_name,rel_name,layer_name
     & ,loc_name,rec_name
      parameter(t_name = 'Time')
      parameter(part_name = 'Particle Number')
      parameter(bin_name = 'Bin Number')
      parameter(rel_name = 'Release Number')
      parameter(layer_name = 'Layer')
      parameter(loc_name = 'Location of Particle')
      parameter(rec_name = 'Record Time')
      
c Dimension IDs

      integer bin_dimid,rel_dimid,t_dimid,part_dimid,layer_dimid
     & ,loc_dimid,rec_dimid
      
c Variable ids

      integer traj_varid,t_varid
      
      integer ndims
      parameter(ndims = 6)
      
      
              retval = nf_open(FILE_NAME, nf_nowrite, ncid)                   !! Open file in read mode
              if (retval .ne. nf_noerr) call handle_err(retval)
              

              retval = nf_inq_varid(ncid,t_name,t_varid)
              if (retval .ne. nf_noerr) then
                print *, 'time varid not found'
                call handle_err(retval)
                end if
              
              retval = nf_inq_dimid(ncid,rel_name,rel_dimid)
              if (retval .ne. nf_noerr) then
                print *, 'rel dimid not found'
               call handle_err(retval)
               end if
              retval = nf_inq_dimid(ncid,rec_name,rec_dimid)
              if (retval .ne. nf_noerr) then
                print *, 'rec dimid not found'
               call handle_err(retval)
               end if
              retval = nf_inq_dimid(ncid,part_name,part_dimid)
              if (retval .ne. nf_noerr) then
                print *, 'particle dimid not found'
               call handle_err(retval)
               end if
              retval = nf_inq_dimid(ncid,bin_name,bin_dimid)
              if (retval .ne. nf_noerr) then
                print *, 'bin dimid not found'
               call handle_err(retval)
               end if

               
              retval = nf_inq_dimlen(ncid,rec_dimid,t_len)
              if (retval .ne. nf_noerr) call handle_err(retval)
              retval = nf_inq_dimlen(ncid,rel_dimid,nrel)
              if (retval .ne. nf_noerr) call handle_err(retval)
              retval = nf_inq_dimlen(ncid,part_dimid,npoints)
              if (retval .ne. nf_noerr) call handle_err(retval)
              retval = nf_inq_dimlen(ncid,bin_dimid,nbins)
              if (retval .ne. nf_noerr) call handle_err(retval)
              
              retval = nf_close(ncid)                   !! close file
              if (retval .ne. nf_noerr) call handle_err(retval)

      end subroutine read_binned_file_dimensions
      
      
      subroutine read_bin_width(file_name,nbins,bin_width)
      
      implicit none
      
      include 'netcdf.inc'
      
      character*(*) file_name
      integer nbins,ncid,retval
      real*8 bin_width(nbins)
      
      character*(*) bin_name
      parameter(bin_name = 'Bin Width')
      
      integer bin_varid, bin_dimid
      
      retval = nf_open(FILE_NAME, nf_nowrite, ncid)                   !! Open file in read mode
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_inq_varid(ncid,bin_name,bin_varid)
      if (retval .ne. nf_noerr) then
      print *, 'bin varid not found'
      call handle_err(retval)
      end if
      
      retval = nf_get_vara_double(ncid,bin_varid,1,nbins,bin_width)   ! Read time first
      if (retval .ne. nf_noerr) call handle_err(retval)
      print*,'bin width read'
      
      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      
      end subroutine read_bin_width
      
       subroutine read_release_bin_width(file_name,nbins,nrel,bin_width)
      
      implicit none
      
      include 'netcdf.inc'
      
      character*(*) file_name
      integer nbins,ncid,retval,nrel
      real*8 bin_width(nbins,nrel)
      integer start(2),count(2)
      
      character*(*) bin_name
      parameter(bin_name = 'Release Bin Width')
      
      integer bin_varid, bin_dimid
      
      retval = nf_open(FILE_NAME, nf_nowrite, ncid)                   !! Open file in read mode
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_inq_varid(ncid,bin_name,bin_varid)
      if (retval .ne. nf_noerr) then
      print *, 'bin varid not found'
      call handle_err(retval)
      end if
      
      start(1) = 1
      start(2) = 1
      count(1) = nbins
      count(2) = nrel
      
      retval = nf_get_vara_double(ncid,bin_varid,start,count,bin_width)   ! Read time first
      if (retval .ne. nf_noerr) call handle_err(retval)
      print*,'bin width read'
      
      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      
      end subroutine read_release_bin_width
    
      end module MOD_traj_NETCDF
