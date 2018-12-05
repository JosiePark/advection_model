c MODULE THAT READS THE DIFFUSIVITY TENSOR

      module mod_diffusion_netcdf
      
      use mod_netcdf_error
      
      contains
      
      subroutine read_diffusivity(diff_file,K)
      
      implicit none
      include 'netcdf.inc'
      
      character*(*) diff_file
      real*8,allocatable,dimension(:,:,:) :: K
      integer ncid,retval,i
      integer diff_varid
      
      character*(*) diff_name
      parameter(diff_name = 'Diffusivity')
      
      integer count(3),start(3),dimid(3),dim_len(3)
      
      retval = nf_open(diff_file, nf_nowrite, ncid)                   !! Open file in read mode
      if (retval .ne. nf_noerr) call handle_err(retval)
      print*,'diff_file found'
      
      
      retval = nf_inq_varid(ncid,diff_name,diff_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      print*,'diffusivity found'
      
      
      
      retval = nf_inq_vardimid(ncid,diff_varid,dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      do i = 1,3
        
        retval = nf_inq_dimlen(ncid,dimid(i),dim_len(i))
        if (retval .ne. nf_noerr) call handle_err(retval)
        
      enddo
      
      print*,'dimension lengths = ',dim_len
      
      do i = 1,3
      
        count(i) = dim_len(i)
        start(i) = 1
      
      enddo
      
      allocate(K(dim_len(1),dim_len(2),dim_len(3)))
      
      print*,shape(K)
      
      
      retval = nf_get_vara_double(ncid,diff_varid,start,count,K)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      end subroutine read_diffusivity
      
c --------------------------------------------------------------------------
c ----------- READ PV MAPPED DIFFUSIVITY ------------------------------------
      
      subroutine read_PV_diffusivity(diff_file,K)
      
      implicit none
      include 'netcdf.inc'
      
      character*(*) diff_file
      real*8,allocatable,dimension(:) :: K
      integer ncid,retval,i
      integer diff_varid
      
      character*(*) diff_name
      parameter(diff_name = 'Diffusivity')
      
      integer count,start,dimid,dim_len
      
      retval = nf_open(diff_file, nf_nowrite, ncid)                   !! Open file in read mode
      if (retval .ne. nf_noerr) call handle_err(retval)
      print*,'diff_file found'
      
      
      retval = nf_inq_varid(ncid,diff_name,diff_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      print*,'diffusivity found'
      
      
      
      retval = nf_inq_vardimid(ncid,diff_varid,dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_inq_dimlen(ncid,dimid,dim_len)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      allocate(K(dim_len))

      
      print*,'dimension lengths = ',dim_len
      
      start = 1
      count = dim_len
      
      retval = nf_get_vara_double(ncid,diff_varid,start,count,K)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      end subroutine read_PV_diffusivity
      
      
c ---------------------------------------------------------------------------
c ---------------------------------------------------------------------------
      
      subroutine create_diffusion_trajfile(file_name,npoints,nbins)
      
      implicit none
      include 'netcdf.inc'
      
      character*(*) file_name
      integer npoints,nbins
      
      character*(*) traj_name,t_name,dim_name,part_name,bin_name
     & ,lvl_name
      parameter(traj_name = 'Trajectories')
      parameter(part_name = 'Particle Number')
      parameter(dim_name = 'Dimension')
      parameter(bin_name = 'Bin')
      parameter(lvl_name = 'Layer')
      parameter(t_name = 'Time')
      
      integer traj_varid,t_varid,lvl_dimid,dim_dimid,part_dimid
     & ,bin_dimid,t_dimid
     
      integer ncid,retval,ndims
      parameter(ndims = 5)
      
      integer dimid(ndims)
     
      retval = nf_create(file_name,nf_clobber,ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_dim(ncid,part_name,npoints,part_dimid) ! defines particle numbers
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_dim(ncid,t_name,nf_unlimited,t_dimid) ! time, allowed to go to infinity
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_dim(ncid,dim_name,2,dim_dimid) ! spatial dimensions for the location
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_dim(ncid,lvl_name,2,lvl_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)  
      
      retval = nf_def_dim(ncid,bin_name,nbins,bin_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      
      
      dimid(1) = dim_dimid
      dimid(2) = part_dimid
      dimid(3) = bin_dimid
      dimid(4) = lvl_dimid
      dimid(5) = t_dimid
      
      retval = nf_def_var(ncid,traj_name,nf_double,ndims,dimid,
     & traj_varid)
      if (retval .ne. nf_noerr) call handle_err(retval) 
      
      retval = nf_def_var(ncid,t_name,nf_double,1,t_dimid,t_varid)
      if (retval .ne. nf_noerr) call handle_err(retval) 
      
      retval = nf_enddef(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      
      retval = nf_close(ncid)                                         !! Close file
      if (retval .ne. nf_noerr) call handle_err(retval)
   
      print *,'!!! Netcdf file created i.e. ', FILE_NAME, '!!!'
        

      end subroutine create_diffusion_trajfile
      
c -------------------------------------------------------------------------------
c ------------------------------------------------------------------------------

      subroutine create_pvdiffusion_trajfile(file_name,npoints,nbins)
      
      implicit none
      include 'netcdf.inc'
      
      character*(*) file_name
      integer npoints,nbins
      
      character*(*) traj_name,t_name,dim_name,part_name,bin_name
     & ,lvl_name
      parameter(traj_name = 'Trajectories')
      parameter(part_name = 'Particle Number')
      parameter(dim_name = 'Dimension')
      parameter(bin_name = 'Bin')
      parameter(lvl_name = 'Layer')
      parameter(t_name = 'Time')
      
      integer traj_varid,t_varid,lvl_dimid,dim_dimid,part_dimid
     & ,bin_dimid,t_dimid
     
      integer ncid,retval,ndims
      parameter(ndims = 3)
      
      integer dimid(ndims)
     
      retval = nf_create(file_name,nf_clobber,ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_dim(ncid,part_name,npoints,part_dimid) ! defines particle numbers
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_dim(ncid,t_name,nf_unlimited,t_dimid) ! time, allowed to go to infinity
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_dim(ncid,bin_name,nbins,bin_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      
      dimid(1) = part_dimid
      dimid(2) = bin_dimid
      dimid(3) = t_dimid
      
      retval = nf_def_var(ncid,traj_name,nf_double,ndims,dimid,
     & traj_varid)
      if (retval .ne. nf_noerr) call handle_err(retval) 
      
      retval = nf_def_var(ncid,t_name,nf_double,1,t_dimid,t_varid)
      if (retval .ne. nf_noerr) call handle_err(retval) 
      
      retval = nf_enddef(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      
      retval = nf_close(ncid)                                         !! Close file
      if (retval .ne. nf_noerr) call handle_err(retval)
   
      print *,'!!! Netcdf file created i.e. ', FILE_NAME, '!!!'
        

      end subroutine create_pvdiffusion_trajfile
      
c ----------------------------------------------------------------------------
c ----------------------------------------------------------------------------
      
      subroutine write_diffusion_trajfile(file_name,npoints,nbins
     & ,x1,y1,x2,y2,time,nrec)
      
      implicit none
      include 'netcdf.inc'
      
      character*(*) file_name
      integer npoints,nbins,nrec,p,b
      
      real*8 x1(npoints,nbins),x2(npoints,nbins)
     & ,y1(npoints,nbins),y2(npoints,nbins),time
      real*8 traj_temp(2,npoints,nbins,2)
     
      character*(*) traj_name,t_name
      parameter(traj_name = 'Trajectories')
      parameter(t_name = 'Time')
      
      integer traj_varid,t_varid
      integer ncid,retval,ndims
      parameter(ndims = 5)
      integer count(ndims),start(ndims)
      
      retval = nf_open(file_name,nf_write,ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_inq_varid(ncid,traj_name,traj_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_inq_varid(ncid,t_name,t_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      do p = 1,npoints
      
        do b = 1,nbins
        
            traj_temp(1,p,b,1) = x1(p,b)
            traj_temp(2,p,b,1) = y1(p,b)
            traj_temp(1,p,b,2) = x2(p,b)
            traj_temp(2,p,b,2) = y2(p,b)
            
        enddo
      
      enddo
      
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
      start(5) = nrec
      count(1) = 2
      count(2) = npoints
      count(3) = nbins
      count(4) = 2
      count(5) = 1
      
      retval = nf_put_vara_double(ncid,traj_varid,start,count,traj_temp)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_put_vara_double(ncid,t_varid,nrec,1,time)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      
      retval = nf_close(ncid)                                         !! Close NetCDF file                
      if (retval .ne. nf_noerr) call handle_err(retval)
        
      print *,'!!! Netcdf file modified i.e. ', FILE_NAME, '!!!'
      
      end subroutine write_diffusion_trajfile
      
c ---------------------------------------------------------------------
c ---------------------------------------------------------------------

      subroutine write_pvdiffusion_trajfile(file_name,npoints,nbins
     & ,y,time,nrec)
      
      implicit none
      include 'netcdf.inc'
      
      character*(*) file_name
      integer npoints,nbins,nrec,p,b
      
      real*8 y(npoints,nbins),time
      real*8 traj_temp(npoints,nbins)
     
      character*(*) traj_name,t_name
      parameter(traj_name = 'Trajectories')
      parameter(t_name = 'Time')
      
      integer traj_varid,t_varid
      integer ncid,retval,ndims
      parameter(ndims = 3)
      integer count(ndims),start(ndims)
      
      retval = nf_open(file_name,nf_write,ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_inq_varid(ncid,traj_name,traj_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_inq_varid(ncid,t_name,t_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      do p = 1,npoints
      
        do b = 1,nbins
        
            traj_temp(p,b) = y(p,b)
            
        enddo
      
      enddo
      
      start(1) = 1
      start(2) = 1
      start(3) = nrec
      count(1) = npoints
      count(2) = nbins
      count(3) = 1
      
      retval = nf_put_vara_double(ncid,traj_varid,start,count,traj_temp)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_put_vara_double(ncid,t_varid,nrec,1,time)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      
      retval = nf_close(ncid)                                         !! Close NetCDF file                
      if (retval .ne. nf_noerr) call handle_err(retval)
        
      print *,'!!! Netcdf file modified i.e. ', FILE_NAME, '!!!'
      
      end subroutine write_pvdiffusion_trajfile
      
      end module mod_diffusion_netcdf
