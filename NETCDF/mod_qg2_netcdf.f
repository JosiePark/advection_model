      module MOD_qg2_NETCDF
      
      use mod_netcdf_error

      contains
      
c ----------------------- MODULE THAT CREATES NETCDF FILE THAT STORES CABARET STREAM FUNCTION 
C ----------------------- AND ENERGIES AND TIME ---------------------------------------------

        subroutine create_netcdf_file(FILE_NAME,Num_x,Num_y,BasinScale)
        implicit none
        include 'netcdf.inc'
            
        integer Num_x, Num_y                                            !! Declare Variables
        real*8 BasinScale 
        character*(*) FILE_NAME
        integer NDIMS, NLVLS, ETYPE
        integer ncid, x_dimid, y_dimid, rec_dimid, lvl_dimid, Et_dimid
        parameter (NDIMS = 4, NLVLS=2, ETYPE=4)  
        character*(*) X_NAME, Y_NAME, REC_NAME, Lvl_NAME, Et_NAME
        parameter (X_NAME='X', Y_NAME='Y', REC_NAME='Time', 
     &              Et_NAME='Energies') 
        parameter (Lvl_NAME = 'Layer')
        integer start(NDIMS), count(NDIMS) 
            
        real*8 X_val(Num_x), Y_val(Num_y), dx, dy                       !! Create variables to hold grid data
        integer x_varid, y_varid, t_varid
        
        character*(*) Psi_name,Ene_NAME 
             !! Assign names to the variaables which  
                       !! will be used in the netcdf file  
        parameter(Psi_NAME='Stream Function')                           !! 'X', 'Y', 'Time' will be use for grid and
                      !! time storage   
        parameter(Ene_NAME='Energy') 

        integer psi_varid, ene_varid, 
     &          dimids(NDIMS), dimid(2)
        
        character*(*) UNITS, X_UNITS, Y_UNITS, Time_UNITS,    !! Assign Units to variables
     &                Psi_UNITS, Ene_UNITS, Ene_INFO
        parameter (UNITS='units')

        parameter(Psi_UNITS = 'non-dimensionalized by U=1 cm/s, L=dx') 

        parameter(Ene_UNITS = 'non-dimensionalized by U=1 cm/s, L=dx')
        parameter(Ene_INFO = 'Potential, Kinetic L1, Kinetic L2, Total') 
        parameter (X_UNITS='cms', Y_UNITS='cms', Time_UNITS='Days')

        
        real, parameter:: START_X=0., START_Y=0.
        integer lvl, x, y, retval                                       !! Loop indics
                                                 
                    
        dx = BasinScale/dfloat(Num_x)
        dy = dx
        do x = 1, Num_x                                                 !! Values for grid points 
            X_val(x) = START_X + dfloat((x-1))*dx
        end do
        do y = 1, Num_y
            Y_val(y) = START_Y + dfloat((y-1))*dy
        end do
        
        retval = nf_create(FILE_NAME, nf_clobber, ncid)                 !! Create netcdf file
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = nf_def_dim(ncid, Lvl_NAME, NLVLS, lvl_dimid)           !! Define Dimensions in NetCDF file
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_dim(ncid, Et_NAME, ETYPE, Et_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_dim(ncid, X_NAME, Num_x, x_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_dim(ncid, Y_NAME, Num_y, y_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_dim(ncid, REC_NAME, NF_UNLIMITED, rec_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = nf_def_var(ncid,X_NAME,NF_DOUBLE,1,x_dimid,x_varid)    !! Create variables in NetCDF file
        if (retval .ne. nf_noerr) call handle_err(retval)               !! '+' or '&' in column 6 always
        retval = nf_def_var(ncid,Y_NAME,NF_DOUBLE,1,y_dimid,y_varid) 
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid,REC_NAME,NF_DOUBLE,1,rec_dimid,t_varid) 
        if (retval .ne. nf_noerr) call handle_err(retval)
            
        retval = nf_put_att_text(ncid, x_varid, UNITS, len(X_UNITS),    !! Assign units
     +        X_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, y_varid, UNITS, len(Y_UNITS), 
     +         Y_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, t_varid, UNITS, len(Time_UNITS), 
     +         Time_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
            
        dimids(1) = x_dimid
        dimids(2) = y_dimid
        dimids(3) = lvl_dimid
        dimids(4) = rec_dimid
        
        retval = nf_def_var(ncid, Psi_NAME, NF_DOUBLE, NDIMS, dimids,   !! Creating Main varialbes
     +        psi_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        
        retval = nf_put_att_text(ncid, psi_varid, UNITS, 
     +       len(Psi_UNITS), Psi_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)

        
        dimid(1) = Et_dimid
        dimid(2) = rec_dimid
        
        retval = nf_def_var(ncid, Ene_NAME, NF_DOUBLE, 2, dimid,        !! Energy Variable which has 2 dimensions
     +        ene_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, ene_varid, UNITS, 
     +       len(Ene_UNITS), Ene_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, ene_varid, 'Energies column 1-4', 
     +       len(Ene_INFO), Ene_INFO)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = nf_enddef(ncid)                                        !! End define mode    
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = nf_put_var_double(ncid, x_varid, X_val)                !! Write data in X, Y
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_var_double(ncid, y_varid, Y_val)
        if (retval .ne. nf_noerr) call handle_err(retval)
            
        retval = nf_close(ncid)                                         !! Close file
        if (retval .ne. nf_noerr) call handle_err(retval)
   
        print *,'!!! Netcdf file created i.e. ', FILE_NAME, '!!!'
        print *,'!!!                                         !!!'

       end subroutine create_netcdf_file
    
c    -------------------------------------------------------------------------
c   READ DATA FROM CABARET NETCDF FILE ------------------------------------
C -------------------------------------------------------------------------

       subroutine read_netcdf(FILE_NAME,psi1,psi2,Num_x,Num_y,read_tim,
     +  new_tim, step_tim)
       
        implicit none
        include 'netcdf.inc'
        
        integer Num_x, Num_y, step_tim
        character*(*) FILE_NAME
        real*8 psi1(Num_x,Num_y), psi2(Num_x,Num_y), read_tim, new_tim
        
        integer NDIMS, NLVLS, ETYPE, i, j
        parameter (NDIMS = 4, ETYPE = 4, NLVLS=2)
        integer ncid, retval, dimids(NDIMS), dimid(2),
     +   x_dimid, y_dimid, rec_dimid, lvl_dimid, Et_dimid, 
     +   start(NDIMS), count(NDIMS), start1(2), count1(2)
     
        character*(*) Psi_NAME, REC_NAME                                       
        parameter(Psi_NAME='Stream Function',REC_NAME='Time')
        integer psi_varid, t_varid, t_len
        real *8 tp_psi(Num_x,Num_y,NLVLS)
        real *8, dimension (:), allocatable :: tp_time
        
        retval = nf_open(FILE_NAME, nf_nowrite, ncid)                   !! Open file in read mode
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = nf_inq_varid(ncid, Psi_NAME, psi_varid)                !! Request varid
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_inq_varid(ncid, REC_NAME, t_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = nf_inq_dimid(ncid, REC_NAME, rec_dimid)                !! Request dimension length
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_inq_dimlen(ncid, rec_dimid, t_len)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        allocate(tp_time(t_len))
        retval = nf_get_vara_double(ncid, t_varid, 1, t_len, tp_time)   ! Read time first
        if (retval .ne. nf_noerr) call handle_err(retval)
        
c        print *, t_len, tp_time(t_len)
        
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
                exit
              end if
            end do
          end if
        else
          new_tim = tp_time(1)
          step_tim = 1
          print *,'Starting simulation from t = ',new_tim,' Days'
        end if
        
        count(1) = Num_x                                                                     
        count(2) = Num_y
        count(3) = 1                                                    !! Start, count tell the program where to  
        count(4) = 1                                                    !! write data in NetCDF file 
        start(1) = 1
        start(2) = 1
        start(3) = 1
        start(4) = step_tim
        
        retval = nf_get_vara_double(ncid, psi_varid, start, count,      !! Read Psi
     +        psi1)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        start(3) = 2
        retval = nf_get_vara_double(ncid, psi_varid, start, count,
     +        psi2)
        if (retval .ne. nf_noerr) call handle_err(retval)        
     
        retval = nf_close(ncid)                                         !! Close file
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        deallocate (tp_time)
        
       end subroutine
       
c    -----------------------------------------------------------------
c   WRITE STREAM FUNCTION AND ENERGY TO CABARET NETCDF FILE ---------
C ---------------------------------------------------------------------

       subroutine write_netcdf(FILE_NAME,psi1,psi2,
     +    epot,ekin1,ekin2,Num_x,Num_y,save_tim,step_tim)
     
        implicit none
        include 'netcdf.inc'
        
        integer Num_x, Num_y, step_tim 
        character*(*) FILE_NAME
        real*8 psi1(Num_x,Num_y), psi2(Num_x,Num_y),
     +     epot, ekin1, ekin2, save_tim
     
        integer NDIMS, NLVLS, ETYPE, i, j
        parameter (NDIMS = 4, ETYPE = 4, NLVLS=2)
        integer ncid, retval, dimids(NDIMS), dimid(2),
     +   x_dimid, y_dimid, rec_dimid, lvl_dimid, Et_dimid, 
     +   start(NDIMS), count(NDIMS), start1(2), count1(2)
     
        character*(*) Psi_NAME,Ene_NAME, REC_NAME 
           
                       
        parameter(Psi_NAME='Stream Function')                           
                        
        parameter(Ene_NAME='Energy', REC_NAME='Time')

        integer psi_varid, ene_varid, t_varid

        
        
        real *8 tp_psi(Num_x,Num_y,NLVLS),    !! Define temporary variables
     +       tp_ene(ETYPE)

     
        retval = nf_open(FILE_NAME, nf_write, ncid)                     !! Open file in write mode
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = nf_inq_varid(ncid, Psi_NAME, psi_varid)                !! Request varid
        if (retval .ne. nf_noerr) call handle_err(retval)

        retval = nf_inq_varid(ncid, Ene_NAME, ene_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_inq_varid(ncid, REC_NAME, t_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
 
        
        
        do i = 1, Num_x                                                 !! Write data in temp. variables
            do j = 1, Num_y
               tp_psi(i,j,1) = psi1(i,j)
               tp_psi(i,j,2) = psi2(i,j)

            end do
        end do 
        
        tp_ene(1) = epot
        tp_ene(2) = ekin1
        tp_ene(3) = ekin2
        tp_ene(4) = epot+ekin1+ekin2
        
        count(1) = Num_x                                                                     
        count(2) = Num_y
        count(3) = NLVLS                                                !! Start, count tell the program where to  
        count(4) = 1                                                    !! write data in NetCDF file 
        start(1) = 1
        start(2) = 1
        start(3) = 1
        start(4) = step_tim
        
        retval = nf_put_vara_double(ncid, psi_varid, start, count,      !! Update Data
     +        tp_psi)
        if (retval .ne. nf_noerr) call handle_err(retval)

        
        count1(1) = ETYPE
        count1(2) = 1
        start1(1) = 1
        start1(2) = step_tim
        
        retval = nf_put_vara_double(ncid, ene_varid, start1, count1,
     +       tp_ene)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_vara_double(ncid, t_varid, step_tim, 1,
     +       save_tim)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = nf_close(ncid)                                         !! Close NetCDF file                
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        print *,'!!! Netcdf file modified i.e. ', FILE_NAME, '!!!'
        print *,'!!!                                         !!!'
      end subroutine write_netcdf
      
c -----------------------------------------------------------------------
c -------------- CREATE NETCDF FILE THAT STORES THE TIME-AVERAGED STREAM FUNCTION---
C ---------- AND THE TIME OVER WHICH THE TIME-AVERAGED STREAM FUNCTION WAS TAKEN -----------
      

       
      subroutine create_ave_file(file_name,ii,jj,basinscale)
       
       implicit none
       include 'netcdf.inc'
       
       character*(*) file_name
       integer ii,jj,ncid,retval,i,j
       real*8 basinscale
       
       integer ndims,nlvls
       integer x_dimid,y_dimid,t_dimid,lvl_dimid
       parameter(ndims = 3, nlvls = 2)
       character*(*) x_name, y_name, t_name,lvl_name
       parameter(x_name = 'X', y_name = 'Y', t_name = 'Time')
       parameter(lvl_name = 'Layer')
       character*(*) psiav_name
       parameter(psiav_name = 'Time-averaged Stream Function')
       
       real*8 x_val(ii), y_val(jj), dx,dy
       integer x_varid,y_varid,t_varid
       
       integer psiav_varid
       
       character*(*) units,x_units,y_units,t_units,psi_units
       parameter(units = 'units')
       parameter(psi_units = 'non-dimensionalized by U=1 cm/s, L=dx') 
        parameter (X_UNITS='cms', Y_UNITS='cms', t_UNITS='Days')
        
        integer dimids(ndims),x,y
        
        real, parameter:: START_X=0., START_Y=0.
        
        
        dx = BasinScale/dfloat(ii)
        dy = dx
        do x = 1, ii                                                 !! Values for grid points 
            X_val(x) = START_X + dfloat((x-1))*dx
        end do
        do y = 1, jj
            Y_val(y) = START_Y + dfloat((y-1))*dy
        end do
        
        retval = nf_create(FILE_NAME, nf_clobber, ncid)                 !! Create netcdf file
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = nf_def_dim(ncid, Lvl_NAME, NLVLS, lvl_dimid)           !! Define Dimensions in NetCDF file
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_dim(ncid, X_NAME, ii, x_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_dim(ncid, Y_NAME, jj, y_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_dim(ncid, t_NAME, 1, t_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = nf_def_var(ncid,X_NAME,NF_DOUBLE,1,x_dimid,x_varid)    !! Create variables in NetCDF file
        if (retval .ne. nf_noerr) call handle_err(retval)               !! '+' or '&' in column 6 always
        retval = nf_def_var(ncid,Y_NAME,NF_DOUBLE,1,y_dimid,y_varid) 
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid,t_NAME,NF_DOUBLE,1,t_dimid,t_varid) 
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = nf_put_att_text(ncid, x_varid, UNITS, len(X_UNITS),    !! Assign units
     +        X_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, y_varid, UNITS, len(Y_UNITS), 
     +         Y_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, t_varid, UNITS, len(T_UNITS), 
     +         T_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        dimids(1) = x_dimid
        dimids(2) = y_dimid
        dimids(3) = lvl_dimid
        
        retval = nf_def_var(ncid, Psiav_NAME, NF_DOUBLE, NDIMS, dimids,   !! Creating Main varialbes
     +        psiav_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        
        retval = nf_put_att_text(ncid, psiav_varid, UNITS, 
     +       len(Psi_UNITS), Psi_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
       
       
        retval = nf_enddef(ncid)                                        !! End define mode    
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = nf_put_var_double(ncid, x_varid, X_val)                !! Write data in X, Y
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_var_double(ncid, y_varid, Y_val)
        if (retval .ne. nf_noerr) call handle_err(retval)
            
        retval = nf_close(ncid)                                         !! Close file
        if (retval .ne. nf_noerr) call handle_err(retval)
   
        print *,'!!! Netcdf file created i.e. ', FILE_NAME, '!!!'
        
       
       
       end subroutine create_ave_file
       
C --------------------------------------------------------------------
C --------- WRITE TIME-AVERAGED STREAM FUNCTION TO FILE ------------
C --------------------------------------------------------------------
       
       subroutine write_ave_file(file_name,ii,jj,psi1_av,psi2_av,time)
       
       implicit none
       
       include 'netcdf.inc'
       
       character*(*) file_name
       integer ii,jj,ncid,retval
       real*8 psi1_av(ii,jj), psi2_av(ii,jj), time
       
       integer ndims, nlvls,i,j
       parameter (ndims = 3, nlvls = 2)
       integer psiav_varid,t_varid
       
       character*(*) psiav_name,t_name
       parameter(psiav_name = 'Time-averaged Stream Function')
       parameter(t_name = 'Time')
       integer count(ndims),start(ndims)
       
       real*8 tp_psi(ii,jj,nlvls)
       
               retval = nf_open(FILE_NAME, nf_write, ncid)                     !! Open file in write mode
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = nf_inq_varid(ncid, Psiav_NAME, psiav_varid)                !! Request varid
        if (retval .ne. nf_noerr) call handle_err(retval)

        retval = nf_inq_varid(ncid, t_NAME, t_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        do i = 1, ii                                               !! Write data in temp. variables
            do j = 1, jj
               tp_psi(i,j,1) = psi1_av(i,j)
               tp_psi(i,j,2) = psi2_av(i,j)

            end do
        end do 
        
        count(1) = ii
        count(2) = jj
        count(3) = nlvls
        start(1) = 1
        start(2) = 1
        start(3) = 1
        
        retval = nf_put_vara_double(ncid,psiav_varid,start,count,tp_psi)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = nf_put_vara_double(ncid,t_varid,1,1,time)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = nf_close(ncid)                                         !! Close NetCDF file                
        if (retval .ne. nf_noerr) call handle_err(retval)
c        print *,'!!! Netcdf file modified i.e. ', FILE_NAME, '!!!'
c        print *,'!!!                                         !!!'
      end subroutine write_ave_file
      
c --------------------------------------------------------------------------
c ----------------- READ TIME-AVERAGED STREAM FUNCTION -----------------------
C --------------------------------------------------------------------------------
      
      subroutine read_ave_file(file_name,ii,jj,psi1_av,psi2_av,time_av)
      
      implicit none
      include 'netcdf.inc'
      
      character*(*) file_name
      integer ii,jj
      real*8 psi1_av(ii,jj), psi2_av(ii,jj),time_av
      
      integer ncid,retval,ndims,nlvls
      parameter(nlvls = 2,ndims=3)
      
      character*(*) psiav_name,t_name
      parameter(psiav_name = 'Time-averaged Stream Function')
      parameter(t_name = 'Time')
      
      integer psiav_varid,time_varid
      integer count(ndims), start(ndims)
      
        retval = nf_open(FILE_NAME, nf_nowrite, ncid)                   !! Open file in read mode
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = nf_inq_varid(ncid, Psiav_NAME, psiav_varid)                !! Request varid
        if (retval .ne. nf_noerr) call handle_err(retval)

        
        count(1) = ii
        count(2) = jj
        count(3) = 1
        start(1) = 1
        start(2) = 1
        start(3) = 1
        
        retval = nf_get_vara_double(ncid, psiav_varid, start, count,      !! Read Psi
     +        psi1_av)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        start(3) = 2
        retval = nf_get_vara_double(ncid, psiav_varid, start, count,
     +        psi2_av)
        if (retval .ne. nf_noerr) call handle_err(retval)   
        
        retval = nf_inq_varid(ncid, t_NAME, time_varid)                !! Request varid
        if (retval .ne. nf_noerr) call handle_err(retval)   
        
        retval = nf_get_vara_double(ncid,time_varid,1,1,time_av)
        if (retval .ne. nf_noerr) call handle_err(retval)  
          
     
        retval = nf_close(ncid)                                         !! Close file
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        print*,'Time-averaged stream function read'
        
        
        
       end subroutine read_ave_file
       
C ------------------------------------------------------------------------------
C ----------- READ TIME FROM CABARET NETCDF FILE --------------------------------
C -------------------------------------------------------------------------
       
        subroutine read_time(file_name,time,t_len)
            !use time_array
            implicit none
            include 'netcdf.inc'
        
            character*(*) file_name
            real*8, dimension (:), allocatable :: time
            integer t_len
        
            integer retval, ncid
            
            character*(*) t_name
            parameter(t_name = 'Time')
            integer t_dimid, t_varid
            
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
              
              retval = nf_inq_dimid(ncid,t_name,t_dimid)
              if (retval .ne. nf_noerr) then
                print *, 'timee dimid not found'
               call handle_err(retval)
               end if
              retval = nf_inq_dimlen(ncid,t_dimid,t_len)
              if (retval .ne. nf_noerr) call handle_err(retval)
              
              
              allocate(time(t_len))
              
              retval = nf_get_vara_double(ncid, t_varid, 1, t_len, time)   ! Read time first
              if (retval .ne. nf_noerr) call handle_err(retval)
              
              retval = nf_close(ncid)
              if (retval .ne. nf_noerr) call handle_err(retval)
              
              !deallocate(tp_time,u,v)
              
              write(*,*) 'time read'
              
            end subroutine read_time
            
c        subroutine psi_av(file_name,ii,jj,psi1_av,psi2_av)
      
c      implicit none
c      include 'netcdf.inc'
      
c      character*(*) file_name
c      integer ii,jj
c      real*8 psi1_av(ii,jj), psi2_av(ii,jj)
      
c      integer ncid,retval,ndims,nlvls
c      parameter(nlvls = 2,ndims=3)
      
c      character*(*) psiav_name
c      parameter(psiav_name = 'Time-averaged Stream Function')
      
c      integer psiav_varid
c      integer count(ndims), start(ndims)
      
c        retval = nf_open(FILE_NAME, nf_nowrite, ncid)                   !! Open file in read mode
c        if (retval .ne. nf_noerr) call handle_err(retval)
        
c        retval = nf_inq_varid(ncid, Psiav_NAME, psiav_varid)                !! Request varid
c        if (retval .ne. nf_noerr) call handle_err(retval)

        
c        count(1) = ii
c        count(2) = jj
c        count(3) = 1
c        start(1) = 1
c        start(2) = 1
c        start(3) = 1
        
c        retval = nf_get_vara_double(ncid, psiav_varid, start, count,      !! Read Psi
c     +        psi1_av)
c        if (retval .ne. nf_noerr) call handle_err(retval)
        
c        start(3) = 2
c        retval = nf_get_vara_double(ncid, psiav_varid, start, count,
c     +        psi2_av)
c        if (retval .ne. nf_noerr) call handle_err(retval)        
     
c        retval = nf_close(ncid)                                         !! Close file
c        if (retval .ne. nf_noerr) call handle_err(retval)
        
c        !print*,'psi_av read'
        
       
c        end subroutine psi_av
      
       

      end module MOD_qg2_NETCDF
