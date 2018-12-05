c Code that saves the location of the jet core as
c approximated from the maximum value of the PV anomaly gradient

      program jetcore
      
      use mod_constants
      use mod_parameter 
      use mod_qg2_netcdf
      use mod_time_array
      use mod_variables
      
      
      
      implicit none
      
      integer t_len,n
      

      call read_time(file_name,time,t_len)
      call create_PV_file(PV_name,ii,jj,basinscale)
      
      allocate(grad(ii,t_len))
      
      do n = 1,t_len
      
      print*, 'time =',time(n)
      
      call read_NETCDF(FILE_NAME,psi1,psi2,ii,jj,time(n),
     +  new_tim, step_tim)
     
        print*, 'psi read'
     
c calculate relative vorticity and vorticity anomaly

      call rel_from_psi(ii,jj,psi1,psi2,rel1,rel2)
      
      call zeta_from_psi_and_rel(ii,jj,S1,S2,psi1,psi2,rel1,rel2
     &                                                    ,z1,z2)
     
        call write_PV(PV_NAME,z1,z2,
     +    ii,jj,time(n),n)
      
     
      enddo
      

      
      end program jetcore
      
       subroutine create_PV_file(FILE_NAME,Num_x,Num_y,BasinScale)
        implicit none
        include 'netcdf.inc'
            
        integer Num_x, Num_y                                            !! Declare Variables
        real*8 BasinScale 
        character*(*) FILE_NAME
        integer NDIMS, NLVLS, ETYPE
        integer ncid, x_dimid, y_dimid, rec_dimid, lvl_dimid, Et_dimid
        parameter (NDIMS = 4, NLVLS=2, ETYPE=4)  
        character*(*) X_NAME, Y_NAME, REC_NAME, Lvl_NAME
        parameter (X_NAME='X', Y_NAME='Y', REC_NAME='Time')
        parameter (Lvl_NAME = 'Layer')
        integer start(NDIMS), count(NDIMS) 
            
        real*8 X_val(Num_x), Y_val(Num_y), dx, dy                       !! Create variables to hold grid data
        integer x_varid, y_varid, t_varid
        
        character*(*) Psi_name
             !! Assign names to the variaables which  
                       !! will be used in the netcdf file  
        parameter(Psi_NAME='PV anomaly')                           !! 'X', 'Y', 'Time' will be use for grid and
                      !! time storage   
        

        integer psi_varid, ene_varid, 
     &          dimids(NDIMS), dimid(2)
        
        character*(*) UNITS, X_UNITS, Y_UNITS, Time_UNITS,    !! Assign Units to variables
     &                Psi_UNITS
        parameter (UNITS='units')

        parameter(Psi_UNITS = 'non-dimensionalized by U=1 cm/s, L=dx') 

 
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

       end
       
              subroutine write_PV(FILE_NAME,psi1,psi2,
     +    Num_x,Num_y,save_tim,step_tim)
     
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
     
        character*(*) Psi_NAME, REC_NAME 
           
                       
        parameter(Psi_NAME='PV anomaly')                           
                        
        parameter( REC_NAME='Time')

        integer psi_varid, ene_varid, t_varid

        
        
        real *8 tp_psi(Num_x,Num_y,NLVLS),    !! Define temporary variables
     +       tp_ene(ETYPE)

     
        retval = nf_open(FILE_NAME, nf_write, ncid)                     !! Open file in write mode
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = nf_inq_varid(ncid, Psi_NAME, psi_varid)                !! Request varid
        if (retval .ne. nf_noerr) call handle_err(retval)


        retval = nf_inq_varid(ncid, REC_NAME, t_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
 
        
        
        do i = 1, Num_x                                                 !! Write data in temp. variables
            do j = 1, Num_y
               tp_psi(i,j,1) = psi1(i,j)
               tp_psi(i,j,2) = psi2(i,j)

            end do
        end do 
        

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


        retval = nf_put_vara_double(ncid, t_varid, step_tim, 1,
     +       save_tim)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = nf_close(ncid)                                         !! Close NetCDF file                
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        print *,'!!! Netcdf file modified i.e. ', FILE_NAME, '!!!'
        print *,'!!!                                         !!!'
      end
      
             subroutine handle_err(errcode)
        implicit none
        include 'netcdf.inc'
        integer errcode

        print *, 'Error: ', nf_strerror(errcode)
        stop 2
       end

      
      
      
