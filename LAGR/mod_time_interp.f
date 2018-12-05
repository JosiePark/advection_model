c Module that contains subroutine to interpolate the time array using
c cubic interpolation to determine the half time step for use in
c RK4 time integration

      module MOD_time_interp
      
      
      
              contains
        
        recursive subroutine interp_time(ii,jj,time_cubic,time_half
     &  ,psi,psi_half)
        implicit none
        
        integer ii,jj,k,i,j
        real*8 psi(ii,jj,4), psi_half(ii,jj), time_cubic(4), time_half
        real*8 li
        
        do i = 1,ii
            do j = 1,jj
                psi_half(i,j)=0
            enddo
        enddo
        
        do i = 1,ii
            do j = 1,jj
             do k = 1,4
                li = lagrange_function(k,time_half,time_cubic)
                psi_half(i,j) = psi_half(i,j) + psi(i,j,k)*li
            enddo
            enddo
        enddo
        
        contains 
        function lagrange_function(k,time_half
     &   ,time_cubic) result(fi)
              
            implicit none
            integer k,i
            real*8 time_half,time_cubic(4)
       
            real*8 fi
       
            fi = 1
       
            if (time_half.eq.time_cubic(k)) then
                fi = 1
            else
                do i = 1,4
                if (i .ne. k) then
                    fi = fi*(time_half-time_cubic(i))/
     &                    (time_cubic(k)-time_cubic(i))
                end if
                enddo
            endif
        
            end function lagrange_function
        end subroutine interp_time
        
            subroutine read_psi(file_name,ii,jj,time_interp,psi1,psi2)
        
        
            implicit none
            
            include 'netcdf.inc'
            
            integer ii,jj, time_interp
            real*8 psi1(ii,jj,4), psi2(ii,jj,4)
            
            character*(*) file_name
            
            integer retval, ncid, nlvls, ndims
            parameter(ndims=4,nlvls=2)
            integer dimids(ndims)
            integer psi_varid, psi_dimid, start(4),count(4)
            
            character*(*) psi_name
            parameter(psi_name='Stream Function')
            
            retval = nf_open(file_name,nf_nowrite,ncid)
            if (retval .ne. nf_noerr) call handle_err(retval)
            
            retval = nf_inq_varid(ncid,Psi_name,psi_varid)
            if (retval .ne. nf_noerr) then 
            print*,'psi_varid not found'
            call handle_err(retval)
            end if
            
            count(1) = ii
            count(2) = jj
            count(3) = 1
            count(4) = 4
            start(1) = 1
            start(2) = 1
            start(3) = 1
            start(4) = time_interp
            
            retval = nf_get_vara_double(ncid,psi_varid,start,count,psi1)
            if (retval .ne. nf_noerr) call handle_err(retval)
            
            start(3) = 2
            retval = nf_get_vara_double(ncid,psi_varid,start,count,psi2)
            if (retval .ne. nf_noerr) call handle_err(retval)
            
            retval = nf_close(ncid)
            if (retval .ne. nf_noerr) call handle_err(retval)
            
            !print*, 'psi read'
            
            
            !return
            
            end subroutine read_psi 
            
c -------------------------------------------------------------------------
c -------------------------------------------------------------------------- 

            subroutine read_eof(file_name,ii,jj,time_interp,psi1,psi2)
        
        
            implicit none
            
            include 'netcdf.inc'
            
            integer ii,jj, time_interp
            real*8 psi1(ii,jj,4), psi2(ii,jj,4)
            
            character*(*) file_name
            
            integer retval, ncid, nlvls, ndims
            parameter(ndims=4,nlvls=2)
            integer dimids(ndims)
            integer psi_varid, psi_dimid, start(4),count(4)
            
            character*(*) psi_name
            parameter(psi_name='Stream Function')
            
            retval = nf_open(file_name,nf_nowrite,ncid)
            if (retval .ne. nf_noerr) call handle_err(retval)
            
            retval = nf_inq_varid(ncid,Psi_name,psi_varid)
            if (retval .ne. nf_noerr) then 
            print*,'psi_varid not found'
            call handle_err(retval)
            end if
            
            count(1) = ii
            count(2) = jj
            count(3) = 1
            count(4) = 4
            start(1) = 1
            start(2) = 1
            start(3) = 1
            start(4) = time_interp
            
            retval = nf_get_vara_double(ncid,psi_varid,start,count,psi1)
            if (retval .ne. nf_noerr) call handle_err(retval)
            
            start(3) = 2
            retval = nf_get_vara_double(ncid,psi_varid,start,count,psi2)
            if (retval .ne. nf_noerr) call handle_err(retval)
            
            retval = nf_close(ncid)
            if (retval .ne. nf_noerr) call handle_err(retval)
            
            !print*, 'psi read'
            
            
            !return
            
            end subroutine read_eof
            
                        subroutine handle_err(errcode)
            implicit none
            include 'netcdf.inc'
            integer errcode

            print *, 'Error: ', nf_strerror(errcode)
            stop 
            end subroutine
      
      end module MOD_time_interp
