c MODULE THAT CONTAINS ERRO HANDLER FOR NETCDF FILES

      module mod_netcdf_error
      
      contains
      
c    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
c    !! This subroutine just checks for errors while writing data in netcdf.  
 
       subroutine handle_err(errcode)
        implicit none
        include 'netcdf.inc'
        integer errcode

        print *, 'Error: ', nf_strerror(errcode)
        stop 2
       end subroutine handle_err
      
      end module mod_netcdf_error
