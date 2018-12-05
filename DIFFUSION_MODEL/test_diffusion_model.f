c CODE THAT TESTS COMPONENTS OF THE DIFFUSION MODEL - I.E INTERPOLATION
C AND DIFFERENTIATION OF THE DIFFUSIVITY COEFFICIENT

      program test_diffusion_model
      
      use mod_diff_input
      use mod_diff_constants
      use mod_diffusivity_functions
      use mod_qg2_netcdf
      use mod_random
      use mod_variables
      use mod_2dcubic
      use mod_diffusion_netcdf
      use mod_1dinterp
      use mod_test_diffusivity_netcdf

      implicit none
      
C READ TIME-AVERAGED STREAM FUNCTION

      call read_ave_file(ave_file,ii,jj,psi1_av,psi2_av,time_av)
      
C CALCULATE COEFFICIENTS FOR 2D-CUBIC SPATIAL INTERPOLATION

      call cubic_coeff_x(ii,jj,psi1_av,a1,b1,c1,d1)
      call cubic_coeff_x(ii,jj,psi2_av,a2,b2,c2,d2)
      
      print*,'cubic_coeff_done'
      
C READ THE DIFFUSIVITY TENSOR

      call read_diffusivity(diff_file,K)
      
      print*, 'diffusivity read'
      
      
c DETERMINE BINS

      do b = 1,nbins+1
      
      bin_corners(b) = (b-1)*(dfloat(jj)/nbins)
      
      enddo

      
c CALCULATE BIN CENTRES

      do b = 1,nbins
      
        bin_centres(b) = .5*(bin_corners(b)+bin_corners(b+1))
      
      enddo
      
c CREATE DIFFUSIVITY FILE

      call create_diffusivity_file(test_file,ii)
      
      
c DETERMINE GRID POINTS 

      do i = 1,ii
        print*, i

            
            
            
      
C SPATIALLY INTERPOLATE TO FIND TIME-AVERAGED VELOCITY AT PARTICLE LOCATION


                  call cubic_poly_x(ii,jj,dfloat(1),dfloat(i)
     &                  ,a1,b1,c1,d1,psi1_x)
     
                  call vel(ii,jj,psi1_x,a1,b1,c1,d1,dfloat(1),dfloat(i)
     &             ,u1,v1)
     

     
c INTERPOLATE THE EDDY-DIFFUSIVITY AT THE PARTICLE LOCATION

                  call diffusivity_interp_1d(bin_centres,nbins,K(1,:,1)
     &                  ,ii,dfloat(i),Kx1)
                  call diffusivity_interp_1d(bin_centres,nbins,K(2,:,1)
     &                  ,ii,dfloat(i),Ky1)
                  call diffusivity_interp_1d(bin_centres,nbins,K(1,:,2)
     &                  ,ii,dfloat(i),Kx2)
                  call diffusivity_interp_1d(bin_centres,nbins,K(2,:,2)
     &                  ,ii,dfloat(i),Ky2)
                  

                  
c APPROXIMATE DERIVATIVE OF THE DIFFUSIVITY - I.E THE DRIFT TERM    
              
                  dKdx1 = 0.
                  dKdx2 = 0.
                  call diffusivity_derivative(bin_centres,nbins
     &                  ,K(2,:,1),ii,dfloat(i),dKdy1)
                  call diffusivity_derivative(bin_centres,nbins
     &                  ,K(2,:,2),ii,dfloat(i),dKdy2)
     
        
C WRITE DIFFUSIVITIES TO FILE

      call write_diffusivity_file(test_file,Kx1,Ky1,Kx2,Ky2,
     & dKdx1,dKdy1,dKdx2,dKdy2,i)


      enddo
      
    

      end program test_diffusion_model
