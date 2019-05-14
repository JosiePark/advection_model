      program variance_derivative
      
      use mod_stochastic_parameters
      use mod_vel_variance_netcdf
      use mod_2dcubic
      
      implicit none
      
      character*(*), parameter :: home_dir = 
     & '/home/clustor2/ma/j/jp1115/DATA/1/'
     
      character*(*), parameter :: sigma_file = 
     &   trim(home_dir) // 'STATS/SIGMA/velocity_variance.nc'
     
      integer ii,jj,ii1,jj1
      parameter(ii = 512, jj = 512)
      parameter(ii1 = 512*2, jj1 = 512*2)
      real*8 tmp1(2,2,ii,jj),tmp2(2,2,ii,jj)
      real*8 sigma11(ii,jj),sigma12(ii,jj),sigma21(ii,jj),sigma22(ii,jj)
      real*8 asigma11(ii,jj),bsigma11(ii,jj)
     & ,csigma11(ii,jj),dsigma11(ii,jj) ! zonal, top layer
      real*8 asigma12(ii,jj),bsigma12(ii,jj)
     & ,csigma12(ii,jj),dsigma12(ii,jj) ! meridional, top layer
      real*8 asigma21(ii,jj),bsigma21(ii,jj)
     & ,csigma21(ii,jj),dsigma21(ii,jj) ! zonal, bottom layer
      real*8 asigma22(ii,jj),bsigma22(ii,jj)
     & ,csigma22(ii,jj),dsigma22(ii,jj) ! meridional, bottom layer
      real*8 adsigma11(ii,jj),bdsigma11(ii,jj)
     & ,cdsigma11(ii,jj),ddsigma11(ii,jj)
      real*8 adsigma12(ii,jj),bdsigma12(ii,jj)
     & ,cdsigma12(ii,jj),ddsigma12(ii,jj)
      real*8 adsigma21(ii,jj),bdsigma21(ii,jj)
     & ,cdsigma21(ii,jj),ddsigma21(ii,jj)
      real*8 adsigma22(ii,jj),bdsigma22(ii,jj)
     & ,cdsigma22(ii,jj),ddsigma22(ii,jj)
      real*8 d_sigma12(ii,jj),d_sigma11(ii,jj)
     & ,d_sigma22(ii,jj),d_sigma21(ii,jj) ! derivatives of sigma
     
      real*8 sigma11_tmp,sigma12_tmp,sigma21_tmp,sigma22_tmp
      real*8 d_sigma11_tmp,d_sigma12_tmp,d_sigma21_tmp,d_sigma22_tmp
      
      real*8 sigma11_interp(ii1,jj1),sigma12_interp(ii1,jj1)
     & ,sigma21_interp(ii1,jj1),sigma22_interp(ii1,jj1)
      real*8 dsigma11_interp(ii1,jj1),dsigma12_interp(ii1,jj1)
     & ,dsigma21_interp(ii1,jj1),dsigma22_interp(ii1,jj1)
     
      real*8 x_c(ii1),y_c(jj1)
     
      integer i,j
     
c READ VELOCITY VARIANCE FROM FILE 

      call read_vel_variance(sigma_file,tmp1,tmp2,ii,jj)
      sigma11 = tmp1(1,1,:,:)
      sigma12 = tmp1(2,2,:,:)
      sigma21 = tmp2(1,1,:,:)
      sigma22 = tmp2(2,2,:,:)
      
c CALCULATE COEFFICIENTS

      call cubic_coeff_x(ii,jj,sigma11,asigma11
     & ,bsigma11,csigma11,dsigma11)
      call cubic_coeff_x(ii,jj,sigma12,asigma12
     & ,bsigma12,csigma12,dsigma12)
      call cubic_coeff_x(ii,jj,sigma21,asigma21
     & ,bsigma21,csigma21,dsigma21)
      call cubic_coeff_x(ii,jj,sigma22,asigma22
     & ,bsigma22,csigma22,dsigma22)
     
C DETERMINE POINTS AT WHICH YOU WANT TO INTERPOLATE

      do i = 1,ii1
        x_c(i) = dfloat(i-1)/4.
        y_c(i) = dfloat(i-1)/4.
      enddo
     
c INTERPOLATE SIGMA AND CALCULATE DERIVATIVE

      do i = 1,ii1
      do j = 1,jj1
      
      call sigma_interpolation(ii,jj,asigma11,bsigma11
     & ,csigma11,dsigma11,1
     & ,x_c(i),y_c(j),sigma11_tmp,d_sigma11_tmp)
      call sigma_interpolation(ii,jj,asigma12,bsigma12
     & ,csigma12,dsigma12,2
     & ,x_c(i),y_c(j),sigma12_tmp,d_sigma12_tmp)
      call sigma_interpolation(ii,jj,asigma21,bsigma21
     & ,csigma21,dsigma21,1
     & ,x_c(i),y_c(j),sigma21_tmp,d_sigma21_tmp)
      call sigma_interpolation(ii,jj,asigma22,bsigma22
     & ,csigma22,dsigma22,2
     & ,x_c(i),y_c(j),sigma22_tmp,d_sigma22_tmp)
     
      sigma11_interp(i,j) = sigma11_tmp
      sigma12_interp(i,j) = sigma12_tmp
      sigma21_interp(i,j) = sigma21_tmp
      sigma22_interp(i,j) = sigma22_tmp
      
      dsigma11_interp(i,j) = d_sigma11_tmp
      dsigma12_interp(i,j) = d_sigma12_tmp
      dsigma21_interp(i,j) = d_sigma21_tmp
      dsigma22_interp(i,j) = d_sigma22_tmp
      
      enddo
      enddo

c WRITE SIGMA TO FILE

      open(1, file = 'sigma.dat',status = 'new')
      
      do i = 1,ii1
      do j = 1,jj1
      write(1,*) sigma11_interp(i,j)
      enddo
      enddo
      
      close(1)

c WRITE DERIVATE OF SIGMA TO FILE

      open(2, file = 'dsigma.dat',status = 'new')
      
      do i = 1,ii1
      do j = 1,jj1
      write(2,*) dsigma11_interp(i,j)
      enddo
      enddo
      
      close(2)
      
      
      end program variance_derivative
