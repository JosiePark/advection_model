      module mod_stochastic_parameters
      
      use mod_2dcubic
      use mod_diffusivity_functions
      
      contains
      
      subroutine lagrangian_velocity(x,dt,npoints,nt,i,u)
      
      implicit none
      
      integer i,t,npoints,nt
      real*8 x(2,npoints,nt),u(npoints,nt-1),dt
      
      do t = 1,nt-1
        u(:,t) = (x(i,:,t+1)-x(i,:,t))/dt
      enddo
      
      end subroutine lagrangian_velocity
      
c -----------------------------------------------------------------------------

      subroutine lagrangian_velocity_variance(u,npoints,nt,i,j,sigma)
      
      implicit none
      
      integer npoints,nt,n,i,j
      real*8 x(2,npoints,nt),sigma(2,2,npoints),u(2,npoints,nt)
      
      do n = 1,npoints
        sigma = DOT_PRODUCT(u(i,n,:),u(j,n,:))
      enddo
      
      end subroutine lagrangian_velocity_variance

c -----------------------------------------------------------------------------

      subroutine autocorrelation(x,dt,npoints,nt,ntau,i,j,R)
      
      implicit none
      
      integer npoints,ntau,nt,i,j,d,m
      real*8 dt
      real*8 x(2,npoints,nt),R(npoints,ntau),u(2,npoints,nt-1)
      real*8 sigma(2,npoints)
      integer tau,n,t
      
      ! calculate lagrangian velocity
      
      do d = 1,2
        call lagrangian_velocity(x,dt,npoints,nt,d,u(d,:,:))
        call lagrangian_velocity_variance(u(d,:,:)
     &   ,npoints,nt,i,j,sigma(d,:))
      enddo
      
      
      ! calculate non-normalised autocorrelation
      
      
      do tau = 0,ntau-1
        do t = 1,(nt-tau-1)
            R(:,tau+1) = R(:,tau+1)+u(i,:,t)*u(j,:,t+tau)
        enddo
      enddo
      
      ! normalise and average R
      
      do n = 1,npoints
        R(n,:) = R(n,:)/sqrt(sigma(i,n)*sigma(j,n))
      enddo
      
      end subroutine autocorrelation
      
c -----------------------------------------------------------------------------

      subroutine random_forcing_1(sigma
     & ,theta,ii,jj,nbins,bin_centres,b)
      
      implicit none
      
      integer ii,jj,k,m,j,i
      real*8 sigma(2,ii,jj),theta(nbins,2)
      real*8 theta_interp(2,jj)
      real*8 b(2,ii,jj)
      integer nbins
      real*8 bin_centres(nbins)
      
      
      ! Interpolate theta meridionally
      
      do k = 1,2
      do j = 1,jj
        call diffusivity_interp_1d(bin_centres,nbins,theta(:,k)
     &   ,ii,dfloat(j-1),theta_interp(j,k))
      enddo
      enddo


      do i = 1,ii
      do j = 1,jj
      do k = 1,2
      
        b(k,i,j) = sqrt(2*sigma(k,i,j)/theta_interp(j,k))
      
      enddo
      enddo
      enddo
      
      end subroutine random_forcing_1
      
c -----------------------------------------------------------------------------

      subroutine variance_derivative(sigma,ii,jj,dsigma)

      implicit none
            
      integer ii,jj,i,j,k
      real*8 sigma(2,ii,jj),dsigma(2,ii,jj)
      
      integer ix,iy,im1,ip1,jm1,jp1
      
      do ix = 1,ii
      do iy = 1,jj
      
            im1=ix-1
            ip1=ix+1
            jm1=iy-1
            jp1=iy+1
            if(im1.eq.0)    im1=ii
            if(ip1.eq.ii+1) ip1=1
            if(jm1.eq.0)    jm1=jj
            if(jp1.eq.jj+1) jp1=1
            
            dsigma(1,ix,iy) = .5*(sigma(1,ip1,iy)-sigma(1,im1,iy))
            dsigma(2,ix,iy) = .5*(sigma(2,ix,jp1)-sigma(2,ix,jm1))


      enddo
      enddo
      
      
      
      end subroutine variance_derivative
      
c -----------------------------------------------------------------------------

      subroutine sigma_interpolation(ii,jj,a,b,c,d,k,x,y,sigma,dsigma)
      
      implicit none
      
      integer ii,jj
      real*8 a(ii,jj),b(ii,jj),c(ii,jj),d(ii,jj)
      real*8 x,y
      real*8 sigma,dsigma
      real*8 sigma_x(4)
      real*8 alpha,beta,gamma,delta
      real*8 d_alpha,d_beta,d_gamma,d_delta
      real*8 ay,ax
      integer yc,xc
      integer i,j,jc(4)
      integer k

      ay = y - int(y)
      yc = int(y) + 1
  
      xc = int(x) + 1
    
      ax = x - int(x)
        
        do i = 1,4
            jc(i) = yc - 2 + i
            if (jc(i) <= 0) then
            jc(i) = jj + jc(i)
            elseif (jc(i) > jj) then
            jc(i) = jc(i) - jj
            endif
 
        enddo
      call cubic_poly_x(ii,jj,x,y,a,b,c,d
     & ,sigma_x)
      call cubic_coeff_y(ii,jj,sigma_x,alpha,beta,gamma,delta)
      
      sigma = alpha + beta*ay + gamma*ay**2 
     &    + delta*ay**3
        
        if (k .eq. 2) then
        dsigma = (beta + 2*gamma*ay + 3*delta*ay**2)
        
        else
        
        d_alpha = b(xc,jc(2))
     &    + 2*c(xc,jc(2))*ax + 3*d(xc,jc(2))*ax**2
        
        d_beta = (-b(xc,jc(1))/3 - b(xc,jc(2))/2
     &    + b(xc,jc(3))-b(xc,jc(4))/6)
     & + 2*(-c(xc,jc(1))/3 - c(xc,jc(2))/2 
     & + c(xc,jc(3)) - c(xc,jc(4))/6)*ax
     & + 3*(-d(xc,jc(1))/3 - d(xc,jc(2))/2
     & + d(xc,jc(3)) - d(xc,jc(4))/6)*ax**2
     
        d_gamma = (b(xc,jc(1))/2 - b(xc,jc(2)) + b(xc,jc(3))/2)
     & + 2*(c(xc,jc(1))/2 - c(xc,jc(2)) + c(xc,jc(3))/2)*ax
     & + 3*(d(xc,jc(1))/2 - d(xc,jc(2)) + d(xc,jc(3))/2)*ax**2
                                                                                                             
        d_delta = (-b(xc,jc(1))/6 + b(xc,jc(2))/2 
     &   - b(xc,jc(3))/2 + b(xc,jc(4))/6)
     & + 2*(-c(xc,jc(1))/6 + c(xc,jc(2))/2
     & - c(xc,jc(3))/2 + c(xc,jc(4))/6)*ax
     & + 3*(-d(xc,jc(1))/6 + d(xc,jc(2))/2
     & - d(xc,jc(3))/2 + d(xc,jc(4))/6)*ax**2

     
        dsigma= (d_alpha + ay*d_beta + d_gamma*ay**2 + d_delta*ay**3)
        
        endif
        
        
      
      end subroutine sigma_interpolation
      

      end module mod_stochastic_parameters
