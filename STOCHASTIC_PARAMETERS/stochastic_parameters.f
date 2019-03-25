      module stochastic_parameters
      
      use mod_2dcubic
      
      contains
      
      subroutine lagrangian_velocity(x,dt,npoints,nt,i,u)
      
      implicit none
      
      integer i,t,npoints,nt
      real*8 x(2,npoints,nt),u(npoints,nt-1)
      
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
      
      end subroutine velocity_variance

c -----------------------------------------------------------------------------

      subroutine autocorrelation(x,dt,npoints,nt,ntau,i,j,R)
      
      implicit none
      
      integer npoints,ntau,nt,i,j,d,j,m
      real*8 dt
      real*8 x(2,npoints,nt),R(npoints,ntau),u(2,npoints,nt-1)
      real*8 sigma(2,2,npoints)
      integer tau,n
      
      ! calculate lagrangian velocity
      
      do d = 1,2
        call lagrangian_velocity(x,dt,npoints,nt,d,u(d,:,:))
      enddo
      
      ! calculate velocity variance
      
      call lagrangian_velocity_variance(u,npoints,nt,i,j,sigma)
      
      ! calculate non-normalised autocorrelation
      
      
      do tau = 0,ntau-1
        do t = 1,(nt-tau-1)
            R(:,tau+1) = R(:,tau+1)+u(i,:,t)*u(j,:,t+tau)
        enddo
      enddo
      
      ! normalise and average R
      
      do n = 1,npoints
        tmp(n,:) = tmp(n,:)/sigma(n)
      enddo
      
      end subroutine autocorrelation
      
c -----------------------------------------------------------------------------

      subroutine random_forcing_1(sigma,asigma,bisgma,csigma,dsigma
     & ,theta,ii,jj,nbins,b)
      
      implicit none
      
      integer ii,jj,k,m,j
      real*8 sigma(2,2,ii,jj),theta(nbins,2,2),b(ii,jj,2,2)
      real*8 asigma(2,2,ii,jj),bsigma(2,2,ii,jj)
     & ,csigma(2,2,ii,jj),dsigma(2,2,ii,jj)
      real*8 theta_interp(jj,2,2)
      
      ! Interpolate theta meridionally
      
      do k = 1,2
      do m = 1,2
      do j = 1,jj
        call diffusivity_interp_1d(bin_centres,nbins,theta(:,k,m)
     &   ,ii,j,theta_interp(j,k,m)
      enddo
      enddo
      enddo
      
      
      ! first calculate b11,b22
      
      
      
      
      
      end subroutine random_forcing_1
      
c -----------------------------------------------------------------------------

      subroutine variance_derivative(sigma,ii,jj,dsigma)
    
c i,j denotes coordinate of velocity variance tensor
c k denotes which component of the derivative you wish to differentiate
c sigma with respect to.

      implicit none
            
      integer ii,jj,i,j,k
      real*8 sigma(2,2,ii,jj),dsigma(2,2,2,ii,jj)
      
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
            
            do i = 1,2
            do j = 1,2
                dsigma(i,j,1,ix,iy) = 
     &                .5*(sigma(i,j,ip1,iy)-sigma(i,j,im1,iy))
                dsigma(i,j,2,ix,iy) = 
     &                .5*(sigma(i,j,ix,jp1)-sigma(i,j,ix,jm1))
            enddo
            enddo


      enddo
      enddo
      
      
      
      end subroutine variance_derivative

      end module stochastic_parameters
