      module stochastic_parameters
      
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

      subroutine velocity_variance(u,npoints,nt,i,j,sigma)
      
      implicit none
      
      integer npoints,nt,n,i,j
      real*8 x(2,npoints,nt),sigma(npoints),u(2,npoints,nt)
      
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
      
      call velocity_variance(u,npoints,nt,i,j,sigma)
      
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

      subroutine drift_correction_1(x_new,x_old,npoints,i,a)
      
      implicit none
      
      integer nt,npoints,i
      real*8 sigma(2,2,npoints),theta(npoints),a1_drift(npoints)
      real*8 x(2,npoints,nt),u(2,npoints,nt-1)
      
      integer d,m,j,k
      
      a = 0.
      
      do d = 1,2
        call lagrangian_velocity(x_new,dt,npoints,nt,d,u(d,:,:))
        call lagrangian_velocity(
      enddo
      
      do m = 1,2
      do j = 1,2
        call velocity_variance(u,npoints,nt,m,j,sigma(m,j,:))
      enddo
      enddo
      
      ! calculate the first term
      do m = 1,2
        a = a + .5*(sigma(i,m)
      enddo
      
      do m = 1,2
      do j = 1,2
      do k = 1,2
      
      enddo
      enddo
      enddo
      
      
      
      
      end subroutine drift_correction_1
      
c -----------------------------------------------------------------------------

      subroutine random_forcing_1
      
      end subroutine random_forcing_1
      
c -----------------------------------------------------------------------------
      
      subroutine markov1_tensor(x,dt,npoints,nt,ntau,i,j,theta)
      
      implicit none
      
      integer npoints,nt,t,n
      real*8 R(ntau),x(2,npoints,nt),theta(npoints)
      
      ! calculate autocorrelation
      
      call autocorrelation(x,dt,npoints,nt,ntau,i,j,R)
      
      ! find when R reaches zero
      
      do n = 1,npoints
      do t = 1,ntau
        if (R(n,t) <=0) then 
            theta(n) = t
            break
        endif
      enddo
      enddo
      
      end subroutine markov1_tensor
      
c -----------------------------------------------------------------------------


      
      
      
    
      
      end module stochastic_parameters
