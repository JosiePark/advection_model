c module that performs one dimensional cubic lagrange interpolation

      module mod_1dinterp
      
      implicit none
      
      contains
      
      subroutine interp_1d(x_c,nx,f_c,x,f)
      
c x_c is the grid on which the gridded data is saved on
c nx is the size of x_c
c f_c is the gridded data
c x is the point at which you wish to perform interpolation
c f is the value of f_c at x
      
      implicit none
      
      integer nx,k,i,j
      
      real*8 x_c(nx), f_c(nx), x, f,l_interp(nx)

      !l_interp = 1.0D0
      f = 0.0D0  

      do i = 1,nx
      l_interp(i) = 1.0D0
      do j = 1,nx
      if( i /= j) then
        l_interp(i) = ((x-x_c(j))/(x_c(i)-x_c(j)))*l_interp(i)
      endif
      enddo
      f = f + l_interp(i)*f_c(i)
      enddo

      end subroutine interp_1d
      
      subroutine linear_interp_1d(x_c,nx,f_c,x,f)
      
      implicit none
      
      integer nx,n
      real*8 x_c(nx),f_c(nx),x,f,xL,xR,fL,fR
      
      ! find neighbouring data points
      
      do n = 1,nx
        if (x_c(n)<x) then
            xL = x_c(n)
            xR = x_c(n+1)
            fL = f_c(n)
            fR = f_c(n+1)
        endif
      enddo
      
      f = fL + (x - xL)*(fR - fL)/(xR-xL)
      
      
      end subroutine linear_interp_1d
      
      end module mod_1dinterp
