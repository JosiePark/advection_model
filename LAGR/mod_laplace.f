c module containing the subroutine necessary for cubic interpolation that doesn't require non- divergence
c it requires the velocity calculated using finite differences
      module mod_laplace
      implicit none
      
      contains
      
              subroutine cubic(ii,jj,psi,x,y
     &    ,u_point,v_point)
     
      implicit none
      
      integer ii,jj,i,j,ic,jc,xc(4),yc(4),im1,ip1,jm1,jp1
      real*8 u(ii,jj),v(ii,jj),x,y,u_point,v_point,ax,ay,ax2,ay2
     & ,ax3,ay3,cff(4),aux1(4),aux2(4),psi(ii,jj)
      real*8 c1,c2
      parameter(c1=1.D0/6.D0,c2=2.D0/6.D0)
      
      ic = int(x) + 1
      jc = int(y) + 1
      !print*,'ic,jc=',ic,jc
      
      ax = x - (ic-1)
      ay = y - (jc-1)
      
      ! periodic boundary conditions
      
      do i = 1,4
        xc(i) = ic +i-2
        yc(i) = jc+i-2
      enddo
      !print*,'xc=',xc
      !print*,'yc=',yc
      
      do i = 1,4
        if (xc(i) .gt. ii) then
          xc(i) = xc(i) - ii

        end if
        if(xc(i) .eq. 0) then
          xc(i) = ii
          endif
        if (yc(i) .gt.jj) then

          yc(i) = yc(i)-jj
          endif
        if (yc(i).eq.0)then
          yc(i) = jj

          endif
      enddo
      
c calculate the velocity using finite differences

       do i = 1,ii
       do j = 1,jj
            im1=i-1
            ip1=i+1
            jm1=j-1
            jp1=j+1
            if(im1.eq.0)    im1=ii
            if(ip1.eq.ii+1) ip1=1
            if(jm1.eq.0)    jm1=jj
            if(jp1.eq.jj+1) jp1=1
            
            u(i,j) = .5*(psi(i,jp1)-psi(i,jm1))

            v(i,j) = -.5*(psi(ip1,j)-psi(im1,j))
            
       enddo
       enddo


            
      
      !print*,'xc,yc=',xc,yc
    
      ax2 = ax**2
      ay2 = ay**2
      ax3 = ax**3
      ay3=ay**3

      cff(1)=  -c2*ay+.5*ay2-c1*ay3
      cff(2)=1.-.5*ay   -ay2+.5*ay3
      cff(3)=      ay+.5*ay2-.5*ay3
      cff(4)=  -c1*ay       +c1*ay3
      
        do j=1,4
         aux1(j)=0.
         aux2(j)=0.
         do i=1,4
            aux1(j)=aux1(j)+cff(i)*u(xc(i),yc(j))
            aux2(j)=aux2(j)+cff(i)*v(yc(i),xc(j))
         enddo
      enddo
      
      cff(1)=  -c2*ax+.5*ax2-c1*ax3
      cff(2)=1.-.5*ax   -ax2+.5*ax3
      cff(3)=      ax+.5*ax2-.5*ax3
      cff(4)=  -c1*ax       +c1*ax3

      u_point= cff(1)*aux1(1)+cff(2)*aux1(2)
     &       +cff(3)*aux1(3)+cff(4)*aux1(4)
      v_point= cff(1)*aux2(1)+cff(2)*aux2(2)
     &       +cff(3)*aux2(3)+cff(4)*aux2(4)
      
      
    
      return
      end subroutine
      
      end module
