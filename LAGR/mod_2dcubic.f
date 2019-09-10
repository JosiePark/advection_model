c Module that contains the subroutines required for spatial 2D cubic interpolation

      module MOD_2Dcubic
      
      contains
      
c --------------------------------------------------------------------

c calculates 1D cubic polynomial coefficients across the ii x jj grid
c each set of cubic polynomial coefficients a,b,c,d are calculated for
c each grid coordinate of y. (i.e. we get ii x jj 1D cubic polynomials defined in the y-direction)   
      
      subroutine cubic_coeff_x(ii,jj,psi,a,b,c,d)
      
c INPUT : ii,jj : grid size
c INPUT : psi(ii,jj) : SNAPSHOT OF THE STREAM FUNCTION
C OUTPUT : a(ii,jj),b(ii,jj),c(ii,jj),d(ii,jj) : polynomial coefficients

      
      
      implicit none
      
      integer ii,jj,i,j,n,ic(4)
      real*8 psi(ii,jj)
     &  ,a(ii,jj),b(ii,jj)
     & ,c(ii,jj),d(ii,jj)
      
      
      
        do i = 1,ii
        ! i denotes the grid point to the left of the interpolation point in x
        do j = 1,jj
        ! j denotes the grid point above the interpolation point in y
        
        ! determine the 4 data points used to calculate the cubic polynomial
        
            do n = 1,4
                ic(n) = i - 2 + n
		        if (ic(n) <= 0) then
		        ic(n) = ii + ic(n)
		        elseif (ic(n) > ii) then
		        ic(n) = ic(n) - ii
		        endif
            enddo

        
            ! constructing coefficients for the 4 1-D polynomials defined on the grid lines 
            ! j =jc(1), jc(2), jc(3), jc(4)
        
                a(i,j) = psi(ic(2),j)
                b(i,j)= -psi(ic(1),j)/3 - psi(ic(2),j)/2 
     &            +psi(ic(3),j) - psi(ic(4),j)/6
                c(i,j) = psi(ic(1),j)/2 - psi(ic(2),j)
     &            + psi(ic(3),j)/2
                d(i,j) = -psi(ic(1),j)/6 + psi(ic(2),j)/2
     &            -psi(ic(3),j)/2 + psi(ic(4),j)/6
            
            
            
        enddo
        enddo
        end subroutine
        
c ---------------------------------------------------------------------------------------
c ---------------------------------------------------------------------------------------

c subroutine that takes the grid ii x jj, the snap shot streamfunction and the x coordinate of the
c interpolation point and returns the cubic polynomials approximated along x across the whole
c y direction in all possible locations of the y coordinate of the interpolation point
c takes a,b,c,d from the previous subroutine as input (it avoids having to continuously recalculate the coefficients)

        subroutine cubic_poly_x(ii,jj,x,y,a,b,c,d,psi_interp_x)
        
c INPUT : ii,jj : grid size
c INPUT : x,y : the interpolation point values
c INPUT : a(ii,jj),b(ii,jj),c(ii,jj),d(ii,jj) : cubic polynomical coefficients
c                                               as calculated using cubic_coeff_x
c OUTPUT : psi_interp_x(4) : returns the evaluated polynomial at the 4 surrounding y coordinates
        
        implicit none
        
        integer ii,jj,i,j,xc,jc(4),yc
        real*8 x, psi_interp_x(4),ax,y
        real*8 a(ii,jj),b(ii,jj)
     & ,c(ii,jj),d(ii,jj)
        
        xc = int(x) + 1 ! determines the grid point to the left of the interpolation point

        ax = x - int(x)
        
        yc = int(y) + 1
        
        do i = 1,4
            jc(i) = yc - 2 + i
            if (jc(i)>jj) then
            jc(i) = jc(i) - jj
            endif
            if(jc(i)<=0) then
            jc(i) = jc(i) + jj
            endif

        enddo
        
        !print*, xc,jc
        
        
        do i = 1,4
            psi_interp_x(i) = a(xc,jc(i)) + b(xc,jc(i))*ax 
     &       + c(xc,jc(i))*ax**2
     & + d(xc,jc(i))*ax**3
        enddo
        
        return
        
        end subroutine
        
        
c ----------------------------------------------------------------------------------------
c ----------------------------------------------------------------------------------------

c subroutine that takes the grid size, psi_x from previous subroutine and the x coordinate of the
c interpolation point and constructs the coefficients for the cubic polynomial
c that approximates the streamfunction on the line x.
        
        subroutine cubic_coeff_y(ii,jj,psi_x,alpha,beta,gamma,delta)
        
c INPUT: ii,jj : grid size
c INPUT: psi_x(4) : output of cubic_poly_x
c OUTPUT : alpha,beta,gamma,delta : cubic polynomial coefficients of the final
c                                   2D cubic polynomial
        
        implicit none
        
        integer ii,jj,i,j,n,jc(4)
        real*8 psi_x(4)
     &        ,alpha,beta,gamma,delta
        
        
            
            alpha = psi_x(2)
            
            beta = - psi_x(1)/3 - psi_x(2)/2 + psi_x(3)
     &            -psi_x(4)/6
     
            gamma = psi_x(1)/2 - psi_x(2) + psi_x(3)/2
            
            delta = - psi_x(1)/6 + psi_x(2)/2
     &             - psi_x(3)/2 + psi_x(4)/6
     
        
        
        return
        
        end subroutine
        
c --------------------------------------------------------------------------------------
c --------------------------------------------------------------------------------------

c subroutine that returns the value of psi approximated at the interpolation point (x,y) called psi_interp

        subroutine cubic_interp(ii,jj,psi_x,y,psi_interp)

c INPUT : ii,jj : grid size
C INPUT : psi_x(4) : output of cubic_poly_x
c INPUT : y : value of interpolation point y
c OUTPUT : psi_interp : psi(x,y) interpolated at x,y
        
        implicit none
        
        integer ii,jj,yc
        real*8 psi_x(4),x,y,psi_interp,ay
        real*8 alpha,beta,gamma,delta
        
        ay = y - int(y)
        yc = int(y) + 1
        
        call cubic_coeff_y(ii,jj,psi_x,alpha,beta,gamma,delta)
        
        psi_interp = alpha + beta*ay + gamma*ay**2 
     &    + delta*ay**3
     
        return
        
        end subroutine
        
c -----------------------------------------------------------------------------------------------
c -----------------------------------------------------------------------------------------------

c subroutine that takes the cubic coefficients a,b,c,d to calculte the interpolated value at x,y
c can be used for any variable defined on the grid ii x jj

        subroutine cubic_interp_full(ii,jj,a,b,c,d,x,y,f)
        
c INPUT: ii,jj : grid size
c INPUT : a(ii,jj),b(ii,jj),c(ii,jj),d(ii,jj) : 1d cubic polynomial coefficient
c                                             : output of cubic_coeff_x
c INPUT : x,y : interpolation point
C ONPUT: f : f = f(x,y), f evaluated at the interpolation point  
        
        implicit none
        
        integer ii,jj
        real*8 a(ii,jj),b(ii,jj),c(ii,jj),d(ii,jj)
        real*8 x,y
        real*8 f,f_x(4)
        
        call cubic_poly_x(ii,jj,x,y,a,b,c,d,f_x)
        call cubic_interp(ii,jj,f_x,y,f)
        
        
        
        end subroutine cubic_interp_full
c -----------------------------------------------------------------------------------------------
c -----------------------------------------------------------------------------------------------

c subroutine that returns the velocities given the coefficients

        subroutine vel(ii,jj,psi_x,a,b,c,d,x,y,u,v) !add a,b,c,d as input when incorporating in final code
        
        implicit none
        
        integer ii,jj,yc,xc,jc(4),i
        real*8 x,y,u,v,psi_x(4),ay,ax
        real*8 a(ii,jj),b(ii,jj)
     &        ,c(ii,jj),d(ii,jj)
        real*8 alpha,beta,gamma,delta
        real*8 d_alpha,d_beta,d_gamma,d_delta
        
        !call cubic_coeff_x(ii,jj,psi,a,b,c,d) ! put outside subroutine i.e in main function at beginning
        
        !call cubic_poly_x(ii,jj,psi,x,a,b,c,d,psi_x)
        
        
        call cubic_coeff_y(ii,jj,psi_x,alpha,beta,gamma,delta) 
        

        
        yc = int(y) + 1
        ay = y - int(y)
  
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
        

        
        
        u = -(beta + 2*gamma*ay + 3*delta*ay**2)
        
        d_alpha = b(xc,jc(2)) + 2*c(xc,jc(2))*ax + 3*d(xc,jc(2))*ax**2
        
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

     
        v = (d_alpha + ay*d_beta + d_gamma*ay**2 + d_delta*ay**3)
        
        return
        
        end subroutine
      
c --------------------------------------------------------------------

c subroutine that calculates cubic polynomial coefficients based on interpolation point location
c using the snapshot of the Potential Vorticity
c returns ii x jj + 6 variables that stores the coefficients for each block for 
c the PV defined along constant y values
c 6 extra grid points are need in the y-direction as the potential vorticity
c is not doubly periodic   

c takes realtive vorticity as input, and calculate full pv using pv  
      
      subroutine pv_cubic_coeff_x(ii,jj,zeta,beta,a,b,c,d)
      
      
      implicit none
      
      integer ii,jj,i,j,n,ic(4),jn
      real*8 pv(ii,-3:jj+3),zeta(ii,jj)
     &  ,a(ii,-3:jj+3),b(ii,-3:jj+3)
     & ,c(ii,-3:jj+3),d(ii,-3:jj+3)
     & ,beta,beta_y(-3:jj+3)
     
C CALCULATE BETA_Y ACROSS THE DOMAIN
      do j =-3,jj+3
        beta_y(j) = beta*dfloat(j-1+jj)
      enddo
C CALCULATE PV ACROSS THE DOMAIN
     
      do j = -3,jj+3
      ! apply doubly periodicity to zeta
        if (j<=0) then
            jn = j+jj
        elseif (j>jj) then
            jn = j-jj
        else
            jn = j
        endif
         
      do i = 1,ii
        PV(i,j) = zeta(i,jn) + beta_y(j)
      enddo
      enddo
      !print*,'PV(1,jj:jj+3) =',PV(1,jj:jj+3)
C CONSTRUCT COEFFICIENTS
      
      
      
        do i = 1,ii
        ! i denotes the grid point to the left of the interpolation point in x
        do j = -3,jj+3
        ! j denotes the grid point above the interpolation point in y
        
        ! determine the 4 data points used to calculate the cubic polynomial
        
            do n = 1,4
                ic(n) = i - 2 + n
		        if (ic(n) <= 0) then
		        ic(n) = ii + ic(n)
		        elseif (ic(n) > ii) then
		        ic(n) = ic(n) - ii
		        endif
            enddo
            
            ! constructing coefficients for the 4 1-D polynomials defined on the grid lines 
            ! j =jc(1), jc(2), jc(3), jc(4)
        
                a(i,j) = PV(ic(2),j)
                b(i,j)= -PV(ic(1),j)/3 - PV(ic(2),j)/2 
     &            +PV(ic(3),j) - PV(ic(4),j)/6
                c(i,j) = PV(ic(1),j)/2 - PV(ic(2),j)
     &            + PV(ic(3),j)/2
                d(i,j) = -PV(ic(1),j)/6 + PV(ic(2),j)/2
     &            -PV(ic(3),j)/2 + PV(ic(4),j)/6
            
            
            
        enddo
        enddo
        
        return
        end subroutine
        
c ---------------------------------------------------------------------------------------
c ---------------------------------------------------------------------------------------

c subroutine that takes the grid ii x jj, the snap shot PV and the x coordinate of the
c interpolation point and returns the cubic polynomials approximated along x across the whole
c y direction in all possible locations of the y coordinate of the interpolation point
c takes a,b,c,d from the previous subroutine as input (it avoids having to continuously recalculate the coefficients)

        subroutine PV_cubic_poly_x(ii,jj,x,y,a,b,c,d,PV_x)
        
        implicit none
        
        integer ii,jj,i,j,xc,jc(4),yc
        real*8 x, PV_x(4),ax,y
        real*8 a(ii,-3:jj+3),b(ii,-3:jj+3)
     & ,c(ii,-3:jj+3),d(ii,-3:jj+3)
        
        xc = int(x) + 1 ! determines the grid point to the left of the interpolation point

        ax = x - int(x)
        
        yc = int(y) + 1
        
        do i = 1,4
            jc(i) = yc - 2 + i
            !if (jc(i)>jj) then
            !jc(i) = jc(i) - jj
            !endif
            !if(jc(i)<=0) then
            !jc(i) = jc(i) + jj
            !endif

        enddo
        
        !print*, xc,jc
        
        
        do i = 1,4
            PV_x(i) = a(xc,jc(i)) + b(xc,jc(i))*ax 
     &       + c(xc,jc(i))*ax**2
     & + d(xc,jc(i))*ax**3
        enddo
        
        end subroutine
        
      end module MOD_2Dcubic
