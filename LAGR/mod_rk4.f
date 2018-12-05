c code that contains subroutine for RK4 method for use with both 
c 2D cubic and bicubic spatial interpolation method

      module mod_rk4
      
      use mod_2dcubic
      use mod_bicubic
      
      contains
      
c-----------------------------------------------------------------
c-----------------------------------------------------------------      
      
      subroutine rk4_2dcubic(ii,jj,x,y,dt,U0,a_old,b_old,c_old
     & ,d_old,a_half,b_half,c_half,d_half,a_new,b_new,c_new,d_new
     & ,x_diff,y_diff)
c INPUT:     
c x, y are the locations of the lagrangian particles
c U0 is the uniform background flow (set to 0 if neglecting)
c ii,jj are the number of grid points in x and y
c a,b,c,d are the coefficients obtained by calling cubic_coeff_x
c old,half,new refer to which time the coefficients are obtained from
c dt is the time step
c OUTPUT :
c x_diff and y_diff are the resulting displacements
      implicit none
      
      integer ii,jj
      real*8 x,y,u_point0,v_point0,u_point1,v_point1,u_point2,v_point2
     & ,u_point3,v_point3,U0,dt,dt05,dt6,xt,yt,x_diff,y_diff
      real*8 a_old(ii,jj),b_old(ii,jj),c_old(ii,jj),d_old(ii,jj)
     & ,a_half(ii,jj),b_half(ii,jj),c_half(ii,jj),d_half(ii,jj)
     & ,a_new(ii,jj),b_new(ii,jj),c_new(ii,jj),d_new(ii,jj) ,psi_x(4) 
      dt05 = dt/2.
      dt6 = dt/6.
      
c ------------- FIRST STEP ----------------------------

      call cubic_poly_x(ii,jj,x,y
     & ,a_old,b_old,c_old,d_old,psi_x)

      call vel(ii,jj,psi_x
     & ,a_old,b_old,c_old,d_old
     & ,x,y,u_point0,v_point0)
     
            xt = x + dt05*(u_point0+U0)
            yt = y + dt05*(v_point0)
        
            if(xt<= 0) then
            xt = xt + dfloat(ii)
            elseif(xt>dfloat(ii)) then
            xt = xt - dfloat(ii)
            endif
            if(yt<=0) then
            yt=yt+dfloat(jj)
            elseif(yt>dfloat(jj)) then
            yt=yt-dfloat(jj)
            endif
 
     
c ------------- SECOND STEP ---------------------------
      call cubic_poly_x(ii,jj,xt,yt
     & ,a_half,b_half,c_half,d_half,psi_x)

      call vel(ii,jj,psi_x
     & ,a_half,b_half,c_half,d_half
     & ,xt,yt,u_point1,v_point1)
     
     
            xt = x + dt05*(u_point1+U0)
            yt = y + dt05*(v_point1)
            
            
            if(xt<= 0) then
            xt = xt + dfloat(ii)
            elseif(xt>dfloat(ii)) then
            xt = xt - dfloat(ii)
            endif
            if(yt<=0) then
            yt=yt+dfloat(jj)
            elseif(yt>dfloat(jj)) then
            yt=yt-dfloat(jj)
            endif
            
c ------------- THIRD STEP ---------------------------
      call cubic_poly_x(ii,jj,xt,yt
     & ,a_half,b_half,c_half,d_half,psi_x)

      call vel(ii,jj,psi_x
     & ,a_half,b_half,c_half,d_half
     & ,xt,yt,u_point2,v_point2)
     
            
            xt = x + dt*(u_point2+U0)
            yt = y + dt*(v_point2)
            


            
            if(xt<= 0) then
            xt = xt + dfloat(ii)
            elseif(xt>dfloat(ii)) then
            xt = xt - dfloat(ii)
            endif
            if(yt<=0) then
            yt=yt+dfloat(jj)
            elseif(yt>dfloat(jj)) then
            yt=yt-dfloat(jj)
            endif
c ------------- FOURTH STEP -------------------------

      call cubic_poly_x(ii,jj,xt,yt
     & ,a_new,b_new,c_new,d_new,psi_x)

      call vel(ii,jj,psi_x
     & ,a_new,b_new,c_new,d_new
     & ,xt,yt,u_point3,v_point3)
     
            x_diff = dt*U0 + dt6*(u_point0 + 2*(u_point1+u_point2)
     & + u_point3)
            y_diff = dt6*(v_point0 + 2*(v_point1+v_point2)
     & +v_point3)

	    if(y_diff .gt. 20) then
         print*,'velocities = ',v_point0,v_point1,v_point2,v_point3
            endif
 
     
      end subroutine
      
c-----------------------------------------------------------------
c-----------------------------------------------------------------
      
      subroutine rk4_bicubic(ii,jj,x,y,dt,U0,M_old,M_half,M_new
     & ,x_diff,y_diff)
c INPUT:     
c x, y are the locations of the lagrangian particles
c U0 is the uniform background flow (set  = 0 if neglecting)
c ii,jj are the number of grid points in x and y
c M is the matrix of coefficients obtained by calling A_matrix
c old,half,new refer to which time the coefficients are obtained from
c dt is the time step
c OUTPUT :
c x_diff and y_diff are the resulting displacements
      implicit none
      
      integer ii,jj
      real*8 x,y,dt,U0,x_diff,y_diff,dt05,dt6,u_point0,v_point0,u_point1
     & ,v_point1,u_point2,v_point2
     & ,u_point3,v_point3
     & ,xt,yt
      real*8 M_old(ii,jj,4,4),M_half(ii,jj,4,4),M_new(ii,jj,4,4)
      
      dt05 = dt/2.
      dt6 = dt/6.
      
c ---------------- FIRST STEP -------------------------------

      call bicubic(ii,jj,M_old,x,y,u_point0,v_point0)
      
     
            xt = x + dt05*(u_point0+U0)
            yt = y + dt05*(v_point0)
            
      
            if(xt<= 0) then
            xt = xt + dfloat(ii)
            elseif(xt>dfloat(ii)) then
            xt = xt - dfloat(ii)
            endif
            if(yt<=0) then
            yt=yt+dfloat(jj)
            elseif(yt>dfloat(jj)) then
            yt=yt-dfloat(jj)
            endif

c ---------------- SECOND STEP -------------------------------

      call bicubic(ii,jj,M_half,xt,yt,u_point1,v_point1)
      
            xt = x + dt05*(u_point1+U0)
            yt = y + dt05*(v_point1)
            
            if(xt<= 0) then
            xt = xt + dfloat(ii)
            elseif(xt>dfloat(ii)) then
            xt = xt - dfloat(ii)
            endif
            if(yt<=0) then
            yt=yt+dfloat(jj)
            elseif(yt>dfloat(jj)) then
            yt=yt-dfloat(jj)
            endif
c ---------------- THIRD STEP -----------------------------
      call bicubic(ii,jj,M_half,xt,yt,u_point2,v_point2)
      

            xt = x + dt*(u_point2+U0)
            yt = y + dt*(v_point2)
            
            if(xt<= 0) then
            xt = xt + dfloat(ii)
            elseif(xt>dfloat(ii)) then
            xt = xt - dfloat(ii)
            endif
            if(yt<=0) then
            yt=yt+dfloat(jj)
            elseif(yt>dfloat(jj)) then
            yt=yt-dfloat(jj)
            endif
c ---------------- FOURTH STEP -----------------------------
      call bicubic(ii,jj,M_new,xt,yt,u_point3,v_point3)
      
            x_diff = dt*U0 + dt6*(u_point0 + 2*(u_point1+u_point2)
     & + u_point3)
            y_diff = dt6*(v_point0 + 2*(v_point1+v_point2)
     & +v_point3)
 
      
      end subroutine
      
      end module mod_rk4
