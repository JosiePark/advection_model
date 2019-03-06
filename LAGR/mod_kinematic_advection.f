      module mod_kinematic_advection
      
      contains
      
      subroutine kinematic_velocity(ii,jj,x,y
     & ,pi,t,a,f,cnd,k,y_centre,u,v)
      
      implicit none
      
      real*8 x,y,t,xnd,time_day
      real*8 u,v
      real*8 a,f,cnd,y_centre,k
      real*8 coeff1,coeff2,coeff3,coeff4
      integer ii,jj
      real*8 pi
      
      xnd = x*2*pi/ii - pi
      time_day = t/86400.
      
      coeff1 = cos(k*(xnd-cnd*time_day))
      coeff2 = tanh(f*(y-y_centre)/jj)
      coeff3 = cosh(f*(y-y_centre)/jj)**2
      coeff4 = sin(k*(xnd-cnd*time_day))
      
      !print*,'coeff1,coeff2,coeff3,coeff4 =',coeff1,coeff2,coeff3,coeff4
      u = 2*a*f*coeff1*coeff2/(coeff3*jj)
      v = k*a*2*pi*coeff4/(coeff3*ii)
      
      !print*,'u,v=',u,v

      
      end subroutine kinematic_velocity
      
c --------------------------------------------------------------------

      subroutine rk4_kinematic(ii,jj,x,y,pi,t,a,f,c,k,y_centre,dt,dt_nd
     & ,x_diff,y_diff)
      
      implicit none
      
      real*8 x,y,xt,yt,x_diff,y_diff
      real*8 t,dt,dt2,a,f,c,y_centre,k
      real*8 u1,v1,u2,v2,u3,v3,u4,v4
      real*8 dt_day,dt_nd
      integer ii,jj
      real*8 pi
      real*8 scale
    
      dt2 = dt/2.
      scale = dfloat(ii)/520.d5
      
      
c perform rk4
      
      call kinematic_velocity(ii,jj,x,y,pi,t
     & ,a,f,c,k,y_centre,u1,v1)
      !print*,'u1,v1=',u1,v1
      
      xt = x + dt2*u1*scale
      yt = y + dt2*v1*scale
      
      call kinematic_velocity(ii,jj,xt,yt,pi,t+dt/2.
     & ,a,f,c,k,y_centre,u2,v2)
      !print*,'u2,v2=',u2,v2
      
      xt = x + dt2*u2*scale
      yt = y + dt2*v2*scale
      
      call kinematic_velocity(ii,jj,xt,yt,pi,t+dt/2.
     & ,a,f,c,k,y_centre,u3,v3)
      !print*,'u3,v3=',u3,v3
      
      xt = x + dt_nd*u3*scale
      yt = y + dt_nd*v3*scale
      
      call kinematic_velocity(ii,jj,xt,yt,pi,t+dt
     & ,a,f,c,k,y_centre,u4,v4)
      !print*,'u4,v4=',u4,v4
      
      x_diff = (dt/(6.))*(u1+2*(u2+u3)+u4)*scale
      y_diff = (dt/(6.))*(v1+2*(v2+v3)+v4)*scale 
      
      end subroutine rk4_kinematic
      
      end module mod_kinematic_advection
