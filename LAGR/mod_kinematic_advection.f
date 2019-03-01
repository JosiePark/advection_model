      module mod_kinematic_advection
      
      contains
      
      subroutine kinematic_velocity(x,y,t,a,f,c,k,y_centre,u,v)
      
      implicit none
      
      real*8 x,y,t
      real*8 u,v
      real*8 a,f,c,y_centre,k
      
      u = 2*A*f*cos(k*(c-k*c*t))
     & *tanh(f*(y-y_centre))/(cosh(f*(y-y_centre))**2)
      v = k*A/(cosh(f*(y-y_centre))**2)*sin(k*(x-k*c*t))
      
      end subroutine kinematic_velocity
      
c --------------------------------------------------------------------

      subroutine rk4_kinematic(x,y,t,a,f,c,k,y_centre,dt,x_diff,y_diff)
      
      implicit none
      
      real*8 x,y,xt,yt,x_diff,y_diff
      real*8 t,dt,dt2,a,f,c,y_centre,k
      real*8 u1,v1,u2,v2,u3,v3,u4,v4
      
      dt2 = dt/2.
      
      
c perform rk4
      
      call kinematic_velocity(x,y,t,a,f,c,k,y_centre,u1,v1)
      
      xt = x + dt2*u1
      yt = y + dt2*v1
      
      call kinematic_velocity(xt,yt,t+dt2,a,f,c,k,y_centre,u2,v2)
      
      xt = x + dt2*u2
      yt = y + dt2*v2
      
      call kinematic_velocity(xt,yt,t+dt2,a,f,c,k,y_centre,u3,v3)
      
      xt = x + dt*u3
      yt = y + dt*v3
      
      call kinematic_velocity(xt,yt,t+dt,a,f,c,k,y_centre,u4,v4)
      
      x_diff = (dt/6.)*(u1+2*(u2+u3)+u4)
      y_diff = (dt/6.)*(v1+2*(v2+v3)+v4) 
      
      end subroutine rk4_kinematic
      
      end module mod_kinematic_advection
