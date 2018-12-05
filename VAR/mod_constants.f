      module MOD_constants
      
      use MOD_qg_input
      
      implicit none
      
c COUNTERS AND ITERATIVES
      
      integer i,j,i1,i2
     & ,k_day,ngpswk,K,IP,KVAR,IM1,JM1,IP1,JP1,cas,step_tim
      
      
C TIME VARIABLES
  
      real*8 UMAX,A1,A2,A3,DL1,DL2,TWEIGHT,TWEIGHT_A
     &,TWEIGHTU,DT_O,DT_XO,BETA_NONDIM_U1,BETA_NONDIM_U2
     &,DT05,DTVISC,TWEIGHT_A1,TWEIGHTU1,T
     &,TIME_AV,TIME_O,TIME_S,u,dt_day 
     
C QG VARIABLES
     
       real *8 new_tim
       real*8 psi1(ii,jj),psi2(ii,jj),phi1(ii,jj),phi2(ii,jj)
     & ,rel1(ii,jj),rel2(ii,jj),z1(ii,jj),z2(ii,jj),bot_topo(ii,jj)
     & ,PVgrad(ii)

c     & ,b1(ii,jj),b2(ii,jj)
     & ,b1(ii/2+1,jj),b2(ii/2+1,jj)!MIKE
     & ,z1_ell(ii,jj),z2_ell(ii,jj)

     & ,beta_y(jj),beta1_y(jj),beta2_y(jj)
     & ,psi1_av(ii,jj),psi2_av(ii,jj)

     & ,pv1(0:ii+1,jj),pv2(0:ii+1,jj)
     & ,q1(0:ii+1,jj),q2(0:ii+1,jj)
     & ,v1(0:ii+1,jj),v2(0:ii+1,jj)

     & ,S1,S2,SS
     & ,theta11,theta12,theta21,theta22
     & ,omega11,omega12,omega21,omega22
     & ,uscale,scale,tscale,cff,cff1,cff2, cff3
     & ,beta_nondim,fo_nondim,visc_nondim,visc_bot_nondim
     & ,time_day,ekin,ekin1,ekin2,epot,etot

      real aux(ii,jj),aux1(ii,jj),aux2(ii,jj),aux3(ii,jj)
     &               ,aux4(ii,jj),aux5(ii,jj),aux6(ii,jj)

     & ,FI1,FI2,U1,U2,PSI_N,PSI_O,DTMAX,V0,D_F,FMIN,FMAX,PSI_M,PSI_P
     & ,FI_M,FI_P

C CABARET VARIABLES /////////////////////////////
C CELL CENTRES
     & ,ZETA_NEW(2,II,JJ),ZETA_OLD(2,II,JJ)
C CELL FACES    
     & ,ZETA_FI(2,II+1,JJ),ZETA_FJ(2,II,JJ+1)
     & ,ZETA_FIN(2,II+1,JJ),ZETA_FJN(2,II,JJ+1)
     & ,UI(2,II,JJ),UJ(2,II,JJ)
     & ,UPI(2,II,JJ),UPJ(2,II,JJ)
     &,UIM(2,II,JJ),UJM(2,II,JJ),UIMO(2,II,JJ),UJMO(2,II,JJ)
     &,FLUX_I(2,II+1,JJ),FLUX_J(2,II,JJ+1)

C CELL EDGES
     &,PSI_EDGE(2,II+1,JJ+1)      
C CELL RESIDUALS
     &, RES(2,II,JJ),RES_VISC(2,II,JJ)
     &, RES_BETA(2,II,JJ),RES_BETA_O(2,II,JJ)

     &,viscosity(ii,jj)
     
      real*8 pi,fo,beta
     & ,slope_x,slope_y,dt

      complex*16 d_ell(jj)

      INTEGER*8 plan1,plan2
      REAL*8 z(ii,jj),phi(ii,jj),phii(ii,jj)
      DOUBLE COMPLEX z11(ii/2+1,jj),phi11(ii/2+1,jj)
      
      
      common /ONE/ psi1,psi2
      common /TWO/ z1,z2,phi1,phi2
      common /THREE/ rel1,rel2
      common /FOUR/ theta11,theta12,theta21,theta22
      common /FIVE/ omega11,omega12,omega21,omega22
      common /EIGHT/ b1,b2,d_ell,z1_ell,z2_ell
      common /NINE/ beta_y,beta1_y,beta2_y,psi1_av,psi2_av
      common /TEN/ pv1,pv2,q1,q2,v1,v2
      common /ELEVEN/ aux,aux1,aux2,aux3,aux4,aux5,aux6
      
      data dt/2160./

      data pi/3.14159265358979323846D0/
      data fo/0.83D-4/
      data beta/2.D-13/
      data slope_x/0.D-5/                                               !! Assign bottom slopes                 
      data slope_y/0.D-5/   

      end module MOD_constants
