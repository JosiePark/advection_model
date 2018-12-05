c
c 2L-QG model in zonally periodic
c
c CHECK FOR DIMENSIONS IN SUBROUTINE solv_ell_mike
c Bottom topography added in beta term in lines 160-170
c Main cycle in lines 470-1000
c Hemant - Addtion of bottom topography and f is defined for that purpose
c Lines 577-600 RES_BETA include topography effects in that.


      program main
      
      
      use MOD_qg_input
      use MOD_qg2_NETCDF
      use MOD_constants
      use MOD_traj_NETCDF
      use MOD_elliptic
      use MOD_variables
      
      IMPLICIT NONE
c      INCLUDE 'fftw_f77.i'     
      include "fftw3.f" 
      

      step_tim = 1                                                      !! Keeps count for saving data
      if(istart.eq.0) then
       call create_netcdf_file(FILE_NAME,ii,jj,basinscale)              !! Create netcdf file   
      end if   
                                                                        !! else give time in days where to start  
      if(istart_ave.eq.0) then
        call create_ave_file(ave_file,ii,jj,basinscale) 
      else
        call read_ave_file(ave_file,ii,jj,psi1_av,psi2_av,time_av)
        psi1_av = psi1_av*time_av
        psi2_av = psi2_av*time_av
      endif

C/// NUMERICAL PARAMETERS FOR OUTPUT

        if(regime.eq.1) then
        visc_bot = 1.D-8
        elseif(regime.eq.2) then
        visc_bot = 25.D-9
        elseif(regime.eq.3) then
        visc_bot = 4.D-8
        elseif(regime.eq.4) then
        visc_bot = 1.D-7
        endif
         

      TIME_DAY=0
      if (istart_ave.eq.0) then
      TIME_AV=0
      endif
c
c--- PB counters
c
      
C NOTE THAT SINCE THE COUNTERS ARE NORMILISED BY THE ORIGINAL DT,
C WHICH MAY BE VERY DIFFERENT TO THE ACTUAL TIME STEP IN THE CODE,
C THE SAVING INTERVAL IN DAYS MAY NOT CORRESPOND TO K_OUT AND K_SAVE
      TIME_O=0.D0
      TIME_S=0.D0
c
c--- nondimensionalization
c
      uscale=1.
      scale=basinscale/dfloat(ii) 
      tscale=scale/uscale

      dt=dt/tscale

      fo_nondim = fo*scale/uscale
      beta_nondim=beta*scale*scale/uscale
      visc_nondim=visc/(scale*uscale)
      visc_bot_nondim=visc_bot*scale/uscale

      do i=1,ii
         do j=1,jj
            viscosity(i,j)=visc_nondim
         enddo
      enddo
      

c--- stratification 
c
      SS=(scale/Rd)**2

      S1=SS/(1.+H1/H2)
      S2=SS-S1
      write(*,*)'S1,S2=',S1,S2
c
c--- beta*y
c
      BETA_NONDIM_U1=BETA_NONDIM+U_0*S1
      BETA_NONDIM_U2=BETA_NONDIM-U_0*S2

!! ---------------- Define Bottom topography (Hemant)--------------------!!
!! h_max is max height of bottom
!! cas=0 for flat bottom; cas=1 for constant slope in x; cas=2 for constant slope in y; !!
!! cas=3 contant slopes in x,y both!!
!! More cases can be created depending upon the problem !!
!! where are q1, pv1 calculations? 
    !    cas = 1         
    !    call bottom_topography(ii,jj,H2,h_max,cas,fo_nondim,bot_topo)

c --- I don't understand this equation ---------
      do j=1,jj
         beta_y(j)=beta_nondim*(dfloat(j-1)-dfloat(jj-1)/2.D0)
         beta1_y(j)=(beta_nondim+S1*U_0)*(dfloat(j-1)-dfloat(jj-1)/2.D0)
         beta2_y(j)=(beta_nondim-S2*U_0)*(dfloat(j-1)-dfloat(jj-1)/2.D0)
      enddo
c
c--- initialization for elliptic solver
c
      theta11=H1/(H1+H2)
      theta12=H2/(H1+H2)
      theta21= 1.
      theta22=-1.

      omega11=1.
      omega12= H2/(H1+H2)
      omega21=1.
      omega22=-H1/(H1+H2)

c MIKE
        do j=1,jj
           do i=1,ii/2+1
               b1(i,j)=    2D0*(
     &         dcos(2D0*pi*(i-1)/ii)+dcos(2D0*pi*(j-1)/jj)-2D0
     &                        )
               b2(i,j)=-SS+2D0*(
     &         dcos(2D0*pi*(i-1)/ii)+dcos(2D0*pi*(j-1)/jj)-2D0
     &                        )
               if(b1(i,j).ne.0D0)then
                 b1(i,j)=1D0/b1(i,j)
                 b2(i,j)=1D0/b2(i,j)
               else
                 b1(i,j)=0D0                 
                 b2(i,j)=0D0
             endif
           enddo
        enddo
c MIKE
       call dfftw_plan_dft_r2c_2d(plan1,ii,jj,z,z11,FFTW_MEASURE)!MIKE real to complex
       call dfftw_plan_dft_c2r_2d(plan2,ii,jj,phi11,phi,FFTW_MEASURE)!MIKE complex to real
 
c
c--- INITIAL CONDITIONS
c
      if (istart_ave.eq.0) then
        do j=1,jj
         do i=1,ii
         
            psi1_av(i,j)=0.
            psi2_av(i,j)=0.
            
         enddo
        enddo
       endif
      
      do j=1,jj
         do i=1,ii

            DO K=1,2
               ZETA_OLD(K,I,J)=0
               ZETA_NEW(K,I,J)=0
               RES_VISC(K,I,J)=0
               RES_BETA(K,I,J)=0
               RES_BETA_O(K,I,J)=0

               UIM(K,I,J)=0
               UJM(K,I,J)=0
               UIMO(K,I,J)=0
               UJMO(K,I,J)=0
               UPI(K,I,J)=0
               UPJ(K,I,J)=0
            END DO
         END DO
      END DO

      DO J=1,JJ
         DO I=1,II+1
            DO K=1,2
               ZETA_FI(K,I,J)=0
               PSI_EDGE(K,I,J)=0
            END DO
         END DO
      END DO

      DO J=1,JJ+1
         DO I=1,II
            DO K=1,2
               ZETA_FJ(K,I,J)=0
                 PSI_EDGE(K,I,J)=0
            END DO
         END DO
      END DO

      DO J=1,JJ+1
         DO I=1,II
            DO K=1,2
               ZETA_FJN(K,I,J)=0
               FLUX_J(K,I,J)=0
            END DO
         END DO
      END DO

      DO J=1,JJ
         DO I=1,II+1
            DO K=1,2
               ZETA_FIN(K,I,J)=0
               FLUX_I(K,I,J)=0
            END DO
         END DO
      END DO

C DEFAULT TIME STEP ////////////

      T=DT
      DL1=DT
      DL2=DT
      TWEIGHT=1.
c
c--- initial conditions
c What is cff??
      if(istart.eq.0)then
        cff=80.
        do j=1,jj
           do i=1,ii
              psi1(i,j)=0.
              psi2(i,j)=0.

              psi1(i,j)=cff*dsin(2.*pi*10.*float(j)/float(jj)) ! set template of jets
           enddo
        enddo
        psi1(50,50)=psi1(50,50)+cff
        psi1(100,100)=psi1(100,100)-cff
        psi1(150,150)=psi1(150,150)+cff
        psi1(200,200)=psi1(200,200)-cff
        psi1(250,250)=psi1(250,250)+cff
        psi1(300,300)=psi1(300,300)-cff
        psi1(350,350)=psi1(350,350)+cff
        psi1(400,400)=psi1(400,400)-cff
        psi1(450,450)=psi1(450,450)+cff
        psi1(500,500)=psi1(500,500)-cff

        call rel_from_psi(ii,jj,psi1,psi2,rel1,rel2)
  
        do j=1,jj
           do i=1,ii
              zeta_OLD(1,i,j)=rel1(i,j)-S1*(psi1(i,j)-psi2(i,j))
              zeta_OLD(2,i,j)=rel2(i,j)-S2*(psi2(i,j)-psi1(i,j))
              zeta_NEW(1,i,j)=zeta_OLD(1,i,j)
              zeta_NEW(2,i,j)=zeta_OLD(2,i,j)
           enddo
        enddo
        
c
c----------------------
c
      elseif(istart.eq.1)then
        write(*,*)'Starting from psi-restart...'

        call read_netcdf(FILE_NAME,psi1,psi2,ii,jj,
     +  read_tim, new_tim, step_tim)
        TIME_DAY = new_tim
        step_tim = step_tim+1
        call rel_from_psi(ii,jj,psi1,psi2,rel1,rel2)

        do j=1,jj
           do i=1,ii
              zeta_OLD(1,i,j)=rel1(i,j)-S1*(psi1(i,j)-psi2(i,j))
              zeta_OLD(2,i,j)=rel2(i,j)-S2*(psi2(i,j)-psi1(i,j))
              zeta_NEW(1,i,j)=zeta_OLD(1,i,j)
              zeta_NEW(2,i,j)=zeta_OLD(2,i,j)
           enddo
        enddo
      endif

C EDGES
      DO I=1,II
         IM1=I-1
         IF(IM1.EQ.0) IM1=II
         DO J=1,JJ
            JM1=J-1
            IF(JM1.EQ.0) JM1=JJ

            PSI_EDGE(1,I,J)=.25D0*(
     &          PSI1(I,J)+PSI1(I,JM1)+PSI1(IM1,J)+PSI1(IM1,JM1)
     &                            )
            PSI_EDGE(2,I,J)=.25D0*(
     &          PSI2(I,J)+PSI2(I,JM1)+PSI2(IM1,J)+PSI2(IM1,JM1)
     &                            )
         END DO
      END DO

      DO J=1,JJ
         DO K=1,2
            PSI_EDGE(K,II+1,J)=PSI_EDGE(K,1,J)
         END DO
      END DO

        DO I=1,II
         DO K=1,2
            PSI_EDGE(K,I,JJ+1)=PSI_EDGE(K,I,1)
         END DO
      END DO

      UMAX=0
      DO K=1,2
         KVAR=1
         IF(K.EQ.2)KVAR=0
         DO J=1,JJ
            DO I=1,II
               UJ(K,I,J)= (PSI_EDGE(K,I+1,J)-PSI_EDGE(K,I,J))
               IF(ABS(UJ(K,I,J)).GT.UMAX) UMAX=ABS(UJ(K,I,J))
               UJMO(K,I,J)=UJ(K,I,J)
            END DO
         END DO

         DO J=1,JJ
            DO I=1,II
               UI(K,I,J)=-(PSI_EDGE(K,I,J+1)-PSI_EDGE(K,I,J))+U_0*KVAR
               IF(ABS(UI(K,I,J)).GT.UMAX) UMAX=ABS(UI(K,I,J))
               UIMO(K,I,J)=UI(K,I,J)
            END DO
         END DO
      END DO
c

c print*,"Just before main cycle",J
c################################# MAIN CYCLE #################################
c
 1001   CONTINUE

         goto 101
         umax=0.
         do j=1,jj
            jm1=j-1
            jp1=j+1
            if(jm1.eq.0) jm1=jj
            if(jp1.eq.jj+1) jp1=1 
            do i=1,ii
               im1=i-1
               ip1=i+1
               if(im1.eq.0) im1=ii
               if(ip1.eq.ii+1) ip1=1

               u=.5*abs(psi1(ip1,j)-psi1(im1,j))
               if(u.gt.umax) umax=u
               u=.5*abs(psi1(i,jp1)-psi1(i,jm1))
               if(u.gt.umax) umax=u
               u=.5*abs(psi2(ip1,j)-psi2(im1,j))
               if(u.gt.umax) umax=u
               u=.5*abs(psi2(i,jp1)-psi2(i,jm1))
               if(u.gt.umax) umax=u
            enddo
         enddo
         write(*,*)'umax=',umax
 101      continue

C DEFINE THE TIME STEP IN DAYS
c      write(*,*)'dt=',dt
      DT_DAY=dt*tscale/86400.
      time_day=time_day+DT_DAY
      k_day=int(time_day+0.001)
      TIME_O=TIME_O+DT_DAY
      TIME_S=TIME_S+DT_DAY

C MAX DT ALLOWABLE
      if(TIME_O.ge.TIME_OUT-.5*dt_day)then
        TIME_O=0.
      endif

      if(TIME_S.ge.TIME_SAVE-.5*dt_day)then
        TIME_S=0.
      endif

      DTMAX=1.E+6
        
      DO J=1,JJ
         DO I=1,II
            DO K=1,2
C I
               cff=CFL/ABS(UI(K,I,J))
               IF(ABS(UI(K,I,J)).GT.0)THEN
                 IF(DTMAX.GT.cff) DTMAX=cff
               ELSE
                 DTMAX=DT
               END IF
C J
               cff=CFL/ABS(UJ(K,I,J))
               IF(ABS(UJ(K,I,J)).GT.0)THEN
                 IF(DTMAX.GT.cff) DTMAX=cff
               ELSE
                 DTMAX=DT
               END IF
            END DO
         END DO
      END DO

C EXTRAPOLATION FROM THE PREVIOUS TIME LEVELS
      DT_XO=DT_O
      DT_O=DT
      DT=MIN(DTMAX,100*T)            
      DL1=DT_XO+DT_O
      DL2=DT_O+DT
      TWEIGHT=DL2/DL1
      TWEIGHTU=DT/DL2
      TWEIGHT_A=0.5*DT/DT_O
      DT05=DT*0.5
      DTVISC=DT*visc_nondim
      TWEIGHT_A1=1.+TWEIGHT_A
      TWEIGHTU1=1.+TWEIGHTU

C //////CABARET PREDICTOR STEP
      DO J=1,JJ
         DO I=1,II
            ZETA_NEW(1,I,J)=ZETA_OLD(1,I,J)+DT05*RES(1,I,J)
            ZETA_NEW(2,I,J)=ZETA_OLD(2,I,J)+DT05*RES(2,I,J)
         END DO
      END DO

C/// TIME INEGRATION OF SOURCE TERMS

      cff1=DT05*BETA_NONDIM_U1
      cff2=DT05*(BETA_NONDIM_U2+slope_y*fo_nondim*scale/H2)             !! beta and effects of bottom slope in y
      cff3=DT05*slope_x*fo_nondim*scale/H2                              !! effects of bottom slope in x
      do j=1,jj
         JP1=J+1
         IF(JP1.EQ.JJ+1) JP1=1
         do i=1,ii
            IP1=i+1
            IF(IP1.EQ.II+1) IP1=1
            RES_BETA(1,I,J)=-cff1*(UJ(1,I,JP1)+UJ(1,I,J))
            RES_BETA(2,I,J)=-cff2*(UJ(2,I,JP1)+UJ(2,I,J))
     &                      -cff3*(UI(2,IP1,J)+UI(2,I,J))               !! Extra term added for slope in x
         enddo
      enddo
C SECOND-ORDER IN TIME INTEGRATION, WHICH SUPRESSES SPURIOUS BETA_Y-WAVES;
C NB: INTERSTINGLY, THE SIMPLE OPERATOR ZETA(N+1)-ZETA(N+1/2) = L(ZETA(N+1/2)), AS SUGGESTED BY VM,
C ALSO LEADS TO A STABLE SOLUTION BUT IT IS CONTAMINATED BY SPURIOUS LARGE-SCALE BETA_Y-WAVES

      DO J=1,JJ
         DO I=1,II
            ZETA_NEW(1,I,J)=ZETA_NEW(1,I,J)
     &        +TWEIGHT_A1*RES_BETA(1,I,J)-TWEIGHT_A*RES_BETA_O(1,I,J)
            ZETA_NEW(2,I,J)=ZETA_NEW(2,I,J)
     &        +TWEIGHT_A1*RES_BETA(2,I,J)-TWEIGHT_A*RES_BETA_O(2,I,J)
         END DO
      END DO

      do j=1,jj
         do i=1,ii
            z1_ell(i,j)= theta11*ZETA_NEW(1,I,J) !barotropic mode
     &                  +theta12*ZETA_NEW(2,I,J)
            z2_ell(i,j)= theta21*ZETA_NEW(1,I,J) !baroclinic mode
     &                  +theta22*ZETA_NEW(2,I,J)
         enddo
      enddo

c      call solv_ell(ii,jj,z1_ell,phi1(1,1),0.D0,b1,d_ell)
c      call solv_ell(ii,jj,z2_ell,phi2(1,1),SS,b2,d_ell)
        phii(:,:)=0D0; !MIKE
        call solv_ell_mike(z1_ell,ii,jj,phii,b1,plan1,plan2) !MIKE
        phi1(:,:)=phii(:,:) !MIKE

        phii(:,:)=0D0; !MIKE
        call solv_ell_mike(z2_ell,ii,jj,phii,b2,plan1,plan2) !MIKE
        phi2(:,:)=phii(:,:) !MIKE
 

c
c--- convert to layers
c
      do j=1,jj
         do i=1,ii
            psi1(i,j)= omega11*phi1(i,j)
     &                +omega12*phi2(i,j)
            psi2(i,j)= omega21*phi1(i,j)
     &                +omega22*phi2(i,j)
         enddo
      enddo

      call rel_from_psi(ii,jj,psi1,psi2,rel1,rel2)

C BULK VISCOSITY

      cff=visc_bot_nondim*DT
      do j=1,jj
         jm1=j-1
         jp1=j+1
         if(jm1.eq.0) jm1=jj
         if(jp1.eq.jj+1) jp1=1
         do i=1,ii
            im1=i-1
            ip1=i+1
            if(im1.eq.0) im1=ii
            if(ip1.eq.ii+1) ip1=1

            RES_VISC(1,I,J)=DTVISC*(-4.*ZETA_NEW(1,i,j)
     &                  +ZETA_NEW(1,im1,j)+ZETA_NEW(1,ip1,j)
     &                  +ZETA_NEW(1,i,jm1)+ZETA_NEW(1,i,jp1)
     &              +S1*(REL1(I,J)-REL2(I,J))
     &                             )
    
            RES_VISC(2,I,J)=DTVISC*(-4.*ZETA_NEW(2,i,j)
     &                  +ZETA_NEW(2,im1,j)+ZETA_NEW(2,ip1,j)
     &                  +ZETA_NEW(2,i,jm1)+ZETA_NEW(2,i,jp1)
     &              +S2*(REL2(I,J)-REL1(I,J))
     &                             )
     &               -cff*rel2(i,j)
          
         enddo
      enddo
        
      DO J=1,JJ
         DO I=1,II
            ZETA_NEW(1,I,J)=ZETA_NEW(1,I,J)+RES_VISC(1,I,J)
            ZETA_NEW(2,I,J)=ZETA_NEW(2,I,J)+RES_VISC(2,I,J)
         END DO
      END DO  

C END OF VISCOSITY

C ///UPDATING CELL-FACE VARIABLES  /////

C FIRST DEFINE EDGE VARIABLES FOR PSI
C PERIODIC BCS IN I: EDGE(1)<->EDGE(II+1),CENTRE (0)<-> CENTRE (II)

      DO J=1,JJ
         JM1=J-1
         IF(JM1.EQ.0) JM1=JJ
         DO I=1,II
            IM1=I-1
            IF(IM1.EQ.0) IM1=II
            IP1=I+1
            IF(IP1.EQ.II+1) IP1=1

            UJM(1,I,J)= .25D0*(
     &           PSI1(IP1,J)+PSI1(IP1,JM1)-PSI1(IM1,J)-PSI1(IM1,JM1)
     &                        )
            UJM(2,I,J)= .25D0*(
     &           PSI2(IP1,J)+PSI2(IP1,JM1)-PSI2(IM1,J)-PSI2(IM1,JM1)
     &                        )
         END DO
      END DO

      DO J=1,JJ
         JM1=J-1
         IF(JM1.EQ.0) JM1=JJ
         JP1=J+1
         IF(JP1.EQ.JJ+1) JP1=1
         DO I=1,II
            IM1=I-1
            IF(IM1.EQ.0) IM1=II

            UIM(1,I,J)=-.25D0*(
     &           PSI1(I,JP1)-PSI1(I,JM1)+PSI1(IM1,JP1)-PSI1(IM1,JM1)
     &                        )+U_0
            UIM(2,I,J)=-.25D0*(
     &           PSI2(I,JP1)-PSI2(I,JM1)+PSI2(IM1,JP1)-PSI2(IM1,JM1)
     &                        )
         END DO
      END DO

      UMAX=0
      DO J=1,JJ
         DO I=1,II
            UJ(1,I,J)=TWEIGHTU1*UJM(1,I,J)-TWEIGHTU*UJMO(1,I,J)      
            IF(ABS(UJM(1,I,J)).GT.UMAX) UMAX=ABS(UJM(1,I,J))
            UJ(2,I,J)=TWEIGHTU1*UJM(2,I,J)-TWEIGHTU*UJMO(2,I,J)    
            IF(ABS(UJM(2,I,J)).GT.UMAX) UMAX=ABS(UJM(2,I,J))
         END DO
      END DO

      DO J=1,JJ
         DO I=1,II
            UI(1,I,J)=TWEIGHTU1*UIM(1,I,J)-TWEIGHTU*UIMO(1,I,J)
            IF(ABS(UIM(1,I,J)).GT.UMAX) UMAX=ABS(UIM(1,I,J))
            UI(2,I,J)=TWEIGHTU1*UIM(2,I,J)-TWEIGHTU*UIMO(2,I,J)
            IF(ABS(UIM(2,I,J)).GT.UMAX) UMAX=ABS(UIM(2,I,J))
         END DO
      END DO


C CALCULATE FLUX VARIABLES OF THE CABARET SCHEME
      
C I DIRECTION
       DO J=1,JJ
          DO I=1,II
             DO K=1,2
                FI1=ZETA_FI(K,I,J)
                IF(UI(K,I,J).GE.0) THEN
C POSITIVE SPEED
                  IF(I.GT.1) THEN
                    FI2=ZETA_FI(K,I-1,J)
                    V0=UI(K,I,J)+UI(K,I-1,J)
                    PSI_N=ZETA_NEW(K,I-1,J)
                    PSI_O=ZETA_OLD(K,I-1,J)
                    D_F=-(PSI_N-PSI_O)*2.-DT05*V0*(FI1-FI2)
                  ELSE
C BCS
                    V0=UI(K,I,J)+UI(K,II,J)
                    FI2=ZETA_NEW(K,II,J)
                    PSI_N=ZETA_NEW(K,II,J)
                    PSI_O=ZETA_OLD(K,II,J)
                    D_F=-(PSI_N-PSI_O)*2.-DT05*V0*(FI1-FI2)
                  END IF
                ELSE
C NEGATIVE SPEED
                  IF(I.LT.II) THEN
                    FI2=ZETA_FI(K,I+1,J)
                    PSI_N=ZETA_NEW(K,I,J)
                    PSI_O=ZETA_OLD(K,I,J)
                    V0=UI(K,I,J)+UI(K,I+1,J)
                    D_F=-(PSI_N-PSI_O)*2.-DT05*V0*(FI2-FI1)
                  ELSE
C BCS
                    FI2=ZETA_FI(K,1,J)
                    PSI_N=ZETA_NEW(K,I,J)
                    PSI_O=ZETA_OLD(K,I,J)
                    V0=UI(K,I,J)+UI(K,1,J)
                    D_F=-(PSI_N-PSI_O)*2.-DT05*V0*(FI2-FI1)
                  END IF
                END IF

C NON-LINEAR CORRECTION TO ENFORCE THE MAXIMUM PRINCIPLE
                FMAX=MAX(FI1,FI2,PSI_N)-D_F
                FMIN=MIN(FI1,FI2,PSI_N)-D_F
                  
                ZETA_FIN(K,I,J)=2.*PSI_N-FI2    
              
                IF(ZETA_FIN(K,I,J).GT.FMAX)ZETA_FIN(K,I,J)=FMAX        
                IF(ZETA_FIN(K,I,J).LT.FMIN)ZETA_FIN(K,I,J)=FMIN
             END DO
C PERIODIC BCS
          END DO
          ZETA_FIN(1,II+1,J)=ZETA_FIN(1,1,J)
          ZETA_FIN(2,II+1,J)=ZETA_FIN(2,1,J)
       END DO

C J DIRECTION
       DO J=1,JJ
          DO I=1,II
             DO K=1,2
                FI1=ZETA_FJ(K,I,J)
                IF(UJ(K,I,J).GE.0) THEN
C POSITIVE SPEED
                  IF(J.GT.1) THEN
                    FI2=ZETA_FJ(K,I,J-1)
                    PSI_N=ZETA_NEW(K,I,J-1)
                    PSI_O=ZETA_OLD(K,I,J-1)
                    V0=UJ(K,I,J)+UJ(K,I,J-1)
                    D_F=-(PSI_N-PSI_O)*2.-DT05*V0*(FI1-FI2)
                  ELSE
C BCS
                    FI2=ZETA_FJ(K,I,JJ)
                    PSI_N=ZETA_NEW(K,I,JJ)
                    PSI_O=ZETA_OLD(K,I,JJ)
                    V0=UJ(K,I,J)+UJ(K,I,JJ)
                    D_F=-(PSI_N-PSI_O)*2.-DT05*V0*(FI1-FI2)

                  END IF    
                ELSE
C NEGATIVE SPEED
                  IF(J.LT.JJ) THEN
                    FI2=ZETA_FJ(K,I,J+1)
                    PSI_N=ZETA_NEW(K,I,J)
                    PSI_O=ZETA_OLD(K,I,J)
                    V0=UJ(K,I,J)+UJ(K,I,J+1)
                    D_F=-(PSI_N-PSI_O)*2.-DT05*V0*(FI2-FI1)
                  ELSE
C BCS
                    FI2=ZETA_FJ(K,I,1)
                    PSI_N=ZETA_NEW(K,I,J)
                    PSI_O=ZETA_OLD(K,I,J)
                    V0=UJ(K,I,J)+UJ(K,I,1)
                    D_F=-(PSI_N-PSI_O)*2.-DT05*V0*(FI2-FI1)

                  END IF
                END IF
              
C NON-LINEAR CORRECTION TO ENFORCE THE MAXIMUM PRINCIPLE
                FMAX=MAX(FI1,FI2,PSI_N)-D_F
                FMIN=MIN(FI1,FI2,PSI_N)-D_F  
                  
                ZETA_FJN(K,I,J)=2.*PSI_N-FI2  

                IF(ZETA_FJN(K,I,J).GT.FMAX)ZETA_FJN(K,I,J)=FMAX        
                IF(ZETA_FJN(K,I,J).LT.FMIN)ZETA_FJN(K,I,J)=FMIN
             END DO
          END DO
       END DO
C PERIODIC BCS
          DO I=1,II
          ZETA_FJN(1,I,JJ+1)=ZETA_FJN(1,I,1)
          ZETA_FJN(2,I,JJ+1)=ZETA_FJN(2,I,1)
          END DO
C UPDATING CELL FACE VARIABLES

       DO J=1,JJ
          DO I=1,II+1
             ZETA_FI(1,I,J)=ZETA_FIN(1,I,J)
             ZETA_FI(2,I,J)=ZETA_FIN(2,I,J)
          END DO
       END DO

       DO J=1,JJ+1
          DO I=1,II
             ZETA_FJ(1,I,J)=ZETA_FJN(1,I,J)
             ZETA_FJ(2,I,J)=ZETA_FJN(2,I,J)
          END DO
       END DO

C COMPUTE FLUXES

       DO J=1,JJ
          DO I=1,II
             FLUX_I(1,I,J)=UI(1,I,J)*ZETA_FI(1,I,J)
             FLUX_I(2,I,J)=UI(2,I,J)*ZETA_FI(2,I,J)
          END DO
          FLUX_I(1,II+1,J)=FLUX_I(1,1,J)
          FLUX_I(2,II+1,J)=FLUX_I(2,1,J)
       END DO

       DO J=1,JJ
          DO I=1,II
             FLUX_J(1,I,J)=UJ(1,I,J)*ZETA_FJ(1,I,J)
             FLUX_J(2,I,J)=UJ(2,I,J)*ZETA_FJ(2,I,J)
          END DO
       END DO
C PERIODIC BCS
       DO I=1,II
          FLUX_J(1,I,JJ+1)=FLUX_J(1,I,1)
          FLUX_J(2,I,JJ+1)=FLUX_J(2,I,1)
       END DO

C// END OF THE CABARET SOLVER; START POSTPROCESSING /////////////////////        
c
c--- time-mean streamfunction
c
        IF(TIME_DAY.GT.MAX_SPIN) THEN
          TIME_AV=TIME_AV+DT
          !print*,'time_av = ',time_av
          !print*,'dt = ,'dt

          do j=1,jj
             do i=1,ii
                psi1_av(i,j)=psi1_av(i,j)+DT*psi1(i,j)
                psi2_av(i,j)=psi2_av(i,j)+DT*psi2(i,j)
             enddo
          enddo
        do j=1,jj
         do i=1,ii
            psi1_av(i,j)=psi1_av(i,j)/TIME_AV
            psi2_av(i,j)=psi2_av(i,j)/TIME_AV
         enddo
      enddo

      call write_ave_file(ave_file,ii,jj,psi1_av,psi2_av,time_av)
      
      psi1_av = psi1_av*time_av
      psi2_av = psi2_av*time_av
      
       !print*,'averaged data written at',time_av

      
      endif
c
c--- diagnostics
c
        if(TIME_O.eq.0)then
           write(*,*)'Time (days) =',TIME_DAY

           call energy(H1,H2,S1,S2,ii,jj,psi1,psi2,ekin1,ekin2,epot)

            etot = epot+ekin1+ekin2
c            write(*,*)
           write(*,'(A5,F12.6,A8,F12.6,A8,F12.6,A8,F12.6)')
     & 'epot=',epot,'; ekin1=',ekin1,'; ekin2=',ekin2,'; Total=',etot
        endif
c
c--- save output
c
        if(TIME_S.eq.0)then
        
c            call zeta_from_psi_and_rel(ii,jj,S1,S2,psi1,psi2,rel1,rel2
c     &                                                    ,z1,z2)  
                                                           
          call write_netcdf(FILE_NAME,psi1,psi2,
     +    epot,ekin1,ekin2,ii,jj,TIME_DAY,step_tim)
          step_tim = step_tim + 1
          
            

        endif

C // UPDATING CENTRE-CELL VARIABLES/////

        DO J=1,JJ
           DO I=1,II
              RES(1,I,J)=-(FLUX_I(1,I+1,J)-FLUX_I(1,I,J)
     &                    +FLUX_J(1,I,J+1)-FLUX_J(1,I,J))
              RES(2,I,J)=-(FLUX_I(2,I+1,J)-FLUX_I(2,I,J)
     &                    +FLUX_J(2,I,J+1)-FLUX_J(2,I,J))
           END DO
        END DO
      
        DO J=1,JJ
           DO I=1,II
              ZETA_NEW(1,I,J)=ZETA_NEW(1,I,J)+DT05*RES(1,I,J)
              ZETA_NEW(2,I,J)=ZETA_NEW(2,I,J)+DT05*RES(2,I,J)
           END DO
        END DO

C UPDATING ZETA_OLD

        DO J=1,JJ
           DO I=1,II
              ZETA_OLD(1,I,J)=ZETA_NEW(1,I,J)
              UIMO(1,I,J)=UIM(1,I,J)
              UJMO(1,I,J)=UJM(1,I,J)

              ZETA_OLD(2,I,J)=ZETA_NEW(2,I,J)
              UIMO(2,I,J)=UIM(2,I,J)
              UJMO(2,I,J)=UJM(2,I,J)

              RES_BETA_O(1,I,J)=RES_BETA(1,I,J)
              RES_BETA_O(2,I,J)=RES_BETA(2,I,J)
           END DO
        END DO

      IF(TIME_DAY.LT.MAX_TIME) GOTO 1001
c
c#################################################
      call dfftw_destroy_plan(plan1) !MIKE
      call dfftw_destroy_plan(plan2) !MIKE

      
      
      
      

      stop
      end program
C
C============================================================================





