      module MOD_variables
      
      
      contains
      
            subroutine rel_from_psi(ii,jj,phi1,phi2,rel1,rel2)
      implicit none
      integer jj,ii,j,i,im1,ip1,jm1,jp1
      real*8 phi1(ii,jj),phi2(ii,jj),rel1(ii,jj),rel2(ii,jj)

      do j=1,jj
         do i=1,ii
            im1=i-1
            ip1=i+1
            jm1=j-1
            jp1=j+1
            if(im1.eq.0)    im1=ii
            if(ip1.eq.ii+1) ip1=1
            if(jm1.eq.0)    jm1=jj
            if(jp1.eq.jj+1) jp1=1

            rel1(i,j)=-4.*phi1(i,j)
     &                        +phi1(im1,j)+phi1(i,jp1)
     &                        +phi1(ip1,j)+phi1(i,jm1)
            rel2(i,j)=-4.*phi2(i,j)
     &                        +phi2(im1,j)+phi2(i,jp1)
     &                        +phi2(ip1,j)+phi2(i,jm1)
         enddo
      enddo

      return
      end subroutine
c
c----------------------------------------------------------------------
c
      subroutine zeta_from_psi_and_rel(ii,jj,S1,S2,psi1,psi2,rel1,rel2
     &                                                    ,zeta1,zeta2)
      implicit none
      integer ii,jj,i,j
      real*8 psi1(ii,jj),psi2(ii,jj),rel1(ii,jj),rel2(ii,jj)
     & ,zeta1(ii,jj),zeta2(ii,jj)
     & ,S1,S2

      do j=1,jj
         do i=1,ii
            zeta1(i,j)=rel1(i,j)
     &                 -S1*(psi1(i,j)-psi2(i,j))
            zeta2(i,j)=rel2(i,j)
     &                 -S2*(psi2(i,j)-psi1(i,j))
         enddo
      enddo

      return
      end subroutine
c
c----------------------------------------------------------------------
c
      subroutine energy(H1,H2,S1,S2,ii,jj,psi1,psi2,ekin1,ekin2,epot)
      implicit none
      integer ii,jj,i,j,im1,ip1,jm1,jp1
      real*8 psi1(ii,jj),psi2(ii,jj)
     & ,H1,H2,S1,S2,ekin1,ekin2,epot,H,fac1,fac2,cff

      H=H1+H2
      ekin1=0.
      ekin2=0.
      epot=0.

      do j=1,jj
         do i=1,ii
            epot=epot+(psi1(i,j)-psi2(i,j))**2
         enddo
      enddo

      do j=1,jj
         do i=1,ii
            im1=i-1
            ip1=i+1
            jm1=j-1
            jp1=j+1
            if(im1.eq.0)    im1=ii
            if(ip1.eq.ii+1) ip1=1
            if(jm1.eq.0)    jm1=jj
            if(jp1.eq.jj+1) jp1=1

           ekin1=ekin1+.25*(
     &   (psi1(im1,j)+psi1(i,j)-psi1(im1,jm1)-psi1(i,jm1))**2           !! What is KE expression
     &  +(psi1(i,jm1)+psi1(i,j)-psi1(im1,jm1)-psi1(im1,j))**2
     &                     )
           ekin2=ekin2+.25*(
     &   (psi2(im1,j)+psi2(i,j)-psi2(im1,jm1)-psi2(i,jm1))**2
     &  +(psi2(i,jm1)+psi2(i,j)-psi2(im1,jm1)-psi2(im1,j))**2
     &                     )
         enddo
      enddo
      cff=0.5/dfloat(ii*jj)
      ekin1=cff*(H1/H)*ekin1
      ekin2=cff*(H2/H)*ekin2
      epot =0.5*cff*((S1*H1+S2*H2)/H)*epot                              ! How to come up this expression for PE    

      return
      end subroutine
      
      subroutine vel_from_psi(ii,jj,psi1,psi2,u1,v1,u2,v2)
      
      implicit none
      
      integer i,j,ii,jj,im1,ip1,jm1,jp1
      real*8 psi1(ii,jj), psi2(ii,jj),u1(ii,jj),u2(ii,jj),v1(ii,jj),
     & v2(ii,jj)      
      
      
      
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
            
            u1(i,j) = -.5*(psi1(i,jp1)-psi1(i,jm1))
            u2(i,j) = -.5*(psi2(i,jp1)-psi2(i,jm1))
            v1(i,j) = .5*(psi1(ip1,j)-psi1(im1,j))
            v2(i,j) = .5*(psi2(ip1,j)-psi2(im1,j))
        
        
        enddo
      enddo
      
      end subroutine vel_from_psi
c -----------------------------------------------------------
c -----------------------------------------------------------
      
      subroutine PVgrad_from_zeta(ii,jj,zeta,PVgrad)

      implicit none
      
      integer ii,jj,i,j,jm1,jp1
      
      real*8 zeta(ii,jj), PVgrad(ii,jj)
      
       do i = 1,ii
        do j = 1,jj

            jm1=j-1
            jp1=j+1

            if(jm1.eq.0)    jm1=jj
            if(jp1.eq.jj+1) jp1=1
            
            PVgrad(i,j) = zeta(i,jp1) - zeta(i,jm1)
            
        enddo
      enddo
      
      end subroutine
      
c ---------------------------------------------------------------
c -----------  FIND ZONALLY AND TEMPORALY AVERAGED FULL PV -----------

      subroutine PV_bar_from_psi(ii,jj,jj1,basinscale,beta
     & ,Rd,H1,H2,U_0,psi1_av,psi2_av,PV_bar)
   
      implicit none
      
      integer ii,jj,coord_range,jj1
      real*8 basinscale,beta,Rd,H1,H2,U_0
      
      real*8 psi1_av(ii,jj),psi2_av(ii,jj),PV_bar(jj1)
      
      real*8 uscale,scale,tscale,SS,S1,S2,beta_nondim,beta_nondim_u1
      
      integer i,j,kk,idx
      
      real*8 rel1_av(ii,jj), rel2_av(ii,jj), zeta1_av(ii,jj)
     & , zeta2_av(ii,jj),beta1_y(jj1)
     & , PV_av(ii,jj1)
     
      parameter(coord_range = 3)
      
      print*,'running PV_bar_from_psi'
      
c set non-dimensional variables

      
      uscale=1.
      scale=basinscale/dfloat(ii) 
      tscale=scale/uscale
      SS=(scale/Rd)**2

      S1=SS/(1.+H1/H2)
      S2=SS-S1
      beta_nondim=beta*scale*scale/uscale
      
      BETA_NONDIM_U1=BETA_NONDIM+U_0*S1
      
      
      
c CALCULATE FULL ZONALLY AVERAGED PV

      call rel_from_psi(ii,jj,psi1_av,psi2_av,rel1_av,rel2_av) ! calculate time-averaged relative vorticity
      print*,'relative vorticity done'
      call zeta_from_psi_and_rel(ii,jj,S1,S2,psi1_av,psi2_av,rel1_av
     & ,rel2_av,zeta1_av,zeta2_av) ! time averaged PV anomaly
       print *,'PV anomaly done'
c CALCULATE ZONAL FULL TIME-AVERAGED PV

      do j = 1,jj1
            beta1_y(j) = beta_nondim_u1*(dfloat(j-1)) ! only do for the top layer
      enddo
      print*,'beta done'
      do j=1,jj
      do kk = 1,coord_range
            idx = j+(kk-1)*jj
            PV_bar(idx) = 0.
            do i = 1,ii
                PV_av(i,idx) = zeta1_av(i,j) + beta1_y(idx)
                PV_bar(idx) = PV_bar(idx) + PV_av(i,idx)
            enddo
            PV_bar(idx) = PV_bar(idx)/dfloat(ii)
      enddo
      enddo
      
      end subroutine PV_bar_from_psi
      
      
      
      end module MOD_variables
      
