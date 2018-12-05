      module MOD_elliptic
      
      !use MOD_constants
      
      contains
      
      
      
      
      
              SUBROUTINE solv_ell_mike(z,ii,jj,phi,b,plan1,plan2)
        IMPLICIT NONE
c        INCLUDE 'fftw_f77.i'  
        INCLUDE 'fftw3.f'

        INTEGER i,j,ii,jj,i1
c        PARAMETER(ii=512,jj=512)
c        PARAMETER(ii=512,jj=ii)
c        PARAMETER(ii=1024,jj=ii)
c        PARAMETER(ii=2048,jj=ii)
c        PARAMETER(ii=4096,jj=ii)
        INTEGER*8 plan1,plan2
        REAL*8 z(ii,jj),phi(ii,jj)
        REAL*8 b(ii/2+1,jj)
        DOUBLE COMPLEX, DIMENSION(:,:),ALLOCATABLE:: z1,phi1,b1
        
c        z1(int(ii/2)+1,jj)
c     & ,phi1(int(ii/2)+1,jj),b1(int(ii/2)+1,jj)
!        common /SOLV_ELL/ z1,phi1,b1
        
        i1 = ii/2+1
        
        ALLOCATE(z1(i1,jj),b1(i1,jj),phi1(i1,jj))


        call dfftw_execute_dft_r2c(plan1,z,z1)

        b1=CMPLX(b)
        phi1(:,:)=b1(:,:)*z1(:,:)

        call dfftw_execute_dft_c2r(plan2,phi1,phi)

        do j=1,jj
           do i=1,ii
              phi(i,j)=phi(i,j)/DBLE(ii*jj)
           enddo
        enddo


        return
        END 
        
c              subroutine rlft3(data,speq,nn1,nn2,nn3,isign)
c      implicit none
c      integer isign,nn1,nn2,nn3
c      complex*16 speq(nn2,nn3),data(nn1/2,nn2,nn3)
cCU    USES fourn
c      integer i1,i2,i3,j1,j2,j3,nn(3)
c      real*8 theta,wi,wpi,wpr,wr,wtemp
c      complex*16 c1,c2,h1,h2,w

c      c1=cmplx(0.5,0.0)
c      c2=cmplx(0.0,-0.5*isign)
c      theta=6.28318530717959d0/dble(isign*nn1)
c      wpr=-2.0d0*sin(0.5d0*theta)**2
c      wpi=sin(theta)
c      nn(1)=nn1/2
c      nn(2)=nn2
c      nn(3)=nn3
c      if(isign.eq.1)then
c        call fourn(data,nn,3,isign)
c        do 12 i3=1,nn3
c          do 11 i2=1,nn2
c            speq(i2,i3)=data(1,i2,i3)
c11        continue
c12      continue
c      endif
c      do 15 i3=1,nn3
c        j3=1
c        if (i3.ne.1) j3=nn3-i3+2
c        wr=1.0d0
c        wi=0.0d0
c        do 14 i1=1,nn1/4+1
c          j1=nn1/2-i1+2
c          do 13 i2=1,nn2
c            j2=1
c            if (i2.ne.1) j2=nn2-i2+2
c            if(i1.eq.1)then
c              h1=c1*(data(1,i2,i3)+conjg(speq(j2,j3)))
c              h2=c2*(data(1,i2,i3)-conjg(speq(j2,j3)))
c              data(1,i2,i3)=h1+h2
c              speq(j2,j3)=conjg(h1-h2)
c            else
c              h1=c1*(data(i1,i2,i3)+conjg(data(j1,j2,j3)))
c              h2=c2*(data(i1,i2,i3)-conjg(data(j1,j2,j3)))
c              data(i1,i2,i3)=h1+w*h2
c              data(j1,j2,j3)=conjg(h1-w*h2)
c            endif
c13        continue
c          wtemp=wr
c          wr=wr*wpr-wi*wpi+wr
c          wi=wi*wpr+wtemp*wpi+wi
c          w=cmplx(sngl(wr),sngl(wi))
c14      continue
c15    continue
c      if(isign.eq.-1)then
c        call fourn(data,nn,3,isign)
c      endif

c      return
c      end
cc
cc-----------------------------------------------------------------------------
cc
c      subroutine fourn(data,nn,ndim,isign)
c      implicit none
c      integer isign,ndim,nn(ndim)
c      real*8 data(*)
c      integer i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,ip3,k1,
c     &k2,n,nprev,nrem,ntot
c      real*8 tempi,tempr,theta,wi,wpi,wpr,wr,wtemp

c      ntot=1
c      do 11 idim=1,ndim
c        ntot=ntot*nn(idim)
c11    continue
c      nprev=1
c      do 18 idim=1,ndim
c        n=nn(idim)
c        nrem=ntot/(n*nprev)
c        ip1=2*nprev
c        ip2=ip1*n
c        ip3=ip2*nrem
c        i2rev=1
c        do 14 i2=1,ip2,ip1
c          if(i2.lt.i2rev)then
c            do 13 i1=i2,i2+ip1-2,2
c              do 12 i3=i1,ip3,ip2
c                i3rev=i2rev+i3-i2
c                tempr=data(i3)
c                tempi=data(i3+1)
c                data(i3)=data(i3rev)
c                data(i3+1)=data(i3rev+1)
c                data(i3rev)=tempr
c                data(i3rev+1)=tempi
c12            continue
c13          continue
c          endif
c          ibit=ip2/2
c1         if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
c            i2rev=i2rev-ibit
c            ibit=ibit/2
c          goto 1
c          endif
c          i2rev=i2rev+ibit
c14      continue
c        ifp1=ip1
c2       if(ifp1.lt.ip2)then
c          ifp2=2*ifp1
c          theta=isign*6.28318530717959d0/(ifp2/ip1)
c          wpr=-2.d0*sin(0.5d0*theta)**2
c          wpi=sin(theta)
c          wr=1.d0
c          wi=0.d0
c          do 17 i3=1,ifp1,ip1
c            do 16 i1=i3,i3+ip1-2,2
c              do 15 i2=i1,ip3,ifp2
c                k1=i2
c                k2=k1+ifp1
c                tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
c                tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
c                data(k2)=data(k1)-tempr
c                data(k2+1)=data(k1+1)-tempi
c                data(k1)=data(k1)+tempr
c                data(k1+1)=data(k1+1)+tempi
c15            continue
c16          continue
c            wtemp=wr
c            wr=wr*wpr-wi*wpi+wr
c            wi=wi*wpr+wtemp*wpi+wi
c17        continue
c          ifp1=ifp2
c        goto 2
c        endif
c        nprev=n*nprev
c18    continue

c      return
c      end
cc
cc-----------------------------------------------------------------------------
cc
c      subroutine realft(data,n,isign)
c        implicit real*8 (a-h,o-z)
c      INTEGER isign,n
c      real*8 data(n)
c      INTEGER i,i1,i2,i3,i4,n2p3

c      theta=3.141592653589793d0/dble(n/2)
c      c1=0.5
c      if (isign.eq.1) then
c        c2=-0.5
c        call four1(data,n/2,+1)
c      else
c        c2=0.5
c        theta=-theta
c      endif
c      wpr=-2.0d0*dsin(0.5d0*theta)**2
c      wpi=dsin(theta)
c      wr=1.0d0+wpr
c      wi=wpi
c      n2p3=n+3
c      do 11 i=2,n/4
c        i1=2*i-1
c        i2=i1+1
c        i3=n2p3-i2
c        i4=i3+1
c        wrs=wr
c        wis=wi
c        h1r=c1*(data(i1)+data(i3))
c        h1i=c1*(data(i2)-data(i4))
c        h2r=-c2*(data(i2)+data(i4))
c        h2i=c2*(data(i1)-data(i3))
c        data(i1)=h1r+wrs*h2r-wis*h2i
c        data(i2)=h1i+wrs*h2i+wis*h2r
c        data(i3)=h1r-wrs*h2r+wis*h2i
c        data(i4)=-h1i+wrs*h2i+wis*h2r
c        wtemp=wr
c        wr=wr*wpr-wi*wpi+wr
c        wi=wi*wpr+wtemp*wpi+wi
c11    continue
c      if (isign.eq.1) then
c        h1r=data(1)
c        data(1)=h1r+data(2)
c        data(2)=h1r-data(2)
c      else
c        h1r=data(1)
c        data(1)=c1*(h1r+data(2))
c        data(2)=c1*(h1r-data(2))
c        call four1(data,n/2,-1)
c      endif

c      return
c      end
cc
cc-----------------------------------------------------------------------------
cc
c      subroutine four1(data,nn,isign)
c        implicit real*8 (a-h,o-z)
c      INTEGER isign,nn
c      real*8 data(2*nn)
c      INTEGER i,istep,j,m,mmax,n
c      n=2*nn
c      j=1
c      do 11 i=1,n,2
c        if(j.gt.i)then
c          tempr=data(j)
c          tempi=data(j+1)
c          data(j)=data(i)
c          data(j+1)=data(i+1)
c          data(i)=tempr
c          data(i+1)=tempi
c        endif
c        m=n/2
c1       if ((m.ge.2).and.(j.gt.m)) then
c          j=j-m
c          m=m/2
c        goto 1
c        endif
c        j=j+m
c11    continue
c      mmax=2
c2     if (n.gt.mmax) then
c        istep=2*mmax
c        theta=6.28318530717959d0/(isign*mmax)
c        wpr=-2.d0*dsin(0.5d0*theta)**2
c        wpi=dsin(theta)
c        wr=1.d0
c        wi=0.d0
c        do 13 m=1,mmax,2
c          do 12 i=m,n,istep
c            j=i+mmax
c            tempr=wr*data(j)-wi*data(j+1)
c            tempi=wr*data(j+1)+wi*data(j)
c            data(j)=data(i)-tempr
c            data(j+1)=data(i+1)-tempi
c            data(i)=data(i)+tempr
c            data(i+1)=data(i+1)+tempi
c12        continue
c          wtemp=wr
c          wr=wr*wpr-wi*wpi+wr
c          wi=wi*wpr+wtemp*wpi+wi
c13      continue
c        mmax=istep
c      goto 2
c      endif

c      return
c      end
        
        end module MOD_elliptic
