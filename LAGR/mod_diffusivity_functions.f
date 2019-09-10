      module mod_diffusivity_functions
      
      use mod_1dinterp
      
      implicit none
    
      contains
      
      subroutine diffusivity_interp_1d(xc,nbins,Kc,ii,x,K)
      
      ! EMPLOYS THE FACT THAT THE DIFFUSIVITY IS DOUBLY PERIODIC TO INTERPOLATE
      ! OTHERIWSE INTERPOLATION SCHEME FAILS WHEN THE POINT IS OUT OF RANGE OF THE BIN CENTRES.
      ! NEED TO ADD CONDITION THAT ENSURES THE DIFFUSIVITY IS NEVER NEGATIVE
      
      ! INPUT : xc, bins centres
      ! INPUT : nbins, number of bins
      ! INPUT : Kc, diffusivity of the bins
      ! INPUT : ii, grid size
      ! INPUT : x, location at which you wish to interpolate
      ! OUTPUT : K, diffusivity at x
      
      implicit none
      
      integer nbins,ii,ic,b
      real*8 xc(nbins),Kc(nbins),x,K,K_interp(4),x_interp(4)
      
      ! account for doubly periodicity
      
      
      if (x.ge.xc(nbins)) then
        x_interp(1:2) = xc(nbins-1:nbins)
        x_interp(3:4) = xc(1:2) + dfloat(ii)
        K_interp(1:2) = Kc(nbins-1:nbins)
        K_interp(3:4) = Kc(1:2)
      else
      
      do b = 1,nbins
        if (x.le.xc(b)) then
            if (b.eq.1) then
            x_interp(1:2) = xc(nbins-1:nbins) - dfloat(ii)
            x_interp(3:4) = xc(1:2)
            K_interp(1:2) = Kc(nbins-1:nbins)
            K_interp(3:4) = Kc(1:2)
            elseif (b .eq. 2) then
            x_interp(1) = xc(nbins) - dfloat(ii)
            x_interp(2:4) = xc(1:3)
            K_interp(1) = Kc(nbins)
            K_interp(2:4) = Kc(1:3)
            elseif (b .eq. nbins) then
            x_interp(1:3) = xc(nbins-2:nbins)
            x_interp(4) = xc(1) + dfloat(ii)
            K_interp(1:3) = Kc(nbins-2:nbins)
            K_interp(4) = Kc(1)
            else
            x_interp = xc(b-2:b+1)
            K_interp = Kc(b-2:b+1)
            endif
        exit
        endif
      enddo
      
      endif
      
      call interp_1d(x_interp,4,K_interp,x,K)
c      print*,'x,k=',x,K
      
      if (K .le. 1.D-3) then
        K =  1.D-3
      endif
      
      end subroutine diffusivity_interp_1d
      
      subroutine diffusivity_derivative(xc,nbins,Kc,ii,x,dKdx)
      
      ! input: Kc is the diffusivity defined at the bins
      ! input: xc is the location of the bin centres
      ! input: nbins is the number of bins
      ! input: x is the point at which you wish to calculate the derivate
      ! ouput: dKdx is the derivate of K at x 
      
      implicit none
      
      integer nbins,ii
      real*8 Kc(nbins),xc(nbins),x,dKdx
      
      real*8 h,x0,K0,K
      
      h = .1
      
      x0 = x+h
      
      call diffusivity_interp_1d(xc,nbins,Kc,ii,x0,K0)
      call diffusivity_interp_1d(xc,nbins,Kc,ii,x,K)
      
      dKdx = (K0-K)/h
      
    
      end subroutine diffusivity_derivative
      
      subroutine parameter_interp_1d(xc,nbins,Kc,ii,x,K)
      
      ! EMPLOYS THE FACT THAT THE DIFFUSIVITY IS DOUBLY PERIODIC TO INTERPOLATE
      ! OTHERIWSE INTERPOLATION SCHEME FAILS WHEN THE POINT IS OUT OF RANGE OF THE BIN CENTRES.
      ! NEED TO ADD CONDITION THAT ENSURES THE DIFFUSIVITY IS NEVER NEGATIVE
      
      ! INPUT : xc, bins centres
      ! INPUT : nbins, number of bins
      ! INPUT : Kc, diffusivity of the bins
      ! INPUT : ii, grid size
      ! INPUT : x, location at which you wish to interpolate
      ! OUTPUT : K, diffusivity at x
      
      implicit none
      
      integer nbins,ii,ic,b
      real*8 xc(nbins),Kc(nbins),x,K,K_interp(4),x_interp(4)
      
      ! account for doubly periodicity
      
      
      if (x.ge.xc(nbins)) then
        x_interp(1:2) = xc(nbins-1:nbins)
        x_interp(3:4) = xc(1:2) + dfloat(ii)
        K_interp(1:2) = Kc(nbins-1:nbins)
        K_interp(3:4) = Kc(1:2)
      else
      
      do b = 1,nbins
        if (x.le.xc(b)) then
            if (b.eq.1) then
            x_interp(1:2) = xc(nbins-1:nbins) - dfloat(ii)
            x_interp(3:4) = xc(1:2)
            K_interp(1:2) = Kc(nbins-1:nbins)
            K_interp(3:4) = Kc(1:2)
            elseif (b .eq. 2) then
            x_interp(1) = xc(nbins) - dfloat(ii)
            x_interp(2:4) = xc(1:3)
            K_interp(1) = Kc(nbins)
            K_interp(2:4) = Kc(1:3)
            elseif (b .eq. nbins) then
            x_interp(1:3) = xc(nbins-2:nbins)
            x_interp(4) = xc(1) + dfloat(ii)
            K_interp(1:3) = Kc(nbins-2:nbins)
            K_interp(4) = Kc(1)
            else
            x_interp = xc(b-2:b+1)
            K_interp = Kc(b-2:b+1)
            endif
        exit
        endif
      enddo
      
      endif
      
      call interp_1d(x_interp,4,K_interp,x,K)
c      print*,'x,k=',x,K
      
      end subroutine parameter_interp_1d
      
      subroutine parameter_derivative(xc,nbins,Kc,ii,x,dKdx)
      
      ! input: Kc is the diffusivity defined at the bins
      ! input: xc is the location of the bin centres
      ! input: nbins is the number of bins
      ! input: x is the point at which you wish to calculate the derivate
      ! ouput: dKdx is the derivate of K at x 
      
      implicit none
      
      integer nbins,ii
      real*8 Kc(nbins),xc(nbins),x,dKdx
      
      real*8 h,x0,K0,K
      
      h = .1
      
      x0 = x+h
      
      call diffusivity_interp_1d(xc,nbins,Kc,ii,x0,K0)
      call diffusivity_interp_1d(xc,nbins,Kc,ii,x,K)
      
      dKdx = (K0-K)/h
      
    
      end subroutine parameter_derivative

      end module mod_diffusivity_functions
