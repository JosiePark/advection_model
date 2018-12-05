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
      

      
c      ic = int(nbins/2.)
      
c      if (x .lt. dfloat(ii)/2.) then
      
c        if (x .le. xc(1)) then
        
c            x_interp(ic+1:nbins) = xc(1:ic)
c            x_interp(1:ic) = xc(ic+1:nbins) - dfloat(ii)
c            K_interp(ic+1:nbins) = Kc(1:ic)
c            K_interp(1:ic) = Kc(ic+1:nbins)
            
c        else
        
c        do b = 1,nbins
        
c        if (x .le. xc(b)) then
    
        
c            x_interp(ic+1:nbins) = xc(b+1:b+ic)
c            x_interp(ic-b+1:ic) = xc(1:b) 
c            x_interp(1:ic-b) = xc(b+ic+1:nbins)-dfloat(ii)
c            K_interp(ic+1:nbins) = Kc(b+1:b+ic)
c            K_interp(ic-b+1:ic) = Kc(1:b) 
c            K_interp(1:ic-b) = Kc(b+ic+1:nbins)
            
c            exit
            
c        endif
        
      
c        enddo
        
c        endif
      
c      else
      
c      if (x.ge.xc(nbins)) then
      
c        x_interp(1:ic) = xc(nbins-(ic+1):nbins)
c        x_interp(ic+1:nbins) = xc(1:ic) + dfloat(ii)
c        K_interp(1:ic) = Kc(nbins-(ic+1):nbins)
c        K_interp(ic+1:nbins) = Kc(1:ic)
      
c      else
      
c          do b = 1,nbins
          
c            if (x .le. xc(b)) then
            
c                x_interp(1:ic) = xc(b-ic:b-1) 
c                x_interp(ic+1:ic+nbins-b+1) = xc(b:nbins)
c                x_interp(ic+2+nbins-b:nbins) = xc(1:b-ic) + dfloat(ii)
c                K_interp(1:ic) = Kc(b-ic:b-1) 
c                K_interp(ic+1:ic+nbins-b+1) = Kc(b:nbins)
c                K_interp(ic+2+nbins-b:nbins) = Kc(1:b-ic)
                
                
c                exit
   
c            endif
          
c        enddo
c      endif
      
c      endif
      
      
      
c      if (x .le. xc(1)) then
      
c        K_interp(1) = Kc(nbins)
c        K_interp(2:nbins) = Kc(1:nbins-1)
c        x_interp(1) = xc(nbins) - dfloat(ii)
c        x_interp(2:nbins) = xc(1:nbins-1)
      
c      elseif (x .ge. xc(nbins)) then
      
c        K_interp(nbins) = Kc(1)
c        K_interp(1:nbins-1) = Kc(2:nbins)
c        x_interp(nbins) = xc(1) + dfloat(ii)
c        x_interp(1:nbins-1) = xc(2:nbins)
      
c      else
        
c        K_interp = Kc
c        x_interp = xc
      
c      endif

      
      call interp_1d(x_interp,4,K_interp,x,K)
      
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

      end module mod_diffusivity_functions
