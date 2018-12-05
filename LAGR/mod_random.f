c module that contains a subroutine to generate random particles

      module MOD_random
      
      ! A module for random number generation from the following distributions:
      !
      !     Distribution                    Function/subroutine name
      !
      !     Normal (Gaussian)               random_normal
      !     Gamma                           random_gamma
      !     Chi-squared                     random_chisq
      !     Exponential                     random_exponential
      !     Weibull                         random_Weibull
      !     Beta                            random_beta
      !     t                               random_t
      !     Multivariate normal             random_mvnorm
      !     Generalized inverse Gaussian    random_inv_gauss
      !     Poisson                         random_Poisson
      !     Binomial                        random_binomial1   *
      !                                     random_binomial2   *
      !     Negative binomial               random_neg_binomial
      !     von Mises                       random_von_Mises
      !     Cauchy                          random_Cauchy
      !
      !  Generate a random ordering of the integers 1 .. N
      !                                     random_order
      !     Initialize (seed) the uniform random number generator for ANY compiler
      !                                     seed_random_number

      !  ** Two functions are provided for the binomial distribution.
      !  If the parameter values remain constant, it is recommended that the
      !  first function is used (random_binomial1).   If one or both of the
      !  parameters change, use the second function (random_binomial2).

      ! The compilers own random number generator, SUBROUTINE RANDOM_NUMBER(r),
      ! is used to provide a source of uniformly distributed random numbers.

      ! N.B. At this stage, only one random number is generated at each call to
      !      one of the functions above.

      ! The module uses the following functions which are included here:
      ! bin_prob to calculate a single binomial probability
      ! lngamma  to calculate the logarithm to base e of the gamma function

      ! Some of the code is adapted from Dagpunar's book:
      !     Dagpunar, J. 'Principles of random variate generation'
      !     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
      !
      ! In most of Dagpunar's routines, there is a test to see whether the value
      ! of one or two floating-point parameters has changed since the last call.
      ! These tests have been replaced by using a logical variable FIRST.
      ! This should be set to .TRUE. on the first call using new values of the
      ! parameters, and .FALSE. if the parameter values are the same as for the
      ! previous call.

      ! Version 1.10, 1 August 1996
      ! Changes from version 1.01
      ! 1. The random_order, random_Poisson & random_binomial routines have been
      !    replaced with more efficient routines.
      ! 2. A routine, seed_random_number, has been added to seed the uniform random
      !    number generator.   This requires input of the required number of seeds
      !    for the particular compiler from a specified I/O unit such as a keyboard.

      !     Author: Alan Miller
      !             CSIRO Division of Mathematics & Statistics
      !             Private Bag 10, Clayton South MDC
      !             Clayton 3169, Victoria, Australia
      !     Phone: (+61) 3 9545-8036      Fax: (+61) 3 9545-8080
      !     e-mail: Alan.Miller @ mel.dms.csiro.au
      !     WWW-page:  http://www.mel.dms.csiro/~alan

      REAL, PRIVATE             :: zero = 0.0, half = 0.5, one = 1.0, two = 2.0
      PRIVATE                   :: integral


      CONTAINS



      REAL FUNCTION random_normal()

      ! Adapted from the following Fortran 77 code
      !      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
      !      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
      !      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.

      !  The function random_normal() returns a normally distributed pseudo-random
      !  number with zero mean and unit variance.

      !  The algorithm uses the ratio of uniforms method of A.J. Kinderman
      !  and J.F. Monahan augmented with quadratic bounding curves.

      IMPLICIT NONE

      !     Local variables
      REAL     :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,           
     &     r1 = 0.27597, r2 = 0.27846, u, v, x, y, q
                  

      !     Generate P = (u,v) uniform in rectangle enclosing acceptance region

      DO
        CALL RANDOM_NUMBER(u)
        CALL RANDOM_NUMBER(v)
        v = 1.7156 * (v - half)

      !     Evaluate the quadratic form
        x = u - s
        y = ABS(v) - t
        q = x**2 + y*(a*y - b*x)

      !     Accept P if inside inner ellipse
        IF (q < r1) EXIT
      !     Reject P if outside outer ellipse
        IF (q > r2) CYCLE
      !     Reject P if outside acceptance region
        IF (v**2 < -4.0*LOG(u)*u**2) EXIT
      END DO

      !     Return ratio of P's coordinates as the normal deviate
      random_normal = v/u
      RETURN

      END FUNCTION random_normal
      
            function ran1(idum)
      parameter (ia=16807,im=2147483647,am=1./im,iq=127773,ir=2836
     & ,ntab=32,ndiv=1+(im-1)/ntab,eps=1.2E-7,rnmx=1.-eps)
      integer iv(ntab)
      save iv,iy
      data iv /ntab*0/, iy /0/

      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=ntab+8,1,-1
          k=idum/iq
          idum=ia*(idum-k*iq)-ir*k
          if (idum.lt.0) idum=idum+im
          if (j.le.ntab) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/iq
      idum=ia*(idum-k*iq)-ir*k
      if (idum.lt.0) idum=idum+im
      j=1+iy/ndiv
      iy=iv(j)
      iv(j)=idum
      ran1=min(am*iy,rnmx)

      return
      end function ran1



      end module MOD_random
