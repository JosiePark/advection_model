c Module that contains subroutines necessary for spatial bicubic interpolation

      module MOD_bicubic
      
      contains
      
c FUNCTION THAT CALCULATES INVERSE
      
      function inv(A) result(Ainv)
          implicit none
      real*8, dimension(:,:), intent(in) :: A
      real*8, dimension(size(A,1),size(A,2)) :: Ainv

      real*8, dimension(size(A,1)) :: work  ! work array for LAPACK
      integer, dimension(size(A,1)) :: ipiv   ! pivot indices
      integer :: n, info

      ! External procedures defined in LAPACK
      external DGETRF
      external DGETRI

      ! Store A in Ainv to prevent it from being overwritten by LAPACK
      Ainv = A
      n = size(A,1)

      ! DGETRF computes an LU factorization of a general M-by-N matrix A
      ! using partial pivoting with row interchanges.
      call DGETRF(n, n, Ainv, n, ipiv, info)

      if (info /= 0) then
         stop 'Matrix is numerically singular!'
      end if

      ! DGETRI computes the inverse of a matrix using the LU factorization
      ! computed by DGETRF.
      call DGETRI(n, Ainv, n, ipiv, work, n, info)

      if (info /= 0) then
         stop 'Matrix inversion failed!'
      end if
        end function inv
        
c MATRIX MULTIPLICATION
      
        function MATMULT(X,Y) result(XY)
        implicit none
        
        real*8, dimension(:,:), intent(in) :: X,Y
        real*8, dimension(size(X,1),size(Y,2)) :: XY
        real*8  sum
        integer :: n,m,p,i,j,k
        
        if (size(X,2) .ne. size(Y,1)) then
        stop  'X and Y mismatch'
        endif
        
        n = size(X,1)
        m = size(X,2)
        p = size(Y,2)
        

        XY = MATMUL(X,Y)
        
        
        
        end function MATMULT
    
C FINDS THE MATRIX OF COEFFICIENTS
        
        subroutine A_matrix(ii,jj,psi,A_mat)
        
c INPUT : ii,jj : grid size
c INPUT : psi(ii,jj) : stream function
c OUTPUT : A_mat(ii,jj,4,4) : 4 X 4 matrix of coefficients stored at each
c                             grid point
        
        implicit none
        
        integer ii,jj,i,j,k,n,m,ic(4),jc(4)

        real*8 A_mat(ii,jj,4,4)
     &   , x_mat(4,4),y_mat(4,4),psi_mat(4,4)
     & ,yc(4),xc(4),x_inv(4,4),y_inv(4,4),psi_grid(4,4)
     & ,psi(ii,jj)
     & ,B(4,4)
        
      
      !print*, psi
      
      do i = 1,ii
      do j = 1,jj
        
        do k = 1,4
            ic(k) = i + k - 2
            jc(k) = j + k - 2
            
                xc(k) = dfloat(ic(k))
                if(ic(k)<=0) then
                ic(i) = ii+ic(k)
                endif
                if(ic(k)>ii) then
                ic(k) = ic(k) - ii
                endif
                
                
                yc(k) = dfloat(jc(k))
                if(jc(k)<=0) then
                jc(k) = jj + jc(k)
                endif
                if(jc(k) > jj) then
                jc(k) = jc(k) - jj
                endif
                
        enddo

        
        do n = 1,4
        do m = 1,4
        
            X_mat(n,m) = xc(n)**(m-1)
            Y_mat(n,m) = yc(m)**(n-1)
        
        enddo
        enddo
        
        
        X_inv = inv(X_mat)
        Y_inv = inv(Y_mat)
        
        do n = 1,4
        do m = 1,4
            psi_grid(m,n) = psi(ic(m),jc(n))
        enddo
        enddo
        

        
        B = matmult(psi_grid,Y_inv)
        A_mat(i,j,:,:) = matmult(X_inv,B)
        
      enddo
      enddo
        
        end subroutine
        
c FINDS VELOCITY
        
        subroutine bicubic(ii,jj,A_mat,x,y,u,v)

C INPUT : ii,jj
c INPUT : A_mat(ii,jj,4,4)
c INPUT : x,y : interpolation point
c OUPUT : u,v : velocities
        
        implicit none
        
        integer ii,jj,i,j,k
        real*8 A_mat(ii,jj,4,4)
     &   ,x,y,u,v,x_col(1,4),y_col(4),diff_y(4)
     & ,diff_x(1,4), XA(1,4), diff_XA(1,4)
        
        do k =1,4
        x_col(1,k) = x**(k-1)
        y_col(k) = y**(k-1)

        enddo
        
      
        diff_y(1) = 0.
        diff_x(1,1) = 0.

      
        do k = 2,4
        diff_y(k) = (k-1)*y**(k-2)
        diff_x(1,k) = (k-1)*x**(k-2)

        enddo
        
        
        XA = matmult(x_col,A_mat(int(x)+1,int(y)+1,:,:))
      diff_XA = matmult(diff_x,A_mat(int(x)+1,int(y)+1,:,:))
      
      u = -dot_product(XA(1,:),diff_y)
      v = dot_product(diff_XA(1,:),y_col)
      
      end subroutine
      
      end module MOD_bicubic
