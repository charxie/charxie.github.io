C  TRANSFORMS I'TH MATRIX OF COMMON/ARRAYS/ TO BOND-ORBITAL BASIS               
C                                                                               
      SUBROUTINE TRANS(N,INP) 
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      INTEGER TOPO                                                              
      COMMON /ARRAYS/A(NBSZR,NBSZR,3) 
      COMMON/EXTRA/T(NBSZR,NBSZR),IHYB(NAZR),TOPO(NAZR,NAZR),
     $T2J(NBSZR)
      logical gaxpy
      dimension t1(n,n)
      
      gaxpy=.true.

      if(gaxpy) then
         do j = 1 , n
            do i = 1 , n
               a(i,j,2)=0.d0
            enddo
         enddo
         do j = 1 , n
            do k = 1 , n
               do i = 1 , n
                  if(k.gt.i) then
                     term=a(k,i,inp)
                  else
                     term=a(i,k,inp)
                  endif
                  a(i,j,2)=a(i,j,2)+term*t(k,j)
               enddo
            enddo
         enddo
         do i = 1 , n
            do j = 1 , n
               a(j,i,inp)=0.d0
               t1(j,i)=t(i,j)
            enddo
         enddo
         do j = 1 , n
            do k = 1 , n
               do i = 1 , j
                  a(i,j,inp)=a(i,j,inp)+t1(i,k)*a(k,j,2)
               enddo
            enddo
         enddo
         do j = 1 , n
            do i = 1 , j-1
               a(j,i,inp)=a(i,j,inp)
            enddo
         enddo
         return
      endif
      
      DO 10 J=1,N                                                               
      DO 10 K=1,N                                                               
      A(J,K,2)=0.0D0                                                            
      DO 10 L=1,N                                                               
      TERM=A(J,L,INP)
      IF(L.GT.J)TERM=A(L,J,INP)
      A(J,K,2)=A(J,K,2)+TERM*T(L,K)                                             
   10 continue                                                                  
      DO 20 J=1,N                                                               
      DO 20 K=1,J                                                               
      A(J,K,INP)=0.0D0
      DO 20 L=1,N
      A(J,K,INP)=A(J,K,INP)+T(L,J)*A(L,K,2)
   20 continue                                                                  

      return                                                                    
      END