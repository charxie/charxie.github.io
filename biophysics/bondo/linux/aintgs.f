      SUBROUTINE AINTGS(X,K)                                                    
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      COMMON/AUXINT/A(17),B(17)                                                 
      A(1) =DEXP(-X)/X                                                          
      do 10  I=1,K                                                              
   10 A(I+1) =(A(I)*dble(I)+DEXP(-X))/X                                       
      K2=K+2
      do 20 I=K2,17
      A(I)=0.0D0
   20 continue
      return
      END