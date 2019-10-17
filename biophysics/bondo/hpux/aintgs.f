      SUBROUTINE AINTGS(X,K)                                                    
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      COMMON/AUXINT/A(17),B(17)                                                 
      A(1) =DEXP(-X)/X                                                          
      DO 10  I=1,K                                                              
   10 A(I+1) =(A(I)*DFLOAT(I)+DEXP(-X))/X                                       
      K2=K+2                                                                    
      DO 20 I=K2,17                                                             
      A(I)=0.0D0                                                                
   20 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
