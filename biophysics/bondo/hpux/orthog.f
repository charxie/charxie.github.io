C  REPLACES INPUT T(N,N) BY THE ORTHOGONAL MATRIX 'CLOSEST' TO T IN             
C  MEAN-SQUARED SENSE (LOWDIN SYMMETRIC ORTHOGONALIZATION).                     
C                                                                               
C  USES B AND LOWER TRIANGLE OF A, OLD T LEFT IN B ON RETURN.                   

      SUBROUTINE ORTHOG(T,NDIM,N) 
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      DIMENSION T(NDIM,NDIM)                                                    
      COMMON/ARRAYS/A(NBSZR,NBSZR),B(NBSZR,NBSZR),D(NBSZR,NBSZR)
      COMMON/GAB/DUM1(3*NBSZR),EIG(NBSZR),
     $DUM2(NAIGAIO2-(4*NBSZR))

      RHO=1.D-6                                                                 
C  FORM OVERLAP MATRIX S = T(TRANSP)*T, STORE IN L.T. OF A...                   
      DO 20 J=1,N                                                               
      DO 20 I=J,N                                                               
      TEMP=0.0D0                                                                
      DO 10 K=1,N                                                               
   10 TEMP=TEMP+T(K,I)*T(K,J)                                                   
   20 A(I,J)=TEMP                                                               
C  DIAGONALIZE OVERLAP MATRIX...                                                
      CALL EIGN(N,RHO)                                                          
C  FORM INVERSE SQUARE ROOT OF S, STORE IN L.T. OF A...                         
      DO 30 K=1,N                                                               
   30 EIG(K)=1.0D0/DSQRT(EIG(K))                                                
      DO 50 J=1,N                                                               
      DO 50 I=J,N                                                               
      TEMP=0.0D0                                                                
      DO 40 K=1,N                                                               
   40 TEMP=TEMP+EIG(K)*B(I,K)*B(J,K)                                            
   50 A(I,J)=TEMP                                                               
C  FORM TNEW = T*S**(-1/2) , STORE IN B...                                      
      DO 70 J=1,N                                                               
      DO 70 I=1,N                                                               
      TEMP=0.0D0                                                                
      DO 60 K=1,N                                                               
      IF(K.GE.J) AKJ=A(K,J)                                                     
      IF(K.LT.J) AKJ=A(J,K)                                                     
   60 TEMP=TEMP+T(I,K)*AKJ                                                      
   70 B(I,J)=TEMP                                                               
C  REPLACE INPUT T BY ORTHOGONALIZED T, STORE OLD T IN B...                     
      DO 80 J=1,N                                                               
      DO 80 I=1,N                                                               
      TEMP=B(I,J)                                                               
      B(I,J)=T(I,J)                                                             
      T(I,J)=TEMP                                                               
   80 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
