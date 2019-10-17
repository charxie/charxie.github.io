C  CALCULATES INTERIOR ANGLE I-J-K, USING ATOMIC COORDS. FROM COMMON/INF        
C                                                                               
      FUNCTION ANGLE(I,J,K) 
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                        
      COMMON/INFO/IJK(NAZR+3),C(NAZR,3),N1,N2,N3,N4

      ANGLE=0.0D0                                                               
      R1=0.0D0                                                                  
      R2=0.0D0                                                                  
      R12=0.0D0                                                                 
      do 10 N=1,3                                                               
      C1=C(I,N)-C(J,N)                                                          
      C2=C(K,N)-C(J,N)                                                          
      R1=R1+C1*C1                                                               
      R2=R2+C2*C2                                                               
      R12=R12+C1*C2                                                             
   10 continue                                                                  
      CTH=R12/DSQRT(R1*R2)                                                      
      IF(DABS(CTH).GE.1.0D0) CTH=CTH/DABS(CTH)                                  
      ANGLE=DACOS(CTH)*180.0D0/3.1415927D0                                      

      return
      END