      SUBROUTINE HARMTR(T,MAXL,E)  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      DIMENSION T(9,9),E(3)        
      COST = E(3)   
      IF((1.D0-COST**2).GT.0.0000000001) GO TO 20 
   10 SINT = 0.D0   
      GO TO 30      
   20 SINT=DSQRT(1.D0-COST**2)     
   30 continue      
      IF(SINT.GT.0.000001D0)  GO TO 50 
   40 COSP = 1.D0   
      SINP = 0.D0   
      GO TO 70      
   50 COSP = E(1)/SINT 
   60 SINP = E(2)/SINT 
   70 continue      
      DO 80 I=1,9   
      DO 80 J=1,9   
   80 T(I,J) = 0.D0 
      T(1,1) =1.D0  
      IF (MAXL.GT.1) GO TO 100     
   90 IF (MAXL.GT.0) GO TO 110     
      GO TO 120     
  100 COS2T = COST**2-SINT**2      
      SIN2T = 2.D0*SINT*COST       
      COS2P = COSP**2-SINP**2      
      SIN2P = 2.D0*SINP*COSP       
C     TRANSFORMATION MATRIX ELEMENTS FOR D FUNCTIONS 
      SQRT3=DSQRT(3.D0) 
      T(5,5) = (3.D0*COST**2-1.D0)/2.D0 
      T(5,6) = -SQRT3   *SIN2T/2.D0 
      T(5,8) = SQRT3   *SINT**2/2.D0 
      T(6,5) = SQRT3   *SIN2T*COSP/2.D0 
      T(6,6) = COS2T*COSP 
      T(6,7) = -COST*SINP 
      T(6,8) =-T(6,5)/SQRT3        
      T(6,9) = SINT*SINP 
      T(7,5) = SQRT3   *SIN2T*SINP/2.D0 
      T(7,6) = COS2T*SINP 
      T(7,7) = COST*COSP 
      T(7,8) = -T(7,5)/SQRT3       
      T(7,9) = -SINT*COSP 
      T(8,5) = SQRT3   *SINT**2*COS2P/2.D0        
      T(8,6) = SIN2T*COS2P/2.D0    
      T(8,7) = -SINT*SIN2P 
      T(8,8) = (1.D0+COST**2)*COS2P/2.D0 
      T(8,9) = -COST*SIN2P 
      T(9,5) = SQRT3   *SINT**2*SIN2P/2.D0        
      T(9,6) = SIN2T*SIN2P/2.D0    
      T(9,7) = SINT*COS2P 
      T(9,8) = (1.D0+COST**2)*SIN2P/2.D0 
      T(9,9) = COST*COS2P 
  110 continue      
C     TRANSFORMATION MATRIX ELEMENTS FOR P FUNCTIONS 
      T(2,2) = COST*COSP 
      T(2,3) = -SINP 
      T(2,4) = SINT*COSP 
      T(3,2) = COST*SINP 
      T(3,3) = COSP 
      T(3,4) = SINT*SINP 
      T(4,2) = -SINT 
      T(4,4) = COST
  120 continue
      return
      END