C  RETURNS COEFFICIENT, P-CHARACTER, AND DIRECTION OF DIRECTED HYBRID           
C  H(1), H(2), H(3), H(4) FROM TRANSFORMATION MATRIX T.                         
C  POW = 0,1,2... FOR S, SP, SP**2,... HYBRID (POW=100 FOR PURE P ORB.)         
C  COEF = COEFFICIENT OF (NORMALIZED) HYBRID                                    
C  THETA, PHI = POLAR AND AZIMUTHAL ANGLES OF DIRECTED HYBRID.                  

      SUBROUTINE HTYPE(H,COEF,POW,THETA,PHI)                                    
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      DIMENSION H(4)                                                            
      DATA EPS/1.D-4/                                                           
      DATA PI/3.14159265357979D0/                                               
      SQ=0.0                                                                    
      DO 10 K=1,4                                                               
   10 SQ=SQ+H(K)*H(K)                                                           
      COEF=DSQRT(SQ)                                                            
      IF(H(1).LT.0.0D0) COEF=-COEF                                              
      DIFF=DABS(COEF-H(1))                                                      
      IF((DABS(COEF).GT.EPS).AND.(DIFF.GT.EPS)) GO TO 20                        
C  ZERO ORBITAL OR PURE S ORBITAL                                               
      POW=0.0                                                                   
      THETA=0.0                                                                 
      PHI=0.0                                                                   
      RETURN                                                                    
C  START WITH NORMALIZED ORBITAL, POSITIVE S COEFFICIENT.                       
   20 DO 30 K=1,4                                                               
   30 H(K)=H(K)/COEF                                                            
      IF(DABS(H(1)).GT.EPS) GO TO 40                                            
C  PURE P ORBITAL...                                                            
      POW=100.0                                                                 
      GO TO 60                                                                  
   40 POW=1.0D0/H(1)/H(1)-1.0D0                                                 
      IF(POW.LE.0.0D0) POW=0.0D0                                                
      IF(POW.GE.99.0D0) POW=99.0D0                                              
      RAD=DSQRT(POW)                                                            
C  FIND COEFFICIENTS OF NORMALIZED P ORBITAL...                                 
      DO 50 K=2,4                                                               
   50 H(K)=H(K)/RAD/H(1)                                                        
C  FIND POLAR ANGLE THETA AND AZIMUTHAL PHI.                                    
   60 CONTINUE                                                                  
      THETA=-1.0                                                                
      SIN2TH=1.0D0-H(4)*H(4)                                                    
      IF(DABS(H(4)).GT.EPS) GO TO 70                                            
      THETA=PI/2.0D0                                                            
   70 IF(SIN2TH.LE.0.0D0) SIN2TH=0.0D0                                          
      SINTH=DSQRT(SIN2TH)                                                       
      IF(DABS(SINTH).GT.EPS) GO TO 80                                           
      THETA=0.0D0                                                               
      IF(H(4).LT.0.0D0) THETA=PI                                                
      PHI=0.0D0                                                                 
      GO TO 110                                                                 
   80 COSPH=H(2)/SINTH                                                          
      SINPH=H(3)/SINTH                                                          
      IF(DABS(COSPH).GT.EPS) GO TO 90                                           
      PHI=PI/2.0D0                                                              
      IF(SINPH.LT.0.0D0) PHI=-PHI                                               
      GO TO 100                                                                 
   90 PHI=DATAN2(SINPH,COSPH)                                                   
  100 IF(THETA.LT.0.0D0) THETA=DATAN2(SINTH,H(4))                               
  110 PHI=PHI*360.D0/(2.0D0*PI)                                                 
      THETA=THETA*360.D0/(2.0D0*PI)                                             
      RETURN                                                                    
      END                                                                       
