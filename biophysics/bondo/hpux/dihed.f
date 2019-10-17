C  CALCULATES DIHEDRAL ANGLE I-J-K-L (IN DEGREES), VIEWED ALONG J TO K          
C  AXIS, CLOCKWISE FROM PLANE OF I-J-K TRIANGLE.                                

      FUNCTION DIHED(I,J,K,L) 
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                        
      DIMENSION PHI(2)                                                          
      COMMON/INFO/N(NAZR+3),C(NAZR,3),N1,N2,N3,N4

C                                                                               
C  FIND ROTATION ANGLES TO ROTATE J-K VECTOR TO Z AXIS...                       
C                                                                               
      X=C(K,1)-C(J,1)                                                           
      Y=C(K,2)-C(J,2)                                                           
      Z=C(K,3)-C(J,3)                                                           
      TH1=DATAN2(X,Z)                                                           
      S1=DSIN(TH1)                                                              
      C1=DCOS(TH1)                                                              
      ZP=X*S1+Z*C1                                                              
      TH2=DATAN2(Y,ZP)                                                          
      S2=DSIN(TH2)                                                              
      C2=DCOS(TH2)                                                              
C                                                                               
C  CALCULATE AZIMUTHAL ANGLES PHI FOR VECTORS I-J AND K-L...                    
C                                                                               
      DO 20 IV=1,2                                                              
      IA=I                                                                      
      IB=J                                                                      
      IF(IV.EQ.1) GO TO 10                                                      
      IA=L                                                                      
      IB=K                                                                      
   10 X=C(IA,1)-C(IB,1)                                                         
      Y=C(IA,2)-C(IB,2)                                                         
      Z=C(IA,3)-C(IB,3)                                                         
      XP=X*C1-Z*S1                                                              
      YP=-X*S1*S2+Y*C2-Z*C1*S2                                                  
   20 PHI(IV)=DATAN2(YP,XP)                                                     
      DIHED=(PHI(2)-PHI(1))*180.0D0/3.1415927D0                                 
      IF(DIHED.LT.0.0D0) DIHED=360.0D0+DIHED                                    

      RETURN                                                                    
      END                                                                       
