      FUNCTION MATCH(I1,I2,I3,J1,J2)                                            
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                        
      DIMENSION I(3),J(2)                                                       
C  MATCH=1 IF (I1,I2,I3) CONTAINS (J1,J2), 0 OTHERWISE                          
      MATCH=0                                                                   
      I(1)=I1                                                                   
      I(2)=I2                                                                   
      I(3)=I3                                                                   
      J(1)=J1                                                                   
      J(2)=J2                                                                   
      DO 20 K=1,2                                                               
      DO 10 L=1,3                                                               
      IF(I(L).EQ.J(K)) GO TO 20                                                 
   10 CONTINUE                                                                  
      RETURN                                                                    
   20 CONTINUE                                                                  
      MATCH=1                                                                   
      RETURN                                                                    
      END                                                                       
