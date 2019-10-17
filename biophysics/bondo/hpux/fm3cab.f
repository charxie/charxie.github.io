C  IDENTIFY PERMUTATION # FOR INPUT BOND ORBITAL                                
C  INPUT BOND IS (A,B,C) = (IM,IL,IS)                                           

      SUBROUTINE FM3CAB(B1,B2,B3,NUM,AB1,AB2,AB3)                               
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                        
      DIMENSION B(3),AB(3),IPERM(3,6)                                           
      DATA IPERM/1,2,3,3,1,2,2,3,1,2,1,3,1,3,2,3,2,1/                           
      B(1)=DABS(B1)                                                             
      B(2)=DABS(B2)                                                             
      B(3)=DABS(B3)                                                             
      IF((B(1).LT.B(2)).OR.(B(2).LT.B(3))) GO TO 10                             
      IP=1                                                                      
      GO TO 100                                                                 
   10 IF((B(3).LT.B(1)).OR.(B(1).LT.B(2))) GO TO 20                             
      IP=2                                                                      
      GO TO 100                                                                 
   20 IF((B(2).LT.B(3)).OR.(B(3).LT.B(1))) GO TO 30                             
      IP=3                                                                      
      GO TO 100                                                                 
   30 IF((B(2).LT.B(1)).OR.(B(1).LT.B(3))) GO TO 40                             
      IP=4                                                                      
      GO TO 100                                                                 
   40 IF((B(1).LT.B(3)).OR.(B(3).LT.B(2))) GO TO 50                             
      IP=5                                                                      
      GO TO 100                                                                 
   50 IF((B(3).LT.B(2)).OR.(B(2).LT.B(1))) GO TO 60                             
      IP=6                                                                      
      GO TO 100                                                                 
   60 CALL ABORT('FM3CAB')                                                      
  100 CONTINUE                                                                  
      B(1)=B1                                                                   
      B(2)=B2                                                                   
      B(3)=B3                                                                   
      IL=IPERM(1,IP)                                                            
      IM=IPERM(2,IP)                                                            
      IS=IPERM(3,IP)                                                            
      IF(NUM.GE.2) GO TO 120                                                    
C  NUM=1, (A*B,-A**2-C**2,B*C)                                                  
      AB(IL)=-B(IM)**2-B(IS)**2                                                 
      AB(IM)=B(IL)*B(IM)                                                        
      AB(IS)=B(IL)*B(IS)                                                        
      GO TO 130                                                                 
C  NUM=2, (C,0,-A)                                                              
  120 CONTINUE                                                                  
      AB(IL)=0.0                                                                
      AB(IM)=B(IS)                                                              
      AB(IS)=-B(IM)                                                             
  130 CONTINUE                                                                  
C  NORMALIZE ANTIBOND                                                           
      T=0.0                                                                     
      DO 140 I=1,3                                                              
      T=T+AB(I)**2                                                              
  140 CONTINUE                                                                  
      T=1.0/SQRT(T)                                                             
      DO 150 I=1,3                                                              
  150 AB(I)=T*AB(I)                                                             
C  RETURN ANTIBOND                                                              
      AB1=AB(1)                                                                 
      AB2=AB(2)                                                                 
      AB3=AB(3)                                                                 
      RETURN                                                                    
      END                                                                       
