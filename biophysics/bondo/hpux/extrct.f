C  EXTRACTS EIGENVECTOR (FROM ARRAY C) OF RANK=IRNK OCCUPANCY TO                
C  RETURN AS THE NATURAL BOND-ORBITAL 'BORB'                                    

      SUBROUTINE EXTRCT(C,EVAL,BORB,OCC,NB,IRNK)                                
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      DIMENSION C(12,12),EVAL(12),BORB(12),IFND(12)                             
      IF(IRNK.GT.NB) CALL ABORT('EXTRCT')                                       
      BOCC=2.1D0                                                                
      OCC=0.0D0                                                                 
      DO 30 IR=1,IRNK                                                           
      IFND(IR)=0                                                                
      IF(IR.GT.1) BOCC=OCC                                                      
      OCC=0.0D0                                                                 
      DO 20 I=1,NB                                                              
      IF(EVAL(I).LT.OCC) GO TO 20                                               
      IF(EVAL(I).GT.BOCC) GO TO 20                                              
      DO 10 K=1,IR                                                              
      IF(IFND(K).EQ.I) GO TO 20                                                 
   10 CONTINUE                                                                  
      OCC=EVAL(I)                                                               
      IOCC=I                                                                    
   20 CONTINUE                                                                  
      IFND(IR)=IOCC                                                             
   30 CONTINUE                                                                  
      DO 40 J=1,NB                                                              
      BORB(J)=C(J,IOCC)                                                         
   40 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       

C  EXTRACTS EIGENVECTOR (FROM ARRAY C) OF RANK=IRNK OCCUPANCY TO                
C  RETURN AS THE NATURAL BOND-ORBITAL 'BORB'                                    

      SUBROUTINE EXTRCT6(C,EVAL,BORB,OCC,NB,IRNK)                                
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      DIMENSION C(24,24),EVAL(24),BORB(24),IFND(24)
      IF(IRNK.GT.NB) CALL ABORT('EXTRCT')                                       
      BOCC=2.1D0                                                                
      OCC=0.0D0                                                                 
      DO 30 IR=1,IRNK                                                           
      IFND(IR)=0                                                                
      IF(IR.GT.1) BOCC=OCC                                                      
      OCC=0.0D0                                                                 
      DO 20 I=1,NB                                                              
      IF(EVAL(I).LT.OCC) GO TO 20                                               
      IF(EVAL(I).GT.BOCC) GO TO 20                                              
      DO 10 K=1,IR                                                              
      IF(IFND(K).EQ.I) GO TO 20                                                 
   10 CONTINUE                                                                  
      OCC=EVAL(I)                                                               
      IOCC=I                                                                    
   20 CONTINUE                                                                  
      IFND(IR)=IOCC                                                             
   30 CONTINUE                                                                  
      DO 40 J=1,NB                                                              
      BORB(J)=C(J,IOCC)                                                         
   40 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       


      SUBROUTINE EXTRCT1(C,EVAL,BORB,OCC,NB,IRNK)                                
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      DIMENSION C(6,6),EVAL(6),BORB(6),IFND(6)                             
      IF(IRNK.GT.NB) CALL ABORT('EXTRCT')                                       
      BOCC=2.1D0                                                                
      OCC=0.0D0                                                                 
      DO 30 IR=1,IRNK                                                           
      IFND(IR)=0                                                                
      IF(IR.GT.1) BOCC=OCC                                                      
      OCC=0.0D0                                                                 
      DO 20 I=1,NB                                                              
      IF(EVAL(I).LT.OCC) GO TO 20                                               
      IF(EVAL(I).GT.BOCC) GO TO 20                                              
      DO 10 K=1,IR                                                              
      IF(IFND(K).EQ.I) GO TO 20                                                 
   10 CONTINUE                                                                  
      OCC=EVAL(I)                                                               
      IOCC=I                                                                    
   20 CONTINUE                                                                  
      IFND(IR)=IOCC                                                             
   30 CONTINUE                                                                  
      DO 40 J=1,NB                                                              
      BORB(J)=C(J,IOCC)                                                         
   40 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
