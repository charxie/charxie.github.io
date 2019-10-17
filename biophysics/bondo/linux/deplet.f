C  REMOVES NATURAL BOND ORBITAL 'BORB' (OCCUPANCY = OCC) FROM DENSITY           
C  MATRIX 'DM'                                                                  

      SUBROUTINE DEPLET(DM,BORB,OCC,NB,IAT1,IAT2,IAT3) 
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                        
      INTEGER CZ,U,UL                                                           
      DIMENSION DM(NBSZR,NBSZR),BORB(12),IAT(3)
      COMMON/INFO1/CZ(NAZR),U(NBSZR),UL(NAZR),LL(NAZR),
     $NELECS,NOCCA,NOCCB  

      IAT(1)=IAT1                                                               
      IAT(2)=IAT2                                                               
      IAT(3)=IAT3                                                               
      NROW=0                                                                    
      NCOL=0                                                                    
      DO 40 I=1,3                                                               
      IA=IAT(I)                                                                 
      IF(IA.EQ.0) GO TO 40                                                      
      IU=UL(IA)                                                                 
      IL=LL(IA)                                                                 
      DO 30 IROW=IL,IU                                                          
      NROW=NROW+1                                                               
      NCOL=0                                                                    
      DO 20 J=1,3                                                               
      JA=IAT(J)                                                                 
      IF(JA.EQ.0) GO TO 20                                                      
      JU=UL(JA)                                                                 
      JL=LL(JA)                                                                 
      DO 10 ICOL=JL,JU                                                          
      NCOL=NCOL+1                                                               
      DM(IROW,ICOL)=DM(IROW,ICOL)-OCC*BORB(NROW)*BORB(NCOL)                     
   10 continue                                                                  
   20 continue                                                                  
   30 continue                                                                  
   40 continue                                                                  
      NB=NROW                                                                   

      return                                                                    
      END                                                                       
