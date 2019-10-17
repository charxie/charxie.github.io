C  ZEROS THE 12X12 MATRIX 'BLK' AND LOADS IN ATOMIC BLOCKS OF DENSITY           
C  MATRIX 'DM' FOR THE ATOMS LISTED IN 'IAT'                                    

      SUBROUTINE LOAD(DM,IAT1,IAT2,IAT3,BLK,NB)
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      INTEGER CZ,U,UL                                                           
      DIMENSION BLK(12,12),DM(NBSZR,NBSZR),IAT(3)
      COMMON/INFO1/CZ(NAZR),U(NBSZR),UL(NAZR),LL(NAZR),
     $NELECS,NOCCA,NOCCB

      IAT(1)=IAT1                                                               
      IAT(2)=IAT2                                                               
      IAT(3)=IAT3                                                               
C  ZERO 'BLK'                                                                   
      DO 10 I=1,12                                                              
      DO 10 J=1,12                                                              
   10 BLK(I,J)=0.0D0                                                            
      NROW=0                                                                    
      NCOL=0                                                                    
      DO 50 I=1,3                                                               
      IA=IAT(I)                                                                 
      IF(IA.EQ.0) GO TO 50                                                      
      IU=UL(IA)                                                                 
      IL=LL(IA)                                                                 
      DO 40 IROW=IL,IU                                                          
      NROW=NROW+1                                                               
      NCOL=0                                                                    
      DO 30 J=1,3                                                               
      JA=IAT(J)                                                                 
      IF(JA.EQ.0) GO TO 30                                                      
      JU=UL(JA)                                                                 
      JL=LL(JA)                                                                 
      DO 20 ICOL=JL,JU                                                          
      NCOL=NCOL+1                                                               
      BLK(NROW,NCOL)=DM(IROW,ICOL)                                              
   20 CONTINUE                                                                  
   30 CONTINUE                                                                  
   40 CONTINUE                                                                  
   50 CONTINUE                                                                  
      NB=NROW                                                                   
      RETURN                                                                    
      END                                                                       

C  ZEROS THE 24X24 MATRIX 'BLK' AND LOADS IN ATOMIC BLOCKS OF DENSITY           
C  MATRIX 'DM' FOR THE ATOMS LISTED IN 'IAT'                                    

      SUBROUTINE LOAD6(DM,IAT1,IAT2,IAT3,IAT4,IAT5,IAT6,BLK,NB)
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      INTEGER CZ,U,UL                                                           
      DIMENSION BLK(24,24),DM(NBSZR,NBSZR),IAT(6)
      COMMON/INFO1/CZ(NAZR),U(NBSZR),UL(NAZR),LL(NAZR),
     $NELECS,NOCCA,NOCCB

      IAT(1)=IAT1                                                               
      IAT(2)=IAT2                                                               
      IAT(3)=IAT3
      IAT(4)=IAT4
      IAT(5)=IAT5
      IAT(6)=IAT6
C  ZERO 'BLK'                                                                   
      DO 10 I=1,24                                                              
      DO 10 J=1,24                                                              
   10 BLK(I,J)=0.0D0                                                            
      NROW=0                                                                    
      NCOL=0                                                                    
      DO 50 I=1,6
      IA=IAT(I)                                                                 
      IF(IA.EQ.0) GO TO 50                                                      
      IU=UL(IA)                                                                 
      IL=LL(IA)                                                                 
      DO 40 IROW=IL,IU                                                          
      NROW=NROW+1                                                               
      NCOL=0                                                                    
      DO 30 J=1,6                                                               
      JA=IAT(J)                                                                 
      IF(JA.EQ.0) GO TO 30                                                      
      JU=UL(JA)                                                                 
      JL=LL(JA)                                                                 
      DO 20 ICOL=JL,JU                                                          
      NCOL=NCOL+1                                                               
      BLK(NROW,NCOL)=DM(IROW,ICOL)                                              
   20 CONTINUE                                                                  
   30 CONTINUE                                                                  
   40 CONTINUE                                                                  
   50 CONTINUE                                                                  
      NB=NROW                                                                   
      RETURN                                                                    
      END