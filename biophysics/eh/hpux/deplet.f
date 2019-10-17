      SUBROUTINE DEPLET(NTIME,BORB,OCC,NB,IAT1,IAT2,IAT3) 
c***********************************************************************
c deplete the density matrix of the contribution from a lone pair
C  REMOVES NATURAL BOND ORBITAL 'BORB' (OCCUPANCY = OCC) FROM DENSITY           
C  MATRIX 'DM'                                                                  
c***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)                                        
      PARAMETER(NATMAX=200,NDMAX=500,NSMP=10)
      INTEGER UL                                                           
      DIMENSION BORB(12),IAT(3)
      COMMON/INFO1/UL(NATMAX),LL(NATMAX),NELECS,NOCCA,NOCCB  
      common /eigt/eigne(ndmax,nsmp),eignvr(ndmax,ndmax,nsmp),
     * 	           eignvi(ndmax,ndmax,nsmp),dm(ndmax,ndmax,nsmp),
     *             pop(ndmax,nsmp)

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
      DM(IROW,ICOL,NTIME)=DM(IROW,ICOL,NTIME)-OCC*BORB(NROW)*BORB(NCOL)                     
   10 CONTINUE                                                                  
   20 CONTINUE                                                                  
   30 CONTINUE                                                                  
   40 CONTINUE                                                                  
      NB=NROW                                                                   

      RETURN                                                                    
      END                                                                       



      SUBROUTINE STASH(BORB,OCC,IBD,IAT1,IAT2,IAT3) 
c*************************************************************************
c seperate bond orbital BORB into polarization coefficients (for storage
c in POL) and hybrids (for storage in Q), keeping track of the number 
c INO(NA) of hybrids accumulated for each atom NA
c*************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NATMAX=200,NDMAX=500)                                        
      INTEGER UL                                                           
      DIMENSION BORB(12),IAT(3),HYB(4)                                          
      COMMON/POL/POL(NDMAX,3),IATHY(NDMAX,3),INO(NATMAX)
      COMMON/INFO1/UL(NATMAX),LL(NATMAX),NELECS,NOCCA,NOCC
      COMMON/Q/Q(4,NDMAX)

      IAT(1)=IAT1                                                               
      IAT(2)=IAT2                                                               
      IAT(3)=IAT3                                                               
      KMAX=0                                                                    
      IOCC=1                                                                    
      IF(OCC.LT.0.1D0) IOCC=-1                                                  
      DO 40 I=1,3                                                               
      IA=IAT(I)                                                                 
      IF(IA.EQ.0) GO TO 40                                                      
      NU=UL(IA)                                                                 
      NL=LL(IA)                                                                 
C  EXTRACT HYBRID FROM BOND ORBITAL FOR ATOM IA                                 
      KMIN=KMAX+1                                                               
      KMAX=KMAX+NU-NL+1                                                         
      MJ=0                                                                      
      DO 10 K=KMIN,KMAX                                                         
      MJ=MJ+1                                                                   
      HYB(MJ)=BORB(K)                                                           
   10 CONTINUE                                                                  
C  EXTRACT POLARIZATION COEFFICIENT, STORE IN 'POL'                             
      PSQ=0.0D0                                                                 
      DO 20 J=1,MJ                                                              
   20 PSQ=PSQ+HYB(J)**2                                                         
      SGN=1.0D0                                                                 
      IF(HYB(1).LT.0.0D0) SGN=-SGN                                              
      P=SGN*SQRT(PSQ)                                                           
      POL(IBD,I)=P                                                              
C  ENTER HYBRID IN APPROPRIATE COLUMN OF Q MATRIX IF 'BORB' OCCUPIED            
      IF(IOCC.LE.0) GO TO 40                                                    
C  ONE MORE HYBRID FOR ATOM IA0                                                 
      INO(IA)=INO(IA)+1                                                         
      NCOL=LL(IA)+INO(IA)-1                                                     
C  PLACE NORMALIZED HYBRID IN APPROPRIATE BLOCK OF Q                            
      NH=NU-NL+1                                                                
      DO 30 NROW=1,NH 
   30 Q(NROW,NCOL)=HYB(NROW)/P

      IATHY(IBD,I)=INO(IA)

   40 CONTINUE

      RETURN                                                                    
      END                                                                       

      SUBROUTINE REFORM(T,NATOM,NDIM)
c*************************************************************************
c Builds the final orthogonal transformation matrix T from polarization
c coefficients (POL) and atomic hybrids(stored in Q, identified in IATHY),
c inserting antibonds as necessary
C  SET UP FINAL FORM OF ORTHOGONAL MATRIX T, USING Q AS TEMPORARY STORAG        
c*************************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NATMAX=200,NDMAX=500)                                        
      INTEGER UL                                                           
      DIMENSION T(NDMAX,NDMAX)                                                    
      DIMENSION S(4,4),EVAL(4),C(4,4),TA(4,4)                                   
      COMMON/Q/Q(4,NDMAX) 
      COMMON/POL/POL(NDMAX,3),IATHY(NDMAX,3),INO(NATMAX)
      COMMON/INFO1/UL(NATMAX),LL(NATMAX),NELECS,NOCCA,NCCB
      COMMON/LBL/LABEL(NDMAX,6),IBX(NDMAX) 
      DATA LLP,LBD,L3C,LSTAR/'LP','BD','3C','*'/
                                      
C  REORDER OCCUPIED BO'S TO PUT LONE PAIRS LAST                                 
      DO 10 NLP=1,NOCCA                                                         
      IF(LABEL(NLP,1).NE.LLP) GO TO 20                                          
   10 CONTINUE                                                                  
   20 NLP=NLP-1                                                                 
      NBD=NOCCA-NLP                                                             
      DO 50 IBD=1,N                                                             
      IF(IBD.GT.NLP) GO TO 30                                                   
C  LONE PAIRS                                                                   
      IBX(IBD)=IBD+NBD                                                          
      GO TO 50                                                                  
C  PAIR BONDS                                                                   
   30 IF(IBD.GT.NOCCA) GO TO 40                                                 
      IBX(IBD)=IBD-NLP                                                          
      GO TO 50                                                                  
C  ANTIBONDS                                                                    
   40 IBX(IBD)=IBD                                                              
   50 CONTINUE                                                                  
C  ZERO ARRAY T...                                                              
      DO 60 I=1,NDIM
      DO 60 J=1,NDIM
   60 T(I,J)=0.0D0                                                              
   70 CONTINUE                                                                  
C  SYMMETRIC ORTHOGONALIZATION OF HYBRIDS0 

      open(file='taohom', unit=10) 
                                    
      DO 170 IA=1,NATOM
      IL=LL(IA)                                                                 
      IU=UL(IA)                                                                 
      NH=IU-IL+1                                                                
Cs      IF(NH.EQ.1) GO TO 170  
Cs Above line is original
Cs
Cs
Cs      PRINT *, 'IA', IA
Cs      PRINT *, 'IL', IL
Cs      PRINT *, 'IU', IU
Cs      PRINT *, 'NH', NH
      IF(NH .EQ.1) THEN
Cs      PRINT *, '...............................'
      Do J=1,NH
      Do I=1,NH
Cs      PRINT *, IL-1+I, IL-1+J, 1.00
      write(10,*)  IL-1+I, IL-1+J, 1.00
      End Do
      End Do
Cs      PRINT *, '...............................'
      GO TO 170
      ELSE
      CONTINUE
      END IF
Cs      
Cs                                                   
C  LOAD IA-BLOCK OF Q INTO TA...                                                
      DO 80 J=1,NH                                                              
      DO 80 I=1,4                                                               
   80 TA(I,J)=Q(I,IL+J-1)                                                       
C  FORM OVERLAP MATRIX S = TA(TRANSP)*TA...                                     
      DO 100 J=1,NH                                                             
      DO 100 I=J,NH                                                             
      TEMP=0.0D0                                                                
      DO 90 K=1,NH                                                              
   90 TEMP=TEMP+TA(K,I)*TA(K,J)                                                 
      S(I,J)=TEMP                                                               
  100 S(J,I)=TEMP                                                               
C  DIAGONALIZE OVERLAP MATRIX...                                                
      CALL JACVEC(NH,S,EVAL,C,4)                                                
C  FORM INVERSE SQUARE ROOT OF S, STORE IN S...                                 
      DO 110 I=1,NH                                                             
  110 EVAL(I)=1.0D0/DSQRT(EVAL(I))                                              
      DO 130 J=1,NH                                                             
      DO 130 I=J,NH                                                             
      TEMP=0.0D0                                                                
      DO 120 K=1,NH                                                             
  120 TEMP=TEMP+EVAL(K)*C(I,K)*C(J,K)                                           
      S(I,J)=TEMP                                                               
  130 S(J,I)=TEMP                                                               
C  FORM NEW TAP=TA*S**(-1/2), STORE IN C...                                     
      DO 150 J=1,NH                                                             
      DO 150 I=1,NH                                                             
      TEMP=0.0D0                                                                
      DO 140 K=1,NH                                                             
  140 TEMP=TEMP+TA(I,K)*S(K,J)                                                  
  150 C(I,J)=TEMP    
Cs
Cs      PRINT *, 'IA', IA
Cs      PRINT *, 'IL', IL
Cs      PRINT *, 'IU', IU
Cs      PRINT *, 'NH', NH
Cs      PRINT *, '...............................'
      Do J=1,NH
      Do I=1,NH
Cs      PRINT *, IL-1+I, IL-1+J, C(I,J)
      write(10,*) IL-1+I, IL-1+J, C(I,J)
      End Do
      End Do
Cs      PRINT *, '..............................'
Cs                                                           
C  REPLACE ORTHOGONALIZED TA IN ARRAY Q...                                      
      DO 160 J=1,NH                                                             
      DO 160 I=1,4                                                              
      Q(I,IL+J-1)=C(I,J)                                                        
  160 CONTINUE                                                                  
  170 CONTINUE
Cs 
      close(unit=10) 
Cs                                                                
C  SYMMETRIC ORTHOGONALIZATION COMPLETE.                                        
C  DEPOSIT FINAL BOND ORBITALS IN MATRIX T                                      
      NBO=0                                                                     
      DO 280 IBD=1,NDIM
      KBD=IBD                                                                   
      IF(LABEL(IBD,2).NE.LSTAR) GO TO 250                                       
C  ANTIBOND ORBITALS0 SEARCH OCCUPIED ORB. LIST TO GET PROPER HYBRIDS...        
C  SEARCH OCCUPIED BOND ORBS. FOR MATCH WITH ANTIBOND ATOMS                     
      DO 240 K=1,NBO                                                            
      DO 180 I=4,6                                                              
      IF(LABEL(K,I).NE.LABEL(IBD,I)) GO TO 240                                  
      IF((LABEL(K,3).LE.0).AND.(LABEL(K,1).EQ.LBD)) GO TO 240                   
C  NEGATIVE IRNK = LABEL(K,3) MEANS BOND ORBITAL WAS ALREADY USED               
  180 CONTINUE                                                                  
C  FOUND MATCH; SET LABEL(K,3)<0 AND RESET POLARIZATION PARAMETERS              
C  FOR ANTIBOND                                                                 
      KBD=K                                                                     
      LABEL(KBD,3)=-LABEL(KBD,3)                                                
      IF(LABEL(IBD,1).NE.L3C) GO TO 230                                         
C  3-CENTER ANTIBONDS FOR (A,B,C)                                               
C      NUM=1                                                                     
C      IF(LABEL(KBD,3).GT.0) NUM=2                                               
C      CALL FM3CAB(POL(KBD,1),POL(KBD,2),POL(KBD,3),NUM,POL(IBD,1),              
C     + POL(IBD,2),POL(IBD,3))                                                   
C      GO TO 250                                                                 
C  2-CENTER ANTIBOND                                                            
  230 POL(IBD,2)=-POL(KBD,1)                                                    
      POL(IBD,1)=POL(KBD,2)                                                     
      GO TO 250                                                                 
  240 CONTINUE                                                                  
C  COULDN'T FIND SUCCESSFUL MATCH...EXIT                                        
      STOP 'ERROR:REFORM'
C  DEPOSIT BOND ORBITALS IN T MATRIX                                            
  250 CONTINUE                                                                  
      DO 270 I=1,3                                                              
      IA=LABEL(IBD,I+3)                                                         
      IF(IA.EQ.0) GO TO 270                                                     
      JL=LL(IA)                                                                 
      JU=UL(IA)                                                                 
      IROW=0                                                                    
      ICOL=JL+IATHY(KBD,I)-1                                                    
      DO 260 J=JL,JU                                                            
      IROW=IROW+1                                                               
  260 T(J,IBX(IBD))=POL(IBD,I)*Q(IROW,ICOL)                                     
  270 CONTINUE                                                                  
      IF(IBD.EQ.KBD) NBO=IBD                                                    
  280 CONTINUE                                                                  
C  RESTORE LABEL(I,3) > 0                                                       
      DO 290 I=1,NDIM
      IF(LABEL(I,3).LT.0) LABEL(I,3)=-LABEL(I,3)                                
  290 CONTINUE                                                                  

      RETURN                                                                    
      END                                                                       
