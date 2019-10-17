C  PRINTS OUT TABLE OF INTERNAL/DIHEDRAL BOND ANGLES (DEGREES) BETWEEN          
C  BONDED ROW ATOMS (AR-BR) AND COLUMN ATOMS (AC-BC). IF AR-BR AND AC-BC        
C  BONDS SHARE A COMMON VERTEX ATOM (SAY, AC=BR), TABLE ENTRY IS THE            
C  INTERIOR ANGLE AR-BR-BC.  IF AR-BR AND AC-BC ARE DIRECTLY BONDED             
C  (BY BR-AC, SAY), ENTRY IS THE DIHEDRAL ANGLE AR-BR-AC-BC, WITH LINE          
C  OF SIGHT FROM R-ATOM TO C-ATOM.  ENTRY IS ZERO OTHERWISE.                    

      SUBROUTINE PRANG(LFN)
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                        
      INTEGER AN,TOPO                                                           
      COMMON/EXTRA/TAB(NBSZR,NBSZR),IHYB(NAZR),TOPO(NAZR,NAZR),
     $T2J(NBSZR)
      COMMON/INFO/NATOMS,ICH,MULT,AN(NAZR),C(NAZR,3),NDUM,NTOT,INO,IP

      DIMENSION LIST(40,2),ISYM(10)                                             
      DATA ISYM/'GH',' H','HE','LI','BE',' B',' C',' N',' O',' F'/              
      DATA DASH1,DASH2/'------','-----'/                                        
      IX(L)=1+AN(LIST(L,1))                                                     
      IY(L)=1+AN(LIST(L,2))                                                     
C                                                                               
C  TABULATE 'LIST' OF BONDED ATOMS FROM ARRAY 'TOPO'.  LIST(N,1) AND            
C  LIST(N,2) ARE THE ATOMS OF THE N'TH BOND.                                    
C                                                                               
      IF(NTOT.LE.3) RETURN                                                      
      LMX=0                                                                     
      NM1=NTOT-1                                                                
      DO 10 I=1,NM1                                                             
      I1=I+1                                                                    
      DO 10 J=I1,NTOT                                                           
      IF(TOPO(I,J).EQ.0) GO TO 10                                               
      LMX=LMX+1                                                                 
      LIST(LMX,1)=I                                                             
      LIST(LMX,2)=J                                                             
   10 CONTINUE                                                                  
C                                                                               
C  FIND ATOM PAIRS I1-J1 AND I2-J2 FOR NEXT PAIR OF BONDS...                    
C                                                                               
      LM1=LMX-1                                                                 
      DO 100 L1=1,LM1                                                           
      LP1=L1+1                                                                  
      DO 100 L2=LP1,LMX                                                         
      TAB(L1,L2)=0.0D0                                                          
      I1=LIST(L1,1)                                                             
      J1=LIST(L1,2)                                                             
      I2=LIST(L2,1)                                                             
      J2=LIST(L2,2)                                                             
C                                                                               
C  CALCULATE INTERIOR ANGLE IF I1-J1 SHARES AN ATOM WITH I2-J2...               
C                                                                               
      IF(I1.NE.I2) GO TO 20                                                     
      TAB(L1,L2)=ANGLE(J1,I1,J2)                                                
      GO TO 90                                                                  
   20 IF(I1.NE.J2) GO TO 30                                                     
      TAB(L1,L2)=ANGLE(J1,I1,I2)                                                
      GO TO 90                                                                  
   30 IF(J1.NE.I2) GO TO 40                                                     
      TAB(L1,L2)=ANGLE(I1,J1,J2)                                                
      GO TO 90                                                                  
   40 IF(J1.NE.J2) GO TO 50                                                     
      TAB(L1,L2)=ANGLE(I1,J1,I2)                                                
      GO TO 90                                                                  
C                                                                               
C  CALCULATE DIHEDRAL ANGLE IF I1-J1 LINKED TO I2-J2 BY A BOND...               
C                                                                               
   50 IF(TOPO(I1,I2).EQ.0) GO TO 60                                             
      TAB(L1,L2)=DIHED(J1,I1,I2,J2)                                             
      GO TO 90                                                                  
   60 IF(TOPO(I1,J2).EQ.0) GO TO 70                                             
      TAB(L1,L2)=DIHED(J1,I1,J2,I2)                                             
      GO TO 90                                                                  
   70 IF(TOPO(J1,I2).EQ.0) GO TO 80                                             
      TAB(L1,L2)=DIHED(I1,J1,I2,J2)                                             
      GO TO 90                                                                  
   80 IF(TOPO(J1,J2).EQ.0) GO TO 90                                             
      TAB(L1,L2)=DIHED(I1,J1,J2,I2)                                             
   90 CONTINUE                                                                  
      TAB(L2,L1)=TAB(L1,L2)                                                     
  100 CONTINUE                                                                  
C                                                                               
C  PRINT OUT BOND-ANGLE TABLE0                                                  
C                                                                               
      WRITE(LFN,1000)                                                           
      NL=1                                                                      
      NU=LMX-1                                                                  
  110 IF((NU-NL).GT.10) NU=NL+10                                                
      WRITE(LFN,1100) (ISYM(IX(L)),LIST(L,1),ISYM(IY(L)),LIST(L,2),  L=N        
     +L,NU)                                                                     
      WRITE(LFN,1200) (DASH1,DASH2,L=NL,NU)                                     
      DO 120 M=NL,LMX                                                           
      MX=M-1                                                                    
      IF(M.GT.NU) MX=NU                                                         
      WRITE(LFN,1300) ISYM(IX(M)),LIST(M,1),ISYM(IY(M)),LIST(M,2),  (TAB        
     +(M,L),L=NL,MX)                                                            
  120 CONTINUE                                                                  
      NL=NU+1                                                                   
      NU=LMX-1                                                                  
      IF(NU.GE.NL) GO TO 110                                                    
      RETURN                                                                    
 1000 FORMAT(//,5X,'TABLE OF INTERIOR/DIHEDRAL BOND ANGLES (DEGREES)0')         
 1100 FORMAT(/,11X,11(2X,A2,I2,'-',A2,I2))                                      
 1200 FORMAT(2X,'-----------',11(A6,A5))                                        
 1300 FORMAT(2X,A2,I2,'-',A2,I2,11(F10.3,1X))                                   
      END                                                                       
