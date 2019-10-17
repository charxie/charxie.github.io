C  CONSTRUCTS ORTHOGONAL MATRIX T FOR TRANSFORMATION FROM AO'S TO               
C  NATURAL HYBRID BOND ORBITALS, USING INPUT DENSITY MATRIX 'DM'.               

      SUBROUTINE NATHYB(DM,T,N,NATOMS)  
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      INTEGER CZ,U,UL 
      character*4 symbol,resname,segname
      integer idres 
      common/iden/idres(nazr),symbol(nazr),resname(nazr),segname(nazr),
     :rx(nazr),ry(nazr),rz(nazr)
      COMMON/LBL/LABEL(NBSZR,6),IBX(NBSZR) 
      COMMON/INFO/NATOMS1,CH,MULT,AN(NAZR),Z(NAZR),Y(NAZR),X(NAZR),
     :N1,NTA,INO1,IPR
      COMMON/INFO1/CZ(NAZR),U(NBSZR),UL(NAZR),LL(NAZR),
     $NELECS,NOCCA,NOCCB
      COMMON/Q/Q(4,NBSZR) 
      COMMON/POL/POL(NBSZR,3),IATHY(NBSZR,3),INO(NAZR)
      DIMENSION DM(NBSZR,NBSZR),T(NBSZR,NBSZR),
     $BLK(12,12),C(12,12),EVAL(12),BORB(12)  
      DIMENSION NAME(3)                                                         
      DATA NAME,ISTAR,IBLNK/'LP','BD','3C','*',' '/                             
      DATA OCCMX/1.9D0/

C  ZERO ARRAYS T, Q, INO, POL, IATHY                                            
      DO 20 I=1,N                                                               
      DO 10 K=1,3                                                               
      POL(I,K)=0.0                                                              
      IATHY(I,K)=0                                                              
   10 Q(K,I)=0.0D0                                                              
      Q(4,I)=0.0D0                                                              
      DO 20 J=1,N                                                               
   20 T(I,J)=0.0D0                                                              
      DO 30 I=1,NAZR
   30 INO(I)=0                                                                  
      IBD=0                                                                     
      NA1=NATOMS+1  

      DO 100 IB1=1,NA1                                                          
      IB=IB1-1                                                                  
      DO 90 IC1=2,NA1                                                           
      IC=IC1-1                                                                  
      IF(IB.NE.0) GO TO 40                                                      

C  LONE PAIRS                                                               

      NCTR=1                                                                    
      IAT1=IC                                                                   
      IAT2=0                                                                    
      IAT3=0 

      GO TO 70                                                          

C  BOND PAIRS                                                                

   40 NCTR=2                                                                    
      IAT1=IB                                                                   
      IAT2=IC                                                                   
      IAT3=0                                                                    
      IF(IAT2.LE.IAT1) GO TO 90

   70 CONTINUE   
      IF(IBD.GE.NOCCA) GO TO 120 

      IF(IAT2 .EQ. 0)  GO TO  71
      IF ((((rx(IAT1)-rx(IAT2))**2.0) +
     :     ((ry(IAT1)-ry(IAT2))**2.0) +
     :     ((rz(IAT1)-rz(IAT2))**2.0) ) .GT. 4.0 ) GO TO 90
   71 CONTINUE
      PRINT *, 'IAT1, IAT2, IAT3', IAT1, IAT2, IAT3

      CALL LOAD(DM,IAT1,IAT2,IAT3,BLK,NB)
      CALL JACVEC(NB,BLK,EVAL,C,12)
      DO 80 IRNK=1,NB                                                           
      CALL EXTRCT(C,EVAL,BORB,OCC,NB,IRNK) 
      odiffa = dabs(occ-occmx)
      if(odiffa.ge.0.15D0) go to 85
      print *, 'occ', occ

      IBD=IBD+1 
      IF(NCTR.EQ.1) THEN
      CALL DEPLET(DM,BORB,OCC,NB,IAT1,IAT2,IAT3)
      ENDIF
      CALL STASH(BORB,OCC,IBD,IAT1,IAT2,IAT3)                                   
      LABEL(IBD,1)=NAME(NCTR)
      LABEL(IBD,2)=IBLNK 
      LABEL(IBD,3)=IRNK             
      LABEL(IBD,4)=IAT1                                                         
      LABEL(IBD,5)=IAT2                                                         
      LABEL(IBD,6)=IAT3 

   80 CONTINUE                                                                  
   85 CONTINUE                                                                  
   90 CONTINUE                                                                  
  100 CONTINUE                                                                  
  120 CONTINUE 

      IBO=IBD                                                                   
      DO 150 I=1,IBO                                                            
      IF(LABEL(I,1).EQ.NAME(1)) GO TO 150                                       
      NAB=1                                                                     
      IF(LABEL(I,1).EQ.NAME(3)) NAB=2                                           
      DO 140 IAB=1,NAB                                                          
      IBD=IBD+1                                                                 
      DO 130 J=1,6                                                              
  130 LABEL(IBD,J)=LABEL(I,J)                                                   
      LABEL(IBD,2)=ISTAR                                                        
  140 CONTINUE                                                                  
  150 CONTINUE                                                                  
      IF(IBD.EQ.N) GO TO 160                                                    
C  MISCOUNTED BOND ORBITALS...EXIT                                              
      WRITE(6,1000) IBD,N                                                       
      CALL ABORT('NATHYB')                                                      
  160 CONTINUE                                                                  

      CALL REFORM(T,NBSZR,N,NATOMS,ipr)

      RETURN                                                                    
 1000 FORMAT(2X,'*** FOUND ',I4,' ORBITALS, BUT DIM. SHOULD BE ',I4,'*')        
      END                                                                       
