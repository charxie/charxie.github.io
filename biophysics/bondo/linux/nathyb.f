C  CONSTRUCTS ORTHOGONAL MATRIX T FOR TRANSFORMATION FROM AO'S TO 
C  NATURAL HYBRID BOND ORBITALS, USING INPUT DENSITY MATRIX 'DM'. 

      SUBROUTINE NATHYB(DM,T,N,NATOMS)  
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)    
      integer CZ,U,UL 
      character*4 symbol,resname,segname
      character*2 bname
      integer idres 
      common/iden/idres(nazr),symbol(nazr),resname(nazr),segname(nazr),
     :rx(nazr),ry(nazr),rz(nazr)
      common/LBL/bname(nbszr,2),label(NBSZR,3),IBX(NBSZR) 
      common/INFO/NATOMS1,CH,MULT,AN(NAZR),Z(NAZR),Y(NAZR),X(NAZR),
     :N1,NTA,INO1,IPR
      common/INFO1/CZ(NAZR),U(NBSZR),UL(NAZR),LL(NAZR),
     :NELECS,NOCCA,NOCCB
      common/Q/Q(4,NBSZR) 
      common/POL/POL(NBSZR,3),IATHY(NBSZR,3),INO(NAZR)
      dimension DM(NBSZR,NBSZR),T(NBSZR,NBSZR)
      dimension BLK(12,12),C(12,12),EVAL(12),BORB(12)  
      character*2 NAME(3)
      character*1 star,blnk
      DATA NAME,STAR,BLNK/'LP','BD','3C','*',' '/
      DATA OCCMX/1.9D0/

C  ZERO ARRAYS T, Q, INO, POL, IATHY         
      do 20 I=1,N              
      do 10 K=1,3              
      POL(I,K)=0.0             
      IATHY(I,K)=0             
   10 Q(K,I)=0.0D0             
      Q(4,I)=0.0D0             
      do 20 J=1,N              
   20 T(I,J)=0.0D0             
      do 30 I=1,NAZR
   30 INO(I)=0          
      IBD=0             
      NA1=NATOMS+1  

      do 100 IB1=1,NA1         
      IB=IB1-1          
      do 90 IC1=2,NA1          
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

   70 continue   
      IF(IBD.GE.NOCCA) GO TO 120 

      IF(IAT2 .EQ. 0)  GO TO  71
      IF ((((rx(IAT1)-rx(IAT2))**2.0) +
     :     ((ry(IAT1)-ry(IAT2))**2.0) +
     :     ((rz(IAT1)-rz(IAT2))**2.0) ) .GT. 4.0 ) GO TO 90
   71 continue
      PRINT *, 'IAT1, IAT2, IAT3', IAT1, IAT2, IAT3

      call LOAD(DM,IAT1,IAT2,IAT3,BLK,NB)
      call JACVEC(NB,BLK,EVAL,C,12)
      do 80 IRNK=1,NB          
      call EXTRCT(C,EVAL,BORB,OCC,NB,IRNK) 
      odiffa = dabs(occ-occmx)
      if(odiffa.ge.0.15D0) go to 85
      print *, 'occ', occ

      IBD=IBD+1 
      IF(NCTR.EQ.1) THEN
      call DEPLET(DM,BORB,OCC,NB,IAT1,IAT2,IAT3)
      ENDIF
      call STASH(BORB,OCC,IBD,IAT1,IAT2,IAT3)       
      bname(IBD,1)=NAME(NCTR)
      bname(IBD,2)=BLNK
      label(IBD,1)=IRNK
      label(IBD,2)=IAT1
      label(IBD,3)=IAT2

   80 continue
   85 continue
   90 continue
  100 continue
  120 continue

      IBO=IBD           
      do 150 I=1,IBO           
      IF(bname(I,1).EQ.NAME(1)) GO TO 150           
      NAB=1             
      IF(bname(I,1).EQ.NAME(3)) NAB=2        
      do 140 IAB=1,NAB         
      IBD=IBD+1
      do 130 J=1,3
  130 label(IBD,J)=label(I,J)         
      bname(IBD,1)=bname(i,1)
      bname(IBD,2)=STAR
  140 continue          
  150 continue          
      IF(IBD.EQ.N) GO TO 160          
C  MISCOUNTED BOND ORBITALS...EXIT           
      WRITE(6,1000) IBD,N             
      call abortBondo('NATHYB')
  160 continue          

      call REFORM(T,NBSZR,N,NATOMS,ipr)

      return            
 1000 FORMAT(2X,'*** FOUND ',I4,' ORBITALS, BUT DIM. SHOULD BE ',I4,'*') 
      END