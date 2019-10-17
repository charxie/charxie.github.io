      SUBROUTINE STASH(BORB,OCC,IBD,IAT1,IAT2,IAT3) 
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)        
      integer CZ,U,UL   
      dimension BORB(12),IAT(3),HYB(4)  
      common/POL/POL(NBSZR,3),IATHY(NBSZR,3),INO(NAZR)
      common/INFO1/CZ(NAZR),U(NBSZR),UL(NAZR),LL(NAZR),
     $NELECS,NOCCA,NOCC
      common/Q/Q(4,NBSZR)

      IAT(1)=IAT1   
      IAT(2)=IAT2   
      IAT(3)=IAT3   
      KMAX=0    
      IOCC=1    
      IF(OCC.LT.0.1D0) IOCC=-1  
      do 40 I=1,3   
      IA=IAT(I) 
      IF(IA.EQ.0) GO TO 40      
      NU=UL(IA) 
      NL=LL(IA) 
C  EXTRACT HYBRID FROM BOND ORBITAL FOR ATOM IA 
      KMIN=KMAX+1   
      KMAX=KMAX+NU-NL+1 
      MJ=0      
      do 10 K=KMIN,KMAX 
      MJ=MJ+1   
      HYB(MJ)=BORB(K)   
   10 continue  
C  EXTRACT POLARIZATION COEFFICIENT, STORE IN 'POL' 
      PSQ=0.0D0 
      do 20 J=1,MJ  
   20 PSQ=PSQ+HYB(J)**2 
      SGN=1.0D0 
      IF(HYB(1).LT.0.0D0) SGN=-SGN  
      P=SGN*SQRT(PSQ)   
      POL(IBD,I)=P  
C  ENTER HYBRID IN APPROPRIATE COLUMN OF Q MATRIX IF 'BORB' OCCUPIED    
      IF(IOCC.LE.0) GO TO 40    
C  ONE MORE HYBRID FOR ATOM IA
      INO(IA)=INO(IA)+1 
      NCOL=LL(IA)+INO(IA)-1     
C  PLACE NORMALIZED HYBRID IN APPROPRIATE BLOCK OF Q    
      NH=NU-NL+1    
      do 30 NROW=1,NH 
   30 Q(NROW,NCOL)=HYB(NROW)/P
      IATHY(IBD,I)=INO(IA)
Cs
Cs      PRINT *, '...............'
Cs      PRINT *, 'I IAT(I) IA INO(IA) IATHY(IBD,I)'
Cs      PRINT *,  I,  IAT(I),  IA, INO(IA),  IATHY(IBD,I)
Cs      PRINT *, '...............'
Cs
   40 continue

      return    
      END