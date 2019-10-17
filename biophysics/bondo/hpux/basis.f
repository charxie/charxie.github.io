      SUBROUTINE BASIS  
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      INTEGER CH,AN,TOPO,ULK                                                    
      INTEGER CZ,U,ULIM,OCCA,OCCB                                               
      INTEGER BANANA(3,2),IB(2),IPI(2,2)                                        
      COMMON/INFO/NA,CH,MULT,AN(NAZR),C(NAZR,3),NO,NTA,N,IPR
      COMMON/EXTRA/T(NBSZR,NBSZR),IHYB(NAZR),TOPO(NAZR,NAZR),
     $T2J(NBSZR)
      COMMON/INFO1/CZ(NAZR),U(NBSZR),ULIM(NAZR),LLIM(NAZR),
     $NELECS,OCCA,OCCB
      COMMON/NHYBRD/IDEN(NBSZR,3)
      DATA BANANA/'B1','B2','B3','B1*','B2*','B3*'/                             
      DATA IB/'BD','BD*'/                                                       
      DATA IPI/'P1','P2','P1*','P2*'/                                           
      DATA IN/'LPR'/                                                            
      DO 10 I=1,NTA                                                             
      IHYB(I)=-1                                                                
      DO 10 J=1,NTA                                                             
      IF(TOPO(I,J).GT.0) IHYB(I)=IHYB(I)+1                                      
   10 CONTINUE                                                                  
      K1=1                                                                      
      K2=2                                                                      
      K3=3                                                                      
      NA1=NA+1                                                                  
C  BONDS AND LONE PAIRS FOR ISGN=1, ANTIBONDS FOR ISGN=2                        
      DO 80 ISGN=1,2                                                            
      DO 70 I=1,NA                                                              
      DO 50 J=I,NA                                                              
      IF(J.EQ.I) GO TO 50                                                       
      ITO=TOPO(I,J)                                                             
      IF(ITO.GE.1) GO TO 30                                                     
C  CHECK FOR BANANA BONDS...                                                    
      IF(NA1.GT.NTA) GO TO 50                                                   
      NIJ=0                                                                     
      DO 20 JG=NA1,NTA                                                          
      IF((TOPO(JG,I).EQ.0).OR.(TOPO(JG,J).EQ.0)) GO TO 20                       
      NIJ=NIJ+1                                                                 
      IDEN(K1,1)=BANANA(NIJ,ISGN)                                               
      IDEN(K1,2)=I                                                              
      IDEN(K1,3)=J                                                              
      K1=K1+1                                                                   
      K2=K2+1                                                                   
      K3=K3+1                                                                   
   20 CONTINUE                                                                  
C  SIGMA AND PI BONDS                                                           
   30 CONTINUE                                                                  
      IDEN(K1,1)=IB(ISGN)                                                       
      IDEN(K1,2)=I                                                              
      IDEN(K1,3)=J                                                              
      IF(ITO.LT.2) GO TO 40                                                     
      IDEN(K2,1)=IPI(1,ISGN)                                                    
      IDEN(K2,2)=I                                                              
      IDEN(K2,3)=J                                                              
      IF(ITO.LT.3) GO TO 40                                                     
      IDEN(K3,1)=IPI(2,ISGN)                                                    
      IDEN(K3,2)=I                                                              
      IDEN(K3,3)=J                                                              
   40 K1=K1+ITO                                                                 
      K2=K2+ITO                                                                 
      K3=K3+ITO                                                                 
   50 CONTINUE                                                                  
      IF(ISGN.EQ.2) GO TO 70                                                    
C  LONE PAIRS                                                                   
      IF(NA1.GT.NTA) GO TO 70                                                   
      DO 60 JG=NA1,NTA                                                          
      IF(TOPO(I,JG).EQ.0) GO TO 60                                              
      IF(IHYB(JG).NE.0) GO TO 60                                                
      IDEN(K1,1)=IN                                                             
      IDEN(K1,2)=I                                                              
      IDEN(K1,3)=JG                                                             
      K1=K1+1                                                                   
      K2=K2+1                                                                   
      K3=K3+1                                                                   
   60 CONTINUE                                                                  
   70 CONTINUE                                                                  
   80 CONTINUE                                                                  
C  COUNT NO = NO. ORBITALS0                                                     
      NO=0                                                                      
      DO 140 I=1,NA                                                             
      LLIM(I)=NO+1                                                              
      K=1                                                                       
      IF(AN(I).LT.11) GO TO 100                                                 
   90 NO=NO+9                                                                   
      CZ(I)=AN(I)-10                                                            
      GO TO 130                                                                 
  100 IF(AN(I).LT.3) GO TO 120                                                  
  110 NO=NO+4                                                                   
      CZ(I)=AN(I)-2                                                             
      GO TO 130                                                                 
  120 NO=NO+1                                                                   
      CZ(I)=AN(I)                                                               
  130 CONTINUE                                                                  
      ULIM(I)=NO                                                                
  140 CONTINUE                                                                  
      N=NO                                                                      
C                                                                               
C  FILL U ARRAY--U(J) IDENTIFIES THE ATOM TO WHICH ORBITAL J IS                 
C  ATTACHED, E.G., ORBITAL 32 IS ATTACHED TO ATOM 7, ETC.                       
      DO 150 K=1,NA                                                             
      LLK=LLIM(K)                                                               
      ULK=ULIM(K)                                                               
      LIM=ULK+1-LLK                                                             
      DO 150 I=1,LIM                                                            
      J=LLK+I-1                                                                 
  150 U(J)=K                                                                    
C  ABOVE LINES COPIED FROM S.R. INTGRL                                          
      DO 170 J=1,NO                                                             
      CALL BOND(J)                                                              
      DO 170 I=1,NO                                                             
      TEMP=0.0D0                                                                
      DO 160 K=1,NO                                                             
  160 TEMP=TEMP+T1(I,K)*T2J(K)                                                  
  170 T(I,J)=TEMP                                                               
C  ORTHOGONALIZE TRANSFORMATION T0                                              
cexp3      CALL ORTHOG(T,80,NO) 
      CALL ORTHOG(T,NBSZR,NO) 
cexp3                                                     
cexp3      IF(IPR.GE.3) CALL ANLYZE(T,80) 
      IF(IPR.GE.3) CALL ANLYZE(T,NBSZR)
cexp3                                           
      RETURN                                                                    
      END                                                                       
