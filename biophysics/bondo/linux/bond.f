      SUBROUTINE BOND(J) 
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      character*3 BANANA(3,2)
      character*3 nden
      integer CH,AN,CZ,U,UL,OCA,OCB,TOPO                                        
      common/INFO/NA,CH,MULT,AN(NAZR),C(NAZR,3),NO,NTA,N,IPR
      common/INFO1/CZ(NAZR),U(NBSZR),UL(NAZR),LL(NAZR),NE,OCA,OCB 
      common/EXTRA/T(NBSZR,NBSZR),IHYB(NAZR),
     :             TOPO(NAZR,NAZR),T2J(NBSZR)
      common/NHYBRD/nden(nbszr),iden(NBSZR,2)
      dimension ELNG(10)
      character LBD*2,LBDS*3,LLPR*3,LP1*2,LP1S*3,LP2*2,LP2S*3
      character LBL*3

      DATA BANANA/'B1 ','B2 ','B3 ','B1*','B2*','B3*'/                             
      DATA ELNG/2.75D0,0.00D0,0.97D0,1.47D0,2.01D0,
     :          2.50D0,3.07D0,3.50D0,4.10D0,0.00D0/
      DATA LBD,LBDS,LLPR,LP1,LP1S,LP2,LP2S/
     :    'BD','BD*','LPR','P1','P1*','P2','P2*'/
     
      DO 10 I=1,NO                                                              
   10 T2J(I)=0.0D0                                                              
      ITYPE=1                                                                   
      ILAST=0                                                                   
      DO 20 K=1,J                                                               
      IF(iden(K,1).LT.ILAST) ITYPE=-1                                           
      ILAST=iden(K,1)
   20 continue                                                                  
      I=iden(J,1)                                                               
      K=iden(J,2)                                                               
      LBL=nden(J)
      IF(K.GT.NA) ITYPE=0                                                       
      NA1=NA+1                                                                  
      IF(TOPO(I,K).GT.0) GO TO 80                                               
C  CHECK FOR BANANA BONDS...                                                    
      ISGN=2-(ITYPE+1)/2                                                        
      IF(ITYPE.EQ.0) GO TO 80                                                   
      DO 30 NUM=1,3                                                             
      IF(LBL.EQ.BANANA(NUM,ISGN)) GO TO 40                                      
   30 continue                                                                  
      GO TO 80                                                                  
   40 continue                                                                  
      NG=0                                                                      
      DO 50 JG=NA1,NTA                                                          
      IF((TOPO(I,JG).NE.0).AND.(TOPO(K,JG).NE.0)) NG=NG+1                       
      IF(NG.EQ.NUM) GO TO 60                                                    
   50 continue                                                                  
      call abortBondo('BANANA')                                                      
   60 continue                                                                  
      NUM=0                                                                     
      NUM2=0                                                                    
      DO 70 JT=1,JG                                                             
      IF(TOPO(I,JT).GT.0) NUM=NUM+1                                             
      IF(TOPO(K,JT).GT.0) NUM2=NUM2+1                                           
   70 continue                                                                  
      GO TO 150                                                                 
C  SIGMA BONDS AND LONE PAIRS                                                   
   80 continue                                                                  
      IF((LBL.NE.LBD).AND.(LBL.NE.LBDS).AND.(LBL.NE.LLPR))GO TO 110             
      NUM=0
      NUM2=0
      DO 100 L=1,NTA
      IF(L.GT.K) GO TO 90                                                       
      IF(TOPO(I,L).GT.0) NUM=NUM+1                                              
   90 IF(L.GT.I) GO TO 100                                                      
      IF(TOPO(K,L).GT.0) NUM2=NUM2+1                                            
  100 continue                                                                  
      GO TO 150                                                                 
C  FIRST PI BOND                                                                
  110 continue                                                                  
      IF((LBL.NE.LP1).AND.(LBL.NE.LP1S)) GO TO 140                              
      NUM=IHYB(I)+2                                                             
      NUM2=IHYB(K)+2                                                            
      DO 130 L=1,NTA                                                            
      IF(L.GE.K) GO TO 120                                                      
      IF(TOPO(I,L).GT.1) NUM=NUM+TOPO(I,L)-1                                    
  120 IF(L.GE.I) GO TO 130                                                      
      IF(TOPO(K,L).GT.1) NUM2=NUM2+TOPO(K,L)-1                                  
  130 continue                                                                  
      GO TO 150                                                                 
C  SECOND PI BOND                                                               
  140 continue                                                                  
      NUM=4                                                                     
      NUM2=4                                                                    
      IF((LBL.NE.LP2).AND.(LBL.NE.LP2S)) call abortBondo('PI ORB')                   
  150 continue                                                                  
      IF(ITYPE.EQ.0) GO TO 190                                                  
      COEF=1.0D0                                                                
      IF(AN(I).EQ.AN(K))GO TO 170                                               
      R=DSQRT((C(I,1)-C(K,1))**2+(C(I,2)-C(K,2))**2+(C(I,3)-C(K,3))**2)         
      DM=R*(1.0D0-DEXP(-.25*(ELNG(AN(I))-ELNG(AN(K)))**2))                      
      IF(ELNG(AN(I)).GT.ELNG(AN(K)))DM=-DM                                      
  160 COEF=DSQRT((R+DM)/(R-DM))                                                 
  170 ANORM=DSQRT(1.0D0+COEF**2)                                                
      IF(ITYPE)180,190,200                                                      
  180 II=LL(I)+NUM-1                                                            
      IK=LL(K)+NUM2-1                                                           
      T2J(II)=COEF/ANORM                                                        
      T2J(IK)=-1.0D0/ANORM                                                      
      GO TO 210                                                                 
  190 II=LL(I)+NUM-1                                                            
      T2J(II)=1.0D0                                                             
      GO TO 210                                                                 
  200 II=LL(I)+NUM-1                                                            
      IK=LL(K)+NUM2-1                                                           
      T2J(II)=1.0D0/ANORM                                                       
      T2J(IK)=COEF/ANORM                                                        
  210 continue                                                                  

      return                                                                    
      END