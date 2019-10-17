C  GENERATES THE I,J ELEMENT OF THE MATRIX WHICH TRANSFORMS                     
C  ATOMIC ORBITALS TO HYBRID ORBITALS                                           

      FUNCTION T1(I,J)  
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      INTEGER CH,AN,CZ,U,UL,OCA,OCB,TOPO                                        
      COMMON /INFO/NA,CH,MULT,AN(NAZR),C(NAZR,3),NO,NTA,N,IPR
      COMMON/INFO1/CZ(NAZR),U(NBSZR),UL(NAZR),LL(NAZR),
     $NE,OCA,OCB
      COMMON/EXTRA/T(NBSZR,NBSZR),IHYB(NAZR),TOPO(NAZR,NAZR),
     $T2J(NBSZR)

C  IF ORBITALS I AND J BELONG TO DIFFERENT ATOMS, T1=0                          
      T1=0.0D0                                                                  
      IF(U(I).NE.U(J))GO TO 300                                                 
      NAT=U(I)                                                                  
C  IF THE HYBRIDIZATION OF THE ATOM TO WHICH I AND J BELONG IS 0,               
C  T1=DELTA(I,J)                                                                
      IF(IHYB(NAT).EQ.0)GO TO 290                                               
C  IP AND JP IDENTIFY THE ELEMENT OF THE 4X4 BLOCK BELONGING TO ATOM            
C  ('NAT') UNDER CONSIDERATION, WITH JP DENOTING THE NUMBER OF THE              
C  HYBRID TO BE CONSTRUCTED, AND IP THE NUMBER OF THE A.O. WHOSE                
C  CONTRIBUTION IS WANTED.                                                      
      IP=I-LL(NAT)+1                                                            
      JP=J-LL(NAT)+1                                                            
      IDIF=IHYB(NAT)+1-JP                                                       
C  IDIF DETERMINES WHETHER THE HYBRID ORBITAL IS SIGMA OR PI-TYPE (OR           
C  WHICH PI-TYPE ORBITAL IT IS)                                                 
      IPM=IP-1                                                                  
      IF(IDIF.LT.0)GO TO 100                                                    
C  SIGMA-TYPE ORBITALS (IDIF.GE.0)                                              
C  SCAN TOPO TO FIND ATOM TO WHICH HYBRID JP IS DIRECTED, THEN FIND             
C  COEFFICIENT T1 ACCORDING TO THE TYPE OF HYBRIDIZATION...                     
      NHYB=0                                                                    
      DO 10 K=1,NTA                                                             
      IF(TOPO(NAT,K).NE.0)NHYB=NHYB+1                                           
      IF(NHYB.EQ.JP)GO TO 20                                                    
   10 continue                                                                  
   20 NCAT=K                                                                    
      R=0.0D0                                                                   
      DO 30 K=1,3                                                               
      R=R+(C(NCAT,K)-C(NAT,K))**2                                               
   30 continue                                                                  
      R=DSQRT(R)                                                                
      IF(IHYB(NAT)-2)40,60,80                                                   
   40 IF(IP.EQ.1)GO TO 50                                                       
      T1=(C(NCAT,IPM)-C(NAT,IPM))*0.70710678120D0/R                             
      GO TO 300                                                                 
   50 T1=0.70710678120D0                                                        
      GO TO 300                                                                 
   60 IF(IP.EQ.1)GO TO 70                                                       
      T1=(C(NCAT,IPM)-C(NAT,IPM))*0.81649658090D0/R                             
      GO TO 300                                                                 
   70 T1=0.57735026920D0                                                        
      GO TO 300                                                                 
   80 IF(IP.EQ.1)GO TO 90                                                       
      T1=(C(NCAT,IPM)-C(NAT,IPM))*0.86602540380D0/R                             
      GO TO 300                                                                 
C  PI-TYPE ORBITALS0                                                            
C     THE CONTRIBUTION OF P-TYPE ORBITALS IS FOUND BY FORMING THE               
C     CROSS-PRODUCT (USING PVEC) OF THE VECTOR ALONG THE BOND AND               
C     A VECTOR IN THE SIGMA PLANE (IF SP**2 TYPE).  THE ELEMENTS OF             
C     THE CROSS-PRODUCT UNIT VECTOR ARE THE DIRECTION COSINES OF THE            
C     RESULTANT VECTOR WITH THE MOLECULAR COORDINATE SYSTEM, I.E.,              
C     THE COEFFICIENTS OF THE ORIGINAL P-ORBITALS (IN THE MOLECULAR             
C     COORDINATE SYSTEM) TO THE P-ORBITAL FOR THE PARTICULAR PI BOND            
C     OF INTEREST.  FOR A TRIPLE BOND, FIRST PI BOND IS PUT PERPEN-             
C     DICULAR TO THE SURROUNDING MOLECULAR SKELETON (IF POSSIBLE),              
C     AND THE SECOND PERPENDICULAR TO THE FIRST.                                
C                                                                               
C  PI-TYPE ORBITALS0  SCAN TOPO TO FIND PI-BONDED ATOM...                       
   90 T1=0.50D0                                                                 
      GO TO 300                                                                 
  100 IF(IP.EQ.1)GO TO 300                                                      
      IDIF=IABS(IDIF)                                                           
      NUM=0                                                                     
      DO 110 K=1,NTA                                                            
      IF(TOPO(NAT,K).GT.1)NUM=NUM+TOPO(NAT,K)-1                                 
      IF(NUM.GE.IDIF)GO TO 120                                                  
  110 continue                                                                  
  120 NCAT=K                                                                    
      IF(IHYB(NAT).EQ.2)GO TO 250                                               
C  IHYB(NAT)=1 (SP TYPE HYBRID); TRIPLE BOND OR CENTER BOND ON ALLENE           
      IF(TOPO(NAT,NCAT).GT.2)GO TO 150                                          
C  ALLENE TYPE (TWO PI ORBITALS)...                                             
      DO 130 K=1,NTA                                                            
      IF(K.EQ.NAT)GO TO 130                                                     
      IF(TOPO(NCAT,K).GT.0)GO TO 140                                            
  130 continue                                                                  
  140 NCATAT=K                                                                  
      T1=PVEC(NAT,NCATAT,NCAT,IPM)                                              
      GO TO 300                                                                 
C  TRIPLE BOND (TWO PI ORBITALS)...                                             
  150 IF(NAT.LT.NCAT)GO TO 160                                                  
      NTEMP=NAT                                                                 
      NAT=NCAT                                                                  
      NCAT=NTEMP                                                                
C  ...SEARCH LEFT FOR NON-LINEAR LINK...                                        
  160 DO 170 K=1,NTA                                                            
      IF(K.EQ.NCAT)GO TO 170                                                    
      IF(TOPO(NAT,K).GT.0)GO TO 180                                             
  170 continue                                                                  
  180 NCAT2=K                                                                   
      IF(IHYB(NCAT2).GE.2)GO TO 210                                             
C  ...SEARCH RIGHT FOR NON-LINEAR LINK...                                       
      DO 190 K=1,NTA                                                            
      IF(K.EQ.NAT)GO TO 190                                                     
      IF(TOPO(NCAT,K).GT.0)GO TO 200                                            
  190 continue                                                                  
  200 NCAT2=K                                                                   
C  ...IF NO NON-LINEAR LINKAGE, USE NCAT2 (=NCAT3) ANYWAY                       
      IF(IHYB(NCAT2).GE.2) GO TO 205                                            
      NCAT3=NCAT2                                                               
      GO TO 235                                                                 
  205 continue                                                                  
  210 DO 220 K=1,NTA                                                            
      IF((K.EQ.NAT).OR.(K.EQ.NCAT))GO TO 220                                    
      IF(TOPO(NCAT2,K).GT.0)GO TO 230                                           
  220 continue                                                                  
  230 NCAT3=K                                                                   
  235 continue                                                                  
      IF(IDIF.EQ.2)GO TO 240                                                    
C  FIRST PI BOND OF TRIPLE BOND...                                              
      T1=PVEC(NAT,NCAT,NCAT3,IPM)                                               
      GO TO 300                                                                 
C  SECOND PI BOND OF TRIPLE BOND...                                             
  240 VX=PVEC(NAT,NCAT,NCAT3,1)                                                 
      VY=PVEC(NAT,NCAT,NCAT3,2)                                                 
      VZ=PVEC(NAT,NCAT,NCAT3,3)                                                 
      UX=C(NCAT2,1)-C(NAT,1)                                                    
      UY=C(NCAT2,2)-C(NAT,2)                                                    
      UZ=C(NCAT2,3)-C(NAT,3)                                                    
      R=DSQRT(UX**2+UY**2+UZ**2)                                                
      UX=UX/R                                                                   
      UY=UY/R                                                                   
      UZ=UZ/R                                                                   
      IF(IPM.EQ.1)T1=UY*VZ-UZ*VY                                                
      IF(IPM.EQ.2)T1=UZ*VX-UX*VZ                                                
      IF(IPM.EQ.3)T1=UX*VY-UY*VX                                                
      GO TO 300                                                                 
C  IHYB(NAT)=2, SP**2 TYPE (SINGLE PI BOND)0                                    
  250 IF((NAT.LT.NCAT).OR.(IHYB(NCAT).EQ.1))GO TO 260                           
      NTEMP=NAT                                                                 
      NAT=NCAT                                                                  
      NCAT=NTEMP                                                                
  260 DO 270 K=1,NTA                                                            
      IF(K.EQ.NCAT)GO TO 270                                                    
      IF(TOPO(NAT,K).GT.0)GO TO 280                                             
  270 continue                                                                  
  280 NCAT2=K                                                                   
      T1=PVEC(NAT,NCAT,NCAT2,IPM)                                               
      GO TO 300                                                                 
  290 IF(I.EQ.J)T1=1.0D0                                                        
  300 continue                                                                  
      return                                                                    
      END                                                                       
