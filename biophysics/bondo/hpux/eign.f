C     RHO= UPPER LIMIT FOR OFF-DIAGONAL ELEMENT                                 
C     NN= SIZE OF MATRIX                                                        
C     A = F MATRIX (ONLY LOWER TRIANGLE IS USED + THIS IS DESTROYED)            
C     EIG = RETURNED EIGENVALUES IN ALGEBRAIC ASCENDING ORDER                   
C     VEC = RETURNED EIGENVECTORS IN COLUMNS                                    

      SUBROUTINE EIGN(NN,RHO) 
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      COMMON/ARRAYS/A(NBSZR,NBSZR),VEC(NBSZR,NBSZR),X(NBSZR,NBSZR)
      COMMON/GAB/GAMMA(NBSZR),BETA(NBSZR),BETASQ(NBSZR),EIG(NBSZR),
     :W(NBSZR),XYZ(NAIGAIO2-(5*NBSZR))
      DIMENSION P(NBSZR),Q(NBSZR) 
      EQUIVALENCE (P(1),BETA(1)),(Q(1),BETA(1))                                 
      DIMENSION IPOSV(NBSZR),IVPOS(NBSZR),IORD(NBSZR) 
      EQUIVALENCE (IPOSV(1),GAMMA(1)),(IVPOS(1),BETA(1)),(IORD(1),BETAS        
     :Q(1))                                                                     

      RHOSQ=RHO*RHO                                                             
      N=NN                                                                      
      IF (N .EQ. 0) GO TO 560                                                   
   10 N1=N-1                                                                    
      N2=N-2                                                                    
      GAMMA(1)=A(1,1)                                                           
      IF(N2) 180,170,20                                                         
   20 DO 160 NR=1,N2                                                            
      B=A(NR+1,NR)                                                              
      S=0.D0                                                                    
      DO 30 I=NR,N2                                                             
   30 S=S+A(I+2,NR)**2                                                          
C     PREPARE FOR POSSIBLE BYPASS OF TRANSFORMATION                             
      A(NR+1,NR)=0.D0                                                           
      IF (S) 150,150,40                                                         
   40 S=S+B*B                                                                   
      SGN=+1.D0                                                                 
      IF (B) 50,60,60                                                           
   50 SGN=-1.D0                                                                 
   60 SQRTS=DSQRT(S)                                                            
      D=SGN/(SQRTS+SQRTS)                                                       
      TEMP=DSQRT(.5D0+B*D)                                                      
      W(NR)=TEMP                                                                
      A(NR+1,NR)=TEMP                                                           
      D=D/TEMP                                                                  
      B=-SGN*SQRTS                                                              
C     D IS FACTOR OF PROPORTIONALITY. NOW COMPUTE AND SAVE W VECTOR.            
C     EXTRA SINGLY SUBSCRIPTED W VECTOR USED FOR SPEED.                         
      DO 70 I=NR,N2                                                             
      TEMP=D*A(I+2,NR)                                                          
      W(I+1)=TEMP                                                               
   70 A(I+2,NR)=TEMP                                                            
C     PREMULTIPLY VECTOR W BY MATRIX A TO OBTAIN P VECTOR.                      
C     SIMULTANEOUSLY ACCUMULATE DOT PRODUCT WP,(THE SCALAR K)                   
      WTAW=0.D0                                                                 
      DO 120 I=NR,N1                                                            
      SUM=0.D0                                                                  
      DO 80 J=NR,I                                                              
   80 SUM=SUM+A(I+1,J+1)*W(J)                                                   
      I1=I+1                                                                    
      IF(N1-I1) 110,90,90                                                       
   90 DO 100 J=I1,N1                                                            
  100 SUM=SUM+A(J+1,I+1)*W(J)                                                   
  110 P(I)=SUM                                                                  
  120 WTAW=WTAW+SUM*W(I)                                                        
C     P VECTOR AND SCALAR K  NOW STORED. NEXT COMPUTE Q VECTOR                  
      DO 130 I=NR,N1                                                            
  130 Q(I)=P(I)-WTAW*W(I)                                                       
C     NOW FORM PAP MATRIX, REQUIRED PART                                        
      DO 140 J=NR,N1                                                            
      QJ=Q(J)                                                                   
      WJ=W(J)                                                                   
      DO 140 I=J,N1                                                             
  140 A(I+1,J+1)=A(I+1,J+1)-2.D0*(W(I)*QJ+WJ*Q(I))                              
  150 BETA(NR)=B                                                                
      BETASQ(NR)=B*B                                                            
  160 GAMMA(NR+1)=A(NR+1,NR+1)                                                  
  170 B=A(N,N-1)                                                                
      BETA(N-1)=B                                                               
      BETASQ(N-1)=B*B                                                           
      GAMMA(N)=A(N,N)                                                           
  180 BETASQ(N)=0.D0                                                            
C     ADJOIN AN IDENTITY MATRIX TO BE POSTMULTIPLIED BY ROTATIONS.              
      DO 200 I=1,N                                                              
      DO 190 J=1,N                                                              
  190 VEC(I,J)=0.D0                                                             
  200 VEC(I,I)=1.D0                                                             
      M=N                                                                       
      SUM=0.D0                                                                  
      NPAS=1                                                                    
      GO TO 330                                                                 
  210 SUM=SUM+SHIFT                                                             
      COSA=1.D0                                                                 
      G=GAMMA(1)-SHIFT                                                          
      PP=G                                                                      
      PPBS=PP*PP+BETASQ(1)                                                      
      PPBR=DSQRT(PPBS)                                                          
      DO 300 J=1,M                                                              
      COSAP=COSA                                                                
      IF (PPBS.GT.1.D-12) GO TO 230                                             
  220 SINA=0.D0                                                                 
      SINA2=0.D0                                                                
      COSA=1.D0                                                                 
      GO TO 270                                                                 
  230 SINA=BETA(J)/PPBR                                                         
      SINA2=BETASQ(J)/PPBS                                                      
      COSA=PP/PPBR                                                              
C     POSTMULTIPLY IDENTITY BY P-TRANSPOSE MATRIX                               
      NT=J+NPAS                                                                 
      IF(NT .LT. N) GO TO 250                                                   
  240 NT=N                                                                      
  250 DO 260 I=1,NT                                                             
      TEMP=COSA*VEC(I,J)+SINA*VEC(I,J+1)                                        
      VEC(I,J+1)=-SINA*VEC(I,J)+COSA*VEC(I,J+1)                                 
  260 VEC(I,J)=TEMP                                                             
  270 DIA=GAMMA(J+1)-SHIFT                                                      
      U=SINA2*(G+DIA)                                                           
      GAMMA(J)=G+U                                                              
      G=DIA-U                                                                   
      PP=DIA*COSA-SINA*COSAP*BETA(J)                                            
      IF(J .NE. M) GO TO 290                                                    
  280 BETA(J)=SINA*PP                                                           
      BETASQ(J)=SINA2*PP*PP                                                     
      GO TO 310                                                                 
  290 PPBS=PP*PP+BETASQ(J+1)                                                    
      PPBR=DSQRT(PPBS)                                                          
      BETA(J)=SINA*PPBR                                                         
  300 BETASQ(J)=SINA2*PPBS                                                      
  310 GAMMA(M+1)=G                                                              
C     TEST FOR CONVERGENCE OF LAST DIAGONAL ELEMENT                             
      NPAS=NPAS+1                                                               
      IF(BETASQ(M) .GT. RHOSQ) GO TO 350                                        
  320 EIG(M+1)=GAMMA(M+1)+SUM                                                   
  330 BETA(M)=0.D0                                                              
      BETASQ(M)=0.D0                                                            
      M=M-1                                                                     
      IF(M .EQ. 0) GO TO 380                                                    
  340 IF(BETASQ(M) .LE. RHOSQ) GO TO 320                                        
C     TAKE ROOT OF CORNER 2 BY 2 NEAREST TO LOWER DIAGONAL IN VALUE             
C     AS ESTIMATE OF EIGENVALUE TO USE FOR SHIFT                                
  350 A2=GAMMA(M+1)                                                             
      R2=0.5D0*A2                                                               
      R1=0.5D0*GAMMA(M)                                                         
      R12=R1+R2                                                                 
      DIF=R1-R2                                                                 
      TEMP=DSQRT(DIF*DIF+BETASQ(M))                                             
      R1=R12+TEMP                                                               
      R2=R12-TEMP                                                               
      DIF=DABS(A2-R1)-DABS(A2-R2)                                               
      IF(DIF .LT. 0.D0) GO TO 370                                               
  360 SHIFT=R2                                                                  
      GO TO 210                                                                 
  370 SHIFT=R1                                                                  
      GO TO 210                                                                 
  380 EIG(1)=GAMMA(1)+SUM                                                       
C     INITIALIZE AUXILIARY TABLES REQUIRED FOR REARRANGING THE VECTORS          
      DO 390 J=1,N                                                              
      IPOSV(J)=J                                                                
      IVPOS(J)=J                                                                
  390 IORD(J)=J                                                                 
C     USE A TRANSPOSITION SORT TO ORDER THE EIGENVALUES                         
      M=N                                                                       
      GO TO 430                                                                 
  400 DO 420 J=1,M                                                              
      IF (EIG(J) .LE. EIG(J+1)) GO TO 420                                       
  410 TEMP=EIG(J)                                                               
      EIG(J)=EIG(J+1)                                                           
      EIG(J+1)=TEMP                                                             
      ITEMP=IORD(J)                                                             
      IORD(J)=IORD(J+1)                                                         
      IORD(J+1)=ITEMP                                                           
  420 CONTINUE                                                                  
  430 M=M-1                                                                     
      IF(M .NE. 0) GO TO 400                                                    
  440 IF(N1 .EQ. 0) GO TO 490                                                   
  450 DO 480 L=1,N1                                                             
      NV=IORD(L)                                                                
      NP=IPOSV(NV)                                                              
      IF(NP .EQ. L) GO TO 480                                                   
  460 LV=IVPOS(L)                                                               
      IVPOS(NP)=LV                                                              
      IPOSV(LV)=NP                                                              
      DO 470 I=1,N                                                              
      TEMP=VEC(I,L)                                                             
      VEC(I,L)=VEC(I,NP)                                                        
  470 VEC(I,NP)=TEMP                                                            
  480 CONTINUE                                                                  
  490 CONTINUE                                                                  
C     BACK TRANSFORM THE VECTORS OF THE TRIPLE DIAGONAL MATRIX                  
      DO 550 NRR=1,N                                                            
      K=N1                                                                      
  500 K=K-1                                                                     
      IF(K .LE. 0) GO TO 540                                                    
  510 SUM=0.D0                                                                  
      DO 520 I=K,N1                                                             
  520 SUM=SUM+VEC(I+1,NRR)*A(I+1,K)                                             
      SUM=SUM+SUM                                                               
      DO 530 I=K,N1                                                             
  530 VEC(I+1,NRR)=VEC(I+1,NRR)-SUM*A(I+1,K)                                    
      GO TO 500                                                                 
  540 CONTINUE                                                                  
  550 CONTINUE                                                                  

  560 RETURN                                                                    
      END                                                                       
