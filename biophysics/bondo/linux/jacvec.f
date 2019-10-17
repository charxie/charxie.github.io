      SUBROUTINE JACVEC(N,A,EIVU,EIVR,NDIM)                                     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      DIMENSION A(NDIM,1),EIVR(NDIM,1),EIVU(1)                                  
      IF(N.GT.1) GO TO 10                                                       
      EIVR(1,1)=1.0D0                                                           
      EIVU(1)=A(1,1)                                                            
      return                                                                    
   10 continue                                                                  
      DO 30 J=1,N                                                               
      DO 20 I=1,N                                                               
   20 EIVR(I,J)=0.0D0                                                           
   30 EIVR(J,J)=1.0D0                                                           
C                                                                               
C FIND THE ABSOLUTELY LARGEST ELEMENT OF A                               
C                                                                               
   40 ATOP=0.0D0                                                                
      DO 60 J=1,N                                                               
      DO 60 I=1,N                                                               
      IF (ATOP- DABS(A(I,J))) 50,60,60                                          
   50 ATOP= DABS(A(I,J))                                                        
   60 continue                                                                  
C                                                                               
C CALCULATE THE STOPPING CRITERION -- DSTOP                              
C                                                                               
   70 AVGF= dble(N*(N-1))*.55D0                                               
      D=0.0D0                                                                   
      DO 80 JJ=2,N                                                              
      DO 80 II=2,JJ                                                             
      S=A(II-1,JJ)/ATOP                                                         
   80 D=S*S+D                                                                   
      DSTOP=(1.D-7)*D                                                           
C                                                                               
C        CALCULATE THE THRESHOLD, THRSH                                         
C                                                                               
      THRSH= DSQRT(D/AVGF)*ATOP                                                 
C                                                                               
C        START A SWEEP                                                          
C                                                                               
   90 IFLAG=0                                                                   
      DO 250 JCOL=2,N                                                           
      JCOL1=JCOL-1                                                              
      DO 250 IROW=1,JCOL1                                                       
      AIJ=A(IROW,JCOL)                                                          
C                                                                               
C        COMPARE THE OFF-DIAGONAL ELEMENT WITH THRSH                            
C                                                                               
      IF ( DABS(AIJ)-THRSH) 250,250,100                                         
  100 AII=A(IROW,IROW)                                                          
      AJJ=A(JCOL,JCOL)                                                          
      S=AJJ-AII                                                                 
C                                                                               
C        CHECK TO SEE IF THE CHOSEN ROTATION IS LESS THAN THE ROUNDING E        
C        IF SO , THEN DO NOT ROTATE.                                            
C                                                                               
      IF ( DABS(AIJ)-1.D-09*DABS(S)) 250,250,110                                
  110 IFLAG=1                                                                   
C                                                                               
C        IF THE ROTATION IS VERY CLOSE TO 45 DEGREES, SET SIN AND COS           
C        TO 1/(ROOT 2).                                                         
C                                                                               
      IF (1.D-10*DABS(AIJ)-DABS(S)) 130,120,120                                 
  120 S=.707106781D0                                                            
      C=S                                                                       
      GO TO 140                                                                 
C                                                                               
C        CALCULATION OF SIN AND COS FOR ROTATION THAT IS NOT VERY CLOSE         
C        TO 45 DEGREES                                                          
C                                                                               
  130 T=AIJ/S                                                                   
      S=0.25D0/ DSQRT(0.25D0+T*T)                                               
C                                                                               
C        COS=C ,  SIN=S                                                         
C                                                                               
      C= DSQRT(0.5D0+S)                                                         
      S=2.D0*T*S/C                                                              
C                                                                               
C        CALCULATION OF THE NEW ELEMENTS OF MATRIX A                            
C                                                                               
  140 DO 150 I=1,IROW                                                           
      T=A(I,IROW)                                                               
      U=A(I,JCOL)                                                               
      A(I,IROW)=C*T-S*U                                                         
  150 A(I,JCOL)=S*T+C*U                                                         
      I2=IROW+2                                                                 
      IF (I2-JCOL) 160,160,180                                                  
  160 continue                                                                  
      DO 170 I=I2,JCOL                                                          
      T=A(I-1,JCOL)                                                             
      U=A(IROW,I-1)                                                             
      A(I-1,JCOL)=S*U+C*T                                                       
  170 A(IROW,I-1)=C*U-S*T                                                       
  180 A(JCOL,JCOL)=S*AIJ+C*AJJ                                                  
      A(IROW,IROW)=C*A(IROW,IROW)-S*(C*AIJ-S*AJJ)                               
      DO 190 J=JCOL,N                                                           
      T=A(IROW,J)                                                               
      U=A(JCOL,J)                                                               
      A(IROW,J)=C*T-S*U                                                         
  190 A(JCOL,J)=S*T+C*U                                                         
C                                                                               
C        ROTATION COMPLETED                                                     
C                                                                               
  200 DO 210 I=1,N                                                              
      T=EIVR(I,IROW)                                                            
      EIVR(I,IROW)=C*T-EIVR(I,JCOL)*S                                           
  210 EIVR(I,JCOL)=S*T+EIVR(I,JCOL)*C                                           
C                                                                               
C        CALCULATE THE NEW NORM D AND COMPARE WITH DSTOP                        
C                                                                               
      S=AIJ/ATOP                                                                
      D=D-S*S                                                                   
      IF (D-DSTOP) 220,240,240                                                  
C                                                                               
C        RECALCULATE DSTOP AND THRSH TO DISCARD ROUNDING ERRORS                 
C                                                                               
  220 D=0.0                                                                     
      DO 230 JJ=2,N                                                             
      DO 230 II=2,JJ                                                            
      S=A(II-1,JJ)/ATOP                                                         
  230 D=S*S+D                                                                   
      DSTOP=(1.D-7)*D                                                           
  240 THRSH= DSQRT(D/AVGF)*ATOP                                                 
  250 continue                                                                  
      IF(IFLAG) 90,260,90                                                       
C                                                                               
C        PLACE EIGENVALUES IN EIVU                                              
C                                                                               
  260 DO 270 J=1,N                                                              
      EIVU(J)=A(J,J)                                                            
  270 continue                                                                  
      return
      END