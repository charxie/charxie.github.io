      SUBROUTINE BINTGS(X,K)                                                    
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
C     FILLS ARRAY OF B-INTEGRALS. NOTE THAT B(I) IS B(I-1) IN THE               
C     USUAL NOTATION                                                            
C     FOR X.GT.3                    EXPONENTIAL FORMULA IS USED                 
C     FOR 2.LT.X.LE.3 AND K.LE.10   EXPONENTIAL FORMULA IS USED                 
C     FOR 2.LT.X.LE.3 AND K.GT.10   15 TERM SERIES IS USED                      
C     FOR 1.LT.X .E.2 AND K.LE.7    EXPONENTIAL FORMULA IS USED                 
C     FOR 1.LT.X.LE.2 AND K.GT.7    12 TERM SERIES IS USED                      
C     FOR .5.LT.X.LE.1 AND K.LE.5   EXPONENTIAL FORMULA IS USED                 
C     FOR .5.LT.X.LE.1 AND K.GT.5    7 TERM SERIES IS USED                      
C     FOR X.LE..5                    6 TERM SERIES IS USED                      
C     ******************************************************************        
      COMMON/AUXINT/A(17),B(17)
      IO=1                                                                      
      K1=K+1                                                                    
      ABSX=DABS(X)                                                              
      IF(ABSX.GT.3.D0) GO TO 120                                                
   10 IF(ABSX.GT.2.D0) GO TO 100                                                
   20 IF(ABSX.GT.1.D0) GO TO 80                                                 
   30 IF(ABSX.GT..5D0) GO TO 60                                                 
   40 IF(ABSX.GT..000001D0) GO TO 50                                            
      GO TO 170                                                                 
   50 LAST=6                                                                    
      GO TO 140                                                                 
   60 IF(K.LE.5) GO TO 120                                                      
   70 LAST=7                                                                    
      GO TO 140                                                                 
   80 IF(K.LE.7) GO TO 120                                                      
   90 LAST=12                                                                   
      GO TO 140                                                                 
  100 IF(K.LE.10) GO TO 120                                                     
  110 LAST=15                                                                   
      GO TO 140                                                                 
C                                                                               
  120 EXPX=DEXP(X)                                                              
      EXPMX=1.D0/EXPX                                                           
      B(1)=(EXPX-EXPMX)/X                                                       
      DO 130 I=1,K                                                              
  130 B(I+1)=(dble(I)*B(I)+(-1.D0)**I*EXPX-EXPMX)/X                           
      GO TO 190                                                                 
  140 LAST1=LAST+1                                                              
      DO 160 IL=IO,K1                                                           
      I=IL-1                                                                    
      Y=0.D0                                                                    
      DO 150 ML=IO,LAST1                                                        
      M=ML-1                                                                    
c      Y=Y+(-X)**M*(1.D0-(-1.D0)**(M+I+1))/(FACT(M)*dble(M+I+1))
      dfact=fact(m)
      Y=Y+(-X)**M*(1.D0-(-1.D0)**(M+I+1))/(dfact*dble(M+I+1))
  150 continue
  160 B(I+1)=Y
      GO TO 190                                                                 
C                                                                               
  170 DO 180 IL=IO,K1                                                           
      I=IL-1                                                                    
  180 B(I+1)=(1.D0-(-1.D0)**(I+1))/dble(I+1)                                  
  190 continue                                                                  
      K2=K1+1                                                                   
      DO 200 I=K2,17                                                            
      B(I)=0.0D0                                                                
  200 continue
  
      return                                                                    
      END