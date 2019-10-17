      FUNCTION PVEC(N1,N2,N3,I) 
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      INTEGER CH,AN                                                             
      COMMON/INFO/NA,CH,MULT,AN(NAZR),C(NAZR,3),NO,NTA,N,IPR

      DIMENSION V12(3),V13(3)                                                   
      R2=0.0D0                                                                  
      R3=0.0D0                                                                  
      DO 10 J=1,3                                                               
      V12(J)=C(N2,J)-C(N1,J)                                                    
      V13(J)=C(N3,J)-C(N1,J)                                                    
      R2=R2+V12(J)**2                                                           
      R3=R3+V13(J)**2                                                           
   10 continue                                                                  
      R2=DSQRT(R2)                                                              
      R3=DSQRT(R3)                                                              
      DOT=0.0D0                                                                 
      DO 20 J=1,3                                                               
      V12(J)=V12(J)/R2                                                          
      V13(J)=V13(J)/R3                                                          
      DOT=DOT+V12(J)*V13(J)                                                     
   20 continue                                                                  
      ARG=1.0D0-DOT**2                                                          
      SINT=DSQRT(ARG)                                                           
      IF(SINT.LE.1.D-4) GO TO 30                                                
      IF(I.EQ.1)PVEC=(V12(2)*V13(3)-V12(3)*V13(2))/SINT                         
      IF(I.EQ.2)PVEC=(V12(3)*V13(1)-V12(1)*V13(3))/SINT                         
      IF(I.EQ.3)PVEC=(V12(1)*V13(2)-V12(2)*V13(1))/SINT                         
      return                                                                    
C  VECTORS NUMERICALLY COLLINEAR; CHOOSE ARBITRARY PERP. DIRECTION              
   30 continue                                                                  
      DO 40 J=1,3                                                               
   40 V13(J)=0.0D0                                                              
      DO 60 J=1,3                                                               
      UV=V12(J)                                                                 
      IF((1.0D0-DABS(UV)).LT.1.D-3) GO TO 60                                    
      V13(J)=1.0D0                                                              
      R=0.0D0                                                                   
      DO 50 K=1,3                                                               
      V13(K)=V13(K)-UV*V12(K)                                                   
   50 R=R+V13(K)**2                                                             
      PVEC=V13(I)/DSQRT(R)                                                      
      return                                                                    
   60 continue                                                                  
      return                                                                    
      END