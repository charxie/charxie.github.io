C  PRINTS OUT TABLE OF INTERNUCLEAR DISTANCES TO FILE = LFN...                  

      SUBROUTINE PRDIST(LFN)
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                        
      INTEGER AN                                                                
      COMMON/ARRAYS/R(NBSZR,NBSZR),BC(2*NBSZRS)
      COMMON/INFO/NATOMS,ICH,MULT,AN(NAZR),C(NAZR,3),NUDM,NTOT,INO,IPR 
      DIMENSION ISYM(10)                                                        
      DATA DASH/'----'/                                                         
      DATA ISYM/'GH',' H','HE','LI','BE',' B',' C',' N',' O',' F'/              

      WRITE(LFN,1000)                                                           
      NL=1                                                                      
      DO 10 I=2,NTOT                                                            
      I1=I-1                                                                    
      DO 10 J=1,I1                                                              
      X=C(J,1)-C(I,1)                                                           
      Y=C(J,2)-C(I,2)                                                           
      Z=C(J,3)-C(I,3)                                                           
   10 R(I,J)=DSQRT(X*X+Y*Y+Z*Z)                                                 
      NU=NTOT-1                                                                 
   20 IF((NU-NL).GT.14) NU=NL+14                                                
      WRITE(LFN,1100) (ISYM(AN(N)+1),N,N=NL,NU)                                 
      WRITE(LFN,1200) (DASH,DASH,N=NL,NU)                                       
      DO 30 N=NL,NTOT                                                           
      NM1=N-1                                                                   
      IF(N.GT.NU) NM1=NU                                                        
      WRITE(LFN,1300) ISYM(AN(N)+1),N,(R(N,M),M=NL,NM1)  
   30 CONTINUE 
      NL=NU+1                                                                   
      NU=NTOT-1                                                                 
      IF(NU.GE.NL) GO TO 20                                                     

      RETURN                                                                    
 1000 FORMAT(//,5X,'TABLE OF INTERNUCLEAR DISTANCES (ANGSTROMS)0')              
 1100 FORMAT(/,11X,16(A2,I2,4X))                                                
 1200 FORMAT(2X,'------',16(A4,A4))                                             
 1300 FORMAT(2X,A2,I2,2X,16F8.4)                                                
      END                                                                       
