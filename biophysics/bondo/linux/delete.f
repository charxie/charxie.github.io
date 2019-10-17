C  DELETES ROW AND COLUMN I OF FOCK MATRIX...                                   
C                                                                               
      SUBROUTINE DELETE(N,NLIST,LIST)  
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      DIMENSION LIST(40)                                                        
      COMMON/ARRAYS/A(NBSZR,NBSZR),B(NBSZR,NBSZR,2) 

      IF(NLIST.LE.0) return                                                     
      DO 30 I=1,NLIST                                                           
      K=LIST(I)                                                                 
      IF((K.LE.0).OR.(K.GT.N)) GO TO 40                                         
      DO 10 J=K,N                                                               
   10 A(J,K)=0.0                                                                
      DO 20 J=1,K                                                               
   20 A(K,J)=0.0                                                                
      A(K,K)=5.0                                                                
   30 continue                                                                  
      return                                                                    
   40 WRITE(6,1000) NLIST,(LIST(I),I=1,NLIST)                                   
      call EXIT                                                                 

      return                                                                    
 1000 FORMAT(1X,'EXIT...ILLEGAL DELETE LIST(',I4,')0',20I3)                     
      END