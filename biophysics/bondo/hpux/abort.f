      SUBROUTINE ABORT(WORD)                                                    
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      PRINT 1000,WORD                                                           
      CALL EXIT                                                                 
      RETURN                                                                    
 1000 FORMAT(//,'***  ERROR EXIT0 ',A6,'  ***')                                 
      END                                                                       
