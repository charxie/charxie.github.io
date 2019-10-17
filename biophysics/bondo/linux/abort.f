      SUBROUTINE abortBondo(WORD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*6 word
      PRINT 1000,WORD
      call EXIT
      return
 1000 FORMAT(//,'***  ERROR EXIT0 ',A6,'  ***')
      END
