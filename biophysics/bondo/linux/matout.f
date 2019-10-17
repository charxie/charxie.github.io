      SUBROUTINE MATOUT(N,MATOP)
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/ARRAYS/A(NBSZR,NBSZR,3) 
      DO 40 M=1,N,11
      K=M+10
      IF (K.LE.N) GO TO 20
   10 K=N
   20 continue
      WRITE(6,1000) (J,J=M,K)
      DO 30 I=1,N
      WRITE(6,1100) I,(A(I,J,MATOP),J=M,K)
   30 continue
      WRITE(6,1200)
   40 continue
      return
 1000 FORMAT(//,7X,11(4X,I2,3X),//)
 1100 FORMAT(1X,I2,4X,50(F9.4))
 1200 FORMAT(//)
      END