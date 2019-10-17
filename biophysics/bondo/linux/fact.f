      FUNCTION FACT(N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PRODT = 1.D0
      IF(N.EQ.0)GO TO 30
   10 DO 20 I=1,N
   20 PRODT=PRODT*dble(I)
   30 FACT=PRODT
      return
      END