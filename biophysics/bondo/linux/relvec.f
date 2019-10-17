      SUBROUTINE RELVEC(R,E,C1,C2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION E(3),C1(3),C2(3)
      X = 0.D0
      DO 10 I=1,3
      E(I) = C2(I)-C1(I)
      X = X+E(I)**2
   10 continue
      R=DSQRT(X)
      DO 40 I=1,3
      IF (R.GT..000001D0) GO TO 30
   20 GO TO 40
   30 E(I)=E(I)/R
   40 continue
      return
      END