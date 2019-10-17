C FETCHES LOWER TRIANGLE OF A(I,J,MOP) IF MOP=1 OR 3.

      FUNCTION ALT(I,J,MOP)
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/ARRAYS/A(NBSZR,NBSZR,3)
      ALT=A(I,J,MOP)
      IF(MOP.EQ.2) return
      IF(I.GT.J) return
      ALT=A(J,I,MOP)
      return
      END