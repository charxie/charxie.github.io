C PROCEDURE FOR CALCULATING REDUCED OVERLAP INTEGRALS
      FUNCTION SS(NN1,LL1,MM,NN2,LL2,ALPHA,BETA)
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER ULIM,ULIM1
      COMMON/ARRAYS/S(NBSZR,NBSZR),Y(9,5,203),Z(17,45),XX(NAIGAIO1)
      COMMON/AUXINT/A(17),B(17)
      
      N1=NN1
      L1=LL1
      M=MM
      N2=NN2
      L2=LL2
      P =(ALPHA + BETA)/2.D0
      PT=(ALPHA - BETA)/2.D0
      X = 0.D0
      M=IABS(M)
C REVERSE QUANTUM NUMBERS IF NECESSARY
      IF((L2.LT.L1).OR.((L2.EQ.L1).AND.(N2.LT.N1))) GO TO 20
   10 GO TO 30
   20 K = N1
      N1= N2
      N2= K
      K= L1
      L1= L2
      L2= K
      PT=-PT
   30 CONTINUE
      K = MOD((N1+N2-L1-L2),2)
C FIND A AND B INTEGRALS
      CALL AINTGS(P,N1+N2)
      CALL BINTGS(PT,N1+N2)
      IF((L1.GT.0).OR.(L2.GT.0)) GO TO 60
C BEGIN SECTION USED FOR OVERLAP INTEGRALS INVOLVING S FUNCTIONS
C FIND Z TABLE NUMBER L
   40 L = (90-17*N1+N1**2-2*N2)/2
      ULIM = N1+N2
      LLIM =0
      ULIM1=ULIM+1
      LLIM1=LLIM+1
      DO 50 I=LLIM1,ULIM1
      IL=I-1
      NNI1=N1+N2-IL+1
   50 X=X+Z(IL+1,L)*A(IL+1)*B(NNI1)/2.D0
      SS=X
      GO TO 80
C BEGIN SECTION USED FOR OVERLAPS INVOLVING NON-S FUNCTIONS
C FIND Y TABLE NUMBER L
   60 L=(5-M)*(24-10*M+M**2)*(83-30*M+3*M**2)/120+
     :  (30-9*L1+L1**2-2*N1)*(28-9*L1+L1**2-2*N1)/8+
     :  (30-9*L2+L2**2-2*N2)/2
      LLIM=0
      LLIM1=LLIM+1
      DO 70 I=LLIM1,9
      IL=I-1
      ULIM=4-MOD(K+IL,2)
      ULIM1=ULIM+1
      DO 70 J=LLIM1,ULIM1
      JL=J-1
      IIII=2*JL+MOD(K+IL,2)+1
   70 X=X+Y(IL+1,JL+1,L)*A(IL+1)*B(IIII)
      SS=X*(FACT(M+1)/8.D0)**2*DSQRT(dble(2*L1+1)*FACT(L1-M)*
     :dble(2*L2+1)*FACT(L2-M)/(4.D0*FACT(L1+M)*FACT(L2+M)))
   80 CONTINUE
   
      return
      END