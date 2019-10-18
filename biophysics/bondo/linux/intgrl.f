C     ATOMIC INTEGRALS FOR CNDO CALCULATIONS
      SUBROUTINE INTGRL
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      double precision MU,NUM,K1,K2
      integer AN,ULIM,ULK,ULL,CZ,U,CHARGE,ANL,ANK,OCCA,OCCB,TOPO
      character option*4,opnclo*3
      integer HUCKEL,CNDO,INDO,CLOSED,OPEN
      common/ARRAYS/S(NBSZR,NBSZR),Y(9,5,203),Z(17,45),XX(NAIGAIO1)
      common/INFO/NATOMS,CHARGE,MULTIP,AN(NAZR),C(NAZR,3),N,NTA,INO,IPR
      common/INFO1/CZ(NAZR),U(NBSZR),ULIM(NAZR),LLIM(NAZR),
     :NELECS,OCCA,OCCB
      common/GAB/XXX(5*NBSZR),GAMMA(NAZR,NAZR),
     :T(9,9),PAIRS(9,9),TEMP(9,9),C1(3),C2(3),
     :YYY(NAIGAIO2-(5*NBSZR)-(NAZR*NAZR)-249)
      common/AUXINT/A(17),B(17)
      common/OPSION/OPTION,OPNCLO,HUCKEL,CNDO,INDO,CLOSED,OPEN
      common/EXTRA/TT(NBSZR,NBSZR),IHYB(NAZR),
     :TOPO(NAZR,NAZR),T2J(NBSZR)
      dimension MU(18),NC(18),LC(9),MC(9),E(3)
      dimension P(NBSZR,NBSZR)
      EQUIVALENCE (P(1,1),Y(1,1,1))

C     DETERMINATION OF SIZE OF AO BASIS IN AND CORE CHARGE CZ
      N=0
      DO 60 I=1,NATOMS
      LLIM(I) = N+1
      K=1
      IF (AN(I).LT.11) GO TO 20
   10 N=N+9
      CZ(I)=AN(I)-10
      GO TO 50
   20 IF (AN(I).LT.3) GO TO 40
   30 N=N+4
      CZ(I) = AN(I)-2
      GO TO 50
   40 N=N+1
      CZ(I)= AN(I)
   50 CONTINUE
      ULIM(I) = N
   60 CONTINUE
      INO=N
C     FILL U ARRAY---U(J) IDENTIFIES THE ATOM TO WHICH ORBITAL J IS
C     ATTACHED E.G. ORBITAL 32 ATTACHED TO  ATOM 7, ETC.
      DO 70 K=1,NATOMS
      LLK = LLIM(K)
      ULK = ULIM(K)
      LIM = ULK+1-LLK
      DO 70 I=1,LIM
      J = LLK+I-1
   70 U(J) = K
C     ASSIGNMENT OF ORBITAL EXPONENTS TO ATOMS BY SLATERS RULES
      MU(2)=1.7D0
      MU(1)=1.2D0
      NC(1)=1
      NC(2)=1
      DO 80 I=3,10
      NC(I)=2
   80 MU(I)=.325D0*dble(I-1)
      DO 90 I=11,18
      NC(I)=3
   90 MU(I)=(.65D0*dble(I)-4.95D0)/3.D0
C     ASSIGNMENT OF ANGULAR MOMENTUM QUANTUM NOS. TO ATOMIC ORBITALS
      LC(1)=0
      LC(2)=1
      LC(3)=1
      LC(4)=1
      LC(5)=2
      LC(6)=2
      LC(7)=2
      LC(8)=2
      LC(9)=2
      MC(1)=0
      MC(2)=1
      MC(3)=-1
      MC(4)=0
      MC(5)=0
      MC(6)=1
      MC(7)=-1
      MC(8)=2
      MC(9)=-2
C     STEP THRU PAIRS OF ATOMS
      DO 400 K=1,NATOMS
      DO 400 L=K,NATOMS
      DO 100 I=1,3
      C1(I) = C(K,I)
  100 C2(I) = C(L,I)
C     CALCULATE UNIT VECTOR ALONG INTERATOM AXIS,E
      call RELVEC(R,E,C1,C2)
      LLK = LLIM(K)
      LLL = LLIM(L)
      ULK = ULIM(K)
      ULL = ULIM(L)
      NORBK=ULK-LLK+1
      NORBL=ULL-LLL+1
      NORBL1=NORBL
      ANK=AN(K)
      ANL=AN(L)
C     LOOP THRU PAIRS OF BASIS FUNCTIONS, ONE ON EACH ATOM
      DO 200 I=1,NORBK
      DO 200 J=1,NORBL
      IF(K.EQ.L) GO TO 160
  110 IF(MC(I).NE.MC(J)) GO TO 150
  120 IF(MC(I).LT.0) GO TO 140
  130 PAIRS(I,J)=DSQRT((MU(ANK)*R)**(2*NC(ANK)+1)*(MU(ANL)*R)**
     :(2*NC(ANL)+1)/(FACT(2*NC(ANK))*FACT(2*NC(ANL))))*(-1.D0)**
     :(LC(J)+MC(J))*
     :SS(NC(ANK),LC(I),MC(I),NC(ANL),LC(J),MU(ANK)*R,MU(ANL)*R)
     
c      write(6,*) i,j,k,l
c      if(an(k).eq.1.and.an(l).eq.6) then
c      write(6,'(10f6.3)') (b(ii),ii=1,17)
c      write(6,*)
c     :ss(NC(ANK),LC(I),MC(I),NC(ANL),LC(J),MU(ANK)*R,MU(ANL)*R)
c      endif
      
      GO TO 190
  140 PAIRS(I,J)=PAIRS(I-1,J-1)
      GO TO 190
  150 PAIRS(I,J)=0.0D0
      GO TO 190
  160 IF (I.EQ.J) GO TO 180
  170 PAIRS(I,J)=0.0D0
      GO TO 190
  180 PAIRS(I,J)=1.0D0
  190 CONTINUE
  200 CONTINUE
  
c      if(an(k).eq.1.and.an(l).eq.6) then
c        write(6,*) r,mu(ank),mu(anl),nc(ank),nc(anl)
c        do i = 1 , norbk
c        write(6,'(10f7.3)') (pairs(i,j),j=1,norbl)
c        enddo
c        stop
c      endif
  
  240 CONTINUE
      LCULK=LC(NORBK)
      LCULL=LC(NORBL1)
      MAXL=MAX0(LCULK,LCULL)
      IF(R.GT.0.000001D0) GO TO 260
  250 GO TO 290
C ROTATE INTEGRALS FROM DIATOMIC BASIS TO MOLECULAR BASIS
  260 call HARMTR(T,MAXL,E)
Csatyam (insert below file intgrlnew form satyam)
      FACTOR=0.585
C     WRITE(6,*)'FACTOR = ',FACTOR
      PAIRS(2,2)=FACTOR*PAIRS(2,2)
      PAIRS(2,3)=FACTOR*PAIRS(2,3)
      PAIRS(3,3)=FACTOR*PAIRS(3,3)
      PAIRS(3,2)=FACTOR*PAIRS(3,2)
Csatyam
      DO 270 I=1,NORBK
      DO 270 J=1,NORBL
      TEMP(I,J) = 0.D0
      DO 270 KK=1,NORBL
      TEMP(I,J) = TEMP(I,J)+T(J,KK)*PAIRS(I,KK)
  270 CONTINUE
      DO 280 I=1,NORBK
      DO 280 J=1,NORBL
      PAIRS(I,J) = 0.D0
      DO 280 KK=1,NORBK
      PAIRS(I,J) = PAIRS(I,J)+T(I,KK)*TEMP(KK,J)
  280 CONTINUE
C     FILL S MATRIX
  290 CONTINUE
      DO 300 I=1,NORBK
      LLKP=LLK+I-1
      DO 300 J=1,NORBL
      LLLP=LLL+J-1
  300 S(LLKP,LLLP)=PAIRS(I,J)
  340 CONTINUE
COMPUTATION OF 1-CENTER COULOMB INTEGRALS OVER SLATER S FUNCTIONS
      N1=NC(ANK)
      N2=NC(ANL)
      K1=MU(ANK)
      K2=MU(ANL)
      IF(K.NE.L) GO TO 370
  350 TERM1 = FACT(2*N1-1)/((2.D0*K2)**(2*N1))
      TERM2 = 0.D0
      LIM = 2*N1
      DO 360 J=1,LIM
      NUM =dble(J)*(2.D0*K1)**(2*N1-J)*FACT(4*N1-J-1)
      DEN = FACT(2*N1-J)*2.D0*dble(N1)*(2.D0*(K1+K2))**(4*N1-J)
      TERM2 = TERM2 + NUM/DEN
  360 CONTINUE
      GO TO 390
COMPUTATION OF 2-CENTER COULOMB INTEGRALS OVER SLATER S FUNCTIONS
  370 TERM1=(R/2.D0)**(2*N2)*SS(0,0,0,2*N2-1,0,0.D0,2.D0*K2*R)
      TERM2 = 0.D0
      LIM = 2*N1
      DO 380 J=1,LIM
  380 TERM2=TERM2+(dble(J)*(2.D0*K1)**(2*N1-J)*(R/2.D0)**(2*N1-J+2*
     :N2))/(FACT(2*N1-J)*2.D0*dble(N1))*SS(2*N1-J,0,0,2*N2-1,0,
     :2.D0*K1*R,2.D0*K2*R)
  390 GAMMA(K,L) = ((2.D0*K2)**(2*N2+1)/FACT(2*N2))*(TERM1-TERM2)
  400 CONTINUE
C SYMMETRIZATION OF OVERLAP AND COULOMB INTEGRAL MATRICES
      DO 410 I=1,N
      DO 410 J=I,N
  410 S(J,I) = S(I,J)
      DO 420 I=1,NATOMS
      DO 420 J=I,NATOMS
  420 GAMMA(J,I) = GAMMA(I,J)
      IF(IPR.LT.5) GO TO 430
      WRITE(6,1000)
      call MATOUT(N,1)
C TRANSFER GAMMA TO 80X80 MATRIX P FOR PRINTING
  430 CONTINUE
      DO 440 I=1,NATOMS
      DO 440 J=1,NATOMS
  440 P(I,J)=GAMMA(I,J)

c      open(37,file='tempq')
c      do j = 1 , n
c         write(37,'(10f8.3)')(s(i,j),i=1,n)
c      enddo
c      close(37)
c      stop

      IF(IPR.LT.5) RETURN
      WRITE(6,1100)
      call MATOUT(NATOMS,2)

      RETURN
 1000 FORMAT(//1X,23HOVERLAP INTEGRAL MATRIX)
 1100 FORMAT(1X,23HCOULOMB INTEGRAL MATRIX)
      END