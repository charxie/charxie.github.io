C EXTENDED HUCKEL THEORY FOR CLOSED SHELLS 
C OVERLAPS ARE IN MATRIX A, COULOMB INTEGRALS (GAMMA) ARE IN MATRIX Q      

      SUBROUTINE HUCKCL  
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
      integer CHARGE,OCCA,OCCB,UL,AN,CZ,U,ULIM,ANI
      character option*4,opnclo*3
      integer HUCKEL,CNDO,INDO,CLOSED,OPEN 
      common/ARRAYS/A(NBSZR,NBSZR),B(NBSZR,NBSZR),D(NBSZR,NBSZR)
      common/INFO/NATOMS,CHARGE,MULTIP,AN(NAZR),C(NAZR,3),N,NTA,INO,IPR
      common/INFO1/CZ(NAZR),U(NBSZR),ULIM(NAZR),LLIM(NAZR),
     :NELECS,OCCA,OCCB
      common/GAB/XXX(5*NBSZR),G(NAZR,NAZR),
     :Q(NBSZR),YYY(NBSZR),ENERGY,XXY(214)
      common/OPSION/OPTION,OPNCLO,HUCKEL,CNDO,INDO,CLOSED,OPEN     
      dimension ENEG(18,3),BETA0(18)     
      dimension G1(18),F2(18)
      
      do 10 I=1,N 
      do 10 J=1,N 
      D(I,J)=A(I,J)  
   10 continue 
      G1(3)=.092012     D0  
      G1(4)=.1407       D0  
      G1(5)=.199265     D0  
      G1(6)=.267708     D0  
      G1(7)=.346029     D0  
      G1(8)=.43423      D0  
      G1(9)=.532305     D0  
      F2(3)=.049865     D0  
      F2(4)=.089125     D0  
      F2(5)=.13041      D0  
      F2(6)=.17372      D0  
      F2(7)=.219055     D0  
      F2(8)=.266415     D0  
      F2(9)=.31580      D0  
      ENEG(1,1)=7.1761  D0  
      ENEG(3,1)=3.1055  D0  
      ENEG(3,2)=1.258   D0  
      ENEG(4,1)=5.94557 D0  
      ENEG(4,2)=2.563   D0  
      ENEG(5,1)=9.59407 D0  
      ENEG(5,2)=4.001   D0  
      ENEG(6,1)=14.051  D0  
      ENEG(6,2)=5.572   D0  
      ENEG(7,1)=19.31637D0  
      ENEG(7,2)=7.275   D0  
      ENEG(8,1)=25.39017D0  
      ENEG(8,2)=9.111   D0  
      ENEG(9,1)=32.2724 D0  
      ENEG(9,2)=11.08   D0  
      ENEG(11,1)=2.804  D0  
      ENEG(11,2)=1.302  D0  
      ENEG(11,3)=0.150  D0  
      ENEG(12,1)=5.1254 D0  
      ENEG(12,2)=2.0516 D0  
      ENEG(12,3)=0.16195D0  
      ENEG(13,1)=7.7706 D0  
      ENEG(13,2)=2.9951 D0  
      ENEG(13,3)=0.22425D0  
      ENEG(14,1)=10.0327D0  
      ENEG(14,2)=4.1325 D0  
      ENEG(14,3)=0.337  D0  
      ENEG(15,1)=14.0327D0  
      ENEG(15,2)=5.4638 D0  
      ENEG(15,3)=0.500  D0  
      ENEG(16,1)=17.6496D0  
      ENEG(16,2)=6.989  D0  
      ENEG(16,3)=0.71325D0  
      ENEG(17,1)=21.5906D0  
      ENEG(17,2)=8.7081 D0  
      ENEG(17,3)=0.97695D0  
      BETA0(1)= -9.     D0  
      BETA0(3)= -9.     D0  
      BETA0(4)= -13.    D0  
      BETA0(5)= -17.    D0  
      BETA0(6)= -21.    D0  
      BETA0(7)= -25.    D0  
      BETA0(8)= -31.    D0  
      BETA0(9)= -39.    D0  
      BETA0(11)=-7.7203 D0  
      BETA0(12)=-9.4471 D0  
      BETA0(13)=-11.3011D0  
      BETA0(14)=-13.065 D0  
      BETA0(15)=-15.070 D0  
      BETA0(16)=-18.150 D0  
      BETA0(17)=-22.330 D0

C FIND NELECS AND FILL H CORE(DIAGONAL) WITH (I+A)/2 

      NELECS=0 
      DO 70 I=1,NATOMS      
      NELECS=NELECS+CZ(I)   
      LL =LLIM(I) 
      UL =ULIM(I) 
      ANI=AN(I)  
      L=0      
      DO 60 J=LL,UL  
      L=L+1    
      IF (L.EQ.1) GO TO 50  
   20 IF (L.LT.5) GO TO 40  
   30 A(J,J)=-ENEG(ANI,3)/27.21D0 
      GO TO 60 
   40 A(J,J)=-ENEG(ANI,2)/27.21D0 
      GO TO 60 
   50 A(J,J)=-ENEG(ANI,1)/27.21D0      
   60 continue 
   70 continue 
      NELECS=NELECS-CHARGE  
      OCCA=NELECS/2  

C FORM HUCKEL HAMILTONIAN IN A (OFF DIAGONAL TWO CENTER TERMS) 

      DO 100 I=2,N 
      K=U(I)   
      L=AN(K)  
      UL=I-1   
      DO 100 J=1,UL  
      KK=U(J)  
      LL=AN(KK)  
      IF ((L.GT.9).OR.(LL.GT.9)) GO TO 90  
   80 A(I,J)=A(I,J)*(BETA0(L)+BETA0(LL))/54.42D0      
      A(J,I)=A(I,J)  
      GO TO 100  
   90 A(I,J)=0.75D0*A(I,J)*(BETA0(L)+BETA0(LL))/54.42D0 
      A(J,I)=A(I,J)  
  100 continue 
      DO 110 I=1,N 
  110 Q(I)=A(I,I) 
      RHO=1.D-6  

      call EIGN(N,RHO)      

C EIGENVECTORS (IN B) ARE CONVERTED INTO DENSITY MATRIX (IN B) 

      DO 150 I=1,N 
      DO 130 J=I,N 
      XXX(J)=0.0D0 
      DO 120 K=1,OCCA       
  120 XXX(J)= XXX(J)+2.D0*B(I,K)*B(J,K)  
  130 continue 
      DO 140 J=I,N 
  140 B(I,J)= XXX(J) 
  150 continue 
      DO 160 I=1,N 
      DO 160 J=I,N 
  160 B(J,I)=B(I,J)  
C ADD V(AB) TO HCORE--CNDO 
      DO  180  I=1,N 
      J=U(I)   
      Q(I)=Q(I)  +0.5D0*G(J,J) 
      DO 170 K=1,NATOMS     
  170 Q(I)=Q(I)-dble(CZ(K))*G(J,K)     
  180 continue 
C EXIT SEGMENT IF ONLY CNDO APPROXIMATIONS ARE DESIRED  
      IF (OPTION.EQ.'CNDO') GO TO 300
      write(6,*) ' INDO *****'
C INDO MODIFICATION (CORRECTION TO U(I,I) )       
  190 DO 290 I=1,NATOMS     
      K=AN(I)  
      J=LLIM(I)  
      IF ((K.GT.1).AND.(K.LT.10)) GO TO 200 
      GO TO 290  
  200 IF (K.LE.3) GO TO 220 
  210 Q(J)=Q(J)   +(dble(CZ(I))-1.5D0)*G1(K)/6.D0   
  220 IF(K.EQ.3) GO TO 260  
  230 IF(K.EQ.4) GO TO 250  
  240 TEMP=G1(K)/3.D0+(dble(CZ(I))-2.5D0)*2.D0*F2(K)/25.D0       
      GO TO 270  
  250 TEMP=G1(K)/4.D0       
      GO TO 270  
  260 TEMP=G1(K)/12.D0      
  270 continue 
      DO 280 L=1,3 
  280 Q(J+L)=Q(J+L)+TEMP    
  290 continue 
  300 continue 
      DO 320 I=1,N 
      DO 310 J=I,N 
  310 A(J,I)=A(I,J)  
  320 A(I,I)=Q(I) 
      IF(IPR.LT.5) return   
      WRITE(6,1000)  
      call SCFOUT(0,1,0)    

      return   
 1000 FORMAT(1X,18H CORE HAMILTONIAN /)
      END