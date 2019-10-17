      SUBROUTINE correcteigen(tau,n,psi,time,nstp,
     :                        tmd,eg,cofeg,ft,cofft)

      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NMAX=500)
      COMPLEX*8 RX(NMAX),VX(NMAX),AX(NMAX),BX(NMAX),CX(NMAX),DX(NMAX)
      COMPLEX*8 PSI(N)
      REAL*4 TMD(NSTP),ft(N,N,NSTP),cofft(N,N,NSTP,3)
      real*4 eg(n,nstp), cofeg(n,nstp,3)
      COMPLEX*8 HXI,CORRX
      COMMON /BLOCK/ RX,VX,AX,BX,CX,DX

      if(n.gt.nmax) stop ' error: n > nmax ! '
      
      DO I = 1 , N
         RX(I)=PSI(I)
      ENDDO

c---> Gear corrector coefficients for a first-order equation
      
      GEAR0 =  95.0 / 288.0 
      GEAR1 =   1.0
      GEAR2 =  25.0 /  24.0
      GEAR3 =  35.0 /  72.0 
      GEAR4 =   5.0 /  48.0
      GEAR5 =   1.0 / 120.0

      C1 = TAU
      C2 = C1 * TAU / 2.0
      C3 = C2 * TAU / 3.0
      C4 = C3 * TAU / 4.0
      C5 = C4 * TAU / 5.0

      CR = GEAR0 * C1
      CV = GEAR1 * C1 / C1
      CA = GEAR2 * C1 / C2
      CB = GEAR3 * C1 / C3
      CC = GEAR4 * C1 / C4
      CD = GEAR5 * C1 / C5
      ITT=(TIME-TMD(1))/(TMD(2)-TMD(1))+1
      DTT=TIME-TMD(ITT) 

      DO I = 1 , N
	 HXI = (0.0,0.0)
	 DO J = 1 , N
	    if(i.eq.j) then
	       EIJ=((COFEG(I,ITT,3)*DTT+COFEG(I,ITT,2))*DTT
     :              +COFEG(I,ITT,1))*DTT+EG(I,ITT)
	    else 
	       EIJ=0.0
	    endif
	    HIJ=((COFFT(J,I,ITT,3)*DTT+COFFT(J,I,ITT,2))*DTT
     :           +COFFT(J,I,ITT,1))*DTT+FT(J,I,ITT)
c            HXI=HXI-((0.0,1.0)*EIJ+HIJ/1.51653)*RX(J)
            HXI=HXI-((0.0,1.0)*EIJ-HIJ)*RX(J)
	 ENDDO
         CORRX = HXI - VX(I)
         RX(I) = RX(I) + SNGL(CR) * CORRX
         VX(I) = VX(I) + SNGL(CV) * CORRX
         AX(I) = AX(I) + SNGL(CA) * CORRX
         BX(I) = BX(I) + SNGL(CB) * CORRX
         CX(I) = CX(I) + SNGL(CC) * CORRX
         DX(I) = DX(I) + SNGL(CD) * CORRX
      ENDDO
      DO I = 1 , N
         PSI(I)=RX(I)
      ENDDO

      RETURN
      END