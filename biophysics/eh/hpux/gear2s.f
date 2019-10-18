c Gear's five-value predictor-corrector method for 
c solving the time-dependent Schroedinger equation
 

      SUBROUTINE predictor2s(tau,n,psi)

c  PREDICTOR ROUTINE 

      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NMAX=2)
      COMPLEX*8 RX(NMAX),VX(NMAX),AX(NMAX),BX(NMAX),CX(NMAX),DX(NMAX)
      COMPLEX*8 PSI(N)
      COMMON /BLOCK2S/ RX,VX,AX,BX,CX,DX

      if(n.gt.nmax) stop ' error: n > nmax ! '

      DO I = 1 , N
         RX(I)=PSI(I)
      ENDDO

      C1 = TAU
      C2 = C1 * TAU / 2.0
      C3 = C2 * TAU / 3.0
      C4 = C3 * TAU / 4.0
      C5 = C4 * TAU / 5.0	      

      DO I = 1 , N

         RX(I)=RX(I)+SNGL(C1)*VX(I)+SNGL(C2)*AX(I)
     : 	            +SNGL(C3)*BX(I)+SNGL(C4)*CX(I)
     :              +SNGL(C5)*DX(I)
         VX(I)=VX(I)+SNGL(C1)*AX(I)+SNGL(C2)*BX(I)
     :	            +SNGL(C3)*CX(I)+SNGL(C4)*DX(I)
         AX(I)=AX(I)+SNGL(C1)*BX(I)+SNGL(C2)*CX(I)
     :              +SNGL(C3)*DX(I)    
         BX(I)=BX(I)+SNGL(C1)*CX(I)+SNGL(C2)*DX(I)
         CX(I)=CX(I)+SNGL(C1)*DX(I)

      ENDDO

      DO I = 1 , N
         PSI(I)=RX(I)
      ENDDO

      RETURN
      END

      SUBROUTINE 
     :corrector2s(tau,n,psi,time,nstp,tmd,vtr,vti,cvr,cvi)

c CORRECTOR ROUTINE

      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NMAX=2,NSMP=10)
      COMPLEX*8 RX(NMAX),VX(NMAX),AX(NMAX),BX(NMAX),CX(NMAX),DX(NMAX)
      COMPLEX*8 PSI(N)
      REAL*4 TMD(NSTP),VTR(N,N,NSTP),VTI(N,N,NSTP)
      REAL*4 CVR(N,N,NSTP,3),CVI(N,N,NSTP,3)
      COMPLEX*8 HXI,CORRX
      COMMON /BLOCK2S/ RX,VX,AX,BX,CX,DX

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
      CS = GEAR1 * C1 / C1
      CA = GEAR2 * C1 / C2
      CB = GEAR3 * C1 / C3
      CC = GEAR4 * C1 / C4
      CD = GEAR5 * C1 / C5

      ITT=(TIME-TMD(1))/(TMD(2)-TMD(1))+1
      DTT=TIME-TMD(ITT) 

      DO I = 1 , N
	 HXI = (0.0,0.0)
	 DO J = 1 , N
	    VIJR=((CVR(I,J,ITT,3)*DTT+CVR(I,J,ITT,2))*DTT
     :            +CVR(I,J,ITT,1))*DTT+VTR(I,J,ITT)
	    VIJI=((CVI(I,J,ITT,3)*DTT+CVI(I,J,ITT,2))*DTT
     :            +CVI(I,J,ITT,1))*DTT+VTI(I,J,ITT)
            HXI=HXI-((0.0,1.0)*VIJR+VIJI)*RX(J)
	 ENDDO
         CORRX = HXI - VX(I)
         RX(I) = RX(I) + SNGL(CR) * CORRX
         VX(I) = VX(I) + SNGL(CS) * CORRX
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
      
      
      SUBROUTINE save2s(n,time)

c SAVE ROUTINE

      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NMAX=2)
      COMPLEX*8 RX(NMAX),VX(NMAX),AX(NMAX),BX(NMAX),CX(NMAX),DX(NMAX)
      COMMON /BLOCK2S/ RX,VX,AX,BX,CX,DX

      OPEN(99,FILE='SAVE2S.CONT',STATUS='UNKNOWN',FORM='UNFORMATTED')

      REWIND 99

      WRITE(99) TIME,N
      DO I = 1 , N
	 WRITE(99) RX(I),VX(I),AX(I),BX(I),CX(I),DX(I)
      ENDDO

      CLOSE(99)	

      RETURN
      END      
      
      SUBROUTINE read2s(n,time0,psi0)

c READ ROUTINE

      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NMAX=2)
      COMPLEX*8 RX(NMAX),VX(NMAX),AX(NMAX),BX(NMAX),CX(NMAX),DX(NMAX)
      COMMON /BLOCK2S/ RX,VX,AX,BX,CX,DX
      COMPLEX*8 PSI0(N)

      OPEN(99,FILE='SAVE2S.CONT',STATUS='OLD',FORM='UNFORMATTED')

      READ(99) TIME0,N0
      
      IF(N0.NE.N) THEN
         WRITE(*,*) N0,N
         STOP ' ERROR: SAVED FILE INCONSISTENT '
      ENDIF
      DO I = 1 , N
	 READ(99) RX(I),VX(I),AX(I),BX(I),CX(I),DX(I)
      ENDDO
      
      DO I = 1 , N
         PSI0(I)=RX(I)
      ENDDO
      
      CLOSE(99)
      
      RETURN
      END