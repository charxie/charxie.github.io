      REAL*8 FUNCTION RANDOM(ISEED)
C-----------------------------------------------------------------------
C     RANDOM NUMBER GENERATOR: UNIFORM DISTRIBUTION (0,1)
C     ISEED: SEED FOR GENERATOR. ON THE FIRST CALL THIS HAS TO
C     HAVE A VALUE IN THE EXCLUSIVE RANGE (1, 2147483647)
C     AND WILL BE REPLACED BY A NEW VALUE TO BE USED IN
C     FOLLOWING CALL.
C
C     REF: Lewis, P.A.W., Goodman, A.S. & Miller, J.M. (1969)
C     "Pseudo-random number generator for the System/360", IBM
C     Systems Journal 8, 136.
C
C     This is a "high-quality" machine independent generator.
C     INTEGERS are supposed to be 32 bits or more.
C     The same algorithm is used as the basic IMSL generator.
C
C      Author: Lennart Nilsson
C
      INTEGER ISEED
      REAL*8 DSEED,DIVIS,DENOM,MULTIP
      DATA  DIVIS/2147483647.D0/
      DATA  DENOM /2147483711.D0/
      DATA  MULTIP/16807.D0/

      IF(ISEED.LE.1) ISEED=314159
      DSEED=MULTIP*ISEED
      DSEED=MOD(DSEED,DIVIS)
      RANDOM=DSEED/DENOM
      ISEED=DSEED

      RETURN
      END

	subroutine interpol(nstp,ndim,tmd,h,cof)
	implicit real*8(a-h,o-z)
	real*4 tmd(nstp),h(ndim,ndim,nstp),cof(ndim,ndim,nstp,3)
	real*4 hij(nstp),cofij(nstp,3)
        common /io/ iin,iout
	
	do i = 1 , ndim
	   do j = 1 , ndim
	      do k = 1 , nstp
	         hij(k)=h(i,j,k)
	      enddo
	      call icsccu(tmd,hij,nstp,cofij,nstp-1,ier)
	      do k = 1 , nstp
	         do l = 1 , 3
	            cof(i,j,k,l)=cofij(k,l)
	         enddo
	      enddo
	   enddo
	enddo

	do l = 1 , 3
	   write(iout,*)
	   write(iout,*) 'l=',l
	   cmax=0.0
	   imax=0
	   jmax=0
	   kmax=0
	   do k = 1 , nstp
	      do i = 1 , ndim
	         do j = 1 , ndim
	            if(abs(cof(i,j,k,l)).ge.cmax)then
	               cmax=abs(cof(i,j,k,l))
	               imax=i
	               jmax=j
	               kmax=k
	            endif
	         enddo
	      enddo
	   enddo
	   write(iout,*)imax,jmax,kmax,cmax
	enddo
	
	return
	end
	

	subroutine storesteps(nstp)
	integer nstp
	
	open(99,file='step',status='unknown',form='formatted')
	write(99,*) nstp
	close(99)
	
	return
	end
	
	subroutine getsteps(nstp)
	integer nstp

	open(99,file='step',status='old',form='formatted')
	read(99,*) nstp
	close(99)
	
	return
	end
	
c calculate the inverse of a matrix

	subroutine matinv(n,a,b)
	implicit real*8(a-h,o-z)
	dimension a(n,n),b(n,n),indx(n)
	
	do i = 1 , n
	   do j = 1 , n
	      b(i,j)=0.d0
	   enddo
	   b(i,i)=1.d0
	enddo
	
	call ludcmp(a,n,n,indx,d)

	do j = 1 , n
	   call lubksb(a,n,n,indx,b(1,j))
	enddo

	return
	end

      subroutine rnum(iseed)
      implicit real*8 (a-h,o-z)
      real*4 ran1
      common /ranum/ a(1000)
      
      do i = 1 , 1000
         a(i)=ran1(iseed)
      enddo    
      
      return
      end

      subroutine slope(x,y,d1y,n)
c
c--> calculates the slope of a curve
c
      implicit real*8 (a-h,o-z)
      dimension x(n),y(n)
      dimension d1y(n)
      do i=1,n-1
         if(x(i+1).ne.x(i)) then
            d1y(i)=(y(i+1)-y(i))/(x(i+1)-x(i))
         else
            d1y(i)=0.d0
         endif
      enddo
      return
      end

      subroutine deriv(x,y,d1y,d2y,n)
c
c--> calculates the first and second derivatives
c
      implicit real*8 (a-h,o-z)
      dimension x(n),y(n)
      dimension d1y(n),d2y(n)
      h=x(2)-x(1)
      d1y(1)=(-y(3)+4.d0*y(2)-3.d0*y(1))/(2.d0*h)
      d2y(1)=(y(3)-2.d0*y(2)+y(1))/(h**2.d0)
      d1y(n)=(3.d0*y(n)-4.d0*y(n-1)+y(n-2))/(2.d0*h)
      d2y(n)=(y(n)-2.d0*y(n-1)+y(n-2))/(h**2.d0)
      do i=2,n-1
         d1y(i)=(y(i+1)-y(i-1))/(2.d0*h)
         d2y(i)=(y(i+1)-2.d0*y(i)+y(i-1))/(h**2.0)
      enddo
      return
      end


      complex function sdot(n,x,y)

c  inner product of two complex vectors S=X*Y

      implicit real*8 (a-h,o-z)
      complex x(n),y(n),x_conj(n),tot
      tot = 0.d0
      incx=1
      incy=1
      xr=0.d0
      xi=0.d0

      do i = 1 , n
         xr=real(x(i))
         xi=aimag(x(i))
	 x_conj(i)=(1.0,0.0)*xr-(0.0,1.0)*xi
      enddo
      
      do i = 1 , n
         ix = 1 + (i-1)*incx
         iy = 1 + (i-1)*incy
         tot = tot + x_conj(ix)*y(iy)
      enddo   
      sdot = tot
      
      return
      end

c find coefficients for interpolating polynomial

      SUBROUTINE POLCOE(X,Y,N,COF)
      implicit real*8(a-h,o-z)
      PARAMETER (NMAX=50)
      DIMENSION X(N),Y(N),COF(N),S(NMAX)
      DO 11 I=1,N
        S(I)=0.
        COF(I)=0.
11    CONTINUE
      S(N)=-X(1)
      DO 13 I=2,N
        DO 12 J=N+1-I,N-1
          S(J)=S(J)-X(I)*S(J+1)
12      CONTINUE
        S(N)=S(N)-X(I)
13    CONTINUE
      DO 16 J=1,N
        PHI=N
        DO 14 K=N-1,1,-1
          PHI=K*S(K+1)+X(J)*PHI
14      CONTINUE
        FF=Y(J)/PHI
        B=1.
        DO 15 K=N,1,-1
          COF(K)=COF(K)+B*FF
          B=S(K)+X(J)*B
15      CONTINUE
16    CONTINUE
      RETURN
      END

c Langrange interpolation

      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
      implicit real*8(a-h,o-z)
      PARAMETER (NMAX=50) 
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N 
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.)PAUSE
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END


C THE LU DECOMPOSITION ROUTINES
C A=LU 
C D,INDX FOR INTERNAL USE

      SUBROUTINE LUDCMP(A,N,NP,INDX,D)

      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMAX=1500,TINY=1.0D-20)
      DIMENSION A(NP,NP),INDX(N),VV(NMAX)

      D=1.
      DO 12 I=1,N
        AAMAX=0.
        DO 11 J=1,N
          IF (DABS(A(I,J)).GT.AAMAX) AAMAX=DABS(A(I,J))
11      CONTINUE
        IF (AAMAX.EQ.0.) PAUSE 'Singular matrix.'
        VV(I)=1./AAMAX
12    CONTINUE
      DO 19 J=1,N
        DO 14 I=1,J-1
          SUM=A(I,J)
          DO 13 K=1,I-1
            SUM=SUM-A(I,K)*A(K,J)
13        CONTINUE
          A(I,J)=SUM
14      CONTINUE
        AAMAX=0.
        DO 16 I=J,N
          SUM=A(I,J)
          DO 15 K=1,J-1
            SUM=SUM-A(I,K)*A(K,J)
15        CONTINUE
          A(I,J)=SUM
          DUM=VV(I)*DABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
16      CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
17        CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(A(J,J).EQ.0.)A(J,J)=TINY
        IF(J.NE.N)THEN
          DUM=1./A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
18        CONTINUE
        ENDIF
19    CONTINUE

      RETURN
      END

C SOLVE THE SET OF N LINEAR EQUATIONS: AX=B;
C B IS INPUT AS THE RHS VECTOR, AND RETURNS THE SOLUTION X;
C A IS THE LU DECOMPOSITION OF THE MATRIX A, GIVEN BY LUDCMP;

      SUBROUTINE LUBKSB(A,N,NP,INDX,B)      
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NMAX=1500)
      DIMENSION A(NP,NP),INDX(NP),B(NP)

      II=0
      DO 12 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO 11 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
11        CONTINUE
        ELSE IF (SUM.NE.0.) THEN
          II=I
        ENDIF
        B(I)=SUM
12    CONTINUE
      DO 14 I=N,1,-1
        SUM=B(I)
        IF(I.LT.N)THEN
          DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
13        CONTINUE
        ENDIF
        B(I)=SUM/A(I,I)
14    CONTINUE

      RETURN
      END

c*********************************************************************
c Random number generator
c*********************************************************************

      FUNCTION RAN1(IDUM)
      DIMENSION R(97)
      PARAMETER (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6)
      PARAMETER (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6)
      PARAMETER (M3=243000,IA3=4561,IC3=51349)
      DATA IFF /0/
      IF (IDUM.LT.0.OR.IFF.EQ.0) THEN
        IFF=1
        IX1=MOD(IC1-IDUM,M1)
        IX1=MOD(IA1*IX1+IC1,M1)
        IX2=MOD(IX1,M2)
        IX1=MOD(IA1*IX1+IC1,M1)
        IX3=MOD(IX1,M3)
        DO 11 J=1,97
          IX1=MOD(IA1*IX1+IC1,M1)
          IX2=MOD(IA2*IX2+IC2,M2)
          R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
11      CONTINUE
        IDUM=1
      ENDIF
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IA2*IX2+IC2,M2)
      IX3=MOD(IA3*IX3+IC3,M3)
      J=1+(97*IX3)/M3
      IF(J.GT.97.OR.J.LT.1)PAUSE
      RAN1=R(J)
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
      RETURN
      END
      
c-----------------------------------------------------------------
c Interpolation routine
c
c   arguments    x      - vector of length nx containing the abscissae
c                           of the nx data points (x(i),y(i)) i=1,...,
c                           nx. (input) x must be ordered so that
c                           x(i) .lt. x(i+1).
c                y      - vector of length nx containing the ordinates
c                           (or function values) of the nx data points.
c                           (input)
c                nx     - number of elements in x and y. (input) nx
c                           must be .ge. 2.
c                c      - spline coefficients. (output) c is an nx-1 by
c                           3 matrix. the value of the spline
c                           approximation at t is
c                           s(t) = ((c(i,3)*d+c(i,2))*d+c(i,1))*d+y(i)
c                           where x(i) .le. t .lt. x(i+1) and
c                           d = t-x(i).
c                ic     - row dimension of matrix c exactly as
c                           specified in the dimension statement in
c                           the calling program. (input)
c                ier    - error parameter. (output)
c                         terminal error
c                           ier = 129, ic is less than nx-1.
c                           ier = 130, nx is less than 2.
c                           ier = 131, input abscissa are not ordered
c                             so that x(1) .lt. x(2) ... .lt. x(nx).
c
c-----------------------------------------------------------------
      subroutine icsccu (x,y,nx,c,ic,ier)
      implicit real*4(a-h,o-z)
      dimension           x(nx),y(nx),c(nx,3)
      dimension           cnx(3)
      common /io/ iin,iout

c                                  first executable statement
      nm1 = nx-1
      ier = 129
      if (ic .lt. nm1) go to 9000
      ier = 130
      if (nx .lt. 2) go to 9000
      ier = 131
      if (nx .eq. 2) go to 45
c                                  compute not-a-knot spline
      do 5 m = 2,nm1
         mm1=m-1
         c(m,2) = x(m)-x(mm1)
         if (c(m,2).le.0.0) go to 9000
         c(m,3) = (y(m)-y(mm1))/c(m,2)
    5 continue
      cnx(2) = x(nx)-x(nm1)
      if (cnx(2).le.0.0) go to 9000
      cnx(3) = (y(nx)-y(nm1))/cnx(2)
      ier = 0
      nm2 = nx-2
      if (nx .gt. 3) go to 10
      c(1,3) = cnx(2)
      c(1,2) = c(2,2)+cnx(2)
      c(1,1) = ((c(2,2)+2.*c(1,2))*c(2,3)*cnx(2)+c(2,2)**2*cnx(3))
     1/c(1,2)
      go to 20
   10 c(1,3) = c(3,2)
      c(1,2) = c(2,2)+c(3,2)
      c(1,1) = ((c(2,2)+2.*c(1,2))*c(2,3)*c(3,2)+c(2,2)**2*c(3,3))
     1/c(1,2)
      do 15 m=2,nm2
         mp1=m+1
         mm1=m-1
         g = -c(mp1,2)/c(mm1,3)
         c(m,1) = g*c(mm1,1)+3.*c(m,2)*c(mp1,3)+3.*c(mp1,2)*c(m,3)
         c(m,3) = g*c(mm1,2)+2.*c(m,2)+2.*c(mp1,2)
   15 continue
   20 g = -cnx(2)/c(nm2,3)
      c(nm1,1) = g*c(nm2,1)+3.*c(nm1,2)*cnx(3)+3.*cnx(2)*c(nm1,3)
      c(nm1,3) = g*c(nm2,2)+2.*c(nm1,2)+2.*cnx(2)
      if (nx.gt.3) go to 25
      cnx(1)=2.*cnx(3)
      cnx(3)=1.
      g=-1./c(nm1,3)
      go to 30
   25 g = c(nm1,2)+cnx(2)
      cnx(1) = ((cnx(2)+2.*g)*cnx(3)*c(nm1,2)+cnx(2)**2*
     1(y(nm1)-y(nx-2))/c(nm1,2))/g
      g = -g/c(nm1,3)
      cnx(3) = c(nm1,2)
   30 cnx(3) = g*c(nm1,2)+cnx(3)
      cnx(1) = (g*c(nm1,1)+cnx(1))/cnx(3)
      c(nm1,1) = (c(nm1,1)-c(nm1,2)*cnx(1))/c(nm1,3)
      do 35 jj=1,nm2
         j = nm1-jj
         c(j,1) = (c(j,1)-c(j,2)*c(j+1,1))/c(j,3)
   35 continue
      do 40 i=2,nm1
         im1 = i-1
         dtau = c(i,2)
         divdf1 = (y(i)-y(im1))/dtau
         divdf3 = c(im1,1)+c(i,1)-2.*divdf1
         c(im1,2) = (divdf1-c(im1,1)-divdf3)/dtau
         c(im1,3) = divdf3/dtau**2
   40 continue
      dtau = cnx(2)
      divdf1 = (y(nx)-y(nm1))/dtau
      divdf3 = c(nm1,1)+cnx(1)-2.*divdf1
      c(nm1,2) = (divdf1-c(nm1,1)-divdf3)/dtau
      c(nm1,3) = divdf3/dtau**2
      go to 9005
   45 if (x(1) .ge. x(2)) go to 9000
      ier = 0
      c(1,1) = (y(2)-y(1))/(x(2)-x(1))
      c(1,2) = 0.0
      c(1,3) = 0.0
      go to 9005
 9000 continue
      write(iout,*)  ' warning:   ier.ne. 0  in icsccu'
      write(iout,*)  '                  ier=',ier
      write(69,*) ' warning:   ier.ne. 0  in icsccu'
      write(69,*) '            ier=',ier
      stop
 9005 return
      end
