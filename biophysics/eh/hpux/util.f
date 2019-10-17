c--------------------------------------------------------------------                                                                               
c to store the number of steps already run
c--------------------------------------------------------------------
	subroutine storesteps(nstp,dotsa)
	integer nstp
	logical dotsa
	
	if(dotsa) then
	   open(99,file='step2s',status='unknown',form='formatted')
	else
	   open(99,file='stepns',status='unknown',form='formatted')
	endif

	write(99,*) nstp
	close(99)
	
	return
	end
	
c--------------------------------------------------------------------                                                                               
c to get the number of steps already run
c--------------------------------------------------------------------
	subroutine getsteps(nstp,dotsa)
	integer nstp
	logical dotsa

        if(dotsa) then
	   open(99,file='step2s',status='old',form='formatted')
	else
	   open(99,file='stepns',status='old',form='formatted')
	endif

	read(99,*) nstp
	close(99)
	
	return
	end
	
c--------------------------------------------------------------------                                                                               
c to convert the atom names into a proper format                                
c--------------------------------------------------------------------
      subroutine getnam(i,string,it)  
      implicit real*8(a-h,o-z)
      character keywrd*320,sline*80,string*80                                   
                                                                               
      common/iounit/iin,iout,istart(40),iungf,iunhl                                         
      common/keywrd/keywrd,sline                                                
                                                                               
      iia=ichar('A')                                                            
      iiz=ichar('Z')                                                            
      if(istart(i).eq.81) then                                                  
         it=1                                                                     
         return                                                                   
      end if                                                                    
      len=istart(I+1)-istart(I)                                                 
      if(len.gt.4) len=4                                                        
      string=sline(istart(I):istart(I+1)-1)                                     
      iis=ichar(string(2:2))                                                    
      if(iis.lt.iia.or.iis.gt.iiz) then                                         
         do k=len,1,-1                                                          
            string(K+1:K+1)=string(K:K)                                               
         enddo
         string(1:1)=' '                                                           
      end if                                                                    
      if(len.lt.4) then                                                         
         do K=len+1,4                                                           
            string(K:K)=' '                                                           
         enddo
      end if                                                                    

      return                                                                    
      end

c--------------------------------------------------------------------
c processes orbital name information                                            
c--------------------------------------------------------------------
      subroutine getorb(itd,isd,ipd,idd,ifd,ip,id,if,nn,text)        
      implicit real*8(a-h,o-z)                                        
      character text*10,keywrd*320,sline*80                                     
                                                                               
      common/iounit/iin,iout,istart(40),iungf,iunhl                                         
      common/keywrd/keywrd,sline                                                

      dimension ip(3),id(5),if(7),text(20)                                      

      itd=index(sline,'T.')                                                     
      isd=index(sline,'S.')                                                     
      ipd=index(sline,'P.')                                                     
      ip(1)=index(sline,'PX')                                                   
      ip(2)=index(sline,'PY')                                                   
      ip(3)=index(sline,'PZ')                                                   
      idd=index(sline,'D.')                                                     
      id(1)=index(sline,'DX2-Y2')                                               
      id(2)=index(sline,'DZ2')                                                  
      id(3)=index(sline,'DXY')                                                  
      id(4)=index(sline,'DXZ')                                                  
      id(5)=index(sline,'DYZ')                                                  
      ifd=index(sline,'F.')                                                     
      if(1)=index(sline,'FZ3')                                                  
      if(2)=index(sline,'FXZ2')                                                 
      if(3)=index(sline,'FYZ2')                                                 
      if(4)=index(sline,'FXYZ')                                                 
      if(5)=index(sline,'FZ(X2-Y2)')                                            
      if(6)=index(sline,'FX(X2-3Y2)')                                           
      if(7)=index(sline,'FY(3X2-Y2)')                                           
      nn=1                                                                      
      if(itd.ne.0) then                                                         
         nn=nn+1                                                                  
         text(nn)='t'                                                             
      end if                                                                    
      if(isd.ne.0) then                                                         
         nn=nn+1                                                                  
         text(nn)='s'                                                             
      end if                                                                    
      if(ipd.ne.0) then                                                         
         nn=nn+1                                                                  
         text(nn)='p'                                                             
      end if                                                                    
      if(ip(1).ne.0) then                                                       
         nn=nn+1                                                                  
         text(nn)='px'                                                            
      end if                                                                    
      if(ip(2).ne.0) then                                                       
         nn=nn+1                                                                  
         text(nn)='py'                                                            
      end if                                                                    
      if(ip(3).ne.0) then                                                       
         nn=nn+1                                                                  
         text(nn)='pz'                                                            
      end if                                                                    
      if(idd.ne.0) then                                                         
         nn=nn+1                                                                  
         text(nn)='d'                                                             
      end if                                                                    
      if(id(1).ne.0) then                                                       
         nn=nn+1                                                                  
         text(nn)='dx2-y2'                                                        
      end if                                                                    
      if(id(2).ne.0) then                                                       
         nn=nn+1                                                                  
         text(nn)='dz2'                                                           
      end if                                                                    
      if(id(3).ne.0) then                                                       
         nn=nn+1                                                                  
         text(nn)='dxy'                                                           
      end if                                                                    
      if(id(4).ne.0) then                                                       
         nn=nn+1                                                                  
         text(nn)='dxz'                                                           
      end if                                                                    
      if(id(5).ne.0) then                                                       
         nn=nn+1                                                                  
         text(nn)='dyz'                                                           
      end if                                                                    
      if(ifd.ne.0) then                                                         
         nn=nn+1                                                                  
         text(nn)='f'                                                             
      end if                                                                    
      if(if(1).ne.0) then                                                       
         nn=nn+1                                                                  
         text(nn)='fz3'                                                           
      end if                                                                    
      if(if(2).ne.0) then                                                       
         nn=nn+1                                                                  
         text(nn)='fxz2'                                                          
      end if                                                                    
      if(if(3).ne.0) then                                                       
         nn=nn+1                                                                  
         text(nn)='fyz2'                                                          
      end if                                                                    
      if(if(4).ne.0) then                                                       
         nn=nn+1                                                                  
         text(nn)='fxyz'                                                          
      end if                                                                    
      if(if(5).ne.0) then                                                       
         nn=nn+1                                                                  
         text(nn)='fz(x2-y2)'                                                     
      end if                                                                    
      if(if(6).ne.0) then                                                       
         nn=nn+1                                                                  
         text(nn)='fx(x2-3y2)'                                                    
      end if                                                                    
      if(if(7).ne.0) then                                                       
         nn=nn+1                                                                  
         text(nn)='fy(3x2-y2)'                                                    
      end if                                                                    
      nn=nn+1                                                                   

      return                                                                    
      end      
c--------------------------------------------------------------------
c reads unformatted input lines                                                 
c--------------------------------------------------------------------
      subroutine readl
      implicit double precision (A-H,O-Z)                                       
      logical leadsp                                                            
      character keywrd*320,sline*80,space*1                                     
     *,nine*1,zero*1,tab*1, comma*1, string*80                                  
      data comma,space,nine,zero/',',' ','9','0'/                               

      common /iounit/ iin,iout,istart(40),iungf,iunhl                                      
      common/keywrd/keywrd,sline                                                

      ilowa = ichar('a')                                                        
      ilowz = ichar('z')                                                        
      icapa = ichar('A')                                                        

   10 read(iin,'(a)',end=60)sline                                               
      if(sline.eq.' ') goto 10
                                                        
c---> a line headed by an exclamation mark is a comment in the input file

      if(index(sline,'!').ne.0) goto 10

c---> turn lower-case letters into upper-case ones
c     so that the letters from the input cards will be capital letters
c     in the output file

      do 20 i=1,80                                                              
         iline=ichar(sline(i:i))                                                
         if(iline.ge.ilowa.and.iline.le.ilowz) then                             
            sline(i:i)=char(iline+icapa-ilowa)                                  
         endif                                                                  
   20 continue                                                                  

c---> put meaningless letters as blanks

      do 30 i=1,80                                                              
   30 if(sline(i:i).lt.space .or.                                               
     * sline(i:i).eq.comma .or.                                                 
     * sline(i:i).gt.char(125)) sline(i:i)=space                                 

c---> initialize istart to interpret blanks as zero's                             

      do 40 i=1,40                                                              
   40 istart(i)=81                                                              

c---> find initial digit of all numbers, check for leading spaces 
c     followed by a character and store in istart                                        

      leadsp=.true.                                                             
      nvalue=0                                                                  
      do 50 i=1,80                                                              
         if (leadsp.and.sline(i:i).ne.space) then                               
            nvalue=nvalue+1                                                     
            istart(nvalue)=i                                                    
         end if                                                                 
         leadsp=(sline(i:i).eq.space)                                           
   50 continue                                                                  

      return                                                                    

   60 sline='END'                                                               
      return                                                                    
      end                                                                       

c--------------------------------------------------------------------
c writes error messages                                                         
c--------------------------------------------------------------------
      subroutine writel(l)
      implicit double precision (A-H,O-Z)                                       
      character*1 star
      character l*80,keywrd*320,SLINE*80                                        

      common /iounit/ iin,iout,istart(40),iungf,iunhl                                     
      common/keywrd/keywrd,sline                                                
      data star/'*'/

      do 10 i=80,1,-1                                                           
      if(l(i:i).ne.' ') then                                                    
       il=i                                                                     
       goto 11                                                                  
      end if                                                                    
   10 continue                                                                  
   11 continue                                                                  

      write(iout,1000)(star,i=1,il)
 1000 format(/1x,80a1)
      write(iout,1001)
 1001 format('  ')
      write(iout,1010)l                                                         
 1010 format(1x,a)                                                              
      write(iout,1000)(star,i=1,il)
      write(iout,1001)
      if(index(keywrd,'FORCE').eq.0) then                                       
       stop                                                                     
      else                                                                      
       return                                                                   
      end if                                                                    

      end

c--------------------------------------------------------------------
c translates ASCII characters into numbers    
c--------------------------------------------------------------------                                                                               
      double precision function reada(a,istart,len)                             
      implicit double precision (a-h,o-z)                                       
      character*1 a(320)                                                        

      nine=ichar('9')                                                           
      izero=ichar('0')                                                          
      minus=ichar('-')                                                          
      idot=ichar('.')                                                           
      ie=ichar('E')                                                             
      idig=0                                                                    
      c1=0                                                                      
      c2=0                                                                      
      ic3=0                                                                     
      one=1.d0                                                                  
      two=1.d0                                                                  
      x = 1.d0                                                                  

      do 10 j=istart,len                                                        
         n=ichar(a(j))                                                          
         m=ichar(a(j+1))                                                        
         if(n.le.nine.and.n.ge.izero .or.n.eq.idot)goto 20                      
         if(n.eq.minus.and.(m.le.nine.and.m.ge.izero                            
     * .or. m.eq.idot)) goto 20                                                 
   10 continue                                                                  
      reada=0.d0                                                                
      return                                                                    

   20 continue                                                                  
      do 30 i=j,len                                                             
         n=ichar(a(i))                                                          
         if(n.le.nine.and.n.ge.izero) then                                      
            idig=idig+1                                                         
            if (idig.gt.10) goto 60                                             
            c1=c1*10+n-izero                                                    
         else if(n.eq.minus.and.i.eq.j) then                                    
            one=-1.d0                                                           
         else if(n.eq.idot) then                                                
            goto 40                                                             
         else                                                                   
            if(n.eq.ie) goto 60                                                 
            goto 80                                                             
         endif                                                                  
   30 continue                                                                  
   40 continue                                                                  
      idig=0                                                                    
      do 50 ii=i+1,len                                                          
         n=ichar(a(ii))                                                         
         if(n.le.nine.and.n.ge.izero) then                                      
            idig=idig+1                                                         
            if (idig.gt.10) goto 60                                             
            c2=c2*10+n-izero                                                    
            x = x /10                                                           
         elseif(n.eq.minus.and.ii.eq.i) then                                    
            x=-x                                                                
         else if(n.eq.ie) then                                                  
            goto 60                                                             
         else                                                                   
            goto 80                                                             
         endif                                                                  
   50 continue                                                                  

   60 idig=0                                                                    
      do 70 k=ii+1,len                                                          
         n=ichar(a(k))                                                          
         if(n.le.nine.and.n.ge.izero) then                                      
            idig=idig+1                                                         
            if (idig.gt.10) goto 80                                             
            ic3=ic3*10+n-izero                                                  
         else if(n.eq.minus.and.k.eq.ii+1) then                                 
            two=-1.D0                                                           
         else                                                                   
            goto 80                                                             
         endif                                                                  
   70 continue                                                                  
                                                                                
c put the pieces together                                                       

   80 continue                                                                  
      coeff=1.0                                                                 
      if(two.eq.-1) then                                                        
      do 90 i=1,ic3                                                             
      coeff=coeff/10.0d0                                                        
   90 continue                                                                  
      else                                                                      
      do 100 i=1,ic3                                                            
      coeff=coeff*10.0d0                                                        
  100 continue                                                                  
      end if                                                                    
      reada=one*(c1+c2*x)*coeff                                                 
      return                                                                    
      end      
      
                                                                                    
c--------------------------------------------------------------------
c vector product
c--------------------------------------------------------------------                                                                               
      subroutine vprod(vp,x,y)                                                  
      implicit double precision(a-h,o-z)                                        
                                                                               
      dimension vp(3),x(3),y(3)                                                 
                                                                               
      vp(1)=x(2)*y(3)-x(3)*y(2)                                                 
      vp(2)=x(3)*y(1)-x(1)*y(3)                                                 
      vp(3)=x(1)*y(2)-x(2)*y(1)                                                 
      return                                                                    

      end        

c--------------------------------------------------------------------
c   the following four functions provide timing information
c--------------------------------------------------------------------
c  timing routines
c
      subroutine initsec()
      implicit real(a-h,o-z)
      real*4 tarray(2),start,etime
      common /seccom/ start
c      start = etime(tarray)
      return
      end
      function seconds()
      implicit real(a-h,o-z)
      real*4 tarray(2),start,etime
      common /seccom/ start
c      t = etime(tarray) - start
c      seconds = t
      return
      end
c
c  return the remaining time left for this job
c  not meaning for this version, just return a large number
c
      subroutine trmain(timlft)
      timlft = 1.e4
      return
      end
c
      subroutine mytime(timestr)
      character*24 timestr,ctime
      integer*4 time
      timestr = ctime(time())
      return
      end
      
c----------------------------------------------------------------------
c dot product of two integer vectors
c----------------------------------------------------------------------      
      subroutine ivdot(n,a,b,c)
      integer n,a(n),b(n),c
      
      do k = 1 , n
         c=c+a(k)*b(k)
      enddo
      
      return  
      end

c----------------------------------------------------------------------
c matrix transformation (symmetric)
c----------------------------------------------------------------------      
      subroutine mtran(n,t,a,b)
      implicit real*8 (a-h,o-z)
      dimension a(n,n),b(n,n),t(n,n)
      
      do i = 1 , n
         do j = 1 , i
            b(i,j)=0.d0
            do k = 1 , n
               do l = 1 , n
                  b(i,j)=b(i,j)+t(i,k)*a(k,l)*t(l,j)
               enddo
            enddo
         enddo
         b(j,i)=b(i,j)
      enddo
      
      return  
      end

c----------------------------------------------------------------------
c  random number generator
c----------------------------------------------------------------------
      subroutine lranst(iseed)
      implicit real (a-h,o-z)
      real*4 rand
      if (iseed.le.0) then
         x = rand(time())
      else
         x = rand(iseed)
      endif
      return
      end
      real function ranl()
      implicit real (a-h,o-z)
      real*4 rand
      ranl = rand(0)
      return
      end

c----------------------------------------------------------------------
C THE LU DECOMPOSITION ROUTINES
C A=LU 
C D,INDX FOR INTERNAL USE
c----------------------------------------------------------------------
      SUBROUTINE LUDCMP(A,N,NP,INDX,D)

      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMAX=500,TINY=1.0D-20)
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

c----------------------------------------------------------------------
C SOLVE THE SET OF N LINEAR EQUATIONS: AX=B;
C B IS INPUT AS THE RHS VECTOR, AND RETURNS THE SOLUTION X;
C A IS THE LU DECOMPOSITION OF THE MATRIX A, GIVEN BY LUDCMP;
c----------------------------------------------------------------------
      SUBROUTINE LUBKSB(A,N,NP,INDX,B)      
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NMAX=500)
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
	
c--------------------------------------------------------------------
c calculate the first and second derivatives
c--------------------------------------------------------------------
      subroutine deriv(x,y,d1y,d2y,n)
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
      
c--------------------------------------------------------------------
c calculate the first derivative
c--------------------------------------------------------------------
      subroutine deriv1(x,y,d1y,n)
      implicit real*8 (a-h,o-z)
      dimension x(n),y(n)
      dimension d1y(n)
      h=x(2)-x(1)
      d1y(1)=(-y(3)+4.d0*y(2)-3.d0*y(1))/(2.d0*h)
      d1y(n)=(3.d0*y(n)-4.d0*y(n-1)+y(n-2))/(2.d0*h)
      do i=2,n-1
         d1y(i)=(y(i+1)-y(i-1))/(2.d0*h)
      enddo
      return
      end      