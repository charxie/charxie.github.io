        program phofcc 
c 
c calculate phonon dispersion using the lattice-inverted embedded-atom 
c method potential 
c for fcc structure 
c 
c written by: 
c                           Qian Xie 
c 
c 
c 
      implicit real*8 (a-h,o-z) 
      real*8 lambda 
      parameter (npt=500) 
      common /debug/ipot,idyn,ifcc 
      common /xyzr/x(200),y(200),z(200),r(200) 
      common /poten/rpt(npt),den(npt),dend1(npt),dend2(npt), 
     +                       pap(npt),papd1(npt),papd2(npt) 
      common /intpl/cd1(npt,3),cd2(npt,3),cp1(npt,3),cp2(npt,3) 
      dimension qx(100),qy(100),qz(100) 
      dimension dyna(3,3),trid(3,3),tman(3,3) 
      dimension b(3),c(3),w(100,3) 
      character cell*3, filen*10, wavec*3 
c 
      call pmina(cell, key,'phofcc.in ' , 'element   ') 
      filen(1:6)='poten.' 
      filen(7:9)=cell 
      open(20,file=filen,status='unknown') 
      call pmina(cell, key,'phofcc.in ' , 'strk type ') 
      open(30,file=cell,status='unknown') 
      read(30,1500) (x(i),y(i),z(i),r(i),i=1,200) 
c 
      open(40,file='matrix',status='unknown') 
      open(69,file='recd',status='unknown') 
c 
      call putin(aeq,amass,wavec)
      call rdpot(lambda,aaa,rhoe)

      call pmina(cell, key,'phofcc.in ' , 'element   ')
      filen(1:2)=cell
      filen(3:6)='egn.'

      call pmina(wavec,  key,'phofcc.in ' , 'wave vecto')

      filen(7:9)=wavec
      open(10,file=filen,status='unknown')

c-->  print headings
      write(10,*) '\    k      w(k,1)      w(k,2)      w(k,3)' 
c
      if(ifcc.eq.1) then 
        write(6 ,1501) (i,x(i),y(i),z(i),r(i),i=1,200) 
        write(69,1501) (i,x(i),y(i),z(i),r(i),i=1,200) 
      endif 
c 
      twopi=8.d0*datan(1.d0) 
      eps=1.0d-10 
      omega=aeq*aeq*aeq/4.d0 
      amass=amass/6.02d0 
      rnne=aeq/dsqrt(2.d0) 
      frhod1=-aaa*rhoe**(1.d0/lambda-1.d0)/lambda 
      frhod2=frhod1*(1.d0/lambda-1.d0)/rhoe 
      do 100 i=1,200 
        x(i)=x(i)*aeq 
        y(i)=y(i)*aeq 
        z(i)=z(i)*aeq 
        r(i)=r(i)*aeq 
100   continue 
c 
      do 1000 k=1,100 
c 
        do 20 m=1,3 
        do 30 n=1,3 
        dyna(m,n)=0.d0 
 30     continue 
 20     continue 
c 
c--> wave vectors for four different directions 
c 
      if(wavec.eq.'00x') then 
        qx(k)=twopi/100.d0/aeq*dfloat(k)*0.0 
        qy(k)=twopi/100.d0/aeq*dfloat(k)*0.0 
        qz(k)=twopi/100.d0/aeq*dfloat(k)*1.0 
      endif 
      if(wavec.eq.'0x1') then 
        qx(k)=twopi/100.d0/aeq*dfloat(k)*0.0 
        qy(k)=twopi/100.d0/aeq*dfloat(k)*1.0 
        qz(k)=twopi/aeq 
      endif 
      if(wavec.eq.'0xx') then 
        qx(k)=twopi/100.d0/aeq*dfloat(k)*0.0 
        qy(k)=twopi/100.d0/aeq*dfloat(k)*1.0 
        qz(k)=twopi/100.d0/aeq*dfloat(k)*1.0 
      endif 
      if(wavec.eq.'xxx') then 
        qx(k)=twopi/100.d0/aeq*dfloat(k)*1.0 
        qy(k)=twopi/100.d0/aeq*dfloat(k)*1.0 
        qz(k)=twopi/100.d0/aeq*dfloat(k)*1.0 
      endif 
c 
      do 200 i=1,200 
      qr=qx(k)*x(i)+qy(k)*y(i)+qz(k)*z(i) 
      if(r(i).lt.rpt(1)) stop ' distance out of potential range' 
      if(r(i).lt.rpt(200)) then 
        itt=(r(i)-rpt(1))/(rpt(2)-rpt(1))+1 
        d=r(i)-rpt(itt) 
        dd1=((cd1(itt,3)*d+cd1(itt,2))*d+cd1(itt,1))*d+dend1(itt) 
        dd2=((cd2(itt,3)*d+cd2(itt,2))*d+cd2(itt,1))*d+dend2(itt) 
        pd1=((cp1(itt,3)*d+cp1(itt,2))*d+cp1(itt,1))*d+papd1(itt) 
        pd2=((cp2(itt,3)*d+cp2(itt,2))*d+cp2(itt,1))*d+papd2(itt) 
      else 
        pd1=0.d0 
        pd2=0.d0 
        dd1=0.d0 
        dd2=0.d0 
      endif 
c 
      d2=pd2+2.d0*frhod1*dd2 
      d1=pd1+2.d0*frhod1*dd1 
c 
        dyna(1,1)=dyna(1,1)+((d2-d1/r(i)) 
     +  *x(i)*x(i)/r(i)/r(i)+d1/r(i))*(1.d0-dcos(qr)) 
        dyna(1,2)=dyna(1,2)+((d2-d1/r(i)) 
     +  *x(i)*y(i)/r(i)/r(i))*(1.d0-dcos(qr)) 
        dyna(1,3)=dyna(1,3)+((d2-d1/r(i)) 
     &  *x(i)*z(i)/r(i)/r(i))*(1.d0-dcos(qr)) 
        dyna(2,2)=dyna(2,2)+((d2-d1/r(i)) 
     &  *y(i)*y(i)/r(i)/r(i)+d1/r(i))*(1.d0-dcos(qr)) 
        dyna(2,3)=dyna(2,3)+((d2-d1/r(i)) 
     &  *y(i)*z(i)/r(i)/r(i))*(1.d0-dcos(qr)) 
        dyna(3,3)=dyna(3,3)+((d2-d1/r(i)) 
     &  *z(i)*z(i)/r(i)/r(i)+d1/r(i))*(1.d0-dcos(qr)) 
200     continue 
c 
        dyna(2,1)=dyna(1,2) 
        dyna(3,1)=dyna(1,3) 
        dyna(3,2)=dyna(2,3) 
c 
        call many(tman,qx(k),qy(k),qz(k)) 
c 
        do 51 m=1,3 
        do 52 n=1,3 
        tman(m,n)=tman(m,n)*frhod2 
        dyna(m,n)=dyna(m,n)+tman(m,n) 
 52     continue 
 51     continue 
c 
c--> write dynamic matrix 
c 
      if(idyn.eq.1) then 
        if(k.eq.1) write(40,*) ' output of dynamic matrix ' 
        write(40,1041) 
        write(40,1044) k,qx(k),qy(k),qz(k) 
        write(40,1042) 
        write(40,1045)((dyna(i,j),j=1,3),i=1,3) 
        write(40,1043) 
        write(40,1045)((tman(i,j),j=1,3),i=1,3) 
      endif 
c 
c--> diagonalization 
c 
        call strq(dyna,3,trid,b,c) 
        call stq(3,b,c,trid,eps,l) 
        call bsort(b,3,1,3) 
c 
c w(q)-THz 
      do 8000 i=1,3 
8000  w(k,i)=dsqrt(b(i)*16.d0/amass)*10.d0/twopi 
      write(6 ,1600) k,w(k,1),w(k,2),w(k,3) 
      write(10,1600) k,w(k,1),w(k,2),w(k,3) 
c 
1000  continue 
c 
1041  format(2x,'q wave vectors') 
1042  format(2x,'dynamic matrix') 
1043  format(2x,'many-body part') 
1044  format(2x,i5,3f15.6) 
1045  format(2x,5x,3f15.6) 
1500  format(4f16.7) 
1501  format(i5,4f16.7) 
1600  format(2x,i5,3f20.15) 
c 
        stop 
        end 
c 
      subroutine many(t,qx,qy,qz) 
c 
c--> contribution of many-body potential 
c 
      implicit real*8 (a-h,o-z) 
      parameter(npt=500) 
      dimension t(3,3) 
      common /poten/rpt(npt),den(npt),dend1(npt),dend2(npt), 
     +                       pap(npt),papd1(npt),papd2(npt) 
      common /intpl/cd1(npt,3),cd2(npt,3),cp1(npt,3),cp2(npt,3) 
      common /xyzr/x(200),y(200),z(200),r(200) 
      tx=0.d0 
      ty=0.d0 
      tz=0.d0 
c 
      do 10 i=1,200 
      qr=qx*x(i)+qy*y(i)+qz*z(i) 
      if(r(i).lt.rpt(1)) stop ' distance out of potential range ' 
      if(r(i).lt.rpt(200)) then 
        itt=(r(i)-rpt(1))/(rpt(2)-rpt(1))+1 
        d=r(i)-rpt(itt) 
        dd1=((cd1(itt,3)*d+cd1(itt,2))*d+cd1(itt,1))*d+dend1(itt) 
      else 
        dd1=0.d0 
      endif 
      tx=tx+dsin(qr)*dd1*x(i)/r(i) 
      ty=ty+dsin(qr)*dd1*y(i)/r(i) 
      tz=tz+dsin(qr)*dd1*z(i)/r(i) 
10    continue 
c 
      t(1,1)=tx*tx 
      t(2,2)=ty*ty 
      t(3,3)=tz*tz 
      t(1,2)=tx*ty 
      t(2,3)=ty*tz 
      t(3,1)=tz*tx 
      t(2,1)=t(1,2) 
      t(3,2)=t(2,3) 
      t(1,3)=t(3,1) 
c 
      return 
      end 
c 
c 
      subroutine rdpot(lambda,aaa,rhoe) 
c 
c-->  input the interatomic potentials 
c 
      implicit real*8 (a-h,o-z) 
      real*8 lambda 
      parameter (npt=500) 
      common /debug/ipot,idyn,ifcc 
      common /poten/rpt(npt),den(npt),dend1(npt),dend2(npt), 
     +                       pap(npt),papd1(npt),papd2(npt) 
      common /intpl/cd1(npt,3),cd2(npt,3),cp1(npt,3),cp2(npt,3) 
      character dummy*1 
c 
      read(20,'(a1)') dummy 
      read(20,'(a1)') dummy 
      read(20,'(a1)') dummy 
      read(20,'(a1)') dummy 
      read(20,'(a1)') dummy 
      read(20,'(a1)') dummy 
      read(20,'(a1)') dummy 
      read(20,'(a1)') dummy 
      read(20,'(a1)') dummy 
      read(20,'(a1)') dummy 
      read(20,1001) dummy,lambda,aaa,rhoe 
      read(20,'(a1)') dummy 
      read(20,'(a1)') dummy 
      read(20,'(a1)') dummy 
      read(20,'(a1)') dummy 

      do 300 ip=1,npt 
300   read(20,1000) rpt(ip),den(ip),pap(ip) 
c 
      call deriv(rpt,den,dend1,dend2,npt) 
      call deriv(rpt,pap,papd1,papd2,npt) 
      call icsccu(rpt,dend1,npt,cd1,npt-1,ier) 
      call icsccu(rpt,dend2,npt,cd2,npt-1,ier) 
      call icsccu(rpt,papd1,npt,cp1,npt-1,ier) 
      call icsccu(rpt,papd2,npt,cp2,npt-1,ier) 
c 
      if(ipot.eq.1) then 
        write(69,*) 
        write(6 ,*) '  check the interatomic potentials ' 
        write(69,*) '  check the interatomic potentials ' 
        write(6 ,*) lambda,aaa,rhoe 
        write(69,*) lambda,aaa,rhoe 
        write(69,*) 
        do 15 ip=1,npt 
        write(6 ,1002) ip,rpt(ip),den(ip),pap(ip) 
        write(69,1002) ip,rpt(ip),den(ip),pap(ip) 
15      continue 
      endif 
c 
1000  format(3x,3f25.20) 
1001  format(3x,a1,3f15.8) 
1002  format(3x,i5,3f20.15) 
c 
      return 
      end 
c 
      subroutine putin(aeq,amass,wavec) 
c----------------------------------------------------------------- 
c  read in the debug parameters and inputs 
c----------------------------------------------------------------- 
c 
      implicit real*8(a-h,o-z) 
      common /debug/ipot,idyn,ifcc 
      character wavec*3 
c 
      call pmina(wavec,  key,'phofcc.in ' , 'wave vecto') 
      call pminr(aeq,    key,'phofcc.in ' , 'lattice co') 
      call pminr(amass,  key,'phofcc.in ' , 'atom mass ') 
c 
c------------debug parameters 
      call pmini(ipot, key,'phofcc.in ' , 'potential ') 
      call pmini(idyn, key,'phofcc.in ' , 'dynamic ma') 
      call pmini(ifcc, key,'phofcc.in ' , 'fcc struct') 
c 
      return 
      end 
c 
c 
      subroutine icsccu (x,y,nx,c,ic,ier) 
c----------------------------------------------------------------- 
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
c 
        implicit real*8(a-h,o-z) 
        dimension           x(500),y(500),c(500,3) 
        dimension           cnx(3) 
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
      write(6,*)  ' warning:   ier.ne. 0  in icsccu' 
      write(6,*)  '            ier=',ier 
      write(69,*) ' warning:   ier.ne. 0  in icsccu' 
      write(69,*) '            ier=',ier 
      stop 
 9005 return 
      end 
c 
        subroutine deriv(x,y,d1y,d2y,n) 
c 
c-->  calculates the first and second derivatives 
c 
        implicit real*8 (a-h,o-z) 
        dimension x(n),y(n) 
        dimension d1y(n),d2y(n) 
        h=x(2)-x(1) 
        d1y(1)=(-y(3)+4.d0*y(2)-3.d0*y(1))/(2.d0*h) 
        d2y(1)=(y(3)-2.d0*y(2)+y(1))/(h**2.d0) 
        d1y(n)=(3.d0*y(n)-4.d0*y(n-1)+y(n-2))/(2.d0*h) 
        d2y(n)=(y(n)-2.d0*y(n-1)+y(n-2))/(h**2.d0) 
        do 1 i=2,n-1 
        d1y(i)=(y(i+1)-y(i-1))/(2.d0*h) 
        d2y(i)=(y(i+1)-2.d0*y(i)+y(i-1))/(h**2.0) 
 1      continue 
        return 
        end 
c 
        subroutine strq(a,n,q,b,c) 
c-------------------------------------------------------------- 
c eigenvalue evaluation for a symmetric matrix 
c subroutine strq : Householder transformation 
c subroutine stq  : diagonalization 
c-------------------------------------------------------------- 
c 
        implicit real*8 (a-h,o-z) 
        dimension a(n,n),q(n,n),b(n),c(n) 
        do 10 i=1,n 
        do 10 j=1,n 
 10     q(i,j)=a(i,j) 
        do 80 i=n,2,-1 
        h=0.0 
        if(i.gt.2) then 
        do 20 k=1,i-1 
 20     h=h+q(i,k)*q(i,k) 
        endif 
        if(h+1.0.eq.1.0) then 
        c(i-1)=0.0 
        if(i.eq.2) c(i-1)=q(i,i-1) 
        b(i)=0.0 
        else 
        c(i-1)=dsqrt(h) 
        if(q(i,i-1).gt.0.0) c(i-1)=-c(i-1) 
        h=h-q(i,i-1)*c(i-1) 
        q(i,i-1)=q(i,i-1)-c(i-1) 
        f=0.0 
        do 50 j=1,i-1 
        q(j,i)=q(i,j)/h 
        g=0.0 
        do 30 k=1,j 
 30     g=g+q(j,k)*q(i,k) 
        if(j+1.le.i-1) then 
        do 40 k=j+1,i-1 
 40     g=g+q(k,j)*q(i,k) 
        endif 
        c(j-1)=g/h 
        f=f+g*q(j,i) 
 50     continue 
        h2=f/(h+h) 
        do 70 j=1,i-1 
        f=q(i,j) 
        g=c(j-1)-h2*f 
        c(j-1)=g 
        do 60 k=1,j 
 60     q(j,k)=q(j,k)-f*c(k-1)-g*q(i,k) 
 70     continue 
        b(i)=h 
        endif 
 80     continue 
        b(1)=0.0 
        do 130 i=1,n 
        if((b(i).ne.0.0).and.(i-1.ge.1)) then 
        do 110 j=1,i-1 
        g=0.0 
        do 90 k=1,i-1 
 90     g=g+q(i,k)*q(k,j) 
        do 100 k=1,i-1 
 100    q(k,j)=q(k,j)-g*q(k,i) 
 110    continue 
        endif 
        b(i)=q(i,i) 
        q(i,i)=1.0 
        if(i-1.ge.1) then 
        do 120 j=1,i-1 
        q(i,j)=0.0 
        q(j,i)=0.0 
 120    continue 
        endif 
 130    continue 
        return 
        end 
c 
        subroutine stq(n,b,c,q,eps,l) 
        implicit real*8 (a-h,o-z) 
        dimension b(n),c(n),q(n,n) 
        c(n)=0.0 
        d=0.0 
        f=0.0 
        do 50 j=1,n 
        it=0 
        h=eps*(dabs(b(j))+dabs(c(j))) 
        if(h.gt.d) d=h 
        m=j-1 
 10     m=m+1 
        if(m.le.n) then 
        if(dabs(c(m)).gt.d) goto 10 
        endif 
        if(m.ne.j) then 
 15     if(it.eq.60) then 
        l=0 
        write(*,18) 
 18     format(1x,'fail!') 
        return 
        endif 
        it=it+1 
        g=b(j) 
        p=(b(j+1)-g)/(2.0*c(j)) 
        r=dsqrt(p*p+1.0) 
        if(p.ge.0.0) then 
        b(j)=c(j)/(p+r) 
        else 
        b(j)=c(j)/(p-r) 
        endif 
        h=g-b(j) 
        do 20 i=j+1,n 
 20     b(i)=b(i)-h 
        f=f+h 
        p=b(m) 
        e=1.0 
        s=0.0 
        do 40 i=m-1,j,-1 
        g=e*c(i) 
        h=e*p 
        if(dabs(p).ge.dabs(c(i))) then 
        e=c(i)/p 
        r=dsqrt(e*e+1.0) 
        c(i+1)=s*p*r 
        s=e/r 
        e=1.0/r 
        else 
        e=p/c(i) 
        r=dsqrt(e*e+1.0) 
        c(i+1)=s*c(i)*r 
        s=1.0/r 
        e=e/r 
        endif 
        p=e*b(i)-s*g 
        b(i+1)=h+s*(e*g+s*b(i)) 
        do 30 k=1,n 
        h=q(k,i+1) 
        q(k,i+1)=s*q(k,i)+e*h 
        q(k,i)=e*q(k,i)-s*h 
 30     continue 
 40     continue 
        c(j)=s*p 
        b(j)=e*p 
        if(dabs(c(j)).gt.d) goto 15 
        endif 
        b(j)=b(j)+f 
 50     continue 
        do 80 i=1,n 
        k=i 
        p=b(i) 
        if(i+1.le.n) then 
        j=1 
 60     j=j+1 
        if(j.le.n) then 
        if(b(j).le.p) then 
        k=j 
        p=b(j) 
        goto 60 
        endif 
        endif 
        endif 
        if(k.ne.i) then 
        b(k)=b(i) 
        b(i)=p 
        do 70 j=1,n 
        p=q(j,i) 
        q(j,i)=q(j,k) 
        q(j,k)=p 
 70     continue 
        endif 
 80     continue 
        l=1 
        return 
        end 
c 
        subroutine bsort(a,n,mm,nn) 
c 
c--> order the eigenfrequencies 
c 
        implicit real*8 (a-h,o-z) 
        dimension a(n) 
        m=nn-mm+1 
 10     if(m.gt.0) then 
        j=m+mm-2 
        m=0 
        do 20 i=mm,j 
        if(a(i).gt.a(i+1)) then 
        d=a(i) 
        a(i)=a(i+1) 
        a(i+1)=d 
        m=i-mm+1 
        endif 
 20     continue 
        goto 10 
        endif 
        return 
        end 
c 

