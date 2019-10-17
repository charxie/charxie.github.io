      program phonon 
c----------------------------------------------------------------- 
c 
c    lattice dynamics  for complex system 
c      control pmtr:  phonon.in 
c      structure   :  cell.*** 
c      pair  ppt   :  see sub opfile 
c 
c----------------------------------------------------------------- 
      implicit real*8(a-h,o-z) 
      parameter (neigd=60,nkptd=3000) 
      real*8 mass(20) 
      common /struk/g(3,3),omega,natm,pos(3,20),ityp(20),ntyp(20), 
     +              ntype,neq(20),vint,area,a(3) 
      dimension qx(nkptd),qy(nkptd),qz(nkptd),wk(nkptd),wvk(3) 
      common /eigv/omcm(neigd),zr(neigd,neigd),zi(neigd,neigd),ne,nblw 
      character*2 elem(20) 
c 
      open(69,file='recd') 
      open( 8,file='band') 
c 
      call rdpmtr 
      call setcel(elem,mass) 
      call radist 
      call getppt(elem) 
      call iwtset(qx,qy,qz,wk,nk) 
c 
      write(8,110)  nk,natm 
      do 30 ik=1,nk 
      call dynmat(mass,qx(ik),qy(ik),qz(ik),wvk) 
      call eigen 
      write(69 ,120)  ik,qx(ik),qy(ik),qz(ik),wk(ik),ne 
      write( 8 ,120)  ik,qx(ik),qy(ik),qz(ik),wk(ik),ne 
      do 20 i=1,ne 
 
      write(16,'(4f15.7)') qx(ik),qy(ik),qz(ik),omcm(i) 
 
      write(6 ,130)  i, omcm(i) 
      write(8 ,130)  i, omcm(i) 
20    write(69,130)  i, omcm(i) 
30    continue 
c 
110   format('  k-point, nk=',i5, '     natom=',i5) 
120   format ('ik,kx,ky,kz,wk,ne = ',i5,4f12.7,i5) 
130   format(5x,i5,6f15.7) 
c 
      stop 
      end 
c 
      subroutine  dynmat(mass,w1,w2,w3,wvk) 
c 
c     calculate the dynamical matrix 
c 
c     input: 
c     mass .... mass of cation and anion in amu. 
c     w1,w2,w3 ...... wavevector in 2*pi/(a(1),a(2),a(3)) 
c 
c     output: 
c     omcm ..... eigenfrequencies in THZ 
c     zr,zi .... eigenvectors for omcm 
c 
      implicit real*8(a-h,o-z) 
      parameter (neigd=60,nkptd=3000) 
      real*8 mass(20) 
      dimension wvk(3) 
      common /debug / ipmt,idat,icel,idyn,ifce,ifun,ieig,ippt, 
     &                iwts,irad,irdp 
      common /struk/g(3,3),omega,natm,pos(3,20),ityp(20),ntyp(20), 
     +              ntype,neq(20),vint,area,a(3) 
      common /finst1/rg(3,3) 
      common /rdist/nsh,rshel(500) 
c--->   task commons for eigenvalue problem 
      common /secl/ ss(neigd,neigd),dd(neigd,neigd),nmat 
      complex dmat,coe1,coe2 
c 
       call pmini(nbor, key, 'phonon.in ', 'neighbor n') 
       if(nbor.gt.10) stop ' too many neighbor' 
c 
       tpi=8.d0*datan(1.d0) 
       wvk(1)=(w1*rg(1,1)+w2*rg(1,2)+w3*rg(1,3))*tpi 
       wvk(2)=(w1*rg(2,1)+w2*rg(2,2)+w3*rg(2,3))*tpi 
       wvk(3)=(w1*rg(3,1)+w2*rg(3,2)+w3*rg(3,3))*tpi 
       write(6 ,1001)   wvk(1),wvk(2),wvk(3) 
       write(69,1001)   wvk(1),wvk(2),wvk(3) 
1001   format (5x,/,'   wavevector in XYZ coor. = ', 3f15.7/) 
c 
       do 5 i=1,3*natm 
       do 5 j=1,3*natm 
       ss(i,j)=0. 
       if(i.eq.j) ss(i,j)=1.d0 
       dd(i,j)=0. 
5      continue 
c 
         coex=1.602177d0/1.66054d0*100.d0 
         coex=coex*100.d0/tpi/tpi 
c 
         if(idyn.ne.0) then 
         write(69,1000) 
1000     format(/, 2x,'       dynamics matrix dmat=') 
         write(69 ,*)     '   i    j        real          imag    ' 
         endif 
c 
         do 10 i=1,3*natm 
         ix=mod(i,3) 
         if(ix.eq.0) ix=3 
         ia=(i-ix)/3+1 
         do 20 j=1,i 
         jx=mod(j,3) 
         if(jx.eq.0) jx=3 
         ja=(j-jx)/3+1 
c 
         ita=ityp(ia) 
         itb=ityp(ja) 
         coe0=coex/dsqrt(mass(ita)*mass(itb)) 
         xo=pos(1,ja)-pos(1,ia) 
         yo=pos(2,ja)-pos(2,ia) 
         zo=pos(3,ja)-pos(3,ia) 
         qr=xo*wvk(1)+yo*wvk(2)+zo*wvk(3) 
         coe1=cexp(cmplx(0.0,1.0)*sngl(qr)) 
         dmat=cmplx(0.0,0.0) 
c 
c----> between different atoms 
         if(ia.ne.ja) then 
           do 30 l1=-4,4 
           do 30 l2=-4,4 
           do 30 l3=-4,4 
c 
           call force(ix,ia,jx,ja,l1,l2,l3,nbor,phij) 
           if(dabs(phij).lt.1.d-20) goto 30 
           rlx=l1*g(1,1)+l2*g(1,2)+l3*g(1,3) 
           rly=l1*g(2,1)+l2*g(2,2)+l3*g(2,3) 
           rlz=l1*g(3,1)+l2*g(3,2)+l3*g(3,3) 
           qr=rlx*wvk(1)+rly*wvk(2)+rlz*wvk(3) 
           coe2=cexp(cmplx(0.0,1.0)*sngl(qr)) 
           dmat=dmat + coe0*phij*coe1*coe2 
30         continue 
c 
         else 
c 
c----> between the same atoms 
            do 50 l1=-4,4 
            do 50 l2=-4,4 
            do 50 l3=-4,4 
              if(l1.eq.0.and.l2.eq.0.and.l3.eq.0) goto 40 
              call force(ix,ia,jx,ia,l1,l2,l3,nbor,phij) 
              if(dabs(phij).lt.1.d-20) goto 40 
              rlx=dble(l1)*g(1,1)+dble(l2)*g(1,2)+dble(l3)*g(1,3) 
              rly=dble(l1)*g(2,1)+dble(l2)*g(2,2)+dble(l3)*g(2,3) 
              rlz=dble(l1)*g(3,1)+dble(l2)*g(3,2)+dble(l3)*g(3,3) 
              qr=rlx*wvk(1)+rly*wvk(2)+rlz*wvk(3) 
              coe2=cexp(cmplx(0.0,1.0)*sngl(qr)) 
              dmat=dmat + coe0*phij*(coe2-1.d0) 
40        continue 
c 
           do 55 kp=1,natm 
              if(kp.eq.ia) goto 55 
              call force(ix,ia,jx,kp,l1,l2,l3,nbor,phij) 
              dmat = dmat - coe0*phij 
55         continue 
50         continue 
c 
          endif 
c 
         dd(j,i)=aimag(dmat) 
         dd(i,j)=real(dmat) 
20       continue 
10       continue 
         nmat=3*natm 
c 
       if(idyn.eq.0) return 
       do 70 i=1,nmat 
       do 70 j=1,j 
       write(69 ,'(i5,i5, 6f15.7)')  i, j, dd(i,j),dd(j,i) 
70     continue 
c 
      return 
      end 
c 
      subroutine eigen 
c 
      implicit real*8(a-h,o-z) 
      parameter (neigd=60,nkptd=3000) 
      logical invs,alle,eonly 
c--->    task commons for eigenvalue problem 
      common /debug / ipmt,idat,icel,idyn,ifce,ifun,ieig,ippt, 
     &                iwts,irad,irdp 
      common /eigv/ omcm(neigd),zr(neigd,neigd),zi(neigd,neigd),ne,nblw 
      common /secl/ s(neigd,neigd),h(neigd,neigd),nmat 
      dimension ee(neigd) 
c 
c::: solve secular equation. 
      invs=.false. 
      alle=.false. 
      eonly=.false. 
      ellow=-1000.d0 
      elup=  1000.d0 
c 
      call diagm(ellow,elup,invs,alle,eonly) 
c 
c     take 'square root' of lambda 
c 
      do 120 i=1,ne 
      ee(i)=dsqrt(dabs(omcm(i))) 
      if(omcm(i)) 100,110,110 
100   omcm(i)=-ee(i) 
      goto 120 
110   omcm(i)= ee(i) 
120   continue 
c 
      return 
      end 
c 
       subroutine force(ix,ia,kx,ka,l1,l2,l3,nbor,phij) 
c----------------------------------------------------------------- 
c 
c   find force matrix element PHI_{ix,kx)(0,ia; L,ka) 
c 
c----------------------------------------------------------------- 
      implicit real*8(a-h,o-z) 
      common /debug / ipmt,idat,icel,idyn,ifce,ifun,ieig,ippt, 
     &                iwts,irad,irdp 
      common /struk/g(3,3),omega,natm,pos(3,20),ityp(20),ntyp(20), 
     +              ntype,neq(20),vint,area,a(3) 
      common/rdist/nsh,rshel(500) 
      dimension rik(3) 
c 
      phij=0.d0 
c 
      rlx=l1*g(1,1)+l2*g(1,2)+l3*g(1,3) 
      rly=l1*g(2,1)+l2*g(2,2)+l3*g(2,3) 
      rlz=l1*g(3,1)+l2*g(3,2)+l3*g(3,3) 
      rik(1)=rlx+pos(1,ka)-pos(1,ia) 
      rik(2)=rly+pos(2,ka)-pos(2,ia) 
      rik(3)=rlz+pos(3,ka)-pos(3,ia) 
      rrik=dsqrt(rik(1)**2+rik(2)**2+rik(3)**2) 
      if(rrik.gt.(rshel(nbor)+0.005d0)) return 
      call ffunc(ix,ia,kx,ka,rrik,func,gunc) 
      phij=-rik(ix)*rik(kx)*func-gunc 
      phij=.5d0*phij 
c 
      return 
      end 
c 
      subroutine ffunc(ix,ia,kx,ka,rik,func,gunc) 
      implicit real*8(a-h,o-z) 
c 
c    find F(ik) and G(ik) 
c 
      common /struk/g(3,3),omega,natm,pos(3,20),ityp(20),ntyp(20), 
     +              ntype,neq(20),vint,area,a(3) 
      common /pot a/ naa ,rada(500) ,pairpa(500)  ,caa(500,3), 
     &               dap(500),dapp(500) 
      common /pot b/ nbb ,radb(500) ,pairpb(500)  ,cbb(500,3), 
     &               dbp(500),dbpp(500) 
      common /potab/ nab ,radab(500) ,pairpab(500) ,cab(500,3), 
     &               dabp(500),dabpp(500) 
      common/rdist/nsh,rshel(500) 
c 
      func=0.d0 
      gunc=0.d0 
      if(ityp(ia).ne.1.or.ityp(ka).ne.1) goto 15 
      do 10 is=1,nsh 
      if(dabs(rik-rshel(is)).lt.1.d-5) then 
      f1p = dap(is) 
      f2p = dapp(is) 
      func=(f2p-f1p/rik)/(rik*rik) 
      if(ix.eq.kx) gunc=f1p/rik 
      return 
      endif 
10    continue 
 
15    if(ntype.gt.2) stop ' ntype error' 
      if(ityp(ia).ne.2.or.ityp(ka).ne.2) goto 25 
      do 20 is=1,nsh 
      if(dabs(rik-rshel(is)).lt.1.d-5) then 
      f1p = dbp(is) 
      f2p = dbpp(is) 
      func=(f2p-f1p/rik)/(rik*rik) 
      if(ix.eq.kx) gunc=f1p/rik 
      return 
      endif 
20    continue 
c 
25    do 30 is=1,nsh 
      if(dabs(rik-rshel(is)).lt.1.d-5) then 
      f1p = dabp(is) 
      f2p = dabpp(is) 
      func=(f2p-f1p/rik)/(rik*rik) 
      if(ix.eq.kx) gunc=f1p/rik 
      return 
      endif 
30    continue 
c 
      return 
      end 
c 
      subroutine getppt(elem) 
c----------------------------------------------------------------- 
c 
c  read in parameters of  the unit cell 
c 
c----------------------------------------------------------------- 
      implicit real*8(a-h,o-z) 
      common /debug / ipmt,idat,icel,idyn,ifce,ifun,ieig,ippt, 
     &                iwts,irad,irdp 
      common /struk/g(3,3),omega,natm,pos(3,20),ityp(20),ntyp(20), 
     +              ntype,neq(20),vint,area,a(3) 
      character elem(20)*2,cell*3, filen*10 
      common /pot a/ naa,rada(500),pairpa(500),caa(500,3), 
     &               dap(500),dapp(500) 
      common /pot b/ nbb,radb(500),pairpb(500),cbb(500,3), 
     &               dbp(500),dbpp(500) 
      common /potab/ nab,radab(500),pairpab(500),cab(500,3), 
     &               dabp(500),dabpp(500) 
      common /rdist/nsh,rshel(500) 
c 
      write(6,*)     '   into getppt' 
      write(69,*)    '   into getppt' 
 
      ndim=500 
      call opfile(elem) 
      if(ntype.eq.1) then 
      call rdpot(10, naa, rada,pairpa) 
      call icsccu(rada,pairpa,naa,caa,naa-1,ier) 
      call dcsevu(rada,pairpa,naa,caa,naa-1, 
     &            rshel,dap,nsh,dapp,nsh,ier,ndim) 
c 
      if(ippt.ne.0) then 
        do 10 ii=1,nsh 
        write(6,*)   ii, rshel(ii),dap(ii),dapp(ii) 
10      write(69,*)  ii, rshel(ii),dap(ii),dapp(ii) 
      endif 
c 
      return 
      endif 
c 
      write(6 ,*)  '  ntype=',ntype,ippt 
      write(69,*)  '  ntype=',ntype,ippt 
c 
      if(ntype.gt.2) stop 'ntype' 
      call rdpot(10, naa, rada ,pairpa) 
      call rdpot(20, nbb, radb ,pairpb) 
      call rdpot(30, nab, radab,pairpab) 
      call icsccu(rada,pairpa,naa,caa,naa-1,ier) 
      call dcsevu(rada,pairpa,naa,caa,naa-1, 
     &            rshel,dap,nsh,dapp,nsh,ier,ndim) 
      if(ippt.ne.0) then 
        write(6 ,*)   '  dap, dapp=' 
        write(69,*)   '  dap, dapp=' 
        do 50 ii=1,nsh 
        write(6,*)   ii, rshel(ii),dap(ii),dapp(ii) 
50      write(69,*)  ii, rshel(ii),dap(ii),dapp(ii) 
      endif 
      call icsccu(radb,pairpb,nbb,cbb,nbb-1,ier) 
      call dcsevu(radb,pairpb,nbb,cbb,nbb-1, 
     &            rshel,dbp,nsh,dbpp,nsh,ier,ndim) 
      if(ippt.ne.0) then 
        write(6 ,*)   '  dbp, dbpp=' 
        write(69,*)   '  dbp, dbpp=' 
        do 60 ii=1,nsh 
        write(6,*)   ii, rshel(ii),dbp(ii),dbpp(ii) 
60      write(69,*)  ii, rshel(ii),dbp(ii),dbpp(ii) 
      endif 
       call icsccu(radab,pairpab,nab,cab,nab-1,ier) 
       call dcsevu(radab,pairpab,nab,cab,nab-1, 
     &            rshel,dabp,nsh,dabpp,nsh,ier,ndim) 
      if(ippt.ne.0) then 
        write(6 ,*)   '  dabp, dabpp=' 
        write(69,*)   '  dabp, dabpp=' 
        do 70 ii=1,nsh 
        write(6,*)   ii, rshel(ii),dabp(ii),dabpp(ii) 
70      write(69,*)  ii, rshel(ii),dabp(ii),dabpp(ii) 
      endif 
c 
      return 
      end 
c 
      subroutine iwtset(q1,q2,q3,wk,nk) 
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
c---  set q points and  weight factor from input file 
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
      implicit real*8 (a-h,o-z) 
      parameter (neigd=60,nkptd=3000) 
      character*10 qfile,alph 
      dimension q1(nkptd),q2(nkptd),q3(nkptd),wk(nkptd) 
c 
c--------------------------------------assigned  q file 
      call pmina(qfile,  key, 'phonon.in ','q file     ') 
      open(1,file=qfile) 
      read(1,1050) alph 
      read (1,1005) kn,scale 
      nk=kn 
      weight=0.0 
      do 25 j=1,kn 
      read  (1,1007) q1(j),q2(j),q3(j),wk(j) 
c--- input data are in unit  0.5*(b1,b2,b3) 
c--- change them to unit (b1,b2,b3) 
      q1(j)=0.5d0*q1(j)/scale 
      q2(j)=0.5d0*q2(j)/scale 
      q3(j)=0.5d0*q3(j)/scale 
c--- weight sum 
      weight=weight+wk(j) 
25    continue 
c 
c---> weight normalization 
      do 30 ik=1,kn 
30    wk(ik)=wk(ik)/weight 
c 
      close(1) 
c 
1005  format (i5,f10.6) 
1007  format (5d10.3,2i5) 
1050  format(a10) 
      return 
      end 
c 
      subroutine radist 
c 
c    find the radial shell length 
c 
      implicit real*8(a-h,o-z) 
      common /debug / ipmt,idat,icel,idyn,ifce,ifun,ieig,ippt, 
     &                iwts,irad,irdp 
      common /struk/g(3,3),omega,natm,pos(3,20),ityp(20),ntyp(20), 
     +              ntype,neq(20),vint,area,a(3) 
      common/rdist/nsh,rshel(500) 
      dimension  rad(10000) 
c 
      if(irad.ne.0) then 
      write(6,*)   '    radist' 
      write(69,*)   '    radist' 
      endif 
c 
      ntot=1 
      do 10 i1=-3,3 
      do 10 i2=-3,3 
      do 10 i3=-3,3 
      xx=i1*g(1,1)+i2*g(1,2)+i3*g(1,3) 
      yy=i1*g(2,1)+i2*g(2,2)+i3*g(2,3) 
      zz=i1*g(3,1)+i2*g(3,2)+i3*g(3,3) 
      do 40 nat=1,natm 
      x=pos(1,nat)+xx 
      y=pos(2,nat)+yy 
      z=pos(3,nat)+zz 
      xyz=dsqrt(x**2+y**2+z**2) 
      if(dabs(xyz).lt.1.d-5) goto 40 
        rad(ntot)=xyz 
        ntot=ntot+1 
40    continue 
10    continue 
      ntot=ntot-1 
c 
c--> re-order the rad() 
c 
      ini=1 
50    if(ini.eq.ntot) goto 55 
      do 60 i=ini+1,ntot 
      if(rad(i).lt.rad(ini)) then 
c exchange   ini   and    i 
         rxx=rad(ini) 
         rad(ini)=rad(i) 
         rad(i)=rxx 
      endif 
60    continue 
      ini=ini+1 
      goto 50 
c 
55    continue 
      if(irad.ne.0) then 
      write(6,*) ' ordered rad(i)=' 
      write(6,*)  ' ntot=',ntot 
      write(6,1000) (rad(i),i=1,ntot) 
      write(69,*)  ' ntot=',ntot 
      write(69,*) ' ordered rad(i)=' 
      write(69,1000) (rad(i),i=1,ntot) 
      endif 
1000  format(5f15.7) 
c 
c---> shell no. 
c 
      nsh=0 
      nstep=0 
70    if(nstep.ge.ntot) goto 90 
      rshel(nsh+1)=rad(nstep+1) 
      nsh=nsh+1 
      if(nsh.eq.10) goto 90 
      istep=nstep 
      do 85 is=istep+1,ntot 
      if(dabs(rshel(nsh)-rad(is)).lt.1.d-5) then 
      nstep=nstep+1 
      endif 
85    continue 
      goto 70 
c 
90    continue 
      if(irad.ne.0) then 
      write(6,*) '      rshel(i)=',nsh 
      write(6,1000)    ( rshel(i),i=1,nsh) 
      write(69,*) '      rshel(i)=',nsh 
      write(69,1000)   ( rshel(i),i=1,nsh) 
      endif 
      return 
      end 
c 
      subroutine rdpmtr 
c----------------------------------------------------------------- 
c  read in debug parameter 
c 
c----------------------------------------------------------------- 
c 
      implicit real*8(a-h,o-z) 
      common /debug / ipmt,idat,icel,idyn,ifce,ifun,ieig,ippt, 
     &                iwts,irad,irdp 
c 
c------------debug parameters 
      call pmini(ipmt, key,'phonon.in ' , 'rdpmtr    ') 
      call pmini(idat, key,'phonon.in ' , 'rdata     ') 
      call pmini(icel, key,'phonon.in ' , 'setcel    ') 
      call pmini(idyn, key,'phonon.in ' , 'dynmat    ') 
      call pmini(ifce, key,'phonon.in ' , 'force     ') 
      call pmini(ifun, key,'phonon.in ' , 'ffunc     ') 
      call pmini(ieig, key,'phonon.in ' , 'eigen     ') 
      call pmini(ippt, key,'phonon.in ' , 'getppt    ') 
      call pmini(iwts, key,'phonon.in ' , 'iwtset    ') 
      call pmini(irad, key,'phonon.in ' , 'radist    ') 
      call pmini(irdp, key,'phonon.in ' , 'rdpot     ') 
c 
      return 
      end 
c 
      subroutine opfile(elem) 
      implicit real*8(a-h,o-z) 
      common /struk/g(3,3),omega,natm,pos(3,20),ityp(20),ntyp(20), 
     +              ntype,neq(20),vint,area,a(3) 
      character elem(20)*2,cell*3, filen*10 
c 
      write(6,*)    '   into opfile' 
      call pmina(cell, key, 'phonon.in ', 'strk type ') 
c 
      if(ntype.eq.1) then 
      filen(1:3)=cell 
      filen(4:6)='pt.' 
      filen(7:9)=elem(1) 
      write(6,*)    '  filen=',filen 
      open(10,file=filen) 
      return 
      endif 
c 
      if(ntype.gt.2)  stop ' ntype' 
      filen(1:2)='pt' 
      filen(3:4)=elem(1) 
      filen(5:6)=elem(2) 
      filen(7:7)='.' 
      filen(8:10)=cell 
      open(30,file=filen) 
c 
      call pmina(cell, key, 'phonon.in ', 'strk tp 1 ') 
c 
      filen(1:3)=cell 
      filen(4:6)='pt.' 
      filen(7:9)=elem(1) 
      open(10,file=filen) 
c 
      call pmina(cell, key, 'phonon.in ', 'strk tp 2 ') 
c 
      filen(1:3)=cell 
      filen(4:6)='pt.' 
      filen(7:9)=elem(2) 
      open(20,file=filen) 
c 
      return 
      end 
c 
        subroutine rdpot(ntap, naa, rad, pairp) 
        implicit real*8(a-h,o-z) 
        dimension rad(500),pairp(500) 
        common /debug / ipmt,idat,icel,idyn,ifce,ifun,ieig,ippt, 
     &                  iwts,irad,irdp 
c 
        rewind ntap 
c 
        read(ntap,'(a1)')  dummy 
        read(ntap,'(a1)')  dummy 
        read(ntap,'(a1)')  dummy 
        read(ntap,'(a1)')  dummy 
        read(ntap,'(a1)')  dummy 
        read(ntap,'(a1)')  dummy 
c 
        naa=0 
        do 10 i=1,500 
        read(ntap, 1001,end=15) i,rad(i),pairp(i) 
        if(irdp.ne.0) then 
          if(mod(i,10).eq.0) write(6 ,1001) i,rad(i),pairp(i) 
                             write(69,1001) i,rad(i),pairp(i) 
        endif 
        naa=naa+1 
10      continue 
15      continue 
c 
1001    format(5x,i5,10E25.14) 
c 
        close(ntap) 
c 
        return 
        end 
c 
      subroutine setcel(elem,mass) 
c----------------------------------------------------------------- 
c 
c  read in parameters of  the unit cell 
c 
c----------------------------------------------------------------- 
      implicit real*8(a-h,o-z) 
      real*8 mass(20) 
      common /debug / ipmt,idat,icel,idyn,ifce,ifun,ieig,ippt, 
     &                iwts,irad,irdp 
      common /struk/g(3,3),omega,natm,pos(3,20),ityp(20),ntyp(20), 
     +              ntype,neq(20),vint,area,a(3) 
      common /finst1/rg(3,3) 
      character elem(20)*2,cell*3, filecell*10 
      call pmina(cell, key, 'phonon.in ', 'strk type ') 
      filecell(1:7)='percel.' 
      filecell(8:10)=cell 
      write(6,*) '    filecell=',filecell 
c 
      call pmini(natm, key, filecell, 'atom No   ') 
      ntype=0 
      do 30 j=1,natm 
        call pminia(jat, key, filecell, 'type(atom ', j) 
        if(ntype.lt.jat) ntype=jat 
30    continue 
      do 50 j=1,ntype 
      nq=0 
      do 40 k=1,natm 
       call pminia(iatom, key, filecell, ' type(atom', k) 
        if(iatom.eq.j) then 
         nq=nq+1 
         ityp(k)=j 
        ntyp(k)=nq 
       endif 
40    continue 
      neq(j)=nq 
50    continue 
c 
c----------------------- element name of atom types 
      do 60 it=1,ntype 
      call pminaa(elem(it), key, 'phonon.in ', 'elem(type ',it) 
      call pminra(mass(it), key, 'phonon.in ', 'mass(type ',it) 
60    continue 
c 
c--------------------------------------lattice 
      na=3 
      do 100 l=1,na 
      call pminra(g(1,l), key, filecell, 'x(lat vec ', l) 
      call pminra(g(2,l), key, filecell, 'y(lat vec ', l) 
      call pminra(g(3,l), key, filecell, 'z(lat vec ', l) 
      a(l)=sqrt(g(1,l)*g(1,l)+g(2,l)*g(2,l)+g(3,l)*g(3,l)) 
100   continue 
c 
      do 200 jat=1,natm 
      j =ityp(jat) 
      ii=ntyp(jat) 
      call pminra(pos(1,jat), key, filecell, 'x(atom     ',jat) 
      call pminra(pos(2,jat), key, filecell, 'y(atom     ',jat) 
      call pminra(pos(3,jat), key, filecell, 'z(atom     ',jat) 
c 
      if(icel.ne.0) write(69,1005) j,ii, (pos(i,jat),i=1,3) 
                    write(6 ,1005) j,ii, (pos(i,jat),i=1,3) 
200    continue 
c 
c---> find the recip. lattice 
      x12=g(2,1)*g(3,2)-g(3,1)*g(2,2) 
      y12=g(3,1)*g(1,2)-g(1,1)*g(3,2) 
      z12=g(1,1)*g(2,2)-g(2,1)*g(1,2) 
      area=dsqrt(x12*x12+y12*y12+z12*z12) 
      omega=g(1,3)*x12+g(2,3)*y12+g(3,3)*z12 
      write(6 ,1003) (a(i),i=1,3),((g(i,j),i=1,3),j=1,3),omega 
      write(69,1003) (a(i),i=1,3),((g(i,j),i=1,3),j=1,3),omega 
      rg(1,1)=(g(2,2)*g(3,3)-g(3,2)*g(2,3))/omega 
      rg(2,1)=(g(3,2)*g(1,3)-g(1,2)*g(3,3))/omega 
      rg(3,1)=(g(1,2)*g(2,3)-g(2,2)*g(1,3))/omega 
 
      rg(1,2)=(g(2,3)*g(3,1)-g(3,3)*g(2,1))/omega 
      rg(2,2)=(g(3,3)*g(1,1)-g(1,3)*g(3,1))/omega 
      rg(3,2)=(g(1,3)*g(2,1)-g(2,3)*g(1,1))/omega 
 
      rg(1,3)=(g(2,1)*g(3,2)-g(3,1)*g(2,2))/omega 
      rg(2,3)=(g(3,1)*g(1,2)-g(1,1)*g(3,2))/omega 
      rg(3,3)=(g(1,1)*g(2,2)-g(2,1)*g(1,2))/omega 
c 
      write(6 ,1002) ((rg(i,j),i=1,3),j=1,3) 
      write(69,1002) ((rg(i,j),i=1,3),j=1,3) 
c 
1002  format(10x,'reciprocal lattice:'/15x,'b1(x,y,z)=',3f12.7/ 
     +       15x,'b2(x,y,z)=',3f12.7/15x,'b3(x,y,z)=',3f12.7) 
1003  format(10x,'length(a.u.) of a1,a2,a3='/15x,3f12.7/ 
     +       10x,'primitive translations:'/15x,'a1(x,y,z)=',3f12.7/ 
     +       15x,'a2(x,y,z)=',3f12.7/15x,'a3(x,y,z)=',3f12.7/ 
     +       10x,'unit cell volume (a.u.**3) is  ',f10.4) 
1005  format(5x,i2,'-th type, ',i2,'-th atoms at(unit):',5x,3f10.6) 
      return 
      end 
