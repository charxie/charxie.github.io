      subroutine stash(borb,occ,ibd,iat1,iat2)

c*************************************************************************
c seperate bond orbital BORB into polarization coefficients (for storage
c in POL) and hybrids (for storage in Q), keeping track of the number 
c INO(NA) of hybrids accumulated for each atom NA
c*************************************************************************

      implicit real*8 (a-h,o-z)
      parameter (natmax=200,ndmax=500)
      integer ul
      dimension borb(8),iat(2),hyb(4)
      common/info1/ul(natmax),ll(natmax),nele,nocc
      common/q/q(4,ndmax)
      common/pol/pol(ndmax,3),iathy(ndmax,3),ino(natmax)
      common/iounit/iin,iout,istart(40),iungf,iunhl

      iat(1)=iat1
      iat(2)=iat2
      kmax=0
      iocc=1
      if(occ.lt.0.1) iocc=-1
      do 40 i = 1 , 2
         ia=iat(i)
         if(ia.eq.0) goto 40
         nu=ul(ia)
         nl=ll(ia)
c...extract hybrid from bond orbital for atom IA
         kmin=kmax+1
         kmax=kmax+nu-nl+1
         mj=0
         do k = kmin, kmax
            mj=mj+1
            hyb(mj)=borb(k)
         enddo
c...extract polarization coefficients, store in POL
         psq=0.0
         do j = 1 , mj
            psq=psq+hyb(j)**2
         enddo
         sgn=1.0
         if(hyb(1).lt.0.0) sgn=-sgn
         p=sgn*dsqrt(psq)
         pol(ibd,i)=p
c...enter hybrid in appropriate column of Q matrix if BORB occupied
         if(iocc.le.0) goto 40
c...one more hybrid for atom IA
         ino(ia)=ino(ia)+1
         ncol=ll(ia)+ino(ia)-1
c...place normalized hybrid in appropriate block of Q
         nh=nu-nl+1
         do nrow = 1 , nh
            q(nrow,ncol)=hyb(nrow)/p
         enddo
         iathy(ibd,i)=ino(ia)
   40 continue

      return
      end

      subroutine reform(t,natom,ndim,nring)

c*************************************************************************
c Builds the final orthogonal transformation matrix T from polarization
c coefficients (POL) and atomic hybrids(stored in Q, identified in IATHY),
c inserting antibonds as necessary
c*************************************************************************

      implicit real*8 (a-h,o-z)
      parameter(natmax=200,ndmax=500,neigd=500)
      integer ul
      character bname*2
      logical invs,alle,eonly,write
      dimension t(ndim,ndim)
      dimension s(4,4),eval(4),c(4,4),ta(4,4)
      common/lbl/bname(ndmax,2),label(ndmax,3),ibx(ndmax)
      common/info1/ul(natmax),ll(natmax),nele,nocc
      common/q/q(4,ndmax)
      common/pol/pol(ndmax,3),iathy(ndmax,3),ino(natmax)
      common/iounit/iin,iout,istart(40),iungf,iunhl
      common/eigv/omcm(neigd),zr(neigd,neigd),zi(neigd,neigd),ne,nblw
      common/secl/s1(neigd,neigd),h1(neigd,neigd),nmat
      common/cring/idc6(6,10),ovar(6,6,10)

      invs=.false.
      alle=.false.
      eonly=.false.
      ellow=-1000.d0
      elup = 1000.d0
      write=.false.
      
c reorder occupied BO's to put lone pairs last
         
      do nlp = 1 , nocc
         if(bname(nlp,1).ne.'LP') goto 20
      enddo
   20 nlp=nlp-1
      nbd=nocc-nlp
      do 50 ibd = 1 , ndim
         if(ibd.gt.nlp) goto 30
         ibx(ibd)=ibd+nbd
         goto 50
   30    if(ibd.gt.nocc) goto 40
         ibx(ibd)=ibd-nlp
         goto 50
   40    ibx(ibd)=ibd
   50 continue
      do i = 1 , ndim
         do j = 1 , i
            t(i,j)=0.0
            t(j,i)=0.0
         enddo
      enddo

c symmetric orthogonalization of the atomic hybrids

      do 170 ia = 1 , natom

         il=ll(ia)
         iu=ul(ia)
         nh=iu-il+1
         if(nh.eq.1) goto 170
c load IA-block of Q into TA...
         do j = 1 , nh
            do i = 1 , 4
               ta(i,j)=q(i,il+j-1)
            enddo
         enddo
c form overlap matrix S = TA(Transp)*TA...
         do j = 1 , nh
            do i = j , nh
               temp=0.0
               do k = 1 , nh
                  temp=temp+ta(k,i)*ta(k,j)               
               enddo
               s(i,j)=temp
               s(j,i)=temp
            enddo
         enddo
c diagonalize the overlap matrix...
         do i = 1 , 4
            do j = 1 , 4
               s1(i,j)=0.0
               if(i.eq.j) s1(i,j)=1.0
            enddo
         enddo
         do i = 1 , 4
            do j = 1 , i
               h1(j,i)=0.0
               h1(i,j)=s(i,j)
            enddo
         enddo
         nmat=4
         call diagm(ellow,elup,invs,alle,eonly)
         do i = 1 , 4
            eval(i)=omcm(i)
         enddo
         do i = 1 , 4
            do j = 1 , 4
               c(i,j)=zr(i,j)
            enddo
         enddo
c form inverse square root of S, store in S...
         do i = 1 , nh
            eval(i)=1.0/dsqrt(eval(i))
         enddo
         do j = 1 , nh
            do i = j , nh
               temp=0.0
               do k = 1 , nh
                  temp=temp+eval(k)*c(i,k)*c(j,k)
               enddo
               s(i,j)=temp
               s(j,i)=temp
            enddo
         enddo
c form new TA=TA*S^{-1/2}, store in C...
         do j = 1 , nh
            do i = 1 , nh
               temp=0.0
               do k = 1 , nh
                  temp=temp+ta(i,k)*s(k,j)
               enddo
               c(i,j)=temp
            enddo
         enddo
c replace orthogonalized TA in array Q...
         do j = 1 , nh
            do i = 1 , 4
               q(i,il+j-1)=c(i,j)
            enddo
         enddo

  170 continue
  
c  check the results of orthogonalization
      if(nring.eq.0) goto 2000
      if(write) write(iout,*)
      if(write) write(iout,*)'Orthogonalized hybrids'
      do iring = 1 , nring
         do iat = 1 , 6
            ilc=ll(idc6(iat,iring))
            iuc=ul(idc6(iat,iring))
            nhc=iuc-ilc+1
            if(write) write(iout,'(10i10)') 
     :                ino(idc6(iat,iring)),idc6(iat,iring)
            do iq = 1 , 4
               if(write) write(iout,'(10f10.5)')
     :              (q(iq,ilc+jq-1),jq=1,ino(idc6(iat,iring)))
            enddo
            do i = 1 , 4
            do j = i , 4
            dotp=0.0
            do k = 1 , 4
            dotp=dotp+q(k,ilc+i-1)*q(k,ilc+j-1)
            enddo
            if(write) 
     :      write(iout,'(1x,a1,i1,a1,i1,a1,f10.5)')'(',i,'*',j,')',dotp
            enddo
            enddo
         enddo
      enddo
2000  continue      

c symmetric orthogonalization complete. deposit final bond orbitals in T

      nbo=0
      do ibd = 1 , ndim
         kbd=ibd
         if(bname(ibd,2)(1:1).ne.'*') goto 250
c search occupied orbitals list to get proper hybrids
c search occupied bond orbitals for match with antibond atoms
c negative irnk = label(k,1) means bond orbital was already used	
         do 240 k = 1 , nbo
            do i = 2 , 3
               if(label(k,i).ne.label(ibd,i)) goto 240
               if(label(k,1).le.0.and.bname(k,1).eq.'BD') goto 240
            enddo
c found match; set label(k,1) < 0 and reset polarization parameters
c for antibond
            kbd=k
            label(kbd,1)=-label(kbd,1)
            pol(ibd,2)=-pol(kbd,1)
            pol(ibd,1)= pol(kbd,2)
            goto 250
  240    continue
c couldnt find successful match...exit
         stop 'error:reform'
c deposit bond orbitals in T matrix
  250    continue
         do 270 i = 1 , 2
            ia=label(ibd,i+1)
            if(ia.eq.0) goto 270
            jl=ll(ia)
            ju=ul(ia)
            irow=0
            icol=jl+iathy(kbd,i)-1
            do j = jl , ju
               irow=irow+1
               t(j,ibx(ibd))=pol(ibd,i)*q(irow,icol)
            enddo
  270    continue
         if(ibd.eq.kbd) nbo=ibd
      enddo

c restore label(i,1) > 0

      do i = 1 , ndim
         if(label(i,1).lt.0) label(i,1)=-label(i,1)
      enddo
      
      return
      end