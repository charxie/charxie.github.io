      subroutine subdiag(natom,ndim,t,nstp,aring,nring)

c***********************************************************************
c constructs orthogonal matrix T for transformation from AO's to
c natural hybrid bond orbitals, using input density matrix DM
c maximum number of aromatic rings allowed = 10
c***********************************************************************      

      implicit real*8 (a-h,o-z)
      parameter(neigd=500,natmax=200,ndmax=500)
      character segname*4,resname*3,bname*2,ac*4,symbol*2
      logical invs,alle,eonly,write,aring,restart,dotsa
      integer ul
      common/ato/ac(natmax),symbol(40)
      common/atom2/exps(40),exps2(40),expp(40),expp2(40),expd(40),              
     *expd2(40),expf(40),expf2(40),cs1(40),cs2(40),cp1(40),cp2(40),             
     *cd1(40),cd2(40),cf1(40),cf2(40),couls(40),coulp(40),could(40),            
     *coulf(40),x(natmax),y(natmax),z(natmax),ires(natmax),
     *resname(natmax),segname(natmax)
      common/optv/rho,fermi,wstart,wend,
     *lb,le,lbr,ler,iovlb,iovle,lcharg,neltr,
     *distl,idistn,ieb,iee,isb,ise,ibb,ibe,ipope,jpope,irest 
      common /eigv/ omcm(neigd),zr(neigd,neigd),zi(neigd,neigd),ne,nblw
      common /secl/ s(neigd,neigd),h(neigd,neigd),nmat
      common/option/restart,dotsa,nstep,unit
      common/lbl/bname(ndmax,2),label(ndmax,3),ibx(ndmax)
      common/lblaro/labar(ndmax,7)
      common/info1/ul(natmax),ll(natmax),nele,nocc
      common/q/q(4,ndmax)
      common/pol/pol(ndmax,3),iathy(ndmax,3),ino(natmax)
      common/iounit/iin,iout,istart(40),iungf,iunhl
      common/cring/idc6(6,10),ovar(6,6,10)
      dimension t(ndim,ndim)
      dimension dblk(8,8),c(8,8),eval(8),borb(8),sblk(8,8),hblk(8,8)
      dimension sss(8,8),hhh(8,8)
      character*2 name(3)
      character*1 star,blnk
      data name,star,blnk/'LP','BD','PB','*',' '/
      data cutoff/2.d0/
      
      invs=.false.
      alle=.false.
      eonly=.false.
      ellow=-1000.d0
      elup=  1000.d0
      write=.false.
      nele=neltr
      if(mod(nele,2).ne.0) stop 'not close-shell system'
      nocc=nele/2
      
c zero arrays
      do i = 1 , ndim
         do k = 1 , 3
            pol(i,k)=0.0
            iathy(i,k)=0
            q(k,i)=0.0
         enddo
         q(4,i)=0.0
         do j = 1 , ndim
            t(i,j)=0.0
         enddo
      enddo
      do i = 1 , natom
         ino(i)=0
      enddo
      ibd=0
      na1=natom+1

      write(iout,*) ' search bonds by diagonalizing 1,2-center sublocks'

      do ib1 = 1 , na1
      ib=ib1-1
      do 90 ic1 = 2 , na1
         ic=ic1-1
         if(ib.ne.0) goto 40

c>>> lone pairs

      nctr=1
      iat1=ic
      iat2=0

      goto 70

c>>> bond pairs

   40 nctr=2
      iat1=ib
      iat2=ic
      if(iat2.le.iat1) goto 90

   70 continue
      if(ibd.ge.nocc) goto 120

      if(iat2.eq.0) goto 71
      if((((x(iat1)-x(iat2))**2.0)+
     *    ((y(iat1)-y(iat2))**2.0)+
     *    ((z(iat1)-z(iat2))**2.0)).gt.cutoff**2) goto 90
   71 continue

      call load(1,iat1,iat2,dblk,sblk,nb)
      call load(2,iat1,iat2,hblk,sblk,nb)

      write(iout,*)'iat1, iat2', iat1, iat2
      
c diagonalize the subblocks

      do i = 1 , nb
         do j = 1 , nb
            s(i,j)=0.0
            if(i.eq.j) s(i,j)=1.0
            h(i,j)=0.0
         enddo
      enddo

      do i = 1 , nb
         do j = 1, i
            h(j,i)=0.0
            h(i,j)=dblk(i,j)
c            s(i,j)=sblk(i,j)
         enddo
      enddo
      
      if(.not.write) goto 256
      write(iout,*) ' subblock of the density matrix'
	      do i = 1 , nb
	         write(iout,'(500f8.4)') (dblk(i,j),j=1,nb)
	      enddo
	      write(iout,*)
      write(iout,*) ' subblock of the hamiltonian matrix'
	      do i = 1 , nb
	         write(iout,'(500f8.4)') (hblk(i,j),j=1,nb)
	      enddo
	      write(iout,*)
      write(iout,*) ' subblock of the overlap matrix'	      
	      do i = 1 , nb
	         write(iout,'(500f8.4)') (sblk(i,j),j=1,nb)
	      enddo
  256 continue	      

      nmat=nb
      call diagm(ellow,elup,invs,alle,eonly)
        
      do k = 1 , nb
         do l = 1 , nb
            sss(k,l)=0.0
            hhh(k,l)=0.0
            do i = 1 , nb
               do j = 1 , nb
                  sss(k,l)=sss(k,l)+zr(i,k)*sblk(i,j)*zr(j,l)
                  hhh(k,l)=hhh(k,l)+zr(i,k)*hblk(i,j)*zr(j,l)
               enddo
            enddo
         enddo
      enddo

c      if(
c     :      (ac(iat1).eq.' CG '.and.ac(iat2).eq.' CD1 ')
c     :  .or.(ac(iat1).eq.' CG '.and.ac(iat2).eq.' CD2 ')
c     :  .or.(ac(iat1).eq.' CD1'.and.ac(iat2).eq.' CE1 ')
c     :  .or.(ac(iat1).eq.' CD2'.and.ac(iat2).eq.' CE2 ')
c     :  .or.(ac(iat1).eq.' CE1'.and.ac(iat2).eq.' CZ  ')
c     :  .or.(ac(iat1).eq.' CE2'.and.ac(iat2).eq.' CZ  ')) then
c            omcm(8)=0.0
c      endif
      emax=0.0
      do i = 1 , nb
         eval(i)=omcm(i)*sss(i,i)
         if(eval(i).gt.emax) emax=eval(i)
      enddo
c      write(iout,'(a5,f8.4)') 'emax=',emax
      
      do i = 1 , nb
         do j = 1 , nb
            c(i,j)=zr(i,j)
         enddo
      enddo
      
      if(.not.write) goto 512
      write(iout,*)' eigenvalue of this subblock'
      write(iout,'(i10,12f8.4)') nb,(omcm(i),i=1,nb)
      write(iout,'(i10,12f8.4)') nb,(sss(i,i),i=1,nb)
      write(iout,'(i10,12f8.4)') nb,(eval(i),i=1,nb)
c      write(iout,*)' eigenvectors'
c      do i = 1 , nb
c         write(iout,'(12f8.4)') (zr(i,j),j=1,nb)
c      enddo
      write(iout,*)
 512  continue

      do irnk = 1 , nb

          call extrct(c,eval,borb,occ,nb,irnk)
          if(occ.gt.2.5) stop 'local occupancy too large!'
          if(iat2.eq.0.and.
     :      (ac(iat1).eq.' N  '.or.ac(iat1).eq.' NT ')) then
             if(2.0-occ.ge.0.8) goto 85
          elseif(
     :          (ac(iat1).eq.' CG '.and.ac(iat2).eq.' CD1 ')
     :      .or.(ac(iat1).eq.' CG '.and.ac(iat2).eq.' CD2 ')
     :      .or.(ac(iat1).eq.' CD1'.and.ac(iat2).eq.' CE1 ')
     :      .or.(ac(iat1).eq.' CD2'.and.ac(iat2).eq.' CE2 ')
     :      .or.(ac(iat1).eq.' CE1'.and.ac(iat2).eq.' CZ  ')
     :      .or.(ac(iat1).eq.' CE2'.and.ac(iat2).eq.' CZ  ')) then
c             if(2.0-occ.ge.0.5.or.occ-2.0.gt.0.5) goto 85
             if(dabs(occ-emax).gt.1.0d-4) goto 85
          else
             if(2.0-occ.ge.0.2.or.occ-2.0.gt.0.5) goto 85
          endif
          write(iout,'(2x,a8,f20.10)') ' occ  = ',occ
          if(write) then
             write(iout,*) ' local overlap matrix SSS'	      
	     do i = 1 , nb
                write(iout,'(500f8.4)') (sss(i,j),j=1,nb)
	     enddo
	     write(iout,*)
             write(iout,*) ' local hamiltonian matrix HHH'	      
	     do i = 1 , nb
                write(iout,'(500f8.4)') (hhh(i,j),j=1,nb)
	     enddo
	     write(iout,*)
	  endif
c          write(iout,'(8f8.4)') (borb(i),i=1,nb)
          ibd=ibd+1
          if(nctr.eq.1) call deplet(borb,occ,nb,iat1,iat2)
	  call stash(borb,occ,ibd,iat1,iat2)
	  bname(ibd,1)(1:2)=name(nctr)
	  bname(ibd,2)(1:1)=blnk
          label(ibd,1)=irnk
	  label(ibd,2)=iat1
	  label(ibd,3)=iat2

      enddo

   85 continue

   90 continue
   
      enddo

  120 continue
  
      write(iout,*) ibd
  
      aring=.false.
      nring=0
      do iat = 1 , natom
         if(resname(iat).eq.'BEN'.or.
     :      resname(iat).eq.'TYR'.or.
     :      resname(iat).eq.'PHE') then
            aring=.true.
            if(ac(iat).eq.' CG ') then 
               nring=nring+1
               idc6(1,nring)=iat
            endif
            if(ac(iat).eq.' CD1') then 
               idc6(2,nring)=iat
            endif
            if(ac(iat).eq.' CD2') then 
               idc6(3,nring)=iat
            endif
            if(ac(iat).eq.' CE1') then 
               idc6(4,nring)=iat
            endif
            if(ac(iat).eq.' CE2') then 
               idc6(5,nring)=iat
            endif
            if(ac(iat).eq.' CZ ') then 
               idc6(6,nring)=iat
            endif
         endif
         if(resname(iat).eq.'TRP') then
            aring=.true.
            if(ac(iat).eq.' CD2') then 
               nring=nring+1
               idc6(1,nring)=iat
            endif
            if(ac(iat).eq.' CE2') then 
               idc6(2,nring)=iat
            endif
            if(ac(iat).eq.' CE3') then 
               idc6(3,nring)=iat
            endif
            if(ac(iat).eq.' CZ2') then 
               idc6(4,nring)=iat
            endif
            if(ac(iat).eq.' CZ3') then 
               idc6(5,nring)=iat
            endif
            if(ac(iat).eq.' CH2') then 
               idc6(6,nring)=iat
            endif
         endif
      enddo
      if(aring) then      
         write(iout,*)
         write(iout,*) 'Handling aromatic rings ...'
         write(iout,*)
      endif
      write(iout,*)'The biomolecule contains',nring,' aromatic rings'
      if(nring.gt.10) stop 'too many aromatic rings'
      if(aring) then
         do iring = 1 , nring
            call fourth(ibd,iring)
            do i = ibd-2,ibd
               bname(i,1)=name(3)
               bname(i,2)=blnk
               labar(i,1)=iring
               labar(i,2)=idc6(1,iring)
               labar(i,3)=idc6(2,iring)
               labar(i,4)=idc6(3,iring)
               labar(i,5)=idc6(4,iring)
               labar(i,6)=idc6(5,iring)
               labar(i,7)=idc6(6,iring)
            enddo
         enddo
      endif
      
      write(iout,*) ibd
      
      ibo=ibd
      do 150 i = 1 , ibo
         if(bname(i,1).eq.name(1)) goto 150
	 ibd=ibd+1
	 do j = 1 , 3
	    label(ibd,j)=label(i,j)
	 enddo
	 if(aring) then
	    do j = 1 , 7
	       labar(ibd,j)=labar(i,j)
	    enddo
	 endif
	 bname(ibd,1)=bname(i,1)
	 bname(ibd,2)=star
  150 continue

      if(ibd.eq.ndim) goto 160
      write(*,*) ibd,ndim
      stop 'Lewis structure error: ibd != ndim'
  160 continue

      call reform(t,natom,ndim,nring)
c      if(aring) call getovar(ndim,nring)
      if(aring) call reformaro(t,ndim,nring)
      call anlyze(t,ndim)

      return
      end


      subroutine load(isel,iat1,iat2,dblk,sblk,nb)

c***********************************************************************      
c zeros the 8x8 matrices and loads in atomic blocks of density matrix
c and overlap matrix for the atoms iat1 and iat2
c***********************************************************************

      implicit real*8 (a-h,o-z)
      parameter(natmax=200,ndmax=500,nsmp=10)
      integer ul
      dimension dblk(8,8),sblk(8,8),iat(2)
      common/info1/ul(natmax),ll(natmax),nele,nocc
      common/inst/ovlp(ndmax,ndmax),hami(ndmax,ndmax),dens(ndmax,ndmax)

      iat(1)=iat1
      iat(2)=iat2
      do i = 1 , 8
         do j = 1 , 8
            dblk(i,j)=0.0
            sblk(i,j)=0.0
         enddo
      enddo
      nrow=0
      ncol=0

      do 50 i = 1 , 2
         ia=iat(i)
         if(ia.eq.0) goto 50
         iu=ul(ia)
         il=ll(ia)
         do irow = il , iu
            nrow=nrow+1
            ncol=0
            do 30 j = 1 , 2
               ja=iat(j)
               if(ja.eq.0) goto 30
               ju=ul(ja)
               jl=ll(ja)
               do icol = jl , ju
                  ncol=ncol+1
                  if(isel.eq.1)dblk(nrow,ncol)=dens(irow,icol)
                  if(isel.eq.2)dblk(nrow,ncol)=hami(irow,icol)
                  sblk(nrow,ncol)=ovlp(irow,icol)
               enddo
   30       continue
         enddo
   50 continue
      nb=nrow

      return
      end


      subroutine ulll(natom,ndim)
      
c***********************************************************************      
c find the orbit index upper and lower bounds ul and ll
c***********************************************************************
      
      implicit real*8 (a-h,o-z)
      parameter(natmax=200,ndmax=500)
      character indt*10,ac*4,symbol*2
      integer ul
      common/info1/ul(natmax),ll(natmax),nele,nocc
      common/ato/ac(natmax),symbol(40)
      common/aoout/indt(1000,2),iatom(natmax),lorb(1000)

      do na = 1 , natom
         ll(na)=ndim
         ul(na)=0
         do k = 1 , ndim
	    if(indt(k,1)(1:2).eq.ac(na)(1:2)
     *	    .and.lorb(k).eq.iatom(na)) then
	       if(k.gt.ul(na)) ul(na)=k
	       if(k.lt.ll(na)) ll(na)=k	     
     	    endif
	 enddo
c	 write(*,'(3i10)') na,ll(na),ul(na)
      enddo
      
      return
      end