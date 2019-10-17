c fix the signs for the natural hybrids to avoid the phase problem
c a comparison array is obtained from a reference structure.
c a sign array which carries the orientational information of the bonds
c is created and compared with the signs stored in the comparison array.

      SUBROUTINE compsignaro(n,ipr)  
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      character*4 symbol,resname,segname
      character*2 bname
      INTEGER CZ,U,UL
      common/INFO1/CZ(NAZR),U(NBSZR),UL(NAZR),LL(NAZR),
     :NELECS,NOCCA,NOCCB
      common/LBL/bname(nbszr,2),label(NBSZR,3),IBX(NBSZR) 
      common/lblaro/labar(nbszr,7)
      common/EXTRA/T(NBSZR,NBSZR),IHYB(NAZR),TOPO(NAZR,NAZR),
     :T2J(NBSZR)
      common/cring/idc6(6,20)
      common/iden/idres(nazr),symbol(nazr),resname(nazr),segname(nazr),
     :rx(nazr),ry(nazr),rz(nazr)
      dimension infosign6(2,6,120)

      npb=0
      nlp=0
      nbd=0
      do i = 1 , n
         if(bname(i,1).eq.'PB') npb=npb+1
         if(bname(i,1).eq.'LP') nlp=nlp+1
         if(bname(i,1).eq.'BD') nbd=nbd+1      
      enddo
      if(npb+nlp+nbd.ne.n.or.mod(npb,6).ne.0.or.mod(nbd,2).ne.0) 
     :   stop 'error in fixsign.f'
      nring=npb/6
      i1beg=nbd/2+1
      i1end=nbd/2+npb/2
      i2beg=i1end+nlp+nbd/2+1
      i2end=n
      
      open(16,file='tsign')
      
      iaro=0
      do i = 1 , n
         if(i.lt.i1beg) then
            read(16,'(i10)') idummy
         elseif(i.ge.i1beg.and.i.le.i1end) then
            iaro=iaro+1
            read(16,'(i10,12x,6(7x,2i5,f8.4))') 
     :       idummy,(infosign6(1,iat,iaro),
     :       infosign6(2,iat,iaro),tmax0,iat=1,6)
         elseif(i.gt.i1end.and.i.lt.i2beg) then
            read(16,'(i10)') idummy
         else
            iaro=iaro+1
            read(16,'(i10,12x,6(7x,2i5,f8.4))') 
     :       idummy,(infosign6(1,iat,iaro),
     :       infosign6(2,iat,iaro),tmax0,iat=1,6)
         endif
      enddo
      if(idummy.ne.n.or.iaro.ne.npb) stop 'error in fixsign.f'
      close(16)
      
c      do ir = 1 , iaro
c         write(6,'(50i5)') (infosign6(1,iat,ir),iat=1,6)
c         write(6,'(50i5)') (infosign6(2,iat,ir),iat=1,6)
c      enddo

      iaro=0
      do i = 1 , n
         do ibd = 1 , n
            if(ibx(ibd).eq.i) goto 105
         enddo
 105     if(bname(ibd,1).eq.'PB') then
           iaro=iaro+1
	   do iat = 1 , 6
	      if(infosign6(1,iat,iaro).ne.0) then
	         tmax6=t(infosign6(1,iat,iaro),i)
                 imax=int(tmax6/(dabs(tmax6)-0.000001))  
                 itemp=imax*infosign6(2,iat,iaro)
                 if(itemp.eq.-1) then
                    do j = ll(idc6(iat,labar(ibd,1))),
     :                     ul(idc6(iat,labar(ibd,1)))
                       t(j,i)=-t(j,i)
                    enddo
                 endif
              endif
	   enddo
	 endif
      enddo

      return                                                                    
      END                                                                       

      SUBROUTINE compsign(n,ipr)  
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*2 bname
      character*4 symbol,resname,segname
      INTEGER CZ,U,UL
      common/INFO1/CZ(NAZR),U(NBSZR),UL(NAZR),LL(NAZR),
     $NELECS,NOCCA,NOCCB
      common/LBL/bname(nbszr,2),label(NBSZR,3),IBX(NBSZR) 
      common/EXTRA/T(NBSZR,NBSZR),IHYB(NAZR),TOPO(NAZR,NAZR),
     $T2J(NBSZR)
      common/iden/idres(nazr),symbol(nazr),resname(nazr),segname(nazr),
     :rx(nazr),ry(nazr),rz(nazr)
      dimension infosign(2,nbszr)

      npb=0
      nlp=0
      nbd=0
      do i = 1 , n
         if(bname(i,1).eq.'PB') npb=npb+1
         if(bname(i,1).eq.'LP') nlp=nlp+1
         if(bname(i,1).eq.'BD') nbd=nbd+1      
      enddo
      if(npb+nlp+nbd.ne.n.or.mod(npb,6).ne.0.or.mod(nbd,2).ne.0) 
     :   stop 'error in fixsign.f'
      nring=npb/6
      i1beg=nbd/2+1
      i1end=nbd/2+npb/2
      i2beg=i1end+nlp+nbd/2+1
      i2end=n
      
      open(16,file='tsign')

      inor=0
      do i = 1 , n
         if(i.lt.i1beg) then
            inor=inor+1
            read(16,'(i10,21x,2i5,f8.4)') 
     :         idummy,(infosign(k,inor),k=1,2),tmax0
	 elseif(i.ge.i1beg.and.i.le.i1end) then
	    read(16,'(i10)') idummy
	 elseif(i.gt.i1end.and.i.lt.i2beg) then
	    inor=inor+1
            read(16,'(i10,21x,2i5,f8.4)') 
     :         idummy,(infosign(k,inor),k=1,2),tmax0
         else
            read(16,'(i10)') idummy
         endif
      enddo
      if(idummy.ne.n.or.inor.ne.nlp+nbd) stop 'error in compsign'
      close(16)

c      do ir = 1 , inor
c         write(6,'(50i5)') infosign(1,ir),infosign(2,ir)
c      enddo

      inor=0
      do i = 1 , n
         do ibd = 1 , n
            if(ibx(ibd).eq.i) goto 105
         enddo
 105     if(bname(ibd,1).ne.'PB') then
            inor=inor+1
	    tmax=t(infosign(1,inor),i)
	    imax=int(tmax/(dabs(tmax)-0.000001))
            itemp=imax*infosign(2,inor)
            if(itemp.eq.-1) then
               do j = 1 , n
                  t(j,i)=-t(j,i)
               enddo
            endif
	 endif
      enddo

      return                                                                    
      END                                                                       

      SUBROUTINE findsign(n,ipr)  
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*4 symbol,resname,segname
      character*2 bname
      INTEGER CZ,U,UL
      common/INFO1/CZ(NAZR),U(NBSZR),UL(NAZR),LL(NAZR),
     :NELECS,NOCCA,NOCCB
      common/LBL/bname(nbszr,2),label(NBSZR,3),IBX(NBSZR) 
      common/lblaro/labar(nbszr,7)
      common/EXTRA/T(NBSZR,NBSZR),IHYB(NAZR),TOPO(NAZR,NAZR),
     :T2J(NBSZR)
      common/cring/idc6(6,20)
      common/iden/idres(nazr),symbol(nazr),resname(nazr),segname(nazr),
     :rx(nazr),ry(nazr),rz(nazr)
      dimension infosign(2,nbszr),infosign6(2,6,20)
      dimension tmax6(6)

      open(16,file='tsign')

      do i = 1 , n
         do ibd = 1 , n
            if(ibx(ibd).eq.i) goto 105
         enddo
 105     if(bname(ibd,1).ne.'PB') then
	    tmax=0.0
	    jmax=0
	    do j = 1 , n
	       if(dabs(t(j,i)).gt.tmax) then
	          tmax=dabs(t(j,i))
	          jmax=j
	       endif
	    enddo
	    infosign(1,i)=jmax
	    if(t(jmax,i).ge.0.d0) then
	       infosign(2,i)= 1
	    else
	       infosign(2,i)=-1
	    endif
	 else
	    do iat = 1 , 6
	    tmax=0.0
	    jmax=0
	    do j = ll(idc6(iat,labar(ibd,1))),ul(idc6(iat,labar(ibd,1)))
	       if(dabs(t(j,i)).gt.tmax) then
	          tmax=dabs(t(j,i))
	          jmax=j
	       endif
	    enddo
	    infosign6(1,iat,labar(ibd,1))=jmax
	    if(t(jmax,i).ge.0.d0) then
	       infosign6(2,iat,labar(ibd,1))= 1
	    else
	       infosign6(2,iat,labar(ibd,1))=-1
	    endif
	    tmax6(iat)=tmax
	    enddo
	 endif
	 if(bname(ibd,1).eq.'LP')
     :      write(16,'(i10,2x,a2,a1,2x,i5,a2,5x,2x,2i5,f8.4)') 
     :       i,bname(ibd,1),bname(ibd,2),
     :       label(ibd,1),symbol(label(ibd,1))(1:2),
     :       (infosign(k,i),k=1,2),tmax
	 if(bname(ibd,1).eq.'BD')
     :      write(16,'(i10,2x,a2,a1,2x,i5,a2,i5,a2,2i5,f8.4)') 
     :       i,bname(ibd,1),bname(ibd,2),
     :       label(ibd,1),symbol(label(ibd,1))(1:2),
     :       label(ibd,2),symbol(label(ibd,2))(1:2),
     :       (infosign(k,i),k=1,2),tmax
	 if(bname(ibd,1).eq.'PB')
     :      write(16,'(i10,2x,a2,a1,i5,2x,6(i5,a2,2i5,f8.4))') 
     :       i,bname(ibd,1),bname(ibd,2),labar(ibd,1),
     :       (idc6(ia,labar(ibd,1)),symbol(idc6(ia,labar(ibd,1)))(1:2),
     :       infosign6(1,ia,labar(ibd,1)),
     :       infosign6(2,ia,labar(ibd,1)),tmax6(ia),
     :       ia=1,6)
      enddo
      
      close(16)
      
      return                                                                    
      END