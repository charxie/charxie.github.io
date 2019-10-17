      subroutine fourth(ibd,iring)

c********************************************************************
c build the fourth vector perpendicular to the sp2 hydrids
c********************************************************************

      implicit real*8 (a-h,o-z)
      parameter(natmax=200,ndmax=500)
      integer ul
      logical write
      character segname*4,resname*3,ac*4,symbol*2
      common/ato/ac(natmax),symbol(40)
      common/atom2/exps(40),exps2(40),expp(40),expp2(40),expd(40),              
     *expd2(40),expf(40),expf2(40),cs1(40),cs2(40),cp1(40),cp2(40),             
     *cd1(40),cd2(40),cf1(40),cf2(40),couls(40),coulp(40),could(40),            
     *coulf(40),x(natmax),y(natmax),z(natmax),ires(natmax),
     *resname(natmax),segname(natmax)
      common/info1/ul(natmax),ll(natmax),nele,nocc
      common/q/q(4,ndmax)
      common/pol/pol(ndmax,3),iathy(ndmax,3),ino(natmax)
      common/iounit/iin,iout,istart(40),iungf,iunhl
      common/cring/idc6(6,10),ovar(6,6,10)
      dimension pre(4,4),abc(3,3),d(3),cba(3,3),indx(3),p(3)

      write=.false.
  
      write(iout,*)
      write(iout,*)
     :   'Build Pz orbitals from three sp2 hybrids of each carbon'
      write(iout,*)

      do 125 iat = 1 , 6
         ilc=ll(idc6(iat,iring))
         iuc=ul(idc6(iat,iring))
         nhc=iuc-ilc+1
      if(write) write(iout,*) ilc,iuc,nhc
      if(write) write(iout,'(10i10)') ino(idc6(iat,iring))
      if(ino(idc6(iat,iring)).ne.3) 
     :                stop 'error in building the 4th orbitals'
      if(write) write(iout,*)
      if(write) write(iout,*)'The original three sp2 hybrids'
      do iq = 1 , 4
         do jq = 1 , ino(idc6(iat,iring))
            pre(jq,iq)=q(iq,ilc+jq-1)
         enddo
         if(write) write(iout,'(10f10.5)')
     :      (q(iq,ilc+jq-1),jq=1,ino(idc6(iat,iring)))
      enddo
      do i = 1 , 3
         do j = 1 , 3
            abc(j,i)=pre(j,i)
            cba(j,i)=0.0
         enddo
         cba(i,i)=1.0
      enddo
      do i = 1 , 3
         d(i)=pre(i,4)
      enddo
      call ludcmp(abc,3,3,indx,ddd)
      do i = 1 , 3
	 call lubksb(abc,3,3,indx,cba(1,i))
      enddo
      do i = 1 , 3
         p(i)=-d(1)*cba(i,1)-d(2)*cba(i,2)-d(3)*cba(i,3)
      enddo
      deno=dsqrt(p(1)*p(1)+p(2)*p(2)+p(3)*p(3)+1.0)
      if(p(2).lt.0.0) deno=-deno
      do i = 1 , 3
         pre(4,i)=p(i)/deno
      enddo      
      pre(4,4)=1.0/deno
      if(write) then
         write(iout,*)
         write(iout,*)
     *     'pz orbitals built from sp2 hybrids (the 4th vector)'
         do j = 1 , 4
            write(iout,'(4f10.5)') (pre(i,j),i=1,4)
         enddo
         write(iout,*)
         write(iout,*)'Check their orthogonalities'
      endif
      sum=0.0
      do i = 1 , 4
         do j = i , 4
            dotp=0.0
            do k = 1 , 4
               dotp=dotp+pre(i,k)*pre(j,k)
            enddo
            if(write) 
     :      write(iout,'(1x,a1,i1,a1,i1,a1,f10.5)')'(',i,'*',j,')',dotp
            sum=sum+0.0
         enddo
      enddo
      if(sum.gt.5.0) stop 'Problem with orthogonalities'
      ino(idc6(iat,iring))=ino(idc6(iat,iring))+1  
c      cmax=0.0
c      do iq = 1 , 4
c         if(dabs(pre(ino(idc6(iat,iring)),iq)).gt.cmax) then
c            cmax=dabs(pre(ino(idc6(iat,iring)),iq))
c            imax=iq
c         endif
c      enddo
      imax=2
      if(pre(ino(idc6(iat,iring)),imax).lt.0.0) then
         do iq = 1 , 4
            pre(ino(idc6(iat,iring)),iq)=-pre(ino(idc6(iat,iring)),iq)
         enddo
      endif

      do iq = 1 , 4
         do jq = 1 , ino(idc6(iat,iring))
            q(iq,ilc+jq-1)=pre(jq,iq)
         enddo
      enddo
125   continue
      ibd=ibd+3
      
      return
      end
      

      subroutine getovar(ndim,nring)

c********************************************************************
c get overlap subblock corresponding to aromatic rings
c********************************************************************
      
      implicit real*8 (a-h,o-z)
      parameter(natmax=200,ndmax=500)
      integer ul
      character bname*2
      common/inst/ovlp(ndmax,ndmax),hami(ndmax,ndmax),dens(ndmax,ndmax)
      common/lbl/bname(ndmax,2),label(ndmax,3),ibx(ndmax)
      common/info1/ul(natmax),ll(natmax),nele,nocc
      common/q/q(4,ndmax)
      common/pol/pol(ndmax,3),iathy(ndmax,3),ino(natmax)
      common/iounit/iin,iout,istart(40),iungf,iunhl
      common/cring/idc6(6,10),ovar(6,6,10)
      dimension index(24),arao(24,24),coef(6,24)

      do iring = 1 , nring

         incr=0
         do iat = 1 , 6
            kl=ll(idc6(iat,iring))
            ku=ul(idc6(iat,iring))
            do k = kl , ku
               incr=incr+1
               index(incr)=k
            enddo
         enddo
         if(incr.ne.24) stop 'error in getovar'
         do iat = 1 , 6
            kl=ll(idc6(iat,iring))
            ku=ul(idc6(iat,iring))
            ieff=0
            do iorb = 1 , 24
               if(index(iorb).ge.kl.and.index(iorb).le.ku) then
                  ieff=ieff+1
                  coef(iat,iorb)=q(ieff,ku)
               else
                  coef(iat,iorb)=0.0
               endif
            enddo
         enddo
         do i = 1 , 24
            do j = 1 , 24
               arao(i,j)=ovlp(index(i),index(j))
            enddo
         enddo
         do i = 1 , 6
            do j = 1 , 6
               ovar(i,j,iring)=0.0
               do k = 1 , 24
                  do l = 1 , 24
                     ovar(i,j,iring)=ovar(i,j,iring)+
     :                    arao(k,l)*coef(i,k)*coef(j,l)
                  enddo
               enddo
            enddo
         enddo
c         write(iout,*)
c         write(iout,'(24f10.5)') ((arao(i,j),i=1,24),j=1,24)
c         write(iout,*)
c         write(iout,'(6f10.5)') ((ovar(i,j,iring),j=1,6),i=1,6)

      enddo
      
      return
      end
      
      subroutine reformaro(t,ndim,nring)

c********************************************************************
c set up the T matrix part for aromatic rings
c********************************************************************
      
      implicit real*8 (a-h,o-z)
      parameter(natmax=200,ndmax=500)
      integer ul
      character bname*2
      common/inst/ovlp(ndmax,ndmax),hami(ndmax,ndmax),dens(ndmax,ndmax)
      common/lbl/bname(ndmax,2),label(ndmax,3),ibx(ndmax)
      common/info1/ul(natmax),ll(natmax),nele,nocc
      common/q/q(4,ndmax)
      common/pol/pol(ndmax,3),iathy(ndmax,3),ino(natmax)
      common/iounit/iin,iout,istart(40),iungf,iunhl
      common/cring/idc6(6,10),ovar(6,6,10)
      common/homo/lumo
      dimension t(ndim,ndim)
      dimension anorm(6),cop(6,6),crs(6,6)
      data cop/1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
     :         2.0, 1.0, 1.0,-1.0,-1.0,-2.0,
     :         0.0,-1.0, 1.0,-1.0, 1.0, 0.0,
     :         0.0, 1.0,-1.0,-1.0, 1.0, 0.0,
     :         2.0,-1.0,-1.0,-1.0,-1.0, 2.0,
     :         1.0,-1.0,-1.0, 1.0, 1.0,-1.0/
      anorm(1)=1.0d0/dsqrt(6.d0)
      anorm(2)=0.5d0/dsqrt(3.d0)
      anorm(3)=0.5d0
      anorm(4)=anorm(3)
      anorm(5)=anorm(2)
      anorm(6)=anorm(1)

c reorder occupied BO's to put lonepairs last and aromatic pi bonds 
c right before lonepairs

      do nlp = 1 , nocc
         if(bname(nlp,1).ne.'LP') goto 20
      enddo
   20 nlp=nlp-1
      npb=0
      do npi = 1 , nocc
         if(bname(npi,1).eq.'PB') npb=npb+1
      enddo
      nbd=nocc-nlp-npb
      
      do ibd = 1 , ndim
         if(ibd.le.nlp) then
            ibx(ibd)=ibd+nbd+npb
         elseif(ibd.gt.nlp.and.ibd.le.nocc) then
            ibx(ibd)=ibd-nlp
         else
            ibx(ibd)=ibd
         endif         
      enddo
      
      lumo=nlp+nbd+npb+1
      
c>>> build the part of the transformation matrix correponding to the PI bonds

      do iring = 1 , nring

      jet=0
      do ibd = nlp+nbd+3*iring-2, nlp+nbd+3*iring
         jet=jet+1
         do iat = 1 , 6
            crs(iat,jet)=cop(iat,jet)*anorm(jet)
         enddo
      enddo
      do ibd = nocc+nbd+3*iring-2, nocc+nbd+3*iring
         jet=jet+1
         do iat = 1 , 6
            crs(iat,jet)=cop(iat,jet)*anorm(jet)
         enddo
      enddo
      
      jet=0
      do ibd = nlp+nbd+3*iring-2, nlp+nbd+3*iring
         jet=jet+1
         do iat = 1 , 6
            jl=ll(idc6(iat,iring))
            ju=ul(idc6(iat,iring))
            irow=0
            icol=jl+ino(idc6(iat,iring))-1
            do j=jl,ju
               irow=irow+1
               t(j,ibx(ibd))=q(irow,icol)*crs(iat,jet)
            enddo
         enddo
      enddo
      do ibd = nocc+nbd+3*iring-2, nocc+nbd+3*iring
         jet=jet+1
         do iat = 1 , 6
            jl=ll(idc6(iat,iring))
            ju=ul(idc6(iat,iring))
            irow=0
            icol=jl+ino(idc6(iat,iring))-1
            do j=jl,ju
               irow=irow+1
               t(j,ibx(ibd))=q(irow,icol)*crs(iat,jet)
            enddo
         enddo
      enddo
      
      
      enddo
      
      return
      end      