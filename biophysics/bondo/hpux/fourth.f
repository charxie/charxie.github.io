c build the fourth vector perpendicular to the sp2 hydrids

      subroutine fourth(ibd,iring,natoms,ipr)
      include 'PARAM.include'
      implicit double precision (a-h,o-z)
      INTEGER CZ,U,UL 
      character*4 symbol,resname,segname
      integer idres 
      common/iden/idres(nazr),symbol(nazr),resname(nazr),segname(nazr),
     :rx(nazr),ry(nazr),rz(nazr)
      COMMON/INFO1/CZ(NAZR),U(NBSZR),UL(NAZR),LL(NAZR),
     :NELECS,NOCCA,NOCCB
      COMMON/Q/Q(4,NBSZR) 
      COMMON/POL/POL(NBSZR,3),IATHY(NBSZR,3),INO(NAZR)
      common/cring/idc6(6,20)
      dimension pre(4,4),abc(3,3),d(3),cba(3,3),indx(3),p(3)
  
      print *
      print *, 'Build Pz orbitals from three sp2 hybrids of each carbon'
      print *

      do 125 iat = 1 , 6
      ilc=ll(idc6(iat,iring))
      iuc=ul(idc6(iat,iring))
      nhc=iuc-ilc+1
      if(ipr.gt.3) write(6,*) ilc,iuc,nhc
      if(ipr.gt.3) write(6,'(10i10)') ino(idc6(iat,iring))
c      if(ino(idc6(iat,iring)).ne.3) 
c     :                stop 'error in building the 4th orbitals'
      if(ipr.gt.3) print *
      if(ipr.gt.3) print *,'The original three sp2 hybrids'
      do iq = 1 , 4
         do jq = 1 , ino(idc6(iat,iring))
            pre(jq,iq)=q(iq,ilc+jq-1)
         enddo
         if(ipr.gt.3) write(6,'(10f10.5)')
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
      if(ipr.gt.3) then
         print *
         print *,'pz orbitals built from sp2 hybrids (the 4th vector)'
         do j = 1 , 4
            write(6,'(4f10.5)') (pre(i,j),i=1,4)
         enddo
         print *
         print *, 'Check their orthogonalities'
      endif
      sum=0.0
      do i = 1 , 4
         do j = i , 4
            dotp=0.0
            do k = 1 , 4
               dotp=dotp+pre(i,k)*pre(j,k)
            enddo
            if(ipr.gt.3) 
     :      write(6,'(1x,a1,i1,a1,i1,a1,f10.5)')'(',i,'*',j,')',dotp
            sum=sum+0.0
         enddo
      enddo
      if(sum.gt.5.0) stop 'Problem with orthogonalities'
      ino(idc6(iat,iring))=ino(idc6(iat,iring))+1  
      do iq = 1 , 4
         do jq = 1 , ino(idc6(iat,iring))
            q(iq,ilc+jq-1)=pre(jq,iq)
         enddo
      enddo
125   continue
      ibd=ibd+3
      
      return
      end