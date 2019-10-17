
      subroutine doublecheck1(iat1,iat2,natoms)
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
      dimension pre(4,4)  
      
      if(iat2.eq.0.and.symbol(iat1)(2:2).eq.'O') then
         print *
         print *,'check the orthogonality of two lone pairs of O'
         print *
         ilc=ll(iat1)
         iuc=ul(iat1)
         nhc=iuc-ilc+1
         write(6,*) ilc,iuc,nhc
         write(6,'(10i10)') ino(iat1)
         if(ino(iat1).ne.2) stop'error in checking OLPs'
         print *
         print *,'The two lone pairs of oxygen'
         do iq = 1 , 4
            do jq = 1 , ino(iat1)
               pre(jq,iq)=q(iq,ilc+jq-1)
            enddo
            write(6,'(10f10.5)')
     :           (q(iq,ilc+jq-1),jq=1,ino(iat1))
         enddo
         dotp=0.0
         do iq = 1 , 4
            dotp=dotp+pre(1,iq)*pre(2,iq)
         enddo
         write(6,*)'dot product = ',dotp
         print *
      endif
      
      return
      end

      subroutine doublecheck2(iat1,natoms)
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
      dimension pre(4,4)  
      
      if(symbol(iat1)(2:2).eq.'O') then
         print *
         print *,'check the orthogonality of O orbitals'
         print *
         ilc=ll(iat1)
         iuc=ul(iat1)
         nhc=iuc-ilc+1
         write(6,*) ilc,iuc,nhc
         write(6,'(10i10)') ino(iat1)
         if(ino(iat1).ne.4) stop'error in checking OLPs'
         print *
         print *,'The two lone pairs of oxygen'
         do iq = 1 , 4
            do jq = 1 , ino(iat1)
               pre(jq,iq)=q(iq,ilc+jq-1)
            enddo
            write(6,'(10f10.5)')
     :           (q(iq,ilc+jq-1),jq=1,ino(iat1))
         enddo
         do m = 1 , 4
            do n = 1 , m
               dotp=0.0
               do iq = 1 , 4
                  dotp=dotp+pre(m,iq)*pre(n,iq)
               enddo
               write(6,'(a1,i1,a1,i1,a2,f8.4)') '(',m,'*',n,')=',dotp
            enddo
         enddo
         print *
      endif
      
      return
      end      