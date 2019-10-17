C  MODIFY ORTHOGONAL MATRIX T FOR TRANSFORMATION FROM AO'S TO               
C  NATURAL HYBRID BOND ORBITALS, USING INPUT FOCK MATRIX 'FM'. 
c  *** for C=O double bond only
c  *** in input file ipr must be set > 3 as to get dm          

      SUBROUTINE cobond(FM,FMAO,DM,T,N,ipr)  
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      parameter(neigd=500)
      real*8 kappa
      character*2 bname
      integer CZ,U,UL
      common/INFO1/CZ(NAZR),U(NBSZR),UL(NAZR),LL(NAZR),NELECS,NOCCA,NOCCB
      common/LBL/bname(nbszr,2),label(NBSZR,3),IBX(NBSZR) 
      common/ARRAYS/A(NBSZR,NBSZR,3)
      common/EXTRA/TNEW(NBSZR,NBSZR),IHYB(NAZR),
     :             TOPO(NAZR,NAZR),T2J(NBSZR)
      common/eigv/omcm(neigd),zr(neigd,neigd),zi(neigd,neigd),nex,nblw
      character*4 symbol,resname,segname
      common/iden/idres(nazr),symbol(nazr),resname(nazr),segname(nazr),
     :rx(nazr),ry(nazr),rz(nazr)
      dimension FM(NBSZR,NBSZR),DM(NBSZR,NBSZR),T(NBSZR,NBSZR)
      dimension FMAO(NBSZR,NBSZR)
      dimension sm(n,n),hm(n,n)
      dimension index_cobond(nbszr),ind1_cobond(nbszr)
      dimension h(2,2),v(2,2),d(2,2),d1(2,2)
      DATA LLP,LBD,L3C/'LP','BD','3C'/

      print *
      print *, 'C=O bond treatment'
      print *, 'total number of bonds and antibonds = ', n

      do i = 1 , n
         do j = 1 , n
            tnew(i,j)=t(i,j)
         enddo
      enddo
     
c      goto 2300
      
      if(ipr.gt.3) call anlyze(t,nbszr)

      num_cobond=0
      
      do ibd = 1 , n
      
         if((symbol(label(ibd,2))(2:2).eq.'C'
     : .and. symbol(label(ibd,3))(2:2).eq.'O')
     : .or. (symbol(label(ibd,2))(2:2).eq.'O'
     : .and. symbol(label(ibd,3))(2:2).eq.'C')
     :     ) then
     
            if(symbol(label(ibd,2))(2:3).eq.'OG'.or.
     :         symbol(label(ibd,3))(2:3).eq.'OG'.or.
     :         symbol(label(ibd,2))(2:3).eq.'OH'.or.
     :         symbol(label(ibd,3))(2:3).eq.'OH'.or.     
     :         symbol(label(ibd,2))(2:3).eq.'OD'.or.
     :         symbol(label(ibd,3))(2:3).eq.'OD'
     :        ) goto 2020            
     
            num_cobond=num_cobond+1
            index_cobond(num_cobond)=ibx(ibd)
            ind1_cobond(num_cobond)=ibd
            
 2020       continue
    
         endif

      enddo
      
      if(ipr.gt.3) then
         print *,'the molecule contains ',num_cobond,' C=O (anti)bonds'
         write(6,'(100i4)') (index_cobond(i),i=1,num_cobond)
         print *
      endif
      
      do i = 1 , num_cobond, 2
      
         if(label(ind1_cobond(i),2).ne.label(ind1_cobond(i+1),2)
     :  .or.label(ind1_cobond(i),3).ne.label(ind1_cobond(i+1),3)) 
     :   stop 'wrong C=O bond found!'
      
	 h(1,1)=fm(index_cobond(i  ),index_cobond(i  ))
	 h(1,2)=fm(index_cobond(i  ),index_cobond(i+1))
	 h(2,1)=fm(index_cobond(i+1),index_cobond(i  ))
	 h(2,2)=fm(index_cobond(i+1),index_cobond(i+1))
	 d(1,1)=dm(index_cobond(i  ),index_cobond(i  ))
	 d(1,2)=dm(index_cobond(i  ),index_cobond(i+1))
	 d(2,1)=dm(index_cobond(i+1),index_cobond(i  ))
	 d(2,2)=dm(index_cobond(i+1),index_cobond(i+1))
	 
	 if(ipr.gt.3) then
	    print *,'print the subblock dm for the',int(i/2)+1,
     :           '-th C=O (anti)bond'
 	    write(6,'(2f10.5)') ((d(m,k),k=1,2),m=1,2)
         endif
	 
	 delta=h(1,1)-h(2,2)
	 hplus=h(1,1)+h(2,2)
	 kappa=dsqrt(delta*delta+4.d0*h(1,2)*h(2,1))
c	 print *
c	 write(6,'(2f10.5)') 0.5*(hplus+kappa),0.5*(hplus-kappa)
	 if(h(1,2).gt.0) then
	    v(1,1)= dsqrt((kappa+delta)/2.0/kappa)
	    v(1,2)= dsqrt((kappa-delta)/2.0/kappa)
	    v(2,1)= dsqrt((kappa-delta)/2.0/kappa)
	    v(2,2)=-dsqrt((kappa+delta)/2.0/kappa)
	 else
	    v(1,1)=-dsqrt((kappa+delta)/2.0/kappa)
	    v(1,2)= dsqrt((kappa-delta)/2.0/kappa)
	    v(2,1)=-dsqrt((kappa-delta)/2.0/kappa)
	    v(2,2)=-dsqrt((kappa+delta)/2.0/kappa)
	 endif
c	 print *
c	 write(6,'(2f10.5)') ((h(m,k),k=1,2),m=1,2)
c	 write(6,'(2f10.5)') ((v(m,k),k=1,2),m=1,2)
	 
	 d1(1,1)=v(1,1)*v(1,1)*d(1,1)+v(1,1)*v(1,2)*d(1,2)
     :          +v(1,2)*v(1,1)*d(2,1)+v(1,2)*v(1,2)*d(2,2)	 
	 d1(1,2)=v(1,1)*v(2,1)*d(1,1)+v(1,1)*v(2,2)*d(1,2)
     :          +v(1,2)*v(2,1)*d(2,1)+v(1,2)*v(2,2)*d(2,2)	 
	 d1(2,1)=v(2,1)*v(1,1)*d(1,1)+v(2,1)*v(1,2)*d(1,2)
     :          +v(2,2)*v(1,1)*d(2,1)+v(2,2)*v(1,2)*d(2,2)	 
	 d1(2,2)=v(2,1)*v(2,1)*d(1,1)+v(2,1)*v(2,2)*d(1,2)
     :          +v(2,2)*v(2,1)*d(2,1)+v(2,2)*v(2,2)*d(2,2)
	 
	 if(ipr.gt.3) print *
	 if(ipr.gt.3) write(6,'(2f10.5)') ((d1(m,k),k=1,2),m=1,2)
	 
	 do k = 1 , n
	    tnew(k,index_cobond(i  ))=v(1,1)*t(k,index_cobond(i  ))
     :                               +v(1,2)*t(k,index_cobond(i+1))
	    tnew(k,index_cobond(i+1))=v(2,1)*t(k,index_cobond(i  ))
     :                               +v(2,2)*t(k,index_cobond(i+1))
	 enddo
	 
	 
      enddo

      call lonepair(fm,dm,n,ipr)
c      call lonepair2(fm,dm,n,ipr)
c      call lonepair1(n,ipr)

      call findsign(n,ipr)
c      call compsign(n,ipr)      

 2300 call anlyze(tnew,nbszr)
c      call anlyze1(tnew,nbszr,5)
c      call anlyze1(tnew,nbszr,6)
c      call anlyze1(tnew,nbszr,58)
c      call anlyze1(tnew,nbszr,59)
     
      do i = 1 , n
         do j = 1 , n
            a(i,j,1)=fmao(i,j)
         enddo
      enddo
      call trans(n,1)
c      do i = 1 , n
c         do j = 1 , i
c            a(j,i,1)=a(i,j,1)
c         enddo
c      enddo
      if(ipr.gt.3) then
         do i = 1 , n
            write(6,'(i4,500f8.4)')i,(a(j,i,1),j=1,n)
         enddo
      endif
      
c    Print out Fock matrix 

c      open(file='fock1', unit=11)
c      open(file='fock2', unit=12)
      open(file='fock' , unit=13)

c      do i = 1 , n
c         if(label(i,2).eq.'*') then
c            ibeg=i
c            goto 5890
c         endif
c      enddo
c 5890 continue
c      Do I=ibeg,N
c         write(11,'(500f10.4)')(27.2*A(I,J,1),J=ibeg,N)
c      End Do
c      close(unit=11)
c      Do I=1,ibeg-1
c         write(12,'(500f10.4)')(27.2*A(I,J,1),J=1,ibeg-1)
c      End Do
c      close(unit=12)
      Do I=1,N
         write(13,'(500d10.4)')(27.2*A(I,J,1),J=1,N)
      End Do
      close(unit=13)
      
c      open(file='fij_lp',unit=14,access='append')
c      write(14,'(100f10.4)') 27.2*a(135,135,1),27.2*a(48,69,1),
c     :                       27.2*a(68,68,1),27.2*a(69,69,1),
c     :                       27.2*a(47,47,1),27.2*a(48,48,1),        
c     :                       27.2*a(68,69,1),27.2*a(47,48,1)
c      close(14)
c      open(file='fij_co',unit=14,access='append')
c      write(14,'(100f10.4)') 27.2*a(18,58,1),27.2*a(19,58,1),
c     :                       27.2*a(18,18,1),27.2*a(19,19,1),
c     :                       27.2*a(58,58,1),27.2*a(59,59,1),
c     :                       27.2*a(58,59,1),27.2*a(18,19,1)
c      close(14)
c      open(file='fij_n1',unit=14,access='append')
c      write(14,'(100f10.4)') 27.2*a(62,62,1),27.2*a(93,93,1),
c     :                       27.2*a(94,94,1),27.2*a(23,23,1),
c     :                       27.2*a(24,24,1),27.2*a(62,23,1),
c     :                       27.2*a(62,24,1),27.2*a(62,93,1),
c     :                       27.2*a(62,94,1)
c      close(14)
     
      if(ipr.le.3) return
      
      print *
      print *,'check the transformation by re-diagonalizing fock matrix'
      print *
      do i = 1 , n
         do j = 1 , i
            hm(i,j)=a(i,j,1)
            hm(j,i)=a(j,i,1)
            sm(i,j)=0.0
            sm(j,i)=0.0
         enddo
         sm(i,i)=1.0
      enddo
      
      call eigen(n,hm,sm)
      write(6,'(5f8.4)')(omcm(i),i=1,n)
      
      print *
      print *
      print *

      return                                                                    
      END