
	subroutine subhdba(natom,ndim)
	
c********************************************************
c Partition Hamiltonian submatrices for donor, acceptor 
c and bridge
c********************************************************

	implicit real*8 (a-h,o-z)
	parameter (natmax=200,ndmax=500)
	character indt*10,ac*4,symbol*2,bname*2,segname*4,resname*3
        common /inst/ ovlp(ndmax,ndmax),hami(ndmax,ndmax),
     *                dens(ndmax,ndmax)
        common /subh/ ovb(ndmax,ndmax),hmb(ndmax,ndmax)
        common/ato/ac(natmax),symbol(40)
        common/atom2/exps(40),exps2(40),expp(40),expp2(40),expd(40),              
     *  expd2(40),expf(40),expf2(40),cs1(40),cs2(40),cp1(40),cp2(40),             
     *  cd1(40),cd2(40),cf1(40),cf2(40),couls(40),coulp(40),could(40),            
     *  coulf(40),x(natmax),y(natmax),z(natmax),ires(natmax),
     *  resname(natmax),segname(natmax)
        common/aoout/indt(1000,2),iatom(natmax),lorb(1000)
        common/iounit/iin,iout,istart(40),iungf,iunhl        
        common/nda/ndon,nacc,idon,iacc,nb,iold(ndmax),jold(ndmax)
        common/lbl/bname(ndmax,2),label(ndmax,3),ibx(ndmax)
	logical write,lcao

	write=.true.
	lcao=.false.
	
	if(.not.lcao) goto 500

	ibegd=0
	iendd=0
	ibega=0
	ienda=0
	
	do i = 1 , ndim
	   if(indt(i,1)(1:2).eq.ac(ndon)(1:2).and.lorb(i).eq.iatom(ndon))
     :        then
	      ibegd=i
	      go to 100
	   endif
	enddo
100	do i = ndim , 1 , -1
	   if(indt(i,1)(1:2).eq.ac(ndon)(1:2).and.lorb(i).eq.iatom(ndon))
     :        then
	      iendd=i
	      go to 200
	   endif
	enddo
200	do i = 1 , ndim
	   if(indt(i,1)(1:2).eq.ac(nacc)(1:2).and.lorb(i).eq.iatom(nacc))
     :	      then
	      ibega=i
	      go to 300
	   endif
	enddo
300	do i = ndim , 1 , -1
	   if(indt(i,1)(1:2).eq.ac(nacc)(1:2).and.lorb(i).eq.iatom(nacc))
     :	      then
	      ienda=i
	      go to 400
	   endif
	enddo
400	continue

	mdon=iendd-ibegd+1
	macc=ienda-ibega+1
	idon=ibegd
	iacc=ibega
	
	return
500     continue

        do i = 1 , ndim

           do ib = 1 , ndim
              if(ibx(ib).eq.i) goto 15
           enddo
15         continue

           if(bname(ib,1).eq.'BD'.and.bname(ib,2)(1:1).ne.'*'.and.(
     :       (ac(label(ib,2))(1:4).eq.' CA '.and.
     :        ac(label(ib,3))(1:4).eq.' HA ').or.
     :       (ac(label(ib,2))(1:4).eq.' CA '.and.
     :        ac(label(ib,3))(1:4).eq.' HA1'))) then
              if(resname(label(ib,2)).eq.'MET') then
                 idon=i
              endif
              if(resname(label(ib,2)).eq.'THR'.and.
     :           ires(label(ib,2)).eq.124) then
                 iacc=i
              endif
           endif

        enddo
        
	if(idon.gt.iacc) then
	   id=iacc
	   ia=idon
	else
	   id=idon
	   ia=iacc
	endif
	mb=0
	do 11 i = 1 , ndim
	   if(i.eq.id.or.i.eq.ia) go to 11
	   mb=mb+1
	   iold(mb)=i
	   nb=0
	   do 12 j = 1 , ndim
	      if(j.eq.id.or.j.eq.ia) go to 12
	      nb=nb+1
	      jold(nb)=j
12	   continue
11      continue
        if(mb.ne.nb) stop 'error mb != nb'

	if(write) then
	   write(iout,*) '   donor state =',idon
	   write(iout,*) 'acceptor state =',iacc
	endif

	return
	end		