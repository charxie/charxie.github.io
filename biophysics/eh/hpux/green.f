	subroutine grn(istep,natom,ndim,e,green)

c***********************************************************************
c calculate the Green function
c note: The array dimension of LU subroutines MUST be the same !
c***********************************************************************	

	implicit real*8 (a-h,o-z)
	parameter(natmax=200,ndmax=500)
	character ac*4,indt*10,symbol*2,segname*4,resname*3
	character donor*8,acceptor*8
        common /inst/ ovlp(ndmax,ndmax),hami(ndmax,ndmax),
     *                dens(ndmax,ndmax)
        common/iounit/iin,iout,istart(40),iungf,iunhl 
        common/ato/ac(natmax),symbol(40)  
        common/aoout/indt(1000,2),iatom(natmax),lorb(1000) 
        common/atom2/exps(40),exps2(40),expp(40),expp2(40),expd(40),
     *  expd2(40),expf(40),expf2(40),cs1(40),cs2(40),cp1(40),cp2(40),
     *  cd1(40),cd2(40),cf1(40),cf2(40),couls(40),coulp(40),could(40),
     *  coulf(40),x(natmax),y(natmax),z(natmax),ires(natmax),
     *  resname(natmax),segname(natmax)
        dimension hmes(ndim,ndim),indx(ndim),hinv(ndim,ndim)
        dimension temp(ndim,ndim)
        dimension green(ndim,ndim)
	logical write
	
	write=.false.
	donor(1:4)=' CA '
	acceptor(1:4)=' CA '
	donor(5:8)='ALAP'
	acceptor(5:8)='ALAP'
	iresd=1
	iresa=8

c>>> select donor and acceptor for the Green's Function

	id=0
	ia=0
	do ie = 1 , natom
	   if(segname(ie).eq.donor(5:8).and.
     :        ac(ie)(1:4).eq.donor(1:4).and.
     :        ires(ie).eq.iresd) 
     :        id=ie
	   if(segname(ie).eq.acceptor(5:8).and.
     :        ac(ie)(1:4).eq.acceptor(1:4).and.
     :        ires(ie).eq.iresa) 
     :        ia=ie	
	enddo	

c	write(*,'(2i5,2a4)') id,ia,segname(id),segname(ia)
	if(id*ia.eq.0) stop 'Warning: No donor or acceptor found !'
	
	do i = 1 , ndim
	   if(indt(i,1)(1:2).eq.ac(id)(1:2).and.
     *        lorb(i).eq.iatom(id)) then
	      ibegd=i
	      go to 100
	   endif
	enddo
100	do i = ndim , 1 , -1
	   if(indt(i,1)(1:2).eq.ac(id)(1:2).and.
     *        lorb(i).eq.iatom(id)) then
	      iendd=i
	      go to 200
	   endif
	enddo
200	do i = 1 , ndim
	   if(indt(i,1)(1:2).eq.ac(ia)(1:2).and.
     *        lorb(i).eq.iatom(ia)) then
	      ibega=i
	      go to 300
	   endif
	enddo
300	do i = ndim , 1 , -1
	   if(indt(i,1)(1:2).eq.ac(ia)(1:2).and.
     *        lorb(i).eq.iatom(ia)) then
	      ienda=i
	      go to 400
	   endif
	enddo
400	continue

c	write(*,*) ibegd,iendd,ibega,ienda

	n1=ibegd
	n2=ibega

	do i = 1 , ndim
	   do j = 1 , ndim
	      hmes(i,j)=e*ovlp(i,j)-hami(i,j)
	      hinv(i,j)=0.d0
	   enddo
	enddo
	
	do i = 1 , ndim
	   do j = 1 , ndim
	      hinv(i,j)=0.d0
	   enddo
	   hinv(i,i)=1.d0
	enddo

	call ludcmp(hmes,ndim,ndim,indx,d)

	do j = 1 , ndim
	   call lubksb(hmes,ndim,ndim,indx,hinv(1,j))
	enddo

	if(write) then
	   write(iout,*) 
	   write(iout,*) ' The inverse of Hamiltonian matrix '
	   write(iout,*)
	   do i = 1 , ndim
	      write(iout,'(500d10.3)') (hinv(i,j),j=1,ndim)
	   enddo
	endif
	
	call mprod(ndim,ovlp,hinv,temp)
	call mprod(ndim,temp,ovlp,green)
	
	write(iout,1000) e,n1,n2,dabs(green(n1,n2))
 1000   format(1x,'E=',f8.4,3x,'G(',i3,',',i3,')=',3e20.10)
 	write(iungf,'(f10.5,e20.10)') real(istep),dabs(green(n1,n2))	

	return
	end

