	subroutine caltda(natom,ndim,e,tda,hdd,haa)

c***********************************************************************
c calculate the TDA
c***********************************************************************	

	implicit real*8 (a-h,o-z)
	parameter(natmax=200,ndmax=500)
	character ac*4,symbol*2,indt*10
        common/ato/ac(natmax),symbol(40)
        common/aoout/indt(1000,2),iatom(natmax),lorb(1000)
        common /inst/ ovlp(ndmax,ndmax),hami(ndmax,ndmax),
     *                dens(ndmax,ndmax)
        common /subh/ ovb(ndmax,ndmax),hmb(ndmax,ndmax)
        common/iounit/iin,iout,istart(40),iungf,iunhl     
        common/nda/ndon,nacc,idon,iacc,nb,iold(ndmax),jold(ndmax)
        dimension hmes(ndim-2,ndim-2),indx(ndim-2),hinv(ndim-2,ndim-2)
        dimension vd(ndim-2),va(ndim-2)
	logical write,test	
	
	test=.false.
	write=.false.
	
	if(nb.ne.ndim-2) stop 'nb != ndim-2'
	
	do i = 1 , ndim
	   do j = 1 , ndim
	      ovb(i,j)=0.0
	      hmb(j,i)=0.0
	   enddo
	enddo

c---> define the bridge subhamiltonian
		   
	if(idon.gt.iacc) then
	   id=iacc
	   ia=idon
	else
	   id=idon
	   ia=iacc
	endif
	do i = 1 , nb
	   do j = 1 , nb
	      ovb(i,j)=ovlp(iold(i),jold(j))
	      hmb(i,j)=hami(iold(i),jold(j))
	   enddo
	enddo
	
	if(write) then      
	   write(iout,*) ' The bridge subH is a ',nb,'*',nb,'matrix .'
	   do i = 1 , nb
	      write(iout,'(500f10.5)') (hmb(i,j),j=1,nb)
	   enddo	   	   	   
	endif

	do i = 1 , nb
	   do j = 1 , nb
	      hmes(i,j)=e*ovb(i,j)-hmb(i,j)
	      hinv(i,j)=0.d0
	   enddo
	enddo
	
	if(write) then
	   do i = 1 , nb
	      write(iout,'(500f10.5)') (hmes(i,j),j=1,nb)
	   enddo
	endif
	
	do i = 1 , nb
	   do j = 1 , nb
	      hinv(i,j)=0.d0
	   enddo
	   hinv(i,i)=1.d0
	enddo

	call ludcmp(hmes,nb,nb,indx,d)

	do j = 1 , nb
	   call lubksb(hmes,nb,nb,indx,hinv(1,j))
	enddo

	if(write) then
	   write(iout,*) 
	   write(iout,*) 'The inverse of bridge Hamiltonian matrix '
	   write(iout,*)
	   do i = 1 , nb
	      write(iout,'(500f8.3)') (hinv(i,j),j=1,nb)
	   enddo
	endif

c---> test the matrix inversion routine	
	
	if(.not.test) go to 300
	do i = 1 , nb
	   do j = 1 , nb
	      hmes(i,j)=hinv(i,j)
	   enddo
	enddo
	do i = 1 , nb
	   do j = 1 , nb
	      hinv(i,j)=0.d0
	   enddo
	   hinv(i,i)=1.d0
	enddo
	call ludcmp(hmes,nb,nb,indx,d)

	do j = 1 , nb
	   call lubksb(hmes,nb,nb,indx,hinv(1,j))
	enddo
	write(iout,*)
	write(iout,*) ' testing matrix inversion routine '
	sum=0.0
	do i = 1 , nb
	   do j = 1 , nb
	      sum=sum+dabs(e*ovb(i,j)-hmb(i,j)-hinv(i,j))
	   enddo
	enddo
	write(iout,*) ' Accumulating error = ', sum
	if(sum.gt.1.d-8) stop ' Error in LU matrix inversion ! '
300     continue


c---> reorder the coupling between D, A and HB

	do i = 1 , nb
	   if(idon.lt.iacc) then
	      vd(i)=hami(id,iold(i))-e*ovlp(id,iold(i))
	      va(i)=hami(ia,iold(i))-e*ovlp(ia,iold(i))
	   else
	      vd(i)=hami(ia,iold(i))-e*ovlp(ia,iold(i))
	      va(i)=hami(id,iold(i))-e*ovlp(id,iold(i))
	   endif
	enddo

	tda=hami(idon,iacc)
	hdd=hami(idon,idon)
	haa=hami(iacc,iacc)
	
c	write(*,'(10f8.3)')(vd(i),i=1,nb)
c	write(*,'(10f8.3)')(va(i),i=1,nb)
	
	do i = 1 , nb
	   do j = 1 , nb
	      tda=tda+vd(i)*hinv(i,j)*va(j)
	      hdd=hdd+vd(i)*hinv(i,j)*vd(j)
	      haa=haa+va(i)*hinv(i,j)*va(j)
	   enddo
	enddo

     	write(iout,'(a5,1x,f10.5,1x,a4,1x,f20.10)')' TDA(',e,') = ',tda
	write(iout,*)
	write(iout,'(a)') ' Two-State Approximation'
	write(iout,'(2f20.10)') hdd,tda
	write(iout,'(2f20.10)') tda,haa
	write(iungf,'(4f20.10)') e,tda,hdd,haa


	return
	end