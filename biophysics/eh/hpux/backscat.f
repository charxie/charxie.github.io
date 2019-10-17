	subroutine backscat(n,ndim,diff)
	
c********************************************************************
c find the two eigenstates responsible for tunneling
c********************************************************************	
	
	implicit real*8 (a-h,o-z)
	parameter(neigd=500,nsmp=10,ndmax=500)
	real*8 homo,lumo
	logical restart,dotsa
        common/iounit/iin,iout,istart(40),iungf,iunhl	
        common /eigv/ omcm(neigd),zr(neigd,neigd),zi(neigd,neigd),ne,nblw
        common /inst/ ovlp(ndmax,ndmax),hami(ndmax,ndmax),
     *                dens(ndmax,ndmax)
        common/nda/ndon,nacc,idon,iacc,nb,iold(ndmax),jold(ndmax)
        common/option/restart,dotsa,nstep,unit
        dimension popda(ndim)
	logical write,lift
	
	write=.true.
	lift=.true.
	open(78,file='donacp')
	write(iout,*) idon,iacc

c lift to the gap
        if(.not.lift) goto 105
c	eta1=4.0
c	eta2=4.0
c	hami(idon,idon)=-10.0
c	hami(iacc,iacc)=-10.04205
        read(78,*) hami(idon,idon),hami(iacc,iacc)
        read(78,*) eta1,eta2
        close(78)

        do i = 1 , ndim
           if(i.ne.idon.and.i.ne.iacc) then
              if(dabs(hami(idon,i)).gt.0.0) then
                 hami(idon,i)=hami(idon,i)/eta1
                 hami(i,idon)=hami(i,idon)/eta1
              endif
              if(dabs(ovlp(idon,i)).gt.0.0) then
                 ovlp(idon,i)=ovlp(idon,i)/eta2
                 ovlp(i,idon)=ovlp(i,idon)/eta2
              endif
	      if(dabs(hami(iacc,i)).gt.0.0) then
	         hami(iacc,i)=hami(iacc,i)/eta1
	         hami(i,iacc)=hami(i,iacc)/eta1
	      endif
	      if(dabs(ovlp(iacc,i)).gt.0.0) then
	         ovlp(iacc,i)=ovlp(iacc,i)/eta2
	         ovlp(i,iacc)=ovlp(i,iacc)/eta2
	      endif
	   endif
	enddo
 105    continue
        if(restart.or.n.gt.1) return

c>>> calculate the eigenenergies and eigenvectors

	call eigen(ndim,homo,lumo)
	
	write(iout,*)'eigenstates with donor and acceptor'
	write(iout,'(i10,f10.5)')(i,omcm(i),i=1,ndim)
	write(iout,*)
	      
	pmax1=0.0
	pmax2=0.0
	do i = 1 , ndim
	   popda(i)=zr(idon,i)**2+zr(iacc,i)**2
	   if(popda(i).gt.pmax1) then
	      pmax1=popda(i)
	      imax1=i
	   endif
	enddo	      
	do i = 1 , ndim
	   if(popda(i).gt.pmax2.and.popda(i).lt.pmax1) then
	      pmax2=popda(i)
	      imax2=i
	   endif
	enddo
	write(*,*) imax1,pmax1,imax2,pmax2
	if(imax2.lt.imax1) then
	   itemp=imax1
	   ptemp=pmax1
	   imax1=imax2
	   pmax1=pmax2
	   imax2=itemp
	   pmax2=ptemp
	endif
        if(imax2-imax1.ne.1) write(*,*) 'eigenstates problem'
        write(iout,*)'eigenvector(',omcm(imax1),')',imax1,pmax1
        write(iout,'(10f8.4)')(zr(j,imax1),j=1,ndim)
        write(iout,*)'eigenvector(',omcm(imax2),')',imax2,pmax2
        write(iout,'(10f8.4)')(zr(j,imax2),j=1,ndim)
        
        diff=dabs(omcm(imax1)-omcm(imax2))
        write(iout,*)'diff=',diff
        
	return
	end
	
