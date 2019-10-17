	subroutine part1(natom,ndim,nstp,tmd,st,ht,st2,ht2)
	
c********************************************************************
c part 1: diagonalize the Hamiltonian time series step by step 
c********************************************************************	
	
	implicit real*8 (a-h,o-z)
	parameter(neigd=500,nsmp=10,ndmax=500)
	character signname*5
	real*8 homo, lumo
	real*4 tmd(nsmp),st(ndmax,ndmax,nsmp),ht(ndmax,ndmax,nsmp)
	real*4 st2(2,2,nsmp),ht2(2,2,nsmp)
	real*4 dm,pop
	logical write,dotda,doegn,dodmx,dogrn,quick,donbo,bspac
	logical lcbo,restart,dotsa,aring,sfind,scomp
        common/iounit/iin,iout,istart(40),iungf,iunhl	
	common /eigt/ dm(ndmax,ndmax,nsmp),pop(ndmax,nsmp)
        common /eigv/ omcm(neigd),zr(neigd,neigd),zi(neigd,neigd),ne,nblw
        common /inst/ ovlp(ndmax,ndmax),hami(ndmax,ndmax),
     *                dens(ndmax,ndmax)
        common /subh/ ovb(ndmax,ndmax),hmb(ndmax,ndmax)
        common/nda/ndon,nacc,idon,iacc,nb,iold(ndmax),jold(ndmax)
        common/option/restart,dotsa,nstep,unit
        dimension green(ndim,ndim),t(ndim,ndim)
	
	pi=4.d0*datan(1.d0)
	
	write(iout,*)
	write(iout,*) '*****  Part 1 called '
	write(iout,*)
	
	write=.false.
	dotda=.false.
	doegn=.true.
	dodmx=.true.
	donbo=.true.
	dogrn=.false.
	quick=.true.
	bspac=.false.
	lcbo=.true.
	scomp=.true.
	sfind=.true.
	if(restart) then
           open(24,file='brhomolumo',access='append')
        else
           open(24,file='brhomolumo',access='sequential')
        endif

	do n = 1 , nstp
	
	   write(iout,*) 'nstp =',n
	
	   do i = 1 , ndim
	      do j = 1 , ndim
	         ovlp(i,j)=st(i,j,n)
	         hami(i,j)=ht(i,j,n)
	      enddo
	   enddo
	   
c>>> calculate the eigenenergies and eigenvectors

	   if(doegn) then

	      call eigen(ndim,homo,lumo)
	      write(iunhl,'(i10,2f10.5)') n, homo, lumo      	         
	      if(nstp.gt.1) go to 103
	      write(iout,*) ' Step of Hamiltonian =', n
	      write(iout,*) ' Overlap matrix'
	      if(write) then
	         do i = 1 , ndim
	         write(iout,'(500f8.4)') (st(i,j,n),j=1,ndim)
	         enddo
	      endif
	      write(iout,*) ' Hamiltonian matrix'
	      if(write) then
	         do i = 1 , ndim
	         write(iout,'(500f8.4)') (ht(i,j,n),j=1,ndim)
	         enddo	   
	      endif
	      write(iout,*)
	      write(iout,*) ' Eigenenergies '
	      do i = 1 , ndim
	         write(iout,'(i10,f10.5)') i,omcm(i)
	         if(dabs(omcm(i)-homo).lt.1.d-6) write(iout,*)
	      enddo
	      write(iout,*)
 103          continue
 
 	   endif
 	   
c>>> calculate the density matrix
 	   
 	   if(dodmx) then
 	      call denmat(natom,ndim,n,st,ht)
c              call diagden(natom,ndim,n,st,ht)
	      if(donbo) then
	         do iden = 1 , ndim
	            do jden = 1 , ndim
	               dens(jden,iden)=dm(jden,iden,n)
	            enddo
	         enddo
	         if(n.eq.1) call ulll(natom,ndim)
 	         call subdiag(natom,ndim,t,n,aring,nring)
                 if(sfind.and.n.eq.1.and..not.restart) then
                    signname='tsign'
                    call findsign(ndim,t,signname)  
                 endif
c                 if(n.eq.4) then
c                    signname='tsig4'
c                    call findsign(ndim,t,signname)  
c                 endif
c                 if(n.eq.5) then
c                    signname='tsig5'
c                    call findsign(ndim,t,signname)  
c                 endif
                 if(scomp) then
 	            call compsign(ndim,t)
c 	            if(aring) call compsignaro(ndim,t)
                 endif
 	         call trans(ndim,t,n,st,ht,aring,nring)
              endif
           endif
           
           call subhdba(natom,ndim)
           if(bspac) call bspace(ndim,n)
	   call backscat(n,ndim,diff)
           if(n.eq.1.and.doegn.and.quick.and.(.not.restart)) then
              call faster(ndim,nstp,diff)
           endif

	   do i = 1 , ndim
	      do j = 1 , ndim
	         st(i,j,n)=ovlp(i,j)
	         ht(i,j,n)=hami(i,j)
	      enddo
	   enddo
	   
c>>> calculate Green's function or TDA
	   if(dotda) then
              ntun=1
	      energy=-10.
	      open(56,file='etun')
	      read(56,*) energy
	      close(56)
	      it=0	   
222	      it=it+1
	      call caltda(natom,ndim,energy,tda,hdd,haa)
	      if(dogrn) call grn(n,natom,ndim,energy,green)
	      energy=energy+0.001
	      if(it.lt.ntun) go to 222
	      if(ntun.gt.1) stop 'TDA as a function of E'
	      if(idon.lt.iacc) then
	         ht2(1,1,n)=hdd
	         ht2(1,2,n)=tda
	         ht2(2,1,n)=tda
	         ht2(2,2,n)=haa
	      else
	         ht2(1,1,n)=haa
	         ht2(1,2,n)=tda
	         ht2(2,1,n)=tda
	         ht2(2,2,n)=hdd
	      endif
	      st2(1,1,n)=st(idon,idon,n)
	      st2(1,2,n)=st(idon,iacc,n)
	      st2(2,1,n)=st(iacc,idon,n)
	      st2(2,2,n)=st(iacc,iacc,n)
	   endif
	   
	enddo 

	close(24)
	 
	return
	end
	
