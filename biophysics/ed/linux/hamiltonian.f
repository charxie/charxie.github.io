c<<<read time-dependent Fock and overlap matrix>>>

	subroutine hamiltonian(tau,nstp,ndim,tmd,h,cof,protein,dotsa,eta)

	implicit real*8 (a-h,o-z)
	parameter (ncut=1)
	character filename1*15,protein*4,filename2*10
	logical static,restart,freeze,fixdna
	logical dotsa
	real*4 tmd(nstp),h(ndim,ndim,nstp),cof(ndim,ndim,nstp,3)
	real*4 hzero(ndim,ndim)
	dimension have(ndim,ndim),hrms(ndim,ndim)
	real*4 hamf(ndim,ndim,0:ncut)
        common/ida/time0,edon,eacc,idon,iacc,i0,j0,restart,nstep
        common /io/ iin,iout
	
	do i = 1 , ndim
	   do j = 1 , ndim
	      hzero(j,i)=0.0
	      have(j,i)=0.d0
	      hrms(j,i)=0.d0
	   enddo
	enddo
	
	do n = 0, ncut
	   do i = 1 , ndim
	      do j = 1 , ndim
	         hamf(j,i,n)=0.0
	      enddo
	   enddo
	enddo
	
	static=.false.
	freeze=.false.
	fixdna=.false.
	imat=4
	jmat=6
	unit=1.51653*1.0
	istart=1
	if(.not.restart) then
	   nstp0=1
	else
	   call getsteps(nstp0,dotsa)
	   nstp0=nstp0-1
	endif
	
	if(.not.restart) then
	   period=unit*(nstp-1)
	else
	   period=unit*(nstp-2)
	endif
	tau=period/nstep
	twopi=8.d0*datan(1.d0)
	comega=0.5*twopi/period
	deltat=unit
	deltao=comega
	nstp1=nstp0+nstp-1
        do i = nstp0 , nstp1
           tmd(i-nstp0+1)=unit*(i-1)
        enddo

	filename1(1:4)='hami'
	filename1(5:10)='/fock.'

	do n = nstp0, nstp1
	   neff=n-nstp0+1
	   if(n.lt.10) then
	      filename1(11:14)='0000'
	      filename1(15:15)=char(48+n)
	   elseif(n.lt.100) then
	      nn1=mod(n,10)
	      nn2=int(n/10)
	      filename1(11:13)='000'
	      filename1(14:14)=char(48+nn2)
	      filename1(15:15)=char(48+nn1)
	   elseif(n.lt.1000) then
	      nn1=int(n/100)
	      nnx=mod(n,100)
	      nn2=int(nnx/10)
	      nn3=mod(nnx,10)
	      filename1(11:12)='00'
	      filename1(13:13)=char(48+nn1)
	      filename1(14:14)=char(48+nn2)
	      filename1(15:15)=char(48+nn3)
	   elseif(n.lt.10000) then
	      nn1=int(n/1000)
	      nnx=mod(n,1000)
	      nn2=int(nnx/100)
	      nny=mod(nnx,100)
	      nn3=int(nny/10)
	      nn4=mod(nny,10)
	      filename1(11:11)='0'
	      filename1(12:12)=char(48+nn1)
	      filename1(13:13)=char(48+nn2)
	      filename1(14:14)=char(48+nn3)
	      filename1(15:15)=char(48+nn4)
	   elseif(n.lt.100000) then
	      nn1=int(n/10000)
	      nnx=mod(n,10000)
	      nn2=int(nnx/1000)
	      nny=mod(nnx,1000)
	      nn3=int(nny/100)
	      nnz=mod(nny,100)
	      nn4=int(nnz/10)
	      nn5=mod(nnz,10)
	      filename1(11:11)=char(48+nn1)
	      filename1(12:12)=char(48+nn2)
	      filename1(13:13)=char(48+nn3)
	      filename1(14:14)=char(48+nn4)
	      filename1(15:15)=char(48+nn5)
	   else
	      STOP ' too many time steps ! '
	   endif
	   write(iout,'(a4,5x,a5,5x,11a)')
     :          'read',filename1(11:15),'Hamiltonian'
       	   open(31,file=filename1,status='old')
	   do i = 1 , ndim
              read(31,'(500d10.4)') (h(i,j,neff),j=1,ndim)
	   enddo
	   close(31)
	enddo
	
	call storesteps(nstp1,dotsa)
	filename2(1:10)='hami.00001'
       	open(31,file=filename2,status='old')
	do i = 1 , ndim
           read(31,'(500d10.4)') (hzero(i,j),j=1,ndim)
	enddo
	close(31)
	
	if(static) then
	   do k = 1 , nstp
	      do i = 1 , ndim
	         do j = 1 , ndim
	            h(j,i,k)=hzero(j,i)
	         enddo
	      enddo
	   enddo
	endif

	if(.not.restart) then
	   diffedon=edon-h(idon,idon,1)
	   diffeacc=eacc-h(iacc,iacc,1)
	   open(73,file='diff')
	   write(73,*) diffedon,diffeacc
	   close(73)
	else
	   open(73,file='diff')
	   read(73,*) diffedon,diffeacc
	   close(73)
	endif
	
	if(freeze) then
	   do i = 1 , ndim
	      do j = 1 , ndim
	         if((i.eq.idon.and.j.eq.idon).or.
     :	            (i.eq.iacc.and.j.eq.iacc)) goto 1220	         
	         do k = 1 , nstp
	            h(j,i,k)=hzero(j,i)
	         enddo
 1220            continue	         
	      enddo
	   enddo
	endif

	do n = 1 , nstp

           if(.not.fixdna) then
	     h(idon,idon,n)=diffedon+h(idon,idon,n)
             h(iacc,iacc,n)=diffeacc+h(iacc,iacc,n)
           else
             h(idon,idon,n)=edon
             h(iacc,iacc,n)=eacc
           endif
	   eta1=eta
	   eta2=eta
	   eta3=eta
	   do i = 1 , ndim
	      if(i.ne.idon.and.i.ne.iacc) then
	         if(abs(h(idon,idon,n)-h(i,i,n)).lt.10.0) then
	            h(idon,i,n)=h(idon,i,n)/eta1
	            h(i,idon,n)=h(i,idon,n)/eta1
	         endif
	         if(abs(h(idon,i,n)).gt.0.5) then
	            h(idon,i,n)=h(idon,i,n)/eta2
	            h(i,idon,n)=h(i,idon,n)/eta2
	         endif
	         if(abs(h(iacc,iacc,n)-h(i,i,n)).lt.10.0) then
	            h(iacc,i,n)=h(iacc,i,n)/eta1
	            h(i,iacc,n)=h(i,iacc,n)/eta1
	         endif
	         if(abs(h(iacc,i,n)).gt.0.5) then
	            h(iacc,i,n)=h(iacc,i,n)/eta2
	            h(i,iacc,n)=h(i,iacc,n)/eta2
	         endif
	      endif
	   enddo
	enddo
	
	return
	end