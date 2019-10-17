   	program BIOED
   	
c********************************************************************
c BIOED:Electronic Dynamics for Electron Transfer in Proteins.
c Version 1.0 (Nov. 1998)
c
c ***     MD trajectory from CHARMM24
c ***     Hamiltonian from EHMACC
c ***     Integrator GEAR,RUKU
c
c********************************************************************   	
   	
	implicit real*8 (a-h,o-z)
c	include 'PARAM.include'
	parameter(neigd=500,nsmp=10,ndmax=500,natmax=200)
	logical static,restart,dotsa,check1,check3,rite
	real*4 tmd(nsmp),st(ndmax,ndmax,nsmp),ht(ndmax,ndmax,nsmp),
     *         ft(ndmax,ndmax,nsmp)
        real*4 sao,trt	
	real*4 st2(2,2,nsmp),ht2(2,2,nsmp)
        common/iounit/iin,iout,istart(40),iungf,iunhl
        common/nda/ndon,nacc,idon,iacc,nb,iold(ndmax),jold(ndmax)
        common/nonadi/sao(ndmax,ndmax,nsmp),trt(ndmax,ndmax,nsmp)
        common/option/restart,dotsa,nstep,unit
        real*4 s0(ndmax,ndmax),h0(ndmax,ndmax)
	
	data ndon/2/
	data nacc/18/
	data tau/0.001d0/
	data nstep/1000/
	data iseed/87654321/
	
	twopi=8.d0*datan(1.d0)
	unit=1.51653*1.0
	rite=.true.
	restart=.false.
	dotsa=.false.
	check1=.false.
	check3=.false.
	static=.true.
	nstp=1

        ipro=5
	iin=7
	iout=6
	iungf=33
	iunhl=34
	
800	write(iout,*) ' restart ? (y/n) '
	read(ipro,'(a)') yesorno
	if(yesorno.eq.'y') then
	   restart=.true.
	elseif(yesorno.eq.'n') then
	   restart=.false.
	else
	   write(iout,*) ' answer [y]es or [n]o '
	   goto 800
	endif

	open(iungf,file='greenfun',status='unknown')
	open(iunhl,file='homolumo',status='unknown')

	write(iout,*)'-------------------- BIOED V1.0 -------------------'
	write(iout,*)'Biophysics Group, University of Cyprus, 1998'
	write(iout,*)'---------------------------------------------------'

	write(iout,*)
	write(iout,*) 'No. of steps = ', nstp
	write(iout,*)
	call readh(tau,natom,ndim,nstp,tmd,st,ht,ft)
	if(.not.check3) call friction(natom,ndim,nstp,tmd,ft)

	call part1(natom,ndim,nstp,tmd,st,ht,st2,ht2)
	
 	if(.not.check3.and..not.dotsa) then
 	   call fribo(ndim,nstp,tmd,ft)
 	else
 	   do n = 1 , nstp
 	      do j = 1 , ndim
 	         do i = 1 , ndim
 	            ft(i,j,n)=0.0
 	         enddo
 	      enddo
 	   enddo
 	endif
	
	if(rite) then
	   i0=6
	   j0=6
	   open(42,file='recd')
	   do n = 1 , nstp
	      write(42,'(i8,500f8.4)') n,st(i0,j0,n),ht(i0,j0,n),
     *           ft(i0,j0,n)	      
	   enddo
	   close(42)
	endif

	if(static) stop 'END OF STATIC COMPUTATION'

	if(check1) then
	   if(.not.restart) then
	      open(51,file='first',form='unformatted')
	      write(51)((st(i,j,1),ht(i,j,1),j=1,ndmax),i=1,ndmax)
	      close(51)
	   endif
	   if(.not.restart) then
	      nbeg=2
	   else
	      nbeg=1
	   endif
	   open(51,file='first',form='unformatted')
	   read(51)((s0(i,j),h0(i,j),j=1,ndmax),i=1,ndmax)
	   close(51)
	   do n = nbeg , nstp
	      do i = 1 , ndim
	         do j = 1 , ndim
	            st(j,i,n)=s0(j,i)
	            ht(j,i,n)=h0(j,i)
	         enddo
	      enddo
	   enddo
        endif
	
	if(.not.static) then 
	   call part2(tau,natom,ndim,nstp,tmd,st,ht,ft,st2,ht2)
	endif


        STOP
	END