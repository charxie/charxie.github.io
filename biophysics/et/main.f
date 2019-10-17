   	program FSRET
   	
c********************************************************************
c few state reduction
c model demo
c********************************************************************   	
   	
	implicit real*8 (a-h,o-z)
	parameter (ndim=9,nred=3,nmax=50)
	logical restart,static
	dimension hami(ndim,ndim)
	dimension h3s(3,3),h2s(2,2)
	dimension dm(ndim,ndim),gf(ndim,ndim)
	common/para/perio(nmax,nmax),omega(nmax,nmax),
     :              ampli(nmax,nmax),delta(nmax,nmax)
	common/eigv/omcm(nmax),zr(nmax,nmax),zi(nmax,nmax),ne,nblw
	common/edyn/tau,nstep,restart
	common/decay/CD,CA,CB(nmax),ED,EA,EM,EB,VD1,VNA,VBB,VBM
	logical findomax

	data tau/0.05d0/
	data nstep/200000/
	data etun/-0.5/
	data ED/0.0/
	data EM/0.0/
	data EA/0.0/
	data EB/4.0/
	data VD1/1.0/
	data VNA/1.0/
	data VBB/1.0/
	data VBM/2.0/

	twopi=8.d0*datan(1.d0)
	restart=.false.
	static=.true.
	findomax=.true.
	etun2s=etun
	etun3s=etun

	open(90,file='input')
	read(90,*) ed,em,ea
	read(90,*) etun2s,etun3s
	read(90,*) ioption1,ioption2,ioption3
	read(90,*) tau,nstep
	close(90)

	call defdynpar(tau,ndim)

	if(restart) then
	   open(11,file='popuns',status='unknown',access='append')
	   open(12,file='popu2s',status='unknown',access='append')
	   open(13,file='popu3s',status='unknown',access='append')
	else
	   open(11,file='popuns',status='unknown',access='sequential')
	   open(12,file='popu2s',status='unknown',access='sequential')
	   open(13,file='popu3s',status='unknown',access='sequential')	
	endif
	
	call hamiltonian0(ndim,hami)
	call eigen(ndim,hami)
	call output(ndim,hami,tau)

        call caltda(ndim,hami,etun2s,tda,tdd,taa)
        h2s(1,1)=tdd
        h2s(1,2)=tda
        h2s(2,1)=tda
        h2s(2,2)=taa
        call hred(ndim,hami,etun3s,nred,h3s)
        idon=1
        imid=int((ndim+1)/2)
        iacc=ndim
        write(*,*) 'test',idon,imid,iacc
        call triple(ndim,hami,etun3s,idon,imid,iacc,h3s)
        
        stop
        

        if(.not.findomax) then
          if(ioption1.eq.1) call nsdyn(ndim,hami,omax)
          if(ioption2.eq.1) call fsdyn(2,h2s,omax,12)
          if(ioption3.eq.1) call fsdyn(3,h3s,omax,13)
        endif

        if(findomax) then
          open(91,file='omaxns')
          open(92,file='omax3s')
          do i = 1 , 300
            em=0.0+0.01*i
            call nsdyn(ndim,hami,omax)
            write(91,'(2f10.5)') em,omax
            call hred(ndim,hami,etun3s,nred,h3s)
            call fsdyn(3,h3s,omax,13)
            write(92,'(2f10.5)') em,omax
          enddo
          close(91)
          close(92)
        endif

        end