c fourier series analysis for time-dependent hamiltonian
c only cosine series allowed, in order to avoid diagonalization problem
c due to the antisymmetry caused by sine series.

	subroutine fourier(ndim,tau,nstep,ncut,hamf)
	implicit real*8 (a-h,o-z)
	dimension time(nstep),hami(ndim,ndim,nstep)
	dimension hamf(ndim,ndim,0:ncut)

	open(91,file='htim')
	open(92,file='hfos')
	open(93,file='hnew')
	
	do nt = 1 , nstep
	   read(91,'(500f10.5)') time(nt),
     :              ((hami(i,j,nt),j=1,ndim),i=1,ndim)
	enddo
	close(91)
	tau=time(2)-time(1)
	timezero=time(1)-tau
	
	do i = 1 , nstep
	   time(i)=time(i)-timezero
	enddo
	period=time(nstep)
	
	twopi=8.d0*datan(1.d0)
	comega=0.5*twopi/period
	deltat=period/nstep
	deltao=comega
	
	do m = 1 , ndim
	   do k = 1 , ndim

	      do i = 0 , ncut
	         hamf(m,k,i)=0.0	         
	         do j = 1 , nstep
	            hamf(m,k,i)=hamf(m,k,i)+hami(m,k,j)
     :                     *dcos(i*deltao*time(j))*deltat
	         enddo
	         hamf(m,k,i)=2.0*hamf(m,k,i)/period
	      enddo
	      
	   enddo
	enddo
	
	do i = 0 , ncut
	   write(92,'(500f10.5)') deltao*i,
     :          ((hamf(m,k,i),m=1,ndim),k=1,ndim)
	enddo
	close(92)
	
	call outham(tau,ndim,nstep,ncut,hami,hamf)
	
	return
	
	do m = 1 , ndim
	   do k = 1 , ndim
	      do i = 1 , nstep
	         hami(m,k,i)=0.5*hamf(m,k,0)
	         do j = 1 , ncut	   
	            hami(m,k,i)=hami(m,k,i)
     :                  +hamf(m,k,j)*dcos(j*deltao*time(i))
	         enddo
	      enddo
	   enddo
	enddo
	
	do nt = 1 , nstep
	   write(93,'(500f10.5)')time(nt),
     :              ((hami(i,j,nt),j=1,ndim),i=1,ndim)
	enddo
	close(93)
	
	return	
        end