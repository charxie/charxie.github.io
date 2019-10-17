c define the localized electronic states
c phi(0,n)  ...... donor state
c phi(i,n)  ...... bridge states
c phi(n,n)  ...... acceptor state
c psi0(n)   ...... initial state
c psi1(n)   ...... expected state

	subroutine init(n,phi,psi0,psi1) 
	  	
	implicit real*8 (a-h,o-z)
	complex phi(n,n),psi0(n),psi1(n)
		
	do i = 1 , n
	   do j = 1 , n
     	      phi(i,j)=(0.0,0.0)
	      if(i.eq.j) phi(i,j)=(1.0,0.0)
	   enddo
	enddo
	
	do i = 1 , n
	   psi0(i)=phi(1,i)
	   psi1(i)=phi(n,i)
	enddo
	
	return
	end


c define parameters

	subroutine defdynpar(tau,n)
	implicit real*8 (a-h,o-z)
	parameter (nmax=50)
	common/ranum/rann(1000)		
	common/para/perio(nmax,nmax),omega(nmax,nmax),
     :              ampli(nmax,nmax),delta(nmax,nmax)

	twopi=8.d0*datan(1.d0)
	iseed=1234
	
	call rnum(iseed)
	
	do i = 1 , n
	   if(i.eq.1) then
              perio(i,i)=10000000.d0
              ampli(i,i)=0.0d0
  	      delta(i,i)=-270.d0
  	   elseif(i.eq.n) then
  	      perio(i,i)=10000.d0
              ampli(i,i)=0.d0
  	      delta(i,i)=270.d0
  	   else
  	      perio(i,i)=10000.0d0*rann(800+i)+10.0
	      ampli(i,i)=0.d0*(1+rann(700+i))
	      delta(i,i)=360.d0*rann(900+i)
	   endif
	enddo
	   
	do i = 1 , n
	   do j = i , n
	      perio(i,j)=1000.d0*rann(900+i)+10.0
	      ampli(i,j)=0.0d0*(1+rann(700+i))
	      delta(i,j)=360.d0*rann(800+i)
	   enddo
	enddo
	
	do i = 1 , n
	   do j = 1 , n
	      if(dabs(perio(i,j)).gt.1.0d-6) then
	         omega(i,j)=twopi/perio(i,j)
	      else
	         omega(i,j)=0.d0
	      endif
	   enddo
	enddo
	
	return
	end