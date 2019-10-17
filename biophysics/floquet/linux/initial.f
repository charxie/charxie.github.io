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


c define the dynamical parameters

	subroutine defdynpar(tau,n,perio,ampli,delta,omega,
     :	                     perhop,amphop,delhop,omehop,    	
     :                       ncase,incromega,dperiod)
	implicit real*8 (a-h,o-z)
	dimension perio(n),ampli(n),delta(n),omega(n)
	dimension perhop(n,n),omehop(n,n),amphop(n,n),delhop(n,n)
	common /ranum/ rann(1000)		
	logical static

	twopi=8.d0*datan(1.d0)
	static=.false.
	
c>>>case 1: The whole bridge vibrates with the same frequency

	if(ncase.eq.1) then

	   do i = 1 , n

	      if(i.eq.1) then
                 perio(i)=10000000.d0
                 ampli(i)=0.0d0
  	         delta(i)=-270.d0
  	      elseif(i.eq.n) then
  	         perio(i)=10000.d0
  	         ampli(i)=0.d0
  	         delta(i)=270.d0
  	      else
  	         perio(i)=incromega*dperiod
	         ampli(i)=3.5d0
	         delta(i)=0.d0
	      endif

	   enddo
	   
	   do i = 1 , n
	      do j = 1 , n
	      
	         if(abs(j-i).le.1) then
	   
	         perhop(i,j)=1.0d0
	         amphop(i,j)=0.0d0
	   	 delhop(i,j)=0.0d0
	   	 
	   	 endif
	   
	      enddo
	   enddo

	endif

	do i = 1 , n
	   if(dabs(perio(i)).gt.1.0d-6) then 
	      omega(i)=twopi/perio(i)
	   else
	      omega(i)=0.d0
	   endif
	enddo
	
	do i = 1 , n
	   do j = 1 , n
	      if(dabs(perhop(i,j)).gt.1.0d-6) then
	         omehop(i,j)=twopi/perhop(i,j)
	      else
	         omehop(i,j)=0.d0
	      endif
	   enddo
	enddo
	
	do i = 1 , n
	   do j = 1 , i
	      if(i.eq.j) then
	         perhop(i,j)=0.d0
	         omehop(i,j)=0.d0
	         delhop(i,j)=0.d0
	         amphop(i,j)=0.d0
	      else
	         perhop(i,j)=perhop(j,i)
	         omehop(i,j)=omehop(j,i)
	         delhop(i,j)=delhop(j,i)
	         amphop(i,j)=amphop(j,i)
	      endif
	   enddo
	enddo
	
	if(static) then
	   do i = 1 , n
	      ampli(i)=0.0
	   enddo
	   do i = 1 , n
	      do j = 1 , n
	         amphop(i,j)=0.0
	      enddo
	   enddo
	endif

	return
	end