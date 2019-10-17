c prepare tabulated time-dependent hamiltonian file

	subroutine hmod(ndim,nstep)
	implicit real*8(a-h,o-z)
	dimension tmd(nstep),hami(ndim,ndim,nstep)
	character filename*17
	logical static
	common /ida/ idon, iacc
	data idon/04/
	data iacc/17/

	static=.false.

        do i = 1 , nstep
c           tmd(i)=1.51653*1.0*i
	   tmd(i)=1.51653*1.0*i
        enddo

	filename(1:13)='../data/fock.'

	do n = 1 , nstep
	   write(*,*) ' reading ',n,' -th Hamiltonian file ... '
	   if(n.lt.10) then
	      filename(14:16)='000'
	      filename(17:17)=char(48+n)
	   elseif(n.lt.100) then
	      nn1=mod(n,10)
	      nn2=int(n/10)
	      filename(14:15)='00'
	      filename(16:16)=char(48+nn2)
	      filename(17:17)=char(48+nn1)
	   elseif(n.lt.1000) then
	      nn1=int(n/100)
	      nnx=mod(n,100)
	      nn2=int(nnx/10)
	      nn3=mod(nnx,10)
	      filename(14:14)='0'
	      filename(15:15)=char(48+nn1)
	      filename(16:16)=char(48+nn2)
	      filename(17:17)=char(48+nn3)
	   elseif(n.lt.10000) then
	      nn1=int(n/1000)
	      nnx=mod(n,1000)
	      nn2=int(nnx/100)
	      nny=mod(nnx,100)
	      nn3=int(nny/10)
	      nn4=mod(nny,10)
	      filename(14:14)=char(48+nn1)
	      filename(15:15)=char(48+nn2)
	      filename(16:16)=char(48+nn3)
	      filename(17:17)=char(48+nn4)
	   else
	      STOP ' too many time steps ! '
	   endif
       	   open(31,file=filename,status='old')
	   do i = 1 , ndim
              read(31,'(500d10.4)') (hami(i,j,n),j=1,ndim)
	   enddo
	   close(31)
	   hami(idon,idon,n)=-1.32
	   hami(iacc,iacc,n)=-1.00
	enddo
	
	if(static) then
	   do nt = 1 , nstep
	      do i = 1 , ndim
	         do j = 1 , ndim
	            hami(i,j,nt)=hami(i,j,1)
	         enddo
	      enddo
	   enddo
	endif

	open(91,file='htim')

	do n = 1 , nstep
	   write(91,'(500f10.5)')
     :            tmd(n),((hami(i,j,n),j=1,ndim),i=1,ndim)
	enddo

	close(91)

	return
	end

	subroutine htime(tau,n,nstep,
     :                   ED,EA,EB,VD1,VBB,VNA,
     :                   ampli,omega,delta,
     :                   amphop,omehop,delhop)

	implicit real*8(a-h,o-z)
	dimension esit(n),hopi(n,n),hami(n,n)
	dimension omega(n),ampli(n),delta(n)
	dimension omehop(n,n),amphop(n,n),delhop(n,n)
	common /ida/ idon, iacc

	twopi=8.d0*datan(1.d0)
	istep=1
	idon=1
	iacc=n
	open(91,file='htim')
	
 10     if(istep.gt.nstep) goto 20	

      	   call hamiltonian(tau,istep,n,ED,EA,EB,VD1,VBB,VNA,
     :	                    ampli,omega,delta,
     :                      amphop,omehop,delhop,     
     :	                    esit,hopi,hami)

	   write(91,'(200f10.5)')istep*tau,((hami(i,j),j=1,n),i=1,n)
        
           istep=istep+1
        
        goto 10

 20     close(91)
        return
	end



c Subroutine to construct the TB Hamiltonian for electron transfer
c i=1:donor state | i=n:acceptor state
c

	subroutine hamiltonian0(ED,EA,EB,VD1,VBB,VNA,n,esit,hopi,hami)

	implicit real*8 (a-h,o-z)
	dimension esit(n),hopi(n,n),hami(n,n)
	common /ranum/ rann(1000)

	do i = 1 , n
	
	   if(i.eq.1) then
  	      esit(i)=ED
  	   elseif(i.eq.n) then
  	      esit(i)=EA
  	   else
  	      esit(i)=EB
  	   endif
  	   
	enddo
	
	call hop0(VD1,VBB,VNA,n,hopi)
	
	do i = 1 , n
	   do j = 1 , n
	   
	      if(j.eq.i) then 
	         hami(i,j)=esit(i)
	      else
	         hami(i,j)=hopi(i,j)
	      endif
	      
	   enddo
	enddo

	return
	end
	
	
	subroutine hop0(VD1,VBB,VNA,n,hopi)

	implicit real*8 (a-h,o-z)
	dimension hopi(n,n)
	common /ranum/ rann(1000)
	common /decay/ CD,CA,CB(100)

	CD=2.0
	CA=2.0
	do i = 1 , n
	   CB(i)=2.0
	enddo

	do i = 1 , n
	   do j = i+1 , n
	   
	      if(abs(j-i).le.1) then

	         if(i.eq.1) then
                 hopi(i,j)=(VD1)*exp(-CD*real(j-i-1))
          	 elseif(j.eq.n) then
	         hopi(i,j)=(VNA)*exp(-CA*real(j-i-1))
	         else
	         hopi(i,j)=(VBB)*exp(-CB(i)*real(j-i-1))
	         endif
	         
	      else
	         
	         hopi(i,j)=0.d0
	         
	      endif
	      
	   enddo
	enddo

	do i = 1 , n
	   do j = 1 , i
	      if(j.eq.i) then 
	         hopi(i,j)=0.0
	      else
	         hopi(i,j)=hopi(j,i)
	      endif
	   enddo
	enddo

	return
	end

c<<<construct time-dependent Hamiltonian>>>


	subroutine hamiltonian(tau,istep,n,ED,EA,EB,VD1,VBB,VNA,
     :                         ampli,omega,delta,
     :                         amphop,omehop,delhop,
     :	                       esit,hopi,hami)

	implicit real*8 (a-h,o-z)
	dimension esit(n),hopi(n,n),hami(n,n)
	dimension ampli(n),omega(n),delta(n)
	dimension amphop(n,n),omehop(n,n),delhop(n,n)
	common /ranum/ rann(1000)
	
	twopi=8.d0*datan(1.d0)
		
c-->   define the time behavior of the energy levels

	do i = 1 , n
	
	   if(i.eq.1) then
  	      esit(i)=ED+ampli(i)*dcos(istep*omega(i)*tau+delta(i))
  	   elseif(i.eq.n) then
  	      esit(i)=EA+ampli(i)*dcos(istep*omega(i)*tau+delta(i))
  	   else
  	      esit(i)=EB+ampli(i)*dcos(istep*omega(i)*tau+delta(i))
  	   endif
  	   
	enddo
	
	call hop(tau,istep,VD1,VBB,VNA,n,amphop,omehop,delhop,
     :	         hopi)

	do i = 1 , n
	   do j = 1 , n
	   
	      if(j.eq.i) then 
	         hami(i,j)=esit(i)
	      else
	         hami(i,j)=hopi(i,j)
	      endif
	      
	   enddo
	enddo

	return
	end
	
	
	subroutine hop(tau,istep,VD1,VBB,VNA,n,amphop,omehop,delhop,
     :	               hopi)

	implicit real*8 (a-h,o-z)
	dimension hopi(n,n)	
	dimension omehop(n,n),amphop(n,n),delhop(n,n)
	common /decay/ CD,CA,CB(100)
	common /ranum/ rann(1000)
		
	do i = 1 , n
	   do j = i+1 , n

	      if(abs(j-i).le.1) then

	         if(i.eq.1) then
        	    hopi(i,j)=VD1
          	 elseif(j.eq.n) then
	            hopi(i,j)=VNA
	         else
	            hopi(i,j)=VBB
	         endif

	      else
	      
	         hopi(i,j)=0.d0
	         
	      endif
	      
	   enddo
	enddo

	do i = 1 , n
	   do j=1,i
	      if(j.eq.i) then 
	         hopi(i,j)=0.0
	      else
	         hopi(i,j)=hopi(j,i)
	      endif
	   enddo
	enddo

	return
	end