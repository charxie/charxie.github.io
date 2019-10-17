c elctron dynamics using predictor-corrector method

	subroutine dyn(tau,n,nstep,restart,static,
     :                 ED,EA,EB,VD1,VBB,VNA,
     :	               ampli,omega,delta,
     :                 amphop,omehop,delhop,
     :                 omax)
     
	implicit real*8 (a-h,o-z)
	dimension esit(n),hopi(n,n),hami(n,n)
	complex psi0(n),psi(n),psi1(n)
	complex phi(n,n),project(n),rho(n,n),flow(n,n),tot
	complex conjugate
	dimension prob(n)
	dimension perio(n),omega(n),ampli(n),delta(n)
	dimension perhop(n,n),omehop(n,n),amphop(n,n),delhop(n,n)
	logical restart,static

        call init(n,phi,psi0,psi1)

	write(*,*)
	write(*,*) ' ------------------------------------------ '
	write(*,*) ' starting predictor-corrector dynamics '
	write(*,*)

	if(restart) then
	   open(11,file='popu',status='unknown',access='append')
	   open(12,file='esit',status='unknown',access='append')
	else
	   open(11,file='popu',status='unknown',access='sequential')
	   open(12,file='esit',status='unknown',access='sequential')
	endif
	
	nwrit=nstep/1000
	omax=0.0

	if(.not.restart) then
	   istep0=0
      	   do i = 1 , n
	      psi(i)=psi0(i)
	   enddo
	else
	   call readrst(n,istep0,psi)
	endif
	istep=istep0+1
       
 20     if(istep-istep0.gt.nstep) go to 30

	if(.not.static) then
      	   call hamiltonian(tau,istep,n,ED,EA,EB,VD1,VBB,VNA,
     :	                    ampli,omega,delta,
     :                      amphop,omehop,delhop,     
     :	                    esit,hopi,hami)
        else
	   call hamiltonian0(ED,EA,EB,VD1,VBB,VNA,n,esit,hopi,hami)
	endif

	if(mod(istep,nwrit).eq.0) 
     :	   write(12,'(50f10.5)') istep*tau/1.51653,
     :                               (esit(i),i=2,n-1)
     
	call predictor(tau,n,hami,psi)
	call corrector(tau,n,hami,psi) 
 
  	istep=istep+1
	
c	if(mod(istep,1000).eq.0) then
c	   do i = 1 , n
c	      write(*,'(100f7.3)') (hami(i,j),j=1,n)
c           enddo
c           write(*,*)
c        endif

	call proj(n,psi,phi,project)
        do i = 1 , n
           prob(i)=(real(project(i)))**2+(aimag(project(i)))**2
        enddo
	pbridge=1.0-prob(1)-prob(n)
	if(pbridge.gt.omax) omax=pbridge
		   
c	call dent(n,project,rho,hami,flow)
c	write(13,'(i10,3x,20f10.5)') istep*tau/1.51653,
c     :             (real(flow(1,i)),i=2,n)

        
	if(mod(istep,nwrit).eq.0) then
	
           sum=0.0
           tot=(0.0,0.0)
           do i = 1 , n
              sum=sum+prob(i)
           enddo
           do i = 1 , n
              do j = 1 , n
                 tot=tot+
     :               project(j)*conjugate(project(i))*hami(i,j)
              enddo
           enddo
 		
c	   write(*,'(a7,3f10.5)') ' step =', istep*tau, sum, real(tot)
	   if(dabs(sum-1.d0).gt.0.05) stop ' divergence error !'
c           write(*,'(100f10.5)') (prob(i),i=1,n)
           
           write(11,'(100f20.10)') istep*tau/1.51653,
     :                             prob(1),pbridge,prob(n)
           

        endif

        go to 20
        
 30	close(11)
 	close(12)

        call save(n,istep)
        
	return	
	end