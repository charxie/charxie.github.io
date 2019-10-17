c*** perform dynamics
c*** integrate the Schroedinger equation

	subroutine dyn(nstp,tau,ndim,tmd,h,cof)
  	implicit real*8 (a-h,o-z)
  	real*4 tmd(nstp),h(ndim,ndim,nstp),cof(ndim,ndim,nstp,3)
  	complex*8 psi0(ndim),psi1(ndim)
	complex*8 project(ndim),psi(ndim),phi(ndim,ndim),tot
        complex*8 rho(ndim,ndim),flow(ndim,ndim)
	dimension prob(ndim)
	logical restart,dotsa,gear,ruku,rescale
        common /ida/ idon,iacc,i0,j0,restart,time0,nstep,edon,eacc
        common /io/ iin,iout

	dotsa=.false.
	gear=.true.
	ruku=.false.
	rescale=.false.

	open(99,file='totp')

        call init(idon,iacc,ndim,phi,psi0,psi1)
        if(.not.restart) then
           time0=0
           do i = 1 , ndim
              psi(i)=psi0(i)
           enddo
        else
           if(gear) call readrst(ndim,time0,psi,dotsa)
           if(ruku) call readrkns(ndim,time0,psi,dotsa)
        endif
	
	write(iout,*)
	write(iout,*) ' ------------------------------------------ '
	write(iout,*) ' starting dynamics '
	write(iout,*)

	nwrit=nstep/(nstp-1)
	istep=1
       
 20     if(istep.gt.nstep) go to 30
 
	time=istep*tau+time0
	itt=(time-tmd(1))/(tmd(2)-tmd(1))+1
	dtt=time-tmd(itt)
c	if(mod(istep,nwrit).eq.0) then
c           hij=((cof(i0,j0,itt,3)*dtt+cof(i0,j0,itt,2))*dtt
c     :          +cof(i0,j0,itt,1))*dtt+h(i0,j0,itt)
c     	   write(12,'(50f20.10)') time/1.51653,hij
c     	endif

        if(gear) then
	   call predictor(tau,ndim,psi)
	   call corrector(tau,ndim,psi,time,nstp,tmd,h,cof) 
	endif
	if(ruku) then
	   call rk4(tau,ndim,psi,time,nstp,tmd,h,cof)
	endif
	
  	istep=istep+1
        
	if(mod(istep,100).eq.0) then
	
c project onto a nonorthogonal basis, disabled for CNDO calculations
c	   call proj(ndim,psi,phi,project)
c	   call dent(ndim,project,rho,hami,flow)
	
           sum=0.0
           tot=(0.0,0.0)
        
           do i = 1 , ndim
c              prob(i)=(real(project(i)))**2+(aimag(project(i)))**2
              prob(i)=real(psi(i))**2+aimag(psi(i))**2
              sum=sum+prob(i)
           enddo
           if(dabs(sum-1.d0).gt.0.1d0) stop ' Total population > 1 !'
c
c rescale the wave function to control error propogation
c error increases with the size of the hamiltonian and the 
c length of the simulation time, in the case of weak coupling,
c rescaling is recommended.
c
           if(dabs(sum-1.d0).gt.0.01d0.and.rescale) then
              temsum=dsqrt(sum)
              do i = 1 , ndim
                 psi(i)=psi(i)/temsum
              enddo 
           endif

c
c calculate the total electronic energy, this number may be
c useful when we talk about k_BT
c
           do i = 1 , ndim
              do j = 1 , ndim
                 hij=((cof(i,j,itt,3)*dtt+cof(i,j,itt,2))*dtt
     :                +cof(i,j,itt,1))*dtt+h(i,j,itt)
c                 tot=tot+project(j)*conjg(project(i))*hij
                 tot=tot+psi(j)*conjg(psi(i))*hij
              enddo
           enddo
           
	   write(iout,'(a5,f20.5,3x,2f10.5)')'time=',time/1.51653,
     :                                       sum,real(tot)
           write(iout,'(5f12.10)') (prob(i),i=1,ndim)
           
           write(11,'(f20.10,2f20.15)')
     :                 time/1.51653,prob(iacc)
c          write(99,'(500f10.5)')
c     :                 time/1.51653,(prob(i),i=1,ndim)
           write(44,'(f20.10,f20.15)')
     :                 time/1.51653,real(tot)
           write(46,'(f20.10,f20.15)') 
     :                 time/1.51653,sum-prob(idon)-prob(iacc)

        endif

        go to 20
        
 30     if(gear) call save(ndim,time,dotsa)
        if(ruku) call saverkns(ndim,time,psi,dotsa)

	close(99)
	return
	end