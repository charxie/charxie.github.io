        subroutine fsdyn(ndim,hami,omax,iunit)

	implicit real*8 (a-h,o-z)
	logical restart,gear,ruku
	complex psi0(ndim),psi(ndim),psi1(ndim)
	complex phi(ndim,ndim),project(ndim),tot
	dimension hami(ndim,ndim),prob(ndim)
	common/edyn/tau,nstep,restart

        restart=.false.
        gear=.false.
        ruku=.true.
        omax=0.0

	call init(ndim,phi,psi0,psi1)

	write(*,*)
	write(*,*) ' ------------------------------------------ '
	write(*,*) ' starting fs dynamics '
	write(*,*)
	
	nwrit=nstep/500

	istep0=0
      	do i = 1 , ndim
	   psi(i)=psi0(i)
        enddo
	istep=istep0+1
       
 20     if(istep-istep0.gt.nstep) go to 30

        if(gear) then
	   call predictor(tau,ndim,hami,psi)
	   call corrector(tau,ndim,hami,psi) 
	endif
	if(ruku) then
	   call rk4(tau,ndim,hami,psi)
	endif
 
  	istep=istep+1
        
	if(mod(istep,nwrit).eq.0) then

	   call proj(ndim,psi,phi,project)
	
           sum=0.0
           tot=(0.0,0.0)
        
           do i = 1 , ndim
              prob(i)=(real(project(i)))**2+(aimag(project(i)))**2
              sum=sum+prob(i)
           enddo

           do i = 1 , ndim
              do j = 1 , ndim
                 tot=tot+
     :               project(j)*conjg(project(i))*hami(i,j)
              enddo
           enddo
 		
	   write(*,*) ' step =', istep, sum, real(tot)
           write(*,'(10f10.5)') (prob(i),i=1,ndim)
           
           if(ndim.eq.2) then
              write(iunit,'(20f20.10)') istep*tau,prob(2)
           endif
           if(ndim.eq.3) then
              write(iunit,'(10f20.10)') istep*tau,prob(2),prob(3)
              if(prob(2).gt.omax) omax=prob(2)
           endif
           
        endif

        go to 20
        
c-------------------------end-------------------------

 30     return
        end