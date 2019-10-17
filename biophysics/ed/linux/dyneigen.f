c*** perform dynamics
c*** integrate the Schroedinger equation

	subroutine dyneigen(tau,nstp,ndim,tmd,eg,ev,ft,neigen,
     :                      emat,einv,cofft)
  	implicit real*8 (a-h,o-z)
  	real*4 tmd(nstp),eg(ndim,nstp),cofeg(ndim,nstp,3)
  	real*4 ev(ndim,ndim,nstp)
  	dimension emat(ndim,ndim),einv(ndim,ndim)
  	real*4 ft(ndim,ndim,nstp),cofft(ndim,ndim,nstp,3)
  	complex*8 psi0(ndim)
	complex*8 psi(ndim)
	dimension prob(ndim)
	logical restart,dotsa,rescale,twostate
        common/ida/time0,edon,eacc,idon,iacc,i0,j0,restart,nstep
        common /io/ iin,iout

	dotsa=.false.
	rescale=.false.
	twostate=.false.
	
        do j = 1 , ndim
           do i = 1 , ndim
              emat(i,j)=0.d0
           enddo
        enddo
        
        do i = 1 , ndim
           prob(i)=0.d0
           psi0(i)=(0.0,0.0)
           psi(i)=(0.0,0.0)
        enddo
	
	do n = 1 , nstp
	   do j = 1 , ndim
	      do i = 1 , ndim
	         cofft(i,j,n,1)=0.0
	         cofft(i,j,n,2)=0.0
	         cofft(i,j,n,3)=0.0
	      enddo
	      cofeg(j,n,1)=0.0
	      cofeg(j,n,2)=0.0
	      cofeg(j,n,3)=0.0
	   enddo
	enddo
	
	call interpol(nstp,ndim,tmd,ft,cofft)
	call interpol2(nstp,ndim,tmd,eg,cofeg)
	
        if(.not.restart) then
           do i = 1 , ndim
              do j = 1 , ndim
                 emat(j,i)=dble(ev(j,i,1))
              enddo
           enddo
           call initeigen(idon,iacc,ndim,psi0,emat,einv)
           time0=0
           do i = 1 , ndim
              psi(i)=psi0(i)
           enddo

c           pacr1=0.0
c           paci1=0.0
c           do mk = 1 , ndim
c              pacr=0.0
c              paci=0.0
c              do i = 1 , ndim
c                 pacr=pacr+ev(mk,i,1)*real(psi(i))
c                 paci=paci+ev(mk,i,1)*aimag(psi(i))
c              enddo
c              pacr1=pacr1+pacr*pacr
c              paci1=paci1+paci*paci
c           enddo
c           write(6,*) pacr1,paci1
c           stop

        else
           call readrst(ndim,time0,psi,dotsa)
        endif
        
	write(iout,*)
	write(iout,*) ' ------------------------------------------ '
	write(iout,*) ' starting dynamics in the eigen space'
	write(iout,*)

	nwrit=nstep/(nstp-1)
	istep=1
       
 20     if(istep.ge.nstep) go to 30
 
	time=istep*tau+time0
	itt=(time-tmd(1))/(tmd(2)-tmd(1))+1
	dtt=time-tmd(itt)

	call predictor(tau,ndim,psi)
	call correcteigen(tau,ndim,psi,time,nstp,
     :                    tmd,eg,cofeg,ft,cofft) 
	
  	istep=istep+1
        
	if(mod(istep,100).eq.0) then
	
           sum=0.0
        
           do i = 1 , ndim
              prob(i)=real(psi(i))**2+aimag(psi(i))**2
              sum=sum+prob(i)
           enddo
           if(dabs(sum-1.d0).gt.0.1d0) stop ' Total population > 1 !'
           if(dabs(sum-1.d0).gt.0.01d0.and.rescale) then
              temsum=dsqrt(sum)
              do i = 1 , ndim
                 psi(i)=psi(i)/temsum
              enddo 
           endif

           write(iout,'(5f12.10)') (prob(i),i=1,ndim)
   	   write(iout,'(a5,f20.5,3x,2f10.5)')'time=',time/1.51653,
     :                                       sum
           write(46,'(f10.5,2f20.15)') time/1.51653,
     :           prob(neigen-1),prob(neigen)
     
           if(.not.twostate) then
              ibeg=1
              iend=ndim
           else
              ibeg=neigen-1
              iend=neigen
           endif
           pacr=0.0
           paci=0.0
           pdor=0.0
           pdoi=0.0
           ix=itt
           do i = ibeg , iend
              pacr=pacr+ev(iacc,i,ix)*real(psi(i))
              paci=paci+ev(iacc,i,ix)*aimag(psi(i))
              pdor=pdor+ev(idon,i,ix)*real(psi(i))
              pdoi=pdoi+ev(idon,i,ix)*aimag(psi(i))
           enddo
           write(iout,*) time/1.51653, 
     :           dmin1(pacr*pacr+paci*paci,pdor*pdor+pdoi*pdoi)
           write(11,'(2f20.15)') time/1.51653,
     :           dmin1(pacr*pacr+paci*paci,pdor*pdor+pdoi*pdoi)
           
        endif

c	if(mod(istep,5000).eq.0) then
c
c           pacr1=0.0
c           paci1=0.0
c           do mk = 1 , ndim
c              pacr=0.0
c              paci=0.0
c              do i = 1 , ndim
c                 pacr=pacr+ev(mk,i,4)*real(psi(i))
c                 paci=paci+ev(mk,i,4)*aimag(psi(i))
c              enddo
c              pacr1=pacr1+pacr*pacr
c              paci1=paci1+paci*paci
c           enddo
c           write(6,*) pacr1,paci1,pacr1+paci1
c           stop
c        endif

        go to 20
        
 30     call save(ndim,time,dotsa)

	return
	end