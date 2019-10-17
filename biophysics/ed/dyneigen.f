c*** perform dynamics
c*** integrate the Schroedinger equation

	subroutine dyneigen(nstp,tau,ndim,tmd,eg,ev,ft,neigen,toh,h,cof)
  	implicit real*8 (a-h,o-z)
	parameter (neigd=248)
  	real*4 tmd(nstp),eg(ndim,nstp),cofeg(ndim,nstp,3)
  	real*4 ev(ndim,ndim,nstp)
  	dimension emat(ndim,ndim)
  	real*4 ft(ndim,ndim,nstp),cofft(ndim,ndim,nstp,3)
  	complex*8 psi0(ndim)
	complex*8 psi(ndim)
	dimension prob(ndim)
	real*4 toh(nstp),h(ndim,ndim,nstp),cof(ndim,ndim,nstp,3)
	dimension htem(ndim,ndim),etem(ndim,ndim)
      	logical restart,dotsa,rescale,twostate,firstTime,coherence
        common /ida/ idon,iacc,i0,j0,restart,time0,nstep,edon,eacc
        common /io/ iin,iout
	common /eigv/ omcm(neigd),zr(neigd,neigd),zi(neigd,neigd),
     :                ne,nblw
	common /msign/ msign(neigd,2)

	dotsa=.false.
	rescale=.false.
	twostate=.true.
	firstTime=.true.
	coherence=.false.
	
	call interpol(nstp,ndim,tmd,ft,cofft)
	call interpol2(nstp,ndim,tmd,eg,cofeg)
	
        if(.not.restart) then
        
           do i = 1 , ndim
              do j = 1 , ndim
                 emat(j,i)=dble(ev(j,i,1))
              enddo
           enddo
           call initeigen(idon,iacc,ndim,psi0,emat)
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
c                pacr=pacr+ev(mk,i,1)*real(psi(i))
c                paci=paci+ev(mk,i,1)*aimag(psi(i))
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
       
 20     if(istep.gt.nstep) go to 30
 
	time=istep*tau+time0
	itt=(time-tmd(1))/(tmd(2)-tmd(1))+1
	dtt=time-tmd(itt)
	itt1=(time-toh(1))/(toh(2)-toh(1))+1
	dtt1=time-toh(itt)

	call predictor(tau,ndim,psi)
	call correcteigen(tau,ndim,psi,time,nstp,
     :                    tmd,eg,cofeg,ft,cofft) 
	
  	istep=istep+1
        
	if(mod(istep,100).eq.0) then
	
c	   goto 1008
	
	   do i = 1 , ndim
	      do j = 1 , ndim
                 htem(j,i)=((cof(j,i,itt1,3)*dtt1+cof(j,i,itt1,2))*dtt1
     :                      +cof(j,i,itt1,1))*dtt1+h(j,i,itt1)
     	         etem(j,i)=0.0	      
              enddo
              etem(i,i)=1.0
	   enddo
	   
           call eigen(-1000.d0,1000.d0,.true.,ndim,htem,etem)
           
           do i = 1 , ndim
              isign=1
              if(msign(i,2)*zr(msign(i,1),i).lt.0.0) isign=-1
              do j = 1 , ndim
                 zr(j,i)=zr(j,i)*isign
              enddo
           enddo

 1008      continue

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
c              pacr=pacr+ev(iacc,i,ix)*real(psi(i))
c              paci=paci+ev(iacc,i,ix)*aimag(psi(i))
c              pdor=pdor+ev(idon,i,ix)*real(psi(i))
c              pdoi=pdoi+ev(idon,i,ix)*aimag(psi(i))
              pacr=pacr+zr(iacc,i)*real(psi(i))
              paci=paci+zr(iacc,i)*aimag(psi(i))
              pdor=pdor+zr(idon,i)*real(psi(i))
              pdoi=pdoi+zr(idon,i)*aimag(psi(i))
           enddo
           if(twostate.and.coherence) then
              do i = 1 , ndim
                 if(i.ne.ibeg.and.i.ne.iend) then
                    cohe1=cohe1+ev(iacc,i,ix)*psi(i)
                 endif
              enddo
              cohe2=cohe1
              cohe1=cohe1*psi(ibeg)*ev(iacc,ibeg,ix)
              cohe2=cohe2*psi(iend)*ev(iacc,iend,ix)
           endif
           if(coherence) then           
              write(iout,*) time/1.51653,pacr*pacr+paci*paci+
     :                      real(cohe1+cohe2)*2.0
              write(11,'(20f20.15)') time/1.51653,pacr*pacr+paci*paci+
     :                               real(cohe1+cohe2)*2.0
           else 
              write(iout,*) time/1.51653, pacr*pacr+paci*paci
              write(11,'(20f20.15)') time/1.51653, pacr*pacr+paci*paci
c             write(11,'(20f20.15)') time/1.51653, pacr*pacr+paci*paci,
c     :                            ev(iacc,neigen-1,ix),psi(neigen-1),
c     :                            ev(iacc,neigen,ix),psi(neigen)
           endif
           
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