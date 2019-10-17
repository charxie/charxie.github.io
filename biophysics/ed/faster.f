c faster calculation of propogation in static case
c the subroutine is also used to check the accuracy of the integrator

	subroutine faster(ndim,tau,tda)

	implicit real*8 (a-h,o-z)
	parameter (neigd=248)
	dimension prob(ndim)
	complex*8 grnf(ndim)
	logical restart
	common /eigv/ omcm(neigd),zr(neigd,neigd),zi(neigd,neigd),
     :                ne,nblw
        common /ida/ idon,iacc,i0,j0,restart,time0,nstep,edon,eacc
        common /io/ iin,iout

	write(iout,*) ' ----------------------------------'
	write(iout,*) ' faster propogation'
	write(iout,*)

	open(45,file='popf')
	open(46,file='palf')

        planck=8.0*atan(1.0)
c        span=planck/dabs(tda)
        span=1000.0*1.51653
	time=0.0
	ntot=1000
	delt=span/ntot
	ico=0
	
 100    continue

	sum=0.0
	do i = 1 , ndim
	   grnf(i)=0.0
	   do k = 1 , ndim
	      grnf(i)=grnf(i)+(zr(i   ,k)+(0.0,1.0)*zi(i   ,k))
     :                       *(zr(idon,k)-(0.0,1.0)*zi(idon,k))
     :                       *cexp(-(0.0,1.0)*cmplx(omcm(k)*time))
	   enddo
	   prob(i)=grnf(i)*conjg(grnf(i))
	   sum=sum+prob(i)
	enddo
	if(dabs(sum-1.d0).gt.0.01) stop ' divergence error !'
	pbridge=sum-prob(idon)-prob(iacc)
        	
	write(46,'(f30.5,f20.15)') time/1.51653,pbridge
	write(45,'(f30.5,2f20.15)') time/1.51653,
     :                         prob(idon),prob(iacc)
	
	time=time+delt
	ico=ico+1
	if(mod(ico,100).eq.0) write(iout,'(i10,f20.15)') ico,sum
	if(ico.le.ntot) goto 100
	
 200	close(45)
 	close(46)
 	
 	return
	end