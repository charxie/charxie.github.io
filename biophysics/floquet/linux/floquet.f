c build the floquet hamiltonian for a multimode bridge

	subroutine flo(n,nflo,ncut,hamf,nstep,tau,omax)

	implicit real*8 (a-h,o-z)
	parameter (neigd=2000)
	dimension hamf(n,n,0:ncut)
	dimension hflo(nflo,nflo)
	dimension prob(n)
	complex conjugate
	complex grnf(n)
	common /eigv/ omcm(neigd),zr(neigd,neigd),zi(neigd,neigd),
     :                ne,nblw
	common /ida/ idon, iacc

	write(*,*) ' ----------------------------------'
	write(*,*) ' performing Floquet analysis'
	write(*,*)

	open(45,file='prob')
	open(46,file='flom')
	open(47,file='bocc')

        twopi=8.d0*datan(1.d0)
	period=nstep*tau
	comega=0.5*twopi/period
	
	do i = 1 , n
	   do j = 1 , n
	      do k = 0 , ncut
	         hamf(i,j,k)=0.5*hamf(i,j,k)
	      enddo
	   enddo
	enddo
	
	leq=(nflo-n)/n
	max=int(leq/2)
	min=max-leq
	
	write(*,*)' frequency cutoff: min = ', min, ' max = ', max
	write(*,*)
	
        do k = 1 , nflo
           do m = 1 , nflo
              k1=mod(k-1,n)+1
              m1=mod(m-1,n)+1
              k2=int((k-1)/n)+min
              m2=int((m-1)/n)+min
              if(m2.eq.0.and.k2.eq.0.and.k1.eq.idon.and.m1.eq.idon) 
     :           m0=m
           
              if(k.eq.m) then
                 hflo(k,m)=hamf(k1,m1,0)+k2*comega
              elseif(k2.eq.m2) then
                 hflo(k,m)=hamf(k1,m1,0)
              else
                 hflo(k,m)=hamf(k1,m1,abs(k2-m2))
              endif
           enddo
	enddo	
	
	do k = 1 , nflo
	   write(46,'(500f10.4)') (hflo(k,m),m=1,nflo)
	enddo
	close(46)	

	call eigen(nflo,hflo)
	
c	write(*,*) ' quasi-energies '
c	write(*,'(100f8.4)') (omcm(i),i=1,nflo)

	time=0.0
	delt=1.0
	omax=0.0

c>>> search for largest bridge occupancy
	
 100    continue

	sum=0.0
	do i = 1 , n
	   grnf(i)=0.0
	   do m = 1 , nflo
	      if(mod(m-1,n)+1.eq.i) then
	         do k = 1 , nflo
	         grnf(i)=grnf(i)+(zr(m ,k)+(0.0,1.0)*zi(m ,k))
     :                          *(zr(m0,k)-(0.0,1.0)*zi(m0,k))
     :                          *cexp(-(0.0,1.0)*cmplx(omcm(k)*time))
	         enddo
	      endif
	   enddo
	   prob(i)=grnf(i)*conjugate(grnf(i))
	   sum=sum+prob(i)
	enddo
	if(dabs(sum-1.d0).gt.0.1) stop ' divergence error !'
	pbridge=sum-prob(idon)-prob(iacc)
        if(pbridge.gt.omax) omax=pbridge
        	
	write(45,'(100f10.5)') time/1.51653,(prob(i),i=1,n)
	write(47,'(100f10.5)') time/1.51653,
     :                         pbridge,prob(idon),prob(iacc)
	
	time=time+delt
	if(time.le.100) goto 100
	
 200	close(45)
 	close(47)
 	
 	return
	end