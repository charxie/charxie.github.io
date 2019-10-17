c faster calculation of propogation in static case
c the subroutine is also used to check the accuracy of the integrator

	subroutine faster(ndim,nstp,diff)

	implicit real*8 (a-h,o-z)
	parameter (ndmax=500,nsmp=10,neigd=500)
	dimension prob(ndim)
	complex*8 grnf(ndim)
        common/iounit/iin,iout,istart(40),iungf,iunhl	
        common/nda/ndon,nacc,idon,iacc,nb,iold(ndmax),jold(ndmax)
        common/inst/ovlp(ndmax,ndmax),hami(ndmax,ndmax),
     *              dens(ndmax,ndmax)
        common/eigv/omcm(neigd),zr(neigd,neigd),zi(neigd,neigd),ne,nblw
        dimension zr1(ndim,ndim),zr2(ndim,ndim)
        dimension zs(ndim,ndim),zt(ndim,ndim)
        dimension temp(ndim,ndim),tinv(ndim,ndim),indx(ndim)

	write(iout,*) ' ----------------------------------'
	write(iout,*) ' faster propogation'
	write(iout,*)

	open(45,file='popf')
	open(46,file='palf')

        planck=8.0*atan(1.0)
c        span=planck/diff
        span=2000.0
	time=0.0
	ntot=1000
	delt=span/ntot
	ico=0
        if(mod(ndim,2).eq.0) then
           nhaf=ndim*0.5
        else
           nhaf=(ndim+1)*0.5
        endif
	   
	do i = 1 , ndim
	   do j = 1 , i
	      temp(i,j)=ovlp(i,j)
	      temp(j,i)=ovlp(j,i)
	      tinv(i,j)=0.0
	      tinv(j,i)=0.0
	   enddo
	   tinv(i,i)=1.0
	enddo
	call ludcmp(temp,ndim,ndim,indx,d)
	do j = 1 , ndim
	   call lubksb(temp,ndim,ndim,indx,tinv(1,j))
	enddo

	do i = 1 , ndim
	   do j = 1 , ndim
	      zs(i,j)=zr(i,j)
	      zt(j,i)=zr(i,j)
	      temp(i,j)=ovlp(i,j)
	   enddo
	enddo
		
c        call mprod(ndim,temp,zs,zr1)
c        call mprod(ndim,zt,temp,zr2)
        call strassen(ndim,nhaf,temp,zs,zr1)
        call strassen(ndim,nhaf,zt,temp,zr2)
        
 100    continue

	sum=0.0
	do i = 1 , ndim
	   grnf(i)=0.0
	   do k = 1 , ndim
	      grnf(i)=grnf(i)+zr1(i,k)*zr2(k,idon)
     :                       *cexp(-(0.0,1.0)*cmplx(omcm(k)*time))
	   enddo
	enddo
	
	do i = 1 , ndim
	   do j = 1 , ndim
	      temp(i,j)=grnf(i)*conjg(grnf(j))
	   enddo
	enddo
	
c	call mprod(ndim,tinv,temp,zt)
	call strassen(ndim,nhaf,tinv,temp,zt)
	
	do i = 1 , ndim
	   prob(i)=zt(i,i)
	   sum=sum+prob(i)
	enddo
	
c	write(iout,'(10f8.4)')(prob(i),i=1,ndim),sum
	
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