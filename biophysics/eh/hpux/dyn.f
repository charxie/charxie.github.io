	subroutine dynamics(tau,ndim,nstp,tmd,st,ht,vtr,vti,cs,ch,cvr,cvi)

c*********************************************************************
c perform electronic dynamics	
c*********************************************************************	   	
   	
	implicit real*8 (a-h,o-z)
	parameter (nsmp=10,ndmax=500)
	character transa,transb
	logical restart,dotsa,rectify,gear,ruku
	real*4 tmd(nsmp),st(ndmax,ndmax,nsmp),ht(ndmax,ndmax,nsmp)
	real*4 vtr(ndmax,ndmax,nsmp),vti(ndmax,ndmax,nsmp)
	complex*8 psi(ndim),psi0(ndim),psi1(ndim),phi(ndim,ndim)
	complex*8 tot
	real*4 cs(ndim,ndim,nstp,3),cvr(ndim,ndim,nstp,3)
	real*4 ch(ndim,ndim,nstp,3),cvi(ndim,ndim,nstp,3)
        common/option/restart,dotsa,nstep,unit
        common/nda/ndon,nacc,idon,iacc,nb,iold(ndmax),jold(ndmax)
        common/iounit/iin,iout,istart(40),iungf,iunhl		
	dimension prob(ndim),dmtr(ndim,ndim),dmti(ndim,ndim)
	dimension temp(ndim,ndim),ovlp(ndim,ndim)
	
	transa='N'
	transb='N'
	twopi=8.d0*datan(1.d0)
        if(mod(ndim,2).eq.0) then
           nhaf=ndim*0.5
        else
           nhaf=(ndim+1)*0.5
        endif
	rectify=.false.
	gear=.false.
	ruku=.true.
	
        call init(idon,iacc,ndim,phi,psi0,psi1)
        if(.not.restart) then
           time0=0
           do i = 1 , ndim
              psi(i)=psi0(i)
           enddo
        else
           if(gear) call readns(ndim,time0,psi)
           if(ruku) call readrkns(ndim,time0,psi)
        endif

	write(iout,*)
	write(iout,*) ' ------------------------------------------ '
	write(iout,*) ' starting dynamics '
	write(iout,*)

	nwrit=nstep/(nstp-1)
	istep=1

 20     if(istep.gt.nstep) go to 30
 
 	time=real(istep)*tau+time0
 	itt=(time-tmd(1))/(tmd(2)-tmd(1))+1
 	dtt=time-tmd(itt)
 	         	      
 	if(gear) then         	      
	   call predictor(tau,ndim,psi)
	   call corrector(tau,ndim,psi,time,nstp,tmd,vtr,vti,cvr,cvi) 
	endif
	if(ruku) then
	   call rk4(tau,ndim,psi,time,nstp,tmd,vtr,vti,cvr,cvi)
	endif
	
  	istep=istep+1
  	
	if(mod(istep,nwrit).eq.0) then
	
           sum=0.0
           tot=(0.0,0.0)
           do i = 1 , ndim
              do j = 1 , ndim
                 sij=((cs(i,j,itt,3)*dtt+cs(i,j,itt,2))*dtt
     :                +cs(i,j,itt,1))*dtt+st(i,j,itt)
                 hij=((ch(i,j,itt,3)*dtt+ch(i,j,itt,2))*dtt
     :                +ch(i,j,itt,1))*dtt+ht(i,j,itt)
                 tot=tot+psi(j)*conjg(psi(i))*hij
                 ovlp(i,j)=sij
              enddo
           enddo
           
           do i = 1 , ndim
              do j = 1 , ndim
                 dmtr(i,j)= real(psi(i)*conjg(psi(j)))
c                 dmti(i,j)=aimag(psi(i)*conjg(psi(j)))
              enddo
           enddo
        
c          call mprod(ndim,ovlp,dmtr,temp)
c          call strassen(ndim,nhaf,ovlp,dmtr,temp)
           call dgemm(transa,transb,ndim,ndim,ndim,1.d0,
     *                ovlp,ndim,dmtr,ndim,0.d0,temp,ndim)         

           do i = 1 , ndim
              prob(i)=temp(i,i)
              sum=sum+prob(i)
           enddo
c           if(dabs(sum-1.0).gt.0.1) stop ' Total population > 1 !'
           if(dabs(sum-1.0).gt.0.05) then
              if(rectify) then
                 do i = 1 , ndim
                    psi(i)=psi(i)/dsqrt(sum)
                 enddo
              endif
           endif

	   write(iout,'(a5,f20.5,3x,2f10.5)')'time=',time/1.51653,
     :                                       sum,real(tot)
           write(iout,'(5f12.10)') (prob(i),i=1,ndim)
           
           write(11,'(f20.10,2f20.15)')
     :                 time/1.51653,prob(idon),prob(iacc)
           write(44,'(f20.10,f20.15)')
     :                 time/1.51653,real(tot)
           write(46,'(f20.10,f20.15)') 
     :                 time/1.51653,sum-prob(idon)-prob(iacc)

        endif

        go to 20
        
 30     if(gear) call savens(ndim,time)
        if(ruku) call saverkns(ndim,time,psi)
         
 	return
	end
	
	subroutine init(idon,iacc,n,phi,psi0,psi1) 
	implicit real*8 (a-h,o-z)
	complex*8 phi(n,n),psi0(n),psi1(n)
		
	do j = 1 , n
	   do i = 1 , n
	      phi(i,j)=(0.0,0.0)
	      if(i.eq.j) phi(i,j)=(1.0,0.0)
	   enddo
	enddo
	do i = 1 , n
	   psi0(i)=phi(idon,i)
	   psi1(i)=phi(iacc,i)
	enddo
	
	return
	end

	