	subroutine dyntsa(tau,nstp,tmd,st2,ht2)

c*********************************************************************
c perform 2sa electronic dynamics	
c*********************************************************************	   	
   	
	implicit real*8 (a-h,o-z)
	parameter (nsmp=10,ndmax=500)
	logical restart,dotsa,gear,ruku
	real*4 tmd(nsmp),st2(2,2,nsmp),ht2(2,2,nsmp),vt2(2,2,nsmp)
	complex*8 psi(2),psi0(2),psi1(2),phi(2,2)
	complex*8 tot
	real*4 cs2(2,2,nstp,3),ch2(2,2,nstp,3),cv2(2,2,nstp,3)
	real*4 cvi(2,2,nstp,3),vti(2,2,nstp)
	real*4 sij(nstp),hij(nstp),cofij(nstp,3),vij(nstp),tij(nstp)
        common/option/restart,dotsa,nstep,unit
        common/nda/ndon,nacc,idon,iacc,nb,iold(ndmax),jold(ndmax)
        common/iounit/iin,iout,istart(40),iungf,iunhl		
	dimension prob(2),dmtr(2,2),hamt(2,2),tinv(2,2)
	dimension temp(2,2),ovlp(2,2),temq(2,2)
	
	twopi=8.d0*datan(1.d0)
	gear=.false.
	ruku=.true.
	
	do i = 1 , 2
	   do j = 1 , 2
	      do k = 1 , nstp
	         vti(i,j,k)=0.0
	         do l = 1 , 3
	            cvi(i,j,k,l)=0.0
	         enddo
	      enddo
	   enddo
	enddo

	do n = 1 , nstp
	   do i = 1 , 2
	      do j = 1 , 2
	         temp(i,j)=st2(i,j,n)
	         hamt(i,j)=ht2(i,j,n)
	         tinv(i,j)=0.0
	      enddo
	      tinv(i,i)=1.0
	   enddo
	   call ludcmp(temp,2,2,indx,d)
	   do j = 1 , 2
	      call lubksb(temp,2,2,indx,tinv(1,j))
	   enddo
	   call mprod(2,tinv,hamt,temp)
	   do i = 1 , 2
	      do j = 1 , 2
	         vt2(i,j,n)=temp(i,j)
	      enddo
	   enddo
        enddo
        
        do k = 1 , nstp
           tij(k)=tmd(k)
        enddo
	do i = 1 , 2
	   do j = 1 , 2
	      do k = 1 , nstp
	         hij(k)=ht2(i,j,k)
	         sij(k)=st2(i,j,k)
	         vij(k)=vt2(i,j,k)
	      enddo
	      call icsccu(tij,hij,nstp,cofij,nstp-1,ier)
	      do k = 1 , nstp
	         do l = 1 , 3
	            ch2(i,j,k,l)=cofij(k,l)
	         enddo
	      enddo
	      call icsccu(tij,sij,nstp,cofij,nstp-1,ier)
	      do k = 1 , nstp
	         do l = 1 , 3
	            cs2(i,j,k,l)=cofij(k,l)
	         enddo
	      enddo
	      call icsccu(tij,vij,nstp,cofij,nstp-1,ier)
	      do k = 1 , nstp
	         do l = 1 , 3
	            cv2(i,j,k,l)=cofij(k,l)
	         enddo
	      enddo
	   enddo
	enddo
	
        call init(1,2,2,phi,psi0,psi1)
        if(.not.restart) then
           time0=0
           do i = 1 , 2
              psi(i)=psi0(i)
           enddo
        else
           if(gear) call read2s(2,time0,psi)
           if(ruku) call readrk2s(2,time0,psi)
        endif

	write(iout,*)
	write(iout,*) ' ------------------------------------------ '
	write(iout,*) ' starting 2S dynamics '
	write(iout,*)

	nwrit=nstep/(nstp-1)
	istep=1

 20     if(istep.gt.nstep) go to 30
 
 	time=real(istep)*tau+time0
 	itt=(time-tmd(1))/(tmd(2)-tmd(1))+1
 	dtt=time-tmd(itt)
 	         	      
        if(gear) then
	   call predictor2s(tau,2,psi)
	   call corrector2s(tau,2,psi,time,nstp,tmd,vt2,vti,cv2,cvi) 
	endif
	if(ruku) then
	   call rk2s4(tau,2,psi,time,nstp,tmd,vt2,vti,cv2,cvi)
	endif
	
  	istep=istep+1
  	
	if(mod(istep,nwrit).eq.0) then
	
           sum=0.0
           tot=(0.0,0.0)
           do i = 1 , 2
              do j = 1 , 2
                 sab=((cs2(i,j,itt,3)*dtt+cs2(i,j,itt,2))*dtt
     :                +cs2(i,j,itt,1))*dtt+st2(i,j,itt)
                 hab=((ch2(i,j,itt,3)*dtt+ch2(i,j,itt,2))*dtt
     :                +ch2(i,j,itt,1))*dtt+ht2(i,j,itt)
                 tot=tot+psi(j)*conjg(psi(i))*hab
                 ovlp(i,j)=sab
                 temq(i,j)=hab
              enddo
           enddo
           
           do i = 1 , 2
              do j = 1 , 2
                 dmtr(i,j)= real(psi(i)*conjg(psi(j)))
              enddo
           enddo
        
           call mprod(2,ovlp,dmtr,temp)

           do i = 1 , 2
              prob(i)=temp(i,i)
              sum=sum+prob(i)
           enddo
           if(dabs(sum-1.0).gt.0.1) stop ' Total population > 1 !'

	   write(iout,'(a5,f20.5,3x,2f10.5)')'time=',time/1.51653,
     :                                       sum,real(tot)
           write(iout,'(5f12.10)') (prob(i),i=1,2)
           
           write(13,'(f20.10,2f20.15)')
     :                 time/1.51653,prob(1),prob(2)

        endif

        go to 20
        
 30     if(gear) call save2s(2,time)
        if(ruku) call saverk2s(2,time,psi)
         
 	return
	end
