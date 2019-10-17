c Time-dependent electronic dynamics in the two-state approximation
c
	subroutine tsa(tau,nstp,tmd,h2s,pdir,pinv)
	implicit real*8 (a-h,o-z)
	real*4 tmd(nstp),h2s(2,2,nstp),cof2s(2,2,nstp,3)
	real*4 hij(nstp),cofij(nstp,3)
	complex*8 psi0(2),psi(2),psi1(2),phi(2,2),project(2),tot
	dimension prob(2)
	dimension pdir(2,2,nstp),pinv(2,2,nstp)
	logical restart,dotsa,gear,ruku
        common /ida/ idon,iacc,i0,j0,restart,time0,nstep,edon,eacc
        common /io/ iin,iout
        
        dotsa=.true.
        gear=.true.
        ruku=.false.

	do i = 1 , 2
	   do j = 1 , 2
	      do k = 1 , nstp
	         hij(k)=h2s(i,j,k)
	      enddo
	      call icsccu(tmd,hij,nstp,cofij,nstp-1,ier)
	      do k = 1 , nstp
	         do l = 1 , 3
	            cof2s(i,j,k,l)=cofij(k,l)
	         enddo
	      enddo
	   enddo
	enddo
	call init(1,2,2,phi,psi0,psi1)
        if(.not.restart) then
           time0=0.0
           do i = 1 , 2
              psi(i)=psi0(i)
           enddo
        else
           if(gear) call readrst(2,time0,psi,dotsa)
           if(ruku) call readrkns(2,time0,psi0,dotsa)
        endif
	nwrit=nstep/(nstp-1)
        istep=1

 20     if(istep.gt.nstep) go to 30
 
	time=istep*tau+time0
	itt=(time-tmd(1))/(tmd(2)-tmd(1))+1
	dtt=time-tmd(itt)
c        if(mod(istep,nwrit).eq.0) then
c           h12=((cof2s(1,2,itt,3)*dtt+cof2s(1,2,itt,2))*dtt
c     :          +cof2s(1,2,itt,1))*dtt+h2s(1,2,itt)
c           write(12,'(50f20.10)') time/1.51653,h12
c        endif
        
        if(gear) then
	   call predictor(tau,2,psi)
	   call corrector(tau,2,psi,time,nstp,tmd,h2s,cof2s) 
	endif
	if(ruku) then
	   call rk4(tau,2,psi,time,nstp,tmd,h2s,cof2s)
	endif
 
  	istep=istep+1
        
	if(mod(istep,100).eq.0) then

	   call proj(2,psi,phi,project)
	   
           sum=0.0
           tot=(0.0,0.0)
        
           do i = 1 , 2
              prob(i)=(real(project(i)))**2+(aimag(project(i)))**2
              sum=sum+prob(i)
           enddo
           if(dabs(sum-1.0).gt.0.1) stop ' Total population > 1 !'

           do i = 1 , 2
              do j = 1 , 2
                 ham=((cof2s(i,j,itt,3)*dtt+cof2s(i,j,itt,2))*dtt
     :                +cof2s(i,j,itt,1))*dtt+h2s(i,j,itt)
                 tot=tot+
     :               project(j)*conjg(project(i))*ham
              enddo
           enddo
 		
	   write(iout,'(a5,f20.5,3x,2f10.5)')'time=',time/1.51653,
     :                                       sum,real(tot)
           write(iout,'(5f12.10)') (prob(i),i=1,2)
           
           write(11,'(f20.10,20f20.15)')
     :          time/1.51653,prob(2),
     :          -pdir(1,2,itt),
     :          -(pinv(1,1,itt)*project(1)+pinv(2,1,itt)*project(2)),
     :          -pdir(2,2,itt),
     :          -(pinv(1,2,itt)*project(1)+pinv(2,2,itt)*project(2))
           write(44,'(10f20.10)')
     :                 time/1.51653,real(tot)

        endif

        go to 20

 30     if(gear) call save(2,time,dotsa)
        if(ruku) call saverkns(2,time,psi,dotsa)
        
	close(11)
	close(44) 
        return
        end