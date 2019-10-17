	subroutine outham(tau,ndim,nstep,ncut,hami,hamf)
	implicit real*8 (a-h,o-z)
	dimension hami(ndim,ndim,nstep)
	dimension hamf(ndim,ndim,0:ncut)
	character*8 filename1,filename2

	twopi=8.d0*datan(1.d0)
	period=nstep*tau
	comega=0.5*twopi/period
	
	filename1(1:6)='h/hami'
	filename2(1:6)='h/hamf'
	
	do i = 1 , ndim
	   if(i.lt.10) then
	      filename1(7:7)='0'
	      filename1(8:8)=char(48+i)
	      filename2(7:7)='0'
	      filename2(8:8)=char(48+i)
	   elseif(i.lt.100) then
	      nn1=mod(i,10)
	      nn2=int(i/10)
	      filename1(7:7)=char(48+nn2)
	      filename1(8:8)=char(48+nn1)
	      filename2(7:7)=char(48+nn2)
	      filename2(8:8)=char(48+nn1)
	   endif	
	   open(31,file=filename1)
	   do nt = 1 , nstep
	      write(31,'(500f10.5)')tau*nt/1.51653,
     :                              (hami(i,j,nt),j=1,ndim)
	   enddo 
	   close(31)
	   open(31,file=filename2)
	   do nt = 0 , ncut
	      write(31,'(500f10.5)')comega*nt*1.51653,
     :                              (hamf(i,j,nt),j=1,ndim)
	   enddo 
	   close(31)
	enddo
	
	return
	end
   	
   	subroutine output(n,hami,tau)
	implicit real*8 (a-h,o-z)
	parameter(neigd=2000)
	dimension hami(n,n)
        common /eigv/ omcm(neigd),zr(neigd,neigd),zi(neigd,neigd),
     :                ne,nblw
	
	write(*,*) ' hamiltonian '
	do i = 1 , n
	   write(*,'(100f10.5)') (hami(i,j),j=1,n)
	enddo
	write(*,*)
	write(*,*) ' eigenenergies '
	write(*,'(100f10.5)') (omcm(i),i=1,n)
	write(*,*)
	
	tda=omcm(2)-omcm(1)
	
	write(*,*) 
	write(*,*) ' energy split '
	write(*,'(f30.20)') tda
	write(*,*)
	
	planck=8.0*atan(1.0)
	span=planck/tda

c... time unit in hbar/eV=1.055/1.6*femtosecond
	
	write(*,*)
	write(*,*) ' Period (fs) '
	write(*,'(f20.1)') span
	write(*,*)

	write(*,*)
	write(*,*) ' Estimated number of steps needed '
	write(*,'(f20.1)') span/tau
	write(*,*)
	
	do i = 1 , n
	   write(*,*) ' eigenvector(',omcm(i),')'
	   write(*,'(100f10.5)') (zr(j,i),j=1,n)
	   write(*,'(100f10.5)') (zi(j,i),j=1,n)
	enddo
	
	return
        end