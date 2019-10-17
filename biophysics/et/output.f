   	
   	subroutine output(n,hami,tau)
	implicit real*8 (a-h,o-z)
	parameter(neigd=50)
	dimension hami(n,n)
        common /eigv/ omcm(neigd),zr(neigd,neigd),zi(neigd,neigd),
     :                ne,nblw
	
	write(*,*) ' hamiltonian '
	do i = 1 , n
	   write(*,'(10f8.4)') (hami(i,j),j=1,n)
	enddo
	write(*,*)
	write(*,*) ' eigenenergies '
	write(*,'(100f10.5)') (omcm(i),i=1,n)
	write(*,*)
	
	tda=omcm(2)-omcm(1)
	
	write(*,*) 
	write(*,*) ' energy split (eV) '
	write(*,'(f30.20)') tda
	write(*,*)
	
	planck=8.0*atan(1.0)
	span=planck/tda
c... time unit in hbar/eV=1.055/1.6*femtosecond
	
	write(*,*)
	write(*,*) ' Period (0.6594 fs) '
	write(*,'(f20.1)') span
	write(*,*)

	write(*,*)
	write(*,*) ' Estimated number of steps needed '
	write(*,'(f20.1)') span/tau
	write(*,*)
	
	do i = 1 , 3
	   write(*,*) ' eigenvector(',omcm(i),')'
	   write(*,'(5f10.5)') (zr(j,i),j=1,n)
c	   write(*,'(5f10.5)') (zi(j,i),j=1,n)
	enddo
	
	return
        end