c calculate the static properties for a model system

	subroutine sta(tau,n,ED,EA,EB,VD1,VBB,VNA)
	implicit real*8 (a-h,o-z)
	parameter (neigd=2000)
	dimension esit(n),hopi(n,n),hami(n,n),hnew(n,n)
	dimension dm(n,n),gf(n,n)
	common /eigv/ omcm(neigd),zr(neigd,neigd),zi(neigd,neigd),
     :                ne,nblw
	data etun/0.0/
	
	call hamiltonian0(ED,EA,EB,VD1,VBB,VNA,n,esit,hopi,hami)
	call eigen(n,hami)
	call output(n,hami,tau)
c	call demat(n,dm)
c	do i = 1 , 100	   
c	   etun=-1.0+0.02*i
c	   call caltda(n,hami,etun,tda)	
c	enddo
c
c>>> calculate the Green function
c	
c	open(34,file='grnf',status='unknown')
c	e=-1.0
c	i=0
c1000	e=e+0.001		
c	i=i+1
c	call green(n,hami,e,gf)
c	write(34,'(36d10.3)') e,gf(4,4)
c	if(i.le.10000) go to 1000
c	call eigenbridge(n,hami,hnew)

	return
	end