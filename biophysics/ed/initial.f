c define the localized electronic states
c phi(idon,n)  ...... donor state
c phi(i,n)  ...... bridge states
c phi(iacc,n)  ...... acceptor state
c psi0(n)   ...... initial state
c psi1(n)   ...... final state

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

        subroutine initeigen(idon,iacc,ndim,psi0,emat)
        implicit real*8 (a-h,o-z)
        dimension emat(ndim,ndim),einv(ndim,ndim)
        complex*8 psi0(ndim)
        
        call matinv(ndim,emat,einv)
        do i = 1 , ndim
           psi0(i)=einv(i,idon)*(1.0,0.0)
        enddo
        
        
        return
        end