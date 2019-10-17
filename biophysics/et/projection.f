c projection psi(t) onto localized states phi

	subroutine proj(n,psi,phi,project)
   	
	implicit real*8 (a-h,o-z)
	complex psi(n),phi(n,n),project(n)

	do i = 1 , n
   	   project(i)=(0.0,0.0)
   	   do j = 1 , n
   	      project(i)=project(i)+phi(i,j)*psi(j)
   	   enddo	
	enddo

	return
	end

	
c time-dependent density matrix

    	subroutine dent(n,project,rho,hami,flow)
    	
    	implicit real*8 (a-h,o-z)
    	complex project(n),rho(n,n),flow(n,n)
    	dimension hami(n,n)

	do i = 1 , n
	   do j = 1 , n
	      rho(i,j)=project(i)*conjg(project(j))
	   enddo
    	enddo
    	
    	do i = 1 , n
    	   do j = 1 , n
    	      flow(i,j)=hami(i,j)*rho(j,i)-rho(i,j)*hami(j,i)
    	      flow(i,j)=-(0.0,1.0)*flow(i,j)
    	   enddo
    	enddo
	
	return
	end
	
c calculate the static density matrix using the eigenvectors

    	subroutine demat(n,dm)
    	
    	implicit real*8 (a-h,o-z)
 	parameter(neigd=50)
	dimension dm(n,n)
        common /eigv/ omcm(neigd),zr(neigd,neigd),zi(neigd,neigd),
     :                ne,nblw
     
     	do k = 1 , n
     	   do l = 1 , n 
     	      dm(k,l)=0.d0
     	      do i = 1 , n
     	         dm(k,l)=dm(k,l)+zr(i,k)*zr(i,l)+zi(i,k)*zi(i,l)
     	      enddo
     	   enddo
     	enddo
     	
     	write(*,*)
     	write(*,*) ' density matrix '
     	write(*,*)
     	
     	do k = 1 , n	
     	   write(*,'(100f10.5)') (dm(k,l),l=1,n)
     	enddo 
     
        return
        end
        
c calculate the static Green's function

    	subroutine green(n,h,e,g)
    	
    	implicit real*8 (a-h,o-z)
 	parameter(nmax=10)
	dimension h(nmax,nmax),g(nmax,nmax),s(nmax,nmax)
	dimension hmes(nmax,nmax),indx(nmax)
	
	do i = 1 , n
	   do j = 1 , n
	      s(i,j)=0.0
	      if(i.eq.j) s(i,j)=1.0
	   enddo
	enddo
     	
	do i = 1 , n
	   do j = 1 , n
	      hmes(i,j)=e*s(i,j)-h(i,j)
	      g(i,j)=0.d0
	   enddo
	   g(i,i)=1.0
	enddo
	
	call ludcmp(hmes,n,n,indx,d)

	do j = 1 , n
	   call lubksb(hmes,n,n,indx,g(1,j))
	enddo

     	
     	write(*,*)
     	write(*,*) ' Green function (', e, ')'
     	write(*,*)
     	
     	do k = 1 , n	
     	   write(*,'(100f10.5)') (g(k,l),l=1,n)
     	enddo 
     
        return
        end        
        