c rebuild the Hamiltonian in the space spanned by the localized
c states of donor and acceptor, and eigenstates of the bridge

    	subroutine eigenbridge(n,h,hnew)
    	
    	implicit real*8 (a-h,o-z)
 	parameter(nmax=10,nm2=nmax-2,neigd=2000)
	dimension h(n,n),hnew(n,n)
	dimension hbr(nm2,nm2),sbr(nm2,nm2)
        common /eigv/ omcm(neigd),zr(neigd,neigd),zi(neigd,neigd),ne,nblw
	
	do i = 1 , n-2
	   do j = 1 , n-2
	      sbr(i,j)=0.0
	      if(i.eq.j) sbr(i,j)=1.0
	   enddo
	enddo
     	
     	do i = 1 , n-2
     	   do j = 1 , n-2
     	      hbr(i,j)=h(i+1,j+1)
     	   enddo
     	enddo
     	
        call eigen(n-2,hbr)

	write(*,*)
	write(*,*) ' eigenenergies '
	write(*,'(100f10.5)') (omcm(i),i=1,n-2)
	write(*,*)
        do i = 1 , n-2
	   write(*,*) ' eigenvector(',omcm(i),')'
	   write(*,'(100f10.5)') (zr(j,i),j=1,n-2)
	   write(*,'(100f10.5)') (zi(j,i),j=1,n-2)
	enddo        
        
        do i = 1 , n-2
           hnew(1,i+1)=0.0
           do j = 1 , n-2
 	      hnew(1,i+1)=hnew(1,i+1)+h(1,j+1)*zr(j,i)
 	   enddo
 	   hnew(i+1,n)=0.0
           do j = 1 , n-2
 	      hnew(i+1,n)=hnew(i+1,n)+h(n,j+1)*zr(j,i)
 	   enddo
        enddo
        hnew(1,1)=h(1,1)
        hnew(1,n)=h(1,n)
        hnew(n,n)=h(n,n)
        
        do i = 1 , n-2
           do j = 1 , n-2
              hnew(i+1,j+1)=0.0
           enddo
           hnew(i+1,i+1)=omcm(i)
        enddo

	do i = 1 , n
	   do j = 1 , i
	      hnew(i,j)=hnew(j,i)
	   enddo
	enddo
        
        write(*,*)
        write(*,*) ' New Hamiltonian '
        write(*,*)
	do i = 1 , n
	   write(*,'(100f10.5)') (hnew(i,j),j=1,n)
	enddo
         

	call eigen(n,hnew)

	write(*,*)
	write(*,*) ' eigenenergies '
	write(*,'(100f10.5)') (omcm(i),i=1,n)
	write(*,*)
c	do i = 1 , n
c	   write(*,*) ' eigenvector(',omcm(i),')'
c	   write(*,'(100f10.5)') (zr(j,i),j=1,n)
c	   write(*,'(100f10.5)') (zi(j,i),j=1,n)
c	enddo        
	
        
        return
        end        
        