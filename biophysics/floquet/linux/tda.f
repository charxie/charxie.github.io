c calculate the static Green's function

    	subroutine caltda(n,h,e,tda)
    	
    	implicit real*8 (a-h,o-z)
 	parameter(nmax=10,nm2=nmax-2)
	dimension h(n,n)
	dimension hbr(nm2,nm2),gbr(nm2,nm2),sbr(nm2,nm2)
	dimension hmes(nm2,nm2),indx(nm2)
	
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
     	
c     	write(*,*) 
c	do i = 1 , n-2
c	   write(*,'(100f10.5)') (hbr(i,j),j=1,n-2)
c	enddo     	
     	
	do i = 1 , n-2
	   do j = 1 , n-2
	      hmes(i,j)=e*sbr(i,j)-hbr(i,j)
	      gbr(i,j)=0.d0
	   enddo
	   gbr(i,i)=1.0
	enddo
	
	call ludcmp(hmes,n-2,n-2,indx,d)

	do j = 1 , n-2
	   call lubksb(hmes,n-2,n-2,indx,gbr(1,j))
	enddo

c     	write(*,*)
c     	write(*,*) ' Green function (', e, ')'
c     	write(*,*)
c    	do k = 1 , n-2
c     	   write(*,'(100f10.5)') (gbr(k,l),l=1,n-2)
c     	enddo 

	tda=0.0

	do i = 1 , n-2
	   do j = 1 , n-2
	      tda=tda+h(1,i+1)*gbr(i,j)*h(j+1,n)
	   enddo
	enddo
     
     	write(*,'(a5,1x,f10.5,1x,a4,1x,f10.5)')' TDA(',e,') = ',tda

     
        return
        end        
        