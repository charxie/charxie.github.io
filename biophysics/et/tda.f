c calculate the static Green's function

    	subroutine caltda(n,h,e,tda,tdd,taa)
    	
    	implicit real*8 (a-h,o-z)
	dimension h(n,n)
	dimension hbr(n-2,n-2),gbr(n-2,n-2),sbr(n-2,n-2)
	dimension hmes(n-2,n-2),indx(n-2)
	
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
     	
     	write(*,*) 
	do i = 1 , n-2
	   write(*,'(100f10.5)') (hbr(i,j),j=1,n-2)
	enddo     	
     	
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

     	write(*,*)
     	write(*,*) ' Green function (', e, ')'
     	write(*,*)
    	do k = 1 , n-2
     	   write(*,'(100f10.5)') (gbr(k,l),l=1,n-2)
     	enddo 

	tda=h(1,n)
	tdd=h(1,1)
	taa=h(n,n)

	do i = 1 , n-2
	   do j = 1 , n-2
	      tda=tda+h(1,i+1)*gbr(i,j)*h(j+1,n)
	      tdd=tdd+h(1,i+1)*gbr(i,j)*h(j+1,1)
	      taa=taa+h(n,i+1)*gbr(i,j)*h(j+1,n)
	   enddo
	enddo
     
        write(*,'(a)') 'two-state reduction'
     	write(*,'(2f10.5)') tdd,tda
     	write(*,'(2f10.5)') tda,taa

     
        return
        end 