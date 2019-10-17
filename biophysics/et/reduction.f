c Reduce the Hamiltonian matrix to effective N-states Hamiltonian(N=3)

	subroutine hred(n,h,e,nred,heff)
	implicit real*8 (a-h,o-z)
 	dimension h(n,n),hbr(n-2,n-2)
	dimension hqq(n-3,n-3),gqq(n-3,n-3),sqq(n-3,n-3)
	dimension hmes(n-3,n-3),indx(n-3)
	dimension hpq(3,n-3),hqp(n-3,3)
	dimension hbrpq(3,n-2),hbrqp(n-2,3)
	dimension hbrpp(2,2),hpp(3,3)
	dimension heff(3,3)
	logical rite
	
	k2=int((n+1)/2)
	rite=.true.
	
	do i = 1 , n-3
	   do j = 1 , n-3
	      sqq(i,j)=0.0
	      if(i.eq.j) sqq(i,j)=1.0
	   enddo
	enddo
	
     	do i = 1 , n-2
     	   do j = 1 , n-2
              hbr(i,j)=h(i+1,j+1)
     	   enddo
     	enddo
     	
     	hbrpp(1,1)=h(1,1)
     	hbrpp(1,2)=h(1,n)
     	hbrpp(2,1)=h(n,1)
     	hbrpp(2,2)=h(n,n)
     	
     	do i = 1 , n-2
     	   hbrpq(1,i)=h(1 ,i+1)
     	   hbrpq(2,i)=h(k2,i+1)
     	   hbrpq(3,i)=h(n ,i+1)
     	   hbrqp(i,1)=h(i+1,1 )
     	   hbrqp(i,2)=h(i+1,k2)
     	   hbrqp(i,3)=h(i+1,n )
     	enddo
     	
     	if(rite) then
     	  write(*,'(100f6.3)') (hbrpq(1,i),i=1,n-2)
     	  write(*,'(100f6.3)') (hbrpq(2,i),i=1,n-2)
	  write(*,'(100f6.3)') (hbrpq(3,i),i=1,n-2)
     	  write(*,'(100f6.3)') (hbrqp(i,1),i=1,n-2)
     	  write(*,'(100f6.3)') (hbrqp(i,2),i=1,n-2)
     	  write(*,'(100f6.3)') (hbrqp(i,3),i=1,n-2)
     	endif
     	
     	kb2=k2-1
     	do i = 1 , n-2
     	   do j = 1 , n-2
     	      if(i.lt.kb2.and.j.lt.kb2) then
                 hqq(i,j)=hbr(i,j)
              endif
     	      if(i.gt.kb2.and.j.lt.kb2) then
                 hqq(i-1,j)=hbr(i,j)
              endif
     	      if(i.lt.kb2.and.j.gt.kb2) then
                 hqq(i,j-1)=hbr(i,j)
              endif
     	      if(i.gt.kb2.and.j.gt.kb2) then
                 hqq(i-1,j-1)=hbr(i,j)
              endif
     	   enddo
     	enddo
      	
	do i = 1 , n-2
	   if(i.lt.kb2) then
	      hpq(1,i)=hbrpq(1,i)
	      hpq(2,i)=hbrpq(2,i)
	      hpq(3,i)=hbrpq(3,i)
	      hqp(i,1)=hbrqp(i,1)
	      hqp(i,2)=hbrqp(i,2)
	      hqp(i,3)=hbrqp(i,3)
	   endif
	   if(i.gt.kb2) then
	      hpq(1,i-1)=hbrpq(1,i)
	      hpq(2,i-1)=hbrpq(2,i)
	      hpq(3,i-1)=hbrpq(3,i)
	      hqp(i-1,1)=hbrqp(i,1)
	      hqp(i-1,2)=hbrqp(i,2)
	      hqp(i-1,3)=hbrqp(i,3)
	   endif
	enddo
	
	hpp(1,1)=h(1 ,1 )
	hpp(1,2)=h(1 ,k2)
	hpp(1,3)=h(1 ,n )
	hpp(2,1)=h(k2,1 )
	hpp(2,2)=h(k2,k2)
	hpp(2,3)=h(k2,n )
	hpp(3,1)=h(n ,1 )
	hpp(3,2)=h(n ,k2)
	hpp(3,3)=h(n ,n )
		
        if(rite) then		
     	  write(*,*) 
	  do i = 1 , n
	    write(*,'(100f6.3)') (h(i,j),j=1,n)
	  enddo     	
	  write(*,*)
	  do i = 1 , n-3
	    write(*,'(100f6.3)') (hqq(i,j),j=1,n-3)
	  enddo
	  write(*,*)     	
     	  write(*,'(100f6.3)') (hpq(1,i),i=1,n-3)
     	  write(*,'(100f6.3)') (hpq(2,i),i=1,n-3)
	  write(*,'(100f6.3)') (hpq(3,i),i=1,n-3)
     	  write(*,'(100f6.3)') (hqp(i,1),i=1,n-3)
     	  write(*,'(100f6.3)') (hqp(i,2),i=1,n-3)
     	  write(*,'(100f6.3)') (hqp(i,3),i=1,n-3)
     	endif
     	
	do i = 1 , n-3
	   do j = 1 , n-3
	      hmes(i,j)=e*sqq(i,j)-hqq(i,j)
	      gqq(i,j)=0.d0
	   enddo
	   gqq(i,i)=1.0
	enddo
	
	call ludcmp(hmes,n-3,n-3,indx,d)

	do j = 1 , n-3
	   call lubksb(hmes,n-3,n-3,indx,gqq(1,j))
	enddo
	
	if(rite) then
     	  write(*,*)
     	  write(*,*) ' Green function (', e, ')'
     	  write(*,*)
    	  do k = 1 , n-3
     	    write(*,'(100f10.5)') (gqq(k,l),l=1,n-3)
     	  enddo      	
     	endif

        do ipp = 1 , 3
           do jpp = 1 , 3
	      hef12=hpp(ipp,jpp)
	      do i = 1 , n-3
	         do j = 1 , n-3
	           hef12=hef12+hpq(ipp,i)*gqq(i,j)*hqp(j,jpp)
	         enddo
	      enddo
	      heff(ipp,jpp)=hef12
	   enddo
	enddo
	
	write(*,*) 'three-state reduction'
	write(*,'(3f10.5)') ((heff(i,j),j=1,3),i=1,3)
	
        return
        end        
                