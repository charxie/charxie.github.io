c Reduce the Hamiltonian matrix to effective N-states Hamiltonian(N=3)
c must be idon < imid < iacc

	subroutine triple(n,h0,e,idon,imid,iacc,heff)
	implicit real*8 (a-h,o-z)
 	dimension h0(n,n),heff(3,3)
 	dimension h1(n-1,n-1),c1(n-1)
 	dimension h2(n-2,n-2),c2(n-2)
 	dimension h3(n-3,n-3),c3(n-3)
 	dimension hmes(n-3,n-3),indx(n-3)
 	dimension gbr(n-3,n-3),cbr(3,n-3)
	logical rite
	
	rite=.true.
	
	call deleteh(idon,n  ,h0,h1)
	call deleteh(imid-1,n-1,h1,h2)
	call deleteh(iacc-2,n-2,h2,h3)
	call deletec(idon,imid,iacc,n,h0,cbr)
	
	if(rite) then
	  do i = 1 , n
	    write(6,'(10f8.4)') (h0(i,j),j=1,n)
	  enddo
	  write(6,*)
	  do i = 1 , n-3
	    write(6,'(10f8.4)') (h3(i,j),j=1,n-3)
	  enddo
	  write(6,*)
	  write(6,'(10f8.4)') (cbr(1,i),i=1,n-3)
	  write(6,'(10f8.4)') (cbr(2,i),i=1,n-3)
	  write(6,'(10f8.4)') (cbr(3,i),i=1,n-3)
	  write(6,*) 
	endif
	
	
	heff(1,1)=h0(idon,idon)
	heff(1,2)=h0(idon,imid)
	heff(1,3)=h0(idon,iacc)
	heff(2,1)=h0(imid,idon)
	heff(2,2)=h0(imid,imid)
	heff(2,3)=h0(imid,iacc)
	heff(3,1)=h0(iacc,idon)
	heff(3,2)=h0(iacc,imid)
	heff(3,3)=h0(iacc,iacc)
	
	if(rite) then
	  do i = 1 , 3
	    write(6,'(3f8.4)') (heff(i,j),j=1,3)
	  enddo
	  write(6,*)
	endif
		
     	
	do i = 1 , n-3
	   do j = 1 , n-3
	      hmes(j,i)=-h3(j,i)
	      gbr(j,i)=0.d0
	   enddo
	   hmes(i,i)=hmes(i,i)+e
	   gbr(i,i)=1.d0
	enddo
	
	call ludcmp(hmes,n-3,n-3,indx,d)

	do j = 1 , n-3
	   call lubksb(hmes,n-3,n-3,indx,gbr(1,j))
	enddo

        if(rite) then
	  do i = 1 , n-3
	    write(6,'(10f8.4)') (gbr(i,j),j=1,n-3)
	  enddo
	  write(6,*)
	endif
	
        do ipp = 1 , 3
           do jpp = 1 , 3
              hef12=0.0
	      do i = 1 , n-3
	         do j = 1 , n-3
	           hef12=hef12+cbr(ipp,i)*gbr(i,j)*cbr(jpp,j)
	         enddo
	      enddo
	      heff(ipp,jpp)=hef12+heff(ipp,jpp)
	   enddo
	enddo
	
	write(*,*) 'three-state reduction'
	write(*,'(3f10.5)') ((heff(i,j),j=1,3),i=1,3)
	
        return
        end        
                