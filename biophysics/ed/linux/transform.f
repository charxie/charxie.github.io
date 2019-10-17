c<<<transform the hamiltonian matrix>>>

	subroutine trans1(ibeg,iend,nstp,ndim,h)

	implicit real*8 (a-h,o-z)
	parameter(neigd=248)
	real*4 h(ndim,ndim,nstp)
	common /eigv/ omcm(neigd),zr(neigd,neigd),zi(neigd,neigd),
     :                ne,nblw
        common /io/ iin,iout
	dimension hsub(2,2)
	dimension ssub(2,2)
	dimension tmat(ndim,ndim)
	dimension hnew(ndim,ndim)
	
	do nt = 1 , nstp
	
	   do i = 1 , 2
	      do j = 1 , 2
	         ssub(i,j)=0.0
	      enddo
	      ssub(i,i)=1.0
	   enddo
	   hsub(1,1)=h(ibeg,ibeg,nt)
	   hsub(1,2)=h(ibeg,iend,nt)
	   hsub(2,1)=h(iend,ibeg,nt)
	   hsub(2,2)=h(iend,iend,nt)

c	   write(iout,*)
c	   write(iout,*)'subblock of hamiltonian matrix'
c	   write(iout,'(2f8.4)') ((hsub(i,j),j=1,2),i=1,2)
c	   write(iout,*)

	   call eigen(-1000.0,1000.0,.true.,2,hsub,ssub)

c	   write(iout,'(2f8.4)') (omcm(i),i=1,2)
c	   write(iout,*)
c	   write(iout,'(2f8.4)') ((zr(i,j),i=1,2),j=1,2)
c	   write(iout,*)

	   do i = 1 , ndim
	      do j = 1 , ndim
	         hnew(i,j)=h(i,j,nt)
	      enddo
	   enddo
	   do i = 1 , ndim
	      hnew(ibeg,i)=zr(1,1)*h(ibeg,i,nt)+zr(2,1)*h(iend,i,nt)
	      hnew(i,ibeg)=hnew(ibeg,i)
	      hnew(iend,i)=zr(1,2)*h(ibeg,i,nt)+zr(2,2)*h(iend,i,nt)
	      hnew(i,iend)=hnew(iend,i)
	   enddo
	   hnew(ibeg,ibeg)=zr(1,1)*zr(1,1)*h(ibeg,ibeg,nt)+
     :	    	           zr(1,1)*zr(2,1)*h(ibeg,iend,nt)+
     :    		   zr(2,1)*zr(1,1)*h(iend,ibeg,nt)+
     :                     zr(2,1)*zr(2,1)*h(iend,iend,nt)
	   hnew(ibeg,iend)=zr(1,1)*zr(1,2)*h(ibeg,ibeg,nt)+
     :	    	           zr(1,1)*zr(2,2)*h(ibeg,iend,nt)+
     :    		   zr(2,1)*zr(1,2)*h(iend,ibeg,nt)+
     :                     zr(2,1)*zr(2,2)*h(iend,iend,nt)
	   hnew(iend,ibeg)=zr(1,2)*zr(1,1)*h(ibeg,ibeg,nt)+
     :	    	           zr(1,2)*zr(2,1)*h(ibeg,iend,nt)+
     :    		   zr(2,2)*zr(1,1)*h(iend,ibeg,nt)+
     :                     zr(2,2)*zr(2,1)*h(iend,iend,nt)
	   hnew(iend,iend)=zr(1,2)*zr(1,2)*h(ibeg,ibeg,nt)+
     :	    	           zr(1,2)*zr(2,2)*h(ibeg,iend,nt)+
     :    		   zr(2,2)*zr(1,2)*h(iend,ibeg,nt)+
     :                     zr(2,2)*zr(2,2)*h(iend,iend,nt)

	   do i = 1 , ndim
	      do j = 1 , ndim
	         h(i,j,nt)=hnew(i,j)
	      enddo
	   enddo

c	   do i = 1 , ndim
c	      write(iout,'(100f8.4)') (hnew(i,j),j=1,ndim)
c	   enddo
c	   write(iout,*)

	enddo
	
	return
	end
