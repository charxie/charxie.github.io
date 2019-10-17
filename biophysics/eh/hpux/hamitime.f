
	subroutine hint(ndim,nstp,tmd,st,ht,ft,vtr,vti,cs,ch,cvr,cvi)

c*******************************************************************
c prepare interpolation coefficients from the time series
c for both the overlap and fock matrices
c and S^{-1}(t)H(t) (vtr), S^{-1}(t)K(t) (vti)
c*******************************************************************

	implicit real*8 (a-h,o-z)
	parameter (ndmax=500,nsmp=10)
	logical rite
	character transa*1,transb*1
	real*4 tmd(nsmp)
	real*4 st(ndmax,ndmax,nsmp),ht(ndmax,ndmax,nsmp)
	real*4 ft(ndmax,ndmax,nsmp)
	real*4 vtr(ndmax,ndmax,nsmp),vti(ndmax,ndmax,nsmp)
	real*4 cs(ndim,ndim,nstp,3),ch(ndim,ndim,nstp,3)
	real*4 cvr(ndim,ndim,nstp,3),cvi(ndim,ndim,nstp,3)
	real*4 tij(nstp),sij(nstp),vij(nstp),hij(nstp),csij(nstp,3),
     :         cvij(nstp,3),chij(nstp,3)
	dimension tinv(ndim,ndim),temp(ndim,ndim),indx(ndim)
	dimension hamt(ndim,ndim),fric(ndim,ndim)
        common/iounit/iin,iout,istart(40),iungf,iunhl		

        alpha=1.d0
        beta=0.d0
        lda=ndim
        ldb=ndim
        ldc=ndim
        transa='N'
        transb='N'
        	
	if(ndim.gt.ndmax) stop ' ndim > ndmax ! '
	if(nstp.gt.nsmp)  stop ' nstp > nsmp ! ' 
	write(iout,*) 'interpolating...'
	
	rite=.true.
        if(mod(ndim,2).eq.0) then
           nhaf=ndim*0.5
        else
           nhaf=(ndim+1)*0.5
        endif
	
	do n = 1 , nstp
	   do j = 1 , ndim
	      do i = 1 , ndim
	         temp(i,j)=st(i,j,n)
	         hamt(i,j)=ht(i,j,n)
	         fric(i,j)=ft(i,j,n)
	         tinv(i,j)=0.0
	      enddo
	      tinv(j,j)=1.0
	   enddo
	   call ludcmp(temp,ndim,ndim,indx,d)
	   do j = 1 , ndim
	      call lubksb(temp,ndim,ndim,indx,tinv(1,j))
	   enddo

c	   call mprod(ndim,tinv,hamt,temp)
c          call strassen(ndim,nhaf,tinv,hamt,temp)
           call dgemm(transa,transb,ndim,ndim,ndim,alpha,
     *                tinv,lda,hamt,ldb,beta,temp,ldc)         

	   do i = 1 , ndim
	      do j = 1 , ndim
	         vtr(i,j,n)=temp(i,j)
	      enddo
	   enddo

c	   call mprod(ndim,tinv,fric,temp)
c	   call strassen(ndim,nhaf,tinv,fric,temp)
           call dgemm(transa,transb,ndim,ndim,ndim,alpha,
     *                tinv,lda,fric,ldb,beta,temp,ldc)         

	   do i = 1 , ndim
	      do j = 1 , ndim
	         vti(i,j,n)=temp(i,j)
	      enddo
	   enddo
	enddo
	
	do k = 1 , nstp
	   tij(k)=tmd(k)
	enddo
	do i = 1 , ndim
	   do j = 1 , ndim

	      do k = 1 , nstp
	         sij(k)=st(i,j,k)
	         hij(k)=ht(i,j,k)
	         do l = 1 , 3
	            csij(k,l)=0.0
	            chij(k,l)=0.0
	         enddo
	      enddo	
              call icsccu(tij,sij,nstp,csij,nstp-1,ier)
              call icsccu(tij,hij,nstp,chij,nstp-1,ier)
	      do k = 1 , nstp
	         do l = 1 , 3
                    cs(i,j,k,l)=csij(k,l)
                    ch(i,j,k,l)=chij(k,l)
                 enddo
              enddo
              
	      do k = 1 , nstp
	         vij(k)=vtr(i,j,k)
	         do l = 1 , 3
	            cvij(k,l)=0.0
	         enddo
	      enddo	
              call icsccu(tij,vij,nstp,cvij,nstp-1,ier)
	      do k = 1 , nstp
	         do l = 1 , 3
                    cvr(i,j,k,l)=cvij(k,l)
                 enddo
              enddo
              
	      do k = 1 , nstp
	         vij(k)=vti(i,j,k)
	         do l = 1 , 3
	            cvij(k,l)=0.0
	         enddo
	      enddo	
              call icsccu(tij,vij,nstp,cvij,nstp-1,ier)
	      do k = 1 , nstp
	         do l = 1 , 3
                    cvi(i,j,k,l)=cvij(k,l)
                 enddo
              enddo
          
           if(i.eq.0.and.j.eq.0) then
              write(iout,*) i,j
              write(iout,'(10d10.4)')(sij(k),k=1,nstp)   
              write(iout,'(10d10.4)')((cs(i,j,k,l),k=1,nstp),l=1,3)
              write(iout,'(10d10.4)')(hij(k),k=1,nstp)   
              write(iout,'(10d10.4)')((ch(i,j,k,l),k=1,nstp),l=1,3)
              write(iout,'(10d10.4)')(vij(k),k=1,nstp)   
              write(iout,'(10d10.4)')((cvr(i,j,k,l),k=1,nstp),l=1,3)
              stop
           endif
              
	   enddo
	enddo
	
	if(rite) then
	   write(iout,*)'S'
	   do l = 1 , 3
	      write(iout,*)
	      write(iout,*) 'l=',l
	      cmax=0.0
	      imax=0
	      jmax=0
	      kmax=0
	      do k = 1 , nstp
	         do i = 1 , ndim
	            do j = 1 , ndim
	               if(abs(cs(i,j,k,l)).ge.cmax)then
	                  cmax=abs(cs(i,j,k,l))
	                  imax=i
	                  jmax=j
	                  kmax=k
	               endif
	            enddo
	         enddo
	      enddo
	      write(iout,*)imax,jmax,kmax,cmax
	   enddo
	   write(iout,*)'H'
	   do l = 1 , 3
	      write(iout,*)
	      write(iout,*) 'l=',l
	      cmax=0.0
	      imax=0
	      jmax=0
	      kmax=0
	      do k = 1 , nstp
	         do i = 1 , ndim
	            do j = 1 , ndim
	               if(abs(ch(i,j,k,l)).ge.cmax)then
	                  cmax=abs(ch(i,j,k,l))
	                  imax=i
	                  jmax=j
	                  kmax=k
	               endif
	            enddo
	         enddo
	      enddo
	      write(iout,*)imax,jmax,kmax,cmax
	   enddo
	endif
	
	return
	end