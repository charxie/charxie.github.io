	subroutine denmat(natom,ndim,n,st,ht)
	
c********************************************************************
c calculate the density matrix and the molecular orbital populations
c********************************************************************	
	
	implicit real*8 (a-h,o-z)
	parameter(nsmp=10,natmax=200,ndmax=500,neigd=500)
	character indt*10,ac*4,symbol*2
	real*4 dm,pop
	real*4 st(ndmax,ndmax,nsmp),ht(ndmax,ndmax,nsmp)
        common /eigv/ omcm(neigd),zr(neigd,neigd),zi(neigd,neigd),ne,nblw
	common /eigt/ dm(ndmax,ndmax,nsmp),pop(ndmax,nsmp)
        common/ato/ac(natmax),symbol(40)
        common/aoout/indt(1000,2),iatom(natmax),lorb(1000)
        common/atom1/ns(40),np(40),nd(40),nf(40),inatom,
     *  isym(natmax),nelem
        common/optv/rho,fermi,wstart,wend,
     *  lb,le,lbr,ler,iovlb,iovle,lcharg,neltr,
     *  distl,idistn,ieb,iee,isb,ise,ibb,ibe,ipope,jpope,irest 
        common/iounit/iin,iout,istart(40),iungf,iunhl            
        dimension gpop(natmax,nsmp)
        logical wrtdm,wrtmp,wrtgp

        wrtdm=.false. 
        wrtmp=.false.
        wrtgp=.false.

	write(iout,*) 'neltr=',neltr

c>>> k,l stand for AOs; i stands for MOs; n stands for time;

  	do k = 1 , ndim
     	   do l = 1 , ndim 
     	      dm(k,l,n)=0.d0
     	      do i = 1 , ndim
	         if(mod(neltr,2).eq.0) then
     	            if(i.le.neltr/2) 
     *	            dm(k,l,n)=dm(k,l,n)+2.0
     *                  *(zr(k,i)-(0.0,1.0)*zi(k,i))
     *                  *(zr(l,i)+(0.0,1.0)*zi(l,i))
                 else
                    if(i.le.int(neltr/2)) then
     	               dm(k,l,n)=dm(k,l,n)
     *                       +2.0*zr(k,i)*zr(l,i)
     *                       +2.0*zi(k,i)*zi(l,i)
                    endif
                    if(i.eq.int(neltr/2)+1) then
     	               dm(k,l,n)=dm(k,l,n)
     *                       +zr(k,i)*zr(l,i)
     *                       +zi(k,i)*zi(l,i)
                    endif
                 endif
     	      enddo
     	   enddo
     	enddo
	        	   
     	tot=0.d0
     	do k = 1 , ndim
     	   do l = 1 , ndim
     	      tot=tot+dm(k,l,n)*st(k,l,n)
     	   enddo
     	enddo
     	   
	if(wrtdm) then
	   write(iout,*)
     	   write(iout,*) ' density matrix '
           write(iout,*)
     	   do k = 1 , ndim	
     	      write(iout,'(500f8.4)') (dm(k,l,n),l=1,ndim)
     	   enddo
     	endif
     	write(iout,'(a)') ' total number of electrons '
     	write(iout,'(f20.10)') tot
     	   
     	if(wrtmp) then
     	   write(iout,*)
           write(iout,*) ' Population at molecular orbitals '
	   write(iout,*) 
	endif    	   
     	   
        do k = 1 , ndim
	   pop(k,n)=0.0
	   do l = 1 , ndim
	      pop(k,n)=pop(k,n)+dm(k,l,n)*st(k,l,n)
	   enddo
	   if(wrtmp) then
	   write(iout,'(i8,2a10,10f8.4)')lorb(k),indt(k,1),indt(k,2),
     *             pop(k,n)
           endif
        enddo
        
c>>> calculate gross population of atoms

	if(wrtgp) then
	   write(iout,*)
	   write(iout,*) ' Gross population at atoms '
	   write(iout,*)
	endif 

	tot=0
	do na = 1 , natom

	   gpop(na,n)=0.0
	      
	   do k = 1 , ndim
	      if(indt(k,1)(1:2).eq.ac(na)(1:2)
     *	      .and.lorb(k).eq.iatom(na)) then
     	         gpop(na,n)=gpop(na,n)+pop(k,n)
     	      endif
     	   enddo
     	   
	   tot=tot+gpop(na,n)     	   
	      
	   if(wrtgp) then
     	      write(iout,'(2x,a2,i3,f10.5)') ac(na),iatom(na),gpop(na,n)
     	   endif
	   
	enddo 
	
	write(iout,*) ' total number of electrons = ',tot

	return
	end
	
	subroutine diagden(natom,ndim,n,st,ht)
	
c********************************************************************
c diagnalize the density matrix 
c********************************************************************	
	
	implicit real*8 (a-h,o-z)
	parameter(nsmp=10,natmax=200,ndmax=500,neigd=500)
        logical invs,alle,eonly
	real*4 dm,pop
	real*4 st(ndmax,ndmax,nsmp),ht(ndmax,ndmax,nsmp)
	common /eigt/ dm(ndmax,ndmax,nsmp),pop(ndmax,nsmp)
        common/iounit/iin,iout,istart(40),iungf,iunhl            
        common/eigv/omcm(neigd),zr(neigd,neigd),zi(neigd,neigd),ne,nblw
        common/secl/s(neigd,neigd),h(neigd,neigd),nmat
        dimension sss(ndmax,ndmax)
        dimension temp(ndim,ndim),temp1(ndim,ndim)
        dimension ddd(ndim,ndim)
        dimension occ(ndim)
        
        invs=.false.
        alle=.false.
        eonly=.false.
        ellow=-1000.d0
        elup=  1000.d0

        do i = 1 , ndim
          do j = 1 , ndim
            s(i,j)=0.0
            if(i.eq.j) s(i,j)=1.0
            h(i,j)=0.0
          enddo
        enddo
        do i = 1 , ndim
          do j = 1, i
            h(j,i)=0.0
            h(i,j)=dm(i,j,n)
            ddd(i,j)=dm(i,j,n)
            ddd(j,i)=dm(j,i,n)
          enddo
        enddo
        nmat=ndim
        
        call diagm(ellow,elup,invs,alle,eonly)
        
        do k = 1 , ndim
           do l = 1 , ndim
              sss(k,l)=0.0
              do i = 1 , ndim
                 do j = 1 , ndim
                    sss(k,l)=sss(k,l)+zr(i,k)*st(i,j,n)*zr(j,l)
                 enddo
              enddo
           enddo
        enddo
        
        sum=0.0
        do i = 1 , ndim
           do j = 1 , ndim
              sum=sum+ddd(i,j)*st(i,j,n)
           enddo
        enddo
        write(iout,*) 
        write(iout,*) ' check normalization'
        write(iout,*) sum
        
        sum=0.0
        do i = 1 , ndim
           occ(i)=omcm(i)*sss(i,i)
           sum=sum+occ(i)
        enddo

	write(iout,*) ' eigenoccupancies '
        write(iout,'(10f6.2)') (occ(i),i=1,ndim)
        write(iout,*) ' total number of electrons'
        write(iout,'(f10.5)') sum
        write(iout,*)

	return
	end