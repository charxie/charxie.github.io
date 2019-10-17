	subroutine part2(tau,natom,ndim,nstp,tmd,st,ht,ft,st2,ht2)
	
c********************************************************************
c part 2: solve the time-dependent Schroedinger equation 
c         by using the interpolated Hamiltonian
c********************************************************************	
	
	implicit real*8 (a-h,o-z)
	parameter(nsmp=10,ndmax=500)
	logical restart,dotsa,rite
	real*4 tmd(nsmp),st(ndmax,ndmax,nsmp),ht(ndmax,ndmax,nsmp),
     *         ft(ndmax,ndmax,nsmp)	
	real*4 vtr(ndmax,ndmax,nsmp),vti(ndmax,ndmax,nsmp)
	real*4 cs(ndim,ndim,nstp,3),ch(ndim,ndim,nstp,3)
	real*4 cvr(ndim,ndim,nstp,3),cvi(ndim,ndim,nstp,3)
	real*4 sao,trt
	real*4 st2(2,2,nsmp),ht2(2,2,nsmp)
        common/iounit/iin,iout,istart(40),iungf,iunhl		
        common/option/restart,dotsa,nstep,unit
        common/nda/ndon,nacc,idon,iacc,nb,iold(ndmax),jold(ndmax)
        common/nonadi/sao(ndmax,ndmax,nsmp),trt(ndmax,ndmax,nsmp)
        dimension dij(nstp),sij(nstp),tij(nstp),hij(nstp)
	
	write(iout,*) 
	write(iout,*) '  Part 2 called '
	write(iout,*)	
	
	rite=.true.
	
	if(restart) then
	   if(.not.dotsa) then
	      open(11,file='popuns',status='unknown',access='append')
	      open(46,file='pallns',status='unknown',access='append')
	   endif
	   if(dotsa) then
	      open(33,file='tdat',status='unknown',access='append')
	      open(13,file='popu2s',status='unknown',access='append')
	   endif
	else
	   if(.not.dotsa) then
	      open(11,file='popuns',status='unknown',access='sequential')
	      open(46,file='pallns',status='unknown',access='sequential')
	   endif
	   if(dotsa) then
	      open(33,file='tdat',status='unknown',access='sequential')
	      open(13,file='popu2s',status='unknown',access='sequential')
	   endif
	endif
	open(44,file='ener')
	
	if(.not.dotsa) then

     	  call hint(ndim,nstp,tmd,st,ht,ft,vtr,vti,cs,ch,cvr,cvi)
	
  	  if(0.ne.0) then
	    do n = 1 , nstp
	      tij(n)=tmd(n)
	    enddo
	    dmax=0.0
	    do i = 1 , ndim
	      do j = 1 , ndim
	         do n = 1 , nstp
	            sij(n)=st(i,j,n)
	            hij(n)=ht(i,j,n)
	         enddo
	         call deriv1(tij,hij,dij,nstp)
	         do n = 1 , nstp
	            if(dabs(dij(n)).gt.dmax) then
	               dmax=abs(dij(n))
	               nmax=n
	               imax=i
	               jmax=j
	            endif
	         enddo
	      enddo
	    enddo
	    write(iout,*)imax,jmax,nmax,dmax
          endif
	
	  call dynamics(tau,ndim,nstp,tmd,st,ht,vtr,vti,cs,ch,cvr,cvi)

	endif
	
	if(dotsa) then
	  if(.not.restart) then
	    do n = 1 , nstp
	       write(33,'(i10,10f20.10)')
     *           n,ht2(1,2,n),ht2(1,1,n)-ht2(2,2,n),st2(1,2,n)	   
	    enddo
	  else
	    call getsteps(nstp1,dotsa)
	    do n = 3 , nstp
	      write(33,'(i10,10f20.10)')nstp1-nstp+n,
     *          ht2(1,2,n),ht2(1,1,n)-ht2(2,2,n),st2(1,2,n)	   
	    enddo
	  endif
	  call dyntsa(tau,nstp,tmd,st2,ht2)
	endif
	
	close(11)
	close(13)
	close(33)
	close(44)
	close(46)
	
	return
	end