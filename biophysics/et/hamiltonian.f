	subroutine hamiltonian0(ndim,hami)

	implicit real*8 (a-h,o-z)
	parameter (nmax=50)
	common/ranum/rann(1000)
	common/decay/CD,CA,CB(nmax),ED,EA,EM,EB,VD1,VNA,VBB,VBM
	dimension hami(ndim,ndim)
	
	CD=2.0
	CA=2.0
	do i = 1 , ndim
	   CB(i)=2
	enddo
	imid=int((ndim+1)/2)

	do i=1,ndim
	   if(i.eq.1) then
  	      hami(i,i)=ED
  	   elseif(i.eq.ndim) then
  	      hami(i,i)=EA
 	   elseif(ndim.gt.3.and.i.eq.imid) then
  	      hami(i,i)=EM
  	   else
  	      hami(i,i)=EB
  	   endif
	enddo
	
	do i = 1 , ndim
	   do j = i+1 , ndim
	      if(abs(j-i).le.1) then
	         if(i.eq.1) then
                   hami(i,j)=(VD1)*exp(-CD*real(j-i-1))
                 elseif(i.eq.imid.or.i.eq.imid-1) then
                   hami(i,j)=(VBM)*exp(-CB(i)*real(j-i-1))
          	 elseif(j.eq.ndim) then
	           hami(i,j)=(VNA)*exp(-CA*real(j-i-1))
	         else
	           hami(i,j)=(VBB)*exp(-CB(i)*real(j-i-1))
	         endif
	      else
	         hami(i,j)=0.d0
	      endif
	      hami(j,i)=hami(i,j)
	   enddo
	enddo
	

	return
	end



c<<<construct time-dependent Hamiltonian>>>

	subroutine hamiltonian(tau,istep,ndim,hami)

	implicit real*8 (a-h,o-z)
	parameter(nmax=50)
	dimension hami(ndim,ndim)
	common/ranum/rann(1000)
	common/decay/CD,CA,CB(nmax),ED,EA,EM,EB,VD1,VNA,VBB,VBM
	common/para/perio(nmax,nmax),omega(nmax,nmax),
     :              ampli(nmax,nmax),delta(nmax,nmax)
	
	twopi=8.d0*datan(1.d0)
		
	do i=1,ndim
	   if(i.eq.1) then
  	      hami(i,i)=ED+ampli(i,i)
     :  	       *dsin(istep*omega(i,i)*tau+delta(i,i))
  	   elseif(i.eq.ndim) then
  	      hami(i,i)=EA+ampli(i,i)
     :     	       *dsin(istep*omega(i,i)*tau+delta(i,i))
  	   else
  	      hami(i,i)=EB+rann(i)*1.0
     :          +ampli(i,i)*dsin(istep*omega(i,i)*tau+delta(i,i))
  	   endif
	enddo
	
	do i = 1 , ndim
	   do j = i+1 , ndim
	      if(abs(j-i).le.1) then
	         if(i.eq.1) then
        	    hami(i,j)=(VD1+ampli(i,j)*dsin(istep*omega(i,j)
     :              *tau+delta(i,j)))*exp(-CD*real(j-i-1))
          	 elseif(j.eq.ndim) then
	            hami(i,j)=(VNA+ampli(i,j)*dsin(istep*omega(i,j)
     :              *tau+delta(i,j)))*exp(-CA*real(j-i-1))
	         else
	            hami(i,j)=(VBB+rann(300+i)
     :              +ampli(i,j)*dsin(istep*omega(i,j)
     :              *tau+delta(i,j)))*exp(-CB(i)*real(j-i-1))
	         endif
	      else
	         hami(i,j)=0.d0
	      endif
	      hami(j,i)=hami(i,j)
	   enddo
	enddo

	return
	end