	subroutine readh(tau,natom,ndim,nstp,tmd,st,ht,ft)

c*******************************************************************
c build the time series of the Hamiltonian matrices
c*******************************************************************

	implicit real*8 (a-h,o-z)
	parameter(nsmp=10,natmax=200,ndmax=500)
	logical sp,pd,restart,dotsa,rite
	integer sorb
        character ac*4,symbol*2,segname*4,resname*3
	real*4 tmd(nsmp),st(ndmax,ndmax,nsmp),ht(ndmax,ndmax,nsmp),
     *         ft(ndmax,ndmax,nsmp)
	real*4 rx,ry,rz
        common/ato/ac(natmax),symbol(40)                                          
        common/sij/s(ndmax,ndmax),h(ndmax,ndmax),c(ndmax),
     *  maxs(natmax),maxp(natmax),maxd(natmax),sorb(2*natmax),
     *  sp(natmax),pd(natmax)
        common/option/restart,dotsa,nstep,unit
        common/iounit/iin,iout,istart(40),iungf,iunhl            	
        common/nda/ndon,nacc,idon,iacc,nb,iold(ndmax),jold(ndmax)
        common/atom2/exps(40),exps2(40),expp(40),expp2(40),expd(40),              
     *  expd2(40),expf(40),expf2(40),cs1(40),cs2(40),cp1(40),cp2(40),             
     *  cd1(40),cd2(40),cf1(40),cf2(40),couls(40),coulp(40),could(40),            
     *  coulf(40),x(natmax),y(natmax),z(natmax),ires(natmax),
     *  resname(natmax),segname(natmax)
        common/traj/rx(natmax,nsmp),ry(natmax,nsmp),rz(natmax,nsmp)
	character filename*15,protein*4
	
	twopi=8.d0*datan(1.d0)
	protein='azur'
	rite=.false.

	if(.not.restart) then
	   nstp0=1
	else
	   call getsteps(nstp0,dotsa)
	   nstp0=nstp0-1
	endif
	if(.not.restart) then
	   period=unit*(nstp-1)
	else
	   period=unit*(nstp-2)
	endif
	tau=period/nstep
	comega=0.5*twopi/period
	deltat=unit
	deltao=comega
	nstp1=nstp0+nstp-1
        do i = nstp0 , nstp1
           tmd(i-nstp0+1)=unit*(i-1)
        enddo

	filename(1:5)='strk/'
	filename(6:9)=protein
	filename(10:10)='.'
	
	do n = nstp0 , nstp1
	   neff=n-nstp0+1
	   write(iout,*) 'reading ',n,' -th pdb file ... '
	   if(n.lt.10) then
	      filename(11:14)='0000'
	      filename(15:15)=char(48+n)
	   elseif(n.lt.100) then
	      nn1=mod(n,10)
	      nn2=int(n/10)
	      filename(11:13)='000'
	      filename(14:14)=char(48+nn2)
	      filename(15:15)=char(48+nn1)
	   elseif(n.lt.1000) then
	      nn1=int(n/100)
	      nnx=mod(n,100)
	      nn2=int(nnx/10)
	      nn3=mod(nnx,10)
	      filename(11:12)='00'
	      filename(13:13)=char(48+nn1)
	      filename(14:14)=char(48+nn2)
	      filename(15:15)=char(48+nn3)
	   elseif(n.lt.10000) then
	      nn1=int(n/1000)
	      nnx=mod(n,1000)
	      nn2=int(nnx/100)
	      nny=mod(nnx,100)
	      nn3=int(nny/10)
	      nn4=mod(nny,10)
	      filename(11:11)='0'
	      filename(12:12)=char(48+nn1)
	      filename(13:13)=char(48+nn2)
	      filename(14:14)=char(48+nn3)
	      filename(15:15)=char(48+nn4)
	   elseif(n.lt.100000) then
	      nn1=int(n/10000)
	      nnx=mod(n,10000)
	      nn2=int(nnx/1000)
	      nny=mod(nnx,1000)
	      nn3=int(nny/100)
	      nnz=mod(nny,100)
	      nn4=int(nnz/10)
	      nn5=mod(nnz,10)
	      filename(11:11)=char(48+nn1)
	      filename(12:12)=char(48+nn2)
	      filename(13:13)=char(48+nn3)
	      filename(14:14)=char(48+nn4)
	      filename(15:15)=char(48+nn5)
	   else
	      STOP ' too many time steps ! '
	   endif
	   open(iin,file=filename,status='old')
	   call getsel
	   call getpos(natom,nstp)
	   do iat = 1 , natom
	      rx(iat,neff)=x(iat)
	      ry(iat,neff)=y(iat)
	      rz(iat,neff)=z(iat)
	   enddo
	   call getpar(natom,ndim,nstp)
	   natom2=natom+natom
	   call smat(natom2,natom,ndim)
	   call hmat(natom,ndim)
	   do i = 1 , ndim
	      do j = 1 , ndim
	         st(i,j,neff)=s(i,j)
	         ht(i,j,neff)=h(i,j)
	      enddo
	   enddo
	   close(iin)	   
	enddo
	call storesteps(nstp1,dotsa)
	
	if(rite) then
	   do n = 1 , nstp
	      write(iout,*) n
	      write(iout,*)
	      do i = 1 , ndim
	         write(iout,'(500f8.4)')(st(i,j,n),j=1,ndim)
	      enddo
	   enddo
	   stop
	endif

	return
	end