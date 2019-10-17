	subroutine outham(index,imat,jmat,tau,ndim,nstep,ncut,hami,hamf)
	implicit real*8 (a-h,o-z)
	real*4 hami(ndim,ndim,nstep)
	real*4 hamf(ndim,ndim,0:ncut)
	character*16 filename1,filename2

	twopi=8.d0*datan(1.d0)
	period=nstep*tau
	comega=0.5*twopi/period
	
	if(index.eq.0) filename1(1:14)='floquet/h/hami'
	if(index.eq.1) filename1(1:14)='floquet/h/hamj'
	filename2(1:14)='floquet/h/hamf'
	
	do i = imat , imat
	   if(i.lt.10) then
	      filename1(15:15)='0'
	      filename1(16:16)=char(48+i)
	      filename2(15:15)='0'
	      filename2(16:16)=char(48+i)
	   elseif(i.lt.100) then
	      nn1=mod(i,10)
	      nn2=int(i/10)
	      filename1(15:15)=char(48+nn2)
	      filename1(16:16)=char(48+nn1)
	      filename2(15:15)=char(48+nn2)
	      filename2(16:16)=char(48+nn1)
	   endif	
	   open(31,file=filename1)
	   do nt = 1 , nstep
	      write(31,'(500f10.5)')tau*nt/1.51653,
     :                              (hami(i,j,nt),j=jmat,jmat)
	   enddo 
	   close(31)
	   open(31,file=filename2)
	   do nt = 0 , ncut
	      write(31,'(500f10.5)')comega*nt*1.51653,
     :                              (hamf(i,j,nt),j=jmat,jmat)
	   enddo 
	   close(31)
	enddo
	
	return
	end
   	
c print out standard output for the static case
c   	
   	subroutine output(n,hami,ovlp,tau,protein,neign)
	implicit real*8 (a-h,o-z)
	parameter(neigd=248)
	dimension hami(n,n),ovlp(n,n)
	logical restart,lcao,hole
	character protein*4
        common /eigv/ omcm(neigd),zr(neigd,neigd),zi(neigd,neigd),
     :                ne,nblw
        common /ida/ idon,iacc,i0,j0,restart,time0,nstep,edon,eacc
        common /io/ iin,iout
	common/orbital/nbd,npb,nlp
        
        lcao=.true.
        hole=.false.
        
        call distmap(n,protein)
        if(mod(n-nlp,2).ne.0) stop 'error with lonepairs!'
        if(lcao) then
           neign=(n+nlp)/2
        else
           if(hole) then
              neign=n-1
           else
              neign=1
           endif
        endif
        
        write(iout,*)'number of lonepairs=',nlp
        write(iout,*)'indices of the two eigenstates=',neign-1,neign
	
	write(iout,*) ' hamiltonian '
c	write(iout,'(4x,500i8)') (i,i=1,n)
c	do i = 1 , n
c	   write(iout,'(i8,500f8.4)') i,(hami(i,j),j=1,n)
c	enddo
	write(iout,*)
	write(iout,*) ' overlap matrix '
c	do i = 1 , n
c	   write(iout,'(500f8.4)') (ovlp(i,j),j=1,n)
c	enddo
	write(iout,*)
	write(iout,*) ' eigenenergies '
	write(iout,'(10f10.4)') (omcm(i),i=1,n)
	write(iout,*)
	
	split=omcm(neign)-omcm(neign-1)
	
	write(iout,*) 
	write(iout,*) ' energy split (eV) '
	write(iout,'(f30.20)') split
	write(iout,*)
	
	planck=8.0*atan(1.0)
	span=planck/split

c... time unit in hbar/eV=1.055/1.6*femtosecond
	
	write(iout,*)
	write(iout,*) ' Period (Unit = 0.6594 fs) '
	write(iout,'(f20.1)') span
	write(iout,*)

	write(iout,*)
	write(iout,*) ' Estimated number of steps needed '
	write(iout,'(f20.1)') span/tau
	write(iout,*)
	
	do i = neign-1 , neign
	   write(iout,*) ' eigenvector(',omcm(i),')'
	   write(iout,'(10i8)') (j,j=1,n)
	   write(iout,'(10f8.4)') (zr(j,i),j=1,n)
c	   write(iout,'(500f8.4)') (zi(j,i),j=1,n)
c	   write(iout,1010) i,idon,zr(idon,i),iacc,zr(iacc,i)
	enddo
 1010   format(i8,8x,i8,'(',f8.4,')',8x,i8,'(',f8.4,')')
	
	open(57,file='vect')
	do i = 1 , n
	   write(57,'(10f10.5)') i,zr(i,neign-1),zr(i,neign)
	enddo
	close(57)
	
	return
        end