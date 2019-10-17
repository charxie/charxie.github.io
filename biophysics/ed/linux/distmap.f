c calculate the distance matrices

	subroutine distmap(ndim,protein)
	implicit real*4 (a-h,o-z)
	character dummy*1,head*3,protein*4
	parameter (natmax=300)
	character filename*8,bondname*2,anti*1
	character*4 symbol,resname
        common/iden/idres(natmax),symbol(natmax),resname(natmax),
     :  x(natmax),y(natmax),z(natmax)
        common/io/iin,iout
	common/orbital/nbd,npb,nlp
	dimension bondname(natmax),anti(natmax),label(natmax,2)
	
	do i = 1 , natmax
	   x(i)=0.0
	   y(i)=0.0
	   z(i)=0.0
	enddo

	filename(1:4)=protein
	filename(5:8)='.pdb'

	open(91,file=filename,status='old')
	read(91,'(a1)') dummy
	read(91,'(a1)') dummy
        do i = 1 , natmax
	   read(91,'(a3,9x,a4,1x,a4,i5,4x,3f8.3,22x)')
     &          head,symbol(i),resname(i),idres(i),x(i),y(i),z(i)
     	   if(head.eq.'END') then 
     	      ntot=i-1
     	      go to 1001
     	   endif
	enddo
 1001   continue
        if(head.ne.'END') stop 'too many atoms!'
	close(91)
	   
	open(92,file='tsign',status='old')
	
	nbd=0
	nlp=0
	npb=0
	do i = 1 , ndim
	   read(92,'(i10,2x,a2,a1)')idummy,bondname(i),anti(i)
	   if(anti(i).eq.' ') then
	      if(bondname(i).eq.'BD') nbd=nbd+1
	      if(bondname(i).eq.'LP') nlp=nlp+1
	      if(bondname(i).eq.'PB') npb=npb+1
	   endif
	enddo
	if(nbd*2+npb*2+nlp.ne.ndim.or.idummy.ne.ndim) 
     :     stop 'error in TSIGN'
	
	close(92)   
	
	return
	end