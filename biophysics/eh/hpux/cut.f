	PROGRAM CUT
c********************************************************
c cut a pdb file
c********************************************************

	implicit real (a-h,o-z)
	character string1*100,string2*100,head*3
	parameter (natmax=500)
	dimension x(natmax),y(natmax),z(natmax)
	character*4 atom(natmax)
	character*3 resi(natmax)
	character*4 segm(natmax)
	dimension ires(natmax)

c>>> read original pdb file

	read(5,'(a100)')  string1
	read(5,'(a100)')  string2		
	
	do i = 1 , natmax
	
	   read(5,'(a3,9x,a4,1x,a3,i6,4x,3f8.3,18x,a4)')
     &          head,atom(i),resi(i),ires(i),x(i),y(i),z(i),segm(i)
     	   if(head.eq.'END') then 
     	      ntot=i-1
     	      go to 1001
     	   endif
	enddo

 1001   continue

	
c>>> create a smaller part

	write(6,'(a100)') string1
	write(6,'(a100)') string2
	isel=0
	do i = 1 , ntot
	   if(segm(i).eq.'PAR1') then
	      isel=isel+1
              write(6,'(a4,i7,1x,a4,1x,a3,i6,4x,3f8.3,18x,a4)')
     &           'ATOM',isel,atom(i),resi(i),ires(i),
     &           x(i),y(i),z(i),segm(i)
           endif
	enddo
	write(6,'(a3)')'END'

	end