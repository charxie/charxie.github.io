	parameter(nres=15)
	dimension idres(nres),rmsmain(nres),rmsside(nres)
	dimension nmain(nres),nside(nres)
	character*1 dummy
	
	open(20,file='index.num')
	do i = 1 , nres
	  read(20,*) idres(i),nmain(i),nside(i)
	enddo
	
	open(10,file='rms.com')
	
	read(10,*) dummy
	read(10,*) dummy
	read(10,*) dummy
	read(10,*) dummy
	do i = 1 , nres
	  read(10,'(i8,2f15.5)') idres(i),rmsmain(i),rmsside(i)
	enddo
	
	write(6,*) 'crystal'

	summain=0.0
	sumside=0.0
	nummain=0
	numside=0
	do k = 1 , 8
	  summain=summain+rmsmain(k)**2*nmain(k)
	  nummain=nummain+nmain(k)
	  sumside=sumside+rmsside(k)**2*nside(k)
	  numside=numside+nside(k)
	enddo
	summain=sqrt(summain/nummain)
	sumside=sqrt(sumside/numside)
	write(6,*) ' Loop 281-288 ', summain, sumside
	
	summain=0.0
	sumside=0.0
	nummain=0
	numside=0
	do k = 9 , 15
	  summain=summain+rmsmain(k)**2*nmain(k)
	  nummain=nummain+nmain(k)
	  sumside=sumside+rmsside(k)**2*nside(k)
	  numside=numside+nside(k)
	enddo
	summain=sqrt(summain/nummain)
	sumside=sqrt(sumside/numside)
	write(6,*) ' Loop 378-384 ', summain, sumside

	read(10,*) dummy
	do i = 1 , nres
	  read(10,'(i8,2f15.5)') idres(i),rmsmain(i),rmsside(i)
	enddo
		
	write(6,*) 'endpoint'

	summain=0.0
	sumside=0.0
	nummain=0
	numside=0
	do k = 1 , 8
	  summain=summain+rmsmain(k)**2*nmain(k)
	  nummain=nummain+nmain(k)
	  sumside=sumside+rmsside(k)**2*nside(k)
	  numside=numside+nside(k)
	enddo
	summain=sqrt(summain/nummain)
	sumside=sqrt(sumside/numside)
	write(6,*) ' Loop 281-288 ', summain, sumside
	
	summain=0.0
	sumside=0.0
	nummain=0
	numside=0
	do k = 9 , 15
	  summain=summain+rmsmain(k)**2*nmain(k)
	  nummain=nummain+nmain(k)
	  sumside=sumside+rmsside(k)**2*nside(k)
	  numside=numside+nside(k)
	enddo
	summain=sqrt(summain/nummain)
	sumside=sqrt(sumside/numside)
	write(6,*) ' Loop 378-384 ', summain, sumside

	read(10,*) dummy
	do i = 1 , nres
	  read(10,'(i8,2f15.5)') idres(i),rmsmain(i),rmsside(i)
	enddo
		
	write(6,*) 'average'

	summain=0.0
	sumside=0.0
	nummain=0
	numside=0
	do k = 1 , 8
	  summain=summain+rmsmain(k)**2*nmain(k)
	  nummain=nummain+nmain(k)
	  sumside=sumside+rmsside(k)**2*nside(k)
	  numside=numside+nside(k)
	enddo
	summain=sqrt(summain/nummain)
	sumside=sqrt(sumside/numside)
	write(6,*) ' Loop 281-288 ', summain, sumside
	
	summain=0.0
	sumside=0.0
	nummain=0
	numside=0
	do k = 9 , 15
	  summain=summain+rmsmain(k)**2*nmain(k)
	  nummain=nummain+nmain(k)
	  sumside=sumside+rmsside(k)**2*nside(k)
	  numside=numside+nside(k)
	enddo
	summain=sqrt(summain/nummain)
	sumside=sqrt(sumside/numside)
	write(6,*) ' Loop 378-384 ', summain, sumside

	
	end