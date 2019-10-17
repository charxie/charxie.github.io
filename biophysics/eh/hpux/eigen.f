
      subroutine eigen(n,homo,lumo)

c *************************************************************************
c Subroutine to get eigenenergies and eigenvectors
c omcm   ...... eigenenergies
c zr,zi  ...... eigenvectors
c *************************************************************************

      implicit real*8 (a-h,o-z)
      parameter (neigd=500,ndmax=500)
      logical invs,alle,eonly
      real*8 homo,lumo
      common /eigv/ omcm(neigd),zr(neigd,neigd),zi(neigd,neigd),ne,nblw
      common /secl/ s(neigd,neigd),h(neigd,neigd),nmat
      common/inst/ovlp(ndmax,ndmax),hami(ndmax,ndmax),dens(ndmax,ndmax)	      
      common/optv/rho,fermi,wstart,wend,
     *lb,le,lbr,ler,iovlb,iovle,lcharg,neltr,
     *distl,idistn,ieb,iee,isb,ise,ibb,ibe,ipope,jpope,irest        

c-->  initialization

      if(n.gt.neigd) stop ' Dimension exceeds limit in eigen.f ! '

      do i = 1 , n
         do j = 1 , n
            s(i,j)=0.0
            if(i.eq.j) s(i,j)=1.0
            h(i,j)=0.0
         enddo
      enddo

      do i = 1 , n
         do j = 1, i
            h(j,i)=0.0
            h(i,j)=hami(i,j)
            s(i,j)=ovlp(i,j)
         enddo
      enddo
      
      nmat=n
      
c--> solve secular equation

      invs=.false.
      alle=.false.
      eonly=.false.
      ellow=-1000.d0
      elup=  1000.d0

      call diagm(ellow,elup,invs,alle,eonly)

      if(mod(neltr,2).eq.0) then
         homo=omcm(neltr/2)
         lumo=omcm(neltr/2+1)
      else
         homo=omcm(int(neltr/2)+1)
         lumo=omcm(int(neltr/2)+2)
      endif	

      return
      end
      
      
	subroutine bspace(ndim,nstp)

c***********************************************************************
c diagonalize the bridge space
c***********************************************************************	

	implicit real*8 (a-h,o-z)
	parameter(natmax=200,ndmax=500,neigd=500)
	character ac*4,symbol*2,indt*10
	logical write,invs,alle,eonly,restart,dotsa	
        common/ato/ac(natmax),symbol(40)
        common/aoout/indt(1000,2),iatom(natmax),lorb(1000)
        common /inst/ ovlp(ndmax,ndmax),hami(ndmax,ndmax),
     *                dens(ndmax,ndmax)
        common /subh/ ovb(ndmax,ndmax),hmb(ndmax,ndmax)
        common/option/restart,dotsa,nstep,unit
        common/iounit/iin,iout,istart(40),iungf,iunhl     
        common/nda/ndon,nacc,idon,iacc,nb,iold(ndmax),jold(ndmax)
        common/eigv/omcm(neigd),zr(neigd,neigd),zi(neigd,neigd),ne,nblw
        common/secl/s(neigd,neigd),h(neigd,neigd),nmat
        common/optv/rho,fermi,wstart,wend,
     *  lb,le,lbr,ler,iovlb,iovle,lcharg,neltr,
     *  distl,idistn,ieb,iee,isb,ise,ibb,ibe,ipope,jpope,irest        
	
        invs=.false.
        alle=.false.
        eonly=.false.
        ellow=-1000.d0
        elup=  1000.d0
	write=.false.
	
	if(nb.ne.ndim-2) stop 'nb != ndim-2'
	
c---> define the bridge subhamiltonian

        write(iout,*) idon,iacc
		   
	if(idon.gt.iacc) then
	   id=iacc
	   ia=idon
	else
	   id=idon
	   ia=iacc
	endif
	do i = 1 , nb
	   do j = 1 , nb
	      ovb(i,j)=ovlp(iold(i),jold(j))
	      hmb(i,j)=hami(iold(i),jold(j))
	   enddo
	enddo
	
        do i = 1 , nb
          do j = 1 , nb
            s(i,j)=0.0
            if(i.eq.j) s(i,j)=1.0
            h(i,j)=0.0
          enddo
        enddo
        do i = 1 , nb
          do j = 1, i
            h(j,i)=0.0
            h(i,j)=hmb(i,j)
            s(i,j)=ovb(i,j)
          enddo
        enddo
c        write(iout,*)
c        do i = 1 , nb
c           write(iout,'(500f8.3)')(hmb(i,j),j=1,nb)
c        enddo
c        write(iout,*)
        
        nmat=nb
        
        call diagm(ellow,elup,invs,alle,eonly)
        
	write(iout,*) 'bridge eigenstates'
        write(iout,'(i10,f10.5)') (i,omcm(i),i=1,nb)
        nelbr=neltr-4
        if(mod(nelbr,2).eq.0) then
           brhomo=omcm(nelbr/2)
           brlumo=omcm(nelbr/2+1)
        else
           brhomo=omcm(int(nelbr/2)+1)
           brlumo=omcm(int(nelbr/2)+2)
        endif
        if(.not.restart) then
           write(24,'(i10,2f10.5)') nstp,brhomo,brlumo
        else
	   call getsteps(nstp1,dotsa)
           if(nstp.ge.3) then
              write(24,'(i10,2f10.5)') nstp1-10+nstp,brhomo,brlumo
           endif
        endif
        
        return
        end
      