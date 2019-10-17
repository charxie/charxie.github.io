c rebuild the Hamiltonian in the space spanned by the localized
c states of donor and acceptor, and eigenstates of the bridge,

    	subroutine bspace(n,h,s,hnew)
    	
    	implicit real*8 (a-h,o-z)
 	parameter(nmax=248,nm2=nmax-2)
	parameter(neigd=248)
        common /eigv/ omcm(neigd),zr(neigd,neigd),zi(neigd,neigd),
     :                ne,nblw
	dimension h(n,n),hnew(n,n),s(n,n)
	dimension hbr(n-2,n-2),gbr(n-2,n-2),sbr(n-2,n-2)
	dimension hdb(n-2),hba(n-2)
	dimension hdbnew(n-2),hbanew(n-2),hbrnew(n-2,n-2)
	dimension hmes(n-2,n-2),indx(n-2)
	logical restart
        common/ida/time0,edon,eacc,id,ia,i0,j0,restart,nstep
        common /io/ iin,iout
        
        if(id.gt.ia) then
           idon=ia
           iacc=id
        endif

c--> form the bridge Hamiltonian matrix

     	do i = 1 , n
     	   do j = 1 , n
     	      if(i.lt.idon.and.
     :	         j.lt.idon) then
                 hbr(i,j)=h(i,j)
              endif
     	      if(i.gt.idon.and.i.lt.iacc.and.
     :           j.lt.idon) then
                 hbr(i-1,j)=h(i,j)
              endif
     	      if(i.gt.iacc.and.
     :           j.lt.idon) then
                 hbr(i-2,j)=h(i,j)
              endif
     	      if(i.lt.idon.and.
     :           j.gt.idon.and.j.lt.iacc) then
                 hbr(i,j-1)=h(i,j)
              endif
     	      if(i.gt.idon.and.i.lt.iacc.and.
     :           j.gt.idon.and.j.lt.iacc) then
                 hbr(i-1,j-1)=h(i,j)
              endif
     	      if(i.gt.iacc.and.
     :           j.gt.idon.and.j.lt.iacc) then
                 hbr(i-2,j-1)=h(i,j)
              endif
     	      if(i.lt.idon.and.
     :           j.gt.iacc) then
                 hbr(i,j-2)=h(i,j)
              endif
     	      if(i.gt.idon.and.i.lt.iacc.and.
     :           j.gt.iacc) then
                 hbr(i-1,j-2)=h(i,j)
              endif
     	      if(i.gt.iacc.and.
     :           j.gt.iacc) then
                 hbr(i-2,j-2)=h(i,j)
              endif
     	   enddo
     	enddo

c--> build the donor-bridge/acceptor-bridge coupling arrays

        do i = 1 , n
           if(i.lt.idon) then
              hdb(i)=h(i,idon)
           endif
           if(i.gt.idon.and.i.lt.iacc) then
              hdb(i-1)=h(i,idon)
           endif
           if(i.gt.iacc) then
              hdb(i-2)=h(i,idon)
           endif
           if(i.lt.idon) then
              hba(i)=h(i,iacc)
           endif
           if(i.gt.idon.and.i.lt.iacc) then
              hba(i-1)=h(i,iacc)
           endif
           if(i.gt.iacc) then
              hba(i-2)=h(i,iacc)
           endif
        enddo        
c	write(iout,'(500f10.5)') (hdb(i),i=1,n-2)
c	write(iout,'(500f10.5)') (hba(i),i=1,n-2)     	

c--> build the bridge overlap matrix

	do i = 1 , n-2
	   do j = 1 , n-2
	      sbr(i,j)=0.0
	      if(i.eq.j) sbr(i,j)=1.0
	   enddo
	enddo

c        do i = 1 , n-2
c          write(iout,'(500f8.4)')(hbr(i,j),j=1,n-2)
c        enddo

	call eigen(-1000.0,1000.0,.true.,n-2,hbr,sbr)
	
	do i = 1 , n-2
	   hdbnew(i)=0.0
	   hbanew(i)=0.0
	   do k = 1 , n-2
	      hdbnew(i)=hdbnew(i)+hdb(k)*zr(k,i)
	      hbanew(i)=hbanew(i)+hba(k)*zr(k,i)
	   enddo
	enddo
c	write(iout,'(500f10.5)') (hdbnew(i),i=1,n-2)
c	write(iout,'(500f10.5)') (hbanew(i),i=1,n-2)     	
	
	do i = 1 , n-2
	   do j = 1 , n-2
	      hbrnew(i,j)=0.0
	   enddo
	   hbrnew(i,i)=omcm(i)
	enddo	
	
	hnew(1,1)=h(idon,idon)
	hnew(1,n)=h(idon,iacc)
	do j = 2 , n-1
	   hnew(1,j)=hdbnew(j-1)
	   hnew(j,1)=hnew(1,j)
	enddo
	hnew(n,1)=h(iacc,idon)
	hnew(n,n)=h(iacc,iacc)
	do j = 2 , n-1
	   hnew(n,j)=hbanew(j-1)
	   hnew(j,n)=hnew(n,j)
	enddo
	do i = 2 , n-1
	   do j = 2 , i
	      hnew(i,j)=hbrnew(i-1,j-1)
	      hnew(j,i)=hnew(i,j)
	   enddo
	enddo	
	
        return
        end