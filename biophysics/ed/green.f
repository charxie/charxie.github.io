c calculate the static Green's function

    	subroutine caltda(n,h,e,tda,hdd,haa)
    	
    	implicit real*8 (a-h,o-z)
	parameter(neigd=248)
        common /eigv/ omcm(neigd),zr(neigd,neigd),zi(neigd,neigd),
     :                ne,nblw
        logical restart
	dimension h(n,n)
	dimension hbr(n-2,n-2),gbr(n-2,n-2),sbr(n-2,n-2)
	dimension hdb(n-2),hba(n-2)
	dimension hmes(n-2,n-2),indx(n-2)
        common /ida/ id,ia,i0,j0,restart,time0,nstep,edon,eacc
        common /io/ iin,iout
        
        if(id.gt.ia) then
           idon=ia
           iacc=id
        else
           idon=id
           iacc=ia
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

c--> form the donor-bridge/acceptor-bridge coupling arrays

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

	do i = 1 , n-2
	   do j = 1 , n-2
	      sbr(i,j)=0.0
	      if(i.eq.j) sbr(i,j)=1.0
	   enddo
	enddo

	do i = 1 , n-2
	   do j = 1 , n-2
	      hmes(i,j)=e*sbr(i,j)-hbr(i,j)
	      gbr(i,j)=0.d0
	   enddo
	   gbr(i,i)=1.0
	enddo
	
	call ludcmp(hmes,n-2,n-2,indx,d)

	do j = 1 , n-2
	   call lubksb(hmes,n-2,n-2,indx,gbr(1,j))
	enddo

	tda=h(idon,iacc)
	hdd=h(idon,idon)
	haa=h(iacc,iacc)

	do i = 1 , n-2
	   do j = 1 , n-2
	      tda=tda+hdb(i)*gbr(i,j)*hba(j)
	      hdd=hdd+hdb(i)*gbr(i,j)*hdb(j)
	      haa=haa+hba(i)*gbr(i,j)*hba(j)
	   enddo
	enddo
     
     	write(iout,'(a5,1x,f10.5,1x,a4,1x,f20.10)')' TDA(',e,') = ',tda
	write(iout,*)
	write(iout,'(a)') ' Two-State Approximation'
	write(iout,'(2f20.10)') h(id,id),h(id,ia)
	write(iout,'(2f20.10)') h(ia,id),h(ia,ia)
	write(iout,'(2f20.10)') hdd,tda
	write(iout,'(2f20.10)') tda,haa

        return
        end        


c Reduce the Hamiltonian matrix to effective N-states Hamiltonian(N=3)

	subroutine hred(n,h,e,heff)
	
	implicit real*8 (a-h,o-z)
	dimension heff(3,3)
	dimension h(n,n)
        common /io/ iin,iout
	
	k2=int(n/2)

	do i = 1 , 3
	   do j = 1 , 3
	      call reduce(n,h,e,k2,i,j,heff(i,j))
	   enddo
	enddo
	
	write(iout,*)
	do i = 1 , 3
	   write(iout,'(100f10.5)') (heff(i,j),j=1,3)
	enddo

	return
	end

    	subroutine reduce(n,h,e,k2,ipp,jpp,hef12)
    	
    	implicit real*8 (a-h,o-z)
 	parameter(nmax=8,nm2=nmax-2,nm3=nmax-3)
	dimension h(n,n),hbr(nm2,nm2)
	dimension hqq(nm3,nm3),gqq(nm3,nm3),sqq(nm3,nm3)
	dimension hmes(nm3,nm3),indx(nm3),hass(nm3,nm3)
	dimension hpq(3,nm3),hqp(nm3,3)
	dimension hbrpq(3,nm2),hbrqp(nm2,3)
	dimension hbrpp(2,2),hpp(3,3)
        common /io/ iin,iout
	
	do i = 1 , n-3
	   do j = 1 , n-3
	      sqq(i,j)=0.0
	      if(i.eq.j) sqq(i,j)=1.0
	   enddo
	enddo
	
     	do i = 1 , n-2
     	   do j = 1 , n-2
              hbr(i,j)=h(i+1,j+1)
     	   enddo
     	enddo
     	
     	hbrpp(1,1)=h(1,1)
     	hbrpp(1,2)=h(1,n)
     	hbrpp(2,1)=h(n,1)
     	hbrpp(2,2)=h(n,n)
     	
     	do i = 1 , n-2
     	   hbrpq(1,i)=h(1 ,i+1)
     	   hbrpq(2,i)=h(k2,i+1)
     	   hbrpq(3,i)=h(n ,i+1)
     	   hbrqp(i,1)=h(i+1,1 )
     	   hbrqp(i,2)=h(i+1,k2)
     	   hbrqp(i,3)=h(i+1,n )
     	enddo
     	
     	write(iout,*) (hbrpq(1,i),i=1,n-2)
     	write(iout,*) (hbrpq(2,i),i=1,n-2)
	write(iout,*) (hbrpq(3,i),i=1,n-2)
     	write(iout,*) (hbrqp(i,1),i=1,n-2)
     	write(iout,*) (hbrqp(i,2),i=1,n-2)
     	write(iout,*) (hbrqp(i,3),i=1,n-2)
     	
     	kb2=k2-1
     	do i = 1 , n-2
     	   do j = 1 , n-2
     	      if(i.lt.kb2.and.j.lt.kb2) then
                 hqq(i,j)=hbr(i,j)
              endif
     	      if(i.gt.kb2.and.j.lt.kb2) then
                 hqq(i-1,j)=hbr(i,j)
              endif
     	      if(i.lt.kb2.and.j.gt.kb2) then
                 hqq(i,j-1)=hbr(i,j)
              endif
     	      if(i.gt.kb2.and.j.gt.kb2) then
                 hqq(i-1,j-1)=hbr(i,j)
              endif
     	   enddo
     	enddo
      	
	do i = 1 , n-2
	   if(i.lt.kb2) then
	      hpq(1,i)=hbrpq(1,i)
	      hpq(2,i)=hbrpq(2,i)
	      hpq(3,i)=hbrpq(3,i)
	      hqp(i,1)=hbrqp(i,1)
	      hqp(i,2)=hbrqp(i,2)
	      hqp(i,3)=hbrqp(i,3)
	   endif
	   if(i.gt.kb2) then
	      hpq(1,i-1)=hbrpq(1,i)
	      hpq(2,i-1)=hbrpq(2,i)
	      hpq(3,i-1)=hbrpq(3,i)
	      hqp(i-1,1)=hbrqp(i,1)
	      hqp(i-1,2)=hbrqp(i,2)
	      hqp(i-1,3)=hbrqp(i,3)
	   endif
	enddo
	
	hpp(1,1)=h(1 ,1 )
	hpp(1,2)=h(1 ,k2)
	hpp(1,3)=h(1 ,n )
	hpp(2,1)=h(k2,1 )
	hpp(2,2)=h(k2,k2)
	hpp(2,3)=h(k2,n )
	hpp(3,1)=h(n ,1 )
	hpp(3,2)=h(n ,k2)
	hpp(3,3)=h(n ,n )
		
     	write(iout,*) 
	do i = 1 , n
	   write(iout,'(100f10.5)') (h(i,j),j=1,n)
	enddo     	
	write(iout,*)
	do i = 1 , n-3
	   write(iout,'(100f10.5)') (hqq(i,j),j=1,n-3)
	enddo     	
     	write(iout,*) (hpq(1,i),i=1,n-3)
     	write(iout,*) (hpq(2,i),i=1,n-3)
	write(iout,*) (hpq(3,i),i=1,n-3)
     	write(iout,*) (hqp(i,1),i=1,n-3)
     	write(iout,*) (hqp(i,2),i=1,n-3)
     	write(iout,*) (hqp(i,3),i=1,n-3)
     	
	do i = 1 , n-3
	   do j = 1 , n-3
	      hmes(i,j)=e*sqq(i,j)-hqq(i,j)
	      gqq(i,j)=0.d0
	      hass(i,j)=0.d0
	   enddo
	   gqq(i,i)=1.0
	   hass(i,i)=1.0
	enddo
	
	call ludcmp(hmes,n-3,n-3,indx,d)

	do j = 1 , n-3
	   call lubksb(hmes,n-3,n-3,indx,gqq(1,j))
	enddo
	
     	write(iout,*)
     	write(iout,*) ' Green function (', e, ')'
     	write(iout,*)
    	do k = 1 , n-3
     	   write(iout,'(100f10.5)') (gqq(k,l),l=1,n-3)
     	enddo 

c	call ludcmp(gqq,n-3,n-3,indx,d)
c	do j = 1 , n-3
c	   call lubksb(gqq,n-3,n-3,indx,hass(1,j))
c	enddo
c       write(iout,*)
c     	write(iout,*) ' Check matrix inversion '
c    	write(iout,*)
c    	do k = 1 , n-3
c    	   write(iout,'(100f10.5)') (hass(k,l),l=1,n-3)
c     	enddo 
     	

	hef12=hpp(ipp,jpp)
	do i = 1 , n-3
	   do j = 1 , n-3
	      hef12=hef12+hpq(ipp,i)*gqq(i,j)*hqp(j,jpp)
	   enddo
	enddo
	
        return
        end