   	program ET
   	
c********************************************************************
c Solving the time-dependent Schroedinger equation for electron transfer
c********************************************************************   	
   	
	implicit real*8 (a-h,o-z)
	parameter (ndim=248,neigd=248,nstp=5,mstp=300)
	real*4 tmd(nstp),h(ndim,ndim,nstp),cof(ndim,ndim,nstp,3)
	real*4 teg(mstp)
	real*4 ft(ndim,ndim,mstp),eg(ndim,mstp),ev(ndim,ndim,mstp)
  	real*4 cofft(ndim,ndim,mstp,3)
	real*4 h2s(2,2,nstp)
	real*8 heig(2,2),voso(2,2),sovo(2,2)
	real*8 tem1(2,2),tem2(2,2)
	real*8 hami(ndim,ndim),ovlp(ndim,ndim)
	real*8 hnew(ndim,ndim)
	complex*8 psi(ndim)
	logical restart,static,quick,dotsa,brgap,both,hole,new2s
	logical alle, eigenspace, remft, turnoff
	common/eigv/omcm(neigd),zr(neigd,neigd),zi(neigd,neigd),ne,nblw
        common/ida/time0,edon,eacc,idon,iacc,i0,j0,restart,nstep
        common /io/ iin,iout
	common/orbital/nbd,npb,nlp
        character yesorno*1,protein*4

	data nstep/10000/
	data nele/2/
	data idon/10/
	data iacc/54/
	data i0/1/
	data j0/1/
	data etun/-9.0/
	data edon/-10.0/
	data eacc/-10.0/
	data iin/5/
	data iout/6/

	protein='azur'
	twopi=8.d0*datan(1.d0)
	restart=.false.
	static=.false.
	quick=.false.
	dotsa=.false.
	brgap=.false.
	new2s=.false.
	alle=.true.
	eigenspace=.true.
	remft=.false.
	turnoff=.false.
        
        both=.true.
        hole=.false.
        
c...initialize arrays (g77 requires explicit initialization)

        do n = 1 , nstp
           tmd(n)=0.0
           do j = 1 , ndim
              do i = 1 , ndim
                 h(i,j,n)=0.0
                 cof(i,j,n,1)=0.0
                 cof(i,j,n,2)=0.0
                 cof(i,j,n,3)=0.0
              enddo
           enddo
           h2s(1,1,n)=0.0
           h2s(1,2,n)=0.0
           h2s(2,1,n)=0.0
           h2s(2,2,n)=0.0
        enddo
        do j = 1 , ndim
           do i = 1 , ndim
              hami(i,j)=0.d0
              ovlp(i,j)=0.d0
              hnew(i,j)=0.d0
           enddo
        enddo
        do m = 1 , mstp
           teg(m)=0.0
           do j = 1 , ndim
              do i = 1 , ndim
                 ft(i,j,m)=0.0
                 ev(i,j,m)=0.0
              enddo
              eg(j,m)=0.0
           enddo
        enddo
        do j = 1 , ndim
           do i = 1 , ndim
              zr(i,j)=0.d0
              zi(i,j)=0.d0
           enddo
           omcm(j)=0.d0
        enddo
        
        do j = 1 , 2
           do i = 1 , 2
              heig(i,j)=0.d0
              voso(i,j)=0.d0
              sovo(i,j)=0.d0
              tem1(i,j)=0.d0
              tem2(i,j)=0.d0
           enddo
        enddo
        
        do i = 1 , ndim
           psi(i)=(0.0,0.0)
        enddo
        
c...read parameters	
	open(73,file='edea')
	read(73,*) idon,iacc
	read(73,*) etun
	read(73,*) edon,eacc
	read(73,*) eta
	close(73)
		
885	write(iout,*) ' restart ? (y/n) '
	read(iin,'(a)') yesorno
	if(yesorno.eq.'y') then
	   restart=.true.
	elseif(yesorno.eq.'n') then
	   restart=.false.
	else
	   write(iout,*) ' answer [y]es or [n]o '
	   goto 885
	endif
	
	if(restart) then
	   if(.not.dotsa.and..not.new2s) then
	      open(11,file='popu',status='unknown',access='append')
	      open(21,file='esit',status='unknown',access='append')
	      open(45,file='bgap',status='unknown',access='append')
	      open(46,file='pall',status='unknown',access='append')
	   endif
	   if(dotsa.or.new2s) then
	      open(11,file='popu2s',status='unknown',access='append')
	      open(12,file='esit2s',status='unknown',access='append')
	      open(13,file='tdat'  ,status='unknown',access='append')
	      open(14,file='hdha'  ,status='unknown',access='append')
	      open(15,file='homo'  ,status='unknown',access='append')
	   endif
	else
	   if(.not.dotsa.and..not.new2s) then
	      open(11,file='popu',status='unknown',access='sequential')
	      open(21,file='esit',status='unknown',access='sequential')
	      open(45,file='bgap',status='unknown',access='sequential')
	      open(46,file='pall',status='unknown',access='sequential')
	   endif
	   if(dotsa.or.new2s) then
	      open(11,file='popu2s',status='unknown',access='sequential')
	      open(12,file='esit2s',status='unknown',access='sequential')
	      open(13,file='tdat'  ,status='unknown',access='sequential')
	      open(14,file='hdha'  ,status='unknown',access='sequential')
	      open(15,file='homo'  ,status='unknown',access='sequential')
	   endif
	endif
	open(42,file='recd')
	open(44,file='ener')

	if(.not.restart) then
	   time0=0
	else
	   if(.not.dotsa) then
	      call readrst(ndim,time0,psi,dotsa)
	   else
	      call readrst(2,time0,psi,dotsa)
	   endif
	endif
	
	call hamiltonian(tau,nstp,ndim,tmd,h,cof,protein,dotsa,eta)
	
	IF(STATIC) THEN
	
           do i = 1 , ndim
              do j = 1 , ndim
	         hami(j,i)=h(j,i,1)
	         ovlp(j,i)=0.0
	      enddo
	      ovlp(i,i)=1.0
	   enddo
	   alle=.true.
	   ellow=-1000.0
	   elup=1000.0
           call eigen(ellow,elup,alle,ndim,hami,ovlp)
	   call output(ndim,hami,ovlp,tau,protein,neign)
	   call caltda(ndim,hami,etun,tda,hdd,haa)
	   if(quick) call faster(ndim,tau,tda)

           goto 896
 895	   write(iout,*) ' continue ? (y/n) '
	   read(iin,'(a)') yesorno
	   if(yesorno.eq.'y') then
	      goto 896
	   elseif(yesorno.eq.'n') then
	      stop
	   else
	      write(iout,*) ' answer [y]es or [n]o '
	      goto 895
	   endif
 896       continue
	   
	ELSE

c>>> do dynamics in the eigenspace

           if(eigenspace) then
              call distmap(ndim,protein)
              if(mod(ndim-nlp,2).ne.0) stop 'error with lonepairs!'
              neign=(ndim+nlp)/2
 	      call interpol(nstp,ndim,tmd,h,cof)
              call kmatrix(nstp,mstp,ndim,tmd,h,cof,
     :                     teg,ft,eg,ev,remft,neign)
     
c turn some part of kmatrix off
              if(turnoff) then
                 do m = 1 , mstp
                    do i = 1 , ndim
                       if(i.eq.neign.or.i.eq.neign-1) goto 3434
              	       ft(neign-1,i,m)=0.0
              	       ft(i,neign-1,m)=0.0
              	       ft(neign,i,m)=0.0
              	       ft(i,neign,m)=0.0
 3434          	       continue
                    enddo
                 enddo
              endif

              call dyneigen(tau,mstp,ndim,teg,eg,ev,ft,neign,
     :                      hami,ovlp,cofft)

              stop

           endif
	
c>>> check the bridge eigenstates

           if(brgap) then

              call distmap(ndim,protein)
              if(mod(ndim-nlp,2).ne.0) stop 'error with lonepairs!'
              if(both) then
                 neign=(ndim+nlp)/2
              else
                 if(hole) then
                    neign=ndim-1
                 else
                    neign=1
                 endif
              endif
        
              do n = 1 , nstp
                 if(.not.restart) then
                    do i = 1 , ndim
                       do j = 1 , ndim
                          hami(j,i)=h(j,i,n)
                          ovlp(j,i)=0.0
                       enddo
                       ovlp(i,i)=1.0
                    enddo
	            call bspace(ndim,hami,ovlp,hnew)
	            write(iout,*) n
	            write(45,'(50f10.5)') 
     *  	         tmd(n)/1.51653,(omcm(i),i=neign-2,neign-1)
	         else
	            if(n.ge.3) then
                       do i = 1 , ndim
                          do j = 1 , ndim
                             hami(j,i)=h(j,i,n)
                             ovlp(j,i)=0.0
                          enddo
                          ovlp(i,i)=1.0
                       enddo
	               call bspace(ndim,hami,ovlp,hnew)
	               write(iout,*) n
	               write(45,'(50f10.5)') 
     *  	            tmd(n)/1.51653,(omcm(i),i=neign-2,neign-1)
	            endif
	         endif
	      enddo
	      close(45)
	   endif

           if(new2s) then

              if(alle) then
                 call distmap(ndim,protein)
                 if(mod(ndim-nlp,2).ne.0) stop 'error with lonepairs!'
                 neign=(ndim+nlp)/2
                 ellow=-1000.0
                 elup=1000.0
              else
                 neign=2
                 ellow=-10.0
                 elup=0.0
              endif
              
              do n = 1 , nstp

                 do i = 1 , ndim
                    do j = 1 , ndim
                       hami(j,i)=h(j,i,n)
                       ovlp(j,i)=0.0
                    enddo
                    ovlp(i,i)=1.0
                 enddo
                 
                 call eigen(ellow,elup,alle,ndim,hami,ovlp)

                 ss1=dsqrt(zr(idon,neign-1)**2+zr(iacc,neign-1)**2)
                 ss2=dsqrt(zr(idon,neign  )**2+zr(iacc,neign  )**2)
                 if(zr(idon,neign-1).gt.0.d0) then
                    voso(1,1) = zr(idon,neign-1)/ss1
                    voso(2,1) = zr(iacc,neign-1)/ss1
                 else
                    voso(1,1) = -zr(idon,neign-1)/ss1
                    voso(2,1) = -zr(iacc,neign-1)/ss1
                 endif

                 sumover1=zr(idon,neign-1)**2*zr(iacc,neign-1)**2
     *                   +zr(idon,neign  )**2*zr(iacc,neign  )**2
                 sumover2=zr(idon,neign-1)**2*zr(iacc,neign-1)**2/ss1**4
     *                   +zr(idon,neign  )**2*zr(iacc,neign  )**2/ss2**4
                 
                 voso(1,2)= voso(2,1)
                 voso(2,2)=-voso(1,1)
                 tem2(1,1)=voso(1,1)
                 tem2(1,2)=voso(1,2)
                 tem2(2,1)=voso(2,1)
                 tem2(2,2)=voso(2,2)
                 write(6,*)
                 write(6,'(2f20.10)')((tem2(i,j),j=1,2),i=1,2)
                 
                 heig(1,1)=omcm(neign-1)
                 heig(2,2)=omcm(neign)
                 heig(1,2)=0.0
                 heig(2,1)=0.0
                 
                 call matinv(2,tem2,sovo)
                 call mprod(2,voso,heig,tem1)
                 call mprod(2,tem1,sovo,tem2)

                 h2s(1,1,n)=tem2(1,1)
                 h2s(1,2,n)=tem2(1,2)
                 h2s(2,1,n)=tem2(2,1)
                 h2s(2,2,n)=tem2(2,2)
                 
                 write(6,'(i10)') n
                 write(6,'(2f20.10)')((tem2(i,j),j=1,2),i=1,2)
                 
	         if(.not.restart) then
	            write(12,'(4f20.10)') tmd(n)/1.51653,
     *	            sumover1,sumover1/sumover2
	            write(13,'(4f20.10)') tmd(n)/1.51653,
     *	            h2s(1,2,n),h2s(1,1,n)-h2s(2,2,n)
	            write(14,'(9f20.10)') tmd(n)/1.51653,
     * 	            h2s(1,1,n),h2s(2,2,n)
                    write(15,'(9f20.10)') tmd(n)/1.51653,
     *              omcm(neign-2),omcm(neign+1)
     	         else
	            if(n.ge.3) then
	            write(12,'(4f20.10)') tmd(n)/1.51653,
     *	            sumover1,sumover1/sumover2
	            write(13,'(4f20.10)') tmd(n)/1.51653,
     *	            h2s(1,2,n),h2s(1,1,n)-h2s(2,2,n)
	            write(14,'(9f20.10)') tmd(n)/1.51653,
     * 	            h2s(1,1,n),h2s(2,2,n)
                    write(15,'(9f20.10)') tmd(n)/1.51653,
     *              omcm(neign-2),omcm(neign+1)
	            endif
	         endif
	         
	      enddo

	      call tsa(tau,nstp,tmd,h2s)
	      stop

	   endif



	   
           if(dotsa.and..not.new2s) then
	      sum=0.0
              do n = 1 , nstp
                 do i = 1 , ndim
                    do j = 1 , ndim
                       hami(j,i)=h(j,i,n)
                    enddo
                 enddo
                 tda=0.0
                 hdd=0.0
                 haa=0.0
	         call caltda(ndim,hami,etun,tda,hdd,haa)
	         if(idon.lt.iacc) then
	            h2s(1,1,n)=hdd
	            h2s(2,2,n)=haa
	            h2s(1,2,n)=tda
	            h2s(2,1,n)=tda
	         else
	            h2s(1,1,n)=haa
	            h2s(2,2,n)=hdd
	            h2s(1,2,n)=tda
	            h2s(2,1,n)=tda
	         endif
	         if(.not.restart) then
	            write(13,'(4f20.10)') tmd(n)/1.51653,
     *	            tda,h2s(1,1,n)-h2s(2,2,n)
	            write(14,'(9f20.10)') tmd(n)/1.51653,
     * 	            h2s(1,1,n),h2s(2,2,n),
     *              h2s(1,1,n)-hami(idon,idon),
     *              h2s(2,2,n)-hami(iacc,iacc)
	         else
	            if(n.ge.3) then
	            write(13,'(4f20.10)') tmd(n)/1.51653,
     *	            tda,h2s(1,1,n)-h2s(2,2,n)
	            write(14,'(9f20.10)') tmd(n)/1.51653,
     * 	            h2s(1,1,n),h2s(2,2,n),
     *              h2s(1,1,n)-hami(idon,idon),
     *              h2s(2,2,n)-hami(iacc,iacc)
	            endif
	         endif
	      enddo
	      close(13)

	      call tsa(tau,nstp,tmd,h2s)

	   endif
	   
	   if(.not.dotsa.and..not.new2s) then
 	      call interpol(nstp,ndim,tmd,h,cof)
	      call dyn(nstp,tau,ndim,tmd,h,cof)
	   endif

        ENDIF

	
 945    close(11)
        close(12)
        close(13)
        close(14)
        close(15)
        close(21)
	close(44)
	close(46)

        STOP
	END