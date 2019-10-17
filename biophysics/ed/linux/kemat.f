        subroutine kemat(nstp,ndim,tmd,h,cof)
	implicit real*8 (a-h,o-z)
	parameter (neigd=248)
	real*4 tmd(nstp),h(ndim,ndim,nstp),cof(ndim,ndim,nstp,3)
	real*4 ft(ndim,ndim,nstp),dft(ndim,ndim,nstp)
	real*4 ftemp(ndim,ndim,nstp)
        dimension cik(nstp),cikd1(nstp),tim1(nstp)
	dimension hami(ndim,ndim),ovlp(ndim,ndim)
	dimension hnew(ndim,ndim),hoto(ndim,ndim)
	dimension msign(ndim,2)
	logical alle,restart
	common /eigv/ omcm(neigd),zr(neigd,neigd),zi(neigd,neigd),
     :                ne,nblw
        common/ida/time0,edon,eacc,idon,iacc,i0,j0,restart,nstep
        common /io/ iin,iout
	common/orbital/nbd,npb,nlp

        delta=(tmd(2)-tmd(1))/nstp*0.1

        open(37,file='friceign')
        open(38,file='diag')
        open(39,file='eigensign')
        open(34,file='dotprod')
        
        alle=.true.
	ellow=-1000.0
	elup=1000.0
	
        do n = 1 , nstp

           tim1(n)=n*delta
           do i = 1 , ndim
              do j = 1 , ndim
                 hami(j,i)=h(j,i,1)+n*delta*
     :                    (h(j,i,2)-h(j,i,1))/(tmd(2)-tmd(1))
                 ovlp(j,i)=0.0
              enddo
              ovlp(i,i)=1.0
           enddo
           call eigen(ellow,elup,alle,ndim,hami,ovlp)
                                 
           if(.not.restart.and.n.eq.1) then
              do i = 1 , ndim
                 jpick=0
                 xpick=0.0
                 do j = 1 , ndim
                    if(dabs(zr(j,i)).gt.xpick) then
                       xpick=dabs(zr(j,i))
                       jpick=j
                    endif
                 enddo
                 msign(i,1)=jpick
                 if(zr(jpick,i).ge.0) then
                    msign(i,2)=1
                 else
                    msign(i,2)=-1
                 endif
                 write(39,'(2i10,f8.4)') msign(i,1),msign(i,2),
     :                                   zr(jpick,i)                    
              enddo
              close(39)
           else
              open(39,file='eigensign',status='old')
              do i = 1 , ndim
                 read(39,'(2i10)') msign(i,1),msign(i,2)
              enddo
              close(39)
              do i = 1 , ndim
                 isign=1
                 if(msign(i,2)*zr(msign(i,1),i).lt.0.0) isign=-1
                 do j = 1 , ndim
                    zr(j,i)=zr(j,i)*isign
                 enddo
                 do k = 1 , ndim
                    dotprod=0.0
                    do j = 1 , ndim
                       dotprod=dotprod+zr(j,i)*ft(j,k,1)
                    enddo
                    if(i.eq.99)
     :                 write(34,'(i10,f10.5)') k,dotprod
                 enddo
              enddo
              write(34,*) n
           endif

           do i = 1 , ndim
              do j = 1 , ndim
                 ft(i,j,n)=zr(i,j)
              enddo
           enddo
                 
        enddo
              
        close(34)
              
        if(mod(ndim,2).eq.0) then
           nhaf=0.5*ndim
        else
           nhaf=(ndim+1)*0.5
        endif

        goto 4343

c...this method does not guarantee that the k matrix is antisymmetric
c...should be deprecated in a formal calculation.
c...however, this method is a good test for how small
c...the time step should be in order to achieve correct derivatives
c...in the hammer-schiffer's approach one cannot know wether the
c...derivatives are true or false. once the derivatives are 
c...miscalculated, the antisymmetry of the k-matrix cannot be
c...obtained.
c...helpful when determining the time step for the eigen dynamics.

        do k = 1 , ndim
           do i = 1 , ndim
              do ntime = 1 , nstp
                 cik(ntime)=ft(i,k,ntime)
              enddo
              call deriv1(tim1,cik,cikd1,nstp)
              do ntime = 1 , nstp
                 dft(i,k,ntime)=cikd1(ntime)
              enddo
           enddo
        enddo

        do n = 1 , nstp
           do i = 1 , ndim
              do j = 1 , ndim
                 hami(j,i)= ft(i,j,n)
                 ovlp(i,j)=dft(i,j,n)
              enddo
           enddo
           call strassen(ndim,nhaf,hami,ovlp,hnew)
           write(6,*) n
           write(37,*) n
           do i = 1 , ndim
              do j = 1 , ndim
                 ft(j,i,n)=hnew(j,i)
              enddo
              write(37,'(500f8.5)')(ft(j,i,n),j=1,ndim)
              write(38,'(i9,f10.6)') i,ft(i,i,n)
           enddo
        enddo
              
        goto 4344
              
c...hammes-schiffer and tully's formula for computing the k-matrix
c...this formula guarantees antisymmetry of the k-matrix automatically.
c...however, as antisymmetry in this case has nth to do with the
c...correctness of the calculation, one has no easy rule to judge
c...the calculation of the derivative. it is a better idea to perform
c...the above scheme first, decide the time step and then use this
c...scheme to ensure the antisymmetry and then the hermitianity of
c...the hamiltonian.
c...this scheme doesn't compute the endpoints, to computer the
c...k-matrix at the two endpoints, we use the first scheme and
c...antisymmetrize the k-matrix to eliminate the small errors.
c...when the time step is small enough, the two approaches give the
c...same results.


 4343   continue

        do k = 1 , ndim
           do i = 1 , ndim
              do ntime = 1 , nstp
                 cik(ntime)=ft(i,k,ntime)
              enddo
              call deriv1(tim1,cik,cikd1,nstp)
              do ntime = 1 , nstp
                 dft(i,k,ntime)=cikd1(ntime)
              enddo
           enddo
        enddo

        do n = 1 , nstp
           if(n.gt.1.and.n.lt.nstp) then
              do i = 1 , ndim
                 do j = 1 , ndim
                    hami(j,i)=ft(i,j,n-1)
                    ovlp(i,j)=ft(i,j,n+1)
                 enddo
              enddo
              call strassen(ndim,nhaf,hami,ovlp,hnew)
              do i = 1 , ndim
                 do j = 1 , ndim
                    hami(j,i)=ft(i,j,n+1)
                    ovlp(i,j)=ft(i,j,n-1)
                 enddo
              enddo
              call strassen(ndim,nhaf,hami,ovlp,hoto)
              do i = 1 , ndim
                 do j = 1 , ndim
                    ftemp(j,i,n)=0.25*(hnew(j,i)-hoto(j,i))/delta
                 enddo
              enddo
           else
              do i = 1 , ndim
                 do j = 1 , ndim
                    hami(j,i)= ft(i,j,n)
                    ovlp(i,j)=dft(i,j,n)
                 enddo
              enddo
              call strassen(ndim,nhaf,hami,ovlp,hnew)
c...this procedure is to antisymmetrize the k-matrix at the endpoints
              do i = 1 , ndim
                 do j = i , ndim
                    ftemp(j,i,n)= hnew(j,i)
                    ftemp(i,j,n)=-hnew(j,i)                          
                    if(i.eq.j) ftemp(i,j,n)=0.0
                 enddo
              enddo
           endif
        enddo

        do n = 1 , nstp
           write(6,*) n
           write(37,*) n
           do i = 1 , ndim
              do j = 1 , ndim
                 ft(j,i,n)=ftemp(j,i,n)
              enddo
              write(37,'(500f8.5)')(ft(j,i,n),j=1,ndim)
              write(38,'(i9,f10.6)') i,ft(i,i,n)
           enddo
        enddo

 4344   continue              

        close(37)
        close(38)
	
        return
        end