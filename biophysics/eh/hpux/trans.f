	subroutine trans(ndim,t,n,st,ht,aring,nring)
	
c********************************************************************
c transform overlap, hamiltonian, and density matrices to the NBO
c basis
c********************************************************************	
	
	implicit real*8 (a-h,o-z)
	parameter(nsmp=10,natmax=200,ndmax=500,neigd=500)
        logical invs,alle,eonly,check,write,fixeh,restart,dotsa,aring
	dimension t(ndim,ndim),t1(ndim,ndim)
	character bname*2,ac*4,symbol*2
        character*1 transa,transb
	real*4 dm,pop,sao,trt
	real*4 st(ndmax,ndmax,nsmp),ht(ndmax,ndmax,nsmp)
        common/ato/ac(natmax),symbol(40)
	common /eigt/ dm(ndmax,ndmax,nsmp),pop(ndmax,nsmp)
        common /inst/ ovlp(ndmax,ndmax),hami(ndmax,ndmax),
     *                dens(ndmax,ndmax)
        common/lbl/bname(ndmax,2),label(ndmax,3),ibx(ndmax)
        common/lblaro/labar(ndmax,7)
        common/eigv/omcm(neigd),zr(neigd,neigd),zi(neigd,neigd),ne,nblw
        common/secl/s(neigd,neigd),h(neigd,neigd),nmat
        common/iounit/iin,iout,istart(40),iungf,iunhl
        common/nonadi/sao(ndmax,ndmax,nsmp),trt(ndmax,ndmax,nsmp)
        common/option/restart,dotsa,nstep,unit
        dimension stold(ndim,ndim),stnew(ndim,ndim)
        dimension htold(ndim,ndim),htnew(ndim,ndim)
        dimension dmold(ndim,ndim),dmnew(ndim,ndim)
        dimension temp(ndim,ndim),rho(ndim,ndim)
        dimension ebond(30),inaro(30),ithread(10),movep(10)

        alpha=1.d0
        beta=0.d0
        lda=ndim
        ldb=ndim
        ldc=ndim
        transa='N'
        transb='N'
        
        check=.false.
        write=.false.
        fixeh=.true.
        if(mod(ndim,2).eq.0) then
           nhaf=ndim*0.5
        else
           nhaf=(ndim+1)*0.5
        endif
        judge=0
        
        if(aring) then
           itemp=0
           do i = 1 , ndim
              do ib = 1 , ndim
                 if(ibx(ib).eq.i) goto 70
              enddo
 70           if(bname(ib,1).eq.'PB'.and.bname(ib,2)(1:1).eq.'*') then
                 itemp=itemp+1
                 inaro(itemp)=i
              endif
           enddo
           itemq=0
           do ir = 1 , itemp
              if(mod(ir-1,3).eq.0) then
                 itemq=itemq+1
                 ithread(itemq)=inaro(ir)
              endif
           enddo
           if(itemq.ne.nring) stop 'error'
        endif
        
        do j = 1 , ndim
           do i = 1 , ndim
              dmold(i,j)=dm(i,j,n)
              stold(i,j)=st(i,j,n)
              htold(i,j)=ht(i,j,n)
              sao(i,j,n)=st(i,j,n)
           enddo
        enddo
        do im = 1 , nring
           movep(im)=2
        enddo
        
 2000   do i = 1 , ndim
           do j = 1 , ndim
              t1(i,j)=t(j,i)
           enddo
        enddo
c        call mprod(ndim,t1,dmold,temp)
c        call mprod(ndim,temp,t ,dmnew)
c        call mprod(ndim,t1,stold,temp)
c        call mprod(ndim,temp,t ,stnew)
c        call mprod(ndim,t1,htold,temp)
c        call mprod(ndim,temp,t ,htnew)
c        call strassen(ndim,nhaf,t1,dmold,temp)
c        call strassen(ndim,nhaf,temp,t ,dmnew)
c        call strassen(ndim,nhaf,t1,stold,temp)
c        call strassen(ndim,nhaf,temp,t ,stnew)
c        call strassen(ndim,nhaf,t1,htold,temp)
c        call strassen(ndim,nhaf,temp,t ,htnew)
        call dgemm(transa,transb,ndim,ndim,ndim,alpha,
     *              t1,lda,dmold,ldb,beta,temp,ldc)         
        call dgemm(transa,transb,ndim,ndim,ndim,alpha,
     *              temp,lda,t,ldb,beta,dmnew,ldc)         
        call dgemm(transa,transb,ndim,ndim,ndim,alpha,
     *              t1,lda,stold,ldb,beta,temp,ldc)         
        call dgemm(transa,transb,ndim,ndim,ndim,alpha,
     *              temp,lda,t,ldb,beta,stnew,ldc)         
        call dgemm(transa,transb,ndim,ndim,ndim,alpha,
     *              t1,lda,htold,ldb,beta,temp,ldc)         
        call dgemm(transa,transb,ndim,ndim,ndim,alpha,
     *              temp,lda,t,ldb,beta,htnew,ldc)         

        if(judge.eq.0) then

c           call mprod(ndim,dmnew,stnew,rho)
c           call strassen(ndim,nhaf,dmnew,stnew,rho)
           call dgemm(transa,transb,ndim,ndim,ndim,alpha,
     *                 dmnew,lda,stnew,ldb,beta,rho,ldc)         

           if(write) then
	      write(iout,*) 'Density matrix'
              do i = 1 , ndim
                 write(iout,'(500f8.3)')(rho(i,j),j=1,ndim)
              enddo
           endif
           if(aring) then
              do iring = 1 , nring
                 rmax=0.0
                 do i = ithread(iring),ithread(iring)+2
                    if(rho(i,i).gt.rmax) rmax=rho(i,i)
                 enddo
                 if(rmax.gt.0.5) then
                    if(movep(iring).ge.4) then
                       write(*,*) iring
                       stop 'fail in forming aromatic pi bond'
                    endif
                    movep(iring)=movep(iring)+1
                    mtemp=movep(iring)
                    write(iout,*) 'rectify',iring,mtemp
                    call autoget(iring,mtemp)
                    call reformaro(t,ndim,nring)
                    goto 2000
                 endif
              enddo
           endif
           
        endif
        
        if(judge.eq.0) then        
           do ibd = 1 , ndim
              sum=0.0
              do j = 1 , ndim
                 t(j,ibd)=t(j,ibd)/dsqrt(stnew(ibd,ibd))
                 sum=sum+t(j,ibd)*t(j,ibd)
              enddo
           enddo
           judge=1
           goto 2000
        endif
        
        if(write) then
           write(iout,*)
	   write(iout,*) 'Hamiltonian matrix'
           do i = 1 , ndim
              write(iout,'(500f8.3)')(htnew(i,j),j=1,ndim)
           enddo
           write(iout,*)
	   write(iout,*) 'Overlap matrix'
           do i = 1 , ndim
              write(iout,'(500f8.3)')(stnew(i,j),j=1,ndim)
           enddo
        endif

        rho_anti_max=0.0
        rho_bond_min=2.0
        ibond=0
        ianti=0
        write(iout,*)'diagonal elements in NBO basis'
        do i = 1 , ndim
           iflag=1
           do ib = 1 , ndim
              if(ibx(ib).eq.i) goto 90
           enddo
 90        if(bname(ib,1).eq.'LP') nctr=1
           if(bname(ib,1).eq.'BD') nctr=2
           if(bname(ib,1).eq.'PB') nctr=6
c fix C=O and aromatic pi antibonds......           
           if(bname(ib,2)(1:1).eq.'*') then
              if(fixeh) then
                 if(bname(ib,1).eq.'PB'.or.
     :             (bname(ib,1).eq.'BD'
     :             .and.htnew(i,i).lt.0.and.
     :             ((ac(label(ib,2))(1:2).eq.' C'.and.
     :               ac(label(ib,3))(1:2).eq.' O').or.
     :              (ac(label(ib,2))(1:2).eq.' O'.and.
     :               ac(label(ib,3))(1:2).eq.' C').or.
     :              (ac(label(ib,2))(1:2).eq.' O'.and.
     :               ac(label(ib,3))(1:2).eq.' H')))
     :              ) then
                    if(ac(label(ib,2))(1:4).eq.' CB '.and.
     :                 ac(label(ib,3))(1:4).eq.' OG1') goto 1234
                    iflag=-1
                    htorg=htnew(i,i)
                    htnew(i,i)=htnew(i,i)+12.0
 1234               continue
                 endif
              endif
           endif
c fix aromatic C=C antibonds......
           if(fixeh.and.aring) then
              if(bname(ib,1).eq.'BD'.and.
     :   ((ac(label(ib,2)).eq.' CG '.and.ac(label(ib,3)).eq.' CD1')
     :.or.(ac(label(ib,2)).eq.' CG '.and.ac(label(ib,3)).eq.' CD2')
     :.or.(ac(label(ib,2)).eq.' CD1'.and.ac(label(ib,3)).eq.' CE1')
     :.or.(ac(label(ib,2)).eq.' CD2'.and.ac(label(ib,3)).eq.' CE2')
     :.or.(ac(label(ib,2)).eq.' CE1'.and.ac(label(ib,3)).eq.' CZ ')
     :.or.(ac(label(ib,2)).eq.' CE2'.and.ac(label(ib,3)).eq.' CZ '))
     :              ) then
                 if(bname(ib,2)(1:1).ne.'*') then
                    ibond=ibond+1
                    ebond(ibond)=htnew(i,i)
                 endif
                 if(bname(ib,2)(1:1).eq.'*') then
                    iflag=-1
                    ianti=ianti+1
                    htorg=htnew(i,i)
                    htnew(i,i)=ebond(ianti)+40.0
                 endif
              endif
           endif
           if(iflag.eq.1) then
              if(nctr.eq.6) then
              write(iout,'(i8,5x,a2,a1,5x,2f10.4,6(a2,i3))')
     :           i,bname(ib,1),bname(ib,2)(1:1),rho(i,i),htnew(i,i),
     :           (ac(labar(ib,nb+1))(1:2),labar(ib,nb+1),nb=1,nctr)
              else
              write(iout,'(i8,5x,a2,a1,5x,2f10.4,2(a2,i3))')
     :           i,bname(ib,1),bname(ib,2)(1:1),rho(i,i),htnew(i,i),
     :           (ac(label(ib,nb+1))(1:2),label(ib,nb+1),nb=1,nctr)
              endif
           else
              if(nctr.eq.6) then
              write(iout,'(i8,5x,a2,a1,5x,2f10.4,6(a2,i3),a3,f10.4)')
     :           i,bname(ib,1),bname(ib,2)(1:1),rho(i,i),htorg,
     :           (ac(labar(ib,nb+1))(1:2),labar(ib,nb+1),nb=1,nctr),
     :           '-->',htnew(i,i)
              else
              write(iout,'(i8,5x,a2,a1,5x,2f10.4,2(a2,i3),a3,f10.4)')
     :           i,bname(ib,1),bname(ib,2)(1:1),rho(i,i),htorg,
     :           (ac(label(ib,nb+1))(1:2),label(ib,nb+1),nb=1,nctr),
     :           '-->',htnew(i,i)
              endif
           endif
           if(bname(ib,2)(1:1).eq.'*') then
              if(rho(ib,ib).gt.rho_anti_max) then
                 rho_anti_max=rho(ib,ib)
                 inx_anti_max=ib
              endif
           else
              if(rho(ib,ib).lt.rho_bond_min) then
                 rho_bond_min=rho(ib,ib)
                 inx_bond_min=ib
              endif    
           endif
        enddo
        write(iout,*)'--------------------------------------------'
        write(iout,'(2x,a23,f8.4,2x,a7,i4)')'minimum bond occupancy=',
     :    rho_bond_min,'at bond',inx_bond_min
        write(iout,'(2x,a23,f8.4,2x,a7,i4)')'maximum anti occupancy=',
     :    rho_anti_max,'at bond',inx_anti_max
        write(iout,*)
        if(rho_anti_max.gt.0.6) stop 'antibond density too big'
        
        do i = 1 , ndim
           do j = 1 , ndim
              hami(i,j)=htnew(i,j)
              ovlp(i,j)=stnew(i,j)
              trt(i,j,n)=t(i,j)
           enddo
        enddo

        if(check) then
           do i = 1 , ndim
              do j = 1 , ndim
                 s(i,j)=0.0
                 if(i.eq.j) s(i,j)=1.0
                 h(i,j)=0.0
              enddo
           enddo
           do i = 1 , ndim
              do j = 1, i
                 h(j,i)=0.0
                 h(i,j)=htnew(i,j)
                 s(i,j)=stnew(i,j)
              enddo
           enddo
           nmat=ndim
           invs=.false.
           alle=.false.
           eonly=.false.
           ellow=-1000.d0
           elup=  1000.d0
           call diagm(ellow,elup,invs,alle,eonly)
           write(iout,*)'check transformation by rediagonalization'
           write(iout,'(i8,f8.4)')(i,omcm(i),i=1,ndim)
           stop
        endif
                        
        return
        end

	subroutine autoget(iring,movep)
	
c********************************************************************
c automatically find the correct orientations for aromatic ring
c********************************************************************	
	
	implicit real*8 (a-h,o-z)
	parameter(nsmp=10,natmax=200,ndmax=500,neigd=500)
	integer ul
        common/q/q(4,ndmax)
        common/info1/ul(natmax),ll(natmax),nele,nocc
        common/cring/idc6(6,10),ovar(6,6,10)
        common/iounit/iin,iout,istart(40),iungf,iunhl

        do iat = 1 , 6
           ilc=ll(idc6(iat,iring))
           iuc=ul(idc6(iat,iring))
           if(q(movep,iuc).lt.0.0) then
              do iorb = 1 , 4
                 q(iorb,iuc)=-q(iorb,iuc)
              enddo
           endif
        enddo
        
        return
        end