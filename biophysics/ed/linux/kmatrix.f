        subroutine kmatrix(nstp,mstp,ndim,tmd,h,cof,
     :                     teg,ft,eg,ev,remft,neigen)
	implicit real*8 (a-h,o-z)
	parameter (neigd=248)
	real*4 tmd(nstp),h(ndim,ndim,nstp),cof(ndim,ndim,nstp,3)
	real*4 teg(mstp)
	real*4 ft(ndim,ndim,mstp),eg(ndim,mstp),ev(ndim,ndim,mstp)
	dimension hami(ndim,ndim),ovlp(ndim,ndim)
	dimension hold(ndim,ndim),hnew(ndim,ndim)
	dimension msign(ndim,2)
	logical alle,restart,forcesign,remft
	common /eigv/ omcm(neigd),zr(neigd,neigd),zi(neigd,neigd),
     :                ne,nblw
        common/ida/time0,edon,eacc,idon,iacc,i0,j0,restart,nstep
        common /io/ iin,iout
	common/orbital/nbd,npb,nlp

	do j = 1 , ndim
	   do i = 1 , ndim
	      hami(i,j)=0.d0
	      ovlp(i,j)=0.d0
	      hold(i,j)=0.d0
	      hnew(i,j)=0.d0
	   enddo
	   msign(j,1)=0
	   msign(j,2)=0
	enddo
	
        delta=(tmd(2)-tmd(1))*0.00001d0
        stepl=(tmd(nstp)-tmd(1))/(mstp-1)
        
        do m = 1 , mstp
           teg(m)=stepl*(m-1)+tmd(1)
        enddo

        open(37,file='friceign')
        
        alle=.true.
        forcesign=.true.
	ellow=-1000.0
	elup=1000.0
	ischeme=1
	
        if(mod(ndim,2).eq.0) then
           nhaf=0.5*ndim
        else
           nhaf=(ndim+1)*0.5
        endif

        do n = 1 , mstp
        
           itt=(teg(n)-tmd(1))/(tmd(2)-tmd(1))+1
           dtt=teg(n)-tmd(itt)
           
           write(6,*) n,itt,dtt

           do i = 1 , ndim
              do j = 1 , ndim
                hami(j,i)=h(j,i,itt)+
     :          ((cof(j,i,itt,3)*dtt+cof(j,i,itt,2))*dtt
     :           +cof(j,i,itt,1))*dtt
                ovlp(j,i)=0.0
              enddo
              ovlp(i,i)=1.0
           enddo

c           do i = 1 , ndim
c              do j = 1 , ndim
c                hami(j,i)=h(j,i,n)
c                ovlp(j,i)=0.0
c              enddo
c              ovlp(i,i)=1.0
c           enddo

           call eigen(ellow,elup,alle,ndim,hami,ovlp)
           if(forcesign.and.n.eq.1) then
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
              enddo
           else
              do i = 1 , ndim
                 isign=1
                 if(msign(i,2)*zr(msign(i,1),i).lt.0.0) isign=-1
                 do j = 1 , ndim
                    zr(j,i)=zr(j,i)*isign
                 enddo
              enddo
           endif
           do i = 1 , ndim
              do j = 1 , ndim
                 ev(j,i,n)=zr(j,i)
              enddo
              eg(i,n)=omcm(i)
           enddo
           if(.not.restart) then
             write(21,'(10f8.4)') teg(n)/1.51653,
     :                           eg(neigen-2,n),eg(neigen-1,n),
     :                           eg(neigen,n),eg(neigen+1,n)
           else
             if(n.ge.3) then
               write(21,'(10f8.4)') teg(n)/1.51653,
     :                             eg(neigen-2,n),eg(neigen-1,n),
     :                             eg(neigen,n),eg(neigen+1,n)
             endif
           endif
           
           f1=zr(idon,neigen-1)**2+zr(idon,neigen)**2
     :       +zr(iacc,neigen-1)**2+zr(iacc,neigen)**2
     
           f2=0.0
           do i = 1 , neigen-2
              f2=f2+zr(idon,i)**2+zr(iacc,i)**2
           enddo
           do i = neigen+1 , ndim
              f2=f2+zr(idon,i)**2+zr(iacc,i)**2
           enddo
           write(45,'(10f8.4)') teg(n)/1.51653,f1,f2
           
        enddo

        if(remft) goto 4030

        do n = 1 , mstp
        
           write(6,*) n

           if(ischeme.eq.1) then

             itt=(teg(n)-tmd(1))/(tmd(2)-tmd(1))+1
             dtt=teg(n)-tmd(itt)+delta

             do i = 1 , ndim
               do j = 1 , ndim
                 hami(j,i)=h(j,i,itt)+
     :           ((cof(j,i,itt,3)*dtt+cof(j,i,itt,2))*dtt
     :            +cof(j,i,itt,1))*dtt
                 ovlp(j,i)=0.0
               enddo
               ovlp(i,i)=1.0
             enddo
             call eigen(ellow,elup,alle,ndim,hami,ovlp)
             if(forcesign) then
               do i = 1 , ndim
                 isign=1
                 if(msign(i,2)*zr(msign(i,1),i).lt.0.0) isign=-1
                 do j = 1 , ndim
                    zr(j,i)=zr(j,i)*isign
                 enddo
               enddo
             endif
             do i = 1 , ndim
               do j = 1 , ndim
                 hami(j,i)=ev(i,j,n)
                 ovlp(i,j)=zr(i,j)
               enddo
             enddo
             call strassen(ndim,nhaf,hami,ovlp,hold)
             do i = 1 , ndim
               do j = 1 , ndim
                 hami(j,i)=zr(i,j)
                 ovlp(i,j)=ev(i,j,n)
               enddo
             enddo
             call strassen(ndim,nhaf,hami,ovlp,hnew)
             do i = 1 , ndim
               do j = 1 , ndim
                 ft(j,i,n)=0.5*(hold(j,i)-hnew(j,i))/delta
               enddo
             enddo
             
             write(6,'(2f15.6)') ft(neigen-1,neigen,n),
     :                           ft(neigen,neigen-1,n)             
             write(37,'(i10,f15.6)') n,ft(neigen-1,neigen,n)
             
           else

             if(n.eq.1) then
               do i = 1 , ndim
                 do j = 1 , ndim
                    hami(j,i)=ev(i,j,n)
                    ovlp(i,j)=ev(i,j,n+1)
                 enddo
               enddo
               call strassen(ndim,nhaf,hami,ovlp,hold)
               do i = 1 , ndim
                 do j = 1 , ndim
                    hami(j,i)=ev(i,j,n+1)
                    ovlp(i,j)=ev(i,j,n)
                 enddo
               enddo
               call strassen(ndim,nhaf,hami,ovlp,hnew)
               do i = 1 , ndim
                 do j = 1 , ndim
                   ft(j,i,n)=0.5*(hold(j,i)-hnew(j,i))/stepl
                 enddo
               enddo
             elseif(n.eq.mstp) then
               do i = 1 , ndim
                 do j = 1 , ndim
                    hami(j,i)=ev(i,j,n-1)
                    ovlp(i,j)=ev(i,j,n)
                 enddo
               enddo
               call strassen(ndim,nhaf,hami,ovlp,hold)
               do i = 1 , ndim
                 do j = 1 , ndim
                    hami(j,i)=ev(i,j,n)
                    ovlp(i,j)=ev(i,j,n-1)
                 enddo
               enddo
               call strassen(ndim,nhaf,hami,ovlp,hnew)
               do i = 1 , ndim
                 do j = 1 , ndim
                   ft(j,i,n)=0.5*(hold(j,i)-hnew(j,i))/stepl
                 enddo
               enddo
             else
               do i = 1 , ndim
                 do j = 1 , ndim
                    hami(j,i)=ev(i,j,n-1)
                    ovlp(i,j)=ev(i,j,n+1)
                 enddo
               enddo
               call strassen(ndim,nhaf,hami,ovlp,hold)
               do i = 1 , ndim
                 do j = 1 , ndim
                    hami(j,i)=ev(i,j,n+1)
                    ovlp(i,j)=ev(i,j,n-1)
                 enddo
               enddo
               call strassen(ndim,nhaf,hami,ovlp,hnew)
               do i = 1 , ndim
                 do j = 1 , ndim
                   ft(j,i,n)=0.25*(hold(j,i)-hnew(j,i))/stepl
                 enddo
               enddo
             endif
             
           endif

        enddo
        
        open(67,file='maxk')
        do n = 1 , mstp
           umax1=0.0
           umax2=0.0
           usum1=0.0
           usum2=0.0
           do i = 1 , ndim
              if(i.eq.neigen-1.or.i.eq.neigen) goto 911
              valtemp1=abs(ft(neigen-1,i,n)/(eg(neigen-1,n)-eg(i,n)))
              valtemp2=abs(ft(neigen,i,n)/(eg(neigen,n)-eg(i,n)))
              usum1=usum1+valtemp1
              usum2=usum2+valtemp2
              if(valtemp1.gt.umax1) then
                 umax1=valtemp1
                 imax1=i
              endif
              if(valtemp2.gt.umax2) then
                 umax2=valtemp2
                 imax2=i
              endif
 911          continue
           enddo
           write(67,'(i8,4f10.5)') n,umax1,umax2,usum1,usum2
        enddo
        close(67)

        goto 4040
           
 4030   continue
        do n = 1 , mstp
          do i = 1 , ndim
            do j = 1 , ndim
               ft(j,i,n)=0.0
            enddo
          enddo
        enddo
           
 4040   continue
                 

        close(37)
	
        return
        end