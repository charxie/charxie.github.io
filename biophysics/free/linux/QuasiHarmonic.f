
        subroutine saveTraj(ntot,ioutR,ioutV)
        implicit real*8 (a-h,o-z)
        parameter(nmax=5000)
        common/block1/rx(nmax),ry(nmax),rz(nmax),
     :                vx(nmax),vy(nmax),vz(nmax),
     :                ax(nmax),ay(nmax),az(nmax),
     :                bx(nmax),by(nmax),bz(nmax),
     :                cx(nmax),cy(nmax),cz(nmax),
     :                fx(nmax),fy(nmax),fz(nmax)

        if(ioutR.ne.0) then
           do i = 1 , ntot
              write(ioutR) rx(i),ry(i),rz(i)
           enddo
        endif
        if(ioutV.ne.0) then
           do i = 1 , ntot
              write(ioutV) vx(i),vy(i),vz(i)
           enddo
        endif

        return
        end
        
        

        subroutine quasi(nprod,ntot,temperature)
        implicit real*8 (a-h,o-z)
        parameter(nmax=5000)
        parameter(neigd=1500)
        logical invs,alle,eonly
        common/box/xbox,ybox,zbox
        common/block1/rx(nmax),ry(nmax),rz(nmax),
     :                vx(nmax),vy(nmax),vz(nmax),
     :                ax(nmax),ay(nmax),az(nmax),
     :                bx(nmax),by(nmax),bz(nmax),
     :                cx(nmax),cy(nmax),cz(nmax),
     :                fx(nmax),fy(nmax),fz(nmax)
        common/eigv/omcm(neigd),zr(neigd,neigd),zi(neigd,neigd),ne,nblw
        common/secl/s(neigd,neigd),h(neigd,neigd),nmat
        dimension xave(ntot),yave(ntot),zave(ntot)
        dimension cov(3*ntot,3*ntot),unit(3*ntot,3*ntot)

        open(11,file='dcd',form='unformatted',status='old')

        ntwo=2*ntot
        nthr=3*ntot

        nmat=nthr
        invs=.false.
        alle=.false.
        eonly=.false.
        ellow=-1000.d0
        elup=  1000.d0
        
        converter=1.38d0/1.6d0*0.000000001d0

        do i = 1 , ntot
           xave(i)=0.d0
           yave(i)=0.d0
           zave(i)=0.d0
        enddo

        do j = 1 , nthr
           do i = 1 , nthr
              cov(i,j)=0.d0
           enddo
        enddo
        
        do n = 1 , nprod
           do i = 1 , ntot
              read(11) rx(i),ry(i),rz(i)
           enddo
           do i = 1 , ntot
              xave(i)=xave(i)+rx(i)
              yave(i)=yave(i)+ry(i)
              zave(i)=zave(i)+rz(i)
           enddo
        enddo
        do i = 1 , ntot
           xave(i)=xave(i)/nprod
           yave(i)=yave(i)/nprod
           zave(i)=zave(i)/nprod
        enddo

        rewind 11

        do n = 1 , nprod
        
           write(6,*) n
        
           do i = 1 , ntot
              read(11) rx(i),ry(i),rz(i)
           enddo
           do i = 1 , ntot
              rx(i)=rx(i)-xave(i)
              ry(i)=ry(i)-yave(i)
              rz(i)=rz(i)-zave(i)
           enddo
           do i = 1 , ntot
              do j = 1 , ntot
                 cov(i     ,j     )=rx(i)*rx(j)+cov(i     ,j     )
                 cov(i     ,j+ntot)=rx(i)*ry(j)+cov(i     ,j+ntot)
                 cov(i     ,j+ntwo)=rx(i)*rz(j)+cov(i     ,j+ntwo)
                 cov(i+ntot,j     )=ry(i)*rx(j)+cov(i+ntot,j     )
                 cov(i+ntot,j+ntot)=ry(i)*ry(j)+cov(i+ntot,j+ntot)
                 cov(i+ntot,j+ntwo)=ry(i)*rz(j)+cov(i+ntot,j+ntwo)
                 cov(i+ntwo,j     )=rz(i)*rx(j)+cov(i+ntwo,j     )
                 cov(i+ntwo,j+ntot)=rz(i)*ry(j)+cov(i+ntwo,j+ntot)
                 cov(i+ntwo,j+ntwo)=rz(i)*rz(j)+cov(i+ntwo,j+ntwo)
              enddo
           enddo
           
        enddo

        write(6,*) ' covariance matrix formed'

        do j = 1 , nthr
           do i = 1 , nthr
              cov(j,i)=cov(j,i)/dble(nprod)
           enddo
        enddo
        
        call matinv(nthr,cov,unit)
        do i = 1 , nthr
           do j = 1 , nthr
              cov(i,j)=unit(i,j)*converter
           enddo
        enddo
        
        write(6,*) ' covariance matrix inverted'

        open(21,file='cova')
        do i = 1 , nthr
           write(21,'(1000f20.5)') (cov(i,j),j=1,nthr)
        enddo
        close(21)

        close(11)

        do m = 1 , nthr
           do n = 1 , nthr
              unit(m,n)=0.d0
              h(m,n)=0.d0
              s(m,n)=0.d0
              if(m.eq.n) unit(m,n)=1.d0
           enddo
        enddo
        
        do i = 1 , nthr
           do j = 1, i
              h(j,i)=0.d0
              h(i,j)=cov(i,j)
              s(i,j)=unit(i,j)
           enddo
        enddo
      
        call diagm(ellow,elup,invs,alle,eonly)

        open(61,file='qdos')  
        do i = 1 , nthr      
           write(61,'(2e20.10)') i,omcm(i)
        enddo
        close(61)

        
        return
        end