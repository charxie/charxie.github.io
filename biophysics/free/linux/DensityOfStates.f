        subroutine freeofT(epot)
        implicit real*8 (a-h,o-z)
        parameter(ndos=100,ntem=100)
        dimension omega(ndos),rho(ndos)
        dimension temp(ntem),free(ntem)
        
        open(11,file='dos.50',status='old')
        sum=0.d0
        do i = 1 , ndos
     	   read(11,'(2f20.10)') omega(i),rho(i)
     	   sum=sum+rho(i)
        enddo
        close(11)
        
        ave=0.d0
        omegadelta=omega(2)-omega(1)
        do i = 1 , ndos
           rho(i)=3.d0*rho(i)/sum/omegadelta
           ave=ave+rho(i)*omega(i)*omegadelta
        enddo
        
        zero=6.626d0/3.2d0*0.001d0*ave
        converter=1.38d0/1.6d0*0.0001d0
        
        do i = 1 , ntem
           temp(i)=i*300.d0+100.d0
           free(i)=epot-zero
           do n = 1 , ndos
              pass=6.626d0*10.d0*omega(n)/(2.d0*1.38d0*temp(i))
              free(i)=free(i)+temp(i)*converter*omegadelta
     *                       *dlog(dexp(pass)-dexp(-pass))*rho(n)
           enddo
        enddo        

        open(12,file='freeharm')
        do i = 1 , ntem
           write(12,'(2f20.10)') temp(i),free(i)
        enddo
        close(12)

        return
        end

        subroutine getdos(epsilon,sigma,rcut,cell)
        implicit real*8 (a-h,o-z)
        parameter(nmax=5000,ndos=100,kpnt=30)
        parameter(ndim=3)
        parameter (neigd=1500)
        logical invs,alle,eonly
        common/box/xbox,ybox,zbox
        common/block1/rx(nmax),ry(nmax),rz(nmax),
     :                vx(nmax),vy(nmax),vz(nmax),
     :                ax(nmax),ay(nmax),az(nmax),
     :                bx(nmax),by(nmax),bz(nmax),
     :                cx(nmax),cy(nmax),cz(nmax),
     :                fx(nmax),fy(nmax),fz(nmax)
        common/hessian/hess(ndim,ndim)
        common/eigv/omcm(neigd),zr(neigd,neigd),zi(neigd,neigd),ne,nblw
        common/secl/s(neigd,neigd),h(neigd,neigd),nmat
        dimension unit(ndim,ndim)
        dimension omega(ndos),rho(ndos)

        twopi=8.d0*datan(1.d0)
        nmat=ndim
        invs=.false.
        alle=.false.
        eonly=.false.
        ellow=-1000.d0
        elup=  1000.d0
        cellinv=twopi/cell/kpnt
        omegamax=500.d0
        omegadelta=omegamax/ndos

        open(91,file='disp')

        do n = 1 , ndos
           omega(n)=n*omegadelta
           rho(n)=0.0
        enddo

        do kx = 1 , kpnt
        do ky = 1 , 1
        do kz = 1 , 1
        
           write(6,'(3i10)') kx,ky,kz

           xv=(kx-1)*cellinv
           yv=(ky-1)*cellinv
           zv=(kz-1)*cellinv

           call buildHessian(epsilon,sigma,rcut,cell,xv,yv,zv)
        
           do m = 1 , ndim
              do n = 1 , ndim
                 unit(m,n)=0.d0
                 h(m,n)=0.d0
                 s(m,n)=0.d0
                 if(m.eq.n) unit(m,n)=1.d0
              enddo
           enddo
        
           do i = 1 , ndim
              do j = 1, i
                 h(j,i)=0.d0
                 h(i,j)=hess(i,j)
                 s(i,j)=unit(i,j)
              enddo
           enddo
      
           call diagm(ellow,elup,invs,alle,eonly)
        
           write(91,'(10f20.10)') xv,(omcm(i),i=1,ndim)
           do i = 1 , ndim
              nget=omcm(i)/omegadelta+1
              if(nget.gt.ndos+1) then
                 write(6,*) nget,omcm(i)
                 stop ' error'
              endif
              rho(nget)=rho(nget)+1.d0
           enddo
           
        enddo
        enddo
        enddo      

c>>> normalized to 3
        sum=0.d0
        do i = 1 , ndos
     	   sum=sum+rho(i)
        enddo
        do i = 1 , ndos
           rho(i)=3.d0*rho(i)/sum
           rho(i)=rho(i)/omegadelta
        enddo

        open(11,file='dos')

c>>> frequency in THz
        converter=dsqrt(1.6d0*6.022d0*10.d0/12.d0)
     	do n = 1 , ndos
     	   write(11,'(2f20.10)') omega(n)*converter,rho(n)
     	enddo

     	close(11)
     	close(91)


        return
        end


        subroutine buildHessian(epsilon,sigma,rcut,cell,
     *                          xvec,yvec,zvec)
        implicit real*8 (a-h,o-z)
        parameter(nmax=5000,ninc=500)
        parameter(ndim=3)
        parameter(amass=1.d0)
        common/box/xbox,ybox,zbox
        common/block1/rx(nmax),ry(nmax),rz(nmax),
     :                vx(nmax),vy(nmax),vz(nmax),
     :                ax(nmax),ay(nmax),az(nmax),
     :                bx(nmax),by(nmax),bz(nmax),
     :                cx(nmax),cy(nmax),cz(nmax),
     :                fx(nmax),fy(nmax),fz(nmax)
        dimension x0(ninc),y0(ninc),z0(ninc),r0(ninc)
        character type*3
        common/hessian/hess(ndim,ndim)
        real*8 lennardJones
        real*8 lennardJonesDerivative1
        real*8 lennardJonesDerivative2

        type='fcc'
        
        open(11,file=type)
        do i = 1 , ninc
           read(11,'(4f16.7)') x0(i),y0(i),z0(i),r0(i)
           x0(i)=x0(i)*cell
           y0(i)=y0(i)*cell
           z0(i)=z0(i)*cell
           r0(i)=r0(i)*cell
        enddo
        close(11)
        
        do m = 1 , 3
           do n = 1 , 3
              hess(m,n)=0.d0
           enddo
        enddo
        
        energy=0.d0

        do i = 1 , ninc

           if(r0(i).le.20.d0) then
              qr=xvec*x0(i)+yvec*y0(i)+zvec*z0(i)
              energy=energy+lennardJones(epsilon,sigma,r0(i))
              temp1=lennardJonesDerivative2(epsilon,sigma,r0(i))
     *             -lennardJonesDerivative1(epsilon,sigma,r0(i))
     *             /r0(i)
              temp2=lennardJonesDerivative1(epsilon,sigma,r0(i))
     *             /r0(i)
              factor=1.d0-dcos(qr)
              hess(1,1)=hess(1,1)+(temp1*x0(i)/r0(i)*x0(i)/r0(i)
     *                            +temp2/r0(i))*factor
              hess(2,2)=hess(2,2)+(temp1*y0(i)/r0(i)*y0(i)/r0(i)
     *                            +temp2/r0(i))*factor
              hess(3,3)=hess(3,3)+(temp1*z0(i)/r0(i)*z0(i)/r0(i)
     *                            +temp2/r0(i))*factor
              hess(1,2)=hess(1,2)+temp1*x0(i)/r0(i)*y0(i)/r0(i)*factor
              hess(1,3)=hess(1,3)+temp1*x0(i)/r0(i)*z0(i)/r0(i)*factor
              hess(2,3)=hess(2,3)+temp1*y0(i)/r0(i)*z0(i)/r0(i)*factor
           endif

        enddo
        
        hess(2,1)=hess(1,2)
        hess(3,1)=hess(1,3)
        hess(3,2)=hess(2,3)
        energy=0.5*energy
        
        do m = 1 , 3
           do n = 1 , 3
              hess(m,n)=hess(m,n)/amass
           enddo
        enddo

c        write(6,*) energy        
c        write(6,'(3e20.5)') ((hess(m,n),m=1,3),n=1,3)
        
        return
        end