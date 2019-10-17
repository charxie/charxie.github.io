c**************************************************************
c Program FREE
c
c This package is designed to investigate novel approaches
c for free energy calculations, based on the Lennard-Jones 
c model.
c
c Special effort is paid on checking the effectiveness of
c the competitive lambda dynamics
c**************************************************************

        implicit real*8 (a-h,o-z)
        parameter(nequi=500000,nprod=1000,nsave=10)
        parameter(nmax=5000,kmax=100000)
        parameter(lmax=20)
        parameter(amass=1.d0)
        common/blockLambda/rl(lmax),vl(lmax),
     :                     al(lmax),bl(lmax),
     :                     cl(lmax),fl(lmax),
     :                     amassLam(lmax)
        common/box/xbox,ybox,zbox
        common/block1/rx(nmax),ry(nmax),rz(nmax),
     :                vx(nmax),vy(nmax),vz(nmax),
     :                ax(nmax),ay(nmax),az(nmax),
     :                bx(nmax),by(nmax),bz(nmax),
     :                cx(nmax),cy(nmax),cz(nmax),
     :                fx(nmax),fy(nmax),fz(nmax)
        common/record/epot(kmax),ekin(kmax),etot(kmax),para(kmax)
        dimension RALambda(lmax)
        real*8 lennardJones
        character rescaled*1
        logical pure
        logical drawLJ,AST,searchGroundState,HARM,saveTrajectory
        logical QUASIHARM,TI,drawSF,Tsallis, DynLambda

        pure=.false.
        DynLambda=.true.
        drawLJ=.false.
        drawSF=.false.
        TI=.false.
        Tsallis=.false.
        AST=.false.
        HARM=.false.
        QUASIHARM=.false.
        searchGroundState=.false.
        saveTrajectory=.false.
        rescaled=' '
        ntype=4
        namb=5

        epsilon=1.0d0
        cell=2.d0
        sigma=cell/dsqrt(2.d0)*0.9155d0
        rcut=20.d0
        temperature=5000.d0
        tolerance=0.01d0*temperature
        delta=0.1d0
        nxc=5
        nyc=5
        nzc=5
        ts=nprod*delta
        
        if(DynLambda) then
           open(16,file='lambda-t')
           open(26,file='running')
        endif
        
        if(drawLJ) then
           open(11,file='poten')
           do i = 1 , 40
              rad=0.02*sigma*i+0.9*sigma
              write(11,*) rad,lennardJones(epsilon,sigma,rad)
           enddo
           close(11)
           stop
        endif
        if(drawSF) then
           open(11,file='switch')
           do i = 1 , 100
              t=i*ts/100.d0
              write(11,*) t,switchingFunction(ntype,ts,t)
           enddo
           close(11)
           stop
        endif
        if(searchGroundState) then
           call bain(epsilon,sigma,rcut)
c           call ground(epsilon,sigma,rcut)
           stop
        endif
        
        converter=2.d0*1.6d0/3.d0/1.38d0*10000.d0

        call lattice(ntot,nxc,nyc,nzc,cell,cell,cell)
        nsolute=ntot/2
        if(.not.pure) call solute(ntot,nsolute,epsilon,sigma)
        if(DynLambda) call initLambda(namb,epsilon,sigma)
        call assignVel(ntot,temperature,energyKinetic)
        if(DynLambda) then
          call forceLambda(ntot,nsolute,namb,rcut,energyPotential,vc,w)
        else
          call force(ntot,epsilon,sigma,rcut,energyPotential,vc,w,pure)
        endif
        write(6,*) ' box size'
        write(6,'(3f10.5)') xbox,ybox,zbox
        write(6,*) ' initial potential energy =',energyPotential,vc
        write(6,*) ' initial kinetic energy =',energyKinetic
        write(6,*) ' initial total energy =',
     *               energyKinetic+energyPotential
     
        getTemp=converter*energyKinetic
        write(6,*) ' desired temperature =',getTemp

        if(HARM) then
c           call getdos(epsilon,sigma,rcut,cell)
           call freeofT(energyPotential)
           stop
        endif

        if(QUASIHARM) then
           call quasi(nprod/10,ntot,temperature)
           stop
        endif

c>>> equilibration should attemp to drive the system
c    into a desired regime where the free energy can be
c    approximately calculated by other methods

        ratio=1.0d0
        if(TI) then
           if(AST) stop ' set AST flag false first!'
           ratio=0.3d0
        endif
        if(AST) then
           if(TI) stop ' set TI flag false first!'
        endif

c initialize the running averages        
        if(DynLambda) then
           do k = 1 , namb
              RALambda(k)=0.0
           enddo
        endif

        do istp = 1 , nequi
           call predictor(ntot,delta)
           if(DynLambda) then
              call predictorLambda(namb,delta)
           endif
           call setPBC(ntot)
           if(DynLambda) then
              call forceLambda(ntot,nsolute,namb,rcut,
     :                         epoten,potc,virial)
           else
              call force(ntot,epsilon,sigma,rcut,
     :                   epoten,potc,virial,pure)
           endif
           if(AST) then
              ratio=switchingFunction(ntype,ts,delta)
              do i = 1 , ntot
                 fx(i)=fx(i)*ratio
                 fy(i)=fy(i)*ratio
                 fz(i)=fz(i)*ratio
              enddo
           endif
           if(TI) then
              do i = 1 , ntot
                 fx(i)=fx(i)*ratio
                 fy(i)=fy(i)*ratio
                 fz(i)=fz(i)*ratio
              enddo
           endif
           call corrector(DynLambda,namb,nsolute,ntot,delta,ekinet)
           if(DynLambda) call correctorLambda(namb,delta)
           
c           if(Tsallis) then
c              epoten=q/(beta*(q-1.d0))*dlog(1.d0-(1.0-q)*beta*
c     :              (epoten+eref))
c           endif

           etotal=epoten+ekinet
           call setPBC(ntot)
           call order(ntot,cell,param)
           if(DynLambda) then
              sumLambda=0.d0
              sumVambda=0.d0
              do k = 1 , namb
                 sumLambda=sumLambda+rl(k)*rl(k)
                 sumVambda=sumVambda+rl(k)*vl(k)
              enddo
              if(mod(istp,10).eq.0) call rescaleLambda(namb)
              do k = 1 , namb
                 RALambda(k)=RALambda(k)+rl(k)*rl(k)
              enddo
           endif
           call rescaleVel(ntot,temperature,tolerance,ekinet,rescaled)
           if(mod(istp,10).eq.0) then
              write(6,'(i6,2f15.5,2x,f15.8,2x,f12.5,f8.0,a1)') 
     *              istp,epoten,ekinet,etotal,
     *              param,converter*ekinet,rescaled
              if(DynLambda) then
                 write(6,'(6x,20f8.3)')(rl(k),k=1,namb),
     *                 sumLambda,sumVambda
                 if(mod(istp,100).eq.0) then
                    write(16,'(i6,20f10.5)')
     *                   istp,(rl(k),k=1,namb)
                    write(26,'(i6,20f10.5)')
     *                   istp,(RALambda(k)/istp,k=1,namb)
                 endif
              endif
              rescaled=' '
           endif
        enddo
        
        if(DynLambda) stop
        
        rescaled=' '
        open(45,file='dcd',form='unformatted')

        parave=0.d0
        do istp = 1 , nprod
           call predictor(ntot,delta)
           call setPBC(ntot)
           call force(ntot,epsilon,sigma,rcut,epoten,potc,virial,pure)
           if(Tsallis) then
              t=istp*delta
              ratio=switchingFunction(ntype,ts,t)
              q=1.5d0
              beta=2.d0
              eref=8.d0
              ratio=ratio*q/(1.d0-(1.d0-q)*beta*(ratio*epoten+eref))
              do i = 1 , ntot
                 fx(i)=fx(i)*ratio
                 fy(i)=fy(i)*ratio
                 fz(i)=fz(i)*ratio
              enddo
           endif
           if(AST) then           
              t=istp*delta
              ratio=switchingFunction(ntype,ts,t)
              do i = 1 , ntot
                 fx(i)=fx(i)*ratio
                 fy(i)=fy(i)*ratio
                 fz(i)=fz(i)*ratio
              enddo
           endif
           if(TI) then
              do i = 1 , ntot
                 fx(i)=fx(i)*ratio
                 fy(i)=fy(i)*ratio
                 fz(i)=fz(i)*ratio
              enddo
           endif
           call corrector(DynLambda,namb,nsolute,ntot,delta,ekinet)
           if(saveTrajectory.and.mod(istp,nsave).eq.0) 
     *            call saveTraj(ntot,45,0)
           etotal=epoten+ekinet
           call setPBC(ntot)
           call order(ntot,cell,param)
           parave=parave+param
           if(mod(istp,nsave).eq.0) then
              epot(istp/nsave)=epoten
              ekin(istp/nsave)=ekinet
              etot(istp/nsave)=etotal
              para(istp/nsave)=param
           endif
           call rescaleVel(ntot,temperature,tolerance,ekinet,rescaled)
           if(mod(istp,10).eq.0) then
              write(6,'(i6,2f15.5,2x,f15.8,2x,f12.5,f8.0,a1)') 
     *              istp,epoten*ratio,ekinet,epoten*ratio+ekinet,
     *              param,converter*ekinet,rescaled
           endif
        enddo
        parave=parave/nprod
        
        write(6,*) ' order parameter =',parave

        open(12,file='fluc')
        averU=0.d0
        do k = 1 , nprod/nsave
           ratio=switchingFunction(ntype,ts,k*nsave*delta)
           write(12,'(f15.5,3f15.5,f15.5,2x,f8.3)') 
     *              temperature/ratio,epot(k),ekin(k),
     *              etot(k),converter*ekin(k),para(k)
           if(TI) averU=averU+epot(k)
        enddo
        close(12)
        if(TI) then
           averU=averU/(nprod/nsave)
           write(6,'(a9,2x,f8.4,2x,a6,2x,f15.8)') 
     *             ' lambda= ',ratio,' <U>= ',averU
        endif
        
        if(AST.or.Tsallis) then
           call work(ntot,delta,nprod/nsave,nsave,temperature,
     *               ntype,energyPotential)
           call entropy(ntot,nprod/nsave)
        endif

        close(45)
        
        end