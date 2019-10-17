c this routine finds the bain path for cubic structure
c
        subroutine bain(epsilon,sigma,rcut)
        implicit real*8 (a-h,o-z)
        parameter(nmax=5000)
        common/box/xbox,ybox,zbox
        common/block1/rx(nmax),ry(nmax),rz(nmax),
     :                vx(nmax),vy(nmax),vz(nmax),
     :                ax(nmax),ay(nmax),az(nmax),
     :                bx(nmax),by(nmax),bz(nmax),
     :                cx(nmax),cy(nmax),cz(nmax),
     :                fx(nmax),fy(nmax),fz(nmax)
        logical TsallisStatistics,pure

        TsallisStatistics=.false.
        open(11,file='bainpath')
        
        nxc=5
        nyc=5
        nzc=5
        pure=.true.

        do n = 1 , 20
           cell=1.2d0*sigma+dble(n*5)*0.005d0*sigma
           emin=100000000.d0
           do m = 1 , 100
              zcell=(0.5+m*0.01)*cell
              call lattice(ntot,nxc,nyc,nzc,cell,cell,zcell)
              call force(ntot,epsilon,sigma,rcut,energy,vc,w,pure)
              if(energy.lt.emin) then
                 emin=energy
                 cmin=zcell/cell
              endif
           enddo
           write(6 ,'(2f10.5)') cmin,emin
           write(11,'(2f10.5)') cmin,emin
        enddo

        close(11)

        return
        end
        
c this routine searches the equilibrium lattice constant
c
        subroutine ground(epsilon,sigma,rcut)
        implicit real*8 (a-h,o-z)
        parameter(nmax=5000)
        common/box/xbox,ybox,zbox
        common/block1/rx(nmax),ry(nmax),rz(nmax),
     :                vx(nmax),vy(nmax),vz(nmax),
     :                ax(nmax),ay(nmax),az(nmax),
     :                bx(nmax),by(nmax),bz(nmax),
     :                cx(nmax),cy(nmax),cz(nmax),
     :                fx(nmax),fy(nmax),fz(nmax)
        logical pure

        open(11,file='eofa')
        
        nxc=5
        nyc=5
        nzc=5
        pure=.true.

        do n = 1 , 100
           cell=1.32d0*sigma+dble(n)*0.005d0*sigma
           call lattice(ntot,nxc,nyc,nzc,cell,cell,cell)
           call force(ntot,epsilon,sigma,rcut,energy,vc,w,pure)
           write(6 ,'(2f10.5)') cell,tsallis(1.1d0,1.d0,energy)
           write(11,'(2f10.5)') cell,tsallis(1.1d0,1.d0,energy)
        enddo

        close(11)

        return
        end
        
c Tsallis effective potential
c        
        real*8 function tsallis(q,beta,epot)
        implicit real*8 (a-h,o-z)
        
        tsallis=q/(q-1.d0)/beta*dlog(1.d0-(1.d0-q)*beta*epot+0.5d0)
        
        return
        end