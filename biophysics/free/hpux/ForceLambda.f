c the force routine for the lambda-dynamics
c
        subroutine forceLambda(ntot,nsolute,namb,rcut,v,vc,w)
        implicit real*8 (a-h,o-z)
        parameter(nmax=5000,lmax=20)
        common/paramLambda/epslam(lmax),siglam(lmax)
        common/blockLambda/rl(lmax),vl(lmax),
     :                     al(lmax),bl(lmax),
     :                     cl(lmax),fl(lmax),
     :                     amassLam(lmax),amassGen
        common/fields/eps(nmax),sig(nmax)
        common/block1/rx(nmax),ry(nmax),rz(nmax),
     :                vx(nmax),vy(nmax),vz(nmax),
     :                ax(nmax),ay(nmax),az(nmax),
     :                bx(nmax),by(nmax),bz(nmax),
     :                cx(nmax),cy(nmax),cz(nmax),
     :                fx(nmax),fy(nmax),fz(nmax)
        common/box/xbox,ybox,zbox
        dimension vadd(namb)

c>>>> convert the force into the unit of A/fs**2 to
c     be directly used in the integrator routine
c     * mass unit=12*1.0E-3/6.022*10E+23 Kg
c     * energy unit=1.6*10E-19 J

        converter=1.6d0*6.022d0/12000.d0
        xbox2 = 0.5d0 * xbox
        ybox2 = 0.5d0 * ybox
        zbox2 = 0.5d0 * zbox
        rcutsq = rcut * rcut

        do i = 1 , ntot
           fx(i) = 0.d0
           fy(i) = 0.d0
           fz(i) = 0.d0
        enddo

        ncut = 0
        v    = 0.d0
        w    = 0.d0

c>>> these two do loops calculate the interatomic forces as usual

        do i = 1, ntot-1

           rxi = rx(i)
           ryi = ry(i)
           rzi = rz(i)
           fxi = fx(i)
           fyi = fy(i)
           fzi = fz(i)

           do j = i + 1, ntot

              rxij = rxi - rx(j)
              ryij = ryi - ry(j)
              rzij = rzi - rz(j)
              if(rxij.gt. xbox2) rxij=rxij-xbox
              if(rxij.le.-xbox2) rxij=rxij+xbox
              if(ryij.gt. ybox2) ryij=ryij-ybox
              if(ryij.le.-ybox2) ryij=ryij+ybox
              if(rzij.gt. zbox2) rzij=rzij-zbox
              if(rzij.le.-zbox2) rzij=rzij+zbox
              rijsq = rxij*rxij + ryij*ryij + rzij*rzij

              if ( rijsq .lt. rcutsq ) then
                 sigsq = sig(i)*sig(j)
                 sr2   = sigsq / rijsq
                 sr6   = sr2 * sr2 * sr2
                 sr12  = sr6 * sr6
                 vij   = (sr12 - sr6)*2.d0*(eps(i)+eps(j))
                 v     = v + vij
                 wij   = vij + sr12*2.d0*(eps(i)+eps(j))
                 w     = w + wij
                 fij   = wij / rijsq
                 fxij  = fij * rxij
                 fyij  = fij * ryij
                 fzij  = fij * rzij
                 fxi   = fxi + fxij
                 fyi   = fyi + fyij
                 fzi   = fzi + fzij
                 fx(j) = fx(j) - fxij
                 fy(j) = fy(j) - fyij
                 fz(j) = fz(j) - fzij
                 ncut  = ncut + 1
              endif

           enddo

           fx(i) = fxi
           fy(i) = fyi
           fz(i) = fzi

        enddo

c        write(6,*) v, (epslam(i),i=1,namb)
        do k = 1 , namb
           vadd(k)=0.d0
        enddo

c>>> in what follows we calculate the interatomic forces with the hybrid
c>>> particle, indexed by 'nsolute'

        do i = 1 , ntot
        
           if(i.ne.nsolute) then
           
              rxij = rx(i) - rx(nsolute)
              ryij = ry(i) - ry(nsolute)
              rzij = rz(i) - rz(nsolute)
              if(rxij.gt. xbox2) rxij=rxij-xbox
              if(rxij.le.-xbox2) rxij=rxij+xbox
              if(ryij.gt. ybox2) ryij=ryij-ybox
              if(ryij.le.-ybox2) ryij=ryij+ybox
              if(rzij.gt. zbox2) rzij=rzij-zbox
              if(rzij.le.-zbox2) rzij=rzij+zbox
              rijsq = rxij*rxij + ryij*ryij + rzij*rzij

              if ( rijsq .lt. rcutsq ) then

                 do k = 1 , namb
                    
                    sigsq = sig(i)*siglam(k)
                    sr2   = sigsq / rijsq
                    sr6   = sr2 * sr2 * sr2
                    sr12  = sr6 * sr6
                    tigsq = sig(i)*sig(i)
                    tr2   = tigsq / rijsq
                    tr6   = tr2 * tr2 * tr2
                    tr12  = tr6 * tr6
                    
                    vij   = ((sr12-sr6)*2.d0*(eps(i)+epslam(k))
     :                      -(tr12-tr6)*4.d0*eps(i))*rl(k)*rl(k)
                    v     = v + vij
                    vadd(k) = vadd(k) + vij/rl(k)/rl(k)
                    wij   = vij + (sr12*2.d0*(eps(i)+epslam(k))
     :                            -tr12*4.d0*eps(i))*rl(k)*rl(k)
                    w     = w + wij
                    
                    fij   = wij / rijsq
                    fxij  = fij * rxij
                    fyij  = fij * ryij
                    fzij  = fij * rzij
                    fx(i) = fx(i) + fxij
                    fy(i) = fy(i) + fyij
                    fz(i) = fz(i) + fzij
                    fx(nsolute)=fx(nsolute)-fxij
                    fy(nsolute)=fy(nsolute)-fyij
                    fz(nsolute)=fz(nsolute)-fzij

                 enddo

              endif

           endif
                      
        enddo
        
c        write(6,*) v
c        vtemp=0.d0
c        do k = 1 , namb
c           vtemp=vtemp+vadd(k)*rl(k)**2
c        enddo
c        write(6,*) (vadd(k),k=1,namb),vtemp

        do i = 1 , ntot
           fx(i) = fx(i) * 6.d0 * converter
           fy(i) = fy(i) * 6.d0 * converter
           fz(i) = fz(i) * 6.d0 * converter
        enddo
        w  = w  * 2.d0
        
        v=v/dble(ntot)
        vc=vc/dble(ntot)
        w=w/dble(ntot)
        
c>>> now evaluate the force for the lambda's
c update the lagrange multiplier first
c 'amassGen' is the mass for the lambda variables

        converter2 = 12000.d0/1.6d0/6.022d0
        sum1=0.d0
        sum2=0.d0
        do k = 1 , namb
          sum2=sum2+rl(k)*rl(k)
          sum1=sum1+2.d0*vadd(k)*rl(k)*rl(k)
     :        -amassGen*vl(k)*vl(k)
     :        -amassLam(k)*rl(k)*rl(k)*converter2*(
     :           vx(nsolute)*vx(nsolute)+
     :           vy(nsolute)*vy(nsolute)+
     :           vz(nsolute)*vz(nsolute))
        enddo
        gamma=sum1/sum2*2.d0

        do k = 1 , namb
           fl(k)=-2.d0*rl(k)*vadd(k)
           fl(k)=fl(k)+amassLam(k)*(
     :           vx(nsolute)*vx(nsolute)+
     :           vy(nsolute)*vy(nsolute)+
     :           vz(nsolute)*vz(nsolute))*rl(k)*converter2+
     :           gamma*rl(k)
        enddo
        

        return
        end