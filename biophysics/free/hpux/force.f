c calculate the forces, return 
c the potential energy 'v', virial function 'w'
c shifted potential 'vc'
c the logical variable 'pure' = false allows different type of particles
c when 'pure' is set true, the force routine is short-circut
c
        subroutine force(ntot,epsilon,sigma,rcut,v,vc,w,pure)
        implicit real*8 (a-h,o-z)
        parameter(nmax=5000)
        common/fields/eps(nmax),sig(nmax)
        common/block1/rx(nmax),ry(nmax),rz(nmax),
     :                vx(nmax),vy(nmax),vz(nmax),
     :                ax(nmax),ay(nmax),az(nmax),
     :                bx(nmax),by(nmax),bz(nmax),
     :                cx(nmax),cy(nmax),cz(nmax),
     :                fx(nmax),fy(nmax),fz(nmax)
        common/box/xbox,ybox,zbox
        logical pure

c>>>> convert the force into the unit of A/fs**2 to
c     be directly used in the integrator routine
c     * mass unit=12*1.0E-3/6.022*10E+23 Kg
c     * energy unit=1.6*10E-19 J

        converter=1.6d0*6.022d0/12000.d0
        xbox2 = 0.5d0 * xbox
        ybox2 = 0.5d0 * ybox
        zbox2 = 0.5d0 * zbox
        rcutsq = rcut * rcut
        if(pure) then
           sigsq  = sigma * sigma
           eps4   = epsilon * 4.d0
           eps24  = epsilon * 24.d0
        endif

        do i = 1 , ntot
           fx(i) = 0.d0
           fy(i) = 0.d0
           fz(i) = 0.d0
        enddo

        ncut = 0
        v    = 0.d0
        w    = 0.d0

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
                 if(pure) then
                    sr2   = sigsq / rijsq
                 else
                    sigsq = sig(i)*sig(j)
                    sr2   = sigsq / rijsq
                 endif
                 sr6   = sr2 * sr2 * sr2
                 sr12  = sr6 * sr6
                 if(pure) then
                    vij   = sr12 - sr6
                 else
                    vij   =(sr12 - sr6)*2.d0*(eps(i)+eps(j))
                 endif
                 v  = v + vij
                 if(pure) then
                    wij   = vij + sr12
                 else
                    wij   = vij + sr12*2.d0*(eps(i)+eps(j))
                 endif
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

        if(pure) then
           sr2 = sigsq / rcutsq
           sr6 = sr2 * sr2 * sr2
           sr12 = sr6 * sr6
           vij = sr12 - sr6
           vc = v - float(ncut)*vij
        endif

        if(pure) then
           do i = 1 , ntot
              fx(i) = fx(i) * eps24 * converter
              fy(i) = fy(i) * eps24 * converter
              fz(i) = fz(i) * eps24 * converter
           enddo
           v  = v  * eps4
           vc = vc * eps4
           w  = w  * eps24 / 3.d0
        else
           do i = 1 , ntot
              fx(i) = fx(i) * 6.d0 * converter
              fy(i) = fy(i) * 6.d0 * converter
              fz(i) = fz(i) * 6.d0 * converter
           enddo
           w  = w  * 2.d0
        endif
        
        v=v/dble(ntot)
        vc=vc/dble(ntot)
        w=w/dble(ntot)

        return
        end
