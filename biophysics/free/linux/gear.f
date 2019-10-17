        subroutine predictor(ntot,delta)
        implicit real*8 (a-h,o-z)
        parameter(nmax=5000)
        common/block1/rx(nmax),ry(nmax),rz(nmax),
     :                vx(nmax),vy(nmax),vz(nmax),
     :                ax(nmax),ay(nmax),az(nmax),
     :                bx(nmax),by(nmax),bz(nmax),
     :                cx(nmax),cy(nmax),cz(nmax),
     :                fx(nmax),fy(nmax),fz(nmax)

        c1 = delta
        c2 = c1*delta/2.d0
        c3 = c2*delta/3.d0
        c4 = c3*delta/4.d0

        do i = 1, ntot
           rx(i) = rx(i)+c1*vx(i)+c2*ax(i)+c3*bx(i)+c4*cx(i)
           ry(i) = ry(i)+c1*vy(i)+c2*ay(i)+c3*by(i)+c4*cy(i)
           rz(i) = rz(i)+c1*vz(i)+c2*az(i)+c3*bz(i)+c4*cz(i)
           vx(i) = vx(i)+c1*ax(i)+c2*bx(i)+c3*cx(i)
           vy(i) = vy(i)+c1*ay(i)+c2*by(i)+c3*cy(i)
           vz(i) = vz(i)+c1*az(i)+c2*bz(i)+c3*cz(i)
           ax(i) = ax(i)+c1*bx(i)+c2*cx(i)
           ay(i) = ay(i)+c1*by(i)+c2*cy(i)
           az(i) = az(i)+c1*bz(i)+c2*cz(i)
           bx(i) = bx(i)+c1*cx(i)
           by(i) = by(i)+c1*cy(i)
           bz(i) = bz(i)+c1*cz(i)
        enddo

        return
        end



        subroutine corrector(DynLambda,namb,nsolute,
     :                       ntot,delta,energyKinetic)
        implicit real*8(a-h,o-z)
        parameter(nmax=5000)
        parameter(lmax=20)
        common/block1/rx(nmax),ry(nmax),rz(nmax),
     :                vx(nmax),vy(nmax),vz(nmax),
     :                ax(nmax),ay(nmax),az(nmax),
     :                bx(nmax),by(nmax),bz(nmax),
     :                cx(nmax),cy(nmax),cz(nmax),
     :                fx(nmax),fy(nmax),fz(nmax)
        common/blockLambda/rl(lmax),vl(lmax),
     :                     al(lmax),bl(lmax),
     :                     cl(lmax),fl(lmax),
     :                     amassLam(lmax)
        logical DynLambda

        amass=1.d0
        if(DynLambda) then
           amassHybrid=0.d0
           do k = 1 , namb
              amassHybrid=amassHybrid+rl(k)*rl(k)*amassLam(k)
           enddo
        endif
        
        gear0 = 19.d0 / 120.d0 
        gear1 =  3.d0 /   4.d0
        gear3 =  1.d0 /   2.d0    
        gear4 =  1.d0 /  12.d0
        c1 = delta
        c2 = c1 * delta / 2.d0
        c3 = c2 * delta / 3.d0
        c4 = c3 * delta / 4.d0
        cr = gear0 * c2
        cv = gear1 * c2 / c1
        cb = gear3 * c2 / c3
        cc = gear4 * c2 / c4

        energyKinetic = 0.d0

        do i = 1, ntot

           axi = fx(i)/amass
           ayi = fy(i)/amass
           azi = fz(i)/amass
           if(DynLambda.and.i.eq.nsolute) then 
              axi = fx(i)/amassHybrid
              ayi = fy(i)/amassHybrid
              azi = fz(i)/amassHybrid
           endif
           corrx = axi - ax(i)
           corry = ayi - ay(i)
           corrz = azi - az(i)

           rx(i) = rx(i) + cr * corrx
           ry(i) = ry(i) + cr * corry
           rz(i) = rz(i) + cr * corrz
           vx(i) = vx(i) + cv * corrx
           vy(i) = vy(i) + cv * corry
           vz(i) = vz(i) + cv * corrz
           ax(i) = axi
           ay(i) = ayi
           az(i) = azi
           bx(i) = bx(i) + cb * corrx
           by(i) = by(i) + cb * corry
           bz(i) = bz(i) + cb * corrz
           cx(i) = cx(i) + cc * corrx
           cy(i) = cy(i) + cc * corry
           cz(i) = cz(i) + cc * corrz

           energyKinetic = energyKinetic + amass * 
     *                  (vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))
           if(DynLambda.and.i.eq.nsolute) then
              energyKinetic = energyKinetic + (amassHybrid-amass)*
     *                  (vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))
           endif

        enddo

        energyKinetic = 0.5*energyKinetic/dble(ntot)
        converter = 12000.d0/1.6d0/6.022d0
        energyKinetic = energyKinetic*converter

c>>> kinetic energy is converted to eV

        return
        end


        subroutine setPBC(ntot)
        implicit real*8(a-h,o-z)
        parameter(nmax=5000)
        common/block1/rx(nmax),ry(nmax),rz(nmax),
     :                vx(nmax),vy(nmax),vz(nmax),
     :                ax(nmax),ay(nmax),az(nmax),
     :                bx(nmax),by(nmax),bz(nmax),
     :                cx(nmax),cy(nmax),cz(nmax),
     :                fx(nmax),fy(nmax),fz(nmax)
        common/box/xbox,ybox,zbox

c        write(6,*) xbox,ybox,zbox

        do i = 1 , ntot

           if(rx(i).lt.0.d0) rx(i)=rx(i)+xbox
           if(rx(i).ge.xbox) rx(i)=rx(i)-xbox
           if(ry(i).lt.0.d0) ry(i)=ry(i)+ybox
           if(ry(i).ge.ybox) ry(i)=ry(i)-ybox
           if(rz(i).lt.0.d0) rz(i)=rz(i)+zbox
           if(rz(i).ge.zbox) rz(i)=rz(i)-zbox

        enddo
         
        return
        end