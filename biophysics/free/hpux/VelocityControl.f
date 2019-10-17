c unit: length in angstrom, time in femtosecond, 
c energy in eletronic volt, mass taken as 12 gram per mole;
c
        subroutine assignVel(ntot,temperature,energyKinetic)
        implicit real*8 (a-h,o-z)
        parameter(nmax=5000)
        parameter(amass=1.d0)
        common/box/xbox,ybox,zbox
        common/block1/rx(nmax),ry(nmax),rz(nmax),
     :                vx(nmax),vy(nmax),vz(nmax),
     :                ax(nmax),ay(nmax),az(nmax),
     :                bx(nmax),by(nmax),bz(nmax),
     :                cx(nmax),cy(nmax),cz(nmax),
     :                fx(nmax),fy(nmax),fz(nmax)
        real*8 gauss

        rtemp = dsqrt(temperature)*0.00026316
        sumvx = 0.0
        sumvy = 0.0
        sumvz = 0.0
        do i = 1 , ntot
           vx(i) = rtemp*gauss()
           vy(i) = rtemp*gauss()
           vz(i) = rtemp*gauss()
           sumvx = sumvx + vx(i)
           sumvy = sumvy + vy(i)
           sumvz = sumvz + vz(i)
        enddo
        sumvx = sumvx / ntot
        sumvy = sumvy / ntot
        sumvz = sumvz / ntot
        do i = 1 , ntot
           vx(i)=vx(i)-sumvx
           vy(i)=vy(i)-sumvy
           vz(i)=vz(i)-sumvz
        enddo
        do i = 1 , ntot
           energyKinetic = energyKinetic + 
     *                     vx(i)*vx(i) + vy(i)*vy(i) + vz(i)*vz(i)
        enddo
        energyKinetic = 0.5*amass*energyKinetic/ntot
        energyKinetic = energyKinetic*12000.d0/1.6d0/6.022d0
        
        return
        end
 
c Random variable from the standard normal distribution; 
c The distribution is a Gaussian function with zero mean and unit 
c variance; 

$HP9000_800 INTRINSICS ON

        real*8 function gauss()
        implicit real*8(a-h,o-z)
        parameter ( A1 = 3.949846138, A3 = 0.252408784 )
        parameter ( A5 = 0.076542912, A7 = 0.008355968 )
        parameter ( A9 = 0.029899776                   )
        real rand

        sum = 0.0
        do i = 1, 12
           sum = sum + rand()
        enddo
        r  = (sum-6.0) / 4.0
        r2 = r * r
        gauss = ((((A9*R2+A7)*R2+A5)*R2+A3)*R2+A1)*R

        return
        end

$HP9000_800 INTRINSICS OFF
        
        subroutine rescaleVel(ntot,temperature,tolerance,ekinet,
     *                        rescaled)
        implicit real*8 (a-h,o-z)
        parameter(nmax=5000)
        parameter(amass=1.0)
        common/box/xbox,ybox,zbox
        common/block1/rx(nmax),ry(nmax),rz(nmax),
     :                vx(nmax),vy(nmax),vz(nmax),
     :                ax(nmax),ay(nmax),az(nmax),
     :                bx(nmax),by(nmax),bz(nmax),
     :                cx(nmax),cy(nmax),cz(nmax),
     :                fx(nmax),fy(nmax),fz(nmax)
        character*1 rescaled
        
        converter=2.d0*1.6d0/3.d0/1.38d0*10000.d0
        getTemp=converter*ekinet
        
        if(dabs(getTemp-temperature).gt.tolerance) then
          factor=dsqrt(temperature/getTemp)
          do i = 1 , ntot
             vx(i)=vx(i)*factor
             vy(i)=vy(i)*factor
             vz(i)=vz(i)*factor
          enddo
          rescaled='*'
        endif
        
        return
        end