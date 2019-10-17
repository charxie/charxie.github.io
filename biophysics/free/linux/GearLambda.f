        subroutine predictorLambda(namb,delta)
        implicit real*8 (a-h,o-z)
        parameter(lmax=20)
        common/blockLambda/rl(lmax),vl(lmax),
     :                     al(lmax),bl(lmax),
     :                     cl(lmax),fl(lmax),
     :                     amassLam(lmax)

        c1 = delta
        c2 = c1*delta/2.d0
        c3 = c2*delta/3.d0
        c4 = c3*delta/4.d0

        do i = 1, namb
           rl(i) = rl(i)+c1*vl(i)+c2*al(i)+c3*bl(i)+c4*cl(i)
           vl(i) = vl(i)+c1*al(i)+c2*bl(i)+c3*cl(i)
           al(i) = al(i)+c1*bl(i)+c2*cl(i)
           bl(i) = bl(i)+c1*cl(i)
        enddo

        return
        end



        subroutine correctorLambda(namb,delta)
c
c amass is the fictitious mass associated with the lambda coordinate
c        
        implicit real*8(a-h,o-z)
        parameter(lmax=20)
        common/blockLambda/rl(lmax),vl(lmax),
     :                     al(lmax),bl(lmax),
     :                     cl(lmax),fl(lmax),
     :                     amassLam(lmax)

        gear0 = 19.d0 / 120.d0 
        gear1 =  3.d0 /   4.d0
        gear3 =  1.d0 /   2.d0    
        gear4 =  1.d0 /  12.d0
        amass=100.d0
        c1 = delta
        c2 = c1 * delta / 2.d0
        c3 = c2 * delta / 3.d0
        c4 = c3 * delta / 4.d0
        cr = gear0 * c2
        cv = gear1 * c2 / c1
        cb = gear3 * c2 / c3
        cc = gear4 * c2 / c4

        do i = 1, namb

           ali = fl(i)/amass
           corr = ali - al(i)

           rl(i) = rl(i) + cr * corr
           vl(i) = vl(i) + cv * corr
           al(i) = ali
           bl(i) = bl(i) + cb * corr
           cl(i) = cl(i) + cc * corr

        enddo

        return
        end
        
        subroutine rescaleLambda(namb)
        implicit real*8(a-h,o-z)
        parameter(lmax=20)
        common/blockLambda/rl(lmax),vl(lmax),
     :                     al(lmax),bl(lmax),
     :                     cl(lmax),fl(lmax),
     :                     amassLam(lmax)

        rsum=0.d0
        vsum=0.d0
        do k = 1 , namb
           rsum=rsum+rl(k)*rl(k)
           vsum=vsum+rl(k)*vl(k)
        enddo
        
        rsum=1.d0/dsqrt(rsum)
        
        do k = 1 , namb
           rl(k)=rl(k)*rsum
           vl(k)=vl(k)*rsum-vsum*rl(k)*rsum*rsum
        enddo

        return
        end