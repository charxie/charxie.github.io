c lambda initializer: let's put all the competitors on the same starting
c line

        subroutine initLambda(namb,epsilon,sigma)
        implicit real*8(a-h,o-z)
        parameter(lmax=20)
        common/paramLambda/epslam(lmax),siglam(lmax)
        common/blockLambda/rl(lmax),vl(lmax),
     :                     al(lmax),bl(lmax),
     :                     cl(lmax),fl(lmax),
     :                     amassLam(lmax)

        do i = 1 , namb
           siglam(i)=sigma
           epslam(i)=epsilon
           rl(i)=1.d0/dsqrt(dble(namb))
           vl(i)=0.d0
           al(i)=0.d0
           bl(i)=0.d0
           cl(i)=0.d0
           fl(i)=0.d0
           amassLam(i)=1.d0+(i-1)*1.d0
        enddo
     
c        write(6,*) (epslam(i),i=1,namb)
     
        return
        end


c create a mixture

        subroutine solute(ntot,nsolute,epsilon,sigma)
        implicit real*8 (a-h,o-z)
        parameter(nmax=5000)
        common/fields/eps(nmax),sig(nmax)
        common/box/xbox,ybox,zbox
        common/block1/rx(nmax),ry(nmax),rz(nmax),
     :                vx(nmax),vy(nmax),vz(nmax),
     :                ax(nmax),ay(nmax),az(nmax),
     :                bx(nmax),by(nmax),bz(nmax),
     :                cx(nmax),cy(nmax),cz(nmax),
     :                fx(nmax),fy(nmax),fz(nmax)

        do i = 1 , ntot
           if(i.eq.nsolute) then
              eps(i)=0.0
              sig(i)=sigma
           else
              eps(i)=epsilon
              sig(i)=sigma
           endif      
        enddo
        
        write(6,*) ' Lennard-Jones mixture created!'

        return
        end

c latticeBuilder
c cell -- lattice constant;
c

        subroutine lattice(ntot,nxc,nyc,nzc,xcell,ycell,zcell)
        implicit real*8 (a-h,o-z)
        parameter(nmax=5000)
        common/box/xbox,ybox,zbox
        common/block1/rx(nmax),ry(nmax),rz(nmax),
     :                vx(nmax),vy(nmax),vz(nmax),
     :                ax(nmax),ay(nmax),az(nmax),
     :                bx(nmax),by(nmax),bz(nmax),
     :                cx(nmax),cy(nmax),cz(nmax),
     :                fx(nmax),fy(nmax),fz(nmax)
        character type*3

        type='fcc'
        xbox=nxc*xcell
        ybox=nyc*ycell
        zbox=nzc*zcell
        
        if(type.eq.'fcc') then      
           ntot = 4*nxc*nyc*nzc
           if(ntot.gt.nmax) stop ' too many atoms '
           do i = 1, ntot
              rx(i) = 0.0
              ry(i) = 0.0
              rz(i) = 0.0           
           enddo
c--> first build the unit cell seed and then grow the entire lattice
           rx(1) =  0.0
           ry(1) =  0.0
           rz(1) =  0.0
           rx(2) =  0.5
           ry(2) =  0.5
           rz(2) =  0.0      
           rx(3) =  0.0
           ry(3) =  0.5
           rz(3) =  0.5
           rx(4) =  0.5
           ry(4) =  0.0
           rz(4) =  0.5
           m = 0
           do iz = 1, nzc
              do iy = 1, nyc
                 do ix = 1, nxc
                    do iref = 1, 4
                       rx(iref+m) = rx(iref) + dble(ix-1)
                       ry(iref+m) = ry(iref) + dble(iy-1)
                       rz(iref+m) = rz(iref) + dble(iz-1)
                    enddo
                    m = m + 4
                 enddo
              enddo
           enddo
        endif
        
        do i = 1 , ntot
           rx(i)=rx(i)*xcell
           ry(i)=ry(i)*ycell
           rz(i)=rz(i)*zcell
        enddo
        
        return
        end