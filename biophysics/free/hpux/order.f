      subroutine order(ntot,cell,param)
      implicit real*8 (a-h,o-z)
      parameter(nmax=5000)
      common/block1/rx(nmax),ry(nmax),rz(nmax),
     :              vx(nmax),vy(nmax),vz(nmax),
     :              ax(nmax),ay(nmax),az(nmax),
     :              bx(nmax),by(nmax),bz(nmax),
     :              cx(nmax),cy(nmax),cz(nmax),
     :              fx(nmax),fy(nmax),fz(nmax)
      common/box/xbox,ybox,zbox

      twopi=8.d0*datan(1.d0)

      xklat = - twopi/cell
      yklat =   twopi/cell
      zklat = - twopi/cell

      sinsum = 0.d0
      cossum = 0.d0

      do i = 1 , ntot
         cossum = cossum + dcos(  xklat * rx(i)
     :                          + yklat * ry(i)
     :                          + zklat * rz(i))
         sinsum = sinsum + dsin(  xklat * rx(i)
     :                          + yklat * ry(i)
     :                          + zklat * rz(i))
      enddo

      cossum = cossum / dble ( ntot )
      sinsum = sinsum / dble ( ntot )
      param  = dsqrt ( cossum*cossum + sinsum*sinsum )
      
      return
      end