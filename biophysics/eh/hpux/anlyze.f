      subroutine anlyze(t,ndim)

c***********************************************************************
c print out details of T matrix
c***********************************************************************

      implicit real*8 (a-h,o-z)
      parameter(natmax=200,ndmax=500)
      integer ul
      logical write
      character nel*2,ac*4,symbol*2,bname*2
      dimension hyb(4),coef(6),pow(6),theta(6),phi(6),nel(6)
      dimension t(ndim,ndim)
      common/ato/ac(natmax),symbol(40)
      common/lbl/bname(ndmax,2),label(ndmax,3),ibx(ndmax)
      common/lblaro/labar(ndmax,7)
      common/info1/ul(natmax),ll(natmax),nele,nocc
      common/iounit/iin,iout,istart(40),iungf,iunhl
      
      write=.false.
      
      if(write) then
         write(iout,*)
         write(iout,*) 'AO-BO transformation matrix'
         write(iout,*)
         write(iout,*) ndim
         do i = 1 , ndim
            write(iout,'(i8,500f8.4)') i,(t(j,i),j=1,ndim)
         enddo
      endif
      write(iout,1400)
      
      do nbond = 1 , ndim
         do ib = 1 , ndim
            if(ibx(ib).eq.nbond) goto 90
         enddo
   90    if(bname(ib,1).eq.'LP') nctr=1
         if(bname(ib,1).eq.'BD') nctr=2
         if(bname(ib,1).eq.'PB') nctr=6
         do n = 1 , nctr
            if(nctr.eq.6) then
               i=labar(ib,n+1)
            else
               i=label(ib,n+1)
            endif
            nel(n)=ac(i)(1:2)
            kl=ll(i)
            ku=ul(i)
            do k = 1 , 4
               hyb(k)=0.0
            enddo
            kh=0
            do k = kl , ku
               kh=kh+1
               hyb(kh)=t(k,nbond)
            enddo
            call htype(hyb,coef(n),pow(n),theta(n),phi(n))
            if(pow(n).eq.100.0) pow(n)=99.99
         enddo
         if(nctr.eq.6) then
            write(iout,1000)
     :        nbond,(bname(ib,j),j=1,2),(nel(n),labar(ib,n+1),
     :        coef(n),pow(n),theta(n),phi(n),n=1,nctr)
         else
            write(iout,1000)
     :        nbond,(bname(ib,j),j=1,2),(nel(n),label(ib,n+1),
     :        coef(n),pow(n),theta(n),phi(n),n=1,nctr)
         endif

      enddo
      write(iout,1500)

 1000 format(1x,i3,'. ',a2,a1,2x,6(a2,i3,' (',f7.4,'*SP',f5.2,';(',             
     + f5.1,',',f6.1,'))',3x))                                                  
 1400 format(6x,120('-'))                                                       
 1500 format(/)                                                                 

      return
      end


      subroutine gett(t,ndim,iwanted)
      
      implicit real*8 (a-h,o-z)
      parameter(ndmax=500,natmax=200)
      integer ul
      character bname*2
      dimension hyb(4),t(ndim,ndim)
      common/lbl/bname(ndmax,2),label(ndmax,3),ibx(ndmax)
      common/info1/ul(natmax),ll(natmax),nele,nocc
      common/iounit/iin,iout,istart(40),iungf,iunhl

      open(15,file='tij',access='append') 
      
      do ib = 1 , ndim
	 if(ibx(ib).eq.iwanted) iwant=ib
      enddo
      if(bname(iwant,1).eq.'LP') nctr=1
      if(bname(iwant,2).eq.'BD') nctr=2
      do ia = 1 , nctr
         kl=ll(label(iwant,ia+1))
         ku=ul(label(iwant,ia+1))
         do k = 1 , 4
            hyb(k)=0.0
         enddo
         kh=0
         do k = kl , ku
            kh=kh+1
            hyb(kh)=t(k,iwanted)
         enddo
      write(15,'(2i5,10f10.5)') iwanted,ia,(hyb(ik),ik=1,kh)
      enddo
      close(15)

      return
      end

c return coefficients, p-characters, and directions of hybrids
c (h(i),i=1,4) from transformation matrix
c pow = 0,1,2,... for s, sp, sp2, ... hybrid (pow=100 for pure p orbital)
c coef = coefficient of hybrid
c theta, phi = polar and azimuthal angles of directed hybrid.

      subroutine htype(h,coef,pow,theta,phi)
      
      implicit real*8 (a-h,o-z)
      dimension h(4)
      data eps/1.d-4/
      data pi/3.14159265357979d0/
      
      sq=0.0
      do k = 1 , 4
         sq=sq+h(k)*h(k)
      enddo
      coef=dsqrt(sq)
      if(h(1).lt.0.0d0) coef=-coef
      diff=dabs(coef-h(1))
      if(dabs(coef).gt.eps.and.diff.gt.eps) goto 20
c zero orbital or pure s orbital      
      pow=0.0
      theta=0.0
      phi=0.0
      return
c start with normalized orbital, positive s coefficient      
   20 do k = 1 , 4
         h(k)=h(k)/coef
      enddo
      if(dabs(h(1)).gt.eps) goto 40
c pure p orbital
      pow=100.0
      goto 60
   40 pow=1.0/h(1)/h(1)-1.0
      if(pow.le.0.0) pow=0.0
      if(pow.ge.99.0) pow=99.0
      rad=dsqrt(pow)
c find coefficients of normalized p orbital
      do k = 2 , 4
         h(k)=h(k)/rad/h(1)
      enddo
   60 continue
      theta=-1.0
      sin2th=1.0-h(4)*h(4)
      if(dabs(h(4)).gt.eps) goto 70
      theta=pi/2.0   
   70 if(sin2th.le.0.0) sin2th=0.0
      sinth=dsqrt(sin2th)
      if(dabs(sinth).gt.eps) goto 80
      theta=0.0
      if(h(4).lt.0.0) theta=pi
      phi=0.0
      goto 110
   80 cosph=h(2)/sinth
      sinph=h(3)/sinth
      if(dabs(cosph).gt.eps) goto 90
      phi=pi/2.0
      if(sinph.lt.0.0) phi=-phi
      goto 100
   90 phi=datan2(sinph,cosph)
  100 if(theta.lt.0.0) theta=datan2(sinth,h(4))
  110 phi=phi*360.d0/(2.d0*pi)
      theta=theta*360.d0/(2.d0*pi)
      
      return
      end