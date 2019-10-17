
      subroutine rk4(tau,n,psi,time,nstp,tmd,h,cof)

c Fourth-order Runge-Kutta method

      implicit real*8 (a-h,o-z)
      complex*8 psi(n),psinew(n)
      real*4 tmd(nstp),h(n,n,nstp),cof(n,n,nstp,3)
      complex*8 ftemp(n),xtemp(n)
      dimension vr(n,n)

      ta2=0.5*tau
      ta3=1.0/3.0
      ta6=1.0/6.0
      dtmd=tmd(2)-tmd(1)
      
      itt0=(time-tau-tmd(1))/dtmd+1
      itt1=(time-ta2-tmd(1))/dtmd+1
      itt2=(time-tmd(1))/dtmd+1
      dtt0=time-tau-tmd(itt0)
      dtt1=time-ta2-tmd(itt1)
      dtt2=time-tmd(itt2)
      
      do i = 1 , n
         psinew(i)=psi(i)
	 ftemp(i)=(0.0,0.0)
	 xtemp(i)=(0.0,0.0)
      enddo

c...f1

      do j = 1 , n
	 do i = 1 , n
	    vr(i,j)=((cof(i,j,itt0,3)*dtt0+cof(i,j,itt0,2))*dtt0
     :            +cof(i,j,itt0,1))*dtt0+h(i,j,itt0)
            ftemp(i)=ftemp(i)-psi(j)*(0.0,1.0)*vr(i,j)
	 enddo
      enddo
      do i = 1 , n
         ftemp(i)=ftemp(i)*tau
         xtemp(i)=psi(i)+0.5*ftemp(i)
         psinew(i)=psinew(i)+ta6*ftemp(i)
         ftemp(i)=(0.0,0.0)
      enddo

c...f2
      
      do j = 1 , n
	 do i = 1 , n
	    vr(i,j)=((cof(i,j,itt1,3)*dtt1+cof(i,j,itt1,2))*dtt1
     :            +cof(i,j,itt1,1))*dtt1+h(i,j,itt1)
            ftemp(i)=ftemp(i)-xtemp(j)*(0.0,1.0)*vr(i,j)
	 enddo
      enddo
      do i = 1 , n
         ftemp(i)=ftemp(i)*tau
         xtemp(i)=psi(i)+0.5*ftemp(i)
         psinew(i)=psinew(i)+ta3*ftemp(i)
         ftemp(i)=(0.0,0.0)
      enddo

c...f3
      
      do j = 1 , n
	 do i = 1 , n
            ftemp(i)=ftemp(i)-xtemp(j)*(0.0,1.0)*vr(i,j)
	 enddo
      enddo
      do i = 1 , n
         ftemp(i)=ftemp(i)*tau
         xtemp(i)=psi(i)+ftemp(i)
         psinew(i)=psinew(i)+ta3*ftemp(i)
         ftemp(i)=(0.0,0.0)
      enddo

c...f4
      
      do j = 1 , n
	 do i = 1 , n
	    vr(i,j)=((cof(i,j,itt2,3)*dtt2+cof(i,j,itt2,2))*dtt2
     :            +cof(i,j,itt2,1))*dtt2+h(i,j,itt2)
            ftemp(i)=ftemp(i)-xtemp(j)*(0.0,1.0)*vr(i,j)
	 enddo
      enddo
      do i = 1 , n
         psi(i)=psinew(i)+ta6*ftemp(i)*tau
      enddo
      
      return
      end
      
      subroutine saverkns(n,time,psi,dotsa)

c save wave function

      implicit real*8 (a-h,o-z)
      complex*8 psi(n)
      logical dotsa
      
      if(.not.dotsa) then
      open(99,file='SAVE.CONT',status='unknown',form='unformatted')
      else
      open(99,file='SAVE2S.CONT',status='unknown',form='unformatted')
      endif

      rewind 99

      write(99) time,n
      do i = 1 , n
	 write(99) psi(i)
      enddo

      close(99)	

      return
      end
      
      subroutine readrkns(n,time0,psi0,dotsa)

c read wave function

      implicit real*8 (a-h,o-z)
      complex*8 psi0(n)
      logical dotsa

      if(.not.dotsa) then
      open(99,file='SAVE.CONT',status='old',form='unformatted')
      else
      open(99,file='SAVE2S.CONT',status='old',form='unformatted')
      endif

      read(99) time0,n0
      
      if(n0.ne.n) then
         write(*,*) n0,n
         stop ' ERROR: SAVED FILE INCONSISTENT '
      endif
      do i = 1 , n
	 read(99) psi0(i)
      enddo
      
      close(99)
      
      return
      end