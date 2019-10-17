
      subroutine rk4(tau,n,hami,psi)

c Fourth-order Runge-Kutta method

      implicit real*8 (a-h,o-z)
      complex psi(n),psinew(n)
      dimension hami(n,n)
      complex ftemp(n),xtemp(n)

      ta2=0.5*tau
      ta3=1.0/3.0
      ta6=1.0/6.0
      
      do i = 1 , n
         psinew(i)=psi(i)
	 ftemp(i)=(0.0,0.0)
	 xtemp(i)=(0.0,0.0)
      enddo

c...f1

      do j = 1 , n
	 do i = 1 , n
            ftemp(i)=ftemp(i)-psi(j)*(0.0,1.0)*hami(i,j)
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
            ftemp(i)=ftemp(i)-xtemp(j)*(0.0,1.0)*hami(i,j)
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
            ftemp(i)=ftemp(i)-xtemp(j)*(0.0,1.0)*hami(i,j)
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
            ftemp(i)=ftemp(i)-xtemp(j)*(0.0,1.0)*hami(i,j)
	 enddo
      enddo
      do i = 1 , n
         psi(i)=psinew(i)+ta6*ftemp(i)*tau
      enddo
      
      return
      end
