      implicit real*8(a-h,o-z)
      parameter(nmax=35)
      dimension x(nmax),y(nmax),z(nmax),sum(nmax),free(nmax)
      dimension temperature(nmax)
      dimension temp(500),fe(500)
      
      open(21,file='free.ana')
      
      u0=-8.51d0
      ak1=1.158*0.0001d0*3.0/3.0
      ak2=0.997d-4
      akb=1.38d0/1.6d0*0.0001d0
      do i = 1 , 500
         temp(i)=(i-1)*100.d0+300.d0
         if(temp(i).lt.27000.d0) then
            fe(i)=u0-(ak1+1.5d0*akb)*temp(i)*dlog(temp(i)/300.d0)
     *             +ak1*(temp(i)-300)
            write(21,*) temp(i),fe(i)
         else
            
         endif
      enddo
      
      close(21)
      
      open(11,file='ti01.dat')
      do i = 1 , nmax
         read(11,*) x(i),y(i),z(i)
         temperature(i)=300.d0/x(i)
      enddo
      close(11)
      
      do n = 1 , nmax-1
         sum(n)=0.d0
c         do i = 1 , n
c            sum(n)=sum(n)+(x(i+1)-x(i))*(y(i+1)+y(i))*0.5d0
c         enddo
         do i = 1 , n
            sum(n)=sum(n)-(temperature(i+1)-temperature(i))
     *            *(y(i+1)/temperature(i+1)**2
     *            +y(i)/temperature(i)**2)*0.5d0
         enddo
      enddo
      
      do n = 2 , nmax
c         free(n)=1.5d0*1.38d0*300.d0/x(n)*dlog(x(n))/1.6d0*0.0001d0
c     *          +(-8.43+sum(n-1))/x(n)
c         write(6,*) 300.d0/x(n),free(n)
         u0=-8.51d0
         free(n)=1.5d0*1.38d0*temperature(n)*dlog(x(n))/1.6d0*0.0001d0
     *         +u0/x(n)
     *         +temperature(n)*sum(n-1)
         write(6,*) temperature(n),free(n)
      enddo
      
      open(14,file='internal.ti')
      do n = 2, nmax
         write(14,*) 300/x(n),y(n)
      enddo
      close(14)
      
      open(13,file='entropy.ti')
      do n = 3, nmax
         write(13,*) 300/x(n),-(free(n)-free(n-1))/(300/x(n)-300/x(n-1))
      enddo
      close(13)
      
      open(12,file='order.par')
      do n = 1 , nmax
         write(12,'(10f20.5)') 300.d0/x(n),z(n)
      enddo
      close(12)
      
      
      end