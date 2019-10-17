
        implicit real*8 (a-h,o-z)
        parameter(nprod=10000)
        dimension free(nprod),temp(nprod),entr(nprod)
        
        open(11,file='entropy')
        open(12,file='free')
        
        do n = 1 , nprod
           read(12,'(f10.2,f20.10)') temp(n),free(n)
        enddo
        
        call slope(temp,free,entr,nprod)
        
        do n = 1 , nprod-1
           if(temp(n).lt.48000.d0)
     *        write(11,'(f10.2,f20.10)') temp(n),-entr(n)
        enddo
        
        close(11)
        close(12)

        end
        
      subroutine slope(x,y,d1y,n)
c
c--> calculates the slope of a curve
c
      implicit real*8 (a-h,o-z)
      dimension x(n),y(n)
      dimension d1y(n)
      do i=1,n-1
         if(x(i+1).ne.x(i)) then
            d1y(i)=(y(i+1)-y(i))/(x(i+1)-x(i))
         else
            d1y(i)=0.d0
         endif
      enddo
      return
      end

        