c routines for matrix multiplication: strassen and gaxpy

       subroutine strassen(ndim,nhaf,a,b,c)
       implicit real*8 (a-h,o-z)
       dimension a(ndim,ndim),b(ndim,ndim),c(ndim,ndim)
       dimension a11(nhaf,nhaf),a12(nhaf,nhaf)
       dimension a21(nhaf,nhaf),a22(nhaf,nhaf)
       dimension b11(nhaf,nhaf),b12(nhaf,nhaf)
       dimension b21(nhaf,nhaf),b22(nhaf,nhaf)
       dimension p(nhaf,nhaf)
       dimension ta(nhaf,nhaf),tb(nhaf,nhaf)
       
       do j = 1 , ndim
          do i = 1 , ndim
             c(i,j)=0.0
          enddo
       enddo
       
       if(mod(ndim,2).eq.0) then
          if(nhaf.ne.ndim/2) stop 'nhaf!=ndim/2'
       else
          if(nhaf.ne.(ndim+1)/2) stop 'nhaf!=(ndim+1)/2'
       endif
       
       do j = 1 , nhaf
          do i = 1 , nhaf
             if(nhaf.eq.ndim/2) then
                a11(i,j)=a(i,j)
                a12(i,j)=a(i,j+nhaf)
                a21(i,j)=a(i+nhaf,j)
                a22(i,j)=a(i+nhaf,j+nhaf)
                b11(i,j)=b(i,j)
                b12(i,j)=b(i,j+nhaf)
                b21(i,j)=b(i+nhaf,j)
                b22(i,j)=b(i+nhaf,j+nhaf)
             else
                a11(i,j)=a(i,j)
                b11(i,j)=b(i,j)
                if(j+nhaf.le.ndim) then
                   a12(i,j)=a(i,j+nhaf)
                   b12(i,j)=b(i,j+nhaf)
                else
                   a12(i,j)=0.0
                   b12(i,j)=0.0
                endif
                if(i+nhaf.le.ndim) then
                   a21(i,j)=a(i+nhaf,j)
                   b21(i,j)=b(i+nhaf,j)
                else
                   a21(i,j)=0.0
                   b21(i,j)=0.0
                endif
                if(i+nhaf.le.ndim.and.j+nhaf.le.ndim) then
                   a22(i,j)=a(i+nhaf,j+nhaf)
                   b22(i,j)=b(i+nhaf,j+nhaf)
                else
                   a22(i,j)=0.0
                   b22(i,j)=0.0
                endif
             endif
          enddo
       enddo
c...p1
       do j = 1 , nhaf
          do i = 1 , nhaf
             ta(i,j)=a11(i,j)+a22(i,j)
             tb(i,j)=b11(i,j)+b22(i,j)
             p(i,j)=0.0
          enddo
       enddo
       do j = 1 , nhaf
          do k = 1 , nhaf
             do i = 1 , nhaf
                p(i,j)=p(i,j)+ta(i,k)*tb(k,j)
             enddo
          enddo
       enddo
       do j = 1 , nhaf
          do i = 1 , nhaf
             c(i,j)=c(i,j)+p(i,j)
             if(i+nhaf.le.ndim.and.j+nhaf.le.ndim) then
                c(i+nhaf,j+nhaf)=c(i+nhaf,j+nhaf)+p(i,j)
             endif
          enddo
       enddo
c...p2
       do j = 1 , nhaf
          do i = 1 , nhaf
             ta(i,j)=a21(i,j)+a22(i,j)
             tb(i,j)=b11(i,j)
             p(i,j)=0.0
          enddo
       enddo
       do j = 1 , nhaf
          do k = 1 , nhaf
             do i = 1 , nhaf
                p(i,j)=p(i,j)+ta(i,k)*tb(k,j)
             enddo
          enddo
       enddo
       do j = 1 , nhaf
          do i = 1 , nhaf
             if(i+nhaf.le.ndim) then
                c(i+nhaf,j)=c(i+nhaf,j)+p(i,j)
             endif
             if(i+nhaf.le.ndim.and.j+nhaf.le.ndim) then
                c(i+nhaf,j+nhaf)=c(i+nhaf,j+nhaf)-p(i,j)
             endif
          enddo
       enddo
c...p3
       do j = 1 , nhaf
          do i = 1 , nhaf
             ta(i,j)=a11(i,j)
             tb(i,j)=b12(i,j)-b22(i,j)
             p(i,j)=0.0
          enddo
       enddo
       do j = 1 , nhaf
          do k = 1 , nhaf
             do i = 1 , nhaf
                p(i,j)=p(i,j)+ta(i,k)*tb(k,j)
             enddo
          enddo
       enddo
       do j = 1 , nhaf
          do i = 1 , nhaf
             if(j+nhaf.le.ndim) then
                c(i,j+nhaf)=c(i,j+nhaf)+p(i,j)
             endif
             if(i+nhaf.le.ndim.and.j+nhaf.le.ndim) then
                c(i+nhaf,j+nhaf)=c(i+nhaf,j+nhaf)+p(i,j)
             endif
          enddo
       enddo
c...p4
       do j = 1 , nhaf
          do i = 1 , nhaf
             ta(i,j)=a22(i,j)
             tb(i,j)=b21(i,j)-b11(i,j)
             p(i,j)=0.0
          enddo
       enddo
       do j = 1 , nhaf
          do k = 1 , nhaf
             do i = 1 , nhaf
                p(i,j)=p(i,j)+ta(i,k)*tb(k,j)
             enddo
          enddo
       enddo
       do j = 1 , nhaf
          do i = 1 , nhaf
             c(i,j)=c(i,j)+p(i,j)
             if(i+nhaf.le.ndim) then
                c(i+nhaf,j)=c(i+nhaf,j)+p(i,j)
             endif
          enddo
       enddo
c...p5
       do j = 1 , nhaf
          do i = 1 , nhaf
             ta(i,j)=a11(i,j)+a12(i,j)
             tb(i,j)=b22(i,j)
             p(i,j)=0.0
          enddo
       enddo
       do j = 1 , nhaf
          do k = 1 , nhaf
             do i = 1 , nhaf
                p(i,j)=p(i,j)+ta(i,k)*tb(k,j)
             enddo
          enddo
       enddo
       do j = 1 , nhaf
          do i = 1 , nhaf
             c(i,j)=c(i,j)-p(i,j)
             if(j+nhaf.le.ndim) then
                c(i,j+nhaf)=c(i,j+nhaf)+p(i,j)
             endif
          enddo
       enddo
c...p6
       do j = 1 , nhaf
          do i = 1 , nhaf
             ta(i,j)=a21(i,j)-a11(i,j)
             tb(i,j)=b11(i,j)+b12(i,j)
             p(i,j)=0.0
          enddo
       enddo
       do j = 1 , nhaf
          do k = 1 , nhaf
             do i = 1 , nhaf
                p(i,j)=p(i,j)+ta(i,k)*tb(k,j)
             enddo
          enddo
       enddo
       do j = 1 , nhaf
          do i = 1 , nhaf
             if(i+nhaf.le.ndim.and.j+nhaf.le.ndim) then
                c(i+nhaf,j+nhaf)=c(i+nhaf,j+nhaf)+p(i,j)
             endif
          enddo
       enddo
c...p7
       do j = 1 , nhaf
          do i = 1 , nhaf
             ta(i,j)=a12(i,j)-a22(i,j)
             tb(i,j)=b21(i,j)+b22(i,j)
             p(i,j)=0.0
          enddo
       enddo
       do j = 1 , nhaf
          do k = 1 , nhaf
             do i = 1 , nhaf
                p(i,j)=p(i,j)+ta(i,k)*tb(k,j)
             enddo
          enddo
       enddo
       do j = 1 , nhaf
          do i = 1 , nhaf
             c(i,j)=c(i,j)+p(i,j)
          enddo
       enddo
          
       return
       end
       
       
       subroutine mprod(ndim,a,b,c)
       implicit real*8 (a-h,o-z)
       dimension c(ndim,ndim)
       dimension a(ndim,ndim),b(ndim,ndim)
       
       do j = 1 , ndim
          do i = 1 , ndim
             c(i,j)=0.0
          enddo
       enddo

       do j = 1 , ndim
          do k = 1 , ndim
             do i = 1 , ndim
                c(i,j)=c(i,j)+a(i,k)*b(k,j)
             enddo
          enddo
       enddo
       
       return
       end