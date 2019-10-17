      implicit real*4 (a-h,o-z)
      parameter(nmax=500,n=100)
      character bondname(nmax)*2,anti(nmax)*1
      character atomname1(nmax)*1,atomname2(nmax)*1
      dimension index1(nmax),index2(nmax)
      dimension polar1(nmax),polar2(nmax)
      dimension pchar1(nmax),pchar2(nmax)
      dimension theta1(nmax),theta2(nmax)
      dimension phian1(nmax),phian2(nmax)
      
      open(11,file='temp/temp_101',status='old')
      open(12,file='temp/pola_101',status='unknown')
      open(13,file='temp/pcha_101',status='unknown')
      open(14,file='temp/orie_101',status='unknown')
       
      do i = 1 , n
      
         read(11,1000) bondname(i),anti(i),
     @                 atomname1(i),index1(i),
     @                 polar1(i),pchar1(i),theta1(i),phian1(i),
     @                 atomname2(i),index2(i),
     @                 polar2(i),pchar2(i),theta2(i),phian2(i)
         write(12,'(i5,2f10.4)') i,polar1(i),polar2(i)
         write(13,'(i5,2f10.4)') i,pchar1(i),pchar2(i)
         write(14,'(i5,4f10.4)') i,theta1(i),phian1(i),
     @                             theta2(i),phian2(i)
      
      enddo
      
      close(11)
      
 1000 FORMAT(6X,A2,A1,2X,A1,I3,2X,F7.4,3X,F5.2,2X,F5.1,1X,F6.1,2X,3X,
     @ A1,I3,2X,F7.4,3X,F5.2,2X,F5.1,1X,F6.1,2X)                                                  

      end