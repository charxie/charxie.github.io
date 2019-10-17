      implicit real*4 (A-H,O-Z)
      parameter (natoms=500)
      character*4 symbol,resname,segname
      character head*4, string*80, title*50
      dimension idres(natoms),symbol(natoms),resname(natoms),
     :segname(natoms),x(natoms),y(natoms),z(natoms)

       open(unit=15,file='proton.pdb')
       
       read(15,'(a50)') title
       read(15,'(a80)') string
       
       Do i=1,natoms
          READ(15,140)head,symbol(i),resname(i),idres(i),
     :                x(i),y(i),z(i),
     :                segname(i)
          if(head.eq.'END ') then
             ntot=i-1
             go to 1001
          endif
       Enddo
       
 1001  continue       
  
       close(15)

       open(16,file='protein.pdb')
      
       write(16,'(a50)') title
       write(16,'(a80)') string
       
       Do i=1,ntot
          write(16,160)'ATOM',i,symbol(i),resname(i),idres(i),
     :                x(i),y(i),z(i),
     :                segname(i)
       Enddo
       write(16,'(a4)') head
      

  140  FORMAT(a4,8x,a4,1x,a4,i5,4x,3f8.3,18x,a4)
  160  FORMAT(a4,i7,1x,a4,1x,a4,i5,4x,3f8.3,18x,a4)


      END                                                                       
