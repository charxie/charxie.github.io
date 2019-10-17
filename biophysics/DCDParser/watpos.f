c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c find water positions in a cube grid
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine watpos(nstp,ntot,cubeLength,numberOfGrids,hits)
      implicit real (a-h,o-z)
      parameter (maxstp=1000,natm=5000,ntil=5,nseg=10)
      character*4 resname,segname,atoname,segtype
      character hdrr*4,title*80
      common/pdb/iatom(natm),rx0(natm),ry0(natm),rz0(natm),
     *           ch(natm),wg(natm),iresi(natm),resname(natm),
     *           segname(natm),atoname(natm),segtype(nseg)
      common/dcd/hdrr,title(ntil),icntrl(20),rx(natm,maxstp),
     *           ry(natm,maxstp),rz(natm,maxstp)
      common/num/num_main,num_side,num_liga,num_wcry,num_wate
      dimension hits(numberOfGrids,numberOfGrids,numberOfGrids)
      dimension sort(100,4)
      logical newFound

      newFound=.false.
      kstep=nstp
      ic=0
      xc=0
      yc=0
      zc=0
      do i = 1 , ntot
        if(segname(i).eq.segtype(1).and.atoname(i).eq.'O8  ') ic=i
c        if(segname(i).eq.segtype(2).and.iresi(i).eq.283.and.
c     *     atoname(i).eq.'OD2 ') ic=i
c       if(segname(i).eq.segtype(2).and.iresi(i).eq.88.and.
c     *     atoname(i).eq.'OE1 ') ic=i
      enddo
      if(ic.eq.0) stop ' no central atom found!'
      do n = 1 , kstep
        xc=xc+rx(ic,n)
        yc=yc+ry(ic,n)
        zc=zc+rz(ic,n)
      enddo
      xc=xc/kstep
      yc=yc/kstep
      zc=zc/kstep
      ibeg=0
      do i = num_liga, ntot
         if(resname(i).eq.'TIP3') then
            ibeg=i
            goto 100
         endif
      enddo    
100   continue
      if(ibeg.eq.0) stop ' no TIP3 water found!'

      if(abs(cubeLength).lt.0.0001) stop ' no cube error'
      if(numberOfGrids.eq.0) stop ' no grid error'
      gridLength=cubeLength/numberOfGrids
      halfCubeLength=0.5*cubeLength
      
      do i = 1 , numberOfGrids
        do j = 1 , numberOfGrids
           do k = 1 , numberOfGrids
              hits(i,j,k)=0.0
           enddo
        enddo
      enddo
      do i = 1 , 100
        sort(i,1)=0.0
        sort(i,2)=0.0
        sort(i,3)=0.0
        sort(i,4)=0.0
      enddo
      
      do n = 1 , kstep
         do i = ibeg , ntot
           if(resname(i).eq.'TIP3'.and.atoname(i).eq.'OH2 '.and.
     *        abs(rx(i,n)-rx(ic,n)).lt.halfCubeLength.and.
     *        abs(ry(i,n)-ry(ic,n)).lt.halfCubeLength.and.
     *        abs(rz(i,n)-rz(ic,n)).lt.halfCubeLength) then
              ix=numberOfGrids/2+(rx(i,n)-rx(ic,n))/gridLength+1
              iy=numberOfGrids/2+(ry(i,n)-ry(ic,n))/gridLength+1
              iz=numberOfGrids/2+(rz(i,n)-rz(ic,n))/gridLength+1
              if(ix*iy*iz.eq.0) stop ' array overflow error'
              hits(ix,iy,iz)=hits(ix,iy,iz)+1.0
           endif
         enddo         
      enddo

      open(83,file='grid')
      open(84,file='sort')
      
      do m = 1 , 100
      do i = 1 , numberOfGrids
        do j = 1 , numberOfGrids
          do k = 1 , numberOfGrids
             if(m.eq.1) then
               if(hits(i,j,k).gt.sort(m,1)) then
                 sort(m,1)=hits(i,j,k)
                 sort(m,2)=i
                 sort(m,3)=j
                 sort(m,4)=k
               endif
             else 
               if(hits(i,j,k).gt.sort(m  ,1).and.
     *            hits(i,j,k).le.sort(m-1,1)) then
                  newFound=.true.
                  do m1 = 1 , m-1
                     if(i.eq.sort(m1,2).and.j.eq.sort(m1,3).and.
     *                  k.eq.sort(m1,4)) then
                        newFound=.false.
                     endif
                  enddo
                  if(newFound) then
                     sort(m,1)=hits(i,j,k)
                     sort(m,2)=i
                     sort(m,3)=j
                     sort(m,4)=k
                  endif
               endif
             endif
          enddo
        enddo
      enddo
      enddo

      do i = 1 , numberOfGrids
        do j = 1 , numberOfGrids
          do k = 1 , numberOfGrids
            if(abs(hits(i,j,k)).gt.0.0001) then
               rrr=sqrt(
     *         ((i-numberOfGrids/2-1)*gridLength)**2+
     *         ((j-numberOfGrids/2-1)*gridLength)**2+
     *         ((k-numberOfGrids/2-1)*gridLength)**2)
               if(rrr.lt.4.0) 
     *         write(83,'(3i5,2x,i6,3x,4f8.3)') i,j,k,hits(i,j,k),
     *         (i-numberOfGrids/2-1)*gridLength+xc,
     *         (j-numberOfGrids/2-1)*gridLength+yc,
     *         (k-numberOfGrids/2-1)*gridLength+zc,rrr
            endif
          enddo
        enddo
      enddo
      
      ktot=0
      do m = 1 , 100
        if(abs(sort(m,1)).gt.0.01) then
          ktot=ktot+sort(m,1)
          write(84,'(2i6,2x,3i5,3x,4f8.3)') m,
     *          sort(m,1),sort(m,2),sort(m,3),sort(m,4),
     *         (sort(m,2)-numberOfGrids/2-1)*gridLength+xc,
     *         (sort(m,3)-numberOfGrids/2-1)*gridLength+yc,
     *         (sort(m,4)-numberOfGrids/2-1)*gridLength+zc,
     *         sqrt(
     *         ((sort(m,2)-numberOfGrids/2-1)*gridLength)**2+
     *         ((sort(m,3)-numberOfGrids/2-1)*gridLength)**2+
     *         ((sort(m,4)-numberOfGrids/2-1)*gridLength)**2)
          write(6,'(f8.3,i8)') 
     *         sqrt(
     *         ((sort(m,2)-numberOfGrids/2-1)*gridLength)**2+
     *         ((sort(m,3)-numberOfGrids/2-1)*gridLength)**2+
     *         ((sort(m,4)-numberOfGrids/2-1)*gridLength)**2),
     *         sort(m,1)
        endif
      enddo
      write(6,*) ' total recorded hits ',ktot
      write(6,*) ' approx. no. of waters in the cube ', ktot/nstp
      
      close(83)
      close(84)
      
      

      return
      end