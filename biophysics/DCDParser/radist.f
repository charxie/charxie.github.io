c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c this subroutine calculates the radial distribution functions of water
c molecules around some selected ligand atoms
c 'atom'  --- name of the selected ligand atom
c 'type'  --- name of the water atoms, OH2 or H*
c 'delta' --- interval of g(r)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine radist(nstp,ntot,atom,type,delta)
      implicit real (a-h,o-z)
      parameter (maxstp=1000,natm=5000,ntil=5,nseg=10)
      character*4 resname,segname,atoname,segtype
      character atom*4,filename*6,type*4
      character hdrr*4,title*80
      common/pdb/iatom(natm),rx0(natm),ry0(natm),rz0(natm),
     *           ch(natm),wg(natm),iresi(natm),resname(natm),
     *           segname(natm),atoname(natm),segtype(nseg)
      common/dcd/hdrr,title(ntil),icntrl(20),rx(natm,maxstp),
     *           ry(natm,maxstp),rz(natm,maxstp)
      common/num/num_main,num_side,num_liga,num_wcry,num_wate
      real rad(100),gor(100)
     
      if(type(1:1).eq.'O') then
         filename(1:2)='O_'
      elseif(type(1:1).eq.'H') then
         filename(1:2)='H_'
      else
         stop ' type must be oxygen or hydrogen of water molecule'
      endif
      
      filename(3:6)=atom
      
      open(91,file=filename)
                 
      icenter=0 
      do i = 1 , ntot
         if(atoname(i).eq.atom) icenter=i
      enddo 
      if(icenter.eq.0) then
         write(6,'(a7,a4,a10)') ' atom: ', atom, ' not found'
      endif
      
      do i = 1 , 100
         rad(i)=delta*float(i)
         gor(i)=0.0
      enddo
      if(rad(100).gt.20.0) then
         write(6,*) ' warning: steplength of g(r) too big'
      endif
     
      do i = 1 , ntot
         if(type.eq.'OH2 ') then
            if(segname(i).eq.'WATE'.and.atoname(i).eq.type) then
               do n = 1 , nstp
                  rrr=(rx(i,n)-rx(icenter,n))*(rx(i,n)-rx(icenter,n))
     *               +(ry(i,n)-ry(icenter,n))*(ry(i,n)-ry(icenter,n))  
     *               +(rz(i,n)-rz(icenter,n))*(rz(i,n)-rz(icenter,n))  
                  if(rrr.lt.rad(100)*rad(100)) then
                     iii=int(sqrt(rrr)/delta)+1
                     gor(iii)=gor(iii)+1.0
                  endif
               enddo
            endif
         endif
         if(type.eq.'H*  ') then
            if(segname(i).eq.'WATE'.and.atoname(i)(1:1).eq.'H') then
               do n = 1 , nstp
                  rrr=(rx(i,n)-rx(icenter,n))*(rx(i,n)-rx(icenter,n))
     *               +(ry(i,n)-ry(icenter,n))*(ry(i,n)-ry(icenter,n))  
     *               +(rz(i,n)-rz(icenter,n))*(rz(i,n)-rz(icenter,n))  
                  if(rrr.lt.rad(100)*rad(100)) then
                     iii=int(sqrt(rrr)/delta)+1
                     gor(iii)=gor(iii)+1.0
                  endif
               enddo
            endif
         endif
      enddo
      
      pi=3.1415926
      volume=4.0*pi*20**3/3.0
      do i = 1 , 100
         gor(i)=gor(i)/float(nstp)
         gor(i)=gor(i)/rad(i)**2
         gor(i)=gor(i)/delta/pi/4.0
         if(type(1:1).eq.'O') gor(i)=gor(i)*volume/(num_wate/3.0)
         if(type(1:1).eq.'H') gor(i)=gor(i)*volume/(2.0*num_wate/3.0)
      enddo

      write(91,'(2f15.5)') (rad(i),gor(i),i=1,100)

      close(91)

      return
      end