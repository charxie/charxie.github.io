c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c this subroutine calculates the root mean square deviations of the
c system from a reference structure
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine rmsdev(nstp,ntot,nref,xref,yref,zref)
      implicit real (a-h,o-z)
      parameter (maxstp=1000,natm=5000,ntil=5,nseg=10)
      character*4 resname,segname,atoname,segtype
      character hdrr*4,title*80
      common/pdb/iatom(natm),rx0(natm),ry0(natm),rz0(natm),
     *           ch(natm),wg(natm),iresi(natm),resname(natm),
     *           segname(natm),atoname(natm),segtype(nseg)
      common/dcd/hdrr,title(ntil),icntrl(20),rx(natm,maxstp),
     *           ry(natm,maxstp),rz(natm,maxstp)
      common/tag/ntag(natm),mtag(natm),ktag(natm)
      common/num/num_main,num_side,num_liga,num_wcry,num_wate
      real rmsd(natm),temp(20)
      real xref(nref),yref(nref),zref(nref)
      
      do i = 1 , ntot
      
         rmsd(i)=0.0
         do n = 1 , nstp
            rmsd(i)=rmsd(i)+(rx(i,n)-xref(i))*(rx(i,n)-xref(i))
            rmsd(i)=rmsd(i)+(ry(i,n)-yref(i))*(ry(i,n)-yref(i))
            rmsd(i)=rmsd(i)+(rz(i,n)-zref(i))*(rz(i,n)-zref(i))
         enddo
         rmsd(i)=rmsd(i)/real(nstp)
         rmsd(i)=sqrt(rmsd(i))
      
      enddo

      do i = 1 , 20
         temp(i)=0.0
      enddo
      
      num1=0
      num2=0
      num3=0
      num4=0
      num5=0
      num6=0
      num7=0
      do i = 1 , ntot
         if(ntag(i).eq.1.and.atoname(i)(1:1).ne.'H') then
            temp(1)=temp(1)+rmsd(i)
            num1=num1+1
         endif
         if(ntag(i).eq.2.and.atoname(i)(1:1).ne.'H') then
            temp(2)=temp(2)+rmsd(i)
            num2=num2+1
         endif
         if(ntag(i).eq.3.and.atoname(i)(1:1).ne.'H') then
            temp(3)=temp(3)+rmsd(i)
            num3=num3+1
         endif
         if(ntag(i).eq.4.and.atoname(i)(1:1).ne.'H') then
            temp(4)=temp(4)+rmsd(i)
            num4=num4+1
         endif
         if(ntag(i).eq.5.and.atoname(i)(1:1).ne.'H') then
            temp(5)=temp(5)+rmsd(i)
            num5=num5+1
         endif
         if(ktag(i).eq.1) then
           if(ntag(i).eq.2.and.atoname(i)(1:1).ne.'H') then
             temp(6)=temp(6)+rmsd(i)
             num6=num6+1
           endif
           if(ntag(i).eq.3.and.atoname(i)(1:1).ne.'H') then
             temp(7)=temp(7)+rmsd(i)
             num7=num7+1
           endif
         endif
      enddo
      temp(1)=temp(1)/num1
      temp(2)=temp(2)/num2
      temp(3)=temp(3)/num3
      temp(4)=temp(4)/num4
      temp(5)=temp(5)/num5
      temp(6)=temp(6)/num6
      temp(7)=temp(7)/num7

      write(6,*) ' RMS deviations from the reference strk'
      write(6,*) ' --------------------------------------'
      write(6,*) ' ligand                    ', temp(1)
      write(6,*) ' main chain                ', temp(2)
      write(6,*) ' side chain                ', temp(3)
      write(6,*) ' xtal water                ', temp(4)
      write(6,*) ' rest water                ', temp(5)
      write(6,*) ' main chain in active zone ', temp(6)
      write(6,*) ' side chain in active zone ', temp(7)
      write(6,*) ' --------------------------------------'
      write(6,*) ' .'
      write(6,*) ' .'
      write(6,*) ' .'
      write(6,*) 

      return
      end