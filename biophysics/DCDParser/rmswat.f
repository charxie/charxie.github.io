c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c this subroutine calculates the root mean square deviations of a
c selected water from a reference structure
c because waters have a different naming convention, they are
c treated differently.
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine rmswat(nstp,ntot,id_wat,seg_wat,nref,
     :                  xref,yref,zref,rms)
      implicit real (a-h,o-z)
      parameter (maxstp=1000,natm=5000,ntil=5,nseg=10)
      character*4 resname,segname,atoname,segtype
      character hdrr*4,title*80
      character seg_wat*4
      common/pdb/iatom(natm),rx0(natm),ry0(natm),rz0(natm),
     *           ch(natm),wg(natm),iresi(natm),resname(natm),
     *           segname(natm),atoname(natm),segtype(nseg)
      common/dcd/hdrr,title(ntil),icntrl(20),rx(natm,maxstp),
     *           ry(natm,maxstp),rz(natm,maxstp)
      common/tag/ntag(natm),mtag(natm),ktag(natm)
      real rmsd(natm)
      real xref(nref),yref(nref),zref(nref)
     
      do i = 1 , ntot
         rmsd(i)=0.0
         do n = 1 , nstp
            rmsd(i)=rmsd(i)+(rx(i,n)-xref(i))*(rx(i,n)-xref(i))
            rmsd(i)=rmsd(i)+(ry(i,n)-yref(i))*(ry(i,n)-yref(i))
            rmsd(i)=rmsd(i)+(rz(i,n)-zref(i))*(rz(i,n)-zref(i))
         enddo
         rmsd(i)=rmsd(i)/real(nstp)
      enddo

      rms=0.0
      num=0

      do i = 1 , ntot
         if(iresi(i).eq.id_wat.and.segname(i).eq.seg_wat.and.
     *      atoname(i)(1:1).ne.'H') then
            if(rmsd(i).gt.100.0)stop ' error in reference coordindates'
            rms=rmsd(i)
            num=num+1
         endif
      enddo
      if(num.ne.1) stop ' failed in finding waters'
      rms=sqrt(rms)

      write(6,*) ' RMS deviations from the reference strk'
      write(6,*) ' --------------------------------------'
      write(6,'(2x,a4,5x,i5)') seg_wat,id_wat
      write(6,'(2x,a4,5x,i5,f10.5)') 'main',num,rms
      write(6,*) ' --------------------------------------'
      write(6,*) ' .'
      write(6,*) ' .'
      write(6,*) ' .'
      write(6,*) 

      return
      end