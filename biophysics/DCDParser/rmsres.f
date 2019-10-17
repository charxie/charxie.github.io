c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c this subroutine calculates the root mean square deviations of a
c selected residue from a reference structure
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine rmsres(nstp,ntot,id_res,name_res,nref,xref,yref,zref,
     *                  num_resmain,num_resside,rms_main,rms_side)
      implicit real (a-h,o-z)
      parameter (maxstp=1000,natm=5000,ntil=5,nseg=10)
      character*4 resname,segname,atoname,segtype
      character hdrr*4,title*80
      character name_res*4
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

      rms_main=0.0
      rms_side=0.0
      num_resmain=0
      num_resside=0

      do i = 1 , ntot
         if(iresi(i).eq.id_res.and.resname(i).eq.name_res.and.
     *      atoname(i)(1:1).ne.'H') then
            if(rmsd(i).gt.100.0)stop ' error in reference coordindates'
            if(ntag(i).eq.2) then
              rms_main=rms_main+rmsd(i)
              num_resmain=num_resmain+1
            endif
            if(ntag(i).eq.3) then
              rms_side=rms_side+rmsd(i)
              num_resside=num_resside+1
            endif
         endif
      enddo
      if(num_resmain.eq.0) stop ' error in finding main chain'
      rms_main=rms_main/num_resmain
      if(num_resside.gt.0) then
        rms_side=rms_side/num_resside
      else
        rms_side=0.0
      endif
      rms_main=sqrt(rms_main)
      rms_side=sqrt(rms_side)

      write(6,*) ' RMS deviations from the reference strk'
      write(6,*) ' --------------------------------------'
      write(6,'(2x,a4,5x,i5)') name_res,id_res
      write(6,'(2x,a4,5x,i5,f10.5)') 'main',num_resmain,rms_main
      write(6,'(2x,a4,5x,i5,f10.5)') 'side',num_resside,rms_side
      write(6,*) ' --------------------------------------'
      write(6,*) ' .'
      write(6,*) ' .'
      write(6,*) ' .'
      write(6,*) 

      return
      end