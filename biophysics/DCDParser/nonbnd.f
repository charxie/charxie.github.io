c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c this subroutine calculates the nonbonded interactions between the
c ligand and the environment
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine nobond(nstp,ntot,iSelected)
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
      dimension vdw(100),ele(100),vdwrms(100),elerms(100)
      
      if(iSelected.ne.0) then
        inum=0
        do i = 1 , ntot
          if(segname(i).eq.segtype(2).and.iresi(i).eq.iSelected) then
            inum=inum+1
            iend=i
          endif
        enddo
        ibeg=iend-inum+1
      else 
        ibeg=num_liga+1
        iend=ntot
      endif
      
      write(6,*) '     No.  Atom     vdw       elec'
      write(6,*) ' -----------------------------------'
      
      totvdw=0.0
      totele=0.0
      do i = 1 , num_liga
        vdw(i)=0.0
        ele(i)=0.0
        vdwrms(i)=0.0
        elerms(i)=0.0
        do j = ibeg, iend
           call dopair(nstp,i,j,v,e,vrms,erms)
           vdw(i)=vdw(i)+v
           ele(i)=ele(i)+e
           vdwrms(i)=vdwrms(i)+vrms
           elerms(i)=elerms(i)+erms
        enddo
        totvdw=totvdw+vdw(i)
        totele=totele+ele(i)
        write(6,'(i8,4x,a4,2f10.5,2f10.3)') 
     :        i,atoname(i),vdw(i),ele(i),sqrt(vdwrms(i)),sqrt(elerms(i))
      enddo
      write(6,*) ' -----------------------------------'
      if(iSelected.eq.0) then
         write(6,*)' interacting with the whole environment'
      else
         write(6,*)' interacting with residue ',iSelected,resname(ibeg)
      endif
      write(6,*) '    total vdw  total elec'
      write(6,'(3x,2f10.5)') totvdw, totele

      return
      end