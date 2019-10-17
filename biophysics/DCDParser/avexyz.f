c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c this subroutine calculates the average coordinates
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine avexyz(nstp,ntot)
      implicit real (a-h,o-z)
      parameter (maxstp=1000,natm=5000,ntil=5,nseg=10)
      character*4 resname,segname,atoname,segtype
      character hdrr*4,title*80
      common/pdb/iatom(natm),rx0(natm),ry0(natm),rz0(natm),
     *           ch(natm),wg(natm),iresi(natm),resname(natm),
     *           segname(natm),atoname(natm),segtype(nseg)
      common/dcd/hdrr,title(ntil),icntrl(20),rx(natm,maxstp),
     *           ry(natm,maxstp),rz(natm,maxstp)
      common/ave/xave(natm),yave(natm),zave(natm)
           
      do i = 1 , ntot
      
         xave(i)=0.0
         yave(i)=0.0
         zave(i)=0.0
         do n = 1 , nstp
            xave(i)=xave(i)+rx(i,n)
            yave(i)=yave(i)+ry(i,n)
            zave(i)=zave(i)+rz(i,n)
         enddo
         xave(i)=xave(i)/real(nstp)
         yave(i)=yave(i)/real(nstp)
         zave(i)=zave(i)/real(nstp)
      
      enddo

      return
      end