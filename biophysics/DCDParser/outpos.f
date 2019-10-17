c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c this subroutine outputs the 3d trajectory of a given atom into a file
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine outpos(nstp,ntot,segnamex,iresix,atonamex)
      implicit real (a-h,o-z)
      parameter (maxstp=1000,natm=5000,ntil=5,nseg=10)
      character*4 resname,segname,atoname,segtype
      character*4 segnamex,atonamex
      character filename*12
      character hdrr*4,title*80
      common/pdb/iatom(natm),rx0(natm),ry0(natm),rz0(natm),
     *           ch(natm),wg(natm),iresi(natm),resname(natm),
     *           segname(natm),atoname(natm),segtype(nseg)
      common/dcd/hdrr,title(ntil),icntrl(20),rx(natm,maxstp),
     *           ry(natm,maxstp),rz(natm,maxstp)
      common/tag/ntag(natm),mtag(natm),ktag(natm)
      common/num/num_main,num_side,num_liga,num_wcry,num_wate


      filename(1:4)=segnamex
      
      if(iresix.lt.10) then
         filename(5:6)='00'
         filename(7:7)=char(48+iresix)
      elseif(iresix.lt.100) then
         nn1=mod(iresix,10)
         nn2=int(iresix/10)
         filename(5:5)='0'
         filename(6:6)=char(48+nn2)
         filename(7:7)=char(48+nn1)
      elseif(iresix.lt.1000) then
         nn1=int(iresix/100)
         nnx=mod(iresix,100)
         nn2=int(nnx/10)
         nn3=mod(nnx,10)
         filename(5:5)=char(48+nn1)
         filename(6:6)=char(48+nn2)
         filename(7:7)=char(48+nn3)
      else
         write(6,*) ' warning:resid overflow!'
      endif
      
      filename(8:8)='.'
            
      filename(9:12)=atonamex      
      
      open(91,file=filename)
                 
      do i = 1 , ntot

         if(segname(i).eq.segnamex.and.iresi(i).eq.iresix.and.
     *      atoname(i).eq.atonamex) then
     
            do n = 1 , nstp
               write(91,'(3f10.5)') rx(i,n),ry(i,n),rz(i,n)
            enddo
            close(91)
            return
          
         endif               
      
      enddo

      return
      end