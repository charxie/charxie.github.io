c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c this subroutine analyzes the distribution of interested protein atoms
c around the ligand, task list:
c
c (1) by default it should find the first 10 nearest amino acids 
c     to the targetted group and automatically single out suspicious 
c     hydrogen bond interactions;
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine prodis(nstp,ntot,target)
      implicit real (a-h,o-z)
      parameter (maxstp=1000,natm=5000,ntil=5,nseg=10,nmax=100)
      character*4 resname,segname,atoname,segtype
      character hdrr*4,title*80
      common/pdb/iatom(natm),rx0(natm),ry0(natm),rz0(natm),
     *           ch(natm),wg(natm),iresi(natm),resname(natm),
     *           segname(natm),atoname(natm),segtype(nseg)
      common/dcd/hdrr,title(ntil),icntrl(20),rx(natm,maxstp),
     *           ry(natm,maxstp),rz(natm,maxstp)
      common/tag/ntag(natm),mtag(natm),ktag(natm)
      real dmin(natm),dmax(natm),dsorted(20)
      integer isorted(20),itarget(nmax)
      character target*4
      
c>>> find out all the atoms of the targetted residue
c>>> excluding those do not seem to form interesting interactions 
c>>> with the exterior

      jadd=0
      do i = 1 , ntot      
         if(resname(i).eq.target
     *      .and.atoname(i)(1:1).ne.'C'
     *      .and.atoname(i).ne.'N1  '.and.atoname(i).ne.'N2  '
     *      .and.atoname(i).ne.'O5  ') then
            jadd=jadd+1
            itarget(jadd)=i
         endif
      enddo
      
      nman=0
      nsid=0
      do i = 1 , ntot
         mtag(i)=0
         if(ntag(i).eq.2) then  
             nman=nman+1
             mtag(i)=1
         endif
         if(ntag(i).eq.3) then
             nsid=nsid+1
             mtag(i)=2
         endif
      enddo

c>>> start to search       

      do j = 1 , jadd
      
         iorigin=itarget(j)

         do i = 1 , ntot
      
            dmin(i)=99.0
            dmax(i)=00.0
            if(mtag(i).eq.0) goto 300
         
            if(mtag(i).eq.1.or.mtag(i).eq.2) then
          
              do n = 1 , nstp
                 rrr=(rx(i,n)-rx(iorigin,n))*(rx(i,n)-rx(iorigin,n))+
     *               (ry(i,n)-ry(iorigin,n))*(ry(i,n)-ry(iorigin,n))+
     *               (rz(i,n)-rz(iorigin,n))*(rz(i,n)-rz(iorigin,n))
                 rrr=sqrt(rrr)
                 if(rrr.lt.dmin(i)) dmin(i)=rrr
                 if(rrr.gt.dmax(i)) dmax(i)=rrr
              enddo
            
            endif
         
300         continue

         enddo
        
c...find the first 10 shortest dmin(i) and store them in dsorted(10), 
c...meanwhile record their indices in isorted(10)
        
         ibeg=1
400      continue
         dsorted(ibeg)=99.0
         do i = 1 , ntot
            if(mtag(i).eq.1.or.mtag(i).eq.2) then
               if(ibeg.eq.1) then
                  if(dsorted(ibeg).gt.dmin(i)) then
                     dsorted(ibeg)=dmin(i)
                     isorted(ibeg)=i
                  endif
               else
                  if(dsorted(ibeg).gt.dmin(i).and.dmin(i).gt.
     *               dsorted(ibeg-1)+0.00001) then
                     dsorted(ibeg)=dmin(i)
                     isorted(ibeg)=i
                  endif
               endif
            endif
         enddo
         ibeg=ibeg+1
         if(ibeg.le.10) goto 400
        
         write(6,*) ' Distribution Analysis of Proteins     '
         write(6,*) ' --------------------------------------'
         write(6,'(1x,a29,a4,1x,a4)')' % first 10 nearest atoms to ' 
     *          ,resname(iorigin),atoname(iorigin)
         write(6,'(5x,a1,3x,a3,2x,a4,2x,a3,3x,a4,4x,a4)') 
     *        'i','res','dmin','seg','type','dmax'
         write(6,*) ' --------------------------------------'
         do i = 1 , 10
            if(mtag(isorted(i)).eq.1) then
               write(6,'(2i6,f6.3,2x,a4,2x,a4,2x,f6.3,2x,a4)')
     *         i,iresi(isorted(i)),dmin(isorted(i)),
     *         resname(isorted(i)),atoname(isorted(i)),
     *         dmax(isorted(i)),'main'
            else
               write(6,'(2i6,f6.3,2x,a4,2x,a4,2x,f6.3,2x,a4)')
     *         i,iresi(isorted(i)),dsorted(i),
     *         resname(isorted(i)),atoname(isorted(i)),
     *         dmax(isorted(i)),'side'
            endif
         enddo
         write(6,*) ' --------------------------------------'
         write(6,*) ' .'
         write(6,*) ' .'
         write(6,*) ' .'
         write(6,*) 
         
      enddo

      return
      end