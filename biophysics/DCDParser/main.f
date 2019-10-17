c ********************************************************************
c * A simple Fortran program to analyze Charmm DCD trajectory file   *
c ********************************************************************

      implicit real (a-h,o-z)
      parameter (maxstp=1000,natm=5000,ntil=5,nseg=10)
      character*4 dummy,resname,segname,atoname,segtype
      character hdrr*4,title*80
      common/pdb/iatom(natm),rx0(natm),ry0(natm),rz0(natm),
     *           ch(natm),wg(natm),iresi(natm),resname(natm),
     *           segname(natm),atoname(natm),segtype(nseg)
      common/dcd/hdrr,title(ntil),icntrl(20),rx(natm,maxstp),
     *           ry(natm,maxstp),rz(natm,maxstp)
      common/tag/ntag(natm),mtag(natm),ktag(natm)
      common/num/num_main,num_side,num_liga,num_wcry,num_wate
      common/ave/xave(natm),yave(natm),zave(natm)
      dimension xref(natm),yref(natm),zref(natm)
      dimension xc(natm),yc(natm),zc(natm)
      character nameOfRes*4
      character endCharacter*3
      dimension indexOfRes(50),nameOfRes(50)
      dimension hits(10,10,10)
      logical getrms,getdis,getpos,getint,dogrid
      
      character*80 dcdFilename,pdbFilename,refFilename

c default values

      iswitch=1
      iSelected=0
      nstp=0
      dogrid=.false.
      getint=.false.
      getrms=.false.
      getdis=.false.
      getpos=.false.
      
c read input file 

      open(77,file='input',status='old')
        read(77,*) pdbFilename
        read(77,*) refFilename
        read(77,*) dcdFilename
        read(77,*) iflag
        if(iflag.eq.1) dogrid=.true.
        read(77,*) iflag
        if(iflag.eq.1) getint=.true.
        read(77,*) iflag
        if(iflag.eq.1) getrms=.true.
        read(77,*) iflag
        if(iflag.eq.1) getdis=.true.
        read(77,*) iflag
        if(iflag.eq.1) getpos=.true.
        read(77,*) iswitch
        read(77,*) iSelected
      close(77)
      
      open(11,file=dcdFilename,
     *        access='sequential',
     *        status='old',form='unformatted',readonly)

      open(12,file=pdbFilename,status='old',readonly)
      open(13,file=refFilename,status='old',readonly)
      
c>>> read the PDB file for the initial structure
c>>> and store the reference coordinates in rx0,ry0,rz0
      
      ntot=0
      ncom=0
      do i = 1 , natm
        read(12,'(a4)') dummy
        if(dummy.eq.'REMA') ncom=ncom+1
        if(dummy.eq.'ATOM') ntot=ntot+1
        if(dummy.eq.'END ') goto 100
      enddo
100   rewind 12
      do i = 1 , ncom
        read(12,'(a4)') dummy
      enddo
      if(dummy.ne.'REMA') stop 'error in reading remarks'
      do i = 1 , ntot
        read(12,'(a4,i7,1x,a4,1x,a4,i5,4x,3f8.3,2f6.2,6x,a4)') 
     *       dummy,iatom(i),atoname(i),resname(i),iresi(i),
     *       rx0(i),ry0(i),rz0(i),ch(i),wg(i),segname(i)
        if(dummy.ne.'ATOM') stop 'error in reading coordinates'
c...assume that resname must contain at least three characters
        if(resname(i)(1:2).eq.'  ') stop 'check residue names'
        if(resname(i)(1:1).eq.' ') then
          resname(i)(1:1)=resname(i)(2:2)
          resname(i)(2:2)=resname(i)(3:3)
          resname(i)(3:3)=resname(i)(4:4)
          resname(i)(4:4)=' '
        endif
c...atom name must be written from the first position        
        do j = 1 , 3
          if(atoname(i)(1:1).eq.' ') then
            atoname(i)(1:1)=atoname(i)(2:2)
            atoname(i)(2:2)=atoname(i)(3:3)
            atoname(i)(3:3)=atoname(i)(4:4)
            atoname(i)(4:4)=' '
          else
            goto 200
          endif        
        enddo
200     continue
      enddo
      rewind 12

c... segment sequence should be ligand+protein+crystallographic water
c... +other water
      itype=1
      segtype(itype)=segname(1)
      do i = 2 , ntot
         if(segname(i).ne.segtype(itype)) then
            itype=itype+1
            segtype(itype)=segname(i)
         endif
      enddo
      radact=15.0
      do i = 1 , ntot
        if(segname(i).eq.segtype(1).and.atoname(i).eq.'C1  ') then
          icenter=i
          goto 1500
        endif
      enddo
1500  continue
      num_main=0.0
      num_side=0.0
      num_liga=0.0
      num_wcry=0.0
      num_wate=0.0
      num_acti=0.0
      do i = 1 , ntot
        if(segname(i).eq.segtype(1)) then
           num_liga=num_liga+1
           ntag(i)=1
        endif
        if(segname(i).eq.segtype(2).and.(
     *     atoname(i).eq.'C   '.or.
     *     atoname(i).eq.'O   '.or.
     *     atoname(i).eq.'N   '.or.
     *     atoname(i).eq.'HN  '.or.
     *     atoname(i).eq.'CA  '.or.
     *     atoname(i).eq.'HA  '.or.
     *     atoname(i).eq.'HA1 '.or.
     *     atoname(i).eq.'HA2 '
     *   ) ) then
           num_main=num_main+1
           ntag(i)=2
         endif
         if(segname(i).eq.segtype(2).and.(
     *      atoname(i).ne.'C   '.and.
     *      atoname(i).ne.'O   '.and.
     *      atoname(i).ne.'N   '.and.
     *      atoname(i).ne.'HN  '.and.
     *      atoname(i).ne.'CA  '.and.
     *      atoname(i).ne.'HA  '.and.
     *      atoname(i).ne.'HA1 '.and.
     *      atoname(i).ne.'HA2 '
     *   ) ) then
            num_side=num_side+1
            ntag(i)=3
         endif
         if(segname(i).eq.segtype(3).or.segname(i).eq.segtype(4)) then
            num_wcry=num_wcry+1
            ntag(i)=4
         endif
         if(segname(i).eq.segtype(5)) then
            num_wate=num_wate+1
            ntag(i)=5
         endif
         rrr=(rx0(i)-rx0(icenter))*(rx0(i)-rx0(icenter))
     *      +(ry0(i)-ry0(icenter))*(ry0(i)-ry0(icenter))
     *      +(rz0(i)-rz0(icenter))*(rz0(i)-rz0(icenter))        
         if(rrr.lt.radact*radact) then
            ktag(i)=1
            num_acti=num_acti+1
         else
            ktag(i)=0
         endif
      enddo
      write(6,*) ' Structure Information     '
      write(6,*) ' --------------------------------------'
      write(6,*) ' segments ',(' ',segtype(i),i=1,itype)
      write(6,*) ' number of total atoms     ', ntot
      write(6,*) ' number of mainchain atoms ', num_main
      write(6,*) ' number of sidechain atoms ', num_side
      write(6,*) ' number of cryst. waters   ', num_wcry
      write(6,*) ' number of other waters    ', num_wate
      write(6,*) ' number of ligand atoms    ', num_liga
      write(6,*) ' number of atoms in active zone ', num_acti
      write(6,*) ' --------------------------------------'
      write(6,*) ' .'
      write(6,*) ' .'
      write(6,*) ' .'
      write(6,*)
      close(12)     

c>>> read reference structure(e.g. xtal) and store the coordinates into
c>>> xc(),yc(),zc()

      ntot1=0
      ncom1=0      
      do i = 1 , natm
        read(13,'(a4)') dummy
        if(dummy.eq.'REMA') ncom1=ncom1+1
        if(dummy.eq.'ATOM') ntot1=ntot1+1
        if(dummy.eq.'END ') goto 300
      enddo
300   rewind 13
      if(ntot1.ne.ntot) stop ' ntot inconsistent in ref and ini pdb'
      do i = 1 , ncom1
        read(13,'(a4)') dummy
      enddo
      if(dummy.ne.'REMA') stop 'error in reading remarks'
      do i = 1 , ntot
        read(13,'(a4,26x,3f8.3)') dummy,xc(i),yc(i),zc(i)
        if(dummy.ne.'ATOM') stop 'error in reading coordinates'
      enddo
      rewind 13
      
c>>> starting to read a DCD trajectory file
     
      read(11) hdrr,icntrl
      read(11) ntitle,(title(i),i=1,ntitle) 
      read(11) natrec
      if(natrec.ne.ntot) stop ' PDB and DCD files conflict'
      if(natrec.gt.natm) stop ' number of atoms overflow!'
      nstp=icntrl(1)
      if(nstp.eq.0) stop ' read no coordinate sets'
      if(nstp.gt.1000) stop ' number of coordinate sets overflow!'
      do n = 1 , nstp
        read(11) (rx(i,n),i=1,natrec)
        read(11) (ry(i,n),i=1,natrec)
        read(11) (rz(i,n),i=1,natrec)
      enddo
      write(6,*) ' Trajectory Information    '
      write(6,*) ' --------------------------------------'
      write(6,*) ' number of total atoms     ', natrec
      write(6,*) ' number of coordinate sets ', nstp
      write(6,*) ' trajectory starts at step ', icntrl(2)
      write(6,*) ' trajectory skip interval  ', icntrl(3)
      write(6,*) ' total time length         ', icntrl(4)
      write(6,*) ' --------------------------------------'
      write(6,*) ' .'
      write(6,*) ' .'
      write(6,*) ' .'
      write(6,*) 
      close(11)
      
c>>> rms deviations from a reference structure

      open(34,file='res.rms',access='sequential')
      open(36,file='rmspro.in' ,access='sequential')
      open(38,file='rmswat.in' ,access='sequential')

      if(getrms) then
         if(iswitch.eq.1) then
            do i = 1 , ntot
               xref(i)=xc(i)
               yref(i)=yc(i)
               zref(i)=zc(i)
            enddo
            write(6,*) ' reference: crystal structure'
         endif
         if(iswitch.eq.2) then
            do i = 1 , ntot
               xref(i)=rx0(i)
               yref(i)=ry0(i)
               zref(i)=rz0(i)
            enddo
            write(6,*) ' reference: initial structure'
         endif
         if(iswitch.eq.3) then
            call avexyz(nstp,ntot)
            do i = 1 , ntot
               xref(i)=xave(i)
               yref(i)=yave(i)
               zref(i)=zave(i)
            enddo
            write(6,*) ' reference: average structure'
         endif
         call rmsdev(nstp,ntot,natm,xref,yref,zref)
         do i = 1 , 100
           read(36,'(a3)') endCharacter
           if(endCharacter.eq.'end') goto 2500
         enddo
2500     nin=i-1      
         rewind 36
         if(nin.eq.0) goto 2600
         do i = 1 , nin
           read(36,'(i3,1x,a3)') indexOfRes(i), nameOfRes(i)
         enddo
         close(36)
         do isel = 1 , nin
            call rmsres(nstp,ntot,indexOfRes(isel),nameOfRes(isel),
     *                  natm,xref,yref,zref,
     *                  num_resmain,num_resside,rms_main,rms_side)
            write(34,'(i8,2f15.5)') indexOfRes(isel),rms_main,rms_side
         enddo
2600     continue
         do i = 1 , 100
           read(38,'(a3)') endCharacter
           if(endCharacter.eq.'end') goto 3500
         enddo
3500     nin=i-1      
         rewind 38
         if(nin.eq.0) goto 3600
         do i = 1 , nin
           read(38,'(i3,1x,a4)') indexOfRes(i), nameOfRes(i)
         enddo
         close(38)
         do isel = 1 , nin
            call rmswat(nstp,ntot,indexOfRes(isel),nameOfRes(isel),
     *                  natm,xref,yref,zref,rms)
            write(34,'(i8,f15.5,5x,a4)') 
     *        indexOfRes(isel),rms,nameOfRes(isel)
         enddo
3600     continue
      endif
      
      close(34)

c>>> concentrate on the motion of waters

      if(getdis) then
         call watdis(nstp,ntot,'HYDA')
         call prodis(nstp,ntot,'HYDA')
      endif

      if(getpos) then
         
c         call outpos(nstp,ntot,'HYDA',  1,'C1  ')
c         call outpos(nstp,ntot,'HYDA',  1,'C7  ')
c         call outpos(nstp,ntot,'HYDA',  1,'O7  ')
c         call outpos(nstp,ntot,'HYDA',  1,'C8  ')
c         call outpos(nstp,ntot,'HYDA',  1,'O8  ')
c         call outpos(nstp,ntot,'HYDA',  1,'C9  ')
c         call outpos(nstp,ntot,'HYDA',  1,'N2  ')
c         call outpos(nstp,ntot,'WATA', 10,'OH2 ')
c         call outpos(nstp,ntot,'WATA', 61,'OH2 ')
c         call outpos(nstp,ntot,'WATA', 66,'OH2 ')
c         call outpos(nstp,ntot,'WATA', 76,'OH2 ')
c         call outpos(nstp,ntot,'WATA',108,'OH2 ')
c         call outpos(nstp,ntot,'WATA',111,'OH2 ')
c         call outpos(nstp,ntot,'WATA',112,'OH2 ')
c         call outpos(nstp,ntot,'WATA',139,'OH2 ')
c         call outpos(nstp,ntot,'WATE', 66,'OH2 ')
c         call outpos(nstp,ntot,'WATE', 89,'OH2 ')
         call outpos(nstp,ntot,'WATA',176,'OH2 ')
         call outpos(nstp,ntot,'WATA',200,'OH2 ')
         call outpos(nstp,ntot,'AGLU',283,'OD1 ')
         call outpos(nstp,ntot,'AGLU', 88,'OE1 ')
         call outpos(nstp,ntot,'AGLU', 88,'OE2 ')
                  
      endif

      if(getint) then
        call nobond(nstp,ntot,iSelected)
      endif

      if(dogrid) then
        call watpos(nstp,ntot,10.0,10,hits)
      endif

c>>> program ends here

      write(6,*) 
      write(6,*) ' Program Ends                          '
      write(6,*) ' --------------------------------------'


      end