      program sasCal
      implicit real (a-h,o-z)
      parameter (natm=5000,nseg=10)
      character*4 dummy,resname,segname,atoname,segtype
      character hdrr*4,title*80
      common/pdb/iatom(natm),rx0(natm),ry0(natm),rz0(natm),
     *           ch(natm),wg(natm),iresi(natm),resname(natm),
     *           segname(natm),atoname(natm),segtype(nseg)

      dimension sascom(natm),saspro(natm),saslig(natm)
      dimension sas(natm)
     

      open(11,file='comp.surf.pdb',status='old',readonly)
      open(12,file='aglu.surf.pdb',status='old',readonly)
      open(13,file='hyda.surf.pdb',status='old',readonly)

      sascom1=0.0
      saspro1=0.0
      saslig1=0.0
      sascomNP=0.0
      sasproNP=0.0
      sasligNP=0.0
      iresOldARG=0
      iresOldGLU=0
      iresOldASP=0
      
c>>> read the SAS of the complex
      
      ntot=0
      ncom=0
      do i = 1 , natm
        read(11,'(a4)') dummy
        if(dummy.eq.'REMA') ncom=ncom+1
        if(dummy.eq.'ATOM') ntot=ntot+1
        if(dummy.eq.'END ') goto 100
      enddo
100   rewind 11

      do i = 1 , ncom
        read(11,'(a4)') dummy
      enddo
      if(dummy.ne.'REMA') stop 'error in reading remarks'
      do i = 1 , ntot
        read(11,'(a4,i7,1x,a4,1x,a4,i5,4x,3f8.3,2f6.2,6x,a4)') 
     *       dummy,iatom(i),atoname(i),resname(i),iresi(i),
     *       rx0(i),ry0(i),rz0(i),ch(i),wg(i),segname(i)
        if(dummy.ne.'ATOM') stop 'error in reading coordinates'
        sascom(i)=wg(i)
        sascom1=sascom1+sascom(i)
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
            goto 211
          endif        
        enddo
211     continue
        if(atoname(i)(1:1).eq.'C'.or.atoname(i)(1:1).eq.'S') then
           sascomNP=sascomNP+sascom(i)
        endif
      enddo
      close(11)


      do i = 1 , ncom
        read(12,'(a4)') dummy
      enddo
      if(dummy.ne.'REMA') stop 'error in reading remarks'
      do i = 1 , ntot
        read(12,'(a4,i7,1x,a4,1x,a4,i5,4x,3f8.3,2f6.2,6x,a4)') 
     *       dummy,iatom(i),atoname(i),resname(i),iresi(i),
     *       rx0(i),ry0(i),rz0(i),ch(i),wg(i),segname(i)
        if(dummy.ne.'ATOM') stop 'error in reading coordinates'
        saspro(i)=wg(i)
        saspro1=saspro1+saspro(i)
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
            goto 212
          endif        
        enddo
212     continue
        if(atoname(i)(1:1).eq.'C'.or.atoname(i)(1:1).eq.'S') then
           sasproNP=sasproNP+saspro(i)
        endif
      enddo
      close(12)

      do i = 1 , ncom
        read(13,'(a4)') dummy
      enddo
      if(dummy.ne.'REMA') stop 'error in reading remarks'
      do i = 1 , ntot
        read(13,'(a4,i7,1x,a4,1x,a4,i5,4x,3f8.3,2f6.2,6x,a4)') 
     *       dummy,iatom(i),atoname(i),resname(i),iresi(i),
     *       rx0(i),ry0(i),rz0(i),ch(i),wg(i),segname(i)
        if(dummy.ne.'ATOM') stop 'error in reading coordinates'
        saslig(i)=wg(i)
        saslig1=saslig1+saslig(i)
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
            goto 213
          endif        
        enddo
213     continue
        if(atoname(i)(1:1).eq.'C'.or.atoname(i)(1:1).eq.'S') then
           sasligNP=sasligNP+saslig(i)
        endif
      enddo
      close(13)
      
      tdSAS=sascom1-saspro1-saslig1
      write(6,*) sascom1, saspro1, saslig1
      write(6,*) tdSAS
      write(6,*)
      tdSASNP=sascomNP-sasproNP-sasligNP
      write(6,*) sascomNP, sasproNP, sasligNP
      write(6,*) tdSASNP
      write(6,*)
      ratio1=tdSASNP/tdSAS
      write(6,*) ratio1
      write(6,*)

      twdSAS=0.0
      twdSASnp=0.0
      tenC=16.0
      tenS=21.0
      tenN0=-6.0
      tenN1=-50.0
      tenO0=-6.0
      tenO1=-24.0
      do i = 1 , ntot
         sas(i)=sascom(i)-saspro(i)-saslig(i)
         if(atoname(i)(1:1).eq.'C') then
            twdSAS=twdSAS+tenC*sas(i)
            twdSASnp=twdSASnp+tenC*sas(i)
         endif
         if(atoname(i)(1:1).eq.'S') then
            twdSAS=twdSAS+tenS*sas(i)
            twdSASnp=twdSASnp+tenS*sas(i)
         endif
         if(atoname(i)(1:1).eq.'N') then
            if(resname(i).eq.'LYS '.and.atoname(i).eq.'NZ  ') then
               twdSAS=twdSAS+tenN1*sas(i)
            elseif(resname(i).eq.'ARG '.and.
     :             atoname(i)(1:2).eq.'NH') then
               if(iresi(i).eq.iresOldARG) then
                  if(abs(saspro(i)).gt.abs(saspro(i-1))) then
                     twdSAS=twdSAS+tenN1*sas(i)
     :                            -tenN1*sas(i-1)
     :                            +tenN0*sas(i-1)
                  else
                     twdSAS=twdSAS+tenN0*sas(i)
                  endif
               else
                  twdSAS=twdSAS+tenN1*sas(i)
                  iresOldARG=iresi(i)
               endif
            elseif(resname(i).eq.'HSP '.and.atoname(i).eq.'ND1 ') then
               twdSAS=twdSAS+tenN1*sas(i)
            else
               twdSAS=twdSAS+tenN0*sas(i)
            endif
         endif
         if(atoname(i)(1:1).eq.'O') then
            if(resname(i).eq.'GLU '.and.atoname(i)(1:2).eq.'OE') then
               if(iresi(i).eq.iresOldGLU) then
                  if(abs(saspro(i)).gt.abs(saspro(i-1))) then
                     twdSAS=twdSAS+tenO1*sas(i)
     :                            -tenO1*sas(i-1)
     :                            +tenO0*sas(i-1)
                  else
                     twdSAS=twdSAS+tenO0*sas(i)
                  endif
               else
                  twdSAS=twdSAS+tenO1*sas(i)
                  iresOldGLU=iresi(i)
               endif
            elseif(resname(i).eq.'ASP '.and.
     :             atoname(i)(1:2).eq.'OD'.and.iresi(i).ne.339) then
               if(iresi(i).eq.iresOldASP) then
                  if(abs(saspro(i)).gt.abs(saspro(i-1))) then
                     twdSAS=twdSAS+tenO1*sas(i)
     :                            -tenO1*sas(i-1)
     :                            +tenO0*sas(i-1)
                  else
                     twdSAS=twdSAS+tenO0*sas(i)
                  endif
               else
                  twdSAS=twdSAS+tenO1*sas(i)
                  iresOldASP=iresi(i)
               endif
            else
               twdSAS=twdSAS+tenO0*sas(i)
            endif
         endif
      enddo
      
      write(6,*) twdSAS, twdSASnp, twdSASnp/twdSAS

      end