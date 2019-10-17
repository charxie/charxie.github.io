c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c this subroutine calculates the nonbonded interactions between a given
c pair of atoms, atom1 must be ligand atom, atom2 must be protein or
c water atom (ensemble averages).
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine dopair(nstp,index1,index2,vdw,ele,vrms,erms)
      implicit real (a-h,o-z)
      parameter (maxstp=1000,natm=5000,ntil=5,nseg=10)
      character*4 resname,segname,atoname,segtype
      character hdrr*4,title*80
      common/pdb/iatom(natm),rx0(natm),ry0(natm),rz0(natm),
     *           ch(natm),wg(natm),iresi(natm),resname(natm),
     *           segname(natm),atoname(natm),segtype(nseg)
      common/dcd/hdrr,title(ntil),icntrl(20),rx(natm,maxstp),
     *           ry(natm,maxstp),rz(natm,maxstp)
      character nameOfAtom1*4,nameOfRes1*4,segResIn1*4
      character nameOfAtom2*4,nameOfRes2*4,segResIn2*4
      dimension vtime(maxstp),etime(maxstp)
     
      cutoffELE=30.0
      cutoffVDW=12.0
      converter=2.0**(5.0/6.0)
      charge1=0.0
      sigma1=0.0
      epsilon1=0.0
      charge2=0.0
      sigma2=0.0
      epsilon2=0.0
      nameOfAtom1=atoname(index1)
      nameOfAtom2=atoname(index2)
      idOfRes1=iresi(index1)
      idOfRes2=iresi(index2)
      nameOfRes1=resname(index1)
      nameOfRes2=resname(index2)
      segResIn1=segname(index1)
      segResIn2=segname(index2)

c>>> find out the nonbonded interaction parameters for atom1
c>>> common atoms

      if(segResIn1.ne.'HYDA'.and.segResIn1.ne.'NHYD'.and.
     *   segResIn1.ne.'OHYD'.and.segResIn1.ne.'MHYD')       
     *   stop ' The ligand must be HYDA,OHYD,NHYD,or MHYD'
      if(nameOfAtom1.eq.'N1  ') then
         if(segResIn1.eq.'OHYD') then
            charge1=-0.0845
         elseif(segResIn1.eq.'NHYD') then
            charge1=-0.111
         elseif(segResIn1.eq.'MHYD') then
            charge1=-0.22
         else
            charge1=-0.43
         endif
         epsilon1=-0.20
         sigma1=1.85
         goto 1000
      endif
      if(nameOfAtom1.eq.'C8  ') then
         charge1=0.63
         epsilon1=-0.11
         sigma1=2.0
         goto 1000
      endif
      if(nameOfAtom1.eq.'O8  '.or.nameOfAtom1.eq.'O7  ') then
         charge1=-0.51
         epsilon1=-0.12
         sigma1=1.7
         goto 1000
      endif
      if(nameOfAtom1.eq.'C1  ') then
         charge1=0.22
         epsilon1=-0.02
         sigma1=2.275
         goto 1000
      endif
      if(nameOfAtom1.eq.'C7  ') then
         charge1=0.51
         epsilon1=-0.11
         sigma1=2.0
         goto 1000
      endif
      if(nameOfAtom1.eq.'N2  ') then
         charge1=-0.47
         epsilon1=-0.2
         sigma1=1.85
         goto 1000
      endif
      if(nameOfAtom1.eq.'HN2 ') then
         charge1=0.31
         epsilon1=-0.046
         sigma1=0.2245
         goto 1000
      endif
      if(nameOfAtom1.eq.'C5  ') then
         charge1=0.25
         epsilon1=-0.02
         sigma1=2.275
         goto 1000
      endif
      if(nameOfAtom1.eq.'H5  ') then
         charge1=0.09
         epsilon1=-0.022
         sigma1=1.32
         goto 1000
      endif
      if(nameOfAtom1.eq.'O5  ') then
         charge1=-0.40
         epsilon1=-0.1521
         sigma1=1.77
         goto 1000
      endif
      if(nameOfAtom1.eq.'C2  '.or.nameOfAtom1.eq.'C3  '.or.
     *   nameOfAtom1.eq.'C4  ') then
         charge1=0.14
         epsilon1=-0.02
         sigma1=2.275
         goto 1000
      endif
      if(nameOfAtom1.eq.'H2  '.or.nameOfAtom1.eq.'H3  '.or.
     *   nameOfAtom1.eq.'H4  '.or.nameOfAtom1.eq.'H61 '.or.
     *   nameOfAtom1.eq.'H62 ') then
         charge1=0.09
         epsilon1=-0.022
         sigma1=1.32
         goto 1000
      endif
      if(nameOfAtom1.eq.'O2  '.or.nameOfAtom1.eq.'O3  '.or.
     *   nameOfAtom1.eq.'O4  '.or.nameOfAtom1.eq.'O6  ') then
         charge1=-0.66
         epsilon1=-0.1521
         sigma1=1.77
         goto 1000
      endif
      if(nameOfAtom1.eq.'HO6 '.or.nameOfAtom1.eq.'HO3 '.or.
     *   nameOfAtom1.eq.'HO4 '.or.nameOfAtom1.eq.'HO2 ') then
         charge1=0.43
         epsilon1=-0.046
         sigma1=0.2245
         goto 1000
      endif
      if(nameOfAtom1.eq.'C6  ') then
         charge1=0.05
         epsilon1=-0.02
         sigma1=2.275
         goto 1000
      endif

c>>> variable atoms      

      if(nameOfAtom1.eq.'HO9 ') then
         charge1=0.43
         epsilon1=-0.046
         sigma1=0.2245
         goto 1000
      endif
      if(nameOfAtom1.eq.'O9  ') then
         charge1=-0.4655
         epsilon1=-0.1521
         sigma1=1.77
         goto 1000
      endif
      if(nameOfAtom1.eq.'HN31'.or.nameOfAtom1.eq.'HN32') then
         charge1=0.31
         epsilon1=-0.046
         sigma1=0.2245
         goto 1000
      endif
      if(nameOfAtom1.eq.'N3  ') then
         charge1=-0.629
         epsilon1=-0.2
         sigma1=1.85
         goto 1000
      endif
      if(nameOfAtom1.eq.'HN1 ') then
         charge1=0.31
         epsilon1=-0.046
         sigma1=0.2245
         goto 1000
      endif
      if(nameOfAtom1.eq.'H91 '.or.nameOfAtom1.eq.'H92 '
     *      .or.nameOfAtom1.eq.'H93 ') then
         charge1=0.09
         epsilon1=-0.022
         sigma1=1.32
         goto 1000
      endif
      if(nameOfAtom1.eq.'C9  ') then
         charge1=-0.17
         epsilon1=-0.08
         sigma1=2.06
         goto 1000
      endif
      
1000  continue

      if(abs(charge1).lt.0.000000001) then
        write(6,'(i8,3a4)') index1,nameOfAtom1,nameOfRes1,segResIn1
        stop 'no charge found for atom1'
      endif
      if(abs(sigma1).lt.0.000000001) then
        write(6,'(i8,3a4)') index1,nameOfAtom1,nameOfRes1,segResIn1
        stop 'no sigma found for atom1'
      endif
      if(abs(epsilon1).lt.0.000000001) then
        write(6,'(i8,3a4)') index1,nameOfAtom1,nameOfRes1,segResIn1
        stop 'no epsilon found for atom1'
      endif
      if(abs(charge1*sigma1*epsilon1).lt.0.0000000001) then
        write(6,'(i8,3a4)') index1,nameOfAtom1,nameOfRes1,segResIn1
        stop ' atom1 must be a ligand atom'
      endif

      sigma1=sigma1*converter

c>>> ------------------------------------------------------------
c>>> find out the nonbonded interaction parameters for atom2
c>>> atom2 can only be an atom of amino acids or tip3 water

c>>> tip3 water
      if(nameOfRes2.eq.'TIP3') then
        if(nameOfAtom2.eq.'H1  '.or.nameOfAtom2.eq.'H2  ') then
          charge2=0.417
          epsilon2=-0.046
          sigma2=0.2245
          goto 2000
        endif
        if(nameOfAtom2.eq.'OH2 ') then
          charge2=-0.834
          epsilon2=-0.1521
          sigma2=1.7682
          goto 2000
        endif
      endif

c>>> The C=O groups are the same for all amino acids

      if(nameOfAtom2.eq.'C   ') then
        charge2=0.51
        epsilon2=-0.11
        sigma2=2.0
        goto 2000
      endif
      if(nameOfAtom2.eq.'O   ') then
        charge2=-0.51
        epsilon2=-0.12
        sigma2=1.70
        goto 2000
      endif

c>>> PROline
  
      if(nameOfRes2.eq.'PRO ') then
        if(nameOfAtom2.eq.'N   ') then
          charge2=-0.29
          epsilon2=-0.20
          sigma2=1.85
          goto 2000
        endif
        if(nameOfAtom2.eq.'CA  ') then
          charge2=0.02
          epsilon2=-0.02
          sigma2=2.275
          goto 2000
        endif
        if(nameOfAtom2.eq.'HA  ') then
          charge2=0.09
          epsilon2=-0.022
          sigma2=1.32
          goto 2000
        endif
        if(nameOfAtom2.eq.'CD  ') then
          charge2=0.00
          epsilon2=-0.055
          sigma2=2.175
          goto 2000
        endif
        if(nameOfAtom2.eq.'HD1 '.or.nameOfAtom2.eq.'HD2 '.or.
     *     nameOfAtom2.eq.'HG1 '.or.nameOfAtom2.eq.'HG2 '.or.
     *     nameOfAtom2.eq.'HB1 '.or.nameOfAtom2.eq.'HB2 ') then
           charge2=0.09
           epsilon2=-0.022
           sigma2=1.32
           goto 2000
        endif
        if(nameOfAtom2.eq.'CB  '.or.nameOfAtom2.eq.'CG  ') then
          charge2=-0.18
          epsilon2=-0.055
          sigma2=2.175
          goto 2000
        endif
      endif

c>>>GLYcine

      if(nameOfRes2.eq.'GLY ') then
        if(nameOfAtom2.eq.'CA  ') then
          charge2=-0.02
          epsilon2=-0.055
          sigma2=2.175
          goto 2000
        endif
        if(nameOfAtom2.eq.'HA1 '.or.nameOfAtom2.eq.'HA2 ') then
          charge2=0.09
          epsilon2=-0.022
          sigma2=1.32
          goto 2000
        endif
        if(nameOfAtom2.eq.'N   ') then
          charge2=-0.47
          epsilon2=-0.20
          sigma2=1.85
          goto 2000
        endif
        if(nameOfAtom2.eq.'HN  ') then
          charge2=0.31
          epsilon2=-0.046
          sigma2=0.2245
          goto 2000
        endif
      endif

c>>> main chain atoms and common side chain atoms     

      if(nameOfAtom2.eq.'N   ') then
        charge2=-0.47
        epsilon2=-0.20
        sigma2=1.85
        goto 2000
      endif
      if(nameOfAtom2.eq.'HN  ') then
        charge2=0.31
        epsilon2=-0.046
        sigma2=0.2245
        goto 2000
      endif
      if(nameOfAtom2.eq.'CA  ') then
        charge2=0.07
        epsilon2=-0.02
        sigma2=2.275
        goto 2000
      endif
      if(nameOfAtom2.eq.'HA  ') then
        charge2=0.09
        epsilon2=-0.022
        sigma2=1.32
        goto 2000
      endif
      if(nameOfAtom2.eq.'HB1 '.or.
     *   nameOfAtom2.eq.'HB2 '.or.
     *   nameOfAtom2.eq.'HB3 ') then
         charge2=0.09
         epsilon2=-0.022
         sigma2=1.32
         goto 2000
      endif

c>>> side chain atoms that are different
      
      if(nameOfRes2.eq.'ALA ') then
        if(nameOfAtom2.eq.'CB  ') then
          charge2=-0.27
          epsilon2=-0.08
          sigma2=2.06
          goto 2000
        endif
      endif

      if(nameOfRes2.eq.'ARG ') then
        if(nameOfAtom2.eq.'CB  '.or.nameOfAtom2.eq.'CG  ') then
          charge2=-0.18
          epsilon2=-0.055
          sigma2=2.175
          goto 2000
        endif
        if(nameOfAtom2.eq.'HG1 '.or.nameOfAtom2.eq.'HG2 '.or.
     *     nameOfAtom2.eq.'HD1 '.or.nameOfAtom2.eq.'HD2 ') then
           charge2=0.09
           epsilon2=-0.022
           sigma2=1.32
           goto 2000
        endif
        if(nameOfAtom2.eq.'CD  ') then
          charge2=0.20
          epsilon2=-0.055
          sigma2=2.175
          goto 2000
        endif
        if(nameOfAtom2.eq.'NE  ') then
          charge2=-0.70
          epsilon2=-0.20
          sigma2=1.85
          goto 2000
        endif
        if(nameOfAtom2.eq.'HE  ') then
          charge2=0.44
          epsilon2=-0.046
          sigma2=0.2245
          goto 2000
        endif
        if(nameOfAtom2.eq.'CZ  ') then
          charge2=0.64
          epsilon2=-0.11
          sigma2=2.0
          goto 2000
        endif
        if(nameOfAtom2.eq.'NH1 '.or.nameOfAtom2.eq.'NH2 ') then
          charge2=-0.80
          epsilon2=-0.20
          sigma2=1.85
          goto 2000
        endif
        if(nameOfAtom2.eq.'HH11'.or.nameOfAtom2.eq.'HH12'.or.
     *     nameOfAtom2.eq.'HH21'.or.nameOfAtom2.eq.'HH22') then
          charge2=0.46
          epsilon2=-0.046
          sigma2=0.2245
          goto 2000
        endif
      endif

      if(nameOfRes2.eq.'ASN ') then
        if(nameOfAtom2.eq.'CB  ') then
          charge2=-0.18
          epsilon2=-0.055
          sigma2=2.175
          goto 2000
        endif
        if(nameOfAtom2.eq.'CG  ') then
          charge2=0.55
          epsilon2=-0.07
          sigma2=2.0
          goto 2000
        endif
        if(nameOfAtom2.eq.'OD1 ') then
          charge2=-0.55
          epsilon2=-0.12
          sigma2=1.70
          goto 2000
        endif
        if(nameOfAtom2.eq.'ND2 ') then
          charge2=-0.62
          epsilon2=-0.20
          sigma2=1.85
          goto 2000
        endif
        if(nameOfAtom2.eq.'HD21') then
          charge2=0.32
          epsilon2=-0.046
          sigma2=0.2245
          goto 2000
        endif
        if(nameOfAtom2.eq.'HD22') then
          charge2=0.30
          epsilon2=-0.046
          sigma2=0.2245
          goto 2000
        endif
      endif

      if(nameOfRes2.eq.'ASP ') then
        if(nameOfAtom2.eq.'CB  ') then
          charge2=-0.28
          epsilon2=-0.055
          sigma2=2.175
          goto 2000
        endif
        if(nameOfAtom2.eq.'CG  ') then
          charge2=0.62
          epsilon2=-0.07
          sigma2=2.0
          goto 2000
        endif
        if(nameOfAtom2.eq.'OD1 '.or.nameOfAtom2.eq.'OD2 ') then
          charge2=-0.76
          epsilon2=-0.12
          sigma2=1.70
          goto 2000
        endif
      endif

      if(nameOfRes2.eq.'CYS ') then
        if(nameOfAtom2.eq.'CB  ') then
          charge2=-0.11
          epsilon2=-0.055
          sigma2=2.175
          goto 2000
        endif
        if(nameOfAtom2.eq.'SG  ') then
          charge2=-0.23
          epsilon2=-0.45
          sigma2=2.0
          goto 2000
        endif
        if(nameOfAtom2.eq.'HG1 ') then
          charge2=0.16
          epsilon2=-0.10
          sigma2=0.45
          goto 2000
        endif
      endif

      if(nameOfRes2.eq.'GLN ') then
        if(nameOfAtom2.eq.'CB  ') then
          charge2=-0.18
          epsilon2=-0.055
          sigma2=2.175
          goto 2000
        endif
        if(nameOfAtom2.eq.'CG  ') then
          charge2=-0.18
          epsilon2=-0.055
          sigma2=2.175
          goto 2000
        endif
        if(nameOfAtom2.eq.'HG1 '.or.nameOfAtom2.eq.'HG2 ') then
          charge2=0.09
          epsilon2=-0.022
          sigma2=1.32
          goto 2000
        endif
        if(nameOfAtom2.eq.'CD  ') then
          charge2=0.55
          epsilon2=-0.07
          sigma2=2.00
          goto 2000
        endif
        if(nameOfAtom2.eq.'OE1 ') then
          charge2=-0.55
          epsilon2=-0.12
          sigma2=1.70
          goto 2000
        endif
        if(nameOfAtom2.eq.'NE2 ') then
          charge2=-0.62
          epsilon2=-0.20
          sigma2=1.85
          goto 2000
        endif
        if(nameOfAtom2.eq.'HE21') then
          charge2=0.32
          epsilon2=-0.046
          sigma2=0.2245
          goto 2000
        endif
        if(nameOfAtom2.eq.'HE22') then
          charge2=0.30
          epsilon2=-0.046
          sigma2=0.2245
          goto 2000
        endif
      endif

      if(nameOfRes2.eq.'GLU ') then
        if(nameOfAtom2.eq.'CB  ') then
          charge2=-0.18
          epsilon2=-0.055
          sigma2=2.175
          goto 2000
        endif
        if(nameOfAtom2.eq.'CG  ') then
          charge2=-0.28
          epsilon2=-0.055
          sigma2=2.175
          goto 2000
        endif
        if(nameOfAtom2.eq.'HG1 '.or.nameOfAtom2.eq.'HG2 ') then
          charge2=0.09
          epsilon2=-0.022
          sigma2=1.32
          goto 2000
        endif
        if(nameOfAtom2.eq.'CD  ') then
          charge2=0.62
          epsilon2=-0.07
          sigma2=2.0
          goto 2000
        endif
        if(nameOfAtom2.eq.'OE1 '.or.nameOfAtom2.eq.'OE2 ') then
          charge2=-0.76
          epsilon2=-0.12
          sigma2=1.70
          goto 2000
        endif
      endif

      if(nameOfRes2.eq.'HSD ') then
        if(nameOfAtom2.eq.'ND1 ') then
          charge2=-0.36
          epsilon2=-0.20
          sigma2=1.85
          goto 2000
        endif
        if(nameOfAtom2.eq.'HD1 ') then
          charge2=0.32
          epsilon2=-0.046
          sigma2=0.2245
          goto 2000
        endif
        if(nameOfAtom2.eq.'CG  ') then
          charge2=-0.05
          epsilon2=-0.05
          sigma2=1.80
          goto 2000
        endif
        if(nameOfAtom2.eq.'CB  ') then
          charge2=-0.09
          epsilon2=-0.055
          sigma2=2.175
          goto 2000
        endif
        if(nameOfAtom2.eq.'NE2 ') then
          charge2=-0.70
          epsilon2=-0.20
          sigma2=1.85
          goto 2000
        endif
        if(nameOfAtom2.eq.'CD2 ') then
          charge2=0.22
          epsilon2=-0.05
          sigma2=1.80
          goto 2000
        endif
        if(nameOfAtom2.eq.'HD2 ') then
          charge2=0.10
          epsilon2=-0.0078
          sigma2=1.468
          goto 2000
        endif
        if(nameOfAtom2.eq.'CE1 ') then
          charge2=0.25
          epsilon2=-0.05
          sigma2=1.80
          goto 2000
        endif
        if(nameOfAtom2.eq.'HE1 ') then
          charge2=0.13
          epsilon2=-0.046
          sigma2=0.90
          goto 2000
        endif
      endif

      if(nameOfRes2.eq.'HSE ') then
        if(nameOfAtom2.eq.'NE2 ') then
          charge2=-0.36
          epsilon2=-0.20
          sigma2=1.85
          goto 2000
        endif
        if(nameOfAtom2.eq.'HE2 ') then
          charge2=0.32
          epsilon2=-0.046
          sigma2=0.2245
          goto 2000
        endif
        if(nameOfAtom2.eq.'CD2 ') then
          charge2=-0.05
          epsilon2=-0.05
          sigma2=1.80
        endif
        if(nameOfAtom2.eq.'HD2 ') then
          charge2=0.09
          epsilon2=-0.0078
          sigma2=1.468
          goto 2000
        endif
        if(nameOfAtom2.eq.'ND1 ') then
          charge2=-0.70
          epsilon2=-0.20
          sigma2=1.85
          goto 2000
        endif
        if(nameOfAtom2.eq.'CG  ') then
          charge2=0.22
          epsilon2=-0.05
          sigma2=1.80
          goto 2000
        endif
        if(nameOfAtom2.eq.'CE1 ') then
          charge2=0.25
          epsilon2=-0.05
          sigma2=1.80
          goto 2000
        endif
        if(nameOfAtom2.eq.'HE1 ') then
          charge2=0.13
          epsilon2=-0.046
          sigma2=0.90
          goto 2000
        endif
        if(nameOfAtom2.eq.'CB  ') then
          charge2=-0.08
          epsilon2=-0.055
          sigma2=2.175
          goto 2000
        endif
      endif
      
      if(nameOfRes2.eq.'HSP ') then
        if(nameOfAtom2.eq.'ND1 '.or.nameOfAtom2.eq.'NE2 ') then
          charge2=-0.51
          epsilon2=-0.20
          sigma2=1.85
          goto 2000
        endif
        if(nameOfAtom2.eq.'HD1 '.or.nameOfAtom2.eq.'HE2 ') then
          charge2=0.44
          epsilon2=-0.046
          sigma2=0.2245
          goto 2000
        endif
        if(nameOfAtom2.eq.'CE1 ') then
          charge2=0.32
          epsilon2=-0.05
          sigma2=1.80
          goto 2000
        endif
        if(nameOfAtom2.eq.'HE1 ') then
          charge2=0.18
          epsilon2=-0.046
          sigma2=0.70
          goto 2000
        endif
        if(nameOfAtom2.eq.'CD2 ') then
          charge2=0.19
          epsilon2=-0.05
          sigma2=1.80
          goto 2000
        endif
        if(nameOfAtom2.eq.'HD2 ') then
          charge2=0.13
          epsilon2=-0.046
          sigma2=0.90
          goto 2000
        endif
        if(nameOfAtom2.eq.'CG  ') then
          charge2=0.19
          epsilon2=-0.05
          sigma2=1.80
          goto 2000
        endif
        if(nameOfAtom2.eq.'CB  ') then
          charge2=-0.05
          epsilon2=-0.055
          sigma2=2.175
          goto 2000
        endif
      endif

      if(nameOfRes2.eq.'ILE ') then
        if(nameOfAtom2.eq.'CB  ') then
          charge2=-0.09
          epsilon2=-0.02
          sigma2=2.275
          goto 2000
        endif
        if(nameOfAtom2.eq.'HB  ') then
          charge2=0.09
          epsilon2=-0.022
          sigma2=1.32
          goto 2000
        endif
        if(nameOfAtom2.eq.'CG2 '.or.nameOfAtom2.eq.'CD  ') then
          charge2=-0.27
          epsilon2=-0.08
          sigma2=2.06
          goto 2000
        endif
        if(nameOfAtom2.eq.'CG1 ') then
          charge2=-0.18
          epsilon2=-0.055
          sigma2=2.175
          goto 2000
        endif
        if(nameOfAtom2.eq.'HG21'.or.nameOfAtom2.eq.'HG22'.or.
     *     nameOfAtom2.eq.'HG23'.or.nameOfAtom2.eq.'HG11'.or.
     *     nameOfAtom2.eq.'HG12'.or.nameOfAtom2.eq.'HD1 '.or.
     *     nameOfAtom2.eq.'HD2 '.or.nameOfAtom2.eq.'HD3 ') then
          charge2=0.09
          epsilon2=-0.022
          sigma2=1.32
          goto 2000
        endif
      endif

      if(nameOfRes2.eq.'LEU ') then
        if(nameOfAtom2.eq.'CB  ') then
          charge2=-0.18
          epsilon2=-0.055
          sigma2=2.175
          goto 2000
        endif
        if(nameOfAtom2.eq.'CG  ') then
          charge2=-0.09
          epsilon2=-0.02
          sigma2=2.275
          goto 2000
        endif
        if(nameOfAtom2.eq.'CD1 '.or.nameOfAtom2.eq.'CD2 ') then
          charge2=-0.27
          epsilon2=-0.08
          sigma2=2.06
          goto 2000
        endif
        if(nameOfAtom2.eq.'HG  '.or.
     *     nameOfAtom2.eq.'HD11'.or.nameOfAtom2.eq.'HD12'.or.
     *     nameOfAtom2.eq.'HD13'.or.nameOfAtom2.eq.'HD21'.or.
     *     nameOfAtom2.eq.'HD22'.or.nameOfAtom2.eq.'HD23') then
          charge2=0.09
          epsilon2=-0.022
          sigma2=1.32
          goto 2000
        endif
      endif

      if(nameOfRes2.eq.'LYS ') then
        if(nameOfAtom2.eq.'CE  ') then
          charge2=0.21
          epsilon2=-0.055
          sigma2=2.175
          goto 2000
        endif
        if(nameOfAtom2.eq.'CD  '.or.nameOfAtom2.eq.'CG  '.or.
     *     nameOfAtom2.eq.'CB  ') then
          charge2=-0.18
          epsilon2=-0.055
          sigma2=2.175
          goto 2000
        endif
        if(nameOfAtom2.eq.'HE1 '.or.nameOfAtom2.eq.'HE2 ') then
          charge2=0.05
          epsilon2=-0.022
          sigma2=1.32
          goto 2000
        endif
        if(nameOfAtom2.eq.'HG1 '.or.nameOfAtom2.eq.'HG2 '.or.
     *     nameOfAtom2.eq.'HD1 '.or.nameOfAtom2.eq.'HD2 ') then
           charge2=0.09
           epsilon2=-0.022
           sigma2=1.32
           goto 2000
        endif
        if(nameOfAtom2.eq.'NZ  ') then
          charge2=-0.30
          epsilon2=-0.20
          sigma2=1.85
          goto 2000
        endif
        if(nameOfAtom2.eq.'HZ1 '.or.nameOfAtom2.eq.'HZ2 '.or.
     *     nameOfAtom2.eq.'HZ3 ') then
           charge2=0.33
           epsilon2=-0.046
           sigma2=0.2245
           goto 2000
        endif
      endif

      if(nameOfRes2.eq.'MET ') then
        if(nameOfAtom2.eq.'CB  ') then
          charge2=-0.18
          epsilon2=-0.055
          sigma2=2.175
          goto 2000
        endif
        if(nameOfAtom2.eq.'CG  ') then
          charge2=-0.14
          epsilon2=-0.055
          sigma2=2.175
          goto 2000
        endif
        if(nameOfAtom2.eq.'HG1 '.or.nameOfAtom2.eq.'HG2 '.or.
     *     nameOfAtom2.eq.'HE1 '.or.nameOfAtom2.eq.'HE2 '.or.
     *     nameOfAtom2.eq.'HE3 ') then
           charge2=0.09
           epsilon2=-0.022
           sigma2=1.32
           goto 2000
        endif
        if(nameOfAtom2.eq.'SD  ') then
          charge2=-0.09
          epsilon2=-0.45
          sigma2=2.0
          goto 2000
        endif
        if(nameOfAtom2.eq.'CE  ') then
          charge2=-0.22
          epsilon2=-0.08
          sigma2=2.06
          goto 2000
        endif
      endif

      if(nameOfRes2.eq.'PHE ') then
        if(nameOfAtom2.eq.'CB  ') then
          charge2=-0.18
          epsilon2=-0.055
          sigma2=2.175
          goto 2000
        endif
        if(nameOfAtom2.eq.'CG  ') then
          charge2=0.0
          epsilon2=-0.07
          sigma2=1.9924
          goto 2000
        endif
        if(nameOfAtom2.eq.'CD1 '.or.nameOfAtom2.eq.'CD2 '.or.
     *     nameOfAtom2.eq.'CE1 '.or.nameOfAtom2.eq.'CE2 '.or.
     *     nameOfAtom2.eq.'CZ  ') then
           charge2=-0.115
           epsilon2=-0.07
           sigma2=1.9924
           goto 2000
        endif
        if(nameOfAtom2.eq.'HD1 '.or.nameOfAtom2.eq.'HD2 '.or.
     *     nameOfAtom2.eq.'HE1 '.or.nameOfAtom2.eq.'HE2 '.or.
     *     nameOfAtom2.eq.'HZ  ') then
           charge2=0.115
           epsilon2=-0.03
           sigma2=1.3582
           goto 2000
        endif
      endif

      if(nameOfRes2.eq.'SER ') then
        if(nameOfAtom2.eq.'CB  ') then
          charge2=0.05
          epsilon2=-0.055
          sigma2=2.175
          goto 2000
        endif
        if(nameOfAtom2.eq.'OG  ') then
          charge2=-0.66
          epsilon2=-0.1521
          sigma2=1.77
          goto 2000
        endif
        if(nameOfAtom2.eq.'HG1 ') then
          charge2=0.43
          epsilon2=-0.046
          sigma2=0.2245
          goto 2000
        endif
      endif

      if(nameOfRes2.eq.'THR ') then
        if(nameOfAtom2.eq.'CB  ') then
          charge2=0.14
          epsilon2=-0.02
          sigma2=2.275
          goto 2000
        endif
        if(nameOfAtom2.eq.'HB  ') then
          charge2=0.09
          epsilon2=-0.022
          sigma2=1.32
          goto 2000
        endif
        if(nameOfAtom2.eq.'OG1 ') then
          charge2=-0.66
          epsilon2=-0.1521
          sigma2=1.77
          goto 2000
        endif
        if(nameOfAtom2.eq.'HG1 ') then
          charge2=0.43
          epsilon2=-0.046
          sigma2=0.2245
          goto 2000
        endif
        if(nameOfAtom2.eq.'CG2 ') then
          charge2=-0.27
          epsilon2=-0.08
          sigma2=2.06
          goto 2000
        endif
        if(nameOfAtom2.eq.'HG21'.or.nameOfAtom2.eq.'HG22'.or.
     *     nameOfAtom2.eq.'HG23') then
           charge2=0.09
           epsilon2=-0.022
           sigma2=1.32
          goto 2000
        endif
      endif

      if(nameOfRes2.eq.'TRP ') then
        if(nameOfAtom2.eq.'CB  ') then
          charge2=-0.18
          epsilon2=-0.055
          sigma2=2.175
          goto 2000
        endif
        if(nameOfAtom2.eq.'CG  ') then
          charge2=-0.03
          epsilon2=-0.07
          sigma2=1.9924
          goto 2000
        endif
        if(nameOfAtom2.eq.'CD1 ') then
          charge2=0.035
          epsilon2=-0.07
          sigma2=1.9924
          goto 2000
        endif
        if(nameOfAtom2.eq.'CD2 ') then
          charge2=-0.02
          epsilon2=-0.09
          sigma2=1.80
          goto 2000
        endif
        if(nameOfAtom2.eq.'HD1 ') then
          charge2=0.115
          epsilon2=-0.03
          sigma2=1.3582
          goto 2000
        endif
        if(nameOfAtom2.eq.'NE1 ') then
          charge2=-0.61
          epsilon2=-0.20
          sigma2=1.85
          goto 2000
        endif
        if(nameOfAtom2.eq.'HE1 ') then
          charge2=0.38
          epsilon2=-0.046
          sigma2=0.2245
          goto 2000
        endif
        if(nameOfAtom2.eq.'CE2 ') then
          charge2=0.13
          epsilon2=-0.09
          sigma2=1.80
          goto 2000
        endif
        if(nameOfAtom2.eq.'CE3 '.or.nameOfAtom2.eq.'CZ2 '.or.
     *     nameOfAtom2.eq.'CZ3 '.or.nameOfAtom2.eq.'CH2 ') then
           charge2=-0.115
           epsilon2=-0.07
           sigma2=1.9924
          goto 2000
        endif
        if(nameOfAtom2.eq.'HE3 '.or.nameOfAtom2.eq.'HZ2 '.or.
     *     nameOfAtom2.eq.'HZ3 '.or.nameOfAtom2.eq.'HH2 ') then
           charge2=0.115
           epsilon2=-0.03
           sigma2=1.3582
          goto 2000
        endif
      endif

      if(nameOfRes2.eq.'TYR ') then
        if(nameOfAtom2.eq.'CB  ') then
          charge2=-0.18
          epsilon2=-0.055
          sigma2=2.175
          goto 2000
        endif
        if(nameOfAtom2.eq.'CG  ') then
          charge2=0.0
          epsilon2=-0.07
          sigma2=1.9924
          goto 2000
        endif
        if(nameOfAtom2.eq.'CZ  ') then
          charge2=0.11
          epsilon2=-0.07
          sigma2=1.9924
          goto 2000
        endif
        if(nameOfAtom2.eq.'CD1 '.or.nameOfAtom2.eq.'CD2 '.or.
     *     nameOfAtom2.eq.'CE1 '.or.nameOfAtom2.eq.'CE2 ') then
           charge2=-0.115
           epsilon2=-0.07
           sigma2=1.9924
          goto 2000
        endif
        if(nameOfAtom2.eq.'HD1 '.or.nameOfAtom2.eq.'HD2 '.or.
     *     nameOfAtom2.eq.'HE1 '.or.nameOfAtom2.eq.'HE2 ') then
           charge2=0.115
           epsilon2=-0.03
           sigma2=1.3582
          goto 2000
        endif
        if(nameOfAtom2.eq.'OH  ') then
          charge2=-0.54
          epsilon2=-0.1521
          sigma2=1.77
          goto 2000
        endif
        if(nameOfAtom2.eq.'HH  ') then
          charge2=0.43
          epsilon2=-0.046
          sigma2=0.2245
          goto 2000
        endif
      endif

      if(nameOfRes2.eq.'VAL ') then
        if(nameOfAtom2.eq.'CB  ') then
          charge2=-0.09
          epsilon2=-0.02
          sigma2=2.275
          goto 2000
        endif
        if(nameOfAtom2.eq.'HB  ') then
          charge2=0.09
          epsilon2=-0.022
          sigma2=1.32
          goto 2000
        endif
        if(nameOfAtom2.eq.'CG1 '.or.nameOfAtom2.eq.'CG2 ') then
          charge2=-0.27
          epsilon2=-0.08
          sigma2=2.06
          goto 2000
        endif
        if(nameOfAtom2.eq.'HG11'.or.nameOfAtom2.eq.'HG12'.or.
     *     nameOfAtom2.eq.'HG13'.or.nameOfAtom2.eq.'HG21'.or.
     *     nameOfAtom2.eq.'HG22'.or.nameOfAtom2.eq.'HG23') then
          charge2=0.09
          epsilon2=-0.022
          sigma2=1.32
          goto 2000
        endif
      endif
      
c>>> patches

      if(nameOfAtom2.eq.'HD2 '.and.nameOfRes2.eq.'ASP ') then
          charge2=0.44
          epsilon2=-0.046
          sigma2=0.2245
          goto 2000
      endif      
      
      if(nameOfRes2.eq.'LSPH') then
        if(nameOfAtom2.eq.'CB  '.or.nameOfAtom2.eq.'CG  '.or.
     *     nameOfAtom2.eq.'CD  ') then
           charge2=-0.18
           epsilon2=-0.055
           sigma2=2.175
          goto 2000
        endif
        if(nameOfAtom2.eq.'HG1 '.or.nameOfAtom2.eq.'HG2 '.or.
     *     nameOfAtom2.eq.'HD1 '.or.nameOfAtom2.eq.'HD2 '.or.
     *     nameOfAtom2.eq.'HE1 '.or.nameOfAtom2.eq.'HE2 '.or.
     *     nameOfAtom2.eq.'H2A1'.or.nameOfAtom2.eq.'H2A2'.or.
     *     nameOfAtom2.eq.'H2A3') then
           charge2=0.09
           epsilon2=-0.022
           sigma2=1.32
          goto 2000
        endif
        if(nameOfAtom2.eq.'CE  ') then
          charge2=0.066
          epsilon2=-0.1562
          sigma2=1.8
          goto 2000
        endif
        if(nameOfAtom2.eq.'NZ  ') then
          charge2=-0.696
          epsilon2=-0.2
          sigma2=1.85
          goto 2000
        endif
        if(nameOfAtom2.eq.'N1  ') then
          charge2=-0.62
          epsilon2=-0.2
          sigma2=1.85
          goto 2000
        endif
        if(nameOfAtom2.eq.'C2  ') then
          charge2=0.31
          epsilon2=-0.18
          sigma2=1.8
          goto 2000
        endif
        if(nameOfAtom2.eq.'C2A ') then
          charge2=-0.27
          epsilon2=-0.055
          sigma2=2.175
          goto 2000
        endif
        if(nameOfAtom2.eq.'C3 ') then
          charge2=0.11
          epsilon2=-0.18
          sigma2=1.8
          goto 2000
        endif
        if(nameOfAtom2.eq.'O3  ') then
          charge2=-0.45
          epsilon2=-0.12
          sigma2=1.7
          goto 2000
        endif
        if(nameOfAtom2.eq.'H3  ') then
          charge2=0.34
          epsilon2=-0.046
          sigma2=0.2245
          goto 2000
        endif
        if(nameOfAtom2.eq.'C4  ') then
          charge2=0.0862
          epsilon2=-0.18
          sigma2=1.8
          goto 2000
        endif
        if(nameOfAtom2.eq.'C4A ') then
          charge2=0.3038
          epsilon2=-0.11
          sigma2=2.1
          goto 2000
        endif
        if(nameOfAtom2.eq.'H4A ') then
          charge2=0.06
          epsilon2=-0.0078
          sigma2=0.87
          goto 2000
        endif
        if(nameOfAtom2.eq.'C5  ') then
          charge2=0.0
          epsilon2=-0.18
          sigma2=1.8
          goto 2000
        endif
        if(nameOfAtom2.eq.'C6  ') then
          charge2=0.195
          epsilon2=-0.18
          sigma2=1.8
          goto 2000
        endif
        if(nameOfAtom2.eq.'H6  ') then
          charge2=0.115
          epsilon2=-0.046
          sigma2=0.9
          goto 2000
        endif
        if(nameOfAtom2.eq.'C5A ') then
          charge2=-0.08
          epsilon2=-0.055
          sigma2=2.175
          goto 2000
        endif
        if(nameOfAtom2.eq.'H5A1'.or.nameOfAtom2.eq.'H5A2') then
          charge2=0.09
          epsilon2=-0.022
          sigma2=1.32
          goto 2000
        endif
        if(nameOfAtom2.eq.'OP4 ') then
          charge2=-0.62
          epsilon2=-0.1521
          sigma2=1.77
          goto 2000
        endif
        if(nameOfAtom2.eq.'P   ') then
          charge2=1.50
          epsilon2=-0.585
          sigma2=2.15
          goto 2000
        endif
        if(nameOfAtom2.eq.'OP1 '.or.nameOfAtom2.eq.'OP2 ') then
          charge2=-0.82
          epsilon2=-0.12
          sigma2=1.7
          goto 2000
        endif
        if(nameOfAtom2.eq.'OP3 ') then
          charge2=-0.68
          epsilon2=-0.1521
          sigma2=1.77
          goto 2000
        endif
        if(nameOfAtom2.eq.'HP  ') then
          charge2=0.34
          epsilon2=-0.046
          sigma2=0.2245
          goto 2000
        endif
      endif

2000  continue

      sigma2=sigma2*converter
      
c      if(abs(charge2).lt.0.000000001) then
c        write(6,'(i8,3a4)') index2,nameOfAtom2,nameOfRes2,segResIn2
c        stop 'no charge found for atom2'
c      endif
      if(abs(sigma2).lt.0.000000001) then
        write(6,'(i8,3a4)') index2,nameOfAtom2,nameOfRes2,segResIn2
        stop 'no sigma found for atom2'
      endif
      if(abs(epsilon2).lt.0.000000001) then
        write(6,'(i8,3a4)') index2,nameOfAtom2,nameOfRes2,segResIn2
        stop 'no epsilon found for atom2'
      endif
c      if(abs(charge2*sigma2*epsilon2).lt.0.0000000001) then
c        write(6,'(i8,3a4)') index2,nameOfAtom2,nameOfRes2,segResIn2
c        stop ' atom1 must be a protein or water atom'
c      endif


c>>> calculate the vdw and electrostatic interactions

      vdw=0.0
      ele=0.0
      do n = 1 , nstp
        rrr=(rx(index1,n)-rx(index2,n))*(rx(index1,n)-rx(index2,n))
     *     +(ry(index1,n)-ry(index2,n))*(ry(index1,n)-ry(index2,n))  
     *     +(rz(index1,n)-rz(index2,n))*(rz(index1,n)-rz(index2,n))  
        if(rrr.lt.cutoffELE*cutoffELE) then
           etime(n)=charge1*charge2/sqrt(rrr)
           etime(n)=etime(n)*332.0636
        else
           etime(n)=0.0 
        endif
        if(rrr.lt.cutoffVDW*cutoffVDW) then
           epsilon12=sqrt(epsilon1*epsilon2)
           sigma12=0.5*(sigma1+sigma2)
           vtime(n)=(sigma12**12/rrr**6-sigma12**6/rrr**3)*4.0*epsilon12
        else
           vtime(n)=0.0
        endif
        vdw=vdw+vtime(n)
        ele=ele+etime(n)
      enddo

      vdw=vdw/nstp
      ele=ele/nstp
      
      vrms=0.0
      erms=0.0
      do n = 1 , nstp
         vrms=vrms+(vtime(n)-vdw)*(vtime(n)-vdw)
         erms=erms+(etime(n)-ele)*(etime(n)-ele)
      enddo
      vrms=vrms/nstp
      erms=erms/nstp

      return
      end