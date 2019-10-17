C  CONSTRUCTS ORTHOGONAL MATRIX T FOR TRANSFORMATION FROM AO'S TO               
C  NATURAL HYBRID BOND ORBITALS, USING INPUT DENSITY MATRIX 'DM'. 
c  *** for aromatic rings              

      SUBROUTINE aromatic(DM,T,N,NATOMS,ipr)  
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)     
      integer CZ,U,UL
      character*4 symbol,resname,segname
      character*2 bname
      integer idres 
      common/iden/idres(nazr),symbol(nazr),resname(nazr),segname(nazr),
     :rx(nazr),ry(nazr),rz(nazr)
      common/LBL/bname(nbszr,2),label(NBSZR,3),IBX(NBSZR) 
      common/lblaro/labar(nbszr,7)
      common/INFO1/CZ(NAZR),U(NBSZR),UL(NAZR),LL(NAZR),
     :NELECS,NOCCA,NOCCB
      common/Q/Q(4,NBSZR) 
      common/POL/POL(NBSZR,3),IATHY(NBSZR,3),INO(NAZR)
      common/cring/idc6(6,20)
      dimension DM(NBSZR,NBSZR),T(NBSZR,NBSZR)
      dimension BLK(12,12),C(12,12),EVAL(12),BORB(12)
      character*2 NAME(4),LLP,LBD,L3C,LPB
      character*1 STAR,BLNK
      logical aring
      DATA NAME,STAR,BLNK/'LP','BD','3C','PB','*',' '/            
      DATA LLP,LBD,L3C,LPB/'LP','BD','3C','PB'/
      DATA OCCMX/1.9D0/

C  ZERO ARRAYS T, Q, INO, POL, IATHY          
      do 20 I=1,N            
      do 10 K=1,3            
      POL(I,K)=0.0           
      IATHY(I,K)=0           
   10 Q(K,I)=0.0D0           
      Q(4,I)=0.0D0           
      do 20 J=1,N            
   20 T(I,J)=0.0D0           
      do 30 I=1,NAZR
   30 INO(I)=0               
      IBD=0 
      NA1=NATOMS+1
      IAT3=0

      do 100 IB1=1,NA1       
      IB=IB1-1               
      do 90 IC1=2,NA1        
      IC=IC1-1               
      IF(IB.NE.0) GO TO 40   

C  LONE PAIRS            

      NCTR=1  
      IAT1=IC 
      IAT2=0  
      GO TO 70       

C  BOND PAIRS             

   40 NCTR=2  
      IAT1=IB 
      IAT2=IC 
      IF(IAT2.LE.IAT1) GO TO 90

   70 continue   
      IF(IBD.GE.NOCCA) GO TO 120 

      IF(IAT2.EQ.0)  GO TO  71
      if((rx(IAT1)-rx(IAT2))**2.0+(ry(IAT1)-ry(IAT2))**2.0+
     :   (rz(IAT1)-rz(IAT2))**2.0.gt.4.0) goto 90
   71 continue
      PRINT *, 'IAT1, IAT2', IAT1, IAT2

      call LOAD(DM,IAT1,IAT2,IAT3,BLK,NB)
      call JACVEC(NB,BLK,EVAL,C,12)

      do 80 IRNK=1,NB        

      call EXTRCT(C,EVAL,BORB,OCC,NB,IRNK) 
      if(occ.le.2.0d0.and.occ.gt.1.8d0) print *, 'occ', occ
      if(occ.gt.2.001d0) then
         write(6,'(f20.10)') occ
	 stop 'Occupany bigger than 2.0'         
      endif
      if(occ.le.1.8d0) go to 85
      
c      if(
c     :   (iat2.eq.0.and.symbol(iat1)(2:2).eq.'N').or.
c     :   (symbol(iat1)(2:2).eq.'C'.and.symbol(iat2)(2:2).eq.'O'
c     :    .and.dabs(borb(1)).lt.0.1d0).or.
c     :   (symbol(iat1)(2:2).eq.'O'.and.symbol(iat2)(2:2).eq.'C'
c     :    .and.dabs(borb(1)).lt.0.1d0)
c     :  ) then
c         borbmax=0.0
c         if(iat2.eq.0) korb=4
c         if(iat2.ne.0) korb=8
c         do iorb = 1 , korb
c            if(dabs(borb(iorb)).gt.borbmax) then
c               borbmax=dabs(borb(iorb))
c               ibmax=iorb 
c            endif
c         enddo
c         if(borb(ibmax).lt.0.d0) then
c            do irev = 1 , korb
c               borb(irev)=-borb(irev)
c            enddo
c         endif
c      else
c         if(borb(1).lt.0.d0) then
c            do irev = 1 , 8
c               borb(irev)=-borb(irev)
c            enddo
c         endif
c      endif
      
      IBD=IBD+1 
      if(nctr.eq.1) call DEPLET(DM,BORB,OCC,NB,IAT1,IAT2,IAT3)
      call STASH(BORB,OCC,IBD,IAT1,IAT2,IAT3) 
      bname(IBD,1)=NAME(NCTR)
      bname(IBD,2)=BLNK 
      label(IBD,1)=IRNK             
      label(IBD,2)=IAT1      
      label(IBD,3)=IAT2      

   80 continue               
   85 continue               
   90 continue               
  100 continue               
  120 continue

      aring=.false.
      nring=0
      do iat = 1 , natoms
         if(resname(iat).eq.'BENZ'.or.
     :      resname(iat).eq.'TYR '.or.
     :      resname(iat).eq.'PHE ') then
            aring=.true.
            if(symbol(iat).eq.' CG ') then 
               nring=nring+1
               idc6(1,nring)=iat
            endif
            if(symbol(iat).eq.' CD1') then 
               idc6(2,nring)=iat
            endif
            if(symbol(iat).eq.' CD2') then 
               idc6(3,nring)=iat
            endif
            if(symbol(iat).eq.' CE1') then 
               idc6(4,nring)=iat
            endif
            if(symbol(iat).eq.' CE2') then 
               idc6(5,nring)=iat
            endif
            if(symbol(iat).eq.' CZ ') then 
               idc6(6,nring)=iat
            endif
         endif
         if(resname(iat).eq.'TRP ') then
            aring=.true.
            if(symbol(iat).eq.' CD2') then 
               nring=nring+1
               idc6(1,nring)=iat
            endif
            if(symbol(iat).eq.' CE2') then 
               idc6(2,nring)=iat
            endif
            if(symbol(iat).eq.' CE3') then 
               idc6(3,nring)=iat
            endif
            if(symbol(iat).eq.' CZ2') then 
               idc6(4,nring)=iat
            endif
            if(symbol(iat).eq.' CZ3') then 
               idc6(5,nring)=iat
            endif
            if(symbol(iat).eq.' CH2') then 
               idc6(6,nring)=iat
            endif
         endif
      enddo
      if(aring) then      
         write(6,*)
         write(6,*) 'Handling aromatic rings ...'
         write(6,*)
      endif
      print *, 'The biomolecule contains',nring,' aromatic rings'
      if(aring) then
         do iring = 1 , nring
            call fourth(ibd,iring,natoms,ipr)
            indextemp=0
            do i = ibd-2,ibd
               indextemp=indextemp+1
               bname(i,1)=name(4)
               bname(i,2)=blnk
               labar(i,1)=iring
               labar(i,2)=idc6(1,iring)
               labar(i,3)=idc6(2,iring)
               labar(i,4)=idc6(3,iring)
               labar(i,5)=idc6(4,iring)
               labar(i,6)=idc6(5,iring)
               labar(i,7)=idc6(6,iring)
            enddo
         enddo
      endif
      
      IBO=IBD 
      do 150 I=1,IBO         
      IF(bname(I,1).EQ.NAME(1)) GO TO 150     
      NAB=1 
      IF(bname(I,1).EQ.NAME(3)) NAB=2         
      do 140 IAB=1,NAB       
      IBD=IBD+1              
      do 130 J=1,3
  130 label(IBD,J)=label(I,J)
      bname(IBD,1)=bname(i,1)
      bname(IBD,2)=STAR
      if(aring) then
         do j = 1 , 7
           labar(ibd,j)=labar(i,j)
         enddo
      endif
  140 continue               
  150 continue               
      IF(IBD.EQ.N) GO TO 160 
C  MISCOUNTED BOND ORBITALS...EXIT            
      WRITE(6,1000) IBD,N    
      call abortBondo('aromat')
  160 continue
  
      call REFORM(nring,T,NBSZR,N,NATOMS,ipr)
      if(aring) then
         call REFORMaro(nring,T,NBSZR,N,NATOMS,ipr)
      endif

      return  
 1000 FORMAT(2X,'*** FOUND ',I4,' ORBITALS, BUT DIM. SHOULD BE ',I4,'*')        
      END