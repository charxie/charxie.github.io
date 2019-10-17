c*****************************************************************************
c Modified Version of BONDO, on Linux
c
c Biophysics Group, University of Cyprus, 2000
c
c Changelog:
c
c * 3-center bond disabled
c * aromatic ring enabled
c * bond tags set
c
c*****************************************************************************
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)     
      include 'PARAM.include'
      character option*4,opnclo*3
      character*3 nden
      integer HUCKEL,CNDO,INDO,CLOSED,OPEN
      integer ORB,EL,AN,CHARGE,CZ,U,ULIM,OCCA,OCCB,TOPO 
      logical WANTC
      common/ARRAYS/TAB(NBSZR,NBSZR),BC(2*NBSZRS)
      common/INFO/NATOMS,CHARGE,MULTIP,AN(NAZR),C(NAZR,3),N,NTA,INO,IPR
      common/PERTBL/EL(18)
      common/ORBIT/ORB(9)
      common/GAB/XYZ(NAIGAIO2)
      common/INFO1/CZ(NAZR),U(NBSZR),ULIM(NAZR),LLIM(NAZR),
     :NELECS,OCCA,OCCB 
      common/OPSION/OPTION,OPNCLO,HUCKEL,CNDO,INDO,CLOSED,OPEN 
      common/AUXINT/A(17),B(17) 
      common/EXTRA/T(NBSZR,NBSZR),IHYB(NAZR),TOPO(NAZR,NAZR),T2J(NBSZR)
      common/NHYBRD/nden(nbszr),IDEN(NBSZR,2)
      common/WANT/IWNAT,IWLCAO,NOTOPO         
      common/HIT/NDEL,LIST(40)

      character*5 WORD(2) 
      integer cn(NAZR), atomnum(NAZR), atonn(NAZR,4)
      character*4 atomtype(NAZR)
      character bondt(NAZR,4)

      DATA ISW/0/
      DATA WORD/'NOMIN','NATUR'/  
C 
C  INPUT BEGINS (READ FROM LFN 7) 
C 
      open (file='input', unit=7)
      READ(7,1500) (AN(I),I=1,20) 
      WRITE(6,1600) (AN(I),I=1,20)  
      READ(7,1700) OPTION,OPNCLO  
      WRITE(6,1800) OPTION,OPNCLO 
      READ(7,1900) NATOMS,CHARGE,MULTIP,NTA,NDEL,IPR,IWNAT,IWLCAO,NOTOPO        
      IF(OPNCLO.EQ.'OPN') IWNAT=0  
      WRITE(6,2000) NATOMS,CHARGE,MULTIP,NTA,IPR,IWNAT,IWLCAO,NOTOPO  
      IF(IWLCAO.GT.0) WRITE(6,1000) OPTION    
      IF(IWLCAO.EQ.0) WRITE(6,1100) OPTION,WORD(IWNAT+1)       
      IF(NDEL.GT.0) READ(7,1200) (LIST(I),I=1,NDEL)  
      WRITE (6,1300) NDEL    
      IF(NDEL.GT.0) WRITE(6,1400) (LIST(I),I=1,NDEL) 

      close(unit=7)

      IF(NOTOPO.GT.0) GO TO 20 

      do j=1,NTA
        do i=1,NTA
          TOPO(i,j)=0
        enddo
      enddo

c READING OF prottn.hin CONTAINING CONNECTIVITY INFO OF REAL+GHOST ATOMD

cx        open(file='prottn.hin', unit=5)
cx        Do i=1,NTA
cx        READ (5,*) atomnum(i), atomtype(i), cn(i),
cx     $         (atonn(i,j), bondt(i,j), j=1,cn(i))
cx        Do k=cn(i)+1,4
cx        atonn(i,k)=0
cx        bondt(i,j)='n'
cx        End Do
cx        End Do
cx        close(unit=5)

       Do i=1,NTA
       Do j=1,cn(i)
c      if single bond topo(i,j)=1
       IF (bondt(i,j) .EQ. 's') THEN
       TOPO(atomnum(i),atonn(i,j))=1
c      if single bond topo(i,j)=2
       ELSEIF (bondt(i,j) .EQ. 'd') THEN
       TOPO(atomnum(i),atonn(i,j))=2
       ELSE
c      if aromatic bond complain
       PRINT *, 'BONDTYPE NOT ALLOWED', bondt(i,j)
       END IF
       End Do
       End Do

   20 IPR0=IPR 

      call SCF 

      STOP

 1000 format(/38X,'>>> STANDARD ',A1,'NDO LCAO-MO CALCULATION <<<'/)  
 1100 format(/35X,'>>> ',A1,'NDO LCBO-MO METHOD ',A5,'AL HYBRIDS <<<'/)        
 1200 format(20I4) 
 1300 format(I4,' B.O.-BASIS FUNCTIONS TO BE DELETED')        
 1400 format(10X,30I4)       
 1500 format(20A3) 
 1600 format(5X,20A3)        
Csexp2 1700 format(A3,1X,A3,1X,I1,4I5) 
 1700 format(A4,1X,A3)  
Csexp2 
 1800 format(5X,A4,1X,A4,1X,A4) 
 1900 format(18I4) 
 2000 format(3X,I4,' ATOMS, CHARGE =',I2,', MULT. =',I2,', TOTAL # ', 
     + 'ATOMS =',I3,', PRINT LEVEL =',I2,', IWNAT =',I2,', IWLCAO =', 
     + I2,', NOTOPO =',I2)   
 2100 format(35I1) 
 2200 format(' ',35I1)       
 2300 format(36X,46('*'))    
 2400 format(36X,'*',44X,'*')   
 2500 format(36X,'*',2X,'REACTION COORDINATE VAR(',I3,') =',F10.4,2X, 
     + '*') 
 2600 format(2X,'GEOMETRY-OPTIMIZED RESULT'/)   
 2700 format(4I5)  
Csexp2 2800 format(12X,'(POWELL RESTART PARAMETERS IAG =',I2,', NTLIM =',I5,         
Csexp2     + ', NTITER =',I5,', LFN =',I5,', OPT. LEVEL =',I2,')')  
      END