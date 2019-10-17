      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      include 'PARAM.include'
      INTEGER  OPTION,OPNCLO,HUCKEL,CNDO,INDO,CLOSED,OPEN                       
      INTEGER ORB,EL,AN,CHARGE,CZ,U,ULIM,OCCA,OCCB,TOPO 
      LOGICAL WANTC                                                             
      COMMON/ARRAYS/TAB(NBSZR,NBSZR),BC(2*NBSZRS)
      COMMON/INFO/NATOMS,CHARGE,MULTIP,AN(NAZR),C(NAZR,3),N,NTA,INO,IPR
      COMMON/PERTBL/EL(18)                                                      
      COMMON/ORBIT/ORB(9)                                                       
      COMMON/GAB/XYZ(NAIGAIO2)
      COMMON/INFO1/CZ(NAZR),U(NBSZR),ULIM(NAZR),LLIM(NAZR),
     $NELECS,OCCA,OCCB 
      COMMON/OPSION/OPTION,OPNCLO,HUCKEL,CNDO,INDO,CLOSED,OPEN                  
      COMMON/AUXINT/A(17),B(17)                                                 
      COMMON/EXTRA/T(NBSZR,NBSZR),IHYB(NAZR),TOPO(NAZR,NAZR),
     $T2J(NBSZR)
      COMMON/NHYBRD/IDEN(NBSZR,3)
      COMMON/WANT/IWNAT,IWLCAO,NOTOPO                                           
      COMMON/HIT/NDEL,LIST(40)                                                  
      DIMENSION WORD(2) 

      INTEGER  cn(NAZR), atomnum(NAZR), atonn(NAZR,4)
      CHARACTER*4 atomtype(NAZR)
      CHARACTER bondt(NAZR,4)

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
      IF(OPNCLO.EQ.OPEN) IWNAT=0                                                
      WRITE(6,2000) NATOMS,CHARGE,MULTIP,NTA,IPR,IWNAT,IWLCAO,NOTOPO            
      IF(IWLCAO.GT.0) WRITE(6,1000) OPTION                                      
      IF(IWLCAO.EQ.0) WRITE(6,1100) OPTION,WORD(IWNAT+1)                        
      IF(NDEL.GT.0) READ(7,1200) (LIST(I),I=1,NDEL)                             
      WRITE (6,1300) NDEL                                                       
      IF(NDEL.GT.0) WRITE(6,1400) (LIST(I),I=1,NDEL)                            

      close(unit=7)

      IF(NOTOPO.GT.0) GO TO 20 


c       ZEROING (INITIALIZING) TOPO MATRIX
        Do i=1,NTA
         Do j=1,NTA
         TOPO(i,j)=0
         End Do
        End Do
c
c       READING OF PROTTN.HIN   CONTAINING
c       CONNECTIVITY INFO OF REAL+GHOST ATOMD
cx         open(file='prottn.hin', unit=5)
cx        Do i=1,NTA
cx        READ (5,*) atomnum(i), atomtype(i), cn(i),
cx     $         (atonn(i,j), bondt(i,j), j=1,cn(i))
cx       Zeroing remaining array elements
cx        Do k=cn(i)+1,4
cx       atonn(i,k)=0
cx        bondt(i,j)='n'
cx        End Do
cx        End Do
cx        close(unit=5)
c
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

      CALL SCF                                                                  

      STOP

 1000 FORMAT(/38X,'>>> STANDARD ',A1,'NDO LCAO-MO CALCULATION <<<'/)            
 1100 FORMAT(/35X,'>>> ',A1,'NDO LCBO-MO METHOD ',A5,'AL HYBRIDS <<<'/)        
 1200 FORMAT(20I4)                                                              
 1300 FORMAT(I4,' B.O.-BASIS FUNCTIONS TO BE DELETED')                         
 1400 FORMAT(10X,30I4)                                                          
 1500 FORMAT(20A3)                                                              
 1600 FORMAT(5X,20A3)                                                           
Csexp2 1700 FORMAT(A3,1X,A3,1X,I1,4I5) 
 1700 FORMAT(A3,1X,A3)  
Csexp2                                             
 1800 FORMAT(5X,A4,1X,A4,1X,A4)                                                 
 1900 FORMAT(18I4)                                                              
 2000 FORMAT(3X,I4,' ATOMS, CHARGE =',I2,', MULT. =',I2,', TOTAL # ',           
     + 'ATOMS =',I3,', PRINT LEVEL =',I2,', IWNAT =',I2,', IWLCAO =',           
     + I2,', NOTOPO =',I2)                                                      
 2100 FORMAT(35I1)                                                              
 2200 FORMAT(' ',35I1)                                                          
 2300 FORMAT(36X,46('*'))                                                       
 2400 FORMAT(36X,'*',44X,'*')                                                   
 2500 FORMAT(36X,'*',2X,'REACTION COORDINATE VAR(',I3,') =',F10.4,2X,           
     + '*')                                                                     
 2600 FORMAT(2X,'GEOMETRY-OPTIMIZED RESULT'/)                                  
 2700 FORMAT(4I5)                                                               
Csexp2 2800 FORMAT(12X,'(POWELL RESTART PARAMETERS IAG =',I2,', NTLIM =',I5,         
Csexp2     + ', NTITER =',I5,', LFN =',I5,', OPT. LEVEL =',I2,')')  
      END                                                                       
