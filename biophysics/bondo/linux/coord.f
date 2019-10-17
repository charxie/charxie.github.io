      SUBROUTINE COORD
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      integer AN,CH,TOPO
      character*4 symbol,resname,segname
      integer idres 
      common/iden/idres(nazr),symbol(nazr),resname(nazr),segname(nazr),
     :rx(nazr),ry(nazr),rz(nazr)
      common/INFO/NATOMS,CH,MULT,AN(NAZR),Z(NAZR),Y(NAZR),X(NAZR),
     :N,NTA,INO,IPR
      common/EXTRA/T(NBSZR,NBSZR),IHYB(NAZR),TOPO(NAZR,NAZR),
     :T2J(NBSZR)
      common/WANT/IWNAT,IWLCAO,NOTOPO
      DATA ISW/0/                                                               

      PRINT *, IWNAT
      PRINT *, NATOMS
      PRINT *, NTA

Cs READING REAL and GHOST PDB FILES AND ASSIGNMENT OF
Cs ATOMIC NUMBERS
Cs AN(NTA) is integer,
Cs IF IWNAT=1 READ ONLY REAL PDB,
Cs IF IWNAT=0 READ BOTH REAL AND GHOST PDB FILES

       open(file='proton.pdb', unit=15)
       read(15,*)
       read(15,*)
       do i=1,NATOMS
          READ(15,140)symbol(i),resname(i),idres(i),x(i),y(i),z(i),
     :                segname(i)
          IF (symbol(i)(2:2) .EQ. 'N') THEN
             AN(i)=7
          ELSEIF(symbol(i)(2:2).EQ.'H'.or.symbol(i)(1:1).EQ.'H')THEN
             AN(i)=1
          ELSEIF(symbol(i)(2:2).EQ.'C'.or.symbol(i)(1:1).EQ.'C')THEN
             AN(i)=6
          ELSEIF(symbol(i)(2:2).EQ.'O'.or.symbol(i)(1:1).EQ.'O')THEN
             AN(i)=8
          ELSEIF(symbol(i)(2:2).EQ.'S'.or.symbol(i)(1:1).EQ.'S')THEN
             AN(i)=16
          ELSEIF(symbol(i)(2:2).EQ.'B')THEN
             AN(i)=5
          ELSE
             PRINT *, 'UNRECOGNIZED ATOM'
             stop
          END IF
       enddo
  140  format(12x,a4,1x,a4,i5,4x,3f8.3,18x,a4)
       close(15)

       IF(IWNAT.NE.1) THEN
          PRINT*, 'OPA'
          open(file='ghost.pdb', unit=15)
          Do i=NATOMS+1, NTA
             READ(15,140) symbol(i),resname(i),idres(i),
     :                    x(i),y(i),z(i),segname(i)
             IF(symbol(i)(1:1) .EQ. 'G') THEN
                AN(i)=0
             ELSE
                PRINT *, 'UNRECOGNIZED ATOM'
                stop
             ENDIF
          EndDo
          close (15)
       ENDIF

      IF(IPR.LT.3) return                                                       
      WRITE(6,1000)                                                             
      DO I=1,NTA                                                            
         WRITE(6,1100)symbol(i)(2:2),AN(I),X(I),Y(I),Z(I),
     :                resname(i),idres(i)
      ENDDO

      do i = 1 , natoms
         rx(i)=x(i)
         ry(i)=y(i)
         rz(i)=z(i)
      enddo
      
      return                                                                    

 1000 format(/1X,'ATOM',2X,'ATOMIC NO. ',2X,'X COORDINATE',3X,
     : 'Y COORDINATE',3X,'Z COORDINATE',2X,'RESIDUE NO.'/)                                                        
 1100 format(2X,A2,6X,I2,6X,F12.7,2(3X,F12.7),5X,A4,I4)

      END