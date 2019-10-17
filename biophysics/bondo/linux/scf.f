      SUBROUTINE SCF
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character option*4,opnclo*3
      integer HUCKEL,CNDO,INDO,CLOSED,OPEN                        
      common/INFO/NATOMS,ICH,MULT,IAN(NAZR),C(NAZR,3),N,NTA,INO,IPR
      common/OPSION/OPTION,OPNCLO,HUCKEL,CNDO,INDO,CLOSED,OPEN                  
      common/WANT/IWNAT,IWLCAO,NOTOPO

      PRINT *, 'IWNAT=', IWNAT

c write a flag file for the external c-shell script file      
      open(99,file='flag_of_bondo')
      write(99,*) 0
      close(99)

c read the PDB file of the molecular structure      
      call COORD  

      PRINT *, 'IWNAT=', IWNAT

      IF(IPR.LT.3) GO TO 10                                                     
      LFN=6                                                                     
   10 continue                                                                  
C                                                                               
C  CONVERT ATOMIC COORDINATES TO ATOMIC UNITS                                   
C                                                                               
      DO 20 I=1,NTA
      DO 20 J=1,3
   20 C(I,J)=C(I,J)/0.529167D0
      IWNAT0=IWNAT
      IWNAT=0
      if((IWLCAO.EQ.0).AND.(NOTOPO.EQ.0)) then
         call BASIS
      endif
      IWNAT=IWNAT0
      call COEFFT
      call INTGRL
      call HUCKCL
      call SCFCLO
      call CPRINT

c when finishing the job change the flag
      open(99,file='flag_of_bondo')
      write(99,*) 1
      close(99)

      return                                                  
      END