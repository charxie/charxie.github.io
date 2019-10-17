      SUBROUTINE SCF
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER OPTION,OPNCLO,HUCKEL,CNDO,INDO,CLOSED,OPEN                        
      COMMON/INFO/NATOMS,ICH,MULT,IAN(NAZR),C(NAZR,3),N,NTA,INO,IPR
      COMMON/OPSION/OPTION,OPNCLO,HUCKEL,CNDO,INDO,CLOSED,OPEN                  
      COMMON/WANT/IWNAT,IWLCAO,NOTOPO
      common/sulfur/isulfur

      PRINT *, 'IWNAT=', IWNAT
      
      open(99,file='flag_of_bondo')
      write(99,*) 0
      close(99)

      CALL COORD  

      PRINT *, 'IWNAT=', IWNAT

      IF(IPR.LT.3) GO TO 10                                                     
      LFN=6                                                                     
   10 CONTINUE                                                                  
C                                                                               
C  CONVERT ATOMIC COORDINATES TO ATOMIC UNITS                                   
C                                                                               
      DO 20 I=1,NTA                                                             
      DO 20 J=1,3                                                               
   20 C(I,J)=C(I,J)/0.529167D0                                                  
      IWNAT0=IWNAT                                                              
      IWNAT=0 
      if((IWLCAO.EQ.0).AND.(NOTOPO.EQ.0)) then
         CALL BASIS
      endif
      IWNAT=IWNAT0 
      CALL COEFFT                                                               
      CALL INTGRL                                                               
      CALL HUCKCL
      if(isulfur.eq.0) CALL SCFCLO
      if(isulfur.eq.1) then
c         iwlcao=1
         CALL SCFCLO
      endif
      CALL CPRINT                                                               

      open(99,file='flag_of_bondo')
      write(99,*) 1
      close(99)

      RETURN                                                                    
      END                                                                       
