      SUBROUTINE EIGOUT(M,K)
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
C                                                                               
C  THIS ROUTINE IS CALLED IN SCFOUT TO PRINT THE EIGENVALUES M TO K             
C                                                                               
      COMMON/GAB/XXX(3*NBSZR),EPSILN(NBSZR),
     $YYY(NAIGAIO2-(4*NBSZR))
      common/homo/lumo

      WRITE(6,1000) (EPSILN(I),I=M,K)
      do i = m , k
         if(i.eq.lumo-1) then
            write(6,1001) epsiln(i)
         endif
         if(i.eq.lumo) then
            write(6,1002) epsiln(i)
         endif
      enddo 

      return                                                                    
 1000 FORMAT(/,15H EIGENVALUES---,20(F9.4),/)
 1001 format(/,15H *** HOMO =    ,F9.4)
 1002 format(  15H *** LUMO =    ,F9.4,/)
      END                                                                      
