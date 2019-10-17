C     CNDO-INDO SCF CLOSED SHELL- PRINTOUT SEGMENT                              

      SUBROUTINE CPRINT 
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      INTEGER  OPTION,OPNCLO,HUCKEL,CNDO,INDO,CLOSED,OPEN                       
      INTEGER CHARGE,AN,U,ULIM,EL,OCCA,OCCB,CZ,ANI                              
      COMMON/ARRAYS/A(NBSZR,NBSZR),B(NBSZR,NBSZR),D(NBSZR,NBSZR) 
      COMMON/GAB/XXX(5*NBSZR),G(NAZR,NAZR),
     $Q(NBSZR),YYY(NBSZR),ENERGY,XXY(214)
      COMMON/INFO/NATOMS,CHARGE,MULTIP,AN(NAZR),C(NAZR,3),N,NTA,INO,IPR 
      COMMON/INFO1/CZ(NAZR),U(NBSZR),ULIM(NAZR),LLIM(NAZR),
     $NELECS,OCCA,OCCB
      COMMON/PERTBL/EL(18)                                                      
      COMMON/OPSION/OPTION,OPNCLO,HUCKEL,CNDO,INDO,CLOSED,OPEN                  
      COMMON/ETOT/ETOT                                                          
      DIMENSION ATENG(18)                                                       
      IF (OPTION.EQ.CNDO) GO TO 10                                              
      ATENG(1)=-0.6387302462  D0                                                
      ATENG(3)=-.2321972405   D0                                                
      ATENG(4)=-1.1219620354  D0                                                
      ATENG(5)=-2.8725750048  D0                                                
      ATENG(6)=-5.9349548261  D0                                                
      ATENG(7)=-10.6731741251 D0                                                
      ATENG(8)=-17.2920850650 D0                                                
      ATENG(9)=-26.2574377875 D0                                                
      GO TO 20                                                                  
   10 CONTINUE                                                                  
      ATENG(1)=-0.6387302462  D0                                                
      ATENG(3)=-.2321972405   D0                                                
      ATENG(4)=-1.1454120355  D0                                                
      ATENG(5)=-2.9774239048  D0                                                
      ATENG(6)=-6.1649936261  D0                                                
      ATENG(7)=-11.0768746252 D0                                                
      ATENG(8)=-18.0819658651 D0                                                
      ATENG(9)=-27.5491302880 D0                                                
      ATENG(11)=-.1977009568  D0                                                
      ATENG(12)=-.8671913833  D0                                                
      ATENG(13)=-2.0364557744 D0                                                
      ATENG(14)=-3.8979034686 D0                                                
      ATENG(15)=-6.7966009163 D0                                                
      ATENG(16)=-10.7658174341D0                                                
      ATENG(17)=-16.0467017940D0                                                
   20 CONTINUE                                                                  
      K=NATOMS-1                                                                
      IF(IPR.LT.6) GO TO 30                                                     
      WRITE(6,1000)                                                             
      CALL SCFOUT(0,2,1)                                                        
   30 CONTINUE                                                                  
      VNN=0.0D0                                                                 
      DO 40 I=1,K                                                               
      L=I+1                                                                     
      DO 40 J=L,NATOMS                                                          
      RAD=DSQRT((C(I,1)-C(J,1))**2+(C(I,2)-C(J,2))**2+(C(I,3)-        
     +C(J,3))**2)                                                               
   40 VNN=VNN+(DFLOAT(CZ(I))*DFLOAT(CZ(J)))/RAD                                 
      IF(IPR.GE.1) WRITE(6,1200) VNN                                            
      ENERGY=ENERGY+VNN                                                         
      ETOT=ENERGY                                                               
      IF(IPR.LT.1) RETURN                                                       
      WRITE(6,1300) ENERGY                                                      
      DO 50 I=1,NATOMS                                                          
      ANI=AN(I)                                                                 
   50 ENERGY=ENERGY-ATENG(ANI)                                                  
      WRITE(6,1400) ENERGY                                                      
      RETURN                                                                    
 1000 FORMAT(' DENSITY MATRIX IN DEORTHOGONALIZED B.O. BASIS')                  
 1200 FORMAT(11X,'NUCL. REPULSION ENERGY',F20.10)                                
 1300 FORMAT(11X,'TOTAL ENERGY',10X,F20.10)                                      
 1400 FORMAT(11X,'BINDING ENERGY',8X,F20.10,' A.U.'/)                            
      END                                                                       
