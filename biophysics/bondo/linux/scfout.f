C  OP = 1, THE EIGENVALUES (IN common/GAB/) ARE PRINTED.                        
C  IF MOP = 1 OR 3, ONLY LOWER TRIANGLE IS PRINTED (SYMMETRIZED)                
C                                                                               
      SUBROUTINE SCFOUT(OP,MOP,IHY)   
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      integer OP,AN,ANII,CZ,U,ORB,ULIM,EL,CHARGE,OCCA,OCCB
      character l3c*2
      character*2 bname
      character*3 nden
      common/GAB/XXX(NAIGAIO2)
      common/INFO/NATOMS,CHARGE,MULTIP,AN(NAZR),C(NAZR,3),NN,NTA,INO,IPR
      common/INFO1/CZ(NAZR),U(NBSZR),ULIM(NAZR),LLIM(NAZR),
     $NELECS,OCCA,OCCB
      common/ORBIT/ORB(9)                                                       
      common/PERTBL/EL(18)                                                      
      common/NHYBRD/nden(nbszr),iden(NBSZR,2)
      common/WANT/IWNAT,IWLCAO,NOTOPO                                           
      common/LBL/bname(nbszr,2),label(NBSZR,3),IBX(NBSZR) 
      dimension A(20)

      N=NN
      IF(IHY.EQ.1)N=INO                                                         
      DO 130 M=1,N,11                                                           
      K=M+10                                                                    
      IF (K.LE.N) GO TO 20                                                      
   10 K=N                                                                       
   20 continue                                                                  
      WRITE(6,1400)                                                             
      IF (OP.EQ.1) GO TO 30                                                     
      GO TO 40                                                                  
   30 call EIGOUT(M,K)                                                          
   40 continue                                                                  
      WRITE(6,1200) (I,I=M,K)                                                   
      DO 120 I=1,N                                                              
      KM=0                                                                      
      DO 45 J=M,K                                                               
      KM=KM+1                                                                   
   45 A(KM)=ALT(I,J,MOP)                                                        
      IF(IHY.EQ.1)GO TO 70                                                      
C  ATOMIC ORBITAL BASIS, IHY=0                                                  
      II=U(I)                                                                   
      ANII=AN(II)                                                               
      L=I-LLIM(II)+1                                                            
   50 WRITE(6,1300) I,II,EL(ANII),ORB(L),(A(J),J=1,KM)                          
      IF (I.EQ.ULIM(II)) GO TO 60                                               
      GO TO 120                                                                 
   60 WRITE(6,1400)                                                             
      GO TO 120                                                                 
C  BOND ORBITAL BASIS, IHY=1                                                    
   70 continue                                                                  
      IF(IWNAT.GT.0) GO TO 80                                                   
      WRITE(6,1500)I,nden(I),iden(I,1),iden(I,2),(A(J),J=1,KM)                
      GO TO 120
C  NATURAL HYBRID B.O. BASIS, IWNAT=1                                           
   80 continue                                                                  
      DO 90 IB=1,N                                                              
      IF(IBX(IB).EQ.I) GO TO 100
   90 continue
  100 continue
      L3C='3C'                                                                  
      IF(bname(IB,1).EQ.L3C) GO TO 110                                          
      WRITE(6,1000) I,(bname(ib,j),j=1,2),(label(IB,J),J=1,3),
     :              (A(J),J=1,KM)
      GO TO 120                                                                 
  110 WRITE(6,1100) I,(label(IB,J),J=1,3),bname(IB,2),(A(J),J=1,KM)             
  120 continue                                                                  
  130 continue                                                                  
      WRITE(6,1400)                                                             
      WRITE(6,1400)                                                             

      return                                                                    
 1000 FORMAT(1X,I3,'. ',A2,A1,I1,'0',I3,'-',I3,20(F9.4))                        
 1100 FORMAT(1X,I3,'. ',I2,'-',I2,'-',I2,A1,'0',20(F9.4))                       
 1200 FORMAT(13X,20I9)                                                          
 1300 FORMAT(1X,I3,I3,A4,1X,A4,20(F9.4))                                        
 1400 FORMAT(1X)                                                                
 1500 FORMAT(1X,I3,2X,A4,1X,I2,'-',I2,20(F9.4))                                 
      END                                                                       
