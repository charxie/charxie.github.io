C  CNDO/INDO CLOSED SHELL SCF SEGMENT                                           
C  GAMMA MATRIX CONTAINED IN G, CORE HAMILTONIAN CONTAINED IN Q AND             
C  UPPER TRIANGLE OF A, AND INITIAL DENSITY MATRIX CONTAINED IN B               
C  OPTIONS   CNDO OR INDO                                                       
c  01/05/98: delete the back-transformation from BO-AO when Z=26
c            to save computational time; 
c
      SUBROUTINE SCFCLO1
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)                                      
      INTEGER  OPTION,OPNCLO,HUCKEL,CNDO,INDO,CLOSED,OPEN                       
      INTEGER CHARGE,OCCA,OCCB,ULIM,U,AN,CZ,Z,TOPO                              
      character*4 symbol,resname,segname
      integer idres 
      common/iden/idres(nazr),symbol(nazr),resname(nazr),segname(nazr),
     :rx(nazr),ry(nazr),rz(nazr)
      COMMON/ARRAYS/A(NBSZR,NBSZR),B(NBSZR,NBSZR),D(NBSZR,NBSZR)
      COMMON/GAB/XXX(5*NBSZR),G(NAZR,NAZR),
     $Q(NBSZR),YYY(NBSZR),ENERGY,XXY(214) 
      COMMON/INFO/NATOMS,CHARGE,MULTIP,AN(NAZR),C(NAZR,3),N,NTA,NO,IPR
      COMMON/INFO1/CZ(NAZR),U(NBSZR),ULIM(NAZR),LLIM(NAZR),
     $NELECS,OCCA,OCCB
      COMMON/OPSION/OPTION,OPNCLO,HUCKEL,CNDO,INDO,CLOSED,OPEN                  
      COMMON/EXTRA/T(NBSZR,NBSZR),IHYB(NAZR),TOPO(NAZR,NAZR),
     $T2J(NBSZR)
      COMMON/WANT/IWNAT,IWLCAO,NOTOPO
      COMMON/LBL/LABEL(NBSZR,6),IBX(NBSZR) 
      dimension fk(nbszr,nbszr),fkao(nbszr,nbszr)

      DIMENSION G1(18),F2(18)                                                   
      COMMON/HIT/NLIST,LIST(40)                                                 
      DIMENSION WORDS(3,2)                                                      
      DATA WORDS/'ATOMIC',' ORBIT','AL','BOND-A','NTIBON','D '/                 
      G1(3)=.092012D0                                                           
      G1(4)=.1407  D0                                                           
      G1(5)=.199265D0                                                           
      G1(6)=.267708D0                                                           
      G1(7)=.346029D0                                                           
      G1(8)=.43423 D0                                                           
      G1(9)=.532305D0                                                           
      F2(3)=.049865D0                                                           
      F2(4)=.089125D0                                                           
      F2(5)=.13041 D0                                                           
      F2(6)=.17372 D0                                                           
      F2(7)=.219055D0                                                           
      F2(8)=.266415D0                                                           
      F2(9)=.31580 D0                                                           
      Z=0                                                                       
      IT=25                                                                     
      RHO=1.D-6                                                                 
      IWBO=1                                                                    
      IF(IWLCAO.GT.0) IWBO=0                                                    
      OLDENG=0.0D0

      PRINT*, 'N=', N
      
C                                                                               
C  GRAND SCF-ITERATIONS LOOP, Z.LT.25                                          
C     Z.LE.25 SCF ITERATIONS BEFORE CONVERGENCE                                
C     Z.EQ.26 B.O. TRANSFORMATION, BASIS-SET TRUNCATION                        
C     Z.EQ.27 DEORTHOGONALIZATION AND EXIT                                     
C                                                                               
   10 CONTINUE                                                                  
                                                              
      IF(Z.NE.25) GO TO 20                                                      
      WRITE(6,1600)                                                             
      CALL EXIT                                                                 
   20 CONTINUE                                                                  
      Z = Z+1
      if(z.eq.27) goto 170
      ENERGY = 0.D0                                                             
C                                                                               
C  CONSTRUCT FOCK MATRIX (IN A) FROM DENSITY MATRIX (IN B)0                     
C  ...TRANSFER CORE HAMILTONIAN TO LOWER TRIANGLE OF A...                       
C                                                                               
      DO 30 I=1,N                                                               
      A(I,I)=Q(I)                                                               
      DO 30 J=I,N                                                               
   30 A(J,I)=A(I,J)                                                             
      DO 40 I=1,N                                                               
      II=U(I)                                                                   
      A(I,I)=A(I,I)-B(I,I)*G(II,II)*0.5D0                                       
      DO 40 K=1,N                                                               
      JJ=U(K)                                                                   
   40 A(I,I)=A(I,I)+B(K,K)*G(II,JJ)                                             
      NM=N-1                                                                    
      DO 50 I=1,NM                                                              
      II=U(I)                                                                   
      LL=I+1                                                                    
      DO 50 J=LL,N                                                              
      JJ=U(J)                                                                   
   50 A(J,I)=A(J,I)-B(J,I)*G(II,JJ)*0.5D0                                       
C  INDO MODIFICATION                                                            
      IF (OPTION.EQ.CNDO) GO TO 100                                             
   60 DO 90 II=1,NATOMS                                                         
      K=AN(II)                                                                  
      I=LLIM(II)                                                                
      IF (K.EQ.1) GO TO 90                                                      
   70 PAA=B(I,I)+B(I+1,I+1)+B(I+2,I+2)+B(I+3,I+3)                               
      A(I,I)=A(I,I)-(PAA-B(I,I))*G1(K)/6.D0                                    
      DO 80 J=1,3                                                               
      A(I+J,I+J)=A(I+J,I+J)-B(I,I)*G1(K)/6.D0-(PAA-B(I,I))*7.D0*F2(K)/5        
     :      0.D0+B(I+J,I+J)*11.D0*F2(K)/50.D0                                   
   80 A(I+J,I)=A(I+J,I)+B(I,I+J)*G1(K)/2.D0                                    
      I1=I+1                                                                    
      I2=I+2                                                                    
      I3=I+3                                                                    
      A(I2,I1)=A(I2,I1)+B(I2,I1)*11.D0*F2(K)/50.D0                              
      A(I3,I1)=A(I3,I1)+B(I3,I1)*11.D0*F2(K)/50.D0                              
      A(I3,I2)=A(I3,I2)+B(I3,I2)*11.D0*F2(K)/50.D0                              
   90 CONTINUE
C                                                                               
C  FOCK MATRIX COMPLETE; EVALUATE AND PRINT OUT ENERGY
C                                                                               
  100 CONTINUE                                                                  
      DO 110 I=1,N                                                              
  110 ENERGY=ENERGY+0.5D0*B(I,I)*(A(I,I)+Q(I))                                  
      DO 120 I=1,NM                                                             
      LL=I+1                                                                    
      DO 120 J=LL,N                                                             
  120 ENERGY=ENERGY+B(I,J)*(A(I,J)+A(J,I))                                      

      IF((NLIST.EQ.0).AND.(Z.EQ.27))GO TO 170                                   
      IF((Z.EQ.27).AND.(IPR.GE.1)) WRITE(6,1200)                                
      IF(IPR.GE.4) WRITE(6,1300) ENERGY                                         
      IF(DABS(ENERGY-OLDENG).GE.0.000001D0) GO TO 170 
C                                                                               
C  ENERGY SATISFIED; SET Z=26 FOR FINAL TWO PASSES
C                                                                               
  130 IF(Z.LE.IT)Z=26                                                           
  140 IF((Z.LT.27).AND.(IPR.GE.1)) WRITE(6,1400)                                
      IF(IPR.GE.1) WRITE(6,1300) ENERGY                                         
      IF(Z.EQ.27)GO TO 170                                                      
C                                                                               
C  Z = 26  TRANSFORM, DELETE, DIAGONALIZE, PRINT B.O. EIGENVECTORS             
C                                                                               
      IF((NLIST.EQ.0).AND.(IPR.LE.1)) GO TO 310                                 

c         do iof = 1 , n
c            do jof = 1 , iof
c               fkao(iof,jof)=a(iof,jof)
c               fkao(jof,iof)=a(iof,jof)
c            enddo
c         enddo

      IF(IWLCAO.GT.0) GO TO 160                                                 
      IF(IWNAT.EQ.0) GO TO 150                                                  

C                                                                               
C  NATURAL HYBRID ORBITAL TRANSFORMATION                                        
C
      write(6,*)	
      write(6,*) 'Start NBO transformation after energy convergence is' 
      write(6,*) 'achieved in the AO basis... '
      write(6,*) 

c      CALL NATHYB(B,T,N,NATOMS)
      call aromatic(B,T,N,NATOMS)

  145 IF(IPR.GE.5) CALL ANLYZE(T,NBSZR)
  150 CALL TRANS(N,1)
  
c>>> stored the upper triangle of the fock matrix in fk(i,j) before a(i,j) 
c>>> is detroyed in the successive diagonalization procedure
      do iof = 1 , n
         do jof = 1 , iof
            fk(iof,jof)=a(iof,jof)
            fk(jof,iof)=a(iof,jof)
         enddo
      enddo

      CALL DELETE(N,NLIST,LIST)                                                 
  160 CALL EIGN(N,RHO)
c      IF(IPR.GE.2) CALL SCFOUT(1,2,IWBO)
      GO TO 180                                                                 

C                                                                               
C  ENERGY NOT SATISFIED; DIAGONALIZE FOCK MATRIX; CONTINUE ITERATIONS           
C                                                                               

  170 CONTINUE                                                                  
      OLDENG=ENERGY                                                             
      IF(Z.LT.27) CALL EIGN(N,RHO)                                              
      IF(Z.EQ.27) THEN
         if(iwlcao.eq.0) then
            do iof = 1 , n
               do jof = 1 , iof
                  a(iof,jof)=fk(iof,jof)
               enddo
            enddo
         endif
         CALL DEORTH(N,IPR,RHO,IWBO,fkao)
      ENDIF
  180 CONTINUE                                                                  
C                                                                               
C  EIGENVECTORS (IN B) ARE CONVERTED INTO DENSITY MATRIX (IN B)                 
C                                                                               
      DO 220 I=1,N                                                              
      DO 200 J=I,N                                                              
      XXX(J)=0.0D0                                                              
      DO 190 K=1,OCCA                                                           
  190 XXX(J)= XXX(J)+B(I,K)*B(J,K)*2.0D0                                        
  200 CONTINUE                                                                  
      DO 210 J=I,N                                                              
  210 B(I,J)= XXX(J)                                                            
  220 CONTINUE                                                                  
      DO 230 I=1,N                                                              
      DO 230 J=I,N                                                              
  230 B(J,I)=B(I,J)  

C                                                                               
C  WRITE OUT B.O. DENSITY MATRIX IF Z=26                                        
C                                                                               
C  ANTIBOND CONTRIBUTION TO ELECTRON DENSITY OF OCCUPIED MO'S                   
      IF(Z.NE.26) GO TO 260                                                     
      IF(IWLCAO.GT.0) GO TO 250                                                 
      SUM=0.0D0                                                                 
      BIG=0.0D0                                                                 
      NBIG=0                                                                    
      NL=OCCA+1                                                                 
      DO 240 NAB=NL,N                                                           
      SUM=SUM+B(NAB,NAB)                                                        
      IF(BIG.GE.B(NAB,NAB)) GO TO 240                                           
      NBIG=NAB                                                                  
      BIG=B(NAB,NAB)                                                            
  240 CONTINUE                                                                  
      IF(IPR.GE.1) WRITE(6,1000) SUM,NBIG,BIG                                   
  250 IF(IPR.LT.4) GO TO 260                                                    
      WRITE(6,1100) (WORDS(I,IWBO+1),I=1,3)                                     
      CALL SCFOUT(1,2,IWBO)                                                     
  260 CONTINUE
  
  300 PRINT *, 'Z =',Z
      IF (Z.LT.27) GO TO 10                                                     
  310 CONTINUE                                                                  
      N=NO                                                                      
      IF(IPR.LT.1) RETURN
      WRITE(6,1301)ENERGY
C                                                                               
C  SUM ORBITAL ENERGIES...                                                      
C                                                                               
      ESIG=0.0D0                                                                
      DO 320 I=1,OCCA                                                           
  320 ESIG=ESIG+2.0D0*XXX((3*NBSZR)+I) 
      WRITE(6,1500)ESIG 

      PRINT *, 'OCCA =', OCCA
      
      RETURN                                                                    
 1000 FORMAT(1X,'TOTAL ANTIBOND DENSITY =',F7.4,', LARGEST CON',                
     + 'TRIBUTION D.M.(',I4,') =',F7.4/)                                        
 1100 FORMAT(1X,'DENSITY MATRIX IN ',2A6,A2,' BASIS')                           
 1200 FORMAT(1X,'ENERGY IN TRUNCATED BOND-ORBITAL BASIS SET0')                  
 1300 FORMAT(/,10X,23H ELECTRONIC ENERGY     ,F13.7)                            
 1301 FORMAT(/,10X,21H ELECTRON ENERGY     ,F13.7)                            
 1400 FORMAT(5X,18H ENERGY SATISFIED )                                          
 1500 FORMAT(11X,'SUM OF ORB. ENERGIES  ',F13.7)                                
 1600 FORMAT(1X,'*** CONVERGENCE FAILURE (25 ITERATIONS)...EXIT.')              
      END                                                                       
