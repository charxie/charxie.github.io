C  CNDO/INDO CLOSED SHELL SCF SEGMENT 
C  GAMMA MATRIX CONTAINED IN G, CORE HAMILTONIAN CONTAINED IN Q AND  
C  UPPER TRIANGLE OF A, AND INITIAL DENSITY MATRIX CONTAINED IN B    
C  OPTIONS   CNDO OR INDO 

      SUBROUTINE SCFCLO
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
      character option*4,opnclo*3
      integer HUCKEL,CNDO,INDO,CLOSED,OPEN 
      integer CHARGE,OCCA,OCCB,ULIM,U,AN,CZ,Z,TOPO 
      character*4 symbol,resname,segname
      character*2 bname
      logical accel
      common/iden/idres(nazr),symbol(nazr),resname(nazr),segname(nazr),
     :rx(nazr),ry(nazr),rz(nazr)
      common/ARRAYS/A(NBSZR,NBSZR),B(NBSZR,NBSZR),D(NBSZR,NBSZR)
      common/GAB/XXX(5*NBSZR),G(NAZR,NAZR),
     :Q(NBSZR),YYY(NBSZR),ENERGY,XXY(214) 
      common/INFO/NATOMS,CHARGE,MULTIP,AN(NAZR),C(NAZR,3),N,NTA,NO,IPR
      common/INFO1/CZ(NAZR),U(NBSZR),ULIM(NAZR),LLIM(NAZR),
     :NELECS,OCCA,OCCB
      common/OPSION/OPTION,OPNCLO,HUCKEL,CNDO,INDO,CLOSED,OPEN       
      common/EXTRA/T(NBSZR,NBSZR),IHYB(NAZR),TOPO(NAZR,NAZR),
     :T2J(NBSZR)
      common/WANT/IWNAT,IWLCAO,NOTOPO
      common/LBL/bname(nbszr,2),label(NBSZR,3),IBX(NBSZR) 
      dimension fk(nbszr,nbszr)
      dimension dm(nbszr,nbszr),fkao(nbszr,nbszr)
      dimension G1(18),F2(18)       
      common/HIT/NLIST,LIST(40)
      character*6 WORDS(3,2)
      DATA WORDS/'ATOMIC',' ORBIT','    AL',
     :           'BOND-A','NTIBON','     D'/
     
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
      IT=45
      it1=it+1
      it2=it+2
      RHO=1.D-6 
      IWBO=1  
      IF(IWLCAO.GT.0) IWBO=0 
      OLDENG=0.0D0
      accel=.true.

      PRINT*, 'N=', N

      if(accel) then
        open(31,file='dmold')
        read(31,*) n0
        if(n0.ne.n) goto 1080
        do idm = 1 , n
          read(31,'(500f10.5)')(b(idm,jdm),jdm=1,n)
        enddo
        close(31)
      endif
 1080 continue

c GRAND SCF-ITERATIONS LOOP, Z.LT.it
c     Z.LE.it  SCF ITERATIONS BEFORE CONVERGENCE 
c     Z.EQ.it1 B.O. TRANSFORMATION, BASIS-SET TRUNCATION  
c     Z.EQ.it2 DEORTHOGONALIZATION AND EXIT    

   10 continue 
      IF(Z.NE.it) GO TO 20 
      WRITE(6,1600)      
      call EXIT 
   20 continue 
      Z = Z+1
      if(z.eq.it2) goto 170
      ENERGY = 0.D0      
C  
C CONSTRUCT FOCK MATRIX (IN A) FROM DENSITY MATRIX (IN B)0 
C ...TRANSFER CORE HAMILTONIAN TO LOWER TRIANGLE OF A... 
C
      do 30 I=1,N 
      A(I,I)=Q(I) 
      do 30 J=I,N 
   30 A(J,I)=A(I,J)      
      do 40 I=1,N 
      II=U(I) 
      A(I,I)=A(I,I)-B(I,I)*G(II,II)*0.5D0      
      do 40 K=1,N 
      JJ=U(K) 
   40 A(I,I)=A(I,I)+B(K,K)*G(II,JJ) 
      NM=N-1  
      do 50 I=1,NM       
      II=U(I) 
      LL=I+1  
      do 50 J=LL,N       
      JJ=U(J) 
   50 A(J,I)=A(J,I)-B(J,I)*G(II,JJ)*0.5D0      
C INDO MODIFICATION     
      IF (OPTION.EQ.'CNDO') GO TO 100
      write(6,*) 'INDO'
   60 do 90 II=1,NATOMS  
      K=AN(II) 
      I=LLIM(II) 
      IF (K.EQ.1) GO TO 90 
   70 PAA=B(I,I)+B(I+1,I+1)+B(I+2,I+2)+B(I+3,I+3) 
      A(I,I)=A(I,I)-(PAA-B(I,I))*G1(K)/6.D0   
      do 80 J=1,3 
      A(I+J,I+J)=A(I+J,I+J)-B(I,I)*G1(K)/6.D0-(PAA-B(I,I))*7.D0*F2(K)/5 
     :      0.D0+B(I+J,I+J)*11.D0*F2(K)/50.D0  
   80 A(I+J,I)=A(I+J,I)+B(I,I+J)*G1(K)/2.D0   
      I1=I+1  
      I2=I+2  
      I3=I+3  
      A(I2,I1)=A(I2,I1)+B(I2,I1)*11.D0*F2(K)/50.D0 
      A(I3,I1)=A(I3,I1)+B(I3,I1)*11.D0*F2(K)/50.D0 
      A(I3,I2)=A(I3,I2)+B(I3,I2)*11.D0*F2(K)/50.D0 
   90 continue
C  
C FOCK MATRIX COMPLETE; EVALUATE AND PRINT OUT ENERGY
C  
  100 continue 
      do 110 I=1,N       
  110 ENERGY=ENERGY+0.5D0*B(I,I)*(A(I,I)+Q(I)) 
      do 120 I=1,NM      
      LL=I+1  
      do 120 J=LL,N      
  120 ENERGY=ENERGY+B(I,J)*(A(I,J)+A(J,I))     

      IF((NLIST.EQ.0).AND.(Z.EQ.it2))GO TO 170
      IF((Z.EQ.it2).AND.(IPR.GE.1)) WRITE(6,1200)
      IF(IPR.GE.3) WRITE(6,1300) ENERGY 
      IF(DABS(ENERGY-OLDENG).GE.0.00005D0) GO TO 170 
C  
C ENERGY SATISFIED; SET Z=it1 FOR FINAL TWO PASSES
C  
  130 IF(Z.LE.IT)Z=it1
  140 IF((Z.LT.it2).AND.(IPR.GE.1)) WRITE(6,1400)
      IF(IPR.GE.1) WRITE(6,1300) ENERGY
      IF(Z.EQ.it2)GO TO 170
C  
C Z = it1 TRANSFORM, DELETE, DIAGONALIZE, PRINT B.O. EIGENVECTORS
C  
      IF((NLIST.EQ.0).AND.(IPR.LE.1)) GO TO 310

c      if(iwlcao.eq.0) then
         do iof = 1 , n
           do jof = 1 , iof
             fkao(iof,jof)=a(iof,jof)
             fkao(jof,iof)=a(iof,jof)
           enddo
         enddo
         open(31,file='dmold')
         rewind 31
         write(31,*) n
         do iof = 1 , n
           write(31,'(500f10.5)')(b(iof,jof),jof=1,n)
         enddo
         close(31)
c      endif
      
      IF(IWLCAO.GT.0) GO TO 160     
      IF(IWNAT.EQ.0) GO TO 150      

C  
C NATURAL HYBRID ORBITAL TRANSFORMATION       
C
      write(6,*)	
      write(6,*) 'Start NBO transformation after energy convergence is' 
      write(6,*) 'achieved in the AO basis... '
      write(6,*) 

c      call NATHYB(B,T,N,NATOMS)
      call aromatic(B,T,N,NATOMS,ipr)
c      call compsignaro(n,ipr)

  145 IF(IPR.GE.5) call ANLYZE(T,NBSZR)
  150 call TRANS(N,1)
  
c>>> stored the upper triangle of the fock matrix in fk(i,j) before a(i,j) 
c>>> is detroyed in the successive diagonalization procedure
      do iof = 1 , n
        do jof = 1 , iof
          fk(iof,jof)=a(iof,jof)
          fk(jof,iof)=a(iof,jof)
        enddo
      enddo

      call DELETE(N,NLIST,LIST)     
  160 call EIGN(N,RHO)
      IF(IPR.GE.5) call SCFOUT(1,2,IWBO)
      GO TO 180 

C  
C ENERGY NOT SATISFIED; DIAGONALIZE FOCK MATRIX; continue ITERATIONS 
C  

  170 continue 
      OLDENG=ENERGY      
      IF(Z.LT.it2) call EIGN(N,RHO)  
      IF(Z.EQ.it2) THEN
      if(iwlcao.eq.0) then
        do iof = 1 , n
          do jof = 1 , iof
            a(iof,jof)=fk(iof,jof)
          enddo
        enddo
      endif
      call DEORTH(N,IPR,RHO,IWBO,fkao)
      ENDIF
  180 continue 
C  
C EIGENVECTORS (IN B) ARE CONVERTED INTO DENSITY MATRIX (IN B)      
C  
      do 220 I=1,N       
      do 200 J=I,N       
      XXX(J)=0.0D0       
      do 190 K=1,OCCA    
  190 XXX(J)= XXX(J)+B(I,K)*B(J,K)
  200 continue 
      do 210 J=I,N       
  210 B(I,J)= XXX(J)*2.d0
  220 continue 
      do 230 I=1,N       
      do 230 J=I,N       
  230 B(J,I)=B(I,J)  
C
C WRITE OUT B.O. DENSITY MATRIX IF Z=it1
C
C  ANTIBOND CONTRIBUTION TO ELECTRON DENSITY OF OCCUPIED MO'S 
      IF(Z.NE.it1) GO TO 260
      IF(IWLCAO.GT.0) GO TO 250
      SUM=0.0D0
      BIG=0.0D0
      NBIG=0
      NL=OCCA+1
      do 240 NAB=NL,N
      SUM=SUM+B(NAB,NAB)
      IF(BIG.GE.B(NAB,NAB)) GO TO 240
      NBIG=NAB
      BIG=B(NAB,NAB)
  240 continue

      if(iwlcao.eq.0) then
        do iof = 1 , n
          do jof = 1 , iof
             dm(iof,jof)=b(iof,jof)
             dm(jof,iof)=b(iof,jof)
          enddo
        enddo
      endif
      
      IF(IPR.GE.1) WRITE(6,1000) SUM,NBIG,BIG  
  250 IF(IPR.LT.4) GO TO 260 
      WRITE(6,1100) (WORDS(I,IWBO+1),I=1,3)    
      call SCFOUT(1,2,IWBO) 
  260 continue
  300 PRINT *, 'Z =',Z
      IF (Z.LT.it2) GO TO 10 
  310 continue 
      N=NO
      IF(IPR.LT.1) return
      WRITE(6,1301)ENERGY
C  
C SUM ORBITAL ENERGIES... 
C  
      ESIG=0.0D0 
      do 320 I=1,OCCA    
  320 ESIG=ESIG+4.0D0*XXX((3*NBSZR)+I) 
      WRITE(6,1500)ESIG 

      PRINT *, 'OCCA =', OCCA

      if(iwlcao.eq.0) then
        call cobond(fk,fkao,dm,t,n,ipr)
      endif

      return
      
 1000 FORMAT(1X,'TOTAL ANTIBOND DENSITY =',F7.4,', LARGEST CON',     
     : 'TRIBUTION D.M.(',I4,') =',F7.4/)       
 1100 FORMAT(1X,'DENSITY MATRIX IN ',2A6,A2,' BASIS')     
 1200 FORMAT(1X,'ENERGY IN TRUNCATED BOND-ORBITAL BASIS SET0')       
 1300 FORMAT(/,10X,23H ELECTRONIC ENERGY     ,F20.10)      
 1301 FORMAT(/,10X,21H ELECTRON ENERGY     ,F20.10)      
 1400 FORMAT(5X,18H ENERGY SATISFIED ) 
 1500 FORMAT(11X,'SUM OF ORB. ENERGIES  ',F20.10) 
 1600 FORMAT(1X,'*** CONVERGENCE FAILURE (45 ITERATIONS)...EXIT.')   
      END