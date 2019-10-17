C  CALLED BY SCFCLO ON FINAL PASS (Z=27) FOR DEORTHOGONALIZATION OF             
C  EIGENVECTORS (C ==> S**(-1/2)*C) AND PRINTING OF FOCK AND OVERLAP            
C  MATRICES (ARRAYS A AND D LEFT SYMMETRIZED).                                  

      SUBROUTINE DEORTH(N,IPR,RHO,IWBO,fkao)
      include 'PARAM.include'
      parameter(neigd=500)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      integer TOPO
      character*2 bname
      common/LBL/bname(NBSZR,2),label(nbszr,3),IBX(NBSZR) 
      common/ARRAYS/A(NBSZR,NBSZR),B(NBSZR,NBSZR),D(NBSZR,NBSZR)
      common/GAB/XXX(5*NBSZR),G(NAZR,NAZR),
     :Q(NBSZR),YYY(NBSZR),ENERGY,XXY(214)
      common/EXTRA/T(NBSZR,NBSZR),IHYB(NAZR),TOPO(NAZR,NAZR),
     :T2J(NBSZR)
      common/HIT/NLIST,LIST(40)                                                 
      common/eigv/omcm(neigd),zr(neigd,neigd),zi(neigd,neigd),ne,nblw
      character*6 WORDS(3,2)                                                      
      DATA WORDS/'ATOMIC',' ORBIT','    AL',
     :           'BOND-A','NTIBON','     D'/
      dimension fkao(nbszr,nbszr),fff(n,n),sss(n,n)
      
      IF(IWBO.EQ.0) GO TO 20                                                    
   10 call TRANS(N,3)                                                           
   20 DO 30 I=1,N                                                               
      DO 30 J=I,N                                                               
      D(I,J)=D(J,I)                                                             
      A(I,J)=A(J,I) 
   30 continue      

c      do i = 1 , n
c         do j = 1 , n
c            sss(i,j)=0.0
c            fff(i,j)=fkao(i,j)
c         enddo
c         sss(i,i)=1.0
c         write(6,'(500f8.4)')(fff(i,j),j=1,n)
c      enddo
c      call eigen(n,fff,sss)
c      write(6,'(5f8.4)')(omcm(i),i=1,n)
c      stop
      
c    Print out Fock and overlap matrices 

      open(file='fock', unit=11)
      open(file='ovlp', unit=12)
      do i = 1 , n
         if(bname(i,2).eq.'*') then
            ibeg=i
            goto 5890
         endif
      enddo
 5890 continue
      if(iwbo.eq.0) ibeg=1
      Do I=ibeg,N
         if(iwbo.eq.0)write(11,'(500d10.4)')(27.2*fkao(I,J),J=ibeg,N)
c         if(iwbo.eq.1)write(11,'(500d10.4)')(27.2*A(I,J),J=ibeg,N)
c         write(12,'(500d10.4)') (D(I,J),J=ibeg,N)
      End Do
      close(unit=11)
      close(unit=12)

      IF(IPR.LT.4) GO TO 40                                                     
      WRITE(6,1000) (WORDS(I,IWBO+1),I=1,3)                                     
      call SCFOUT(0,1,IWBO)                                                     
      WRITE(6,1100) (WORDS(I,IWBO+1),I=1,3)                                     
c      call SCFOUT(0,3,IWBO)                                                     
   40 continue                                                                  
C                                                                               
C  DEORTHOGONALIZATION CURRENTLY BYPASSED IF IPR.LT.6...                        
C                                                                               
      IF(IPR.LT.6) return                                                       
C                                                                               
C  STORE OVERLAP MATRIX IN A, FOCK MATRIX IN D...                               
C                                                                               
      DO 50 I=1,N                                                               
      DO 50 J=I,N                                                               
      B(I,J)=A(I,J)                                                             
      A(I,J)=D(I,J)                                                             
      D(I,J)=B(I,J)                                                             
      A(J,I)=A(I,J)                                                             
      D(J,I)=D(I,J)                                                             
   50 continue                                                                  
C                                                                               
C  DIAGONALIZE S, FORM EIGENVALUES OF S**(-1/2) IN T2J...                       
C                                                                               
      call EIGN(N,RHO)                                                          
      DO 60 I=1,N                                                               
      IP=I+240                                                                  
      T2J(I)=1.0D0/DSQRT(XXX(IP))                                               
C                                                                               
C  STORE EIGENVECTORS (U) OF S IN D, PLACE FOCK MATRIX BACK IN A...             
C                                                                               
      DO 60 J=1,N                                                               
      A(I,J)=D(I,J)                                                             
      D(I,J)=B(I,J)                                                             
   60 continue                                                                  
C                                                                               
C  ...FOR TRUNCATION AND DIAGONALIZATION (EIGENVECTORS C IN B).                 
C                                                                               
      call DELETE(N,NLIST,LIST)                                                 
      call EIGN(N,RHO)                                                          
C                                                                               
C  TO FORM S**(-1/2)*C, FIRST STORE ES**(-1/2)*U*C IN B...                      
C                                                                               
      DO 70 I=1,N                                                               
      DO 70 J=1,N                                                               
      A(I,J)=0.0D0                                                              
      DO 70 K=1,N                                                               
      A(I,J)=A(I,J)+D(K,I)*B(K,J)                                               
   70 continue                                                                  
      DO 80 I=1,N                                                               
      DO 80 J=1,N                                                               
      B(I,J)=T2J(I)*A(I,J)                                                      
   80 continue                                                                  
C                                                                               
C  ...THEN MULTIPLY UT*B = S**(-1/2)*C (DEORTHOGONALIZED E.V.) IN A...          
C                                                                               
      DO 90 I=1,N                                                               
      DO 90 J=1,N                                                               
      A(I,J)=0.0D0                                                              
      DO 90 K=1,N                                                               
      A(I,J)=A(I,J)+D(I,K)*B(K,J)                                               
   90 continue                                                                  
C                                                                               
C  ...AND STORE RESULT IN B.                                                    
C                                                                               
      DO 100 I=1,N                                                              
      DO 100 J=1,N                                                              
      B(I,J)=A(I,J)                                                             
  100 continue                                                                  

      return                                                                    
1000  FORMAT(1X,'FOCK MATRIX IN ',2A6,A2,' BASIS')                              
1100  FORMAT(1X,'OVERLAP MATRIX IN ',2A6,A2,' BASIS')                           
      END                                                                       
