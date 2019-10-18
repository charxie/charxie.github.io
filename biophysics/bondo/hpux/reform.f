C  SET UP FINAL FORM OF ORTHOGONAL MATRIX T, USING Q AS TEMPORARY STORAGE        

      SUBROUTINE REFORM(nring,T,NDIM,N,NATOMS,ipr)
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                        
      INTEGER CZ,U,UL                                                           
      DIMENSION T(NDIM,NDIM)                                                    
      DIMENSION S(4,4),EVAL(4),C(4,4),TA(4,4)                                   
      COMMON/Q/Q(4,NBSZR) 
      COMMON/POL/POL(NBSZR,3),IATHY(NBSZR,3),INO(NAZR)
      COMMON/INFO1/CZ(NAZR),U(NBSZR),UL(NAZR),LL(NAZR),
     $NELECS,NOCCA,NCCB
      COMMON/LBL/LABEL(NBSZR,6),IBX(NBSZR) 
      common/cring/idc6(6,20)
      common/homo/lumo
      DATA LLP,LBD,L3C,LPB,LSTAR/'LP','BD','3C','PB','*'/

C  REORDER OCCUPIED BO'S TO PUT LONE PAIRS LAST                                 

      DO 10 NLP=1,NOCCA                                                         
      IF(LABEL(NLP,1).NE.LLP) GO TO 20                                          
   10 CONTINUE                                                                  
   20 NLP=NLP-1                                                                 
      NBD=NOCCA-NLP                                                             
      DO 50 IBD=1,N                                                             
      IF(IBD.GT.NLP) GO TO 30                                                   
C  LONE PAIRS                                                                   
      IBX(IBD)=IBD+NBD                                                          
      GO TO 50                                                                  
C  PAIR BONDS                                                                   
   30 IF(IBD.GT.NOCCA) GO TO 40                                                 
      IBX(IBD)=IBD-NLP                                                          
      GO TO 50                                                                  
C  ANTIBONDS                                                                    
   40 IBX(IBD)=IBD                                                              
   50 CONTINUE                                                                  

      lumo=nbd+nlp+1

C  ZERO ARRAY T...                                                              
      DO 60 I=1,N                                                               
      DO 60 J=1,N                                                               
   60 T(I,J)=0.0D0                                                              
   70 CONTINUE                                                                  
C  SYMMETRIC ORTHOGONALIZATION OF HYBRIDS
      DO 170 IA=1,NATOMS                                                        
      IL=LL(IA)                                                                 
      IU=UL(IA)                                                                 
      NH=IU-IL+1                                                                
      IF(NH.EQ.1) GO TO 170  
C  LOAD IA-BLOCK OF Q INTO TA...                                                
      DO 80 J=1,NH                                                              
      DO 80 I=1,4                                                               
   80 TA(I,J)=Q(I,IL+J-1)                                                       
C  FORM OVERLAP MATRIX S = TA(TRANSP)*TA...                                     
      DO 100 J=1,NH                                                             
      DO 100 I=J,NH                                                             
      TEMP=0.0D0                                                                
      DO 90 K=1,NH                                                              
   90 TEMP=TEMP+TA(K,I)*TA(K,J)                                                 
      S(I,J)=TEMP                                                               
  100 S(J,I)=TEMP                                                               
C  DIAGONALIZE OVERLAP MATRIX...                                                
      CALL JACVEC(NH,S,EVAL,C,4)                                                
C  FORM INVERSE SQUARE ROOT OF S, STORE IN S...                                 
      DO 110 I=1,NH                                                             
  110 EVAL(I)=1.0D0/DSQRT(EVAL(I))                                              
      DO 130 J=1,NH                                                             
      DO 130 I=J,NH                                                             
      TEMP=0.0D0                                                                
      DO 120 K=1,NH                                                             
  120 TEMP=TEMP+EVAL(K)*C(I,K)*C(J,K)                                           
      S(I,J)=TEMP                                                               
  130 S(J,I)=TEMP                                                               
C  FORM NEW TAP=TA*S**(-1/2), STORE IN C...                                     
      DO 150 J=1,NH                                                             
      DO 150 I=1,NH                                                             
      TEMP=0.0D0                                                                
      DO 140 K=1,NH                                                             
  140 TEMP=TEMP+TA(I,K)*S(K,J)                                                  
  150 C(I,J)=TEMP    
C  REPLACE ORTHOGONALIZED TA IN ARRAY Q...                                      
      DO 160 J=1,NH                                                             
      DO 160 I=1,4                                                              
      Q(I,IL+J-1)=C(I,J)                                                        
  160 CONTINUE                                                                  
  170 CONTINUE
C  SYMMETRIC ORTHOGONALIZATION COMPLETE.                                        
c  check the results of orthogonalization
      if(nring.eq.0) goto 2000
      if(ipr.gt.3) print *
      if(ipr.gt.3) print *, 'Orthogonalized hybrids'
      do iring = 1 , nring
         do iat = 1 , 6
            ilc=ll(idc6(iat,iring))
            iuc=ul(idc6(iat,iring))
            nhc=iuc-ilc+1
            if(ipr.gt.3) write(6,'(10i10)') ino(idc6(iat,iring)) 
            do iq = 1 , 4
               if(ipr.gt.3) write(6,'(10f10.5)')
     :              (q(iq,ilc+jq-1),jq=1,ino(idc6(iat,iring)))
            enddo
            do i = 1 , 4
            do j = i , 4
            dotp=0.0
            do k = 1 , 4
            dotp=dotp+q(k,ilc+i-1)*q(k,ilc+j-1)
            enddo
            if(ipr.gt.3) 
     :      write(6,'(1x,a1,i1,a1,i1,a1,f10.5)')'(',i,'*',j,')',dotp
            enddo
            enddo
         enddo
      enddo
2000  continue      
C  DEPOSIT FINAL BOND ORBITALS IN MATRIX T                                      
      NBO=0                                                                     
      DO 280 IBD=1,N                                                            
      KBD=IBD
      IF(LABEL(IBD,2).NE.LSTAR) GO TO 250
C  ANTIBOND ORBITALS0 SEARCH OCCUPIED ORB. LIST TO GET PROPER HYBRIDS...        
C  SEARCH OCCUPIED BOND ORBS. FOR MATCH WITH ANTIBOND ATOMS                     
      DO 240 K=1,NBO                                                            
      DO 180 I=4,6                                                              
      IF(LABEL(K,I).NE.LABEL(IBD,I)) GO TO 240                                  
      IF((LABEL(K,3).LE.0).AND.(LABEL(K,1).EQ.LBD)) GO TO 240                   
C  NEGATIVE IRNK = LABEL(K,3) MEANS BOND ORBITAL WAS ALREADY USED               
  180 CONTINUE                                                                  
C  FOUND MATCH; SET LABEL(K,3)<0 AND RESET POLARIZATION PARAMETERS              
C  FOR ANTIBOND                                                                 
      KBD=K                                                                     
      LABEL(KBD,3)=-LABEL(KBD,3)                                                
      IF(LABEL(IBD,1).NE.L3C) GO TO 230                                         
C  3-CENTER ANTIBONDS FOR (A,B,C)                                               
      NUM=1                                                                     
      IF(LABEL(KBD,3).GT.0) NUM=2                                               
      CALL FM3CAB(POL(KBD,1),POL(KBD,2),POL(KBD,3),NUM,POL(IBD,1),              
     + POL(IBD,2),POL(IBD,3))                                                   
      GO TO 250                                                                 
C  2-CENTER ANTIBOND                                                            
  230 POL(IBD,2)=-POL(KBD,1)                                                    
      POL(IBD,1)=POL(KBD,2)                                                     
      GO TO 250                                                                 
  240 CONTINUE                                                                  
C  COULDN'T FIND SUCCESSFUL MATCH...EXIT                                        
      CALL ABORT('REFORM')                                                      
C  DEPOSIT BOND ORBITALS IN T MATRIX                                            
  250 CONTINUE                                                                  
      DO 270 I=1,3                                                              
      IA=LABEL(IBD,I+3)                                                         
      IF(IA.EQ.0) GO TO 270                                                     
      JL=LL(IA)                                                                 
      JU=UL(IA)                                                                 
      IROW=0                                                                    
      ICOL=JL+IATHY(KBD,I)-1                                                    
      DO 260 J=JL,JU                                                            
      IROW=IROW+1                                                               
      T(J,IBX(IBD))=POL(IBD,I)*Q(IROW,ICOL)
  260 continue
  270 CONTINUE                                                                  
      IF(IBD.EQ.KBD) NBO=IBD                                                    
  280 CONTINUE
      
C  RESTORE LABEL(I,3) > 0                                                       
      DO 290 I=1,N                                                              
      IF(LABEL(I,3).LT.0) LABEL(I,3)=-LABEL(I,3)                                
  290 CONTINUE
  
      RETURN                                                                    
      END                                                                       