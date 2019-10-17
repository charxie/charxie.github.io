C  PRINTS OUT DETAILS OF BOND-ORBITAL TRANSFORMATION FROM MATRIX T.             

      SUBROUTINE ANLYZE(T,NDIM) 
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      INTEGER CH,AN,CZ,U,UL,OCA,OCB                                             
      DIMENSION HYB(4),COEF(3),POW(3),THETA(3),PHI(3),CONJ(2),ORB(3),  N        
     + EL(3),TYPE(3),NAME(18),T(NDIM,NDIM),WORD(2)                              
      COMMON/INFO/NA,CH,MULT,AN(NAZR),C(NAZR,3),NORB,NTA,NBO,IPR
      COMMON/INFO1/CZ(NAZR),U(NBSZR),UL(NAZR),LL(NAZR),NE,OCA,OCB
      COMMON/NHYBRD/IDEN(NBSZR,3)
      COMMON/LBL/LABEL(NBSZR,6),IBX(NBSZR) 
      COMMON/WANT/IWNAT,IWLCAO,NOTOPO                                           
      DATA NAME/'H','HE','LI','BE','B','C','N','O','F','NE',                    
     + 'NA','MG','AL','SI','P','S','CL','AR'/                                   
      DATA TYPE/'SP',' P',' S'/                                                 
      DATA CONJ/' AND','    '/                                                  
      DATA WORD/'NOMIN','NATUR'/                                                
      DATA LLP,LBD,L3C/'LP','BD','3C'/
      
      if(ipr.gt.3) then
         print *
         print *, 'AO-BO transformation matrix'
         print *
         write(6,*) norb
         do i = 1 , norb
            write(6,'(i8,500f8.4)') i,(t(j,i),j=1,norb)
         enddo
      endif
     
      WRITE(6,1100) WORD(IWNAT+1)                                               
      IF(IWNAT.GT.0) GO TO 70                                                   
      WRITE(6,1300)                                                             
      WRITE(6,1400)                                                             
      DO 60 NBOND=1,NORB                                                        
      DO 40 N=1,2                                                               
      I=IDEN(NBOND,N+1)                                                         
      IF(AN(I).GT.0) GO TO 10                                                   
      NEL(N)='GH'                                                               
      POW(N)=0.0                                                                
      COEF(N)=0.0                                                               
      THETA(N)=0.0                                                              
      PHI(N)=0.0                                                                
      ORB(N)='ST'                                                               
      GO TO 40                                                                  
   10 NEL(N)=NAME(AN(I))                                                        
      KL=LL(I)                                                                  
      KU=UL(I)                                                                  
      DO 20 K=1,4                                                               
   20 HYB(K)=0.0                                                                
      KH=0                                                                      
      DO 30 K=KL,KU                                                             
      KH=KH+1   
      PRINT *, 'K NBOND', K, NBOND
   30 HYB(KH)=T(K,NBOND)                                                        
      CALL HTYPE(HYB,COEF(N),POW(N),THETA(N),PHI(N))                            
      ORB(N)='SP'                                                               
      IF(POW(N).EQ.100.0) ORB(N)=' P'                                           
      IF(POW(N).EQ.0.0) ORB(N)=' S'                                             
      IF(POW(N).EQ.100.0) POW(N)=0.0                                            
   40 CONTINUE                                                                  
C  FIND ANGLES FOR LINE OF CENTERS CONNECTING THE ATOMS.                        
      I1=IDEN(NBOND,2)                                                          
      I2=IDEN(NBOND,3)                                                          
      HYB(1)=0.0                                                                
      DO 50 K=2,4                                                               
      HYB(K)=C(I2,K-1)-C(I1,K-1)                                                
   50 CONTINUE                                                                  
      CALL HTYPE(HYB,DUM1,DUM2,TH,PH)                                           
      PRINT 1200,NBOND,(IDEN(NBOND,I),I=1,3),(COEF(I),NEL(I),ORB(I),  PO        
     + W(I),THETA(I),PHI(I),CONJ(I),I=1,2),NEL(1),NEL(2),TH,PH  
   60 CONTINUE 
C     GO TO 130
      GO TO 530   
   70 CONTINUE                                                                  
C  NATURAL HYBRIDS...                                                           
      WRITE(6,1400)                                                             
      DO 130 NBOND=1,NORB                                                       
      DO 80 IB=1,NORB                                                           
      IF(IBX(IB).EQ.NBOND) GO TO 90                                             
   80 CONTINUE                                                                  
   90 LBL=LABEL(IB,1)                                                           
      IF(LBL.EQ.LLP) NCTR=1                                                     
      IF(LBL.EQ.LBD) NCTR=2                                                     
      IF(LBL.EQ.L3C) NCTR=3                                                     
      DO 120 N=1,NCTR                                                           
      I=LABEL(IB,N+3)                                                           
      NEL(N)=NAME(AN(I))                                                        
      KL=LL(I)                                                                  
      KU=UL(I)                                                                  
      DO 100 K=1,4                                                              
  100 HYB(K)=0.0                                                                
      KH=0                                                                      
      DO 110 K=KL,KU                                                            
      KH=KH+1                                                                   
  110 HYB(KH)=T(K,NBOND)                                                        
      CALL HTYPE(HYB,COEF(N),POW(N),THETA(N),PHI(N))                            
      IF(POW(N).EQ.100.0) POW(N)=99.99                                          
  120 CONTINUE                                                                  
      WRITE(6,1000) NBOND,(LABEL(IB,J),J=1,2),(NEL(N),LABEL(IB,N+3),   C        
     +OEF(N),POW(N),THETA(N),PHI(N),N=1,NCTR)                                   
  130 CONTINUE     
  530 CONTINUE
      WRITE(6,1500) 
      RETURN                                                                    
 1000 FORMAT(1X,I3,'. ',A2,A1,2X,3(A1,I3,' (',F7.4,'*SP',F5.2,';(',             
     + F5.1,',',F6.1,'))',3X))                                                  
 1100 FORMAT(//5X,'BOND ORBITAL TRANSFORMATION0 ',A5,'AL HYBRIDS'/)             
 1200 FORMAT(5X,I3,'. ',A3,I2,'-',I3,'0 ',2(F8.4,'*',A2,A2,                     
     + F5.2,' (',F7.3,',',F8.3,')',A4),A2,'-',A2,' (',F7.3,','F8.3,')')         
 1300 FORMAT(7X,'BOND ORBITAL   COEFF.   HYBRID',6X,'(DIRECTION)',9X,           
     + 'COEFF.   HYBRID',6X,'(DIRECTION)',10X,'LINE OF CENTERS')                
 1400 FORMAT(6X,120('-'))                                                       
 1500 FORMAT(/)                                                                 
      END                                                                       


C  PRINTS OUT DETAILS OF BOND-ORBITAL TRANSFORMATION FROM MATRIX T.             

      SUBROUTINE ANLYZE1(T,NDIM,iwanted) 
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      INTEGER CH,AN,CZ,U,UL,OCA,OCB                                             
      DIMENSION HYB(4),T(NDIM,NDIM)
      COMMON/INFO/NA,CH,MULT,AN(NAZR),C(NAZR,3),NORB,NTA,NBO,IPR
      COMMON/INFO1/CZ(NAZR),U(NBSZR),UL(NAZR),LL(NAZR),NE,OCA,OCB
      COMMON/LBL/LABEL(NBSZR,6),IBX(NBSZR) 
      DATA LLP,LBD,L3C/'LP','BD','3C'/
      
      open(15,file='tij',access='append') 
      
      do ib = 1 , norb
	 if(ibx(ib).eq.iwanted) iwant=ib
      enddo
      lbl=label(iwant,1)
      IF(LBL.EQ.LLP) NCTR=1                                                     
      IF(LBL.EQ.LBD) NCTR=2                                                     
      IF(LBL.EQ.L3C) NCTR=3                                                     
      DO 121 ia=1,NCTR                                                           
      I=LABEL(iwant,ia+3)
      KL=LL(I)
      KU=UL(I)
      DO 101 K=1,4                                                              
  101 HYB(K)=0.0                                                                
      KH=0                                                                      
      DO 111 K=KL,KU                                                            
      KH=KH+1                                                                   
  111 HYB(KH)=T(K,iwanted)
      write(15,'(2i5,10f10.5)') iwanted,ia,(hyb(ik),ik=1,kh)
  121 CONTINUE 
      close(15)

      RETURN                                                                    
      END                                                                       
