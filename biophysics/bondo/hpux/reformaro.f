C  SET UP FINAL FORM OF ORTHOGONAL MATRIX T, USING Q AS TEMPORARY STORAGE
c  *** for aromatic rings

      SUBROUTINE REFORMaro(nring,T,NDIM,N,NATOMS,ipr)
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER CZ,U,UL
      character*4 symbol,resname,segname
      integer idres 
      common/iden/idres(nazr),symbol(nazr),resname(nazr),segname(nazr),
     :rx(nazr),ry(nazr),rz(nazr)
      DIMENSION T(NDIM,NDIM)
      COMMON/ARRAYS/A(NBSZR,NBSZR),B(NBSZR,NBSZR),D(NBSZR,NBSZR)
      COMMON/Q/Q(4,NBSZR)
      COMMON/POL/POL(NBSZR,3),IATHY(NBSZR,3),INO(NAZR)
      COMMON/INFO1/CZ(NAZR),U(NBSZR),UL(NAZR),LL(NAZR),
     $NELECS,NOCCA,NCCB
      COMMON/LBL/LABEL(NBSZR,6),IBX(NBSZR)
      common/cring/idc6(6,20)
      common/homo/lumo
      DATA LLP,LBD,L3C,LPB,LSTAR/'LP','BD','3C','PB','*'/
      dimension anorm(6),cop(6,6)
      data cop/1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
     :         2.0, 1.0, 1.0,-1.0,-1.0,-2.0,
     :         0.0,-1.0, 1.0,-1.0, 1.0, 0.0,
     :         0.0, 1.0,-1.0,-1.0, 1.0, 0.0,
     :         2.0,-1.0,-1.0,-1.0,-1.0, 2.0,
     :         1.0,-1.0,-1.0, 1.0, 1.0,-1.0/
      anorm(1)=1.0d0/dsqrt(6.d0)
      anorm(2)=0.5d0/dsqrt(3.d0)
      anorm(3)=0.5d0
      anorm(4)=anorm(3)
      anorm(5)=anorm(2)
      anorm(6)=anorm(1)
       
C  REORDER OCCUPIED BO'S TO PUT LONE PAIRS LAST and pi bonds second last

      DO 10 NLP=1,NOCCA
      IF(LABEL(NLP,1).NE.LLP) GO TO 20
   10 CONTINUE
   20 NLP=NLP-1
      NPB=0
      DO NPI=1,NOCCA
      IF(LABEL(NPI,1).EQ.LPB) NPB=NPB+1
      ENDDO
      NBD=NOCCA-NLP-NPB
      
      do ibd = 1 , n
         if(ibd.le.nlp) then
            ibx(ibd)=ibd+nbd+npb
         elseif(ibd.gt.nlp.and.ibd.le.nocca) then
            ibx(ibd)=ibd-nlp
         else
            ibx(ibd)=ibd
         endif         
      enddo
      
      lumo=nlp+nbd+npb+1
      
c>>> build the part of the transformation matrix correponding to the PI bonds

      do iring = 1 , nring

      jet=0
      do ibd = nlp+nbd+3*iring-2, nlp+nbd+3*iring
         jet=jet+1
         do iat = 1 , 6
            jl=ll(idc6(iat,iring))
            ju=ul(idc6(iat,iring))
            irow=0
            icol=jl+ino(idc6(iat,iring))-1
            do j=jl,ju
               irow=irow+1
               t(j,ibx(ibd))=q(irow,icol)*cop(iat,jet)*anorm(jet)
            enddo
         enddo
      enddo
      do ibd = nocca+nbd+3*iring-2, nocca+nbd+3*iring
         jet=jet+1
         do iat = 1 , 6
            jl=ll(idc6(iat,iring))
            ju=ul(idc6(iat,iring))
            irow=0
            icol=jl+ino(idc6(iat,iring))-1
            do j=jl,ju
               irow=irow+1
               t(j,ibx(ibd))=q(irow,icol)*cop(iat,jet)*anorm(jet)
            enddo
         enddo
      enddo
      
      enddo
      
c>>> j --- AO ; ibx(i) --- BO
  
  281 if(ipr.gt.3) then
         print *, 'TM(AO<--->BO)'
         write(6,'(8x,100i8)') (j,j=1,n)
         do i = 1 , n
            write(6,'(i8,500f8.4)')ibx(i),(t(j,ibx(i)),j=1,n)
         enddo
         print *
      endif

      RETURN                                                                    
      END