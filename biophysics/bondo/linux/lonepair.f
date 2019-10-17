C  MODIFY ORTHOGONAL MATRIX T FOR TRANSFORMATION FROM AO'S TO               
C  NATURAL HYBRID BOND ORBITALS, USING INPUT FOCK MATRIX 'FM'. 
c  *** for O lone pairs only
c  *** in input file ipr must be set > 3 as to get dm          

      SUBROUTINE lonepair(FM,DM,N,ipr)  
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      real*8 kappa
      character*2 bname
      COMMON/LBL/bname(nbszr,2),label(NBSZR,3),IBX(NBSZR) 
      COMMON/EXTRA/T(NBSZR,NBSZR),IHYB(NAZR),TOPO(NAZR,NAZR),
     :T2J(NBSZR)
      DIMENSION FM(NBSZR,NBSZR),DM(NBSZR,NBSZR)
      character*4 symbol,resname,segname
      integer idres 
      common/iden/idres(nazr),symbol(nazr),resname(nazr),segname(nazr),
     :rx(nazr),ry(nazr),rz(nazr)
      dimension tnew(nbszr,nbszr)
      dimension hm(n,n)
      dimension index_lonepair(nbszr),ind1_lonepair(nbszr)
      dimension h(2,2),v(2,2),d(2,2),d1(2,2)

      do i = 1 , n
         do j = 1 , n
            tnew(i,j)=t(i,j)
         enddo
      enddo
      
      num_lonepair=0
      
      do ibd = 1 , n
      
         if(symbol(label(ibd,2))(2:2).eq.'O'
     :      .and.bname(ibd,1).eq.'LP') then
     
            num_lonepair=num_lonepair+1
            index_lonepair(num_lonepair)=ibx(ibd)
            ind1_lonepair(num_lonepair)=ibd
            
         endif
    
      enddo
      
      if(ipr.gt.3) then
         print *,'the molecule contains ',num_lonepair,' O lonepair'
         write(6,'(100i4)') (index_lonepair(i),i=1,num_lonepair)
         print *
      endif
      
      do i = 1 , num_lonepair, 2

         if(bname(ind1_lonepair(i),2).ne.bname(ind1_lonepair(i+1),2))
     :   stop 'wrong O lone pair found!'
      
	 h(1,1)=fm(index_lonepair(i  ),index_lonepair(i  ))
	 h(1,2)=fm(index_lonepair(i  ),index_lonepair(i+1))
	 h(2,1)=fm(index_lonepair(i+1),index_lonepair(i  ))
	 h(2,2)=fm(index_lonepair(i+1),index_lonepair(i+1))
	 d(1,1)=dm(index_lonepair(i  ),index_lonepair(i  ))
	 d(1,2)=dm(index_lonepair(i  ),index_lonepair(i+1))
	 d(2,1)=dm(index_lonepair(i+1),index_lonepair(i  ))
	 d(2,2)=dm(index_lonepair(i+1),index_lonepair(i+1))
	 
	 if(ipr.gt.3) then
	    print *,'print the subblock dm for the',int(i/2)+1,
     :           '-th O lone pair'
 	    write(6,'(2f10.5)') ((d(m,k),k=1,2),m=1,2)
         endif
	 
	 delta=h(1,1)-h(2,2)
	 hplus=h(1,1)+h(2,2)
	 kappa=dsqrt(delta*delta+4.d0*h(1,2)*h(2,1))
c	 print *
c	 write(6,'(2f10.5)') 0.5*(hplus+kappa),0.5*(hplus-kappa)
	 if(h(1,2).gt.0) then
	    v(1,1)= dsqrt((kappa+delta)/2.0/kappa)
	    v(1,2)= dsqrt((kappa-delta)/2.0/kappa)
	    v(2,1)= dsqrt((kappa-delta)/2.0/kappa)
	    v(2,2)=-dsqrt((kappa+delta)/2.0/kappa)
	 else
	    v(1,1)=-dsqrt((kappa+delta)/2.0/kappa)
	    v(1,2)= dsqrt((kappa-delta)/2.0/kappa)
	    v(2,1)=-dsqrt((kappa-delta)/2.0/kappa)
	    v(2,2)=-dsqrt((kappa+delta)/2.0/kappa)
	 endif
c	 print *
c	 write(6,'(2f10.5)') ((h(m,k),k=1,2),m=1,2)
c	 write(6,'(2f10.5)') ((v(m,k),k=1,2),m=1,2)
	 
	 d1(1,1)=v(1,1)*v(1,1)*d(1,1)+v(1,1)*v(1,2)*d(1,2)
     :          +v(1,2)*v(1,1)*d(2,1)+v(1,2)*v(1,2)*d(2,2)	 
	 d1(1,2)=v(1,1)*v(2,1)*d(1,1)+v(1,1)*v(2,2)*d(1,2)
     :          +v(1,2)*v(2,1)*d(2,1)+v(1,2)*v(2,2)*d(2,2)	 
	 d1(2,1)=v(2,1)*v(1,1)*d(1,1)+v(2,1)*v(1,2)*d(1,2)
     :          +v(2,2)*v(1,1)*d(2,1)+v(2,2)*v(1,2)*d(2,2)	 
	 d1(2,2)=v(2,1)*v(2,1)*d(1,1)+v(2,1)*v(2,2)*d(1,2)
     :          +v(2,2)*v(2,1)*d(2,1)+v(2,2)*v(2,2)*d(2,2)
	 
	 if(ipr.gt.3) print *
	 if(ipr.gt.3) write(6,'(2f10.5)') ((d1(m,k),k=1,2),m=1,2)
	 
	 do k = 1 , n
	    tnew(k,index_lonepair(i  ))=v(1,1)*t(k,index_lonepair(i  ))
     :                               +v(1,2)*t(k,index_lonepair(i+1))
	    tnew(k,index_lonepair(i+1))=v(2,1)*t(k,index_lonepair(i  ))
     :                               +v(2,2)*t(k,index_lonepair(i+1))
	 enddo
	 
	 
      enddo
      
      do i = 1 , n
         do j = 1 , i
            t(i,j)=tnew(i,j)
            t(j,i)=tnew(j,i)
         enddo
      enddo

      return
      END                                                                       

      SUBROUTINE lonepair1(n,ipr)  
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*2 bname
      COMMON/LBL/bname(nbszr,2),label(NBSZR,3),IBX(NBSZR) 
      COMMON/EXTRA/T(NBSZR,NBSZR),IHYB(NAZR),TOPO(NAZR,NAZR),
     $T2J(NBSZR)
      character*4 symbol,resname,segname
      integer idres 
      common/iden/idres(nazr),symbol(nazr),resname(nazr),segname(nazr),
     :rx(nazr),ry(nazr),rz(nazr)
      dimension tnew(nbszr,nbszr)
      dimension index_lonepair(nbszr),ind1_lonepair(nbszr)
      dimension v(2,2)

      v(1,1)= 1.0/dsqrt(2.d0)
      v(1,2)= 1.0/dsqrt(2.d0)
      v(2,1)= 1.0/dsqrt(2.d0)
      v(2,2)=-1.0/dsqrt(2.d0)

      do i = 1 , n
         do j = 1 , n
            tnew(i,j)=t(i,j)
         enddo
      enddo
      
      num_lonepair=0
      
      do ibd = 1 , n
      
         if(symbol(label(ibd,2))(2:2).eq.'O'
     :      .and.bname(ibd,1).eq.'LP') then
     
            num_lonepair=num_lonepair+1
            index_lonepair(num_lonepair)=ibx(ibd)
            ind1_lonepair(num_lonepair)=ibd
            
         endif
    
      enddo
      
      if(ipr.gt.3) then
         print *,'the molecule contains ',num_lonepair,' O lonepair'
         write(6,'(100i4)') (index_lonepair(i),i=1,num_lonepair)
         print *
      endif
      
      do i = 1 , num_lonepair, 2

         if(label(ind1_lonepair(i),2).ne.label(ind1_lonepair(i+1),2))
     :   stop 'wrong O lone pair found!'
      
	 do k = 1 , n
	    tnew(k,index_lonepair(i  ))=v(1,1)*t(k,index_lonepair(i  ))
     :                               +v(1,2)*t(k,index_lonepair(i+1))
	    tnew(k,index_lonepair(i+1))=v(2,1)*t(k,index_lonepair(i  ))
     :                               +v(2,2)*t(k,index_lonepair(i+1))
	 enddo
	 
      enddo
      
      do i = 1 , n
         do j = 1 , i
            t(i,j)=tnew(i,j)
            t(j,i)=tnew(j,i)
         enddo
      enddo

      return
      END                                                                       



      SUBROUTINE lonepair2(FM,DM,N,ipr)  
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      real*8 kappa
      character*2 bname
      COMMON/LBL/bname(nbszr,2),label(NBSZR,3),IBX(NBSZR)
      COMMON/EXTRA/T(NBSZR,NBSZR),IHYB(NAZR),TOPO(NAZR,NAZR),
     $T2J(NBSZR)
      DIMENSION FM(NBSZR,NBSZR),DM(NBSZR,NBSZR)
      character*4 symbol,resname,segname
      integer idres 
      common/iden/idres(nazr),symbol(nazr),resname(nazr),segname(nazr),
     :rx(nazr),ry(nazr),rz(nazr)
      dimension tnew(nbszr,nbszr)
      dimension hm(n,n)
      dimension index_lonepair(nbszr),ind1_lonepair(nbszr)
      dimension h(2,2),v(2,2),d(2,2),d1(2,2)

      do i = 1 , n
         do j = 1 , n
            tnew(i,j)=t(i,j)
         enddo
      enddo
      
      num_lonepair=0
      
      do ibd = 1 , n
      
         if((symbol(label(ibd,2))(2:3).eq.'OD'
     :      .or.symbol(label(ibd,2))(2:3).eq.'OG')
     :      .and.bname(ibd,1).eq.'LP') then
     
            num_lonepair=num_lonepair+1
            index_lonepair(num_lonepair)=ibx(ibd)
            ind1_lonepair(num_lonepair)=ibd
            
         endif
    
      enddo
      
      if(ipr.gt.3) then
         print *,'the molecule contains ',num_lonepair,' O lonepair'
         write(6,'(100i4)') (index_lonepair(i),i=1,num_lonepair)
         print *
      endif
      
      do i = 1 , num_lonepair, 2

         if(label(ind1_lonepair(i),2).ne.label(ind1_lonepair(i+1),2))
     :   stop 'wrong O lone pair found!'
      
	 h(1,1)=fm(index_lonepair(i  ),index_lonepair(i  ))
	 h(1,2)=fm(index_lonepair(i  ),index_lonepair(i+1))
	 h(2,1)=fm(index_lonepair(i+1),index_lonepair(i  ))
	 h(2,2)=fm(index_lonepair(i+1),index_lonepair(i+1))
	 d(1,1)=dm(index_lonepair(i  ),index_lonepair(i  ))
	 d(1,2)=dm(index_lonepair(i  ),index_lonepair(i+1))
	 d(2,1)=dm(index_lonepair(i+1),index_lonepair(i  ))
	 d(2,2)=dm(index_lonepair(i+1),index_lonepair(i+1))
	 
	 if(ipr.gt.3) then
	    print *,'print the subblock dm for the',int(i/2)+1,
     :           '-th O lone pair'
 	    write(6,'(2f10.5)') ((d(m,k),k=1,2),m=1,2)
         endif
	 
	 delta=h(1,1)-h(2,2)
	 hplus=h(1,1)+h(2,2)
	 kappa=dsqrt(delta*delta+4.d0*h(1,2)*h(2,1))
c	 print *
c	 write(6,'(2f10.5)') 0.5*(hplus+kappa),0.5*(hplus-kappa)
	 if(h(1,2).gt.0) then
	    v(1,1)= dsqrt((kappa+delta)/2.0/kappa)
	    v(1,2)= dsqrt((kappa-delta)/2.0/kappa)
	    v(2,1)= dsqrt((kappa-delta)/2.0/kappa)
	    v(2,2)=-dsqrt((kappa+delta)/2.0/kappa)
	 else
	    v(1,1)=-dsqrt((kappa+delta)/2.0/kappa)
	    v(1,2)= dsqrt((kappa-delta)/2.0/kappa)
	    v(2,1)=-dsqrt((kappa-delta)/2.0/kappa)
	    v(2,2)=-dsqrt((kappa+delta)/2.0/kappa)
	 endif
c	 print *
c	 write(6,'(2f10.5)') ((h(m,k),k=1,2),m=1,2)
c	 write(6,'(2f10.5)') ((v(m,k),k=1,2),m=1,2)
	 
	 d1(1,1)=v(1,1)*v(1,1)*d(1,1)+v(1,1)*v(1,2)*d(1,2)
     :          +v(1,2)*v(1,1)*d(2,1)+v(1,2)*v(1,2)*d(2,2)	 
	 d1(1,2)=v(1,1)*v(2,1)*d(1,1)+v(1,1)*v(2,2)*d(1,2)
     :          +v(1,2)*v(2,1)*d(2,1)+v(1,2)*v(2,2)*d(2,2)	 
	 d1(2,1)=v(2,1)*v(1,1)*d(1,1)+v(2,1)*v(1,2)*d(1,2)
     :          +v(2,2)*v(1,1)*d(2,1)+v(2,2)*v(1,2)*d(2,2)	 
	 d1(2,2)=v(2,1)*v(2,1)*d(1,1)+v(2,1)*v(2,2)*d(1,2)
     :          +v(2,2)*v(2,1)*d(2,1)+v(2,2)*v(2,2)*d(2,2)
	 
	 if(ipr.gt.3) print *
	 if(ipr.gt.3) write(6,'(2f10.5)') ((d1(m,k),k=1,2),m=1,2)
	 
	 do k = 1 , n
	    tnew(k,index_lonepair(i  ))=v(1,1)*t(k,index_lonepair(i  ))
     :                               +v(1,2)*t(k,index_lonepair(i+1))
	    tnew(k,index_lonepair(i+1))=v(2,1)*t(k,index_lonepair(i  ))
     :                               +v(2,2)*t(k,index_lonepair(i+1))
	 enddo
	 
	 
      enddo
      
      do i = 1 , n
         do j = 1 , i
            t(i,j)=tnew(i,j)
            t(j,i)=tnew(j,i)
         enddo
      enddo

      return
      END