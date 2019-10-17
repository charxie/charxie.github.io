      subroutine smat(natom2,natom,ndim)            

c **********************************************************************
c to calculate the Sij matrix                                                   
c **********************************************************************

      implicit real*8(a-h,o-z)
      logical sp,pd
      integer sorb
      character symbol*2,ac*4,keywrd*320,sline*80,segname*4,resname*3
      parameter (natmax=200,ndmax=500)
      common/keywrd/keywrd,sline
      common/dex/add(3),nb,ne,isub,nsub
      common/ato/ac(natmax),symbol(40)
      common/atom1/ns(40),np(40),nd(40),nf(40),inatom,                          
     *isym(natmax),nelem                                                        
      common/atom2/exps(40),exps2(40),expp(40),expp2(40),expd(40),
     *expd2(40),expf(40),expf2(40),cs1(40),cs2(40),cp1(40),cp2(40),
     *cd1(40),cd2(40),cf1(40),cf2(40),couls(40),coulp(40),could(40),
     *coulf(40),x(natmax),y(natmax),z(natmax),ires(natmax),
     *resname(natmax),segname(natmax) 
      common/optv/rho,fermi,wstart,wend,
     *lb,le,lbr,ler,iovlb,iovle,icharg,nelec,
     *distl,idistn,ieb,iee,isb,ise,ihb,ihe,ipope,jpope,irest
      common/sij/s(ndmax,ndmax),h(ndmax,ndmax),c(ndmax),
     *maxs(natmax),maxp(natmax),maxd(natmax),sorb(2*natmax),
     *sp(natmax),pd(natmax)
      data rho/100.0/	

      do i=1,ndim
         do j=1,i
            s(i,j)=0.0
            s(j,i)=0.0
         enddo
      enddo
      iorb=1                                                                    
c                                                                               
c prepare the storing area                                                      
c                                                                               
      do 20 i=1,natom                                                           
      keyi=isym(i)                                                              
      sorb(i)=iorb                                                              
      sorb(i+natom)=iorb+ndim                                                   
      if(ns(keyi).ne.0) iorb=iorb+1                                             
      maxd(i)=nd(keyi)                                                          
      if(np(keyi).ne.0) iorb=iorb+3                                             
      if(nd(keyi).ne.0) iorb=iorb+5                                             
      if(nf(keyi).ne.0) iorb=iorb+7                                             
      maxp(i)=np(keyi)                                                          
      maxs(i)=ns(keyi)                                                          
      sp(i)=exps(keyi).eq.expp(keyi)                                            
      pd(i)=expp(keyi).eq.expd(keyi)                                            
      if(pd(i)) maxp(i)=max0(maxp(i),maxd(i))                                   
      if(sp(i)) maxs(i)=max0(ns(keyi),maxp(i))                                  
   20 continue                                                                  
c                                                                               
c prepare the file to store s(i,j)                                              
c                                                                               
      if(natom.eq.1) go to 30                                                   
      nb=2                                                                      
      ne=natom                                                                  
      isub=0                                                                    
      nsub=0                                                                    
      add(1)=0.0                                                                
      add(2)=0.0                                                                
      add(3)=0.0                                                                
      call movlap(natom2,natom,ndim)                
   30 do 40 i=1,ndim                                                            
      s(i,i)=1.0                                                                
      do 40 j=1,i                                                               
      s(i,j)=s(j,i)                                                             
   40 continue                                                                  

      return                                                                    
      end                 

      
      subroutine movlap(natom2,natom,ndim)
      
c***************************************************************************
c calculate the overlap matrix
c***************************************************************************

      implicit double precision(a-h,o-z)                                        
      logical sp,pd
      integer sorb                                                              
      character symbol*2,ac*4,segname*4,resname*3
      parameter (natmax=200,ndmax=500)
      common/optv/rho,fermi,wstart,wend,
     *ilb,ile,lbr,ler,iovlb,iovle,icharg,nelec,
     *distl,idistn,ieb,iee,isb,ise,iibb,iibe,ipope,jpope,irest
      common/dex/add(3),nb,ne,isub,nsub                                         
      common/ato/ac(natmax),symbol(40)                                          
      common/atom1/ns(40),np(40),nd(40),nf(40),inatom,                          
     *isym(natmax),nelem                                                        
      common/atom2/exps(40),exps2(40),expp(40),expp2(40),expd(40),              
     *expd2(40),expf(40),expf2(40),cs1(40),cs2(40),cp1(40),cp2(40),             
     *cd1(40),cd2(40),cf1(40),cf2(40),couls(40),coulp(40),could(40),            
     *coulf(40),x(natmax),y(natmax),z(natmax),ires(natmax),
     *resname(natmax),segname(natmax) 
      common/loclap/sk1,sk2,r,l1,l2,m,n1,n2,max                                 
      common/sij/s(ndmax,ndmax),h(ndmax,ndmax),c(ndmax),
     *maxs(natmax),maxp(natmax),maxd(natmax),sorb(2*natmax),
     *sp(natmax),pd(natmax)
      dimension ptr(9),dtr(25),ftr(49)                                          

      equivalence (ptr(3),ca),(ptr(8),cb)                                       
      data sqrt3/1.732050807568877/                                             
      data aui/1.889644746/
      data ptr(9)/0.0/,dtr(12)/0.0/,dtr(22)/0.0/
      data sqrt6/2.449489743/,sqrt10/3.16227766/                                
      data sqrt15/3.872983346/                                                  
      data ftr(15)/0.0/,ftr(22)/0.0/,ftr(43)/0.0/                               

      if(nb.gt.natom) go to 10                                                  
      assign 30 to ibrnch                                                       
      go to 20                                                                  
   10 assign 40 to ibrnch                                                       
      im1=natom                                                                 
   20 continue                                                                  

      do 420 i=nb,ne                                                            

      ico=i-isub                                                                
      keyi=isym(ico)                                                            
      iorbs=sorb(i)                                                             
      jnd=iorbs-nsub                                                            

      go to ibrnch,(30,40)                                                      
   30 im1=i-1                                                                   
   40 continue                                                                  

      do 410 j=1,im1                                                            

      keyj=isym(j)                                                              
      jorbs=sorb(j)                                                             

      delx=x(ico)-x(j)+add(1)                                                   
      dely=y(ico)-y(j)+add(2)                                                   
      delz=z(ico)-z(j)+add(3)                                                   

      rt2=delx*delx+dely*dely                                                   
      r=dsqrt(rt2+delz*delz)                                                    
      if(r.gt.rho) r=30.0                                                       
      if(rt2.gt.1.0d-10) go to 50                                               
      cb=1.0                                                                    
      sb=0.0                                                                    
      sa=0.0                                                                    
      goto 60                                                                   
   50 t=dsqrt(rt2)                                                              
      cb=delx/t                                                                 
      sb=dely/t                                                                 
      sa=t/r                                                                    
   60 ca=delz/r                                                                 

c                                                                               
c  the transformation matrices are calculated explicitly.                       
c  ptr is the matrix for projecting the x,y,z orbitals                          
c  onto the local system. the elements are ordered so that first                
c  x then y then z is projected onto the z' axis (sigma).                       
c  then the 3 are projected onto the x' axis and then the y' (pi).              
c  the d orbitals are handled similarly. the order of projection                
c  is x2-y2,z2,xy,xz,yz first onto z2'(sigma)and then onto xz' and              
c  yz'(pi). finally the 5 orbitals are projected onto x'2-y'2 and               
c  then xy' (delta).                                                            
c     for the f-orbitals, the ordering of the orbitals is the same              
c     as the order of the projections:                                          
c     z3 (sigma), xz2 & yz2 (pi), xyz & z(x2-y2) (delta),                       
c     x(x2-3y2) & y(3x2-y2) (phi)                                               
c                                                                               
c  those ptr and dtr which are zero are initialized in a data state-            
c  ment.  ca and cb have been equivalenced to ptr(3) and ptr(8)                 
c     respectively to save time.                                                
c                                                                               
      ptr(1)=sa*cb                                                              
      ptr(2)=sa*sb                                                              
c ... ptr(3)=ca                                                                 
      ptr(4)=ca*cb                                                              
      ptr(5)=ca*sb                                                              
      ptr(6)=-sa                                                                
      ptr(7)=-sb                                                                
c ... ptr(8)=cb                                                                 
c ... ptr(9)=0.0
      if(nd(keyi)+nd(keyj).eq.0) goto 80 
      ca2=ca*ca                                                                 
      sa2=sa*sa                                                                 
      cb2=cb*cb                                                                 
      sb2=sb*sb                                                                 
      cbsb=cb*sb                                                                
      casa=ca*sa                                                                
      cb2sb2=cb2-sb2                                                            
      dtr(1)=sqrt3*0.5*sa2*cb2sb2                                               
      dtr(2)=1.0-1.5*sa2                                                        
      dtr(3)=sqrt3*cbsb*sa2                                                     
      dtr(4)=sqrt3*casa*cb                                                      
      dtr(5)=sqrt3*casa*sb                                                      
      dtr(6)=casa*cb2sb2                                                        
      dtr(7)=-sqrt3*casa                                                        
      dtr(8)=2.0*casa*cbsb                                                      
      dtr(9)=cb*(ca2-sa2)                                                       
      dtr(10)=sb*(ca2-sa2)                                                      
      dtr(11)=-2.0*sa*cbsb                                                      
c ... dtr(12)=0.0                                                               
      dtr(13)=sa* cb2sb2                                                        
      dtr(14)=-ptr(5)                                                           
      dtr(15)=ptr(4)                                                            
      if(nd(keyi)*nd(keyj).eq.0) goto 70                                        
      dtr(16)=0.5*(1.0+ca2)*cb2sb2                                              
      dtr(17)=0.5*sqrt3*sa2                                                     
      dtr(18)=cbsb*(1.0+ca2)                                                    
      dtr(19)=-casa*cb                                                          
      dtr(20)=-casa*sb                                                          
      dtr(21)=-2.0*ca*cbsb                                                      
c ... dtr(22)=0.0                                                               
      dtr(23)=ca*cb2sb2                                                         
      dtr(24)=ptr(2)                                                            
      dtr(25)=-ptr(1)                                                           
   70 if (nf(keyi)+nf(keyj).eq.0) goto 80                                       
      s2b=2.0*sb*cb                                                             
      sa3=sa2*sa                                                                
      c2b=(cb2-sb2)                                                             
      s3b=(c2b*sb + s2b*cb)                                                     
      c3b=(c2b*cb - s2b*sb)                                                     
      s2a=2.0*sa*ca                                                             
      ftr(1)=0.5*ca*(5.0*ca2 -3.0)                                              
      ftr(2)=sqrt6*0.25*cb*sa*(5.0*ca2 -1.0)                                    
      ftr(3)=sqrt6*0.25*sb*sa*(5.0*ca2 -1.0)                                    
      ftr(4)=sqrt15*0.5*s2b*ca*sa2                                              
      ftr(5)=sqrt15*0.5*c2b*ca*sa2                                              
      ftr(6)=sqrt10*0.25*c3b*sa3                                                
      ftr(7)=sqrt10*0.25*s3b*sa3                                                
      ftr(8)=-sqrt6*0.25*sa*(5.0*ca2 -1.0)                                      
      ftr(9)=0.25*cb*ca*(15.0*ca2 -11.0)                                        
      ftr(10)=0.25*sb*ca*(15.0*ca2 -11.0)                                       
      ftr(11)=sqrt10*0.25*s2b*sa*(3.0*ca2 - 1.0)                                
      ftr(12)=sqrt10*0.25*c2b*sa*(3.0*ca2 - 1.0)                                
      ftr(13)=sqrt15*0.25*c3b*ca*sa2                                            
      ftr(14)=sqrt15*0.25*s3b*ca*sa2                                            
c     ftr(15)=0.0                                                               
      ftr(16)=-0.25*sb*(5.0*ca2 -1.0)                                           
      ftr(17)=0.25*cb*(5.0*ca2 -1.0)                                            
      ftr(18)=sqrt10*0.25*c2b*s2a                                               
      ftr(19)=-sqrt10*s2b*s2a*0.25                                              
      ftr(20)=-sqrt15*0.25*s3b*sa2                                              
      ftr(21)=sqrt15*0.25*c3b*sa2                                               
      if (nd(keyi)*nd(keyj).eq.0) goto 80                                       
      c2a=ca2-sa2                                                               
c     ftr(22)=0.0                                                               
      ftr(23)=sqrt10*0.5*sb*ca*sa                                               
      ftr(24)=-sqrt10*0.5*cb*ca*sa                                              
      ftr(25)=c2b*c2a                                                           
      ftr(26)=-s2b*c2a                                                          
      ftr(27)=-sqrt6*0.25*s3b*s2a                                               
      ftr(28)=sqrt6*0.25*c3b*s2a                                                
      ftr(29)=sqrt15*0.5*ca*sa2                                                 
      ftr(30)=sqrt10*0.25*cb*sa*(1.0 -3.0*ca2)                                  
      ftr(31)=sqrt10*0.25*sb*sa*(1.0 -3.0*ca2)                                  
      ftr(32)=0.5*s2b*ca*(3.0*ca2 -1.0)                                         
      ftr(33)=0.5*c2b*ca*(3.0*ca2 -1.0)                                         
      ftr(34)=sqrt6*0.25*c3b*sa*(1.0 + ca2)                                     
      ftr(35)=sqrt6*0.25*s3b*sa*(1.0 + ca2)                                     
      if (nf(keyi)*nf(keyj).eq.0) goto 80                                       
      ftr(36)=-sqrt10*0.25*sa3                                                  
      ftr(37)=sqrt15*0.25*cb*ca*sa2                                             
      ftr(38)=sqrt15*0.25*sb*ca*sa2                                             
      ftr(39)=-sqrt6*0.25*s2b*sa*(1.0 + ca2)                                    
      ftr(40)=-sqrt6*0.25*c2b*sa*(1.0 + ca2)                                    
      ftr(41)=0.25*c3b*ca*(3.0 + ca2)                                           
      ftr(42)=0.25*s3b*ca*(3.0 + ca2)                                           
c     ftr(43)=0.0                                                               
      ftr(44)=-sqrt15*0.25*sb*sa2                                               
      ftr(45)=sqrt15*0.25*cb*sa2                                                
      ftr(46)=-sqrt6*0.25*c2b*s2a                                               
      ftr(47)=sqrt6*0.25*s2b*s2a                                                
      ftr(48)=-0.25*s3b*(1.0 +3.0*ca2)                                          
      ftr(49)=0.25*c3b*(1.0 + 3.0*ca2)                                          
   80 r=r*aui
c
c     (s(i):s(j)).                                                              
c                                                                               
      n2=ns(keyj)
      n1=ns(keyi)
      if(n1.eq.0.or.n2.eq.0) goto 90
      call mov(sigma,pi,delta,phi,ico,j,r,n1,n2,0,0,1)                          
      s(jorbs,jnd)=sigma
c                                                                               
c     (p(i):s(j)).                                                              
c                                                                               
   90 n1=np(keyi)                                                               
      if(n1.eq.0.or.n2.eq.0) goto 110                                           
      call mov(sigma,pi,delta,phi,ico,j,r,n1,n2,1,0,1)                          
      sigma=-sigma                                                              
      do ic=1,3                                                             
         s(jorbs,jnd+ic)=ptr(ic)*sigma
      enddo
c                                                                               
c     (p(i):p(j)).                                                              
c                                                                               
  110 n2=np(keyj)                                                               
      if(n1.eq.0.or.n2.eq.0) goto 130                                           
      call mov(sigma,pi,delta,phi,ico,j,r,n1,n2,1,1,2)                          
      sigma=-sigma                                                              
      do jc=1,3                                                             
         do ic=jc,3                                                            
            s(jorbs+jc,jnd+ic)=ptr(jc)*ptr(ic)*sigma + (ptr(jc+3)*                    
     *      ptr(ic+3)+ptr(jc+6)*ptr(ic+6))*pi                                        
            s(jorbs+ic,jnd+jc)=s(jorbs+jc,jnd+ic)
	 enddo
      enddo
c                                                                               
c     (s(i):p(j)).                                                              
c                                                                               
  130 n1=ns(keyi)                                                               
      if(n1.eq.0.or.n2.eq.0) goto 150                                           
      call mov(sigma,pi,delta,phi,ico,j,r,n1,n2,0,1,1)                          
      do ic=1,3                                                             
         s(jorbs+ic,jnd)=ptr(ic)*sigma
      enddo
c                                                                               
c     (s(i):d(j)).                                                              
c                                                                               
  150 n2=nd(keyj)                                                               
      if(n1.eq.0.or.n2.eq.0) goto 180                                           
      call mov(sigma,pi,delta,phi,ico,j,r,n1,n2,0,2,1)                          
  160 jd=jorbs+3                                                                
      do 170 ic=1,5                                                             
  170 s(jd+ic,jnd)=dtr(ic)*sigma                                                
c                                                                               
c     (p(i):d(j)).                                                              
c                                                                               
  180 n1=np(keyi)                                                               
      if(n1.eq.0.or.n2.eq.0) goto 200                                           
      call mov(sigma,pi,delta,phi,ico,j,r,n1,n2,1,2,2)                          
      sigma=-sigma                                                              
      do 190 ic=1,3                                                             
      do 190 jc=1,5                                                             
  190 s(jd+jc,jnd+ic)=dtr(jc)*ptr(ic)*sigma+(dtr(jc+5)*ptr(ic+3)                
     * +dtr(jc+10)*ptr(ic+6))*pi                                                
c                                                                               
c     (d(i):s(j)).                                                              
c                                                                               
  200 n1=nd(keyi)                                                               
      n2=ns(keyj)                                                               
      if(n1.eq.0.or.n2.eq.0) goto 230                                           
      call mov(sigma,pi,delta,phi,ico,j,r,n1,n2,2,0,1)                          
  210 id=jnd+3                                                                  
      do 220 ic=1,5                                                             
  220 s(jorbs,id+ic)=dtr(ic)*sigma                                              
c                                                                               
c     (d(i):p(j)).                                                              
c                                                                               
  230 n2=np(keyj)                                                               
      if(n1.eq.0.or.n2.eq.0) goto 250                                           
      call mov(sigma,pi,delta,phi,ico,j,r,n1,n2,2,1,2)                          
      pi=-pi                                                                    
      id=jnd+3                                                                  
      do 240 jc=1,3                                                             
      do 240 ic=1,5                                                             
  240 s(jorbs+jc,id+ic)=ptr(jc)*dtr(ic)*sigma+(ptr(jc+3)*dtr(ic+5)              
     * +ptr(jc+6)*dtr(ic+10))*pi                                                
c                                                                               
c     (d(i):d(j)).                                                              
c                                                                               
  250 n2=nd(keyj)                                                               
      if(n1.eq.0.or.n2.eq.0) goto 270                                           
      call mov(sigma,pi,delta,phi,ico,j,r,n1,n2,2,2,3)                          
      pi=-pi                                                                    
      jd=jorbs+3                                                                
      do 260 ic=1,5                                                             
      id=jnd+3                                                                  
      do 260 jc=1,5                                                             
      s(jd+jc,id+ic) =dtr(ic)*dtr(jc)*sigma+(dtr(ic+5)*dtr(jc+5)                
     * +dtr(ic+10)*dtr(jc+10))*pi+(dtr(ic+15)*dtr(jc+15)                        
     * +dtr(ic+20)*dtr(jc+20))*delta                                            
  260 s(jd+ic,id+jc)=s(jd+jc,id+ic)                                             
c                                                                               
c     (s(i):f(j))                                                               
c                                                                               
  270 n1=ns(keyi)                                                               
      n2=nf(keyj)                                                               
      if (n1.eq.0.or.n2.eq.0) go to 290                                         
      call mov(sigma,pi,delta,phi,ico,j,r,n1,n2,0,3,1)                          
      jf=jorbs+8                                                                
      do 280 ic=1,7                                                             
  280 s(jf+ic,jnd)=ftr(ic)*sigma                                                
c                                                                               
c     (p(i):f(j))                                                               
c                                                                               
  290 n1=np(keyi)                                                               
      if (n1.eq.0.or.n2.eq.0) go to 310                                         
      call mov(sigma,pi,delta,phi,ico,j,r,n1,n2,1,3,2)                          
      sigma=-sigma                                                              
      jf=jorbs+8                                                                
      do 300 ic=1,3                                                             
      do 300 jc=1,7                                                             
  300 s(jf+jc,jnd+ic)=ftr(jc)*ptr(ic)*sigma                                     
     * +(ftr(jc+7)*ptr(ic+3) + ftr(jc+14)*ptr(ic+6))*pi                         
c                                                                               
c     (f(i):s(j))                                                               
c                                                                               
  310 n1=nf(keyi)                                                               
      n2=ns(keyj)                                                               
      if (n1.eq.0.or.n2.eq.0) go to 330                                         
      call mov(sigma,pi,delta,phi,ico,j,r,n1,n2,3,0,1)                          
      sigma=-sigma                                                              
      if=jnd+8                                                                  
      do 320 ic=1,7                                                             
  320 s(jorbs,if+ic)=ftr(ic)*sigma                                              
c                                                                               
c     (f(i):p(j))                                                               
c                                                                               
  330 n2=np(keyj)                                                               
      if (n1.eq.0.or.n2.eq.0) go to 350                                         
      call mov(sigma,pi,delta,phi,ico,j,r,n1,n2,3,1,2)                          
      sigma=-sigma                                                              
      if=jnd+8                                                                  
      do 340 jc=1,3                                                             
      do 340 ic=1,7                                                             
  340 s(jorbs+jc,if+ic)=ptr(jc)*ftr(ic)*sigma                                   
     * +(ptr(jc+3)*ftr(ic+7) + ptr(jc+6)*ftr(ic+14))*pi                         
c                                                                               
c     (f(i):d(j))                                                               
c                                                                               
  350 n2=nd(keyj)                                                               
      if (n1.eq.0.or.n2.eq.0) go to 370                                         
      call mov(sigma,pi,delta,phi,ico,j,r,n1,n2,3,2,3)                          
      sigma=-sigma                                                              
      delta=-delta                                                              
      jd=jorbs+3                                                                
      if=jnd+8                                                                  
      do 360 jc=1,5                                                             
      do 360 ic=1,7                                                             
  360 s(jd+jc,if+ic)=ftr(ic)*dtr(jc)*sigma                                      
     * +(ftr(ic+7)*dtr(jc+5) + ftr(ic+14)*dtr(jc+10))*pi                        
     * +(ftr(ic+21)*dtr(jc+15) + ftr(ic+28)*dtr(jc+20))*delta                   
c                                                                               
c     (f(i):f(j))                                                               
c                                                                               
  370 n2=nf(keyj)                                                               
      if (n1.eq.0.or.n2.eq.0) go to 390                                         
      call mov(sigma,pi,delta,phi,ico,j,r,n1,n2,3,3,4)                          
      sigma=-sigma                                                              
      delta=-delta                                                              
      jf=jorbs+8                                                                
      if=jnd+8                                                                  
      do 380 jc=1,7                                                             
      do 380 ic=1,7                                                             
  380 s(jf+jc,if+ic)=ftr(ic)*ftr(jc)*sigma                                      
     * +(ftr(ic+7)*ftr(jc+7) + ftr(ic+14)*ftr(jc+14))*pi                        
     * +(ftr(ic+21)*ftr(jc+21) + ftr(ic+28)*ftr(jc+28))*delta                   
     * +(ftr(ic+35)*ftr(jc+35) + ftr(ic+42)*ftr(jc+42))*phi                     
c                                                                               
c     (d(i):f(j))                                                               
c                                                                               
  390 n1=nd(keyi)                                                               
      if (n1.eq.0.or.n2.eq.0) go to 410                                         
      call mov(sigma,pi,delta,phi,ico,j,r,n1,n2,2,3,3)                          
      pi=-pi                                                                    
      id=jnd+3                                                                  
      jf=jorbs+8                                                                
      do 400 jc=1,7                                                             
      do 400 ic=1,5                                                             
  400 s(jf+jc,id+ic)=dtr(ic)*ftr(jc)*sigma                                      
     * +(dtr(ic+5)*ftr(jc+7) + dtr(ic+10)*ftr(jc+14))*pi                        
     * +(dtr(ic+15)*ftr(jc+21) + dtr(ic+20)*ftr(jc+28))*delta                   
  410 continue                                                                  

  420 continue                                                                  
      if (nb.le.natom) return                                                   
      do 440 i=1,ndim                                                           
      if (i.eq.ndim) goto 440                                                   
      l=i+1                                                                     
      do 430 j=l,ndim                                                           
      s(i,j)=s(i,j)-s(j,i)                                                      
  430 s(j,i)=s(i,j)+2.0*s(j,i)                                                  
  440 s(i,i)=2.0*s(i,i)                                                         

      return                                                                    
      end                                                             
      

      subroutine mov(sigma,pi,delta,phi,i,j,rr,na,nb,la,lb,nn)                  
c                                                                               
c  ********************************************************************         
c  *                                                                  *         
c  *     subroutine mov          called from movlap                   *         
c  *                                                                  *         
c  *                                                                  *         
c  *       mov      subroutine to calculate the overlap               *         
c  *                component independent of the angle between        *         
c  *                the atoms                                         *         
c  *                                                                  *         
c  ********************************************************************         
c                                                                               
c                                                                               
      implicit real*8(a-h,o-z)
      character ac*4,symbol*2,segname*4,resname*3
      parameter (natmax=200)                                                    
c                                                                               
c  exps,expp,expd,expda,expf,expfa have been combined into oexp                 
c  c1,c2,c3,c4                     have been combined into cdf                  
c                                                                               
      common/ato/ac(natmax),symbol(40)                                          
      common/atom1/ns(40),np(40),nd(40),nf(40),natom,                           
     *isym(natmax),nelem                                                        
      common/atom2/oexp(40,8),                                                  
     * cdf(40,8),couls(40),coulp(40),could(40),                                 
     * coulf(40),x(natmax),y(natmax),z(natmax),ires(natmax),
     * resname(natmax),segname(natmax) 
      common/loclap/sk1,sk2,r,l1,l2,m,n1,n2,max                                 
      dimension rll(4),a(30),b(30)                                              
      data oldsk1/-1./,oldsk2/-1./,oldr/-1./                                    

      sigma=0.                                                                  
      pi=0.                                                                     
      delta=0.                                                                  
      phi=0.                                                                    

      ia=1                                                                      
      ib=1                                                                      
      l1=la                                                                     
      l2=lb                                                                     
      n1=na                                                                     
      n2=nb                                                                     
      r=rr                                                                      
      max=na+nb                                                                 
c                                                                               
c  if double zeta orbitals used for d or f                                      
c                                                                               
      keyi=isym(i)                                                              
      keyj=isym(j)                                                              

      idbl2=(la+1)*2                                                            
      if(cdf(keyi,idbl2).ne.0.0) ia=2                                           
      idbl2=(lb+1)*2                                                            
      if(cdf(keyj,idbl2).ne.0.0) ib=2                                           
c                                                                               
c  loop over orbitals                                                           
c                                                                               
      do 30 ik=1,ia                                                             
      do 30 il=1,ib                                                             
      lc=2*la+ik                                                                
      sk1=oexp(keyi,lc)                                                         
      ld=2*lb+il                                                                
      sk2=oexp(keyj,ld)                                                         

C      if(oldsk1.eq.sk1.and.oldsk2.eq.sk2.and.oldr.eq.r) goto 10                 
      call abfns(a,b)
      oldsk1=sk1                                                                
      oldsk2=sk2                                                                
      oldr=r                                                                    
   10 do 20 ij=1,nn                                                             
      m=ij-1                                                                    
      call lovlap(rll(ij),a,b)                                                  
   20 continue                                                                  

      xx=cdf(keyi,lc)                                                           
      yy=cdf(keyj,ld)                                                           

      sigma=sigma+xx*yy*rll(1)
      pi=pi+xx*yy*rll(2)                                                        
      delta=delta+xx*yy*rll(3)                                                  
      phi=phi+xx*yy*rll(4)                                                      
   30 continue                                                                  

      return                                                                    
      end                                                                       


      subroutine abfns(a,b)                                                     
c                                                                               
c  ********************************************************************         
c  *                                                                  *         
c  *     subroutine abfns        called from movlap                   *         
c  *                                                                  *         
c  *                                                                  *         
c  *       abfns    subroutine to calculate the ab functions          *         
c  *                                                                  *         
c  ********************************************************************         
c                                                                               

      implicit real*8(a-h,o-z)
c     this only works for principal quantum # < or = 7                          

      common/loclap/sk1,sk2,rr,l1,l2,m,n1,n2,maxcal                             
      dimension a(30),b(30)                                                     

      j=maxcal+1                                                                
      rho1=0.5*(sk1+sk2)*rr                                                     
      rho2=0.5*(sk1-sk2)*rr                                                     
      if(dabs(rho1).gt.165.0) go to 240                                         
      if(dabs(rho2).gt.165.0) go to 240                                         
      c=dexp(-rho1)                                                             
      a(1)=c/rho1                                                               
      do 10 i=2,j                                                               
   10 a(i)=(dfloat(i-1)*a(i-1)+c)/rho1                                          
      ix=j                                                                      
      ir=dabs(2.0*rho2)                                                         
      is=min0(ir+1,19)                                                          
      if(rho2) 20,210,20                                                        
   20 d=dexp(rho2)                                                              
      h=1.0/d                                                                   
c                                                                               
c  if the value of rho2 is too small the sinh must be obtained                  
c  by summing the infinite series rather than by addition of                    
c  two exponentials.                                                            
c                                                                               
      r=d-h                                                                     
      if(dabs(r)-0.1) 30,60,60                                                  
   30 ra=rho2                                                                   
      rho22=rho2*rho2                                                           
      t=rho2                                                                    
      do 40 i=2,50,2                                                            
      t=t*rho22/dfloat(i*i+i)                                                   
      ra=ra+t                                                                   
      if(t.lt.1.d-30) go to 50                                                  
   40 continue                                                                  
   50 continue                                                                  
      r=ra+ra                                                                   
c                                                                               
c  as many successive b-functions are generated from b-0 by the                 
c  recurrence formulae as accuracy will permit.                                 
c                                                                               
   60 b(1)=r/rho2                                                               
      do 200 i=2,ix,is                                                          
      if(ir.eq.0) go to 120                                                     
   70 il=is-1                                                                   
      if(1.gt.il) go to 110                                                     
      do 100 j=1,il                                                             
      k=i+j-1                                                                   
      if((-1)**k) 80,80,90                                                      
   80 b(k)=(r+dfloat(k-1)*b(k-1))/rho2                                          
      go to 100                                                                 
   90 b(k)=-(d+h-dfloat(k-1)*b(k-1))/rho2                                       
  100 continue                                                                  
  110 continue                                                                  
  120 in=i+is-1                                                                 
c                                                                               
c  after the recurrence formulae have been applied an appropriate               
c  number of times the next b-function is obtained by summation                 
c  of the infinite series.                                                      
c                                                                               
      if(in-ix) 130,130,230                                                     
  130 if((-1)**in) 170,170,140                                                  
  140 tr=rho2                                                                   
  150 b(in)=-2.0*tr/dfloat(in+1)                                                
      do 160 j=1,500                                                            
      tr=tr*rho2*rho2/dfloat((2*j)*(2*j+1))                                     
      if(dabs(tr/b(in))-1.0d-7 ) 200,200,160                                    
  160 b(in)=b(in)-2.0*tr/dfloat(in+1+2*j)                                       
  170 tr=1.                                                                     
  180 b(in)=2.0*tr/dfloat(in)                                                   

      do 190 j=1,500                                                            
      tr=tr*rho2*rho2/dfloat((2*j)*(2*j-1))                                     
      if(dabs(tr/b(in))-1.0d-7 ) 200,200,190                                    
  190 b(in)=b(in)+2.0*tr/dfloat(in+2*j)                                         
  200 continue                                                                  
c                                                                               
c  if the argument is zero a separate formula must be used.                     
c                                                                               
      go to 230                                                                 
  210 do 220 i=1,ix,2                                                           
      b(i)=2.0/dfloat(i)                                                        
  220 b(i+1)=0.0                                                                
  230 return                                                                    
  240 do 250 i=1,20                                                             
      a(i)=0.0                                                                  
  250 b(i)=0.0                                                                  
      go to 230                                                                 

      end                                                                       


      subroutine lovlap (strad,a,b)                                             
c                                                                               
c  ********************************************************************         
c  *                                                                  *         
c  *     subroutine lovlap       called from movlap                   *         
c  *                                                                  *         
c  *                                                                  *         
c  *       lovlap   subroutine to calculate the overlap               *         
c  *                component independent of the angle between        *         
c  *                the atoms                                         *         
c  *                                                                  *         
c  ********************************************************************         
c                                                                               
c                                                                               
      implicit real*8(a-h,o-z)                                        

      common/loclap/sk1,sk2,r,l1,l2,m1,n1,n2,max                                

      dimension a(30),b(30)                                                     
      dimension fact(25)                                                        
      dimension bincoe(8,8)                                                     
      data bincoe/8*1.0,   0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,                     
     * 2*0.0,1.0,3.0,6.0,10.0,15.0,21.0,3*0.0,1.0,4.0,                          
     * 10.0,20.0,35.0,4*0.0,1.0,5.0,15.0,35.0,5*0.0,                            
     * 1.0,6.0,21.0,6*0.0,1.0,7.0,7*0.0,1.0/                                    

      maxxa=1                                                                   
      maxxb=1                                                                   
      maxxc=1                                                                   
      fact(1)=1.0                                                               

      do 10 i=2,25                                                              
   10 fact(i)=fact(i-1)*dfloat(i-1)                                             

      m2=m1                                                                     
      strad=0.0                                                                 
      rhoa=r*sk1                                                                
      rhob=r*sk2                                                                
      rhoap=rhoa**n1                                                            
      rhoap=rhoap*rhoap                                                         
      rhoap=rhoap*rhoa                                                          
      rhoab=rhob**n2                                                            
      rhoab=rhoab*rhoab                                                         
      rhoab=rhoab*rhob                                                          
      rhopo=rhoap*rhoab                                                         
      terma=0.5**(l1+l2+1)*dsqrt(dfloat((l1+l1+1)*(l2+l2+1))*                   
     * fact(l1-m1+1)*fact(l2-m1+1)/(fact(n1+n1+1)*fact(n2+n2+1)*                
     * fact(l1+m1+1)*fact(l2+m1+1))*rhopo)                                      
      jend=1+((l1-m1)/2)                                                        
      kend=1+((l2-m2)/2)                                                        
      ieb=m1+1                                                                  

      do 30 j=1,jend                                                            
      ju=j-1                                                                    
      iab=n1-l1+ju+ju+1                                                         
      icb=l1-m1-ju-ju+1                                                         
      con1=fact(l1+l1-ju-ju+1)/(fact(l1-m1-ju-ju+1)*fact(ju+1)*                 
     * fact(l1-ju+1))                                                           
      do 30 k=1,kend                                                            
      ku=k-1                                                                    
      con12=con1*fact(l2+l2-ku-ku+1)/(fact(l2-m2-ku-ku+1)*fact(ku+1)*           
     * fact(l2-ku+1))                                                           
      iev=ju+ku+l2                                                              
      if(2*(iev/2).ne.iev) con12=-con12                                         
      ibb=n2-l2+ku+ku+1                                                         
      idb=l2-m2-ku-ku+1                                                         
      value=0.0                                                                 
      do 20 i6=1,ieb                                                            
      do 20 i5=1,ieb                                                            
      value1=bincoe(ieb,i6)*bincoe(ieb,i5)                                      
      iev=i5+i6                                                                 
      if(2*(iev/2).ne.iev) value1=-value1                                       
      do 20 i4=1,idb                                                            
      value1=-value1                                                            
      value2=bincoe(idb,i4)*value1                                              
      do 20 i3=1,icb                                                            
      value3=bincoe(icb,i3)*value2                                              
      do 20 i2=1,ibb                                                            
      value3=-value3                                                            
      value4=bincoe(ibb,i2)*value3                                              
      do 20 i1=1,iab                                                            
      term=value4*bincoe(iab,i1)                                                
      ir=i1+i2+ieb+ieb-i6-i6-i3+idb-i4+icb-1                                    
      ip=iab-i1+ibb-i2+ieb+ieb-i5-i5+icb-i3+idb-i4+1                            
   20 value=value+a(ip)*b(ir)*term                                              
   30 strad=strad+value*con12                                                   

      strad=strad*terma                                                         

      return                                                                    
      end