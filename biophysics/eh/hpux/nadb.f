      subroutine difov(natom,ndim,mtraj,ntraj)
      
c***************************************************************************
c calculate advanced or retarded overlap matrix
c***************************************************************************

      implicit double precision(a-h,o-z)                                        
      parameter (natmax=200,ndmax=500,nsmp=10)
      logical sp,pd                                                       
      integer sorb                                                              
      character symbol*2,ac*4,segname*4,resname*3
      real*4 rx,ry,rz
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
      common/traj/rx(natmax,nsmp),ry(natmax,nsmp),rz(natmax,nsmp)
      dimension ptr(9),dtr(25),ftr(49)                                          
      equivalence (ptr(3),ca),(ptr(8),cb)                                       
      data sqrt3/1.732050807568877/                                             
      data aui/1.889644746/
      data ptr(9)/0.0/,dtr(12)/0.0/,dtr(22)/0.0/
      data sqrt6/2.449489743/,sqrt10/3.16227766/                                
      data sqrt15/3.872983346/                                                  
      data ftr(15)/0.0/,ftr(22)/0.0/,ftr(43)/0.0/                               
      
      do i=1,ndim
         do j=1,ndim
            s(i,j)=0.0
         enddo
      enddo
      if(natom.eq.1) go to 3000
      isub=0
      nsub=0
      nb=1
      ne=natom

      do 420 i=nb,ne

      ico=i-isub                                                                
      keyi=isym(ico)                                                            
      iorbs=sorb(i)                                                             
      jnd=iorbs-nsub

      do 410 j=nb,ne
      if(j.eq.i.and.mtraj.eq.ntraj) goto 410

      keyj=isym(j)                                                              
      jorbs=sorb(j)

  345 delx=rx(ico,mtraj)-rx(j,ntraj)
      dely=ry(ico,mtraj)-ry(j,ntraj)
      delz=rz(ico,mtraj)-rz(j,ntraj)

      rt2=delx*delx+dely*dely
      r=dsqrt(rt2+delz*delz)
      if(r.lt.1.d-6) then
         write(*,*) ico,mtraj,j,ntraj,r,' warning: pair too close'
c shift a little bit to avoid numerical error
         rx(j,ntraj)=rx(j,ntraj)+0.0001
         ry(j,ntraj)=ry(j,ntraj)+0.0001
         rz(j,ntraj)=rz(j,ntraj)+0.0001
         goto 345
      endif
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

  150 continue
  410 continue
  420 continue
      if (nb.le.natom) goto 3000
      do 440 i=1,ndim                                                           
      if (i.eq.ndim) goto 440                                                   
      l=i+1                                                                     
      do 430 j=l,ndim                                                           
      s(i,j)=s(i,j)-s(j,i)                                                      
  430 s(j,i)=s(i,j)+2.0*s(j,i)                                                  
  440 s(i,i)=2.0*s(i,i)                                                         

 3000 if(mtraj.eq.ntraj) then
         do i=1,ndim
            s(i,i)=1.0
         enddo            
      endif

      return
      end
