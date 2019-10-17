
      subroutine hmat(natom,ndim)                                         

c *************************************************************************
c calculates the Hij matrix                                                     
c *************************************************************************

      implicit real*8(a-h,o-z)                                        
      parameter (natmax=200,ndmax=500)                                                    
      integer sorb
      character ac*4,symbol*2,keywrd*320,sline*80
      character segname*4,resname*3
      common/ato/ac(natmax),symbol(40)                                          
      common/atom1/ns(40),np(40),nd(40),nf(40),inatom,                          
     *isym(natmax),nelem                                                        
      common/atom2/oexp(40,8),cdf(40,8),couls(40),coulp(40),could(40),          
     *coulf(40),x(natmax),y(natmax),z(natmax),ires(natmax),
     *resname(natmax),segname(natmax)                                  
      common/keywrd/keywrd,sline                                                
      common/optv/rho,fermi,wstart,wend,
     *lb,le,lbr,ler,iovlb,iovle,icharg,nelec,
     *distl,idistn,ieb,iee,isb,ise,ihb,ihe,ipope,jpope,irest
      common/sij/s(ndmax,ndmax),h(ndmax,ndmax),c(ndmax),
     *maxs(natmax),maxp(natmax),maxd(natmax),sorb(2*natmax),
     *sp(natmax),pd(natmax)
      
      ih=1                                                                      
c                                                                               
c Hii                                                                           
c                                                                               
      do 70 i=1,natom                                                           
      keyi=isym(i)                                                              
      if(ns(keyi).eq.0) goto 10                                                 
      c(ih)=couls(keyi)                                                         
      ih=ih+1                                                                   
   10 if(np(keyi).eq.0) goto 30                                                 

      hh=coulp(keyi)                                                            
      do 20 jh=1,3                                                              
      c(ih)=hh                                                                  
      ih=ih+1                                                                   
   20 continue                                                                  
   30 if(nd(keyi).eq.0) goto 50                                                 

      hh=could(keyi)                                                            
      do 40 jh=1,5                                                              
      c(ih)=hh                                                                  
      ih=ih+1                                                                   
   40 continue                                                                  
   50 if(nf(keyi).eq.0) goto 70                                                 

      hh=coulf(keyi)                                                            
      do 60 jh=1,7                                                              
      c(ih)=hh                                                                  
      ih=ih+1                                                                   
   60 continue                                                                  
   70 continue                                                                  
c                                                                               
c nhij                                                                          
c                                                                               
      nhij=1                                                                    
      if(index(keywrd,'UWHIJ').ne.0) nhij=0                                     
      if(nhij.gt.0) goto 80                                                     
      assign 110 to ibr1                                                        
      goto 90                                                                   
   80 continue                                                                  
      assign 100 to ibr1                                                        
   90 continue                                                                  

      con=0.875
      const=con
c                                                                               
c hij                                                                           
c                                                                               
      h(1,1)=c(1)                                                               
      do i=2,ndim                                                           
         h(i,i)=c(i)
         i1=i-1                                                                    
         do j=1,i1                                                             
            hh=c(i)+c(j)                                                              
            goto ibr1,(110,100)                                                       
  100       et=(c(i)-c(j))/hh                                                         
            et=et*et                                                                  
            const=con+et/2.0+et*et*(0.5-con)                                          
  110       h(i,j)=const*hh*s(i,j)                                                    
            h(j,i)=h(i,j)                                                             
         enddo                                                                  
      enddo                                                                  

      return
      end