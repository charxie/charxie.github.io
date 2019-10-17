      subroutine getpar(natom,ndimen,nstp)

c **********************************************************************
c read the atomic parameters
c the number of valence electrons, and the charge                               
c **********************************************************************

      implicit real*8 (a-h,o-z)                                        
      character ac*4,keywrd*320,symbol*2,sline*80,ac1*4,                        
     *symbd*2,elemt*2,tex1*11,text*2,elems*2,elem*7,textf*10,                   
     *tex2*11,line*80,indt*10,segname*4,resname*3
      parameter (natmax=200)                                                    
      common/ato/ac(natmax),symbol(40)                                          
      common/atom1/ns(40),np(40),nd(40),nf(40),inatom,                           
     *isym(natmax),nelem
      common/atom2/exps(40),exps2(40),expp(40),expp2(40),expd(40),              
     *expd2(40),expf(40),expf2(40),cs1(40),cs2(40),cp1(40),cp2(40),             
     *cd1(40),cd2(40),cf1(40),cf2(40),couls(40),coulp(40),could(40),            
     *coulf(40),x(natmax),y(natmax),z(natmax),ires(natmax),
     *resname(natmax),segname(natmax) 
      common/iounit/iin,iout,istart(40),iungf,iunhl
      common/keywrd/keywrd,sline                                                
      common/bdata/nsd(103),npd(103),ndd(103),nfd(103),                         
     1coulsd(103),expss(103),expsd1(103),expsd2(103),csd1(103),                 
     *csd2(103),                                                                
     2coulpd(103),expps(103),exppd1(103),exppd2(103),cpd1(103),                 
     *cpd2(103),                                                                
     3couldd(103),expds(103),expdd1(103),expdd2(103),cdd1(103),                 
     *cdd2(103),                                                                
     4coulfd(103),expfs(103),expfd1(103),expfd2(103),cfd1(103),                 
     *cfd2(103)                                                                 
      common/delet/ind(1000),ndim,ival,norbe(natmax),ndel                       
      common/bdata1/iveld(103)                                                  
      common/bdata2/symbd(103)                                                  
      common/optv/rho,fermi,wstart,wend,
     *lb,le,lbr,ler,iovlb,iovle,lcharg,nelec,
     *distl,idistn,ieb,iee,isb,ise,ibb,ibe,ipope,jpope,irest
      common/aoout/indt(1000,2),iatom(natmax),lorb(1000)
      common/temp1/norbef(natmax)
      dimension kelem(40)
      dimension idzeta(40,4),iparam(40),text(4),elems(40),                      
     *elem(40),textf(20),idummy(1000)                                           

      ndel=0
      do i = 1 , natom
         ac(i)(4:4)=ac(i)(3:3)
         ac(i)(3:3)=ac(i)(2:2)
         ac(i)(2:2)=ac(i)(1:1)
         ac(i)(1:1)=' '
      enddo
c                                                                               
c search of elements                                                            
c                                                                               
      nelem=1                                                                   
      ac1=ac(1)                                                                 
      symbol(1)=ac1(1:2)
      isym(1)=1                                                                 

      do i=2,natom                                                           
         ac1=ac(i)                                                                 
         ifound=0                                                                  
         do j=1,nelem                                                           
            if(ac1(1:2).eq.symbol(j)) then                                            
               isym(i)=j                                                                
               ifound=1                                                                 
               goto 20                                                                  
            end if                                                                    
         enddo
   20    if(ifound.eq.0) then                                                      
            if(nelem.eq.38) then                                                      
               sline='number of elements limited to 39'                                  
               call writel(sline)                                                        
               STOP
            else
               nelem=nelem+1                                                            
               symbol(nelem)=ac1(1:2)                                                   
               isym(i)=nelem                                                            
            end if                                                                    
         end if                                                                    
      enddo
c                                                                               
c DZETA cards                                                                   
c                                                                               
      ibl=0                                                                     
      do i=1,nelem                                                           
         do j=1,4                                                               
            idzeta(i,j)=0                                                             
         enddo
         iparam(i)=0                                                               
      enddo
     
   60 if(index(sline,'DZETA').eq.0) goto 90                                     
      if(ibl.eq.0) write(iout,1010)                                              
      ibl=1                                                                     
      tex1='DZETA   '                                                           
      is=index(sline,'S.')                                                      
      ip=index(sline,'P.')                                                      
      id=index(sline,'D.')                                                      
      if=index(sline,'F.')                                                      
      nn=1                                                                      
      if(is.ne.0) then                                                          
         nn=nn+1                                                                  
         text(1)=' s'                                                             
      else                                                                      
         text(1)=' '                                                              
      end if                                                                    
      if(ip.ne.0) then                                                          
         nn=nn+1                                                                  
         text(2)=' p'                                                             
      else                                                                      
         text(2)='  '                                                             
      end if                                                                    
      if(id.ne.0) then                                                          
         nn=nn+1                                                                  
         text(3)=' d'                                                             
      else                                                                      
         text(3)='  '                                                             
      end if                                                                    
      if(if.ne.0) then                                                          
         nn=nn+1                                                                  
         text(4)=' f'                                                             
      else                                                                      
         text(4)='  '                                                             
      end if                                                                    
      nn=nn+2                                                                   
      k=0                                                                       
      do 80 i=nn,40                                                             
         ipos=istart(i)                                                            
         if(ipos.eq.81) goto 80                                                    
         elemt=sline(ipos:ipos+1)                                                  
         if(elemt(2:2).eq.' ') then                                                
            elemt(2:2)=elemt(1:1)                                                    
            elemt(1:1)=' '                                                           
         end if                                                                    
         do j=1,nelem                                                           
            if(elemt.eq.symbol(j)) then                                               
               k=k+1                                                                     
               elems(k)=elemt                                                            
               if(is.ne.0) idzeta(j,1)=1                                                 
               if(ip.ne.0) idzeta(j,2)=1                                                 
               if(id.ne.0) idzeta(j,3)=1                                                 
               if(if.ne.0) idzeta(j,4)=1                                                 
            end if                                                                    
	 enddo
   80 continue
      write(iout,1000)tex1,(text(l),l=1,4),(elems(l),l=1,k)                     
 1000 format(2x,a8,4a2,' of ',20(1x,a2))                                        
      call readl                                                                
      goto 60                                                                   
c                                                                               
c STO cards: This allows the user to change STO parameters
c                                                                               
   90 itest=1                                                                   
      write(iout,1010)                                                          
 1010 format(' ')                                                               
  100 if(index(sline,'STO').eq.0) goto 220                                      
      ipos=istart(2)                                                            
      elemt=sline(ipos:ipos+1)                                                  
      if(elemt(2:2).eq.' ') then                                                
         elemt(2:2)=elemt(1:1)                                                    
         elemt(1:1)=' '                                                           
      end if                                                                    
      ifound=0                                                                  
      do 110 i=1,nelem                                                          
         if(elemt.eq.symbol(i)) then                                               
            ifound=1                                                                 
            numb=i                                                                   
            iparam(i)=1                                                              
         end if                                                                    
  110 continue                                                                  
      if(ifound.eq.0) numb=40                                                   
      if(ifound.eq.0) goto 121
      do 120 j=1,103                                                            
      if(symbol(numb).ne.symbd(j)) goto 120                                     
      numb1=j                                                                   
  120 continue
  121 continue
      ns(numb)=0                                                                
      np(numb)=0                                                                
      nd(numb)=0                                                                
      nf(numb)=0                                                                
      ii=2                                                                      
      if(index(sline,'S.').ne.0) then                                           
         ii=ii+1                                                                  
         ns(numb)=reada(sline,istart(ii),80)                                      
      end if                                                                    
      if(index(sline,'P.').ne.0) then                                           
         ii=ii+1                                                                  
         np(numb)=reada(sline,istart(ii),80)                                      
      end if                                                                    
      if(index(sline,'D.').ne.0) then                                           
         ii=ii+1                                                                  
         nd(numb)=reada(sline,istart(ii),80)                                      
      end if                                                                    
      if(index(sline,'F.').ne.0) then                                           
         ii=ii+1                                                                  
         nf(numb)=reada(sline,istart(ii),80)                                      
      end if                                                                    
  130 call readl                                                                

      if(ns(numb).eq.0) goto 150                                                
      if(ifound.eq.0) goto 140                                                  
      if(ns(numb).ne.nsd(numb1)) itest=0                                        
      couls(numb)=reada(sline,istart(1),80)                                     
      if(couls(numb).eq.0.0) itest=0                                            
      if(istart(2).le.80) exps(numb)=reada(sline,istart(2),80)                  
      if(exps(numb).eq.0.0) itest=0                                             
      if(istart(3).le.80) cs1(numb)=reada(sline,istart(3),80)                   
      if(cs1(numb).eq.0.0) cs1(numb)=1.0                                        
      if(cs1(numb).lt.0.9999) idzeta(numb,1)=1                                  
      if(istart(4).le.80) exps2(numb)=reada(sline,istart(4),80)                 
      if(istart(5).le.80) cs2(numb)=reada(sline,istart(5),80)                   
      if(cs2(numb).gt.0.0001) idzeta(numb,1)=1                                  
  140 call readl                                                                

  150 if(np(numb).eq.0) goto 170                                                
      if(ifound.eq.0) goto 160                                                  
      if(np(numb).ne.npd(numb1)) itest=0                                        
      coulp(numb)=reada(sline,istart(1),80)                                     
      if(coulp(numb).eq.0.0) itest=0                                            
      if(istart(2).le.80) expp(numb)=reada(sline,istart(2),80)                  
      if(expp(numb).eq.0.0) itest=0                                             
      if(istart(3).le.80) cp1(numb)=reada(sline,istart(3),80)                   
      if(cp1(numb).eq.0.0) cp1(numb)=1.0                                        
      if(cp1(numb).lt.0.9999) idzeta(numb,2)=1                                  
      if(istart(4).le.80) expp2(numb)=reada(sline,istart(4),80)                 
      if(istart(5).le.80) cp2(numb)=reada(sline,istart(5),80)                   
      if(cp2(numb).gt.0.0001) idzeta(numb,2)=1                                  
  160 call readl                                                                

  170 if(nd(numb).eq.0) goto 190                                                
      if(ifound.eq.0) goto 180                                                  
      if(nd(numb).ne.ndd(numb1)) itest=0                                        
      could(numb)=reada(sline,istart(1),80)                                     
      if(could(numb).eq.0.0) itest=0                                            
      if(istart(2).le.80) expd(numb)=reada(sline,istart(2),80)                  
      if(expd(numb).eq.0.0) itest=0                                             
      if(istart(3).le.80) cd1(numb)=reada(sline,istart(3),80)                   
      if(cd1(numb).eq.0.0) cd1(numb)=1.0                                        
      if(cd1(numb).lt.0.9999) idzeta(numb,3)=1                                  
      if(istart(4).le.80) expd2(numb)=reada(sline,istart(4),80)                 
      if(istart(5).le.80) cd2(numb)=reada(sline,istart(5),80)                   
      if(cd2(numb).gt.0.0001) idzeta(numb,3)=1                                  
  180 call readl                                                                

  190 if(nf(numb).eq.0) goto 210                                                
      if(ifound.eq.0) goto 200                                                  
      if(nf(numb).ne.nfd(numb1)) itest=0                                        
      coulf(numb)=reada(sline,istart(1),80)                                     
      if(coulf(numb).eq.0.0) itest=0                                            
      if(istart(2).le.80) expf(numb)=reada(sline,istart(2),80)                  
      if(expf(numb).eq.0.0) itest=0                                             
      if(istart(3).le.80) cf1(numb)=reada(sline,istart(3),80)                   
      if(cf1(numb).eq.0.0) cf1(numb)=1.0                                        
      if(cf1(numb).lt.0.9999) idzeta(numb,4)=1                                  
      if(istart(4).le.80) expf2(numb)=reada(sline,istart(4),80)                 
      if(istart(5).le.80) cf2(numb)=reada(sline,istart(5),80)                   
      if(cf2(numb).gt.0.0001) idzeta(numb,4)=1                                  
  200 call readl                                                                

  210 goto 100                                                                  
c                                                                               
c search in the data block file: this is default
c iparam(i)=1 --- the parameters for this atom is read from a STO card
c before, or is given in the following routine
c                                                                               
  220 do 320 i=1,nelem                                                          
      if(iparam(i).eq.1) goto 320                                               
      ifound=0                                                                  
      do 230 j=1,103                                                            
      if(symbol(i).eq.symbd(j)) goto 240                                        
  230 continue                                                                  
      goto 320                                                                  
  240 iparam(i)=1                                                               
      ns(i)=nsd(j)                                                              
      np(i)=npd(j)                                                              
      nd(i)=ndd(j)                                                              
      nf(i)=nfd(j)                                                              
      if(ns(i).eq.0.and.np(i).eq.0.and.                                         
     *nd(i).eq.0.and.nf(i).eq.0) itest=0                                        

      if(ns(i).eq.0.0) goto 260                                                 
      couls(i)=coulsd(j)                                                        
      if(couls(i).eq.0.0) itest=0                                               
      if(idzeta(i,1).gt.0.9999) goto 250                                        
      exps(i)=expss(j)                                                          
      if(exps(i).eq.0.0) then                                                   
       idzeta(i,1)=1                                                            
       write(iout,1020)symbol(i)                                                
 1020  format(2x,'single zeta not available for the s',                         
     *' orbital  of ',a2)                                                       
       goto 250                                                                 
      end if                                                                    
      exps2(i)=0.0                                                              
      cs1(i)=1.0                                                                
      cs2(i)=0.0                                                                
      goto 260                                                                  
  250 exps(i)=expsd1(j)                                                         
      if(exps(i).eq.0.0) itest=0                                                
      exps2(i)=expsd2(j)                                                        
      cs1(i)=csd1(j)                                                            
      if(cs1(i).eq.0.0) itest=0                                                 
      cs2(i)=csd2(j)                                                            

  260 if(np(i).eq.0.0) goto 280                                                 
      coulp(i)=coulpd(j)                                                        
      if(coulp(i).eq.0.0) itest=0                                               
      if(idzeta(i,2).gt.0.9999) goto 270                                        
      expp(i)=expps(j)                                                          
      if(expp(i).eq.0.0) then                                                   
       idzeta(i,2)=1                                                            
       write(iout,1030)symbol(i)                                                
 1030  format(2x,'single zeta not available for the p',                         
     *' orbitals of ',a2)                                                       
       goto 270                                                                 
      end if                                                                    
      expp2(i)=0.0                                                              
      cp1(i)=1.0                                                                
      cp2(i)=0.0                                                                
      goto 280                                                                  
  270 expp(i)=exppd1(j)                                                         
      if(expp(i).eq.0.0) itest=0                                                
      expp2(i)=exppd2(j)                                                        
      cp1(i)=cpd1(j)                                                            
      if(cp1(i).eq.0.0) itest=0                                                 
      cp2(i)=cpd2(j)                                                            

  280 if(nd(i).eq.0.0) goto 300                                                 
      could(i)=couldd(j)                                                        
      if(could(i).eq.0.0) itest=0                                               
      if(idzeta(i,3).gt.0.9999) goto 290                                        
      expd(i)=expds(j)                                                          
      if(expd(i).eq.0.0) then                                                   
       idzeta(i,3)=1                                                            
       write(iout,1040)symbol(i)                                                
 1040  format(2x,'single zeta not available for the d',                         
     *' orbitals of ',a2)                                                       
       goto 290                                                                 
      end if                                                                    
      expd2(i)=0.0                                                              
      cd1(i)=1.0                                                                
      cd2(i)=0.0                                                                
      goto 300                                                                  
  290 expd(i)=expdd1(j)                                                         
      if(expd(i).eq.0.0) itest=0                                                
      expd2(i)=expdd2(j)                                                        
      cd1(i)=cdd1(j)                                                            
      if(cd1(i).eq.0.0) itest=0                                                 
      cd2(i)=cdd2(j)                                                            

  300 if(nf(i).eq.0.0) goto 320                                                 
      coulf(i)=coulfd(j)                                                        
      if(coulf(i).eq.0.0) itest=0                                               
      if(idzeta(i,4).gt.0.9999) goto 310                                        
      expf(i)=expfs(j)                                                          
      if(expf(i).eq.0.0) then                                                   
       idzeta(i,4)=1                                                            
       write(iout,1050)symbol(i)                                                
 1050  format(2x,'single zeta not available for the f',                         
     *' orbitals of ',a2)                                                       
       goto 310                                                                 
      end if                                                                    
      expf2(i)=0.0                                                              
      cf1(i)=1.0                                                                
      cf2(i)=0.0                                                                
      goto 320                                                                  
  310 expf(i)=expfd1(j)                                                         
      if(expf(i).eq.0.0) itest=0                                                
      expf2(i)=expfd2(j)                                                        
      cf1(i)=cfd1(j)                                                            
      if(cf1(i).eq.0.0) itest=0                                                 
      cf2(i)=cfd2(j)                                                            

  320 continue                                                                  
c                                                                               
c normalization                                                                 
c                                                                               
      do 360 i=1,nelem                                                          
      if(cs2(i).eq.0.0) goto 330                                                
      s=(4.0*exps(i)*exps2(i)/(exps(i)+exps2(i))**2)**(ns(i)+0.5)               
      s=1.0/dsqrt(cs1(i)**2+cs2(i)**2+2.0*s*cs1(i)*cs2(i))                      
      cs1(i)=s*cs1(i)                                                           
      cs2(i)=s*cs2(i)                                                           
  330 if(cp2(i).eq.0.0) goto 340                                                
      s=(4.0*expp(i)*expp2(i)/(expp(i)+expp2(i))**2)**(np(i)+0.5)               
      s=1.0/dsqrt(cp1(i)**2+cp2(i)**2+2.0*s*cp1(i)*cp2(i))                      
      cp1(i)=s*cp1(i)                                                           
      cp2(i)=s*cp2(i)                                                           
  340 if(cd2(i).eq.0.0) goto 350                                                
      s=(4.0*expd(i)*expd2(i)/(expd(i)+expd2(i))**2)**(nd(i)+0.5)               
      s=1.0/dsqrt(cd1(i)**2+cd2(i)**2+2.0*s*cd1(i)*cd2(i))                      
      cd1(i)=s*cd1(i)                                                           
      cd2(i)=s*cd2(i)                                                           
  350 if(cf2(i).eq.0.0) goto 360                                                
      s=(4.0*expf(i)*expf2(i)/(expf(i)+expf2(i))**2)**(nf(i)+0.5)               
      s=1.0/dsqrt(cf1(i)**2+cf2(i)**2+2.0*s*cf1(i)*cf2(i))                      
      cf1(i)=s*cf1(i)                                                           
      cf2(i)=s*cf2(i)                                                           
  360 continue                                                                  
c                                                                               
c output                                                                        
c
      if(nstp.gt.1) go to 4444
      do 370 i=1,nelem                                                          
      write(iout,1060)i,symbol(i)                                               
 1060 format(2x,'element #',i2,'  : ',a2)                                       
      if(ns(i).ne.0) then                                                       
      if(idzeta(i,1).eq.0) then                                                 
        write(iout,1070) ns(i),couls(i),exps(i),cs1(i)                            
      else                                                                      
        write(iout,1080) ns(i),couls(i),exps(i),cs1(i),exps2(i),cs2(i)            
      end if                                                                    
      end if                                                                    
 1070 format(2x,i2,'s      Hii = ',f10.4,'  zeta1 = ',f10.5,                    
     *'  c1 = ',f10.5)                                                          
 1080 format(2x,i2,'s      Hii = ',f10.4,'  zeta1 = ',f10.5,                    
     *'  c1 = ',f10.5,/,27x,'  zeta2 = ',f10.5,'  c2 = ',f10.5)                 
      if(np(i).ne.0) then                                                       
      if(idzeta(i,2).eq.0) then                                                 
        write(iout,1090) np(i),coulp(i),expp(i),cp1(i)                            
      else                                                                      
        write(iout,1100) np(i),coulp(i),expp(i),cp1(i),expp2(i),cp2(i)            
      end if                                                                    
      end if                                                                    
 1090 format(2x,i2,'p      Hii = ',f10.4,'  zeta1 = ',f10.5,                    
     *'  c1 = ',f10.5)                                                          
 1100 format(2x,i2,'p      Hii = ',f10.4,'  zeta1 = ',f10.5,                    
     *'  c1 = ',f10.5,/,27x,'  zeta2 = ',f10.5,'  c2 = ',f10.5)                 
      if(nd(i).ne.0) then                                                       
      if(idzeta(i,3).eq.0) then                                                 
        write(iout,1110) nd(i),could(i),expd(i),cd1(i)                            
      else                                                                      
        write(iout,1120) nd(i),could(i),expd(i),cd1(i),expd2(i),cd2(i)            
      end if                                                                    
      end if                                                                    
 1110 format(2x,i2,'d      Hii = ',f10.4,'  zeta1 = ',f10.5,                    
     *'  c1 = ',f10.5)                                                          
 1120 format(2x,i2,'d      Hii = ',f10.4,'  zeta1 = ',f10.5,                    
     *'  c1 = ',f10.5,/27x,'  zeta2 = ',f10.5,'  c2 = ',f10.5)                  
      if(nf(i).ne.0) then                                                       
      if(idzeta(i,4).eq.0) then                                                 
        write(iout,1130) nf(i),coulf(i),expf(i),cf1(i)                            
      else                                                                      
        write(iout,1140) nf(i),coulf(i),expf(i),cf1(i),expf2(i),cf2(i)            
      end if                                                                    
      end if                                                                    
 1130 format(2x,i2,'f      Hii = ',f10.4,'  zeta1 = ',f10.5,                    
     *'  c1 = ',f10.5)                                                          
 1140 format(2x,i2,'f      Hii = ',f10.4,'  zeta1 = ',f10.5,                    
     *'  c1 = ',f10.5,/,27x,'  zeta2 = ',f10.5,'  c2 = ',f10.5)                 
  370 continue                                                                  
 4444 continue
 
      if(itest.eq.0) then                                                       
        line='suspicious STO parameters, if ok use the FORCE option'              
        call writel(line)                                                         
      end if                                                                    
c                                                                               
c type and number of orbitals
c
c iatom ...... label for an atom of species i
c lorb  ...... label for an orbital
c
      do k = 1 , nelem
         kelem(k)=0
         do i = 1 , natom
            if(ac(i)(1:2).eq.symbol(k)) then 
               kelem(k)=kelem(k)+1
               iatom(i)=kelem(k)
            endif
         enddo
      enddo
      ndim=0
      do i = 1 , natom
         do k = 1, nelem
         if(isym(i).eq.k) then
         if(ns(k).ne.0) then
            ndim=ndim+1
            lorb(ndim)=iatom(i)
         endif
         if(np(k).ne.0) then
            do ii=1,3
               ndim=ndim+1
               lorb(ndim)=iatom(i)
            enddo
         endif
         if(nd(k).ne.0) then
            do ii=1,5
               ndim=ndim+1
               lorb(ndim)=iatom(i)
            enddo
         endif
         if(nf(k).ne.0) then
            do ii=1,7
               ndim=ndim+1
               lorb(ndim)=iatom(i)
            enddo
         endif
         endif
         enddo                                             
      enddo
      
      ndim=0                                                                    
      do 390 i=1,natom                                                          
      n=0
      do 380 j=1,nelem                                                          
      if(isym(i).ne.j) goto 380                                                 
      if(ns(j).ne.0) then                                                       
         if(ndim.ge.1000) goto 400                                                 
         indt(ndim+1,1)=ac(i)                                                      
         indt(ndim+1,2)='s'                                                        
         ndim=ndim+1                                                               
         n=n+1                                                                     
      end if                                                                    
      if(np(j).ne.0) then                                                       
         if(ndim.ge.998) goto 400                                                  
         indt(ndim+1,1)=ac(i)                                                      
         indt(ndim+1,2)='px'                                                       
         indt(ndim+2,1)=ac(i)                                                      
         indt(ndim+2,2)='py'                                                       
         indt(ndim+3,1)=ac(i)                                                      
         indt(ndim+3,2)='pz'                                                       
         ndim=ndim+3                                                               
         n=n+3                                                                     
      end if                                                                    
      if(nd(j).ne.0) then                                                       
         if(ndim.ge.996) goto 400                                                  
         indt(ndim+1,1)=ac(i)                                                      
         indt(ndim+1,2)='dx2-y2'                                                   
         indt(ndim+2,1)=ac(i)                                                      
         indt(ndim+2,2)='dz2'                                                      
         indt(ndim+3,1)=ac(i)                                                      
         indt(ndim+3,2)='dxy'                                                      
         indt(ndim+4,1)=ac(i)                                                      
         indt(ndim+4,2)='dxz'                                                      
         indt(ndim+5,1)=ac(i)                                                      
         indt(ndim+5,2)='dyz'                                                      
         ndim=ndim+5                                                               
         n=n+5                                                                     
      end if                                                                    
      if(nf(j).ne.0) then                                                       
         if(ndim.ge.994) goto 400                                                  
         indt(ndim+1,1)=ac(i)                                                      
         indt(ndim+1,2)='fz3'                                                      
         indt(ndim+2,1)=ac(i)                                                      
         indt(ndim+2,2)='fxz2'                                                     
         indt(ndim+3,1)=ac(i)                                                      
         indt(ndim+3,2)='fyz2'                                                     
         indt(ndim+4,1)=ac(i)                                                      
         indt(ndim+4,2)='fxyz'                                                     
         indt(ndim+5,1)=ac(i)                                                      
         indt(ndim+5,2)='fz(x2-y2)'                                                
         indt(ndim+6,1)=ac(i)                                                      
         indt(ndim+6,2)='fx(x2-3y2)'                                               
         indt(ndim+7,1)=ac(i)                                                      
         indt(ndim+7,2)='fy(3x2-y2)'                                               
         ndim=ndim+7                                                               
         n=n+7                                                                     
      end if                                                                    
      norbef(i)=n                                                               
      norbe(i)=n                                                                
  380 continue                                                                  
  390 continue                                                                  
      
      ndimen=ndim                                                                         
      write(iout,1150)ndim                                                      
 1150 format(1x,'Number of orbitals before restriction = ',i4)                     
      goto 410                                                                  
                                                                               
  400 line='number of atomic orbitals limited to 1000'                          
      call writel(line)                                                         
       STOP
c
c deletion                                                                      
c                                                                               
  410 do 420 i=1,ndim                                                           
      ind(i)=0                                                                  
      idummy(i)=0                                                               
  420 continue                                                                  

      if(index(sline,'AODELETE').eq.0) goto 460                                 
  430 if(index(sline,'AODELETE').eq.0) goto 440                                 
      call search(textf,nn,elem,kk,norbef,idummy,norbe,ind)                     
      call readl                                                                
      textf(1)=' AODELETE '                                                     
      tex1='   '                                                                
      tex2='     of    '                                                        
      write(iout,1160)textf(1),tex1,(textf(l),l=2,nn),tex2,                     
     *(elem(l),l=1,kk)                                                          
 1160 format(2x,130a)                                                           
      goto 430                                                                  

  440 do 450 i=1,ndim                                                           
      if(ind(i).eq.1) ndel=ndel+1                                               
  450 continue                                                                  
      write(iout,1170)ndel                                                      
 1170 format(2x,'# of orbitals deleted = ',i3)                                  

  460 ikeep=0                                                                   
  470 if(index(sline,'AOKEEP').eq.0) goto 480                                   
      call search(textf,nn,elem,kk,norbef,idummy,norbe,ind)                     
      call readl                                                                
      textf(1)=' AOKEEP'                                                        
      tex1='   '                                                                
      tex2='     of    '                                                        
      write(iout,1160)textf(1),tex1,(textf(l),l=2,nn),tex2,                     
     *(elem(l),l=1,kk)                                                          
      ikeep=1                                                                   
      goto 470                                                                  

  480 if(ikeep.eq.0) goto 500                                                   
      nkeep=0                                                                   
      do 490 i=1,ndim                                                           
      if(ind(i).eq.0) then                                                      
         ind(i)=1                                                                 
      else                                                                      
         ind(i)=0                                                                 
         nkeep=nkeep+1                                                            
      end if                                                                    
  490 continue                                                                  
      write(iout,1180)nkeep                                                     
 1180 format(2x,'# of orbitals kept = ',i3)                                     
      ndel=ndim-nkeep                                                           
                                                                                
c                                                                               
c modification of norbe and indt                                                
c                                                                               
  500 j1=1                                                                      
      ii=0                                                                      
      do 520 i=1,natom                                                          
      j2=j1+norbe(i)-1                                                          
      j11=j1                                                                    
      do 510 j=j11,j2                                                           
      j1=j1+1                                                                   
      if(ind(j).eq.1) then                                                      
       norbe(i)=norbe(i)-1                                                      
      else                                                                      
       ii=ii+1                                                                  
       indt(ii,1)=indt(j,1)                                                     
       indt(ii,2)=indt(j,2)                                                     
      end if                                                                    
  510 continue                                                                  
  520 continue                                                                  

c      if(ndel.ne.0) then                                                        
c         write(iout,1190)ii                                                       
c 1190    format(2x,'# of orbitals after restriction  = ',i4)                       
c      else                                                                      
c         write(iout,1200)                                                         
c 1200    format(2x,'no restriction')                                               
c      end if                                                                    
c                                                                               
c # of valence electrons per unit cell and charge                               
c                                                                               
      icharg=0                                                                  
      ii1=index(keywrd,'CHARGE=')                                               
      if(ii1.ne.0) icharg=lcharg                                                

      ival=nelec                                                                
      ii2=index(keywrd,'NELEC=')                                                

      ivalc=0                                                                   
      do 560 k=1,nelem                                                          
      do 530 l=1,103                                                            
      if(symbol(k).eq.symbd(l)) goto 540                                        
  530 continue                                                                  
  540 do 550 m=1,natom                                                          
      ac1=ac(m)                                                                 
      if(ac1(1:2).eq.symbol(k)) ivalc=ivalc+iveld(l)                            
  550 continue                                                                  
  560 continue                                                                  
      ivalc=ivalc-icharg                                                        

      if(ii2.eq.0) nelec=ivalc                                                  
      if(ii2.eq.0) ival=ivalc                                                   

c  570 write(iout,1210)nelec                                                     
c 1210 format(/2x,'# of valence electrons per unit cell = ',i4)                  
      if(ii1.ne.0) write(iout,1220)icharg                                       
 1220 format(2x,'charge of the system = ',i4)                                   
      write(iout,1010)                                                          

      if(ii2.ne.0) then
         if(ival.ne.ivalc) then
           line='# of calculated electrons different'//
     *       ' from input, if ok use the FORCE option'
           call writel(line)
         end if
      end if

      return                                                                    
      end 

      subroutine search(text,nn,elem,kk,iorbf,indf,iorb,ixx)                    

c **********************************************************************
c to classify and assemble the set of atoms and that of atomic                  
c orbitals for various calculations                                             
c **********************************************************************

      implicit real*8(a-h,o-z)                                        
      character text*10,string*80,elemt*7,elem*7,ac1*4,ac*4,                    
     *symbol*2,sline*80                                                         
      dimension text(20),elem(40),ip(3),id(5),if(7)                             
      parameter (natmax=200)                                                    

      common/ato/ac(natmax),symbol(40)                                          
      common/atom1/ns(40),np(40),nd(40),nf(40),natom,                           
     *isym(natmax),nelem                                                        
      dimension iorbf(natmax),indf(1000),iorb(natmax),ixx(1000)                 

      call getorb(itd,isd,ipd,idd,ifd,ip,id,if,nn,text)                         
      kk=0                                                                      
      nnp1=nn+1                                                                 
      do 100 i=nnp1,40                                                          
c get the name                                                                  
      it=0                                                                      
      call getnam(i,string,it)                                                  
      if(it.eq.1) goto 100                                                      
      elemt=string(1:7)                                                         
      kk=kk+1                                                                   
      elem(kk)=elemt                                                            
      idorb=1                                                                   
      ijk=0                                                                     
      do 90 j=1,natom                                                           
      idorb=idorb+iorbf(j)                                                      
      ijk=ijk+iorb(j)                                                           
      ac1=ac(j)                                                                 
      if(string(1:2).eq.'XX') elemt=ac(j)                                       
      if(string(4:5).eq.'XX') elemt(4:5)=ac1(4:5)                               
      do 10 k=1,7                                                               
      if(k.gt.2) then                                                           
      if(elemt(k:k).ne.ac1(k:k).and.elemt(k:k).ne.' ') goto 90                  
      else                                                                      
      if(elemt(k:k).ne.ac1(k:k)) goto 90                                        
      end if                                                                    
   10 continue                                                                  
      idorb=idorb-iorbf(j)                                                      
      ijk=ijk-iorb(j)                                                           
      do 20 k=1,nelem                                                           
      if(symbol(k).ne.elemt(1:2)) goto 20                                       
      numb=k                                                                    
   20 continue                                                                  
      ia=iorbf(j)                                                               
      if(ns(numb).eq.0) goto 30                                                 
      if(indf(idorb).eq.0) then                                                 
       ijk=ijk+1                                                                
       if(itd.ne.0) ixx(ijk)=1                                                  
       if(isd.ne.0) ixx(ijk)=1                                                  
      end if                                                                    
      idorb=idorb+1                                                             
   30 if(np(numb).eq.0) goto 50                                                 
      do 40 ik=1,3                                                              
      if(indf(idorb).eq.0) then                                                 
      ijk=ijk+1                                                                 
      if(itd.ne.0) ixx(ijk)=1                                                   
      if(ipd.ne.0) ixx(ijk)=1                                                   
      if(ip(ik).ne.0) ixx(ijk)=1                                                
      end if                                                                    
      idorb=idorb+1                                                             
   40 continue                                                                  
   50 if(nd(numb).eq.0) goto 70                                                 
      do 60 ik=1,5                                                              
      if(indf(idorb).eq.0) then                                                 
      ijk=ijk+1                                                                 
      if(itd.ne.0) ixx(ijk)=1                                                   
      if(idd.ne.0) ixx(ijk)=1                                                   
      if(id(ik).ne.0) ixx(ijk)=1                                                
      end if                                                                    
      idorb=idorb+1                                                             
   60 continue                                                                  
   70 if(nf(numb).eq.0) goto 90                                                 
      do 80 ik=1,7                                                              
      if(indf(idorb).eq.0) then                                                 
      ijk=ijk+1                                                                 
      if(itd.ne.0) ixx(ijk)=1                                                   
      if(ifd.ne.0) ixx(ijk)=1                                                   
      if(if(ik).ne.0) ixx(ijk)=1                                                
      end if                                                                    
      idorb=idorb+1                                                             
   80 continue                                                                  
   90 continue                                                                  
  100 continue                                                                  

      nn=nn-1                                                                   

      return                                                                    
      end