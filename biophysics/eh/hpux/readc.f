      subroutine getpos(natom,nstp)
      
c **********************************************************************        
c to read atomic positions from the pdb files            
c **********************************************************************

      implicit real*8(A-H,O-Z)                                        
      character string*80
      character ac*4,keywrd*320,sline*80,symbol*2,segname*4,resname*3
      parameter (natmax=200)                                                    
      common/iounit/iin,iout,istart(40),iungf,iunhl                                               
      common/keywrd/keywrd,sline                                                
      common/ato/ac(natmax),symbol(40)                                          
      common/atom1/ns(40),np(40),nd(40),nf(40),inatom,                           
     *isym(natmax),nelem                                                        
      common/atom2/exps(40),exps2(40),expp(40),expp2(40),expd(40),              
     *expd2(40),expf(40),expf2(40),cs1(40),cs2(40),cp1(40),cp2(40),             
     *cd1(40),cd2(40),cf1(40),cf2(40),couls(40),coulp(40),could(40),            
     *coulf(40),x(natmax),y(natmax),z(natmax),ires(natmax),
     *resname(natmax),segname(natmax)

      iia=ichar('A')                                                            
      iiz=ichar('Z')                                                            
      nz=0
      call readl
      if(index(sline,'ATOM').eq.0) then                                          
         sline='card ATOM expected'                                                
         if(nstp.eq.1) call writel(sline)                                                       
      end if                                                                    

c>>> start reading x,y,z coordinates

  10  continue
      if(index(sline,'ATOM').eq.0) goto 50
      nz=nz+1                                                                   
      if(nz.gt.natmax) then                                                     
         sline='number of atoms limited to 200'                                    
         if(nstp.eq.1) call writel(sline)
         stop
      end if

      len=istart(4)-istart(3)                                                   
      if(len.gt.4) len=4                                                        
      string=sline(istart(3):istart(4)-1)                                       
      is=ichar(string(1:1))
      if(is.lt.iia.or.is.gt.iiz) then                                           
         do i=len,1,-1
            string(i:i)=string(i+1:i+1)                                               
   	 enddo
         string(len:len)=' '                                                           
      end if                                                                    
      if(len.lt.4) then                                                         
         do i=len+1,4                                                           
            string(i:i)=' '                                                           
   	 enddo
      end if

      ac(nz)      = string(1:4)
      resname(nz) = sline(istart(4):istart(4)+2)
      ires(nz)    = reada(sline,istart(5),80)
      x(nz)       = reada(sline,istart(6),80)                                     
      y(nz)       = reada(sline,istart(7),80)                                     
      z(nz)       = reada(sline,istart(8),80)
      segname(nz) = sline(istart(9):istart(9)+3)
      call readl                                                                
      goto 10                                                                   

   50 continue
      natom=nz

      write(iout,*) 'Number of atoms = ',natom
      if(nstp.gt.1) return
      
      do i = 1 , natom
         write(iout,'(2x,a1,4x,i3,6x,f12.7,2(3x,f12.7),5x,a3,2x,a4)')
     :        ac(i)(1:1),i,x(i),y(i),z(i),resname(i),segname(i)         
      enddo 

      return                                                                    
      end                                                                       

      subroutine getsel                                                         

c **********************************************************************                                                                               
c to read a title and all key words                                             
c **********************************************************************                                                                               

      implicit real*8(A-H, O-Z)                                       
      character keywrd*320,sline*80,title*150                                   
      character space*1,space2*2                           
      data space,space2/' ','  '/                                               

      common/iounit/iin,iout,istart(40),iungf,iunhl                                        
      common/keywrd/keywrd,sline                                                

      title=' '                                                                 
      call readl                                                                
      if(index(sline,'REMARK   TITLE').ne.0) then                                        
         title=sline                                                               
         title(1:15)='TRAJECTORY FR: '                                                        
         call readl                                                                
      end if                                                                    

      keywrd=' '                                                                
      ip2=1                                                                     
   10 if(index(sline,'REMARK  KEYWRD').eq.0) goto 40                                    
      do 20 i=80,istart(2),-1                                                   
         if(sline(i:i).ne.space) goto 30                                           
   20 continue                                                                  
   30 ip1=i+1                                                                   
      keywrd=keywrd(1:ip2)//sline(istart(2):i)//' '                             
      ip2=ip2+i+2-istart(2)                                                     
      if(ip2.ge.320) then                                                       
         sline='too many key words'                                               
         call writel(sline)                                                       
         return                                                                   
      end if                                                                    
      call readl                                                                
      goto 10                                                                   
   40 continue                                                                  
c                                                                               
c output                                                                        
c                                                                               
c      write(iout,1000) title                                                     
 1000 format(/' ',a//)                                                         

      return                                                                    
      end               
