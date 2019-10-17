      implicit real*8(a-h,o-z)
      parameter(ndim=994)
      character filename*10
      dimension tmd(ndim), tdat(ndim), hdat(ndim), hdat1(ndim)
      dimension gamma(ndim), prob(ndim)
      character file1*21, file2*21
      
      file1(1:21)='laze/tdat.10-50.eta=1'
      file2(1:21)='laze/laze.10-50.eta=1'
      
      open(10,file=file1,status='old')
      open(30,file=file2)
      
      filename(1:5)='fock.'
      
      do n = 1 , ndim

         read(10,'(3f20.10)') time, tda, hddMinushaa

         tmd(n)=time
         tdat(n)=tda
         hdat(n)=hddMinushaa
         
         goto 1010

         if(n.lt.10) then
            filename(6:9)='0000'
            filename(10:10)=char(48+n)
         elseif(n.lt.100) then
            nn1=mod(n,10)
            nn2=int(n/10)
            filename(6:8)='000'
            filename(9:9)=char(48+nn2)
            filename(10:10)=char(48+nn1)
         elseif(n.lt.1000) then
            nn1=int(n/100)
            nnx=mod(n,100)
            nn2=int(nnx/10)
            nn3=mod(nnx,10)
            filename(6:7)='00'
            filename(8:8)=char(48+nn1)
            filename(9:9)=char(48+nn2)
            filename(10:10)=char(48+nn3)
         elseif(n.lt.10000) then
            nn1=int(n/1000)
            nnx=mod(n,1000)
            nn2=int(nnx/100)
            nny=mod(nnx,100)
            nn3=int(nny/10)
            nn4=mod(nny,10)
            filename(6:6)='0'
            filename(7:7)=char(48+nn1)
            filename(8:8)=char(48+nn2)
            filename(9:9)=char(48+nn3)
            filename(10:10)=char(48+nn4)
         elseif(n.lt.100000) then
            nn1=int(n/10000)
            nnx=mod(n,10000)
            nn2=int(nnx/1000)
            nny=mod(nnx,1000)
            nn3=int(nny/100)
            nnz=mod(nny,100)
            nn4=int(nnz/10)
            nn5=mod(nnz,10)
            filename(6:6)=char(48+nn1)
            filename(7:7)=char(48+nn2)
            filename(8:8)=char(48+nn3)
            filename(9:9)=char(48+nn4)
            filename(10:10)=char(48+nn5)
         else
            STOP ' too many time steps ! '
         endif
         
         open(20,file=filename)
         write(20,'(2d10.4)') hddMinushaa, tda
         write(20,'(2d10.4)') tda        , 0
         close(20)
         
1010     continue         

      enddo
      
      call deriv1(tmd,hdat,hdat1,ndim)
      

      write(30,'(a1,5x,a21)') 'A',file2
      
      write(30,'(2a)') '    T*       TDA       dHDA/dT     LZTIME',
     :'     RABITIME     LZ/RABI     LZ PROB'
      write(30,'(2a)') '  (fs)       (eV)      (eV/fs)      (fs)',
     :'        (fs)       (none)       (none)'
      
      aver_tda=0.0
      aver_hda=0.0
      aver_zen=0.0
      aver_rab=0.0
      aver_rat=0.0
      aver_pro=0.0
      
      rmsd_tda=0.0
      rmsd_hda=0.0
      rmsd_zen=0.0
      rmsd_rab=0.0
      rmsd_rat=0.0
      rmsd_pro=0.0
      
      num=0
      
      do i = 1 , ndim-1
         if(hdat(i)*hdat(i+1).lt.0.0) then
       	   num=num+1
           rabi=6.626/1.6/dabs(tdat(i))
           zener=dabs(tdat(i)/hdat1(i))
	   prob(i)=1.0-dexp(-6.28**2*zener/rabi)
           write(30,'(i6,5e12.3,e12.3)') 
     :              tmd(i),tdat(i),hdat1(i),
     :              zener,rabi,zener/rabi,
     :              prob(i)
           aver_tda=aver_tda+dabs(tdat(i)*tdat(i))
           aver_hda=aver_hda+dabs(hdat1(i))
           aver_zen=aver_zen+dabs(zener)
           aver_rab=aver_rab+dabs(rabi)
           aver_rat=aver_rat+dabs(zener/rabi)
           aver_pro=aver_pro+dabs(prob(i))
         endif
      enddo
      aver_tda=dsqrt(aver_tda/num)
      aver_hda=aver_hda/num
      aver_zen=aver_zen/num
      aver_rab=aver_rab/num
      aver_rat=aver_rat/num
      aver_pro=aver_pro/num
      write(30,*)' -----------------------------------------------'
      write(30,'(i6,6e12.3)') 
     :      num,
     :      aver_tda,aver_hda,
     :      aver_zen,aver_rab,
     :      aver_rat,aver_pro
            


      do i = 1 , ndim-1
         if(hdat(i)*hdat(i+1).lt.0.0) then
       	   num=num+1
           rabi=6.626/1.6/dabs(tdat(i))
           zener=dabs(tdat(i)/hdat1(i))
           rmsd_tda=rmsd_tda+(dabs(tdat(i))-aver_tda)**2
           rmsd_hda=rmsd_hda+(dabs(hdat1(i))-aver_hda)**2
           rmsd_zen=rmsd_zen+(dabs(zener)-aver_zen)**2
           rmsd_rab=rmsd_rab+(dabs(rabi)-aver_rab)**2
           rmsd_rat=rmsd_rat+(dabs(zener/rabi)-aver_rat)**2
           rmsd_pro=rmsd_pro+(dabs(prob(i))-aver_pro)**2
         endif
      enddo
      rmsd_tda=dsqrt(rmsd_tda/num)
      rmsd_hda=dsqrt(rmsd_hda/num)
      rmsd_zen=dsqrt(rmsd_zen/num)
      rmsd_rab=dsqrt(rmsd_rab/num)
      rmsd_rat=dsqrt(rmsd_rat/num)
      rmsd_pro=dsqrt(rmsd_pro/num)
      write(30,'(6x,6e12.3)') 
     :      rmsd_tda,rmsd_hda,
     :      rmsd_zen,rmsd_rab,
     :      rmsd_rat,rmsd_pro



      end
      
      subroutine deriv1(x,y,d1y,n)
      implicit real*8 (a-h,o-z)
      dimension x(n),y(n)
      dimension d1y(n)
      h=x(2)-x(1)
      d1y(1)=(-y(3)+4.d0*y(2)-3.d0*y(1))/(2.d0*h)
      d1y(n)=(3.d0*y(n)-4.d0*y(n-1)+y(n-2))/(2.d0*h)
      do i=2,n-1
         d1y(i)=(y(i+1)-y(i-1))/(2.d0*h)
      enddo
      return
      end      