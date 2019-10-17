      program aqvist
      implicit real(a-h,o-z)
      parameter(nmax=50)
      dimension vdwProHydan(nmax),eleProHydan(nmax)
      dimension vdwSolHydan(nmax),eleSolHydan(nmax)
      dimension vdwProNhydan(nmax),eleProNhydan(nmax)
      dimension vdwSolNhydan(nmax),eleSolNhydan(nmax)
      dimension vdwProMhydan(nmax),eleProMhydan(nmax)
      dimension vdwSolMhydan(nmax),eleSolMhydan(nmax)
      character endSign*6

c  read interaction energies
      
      open(11,file='proint.hydan')
      open(12,file='watint.hydan')
      open(21,file='proint.nhydan')
      open(22,file='watint.nhydan')
      open(31,file='proint.mhydan')
      open(32,file='watint.mhydan')

      ncom=0
      
      write(6,*) ' Ligands      No.    vdw       ele'
      write(6,*) '-------------------------------------'
      
      read(11,*)
      read(11,*)
      do i = 1 , nmax
         read(11,'(a6)') endSign
         if(endSign.eq.'  ----') then
            nend=i-1
            goto 110
         endif
      enddo
 110  rewind 11
      read(11,*)
      read(11,*)
      do i = 1 , nend
         read(11,'(16x,2f10.5)') vdwProHydan(i),eleProHydan(i)
      enddo
      close(11)
      vtotProHydan=0.0
      etotProHydan=0.0
      incl=nend-ncom
      do i = 1 , incl
         vtotProHydan=vtotProHydan+vdwProHydan(i)
         etotProHydan=etotProHydan+eleProHydan(i)
      enddo
      write(6,'(a11,i6,2f10.4)') 
     :      ' Hydan(Pro)',incl,vtotProHydan,etotProHydan
      read(12,*)
      read(12,*)
      do i = 1 , nend
         read(12,'(16x,2f10.5)') vdwSolHydan(i),eleSolHydan(i)
      enddo
      close(12)
      vtotSolHydan=0.0
      etotSolHydan=0.0
      do i = 1 , incl
         vtotSolHydan=vtotSolHydan+vdwSolHydan(i)
         etotSolHydan=etotSolHydan+eleSolHydan(i)
      enddo
      write(6,'(a11,i6,2f10.4)') 
     :      ' Hydan(Sol)',incl,vtotSolHydan,etotSolHydan

      
      read(21,*)
      read(21,*)
      do i = 1 , nmax
         read(21,'(a6)') endSign
         if(endSign.eq.'  ----') then
            nend=i-1
            goto 210
         endif
      enddo
 210  rewind 21
      read(21,*)
      read(21,*)
      do i = 1 , nend
         read(21,'(16x,2f10.5)') vdwProNhydan(i),eleProNhydan(i)
      enddo
      close(21)
      vtotProNhydan=0.0
      etotProNhydan=0.0
      incl=nend-ncom
      do i = 1 , incl
         vtotProNhydan=vtotProNhydan+vdwProNhydan(i)
         etotProNhydan=etotProNhydan+eleProNhydan(i)
      enddo
      write(6,'(a12,i5,2f10.4)') 
     :      ' Nhydan(Pro)',incl,vtotProNhydan,etotProNhydan
      read(22,*)
      read(22,*)
      do i = 1 , nend
         read(22,'(16x,2f10.5)') vdwSolNhydan(i),eleSolNhydan(i)
      enddo
      close(22)
      vtotSolNhydan=0.0
      etotSolNhydan=0.0
      do i = 1 , incl
         vtotSolNhydan=vtotSolNhydan+vdwSolNhydan(i)
         etotSolNhydan=etotSolNhydan+eleSolNhydan(i)
      enddo
      write(6,'(a12,i5,2f10.4)') 
     :      ' Nhydan(Sol)',incl,vtotSolNhydan,etotSolNhydan

      
      read(31,*)
      read(31,*)
      do i = 1 , nmax
         read(31,'(a6)') endSign
         if(endSign.eq.'  ----') then
            nend=i-1
            goto 310
         endif
      enddo
 310  rewind 31
      read(31,*)
      read(31,*)
      do i = 1 , nend
         read(31,'(16x,2f10.5)') vdwProMhydan(i),eleProMhydan(i)
      enddo
      close(31)
      vtotProMhydan=0.0
      etotProMhydan=0.0
      incl=nend-ncom
      do i = 1 , incl
         vtotProMhydan=vtotProMhydan+vdwProMhydan(i)
         etotProMhydan=etotProMhydan+eleProMhydan(i)
      enddo
      write(6,'(a12,i5,2f10.4)') 
     :      ' Mhydan(Pro)',incl,vtotProMhydan,etotProMhydan
      read(32,*)
      read(32,*)
      do i = 1 , nend
         read(32,'(16x,2f10.5)') vdwSolMhydan(i),eleSolMhydan(i)
      enddo
      close(32)
      vtotSolMhydan=0.0
      etotSolMhydan=0.0
      do i = 1 , incl
         vtotSolMhydan=vtotSolMhydan+vdwSolMhydan(i)
         etotSolMhydan=etotSolMhydan+eleSolMhydan(i)
      enddo
      write(6,'(a12,i5,2f10.4)') 
     :      ' Mhydan(Sol)',incl,vtotSolMhydan,etotSolMhydan

      write(6,*) '-------------------------------------'

c Aqvist      
      alpha=0.5
c      betaHydan=0.16
c      betaMhydan=0.16
c      betaNhydan=0.16
      
c NDR correlation
c      slope=(0.84-0.63)/(0.85-0.10)
c      betaHydan =(0.34-0.63)/slope+0.1
c      betaMhydan=(0.45-0.63)/slope+0.1
c      betaNhydan=(0.33-0.63)/slope+0.1
      
c WNDR correlation
      slope=(1.06-1.44)/(0.87-0.13)
      betaHydan =(6.17-0.13)/slope+1.44
      betaMhydan=(2.23-0.13)/slope+1.44
      betaNhydan=(8.54-0.13)/slope+1.44
      
      beta=betaHydan
      dGHydan=alpha*(etotProHydan-etotSolHydan)
     :        +beta*(vtotProHydan-vtotSolHydan)
      beta=betaNhydan
      dGNhydan=alpha*(etotProNhydan-etotSolNhydan)
     :         +beta*(vtotProNhydan-vtotSolNhydan)
      beta=betaMhydan
      dGMhydan=alpha*(etotProMhydan-etotSolMhydan)
     :         +beta*(vtotProMhydan-vtotSolMhydan)
      
      write(6,*) ' Aqvist alpha=',alpha,' beta=',beta
      write(6,*) '   dGHydan     dGNhydan     dGMhydan'
      write(6,*) '-------------------------------------'
      write(6,'(3f12.4)') dGHydan,dGNhydan,dGMhydan
      write(6,*) dGMhydan-dGHydan, dGNhydan-dGHydan

      open(99,file='curve.beta')

      do i = 1 , 100
      beta=-0.50+i*0.01
      dGHydan=alpha*(etotProHydan-etotSolHydan)
     :        +beta*(vtotProHydan-vtotSolHydan)
      dGNhydan=alpha*(etotProNhydan-etotSolNhydan)
     :         +beta*(vtotProNhydan-vtotSolNhydan)
      dGMhydan=alpha*(etotProMhydan-etotSolMhydan)
     :         +beta*(vtotProMhydan-vtotSolMhydan)
      write(99,*) beta, dGHydan+3.44, dGMhydan-0.11, dGNhydan+1.15
      enddo
      

      
      end