      subroutine friction(natom,ndim,nstp,tmd,ft)            

c **********************************************************************
c to calculate electronic friction matrix
c **********************************************************************

      implicit real*8(a-h,o-z)
      parameter (natmax=200,ndmax=500,nsmp=10)
      logical sp,pd
      integer sorb,ul
      logical rite,check,polish
      common/sij/s(ndmax,ndmax),h(ndmax,ndmax),c(ndmax),
     *maxs(natmax),maxp(natmax),maxd(natmax),sorb(2*natmax),
     *sp(natmax),pd(natmax)
      common/info1/ul(natmax),ll(natmax),nele,nocc
      common/iounit/iin,iout,istart(40),iungf,iunhl
      real*4 tmd(nsmp),ft(ndmax,ndmax,nsmp)
      
      rite=.false.
      check=.false.
      polish=.true.
      delta=tmd(2)-tmd(1)
      
      if(check) then
         do ntime = 1 , nstp
            call difov(natom,ndim,ntime,ntime)
            do i = 1 , ndim
               do j = 1 , ndim
                  ft(j,i,ntime)=s(j,i)
               enddo
            enddo
         enddo
         goto 3000
      endif
      
      do ntime = 1 , nstp 
      
         do i = 1 , ndim
            do j = 1 , ndim
               ft(j,i,ntime)=0.0
            enddo
         enddo
         
         if(ntime.gt.1.and.ntime.lt.nstp) then
            call difov(natom,ndim,ntime,ntime+1)
            do i = 1 , ndim
               do j = 1 , ndim
                  ft(j,i,ntime)=ft(j,i,ntime)+s(j,i)
               enddo
            enddo
            call difov(natom,ndim,ntime,ntime-1)
            do i = 1 , ndim
               do j = 1 , ndim
                  ft(j,i,ntime)=ft(j,i,ntime)-s(j,i)
               enddo
            enddo
            do i = 1 , ndim
               do j = 1 , ndim
                  ft(j,i,ntime)=ft(j,i,ntime)/(2.0*delta)
               enddo
            enddo
         endif

         if(ntime.eq.1) then
            call difov(natom,ndim,1,3)
            do i = 1 , ndim
               do j = 1 , ndim
                  ft(j,i,ntime)=ft(j,i,ntime)-s(j,i)
               enddo
            enddo
            call difov(natom,ndim,1,2)
            do i = 1 , ndim
               do j = 1 , ndim
                  ft(j,i,ntime)=ft(j,i,ntime)+4.0*s(j,i)
               enddo
            enddo
            call difov(natom,ndim,1,1)
            do i = 1 , ndim
               do j = 1 , ndim
                  ft(j,i,ntime)=ft(j,i,ntime)-3.0*s(j,i)
               enddo
            enddo
            do i = 1 , ndim
               do j = 1 , ndim
                  ft(j,i,ntime)=ft(j,i,ntime)/(2.0*delta)
               enddo
            enddo
         endif

         if(ntime.eq.nstp) then
            call difov(natom,ndim,nstp,nstp)
            do i = 1 , ndim
               do j = 1 , ndim
                  ft(j,i,ntime)=ft(j,i,ntime)+3.0*s(j,i)
               enddo
            enddo
            call difov(natom,ndim,nstp,nstp-1)
            do i = 1 , ndim
               do j = 1 , ndim
                  ft(j,i,ntime)=ft(j,i,ntime)-4.0*s(j,i)
               enddo
            enddo
            call difov(natom,ndim,nstp,nstp-2)
            do i = 1 , ndim
               do j = 1 , ndim
                  ft(j,i,ntime)=ft(j,i,ntime)+s(j,i)
               enddo
            enddo
            do i = 1 , ndim
               do j = 1 , ndim
                  ft(j,i,ntime)=ft(j,i,ntime)/(2.0*delta)
               enddo
            enddo
         endif

      enddo
      
      if(polish) then
         call ulll(natom,ndim)
         do ntime = 1 , nstp
            do i = 1 , ndim
               ft(i,i,ntime)=0.0
            enddo
            do iat = 1 , natom
               if(ul(iat)-ll(iat).eq.3) then
                  ipx=ll(iat)+1
                  ipy=ll(iat)+2
                  ipz=ll(iat)+3
                  ft(ipx,ipy,ntime)=0.0
                  ft(ipy,ipz,ntime)=0.0
                  ft(ipz,ipx,ntime)=0.0
                  ft(ipx,ipz,ntime)=0.0
                  ft(ipz,ipy,ntime)=0.0
                  ft(ipy,ipx,ntime)=0.0
               endif
            enddo
         enddo
      endif
      
3000  if(rite) then
         do ntime = 1 , nstp
            write(iout,*) ntime
            write(iout,*)
            do i = 1 , ndim
               write(iout,'(500f10.8)')(ft(i,j,ntime),j=1,ndim)
            enddo
         enddo
         stop
      endif

      return
      end
      
      subroutine fribo(ndim,nstp,tmd,ft)
	
c********************************************************************
c calculate nonadiabatic term in BO basis
c <a name="fribo">
c********************************************************************	
	
      implicit real*8 (a-h,o-z)
      parameter(nsmp=10,natmax=200,ndmax=500,neigd=500)
      real*4 sao,trt
      real*4 ft(ndmax,ndmax,nsmp),tmd(nsmp)
      logical rite,polish
      character*1 transa,transb
      common/iounit/iin,iout,istart(40),iungf,iunhl
      common/nonadi/sao(ndmax,ndmax,nsmp),trt(ndmax,ndmax,nsmp)
      dimension cik(nstp),cikd1(nstp),tim(nstp)
      dimension td1(ndim,ndim,nstp)
      dimension tran(ndim,ndim),trat(ndim,ndim),trad(ndim,ndim)
      dimension over(ndim,ndim),fric(ndim,ndim),temp(ndim,ndim)

      alpha=1.d0
      beta=0.d0
      lda=ndim
      ldb=ndim
      ldc=ndim
      transa='N'
      transb='N'
        	
      rite=.false.
      polish=.true.
      if(mod(ndim,2).eq.0) then
         nhaf=ndim*0.5
      else
         nhaf=(ndim+1)*0.5
      endif

c>>> k,l...indices for BO; i,j...indices for AO

      do n = 1 , nstp
         tim(n)=tmd(n)
      enddo
      do k = 1 , ndim
         do i = 1 , ndim
            do n = 1 , nstp
               cik(n)=trt(i,k,n)
            enddo
            call deriv1(tim,cik,cikd1,nstp)
            do n = 1 , nstp
               td1(i,k,n)=cikd1(n)
            enddo
         enddo
      enddo
        
      dmax=0.0
      do n = 1 , nstp
         do k = 1 , ndim
            do i = 1 , ndim
               if(dabs(td1(i,k,n)).gt.dmax)then
                  dmax=dabs(td1(i,k,n))
                  nmax=n
                  kmax=k
                  imax=i
               endif
            enddo
         enddo
      enddo
      write(iout,*)'largest derivative of T-matrix element'
      write(iout,'(3i10,f10.5)') nmax,kmax,imax,dmax
      if(dmax.gt.2.d0)stop'check if the above derivative is too large'

      write(iout,*)
      write(iout,*)'transforming nonadiabatic term into BO basis...'

      do ntime = 1 , nstp
         write(iout,*) ntime
         do j = 1 , ndim
            do i = 1 , ndim
               tran(i,j)=trt(i,j,ntime)
               trat(i,j)=trt(j,i,ntime)
               trad(i,j)=td1(i,j,ntime)
               over(i,j)=sao(i,j,ntime)
               fric(i,j)= ft(i,j,ntime)
            enddo
         enddo
c         call mprod(ndim,trat,over,temp)
c         call mprod(ndim,temp,trad,over)
c         call mprod(ndim,trat,fric,temp)
c         call mprod(ndim,temp,tran,fric)
c         call strassen(ndim,nhaf,trat,over,temp)
c         call strassen(ndim,nhaf,temp,trad,over)
c         call strassen(ndim,nhaf,trat,fric,temp)
c         call strassen(ndim,nhaf,temp,tran,fric)
          call dgemm(transa,transb,ndim,ndim,ndim,alpha,
     *               trat,lda,over,ldb,beta,temp,ldc)         
          call dgemm(transa,transb,ndim,ndim,ndim,alpha,
     *               temp,lda,trad,ldb,beta,over,ldc)         
          call dgemm(transa,transb,ndim,ndim,ndim,alpha,
     *               trat,lda,fric,ldb,beta,temp,ldc)         
          call dgemm(transa,transb,ndim,ndim,ndim,alpha,
     *               temp,lda,tran,ldb,beta,fric,ldc)         

         do j = 1 , ndim
            do i = 1 , ndim
               ft(i,j,ntime)=fric(i,j)+over(i,j)
            enddo
            if(polish) ft(j,j,ntime)=0.0
         enddo
      enddo

      dmax=0.0
      do n = 1 , nstp
         do k = 1 , ndim
            do i = 1 , ndim
               if(abs(ft(i,k,n)).gt.dmax)then
                  dmax=abs(ft(i,k,n))
                  nmax=n
                  kmax=k
                  imax=i
               endif
            enddo
         enddo
      enddo
      write(iout,*)'largest nonadiabatic force'
      write(iout,'(3i10,f10.5)') nmax,kmax,imax,dmax
      if(dmax.gt.1.d0) stop 'check F-term'
      if(rite) then   
         do n = 1 , nstp
            write(iout,*) n
            do j = 1 , ndim
c               write(iout,'(500f10.5)')(ft(i,j,n),i=1,ndim)
               write(iout,'(500f10.5)')ft(j,j,n)
            enddo
            write(iout,*)
         enddo
         stop
      endif
     
      return
      end