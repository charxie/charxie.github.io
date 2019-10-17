c Subroutine to get eigenenergies and eigenvectors
c omcm   ...... eigenenergies
c zr,zi  ...... eigenvectors

      subroutine eigen(n,hami)

      implicit real*8 (a-h,o-z)
      parameter (neigd=2000)
      logical invs,alle,eonly
      common /eigv/ omcm(neigd),zr(neigd,neigd),zi(neigd,neigd),ne,nblw
      common /secl/ s(neigd,neigd),h(neigd,neigd),nmat
      dimension ee(neigd)
      dimension hami(n,n)

      do i = 1 , n
         do j = 1 , n
            s(i,j)=0.0
            if(i.eq.j) s(i,j)=1.0
            h(i,j)=0.0
         enddo
      enddo

      do i = 1 , n
         do j = 1, i
            h(j,i)=0.0
            h(i,j)=hami(i,j)
         enddo
      enddo
      
      nmat=n
      
c--> solve secular equation

      invs=.false.
      alle=.false.
      eonly=.false.
      ellow=-1000.d0
      elup=  1000.d0

      call diagm(ellow,elup,invs,alle,eonly)
      
      return
      end