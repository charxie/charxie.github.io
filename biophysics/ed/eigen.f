c Subroutine to solve the Roothaan equation
c omcm   ...... eigenenergies
c zr,zi  ...... eigenvectors

      subroutine eigen(ellow,elup,alle,n,hami,ovlp)

      implicit real*8 (a-h,o-z)
      parameter (neigd=248)
      logical invs,alle,eonly
      common /eigv/ omcm(neigd),zr(neigd,neigd),zi(neigd,neigd),ne,nblw
      common /secl/ s(neigd,neigd),h(neigd,neigd),nmat
      dimension ee(neigd)
      dimension hami(n,n),ovlp(n,n)

      if(n.gt.neigd) stop ' Dimension exceeds limit in eigen.f ! '

c--->    initialization

      do j = 1 , n
         do i = 1 , n
            s(i,j)=0.0
            if(i.eq.j) s(i,j)=1.0
            h(i,j)=0.0
         enddo
      enddo

      do i = 1 , n
         do j = 1, i
            h(j,i)=0.0
            h(i,j)=hami(i,j)
            s(i,j)=ovlp(i,j)
         enddo
      enddo
      
      nmat=n
      invs=.false.
      eonly=.false.

      call diagm(ellow,elup,invs,alle,eonly)
      

      return
      end