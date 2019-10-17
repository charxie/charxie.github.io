
      subroutine extrct(c,eval,borb,occ,nb,irnk)

c***********************************************************************
c fetch successive eigenvectors, checking for double occupancy, if
c successful, the prospective bond orbital BORB is given an appropriate
c LABEL
c***********************************************************************

      implicit real*8 (a-h,o-z)
      dimension c(8,8),eval(8),borb(8),ifnd(8)

      if(irnk.gt.nb) stop ' extrct error'
      bocc=3.0
      occ=0.0
      do ir = 1 , irnk
         ifnd(ir)=0
         if(ir.gt.1) bocc=occ
         occ=0.0
         do 20 i = 1 , nb
            if(eval(i).lt.occ) goto 20
            if(eval(i).gt.bocc) goto 20
            do k = 1 , ir
	       if(ifnd(k).eq.i) goto 20
	    enddo
	    occ=eval(i)
	    iocc=i
   20    continue
         ifnd(ir)=iocc
      enddo
      do j = 1 , nb
         borb(j)=c(j,iocc)
      enddo

      return
      end


      subroutine deplet(borb,occ,nb,iat1,iat2)

c***********************************************************************
c deplete the density matrix of the contribution from a lone pair
c***********************************************************************

      implicit real*8 (a-h,o-z)
      parameter(natmax=200,ndmax=500,nsmp=10)
      integer ul
      dimension borb(8),iat(2)
      common/info1/ul(natmax),ll(natmax),nele,nocc
      common /inst/ ovlp(ndmax,ndmax),hami(ndmax,ndmax),
     *              dens(ndmax,ndmax)

      iat(1)=iat1
      iat(2)=iat2
      nrow=0
      ncol=0

      do 40 i = 1 , 2
         ia=iat(i)
         if(ia.eq.0) goto 40
         iu=ul(ia)
         il=ll(ia)
         do irow = il , iu
            nrow=nrow+1
            ncol=0
            do 20 j = 1 , 2
               ja=iat(j)
               if(ja.eq.0) goto 20
               ju=ul(ja)
               jl=ll(ja)
	       do icol = jl , ju
	         ncol=ncol+1
	         dens(irow,icol)=dens(irow,icol)-occ*borb(nrow)*borb(ncol)
	       enddo
   20       continue
         enddo
   40 continue
      nb=nrow

      return
      end

      