c delete an orbital and reform hamiltonian

	subroutine deleteh(index,n,hold,hnew)
	implicit real*8 (a-h,o-z)
 	dimension hold(n,n),hnew(n-1,n-1)
	
	do j = 1 , n
	  do i = 1 , n
	     if(i.lt.index.and.j.lt.index) hnew(i,j)=hold(i,j)
	     if(i.gt.index.and.j.lt.index) hnew(i-1,j)=hold(i,j)
	     if(i.lt.index.and.j.gt.index) hnew(i,j-1)=hold(i,j)
	     if(i.gt.index.and.j.gt.index) hnew(i-1,j-1)=hold(i,j)
	  enddo
	enddo
        
        return
        end        


c delete orbitals and reform couplings

	subroutine deletec(idon,imid,iacc,n,hold,cbr)
	implicit real*8 (a-h,o-z)
 	dimension hold(n,n),cbr(3,n-3)
	
        do i = 1 , n
          if(i.lt.idon) cbr(1,i)=hold(i,idon)
          if(i.gt.idon.and.i.lt.imid) cbr(1,i-1)=hold(i,idon)
          if(i.gt.imid.and.i.lt.iacc) cbr(1,i-2)=hold(i,idon)
          if(i.gt.iacc) cbr(1,i-3)=hold(i,idon)
        enddo
        
        do i = 1 , n
          if(i.lt.idon) cbr(2,i)=hold(i,imid)
          if(i.gt.idon.and.i.lt.imid) cbr(2,i-1)=hold(i,imid)
          if(i.gt.imid.and.i.lt.iacc) cbr(2,i-2)=hold(i,imid)
          if(i.gt.iacc) cbr(2,i-3)=hold(i,imid)
        enddo
        
        do i = 1 , n
          if(i.lt.idon) cbr(3,i)=hold(i,iacc)
          if(i.gt.idon.and.i.lt.imid) cbr(3,i-1)=hold(i,iacc)
          if(i.gt.imid.and.i.lt.iacc) cbr(3,i-2)=hold(i,iacc)
          if(i.gt.iacc) cbr(3,i-3)=hold(i,iacc)
        enddo

        
        return
        end        
                                