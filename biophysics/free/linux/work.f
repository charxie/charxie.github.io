        subroutine work(ntot,delta,nprod,nsave,temperature,ntype,free0)
        implicit real*8 (a-h,o-z)
        parameter(kmax=100000)
        common/record/epot(kmax),ekin(kmax),etot(kmax),para(kmax)
        dimension workFunction(nprod)
        
        open(11,file='work')
        open(12,file='free')
        
        delta1=delta*nsave
        ts=nprod*delta1
        ratio0=switchingFunction(ntype,ts,delta1)
        do n = 1 , nprod
           workFunction(n)=0.d0
           ratio=switchingFunction(ntype,ts,n*delta1)
           do i = 1 , n
              t=i*delta1
              workFunction(n)=workFunction(n)+epot(i)
     *            *delta1*switchingFunctionDerivative(ntype,ts,t)
           enddo
           write(11,'(f10.5,f20.10)') ratio,workFunction(n)
           temp=temperature/ratio
           converter=1.38d0*1.5d0/1.6d0*0.0001d0
           free=(free0+workFunction(n))*ratio0/ratio
     *         +1.5d0*1.38d0/1.6d0*0.0001d0*dlog(ratio/ratio0)*temp
           write(12,'(f10.2,f20.10)') temp,free
        enddo
        
        close(11)
        close(12)

        return
        end
        
        
        subroutine entropy(ntot,nprod)
        implicit real*8 (a-h,o-z)
        dimension free(nprod),temp(nprod),entr(nprod)
        
        open(11,file='entropy')
        open(12,file='free')
        
        do n = 1 , nprod
           read(12,'(f10.2,f20.10)') temp(n),free(n)
        enddo
        
        call slope(temp,free,entr,nprod)
        
        do n = 1 , nprod
           write(11,'(f10.2,f20.10)') temp(n),-entr(n)
        enddo
        
        close(11)
        close(12)

        return
        end