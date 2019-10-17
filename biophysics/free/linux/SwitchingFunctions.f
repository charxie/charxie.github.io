        real*8 function switchingFunction(ntype,ts,t)
        implicit real*8 (a-h,o-z)
        
        offset=0.01d0*ts
        tau=t/ts
        if(ntype.eq.1) switchingFunction=(t+offset)/ts
        if(ntype.eq.2) switchingFunction=1.d0-t/ts+offset/ts
        if(ntype.eq.3) switchingFunction=1.01d0-tau**5*
     *                 (70*tau**4-315*tau**3+540*tau**2-420*tau+126)
        if(ntype.eq.4) switchingFunction=(1.d0-tau)**5+0.01d0
        
        return
        end

        real*8 function switchingFunctionDerivative(ntype,ts,t)
        implicit real*8 (a-h,o-z)
        
        tau=t/ts
        if(ntype.eq.1) switchingFunctionDerivative=1.d0/ts
        if(ntype.eq.2) switchingFunctionDerivative=-1.d0/ts
        if(ntype.eq.3) switchingFunctionDerivative=1.d0/ts*
     *                 (-9*70*tau**8+8*315*tau**7-7*540*tau**6
     *                 +6*420*tau**5-5*126*tau**4)
        if(ntype.eq.4) switchingFunctionDerivative=-5.d0/ts*
     *                 (1.d0-tau)**4
        
        return
        end
