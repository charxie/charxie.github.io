        
        real*8 function lennardJones(epsilon,sigma,r)
        implicit real*8 (a-h, o-z)

        sr=sigma/r
        sr6=sr**6
        sr12=sr**12
        lennardJones=4.0*epsilon*(sr12-sr6)

        return
        end
        
        real*8 function lennardJonesDerivative1(epsilon,sigma,r)
        implicit real*8 (a-h, o-z)

        sr=sigma/r
        sr6=sr**6
        sr12=sr**12
        eps24=24.d0*epsilon
        lennardJonesDerivative1=eps24*(-2.d0*sr12+sr6)/r

        return
        end
        
        real*8 function lennardJonesDerivative2(epsilon,sigma,r)
        implicit real*8 (a-h, o-z)

        sr=sigma/r
        sr6=sr**6
        sr12=sr**12
        eps24=24.d0*epsilon
        lennardJonesDerivative2=eps24*(26.d0*sr12-7.d0*sr6)/(r*r)

        return
        end                