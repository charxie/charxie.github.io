FLAGS =  -r8 -i4 -fast -O4

OBJ = main.o LatticeBuilder.o VelocityControl.o gear.o force.o \
      order.o GroundState.o SwitchingFunctions.o util.o work.o \
      diage.o LennardJones.o DensityOfStates.o QuasiHarmonic.o \
      ForceLambda.o GearLambda.o
      
md:	$(OBJ) makefile
	f77 $(FFLAGS) -o md $(OBJ) 

main.o:               main.f               makefile
LatticeBuilder.o:     LatticeBuilder.f     makefile
VelocityControl.o:    VelocityControl.f    makefile
gear.o:               gear.f               makefile
force.o:              force.f              makefile
order.o:              order.f              makefile
GroundState.o:        GroundState.f        makefile
SwitchingFunctions.o: SwitchingFunctions.f makefile
util.o:               util.f               makefile
work.o:               work.f               makefile
diage.o:              diage.f              makefile
LennardJones.o:       LennardJones.f       makefile
DensityOfStates.o:    DensityOfStates.f    makefile
QuasiHarmonic.o:      QuasiHarmonic.f      makefile
ForceLambda.o:        ForceLambda.f        makefile
GearLambda.o:         GearLambda.f         makefile