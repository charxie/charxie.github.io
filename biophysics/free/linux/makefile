FFLAGS = -fno-f2c -fno-backslash -O2 -fugly-complex -fno-globals -Wno-globals
DEFS   = -DGNU
FC := g77 $(FFLAGS) $(DEFS)
LD := g77 $(FFLAGS)

OBJS =  main.o LatticeBuilder.o VelocityControl.o gear.o force.o \
        order.o GroundState.o SwitchingFunctions.o util.o work.o \
        diage.o LennardJones.o DensityOfStates.o QuasiHarmonic.o \
        ForceLambda.o GearLambda.o


mdyn : $(OBJS)
	$(LD) -o mdyn $(OBJS)

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