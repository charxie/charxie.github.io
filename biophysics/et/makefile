FLAGS = -r8 -i4 -fast -O4

OBJ = main.o hamiltonian.o diage.o eigen.o output.o \
      initial.o util.o gear.o ruku.o projection.o tda.o \
      bspace.o reduction.o nsdyn.o fsdyn.o delete.o triple.o
et: 		$(OBJ) makefile
	f77 $(FFLAGS) -o et $(OBJ)

main.o:	        	 main.f makefile

hamiltonian.o:           hamiltonian.f makefile

diage.o:                 diage.f makefile

eigen.o:                 eigen.f makefile

output.o:                output.f makefile

initial.o:               initial.f makefile

util.o:                  util.f makefile

gear.o:                  gear.f makefile

ruku.o:                  ruku.f makefile

projection.o:            projection.f makefile

tda.o:                   tda.f makefile

bspace.o:                bspace.f makefile

reduction.o:             reduction.f makefile

nsdyn.o:                 nsdyn.f makefile

fsdyn.o:                 fsdyn.f makefile

delete.o:                delete.f makefile

triple.o:                triple.f makefile