FLAGS = -r8 -i4 -fast -O4

OBJ = main.o hamiltonian.o diage.o eigen.o output.o statics.o dynamics.o \
      initial.o util.o gear.o projection.o tda.o floquet.o \
      bspace.o fourier.o

et: 		$(OBJ) makefile
	f77 $(FFLAGS) -o et $(OBJ)

main.o:	        	 main.f makefile

hamiltonian.o:           hamiltonian.f makefile

diage.o:                 diage.f makefile

eigen.o:                 eigen.f makefile

output.o:                output.f makefile

statics.o:               statics.f makefile

dynamics.o:              dynamics.f makefile

initial.o:               initial.f makefile

util.o:                  util.f makefile

gear.o:                  gear.f makefile

projection.o:            projection.f makefile

tda.o:                   tda.f makefile

floquet.o:               floquet.f makefile

bspace.o:                bspace.f makefile

fourier.o:               fourier.f makefile