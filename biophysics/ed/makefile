FLAGS = -r8 -i4 -fast -O4

OBJ = main.o hamiltonian.o diage.o eigen.o output.o dynamics.o \
      initial.o util.o gear.o ruku.o denmat.o green.o bspace.o tsa.o \
      distmap.o faster.o kemat.o matmat.o kmatrix.o dyneigen.o \
      geareigen.o

et: 		$(OBJ) makefile
	f77 $(FFLAGS) -o et $(OBJ)

main.o:	        	 main.f makefile
hamiltonian.o:           hamiltonian.f makefile
diage.o:                 diage.f makefile
eigen.o:                 eigen.f makefile
output.o:                output.f makefile
dynamics.o:              dynamics.f makefile
initial.o:               initial.f makefile
util.o:                  util.f makefile
gear.o:                  gear.f makefile
ruku.o:                  ruku.f makefile
denmat.o:                denmat.f makefile
green.o:                 green.f makefile
bspace.o:                bspace.f makefile
tsa.o:                   tsa.f makefile
distmap.o:               distmap.f makefile
faster.o:                faster.f makefile
kemat.o:                 kemat.f makefile
matmat.o:                matmat.f makefile
kmatrix.o:               kmatrix.f makefile
dyneigen.o:              dyneigen.f makefile
geareigen.o:             geareigen.f makefile