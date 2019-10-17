FLAGS = -r8 -i4 -fast -O4

OBJ = main.o rmsdev.o watdis.o prodis.o outpos.o rmsres.o avexyz.o \
      rmswat.o dopair.o nonbnd.o watpos.o radist.o
      
dcd:	$(OBJ) makefile
	f77 $(FFLAGS) -o dcd $(OBJ) 

main.o:        main.f   makefile
rmsdev.o:      rmsdev.f makefile
watdis.o:      watdis.f makefile
prodis.o:      prodis.f makefile
outpos.o:      outpos.f makefile
rmsres.o:      rmsres.f makefile
avexyz.o:      avexyz.f makefile
rmswat.o:      rmswat.f makefile
dopair.o:      dopair.f makefile
nonbnd.o:      nonbnd.f makefile
watpos.o:      watpos.f makefile
radist.o:      radist.f makefile