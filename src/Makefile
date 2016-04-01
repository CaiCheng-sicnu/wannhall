include ../make.sys

LIBS = ../modules/mod.a

wannhall: main.f90 delta.o $(LIBS)
	$(F90) $(F90FLAGS) -fpp -o wannhall.x main.f90 delta.o $(LIBS) -I$(MKLINCLUDE) -I../modules -I. $(LAPACK95LIBS) $(MKLLIBS)

delta.o: delta.f90
	$(F90) $(F90FLAGS) -fpp -c delta.f90 -I$(MKLINCLUDE) -I../modules -I.

$(LIBS):
	cd ../modules
	make mod.a
	cd ..

clean:
	rm -rf *.o *.mod *.x
