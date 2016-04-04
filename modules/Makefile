include ../make.sys

objects = constants.o para.o input.o bxsfdata.o wanndata.o

mod.a: $(objects)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean:
	rm -rf *.o *.mod *.a
