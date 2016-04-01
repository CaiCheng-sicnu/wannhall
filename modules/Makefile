include ../make.sys

objects = constants.o para.o wanndata.o input.o bxsfdata.o

mod.a: $(objects)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean:
	rm -rf *.o *.mod *.a
