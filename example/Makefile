TARGET = \
		mmio.o \
		example_utility.o \
		example.out

all: $(TARGET)

include ../Makefile.inc

OBJECTS = \
		../src/kep.o \
 		../src/kep_lanczos.o  \
		../src/kep_silanczos.o \
 		../src/xorshift128plus.o \
 		../src/kep_utility.o \
 		./mmio.o \
 		./example_utility.o \
 		-I../src

mmio.o : mmio.f
	$(FC) $(FFLAGS) $(INCLUDES) -c $< -o $@

example_utility.o : example_utility.f90 
	$(FC) $(FFLAGS) $(INCLUDES) -c $< -o $@

example.out : example.f90
	$(FC) $(FFLAGS) $(INCLUDES) $< -o $@ $(OBJECTS) $(LIBS)

.PHONY: clean
clean:
	rm -f $(TARGET) mod_example_readmtx.mod
