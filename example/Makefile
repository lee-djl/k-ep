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
 		../src/kep_utility.o \
 		./mmio.o \
 		./example_utility.o \
 		-I./

mmio.o : mmio.f
	$(FC) $(FFLAGS) $(INCLUDE) -c $< 

example_utility.o : example_utility.f90 
	$(FC) $(FFLAGS) $(INCLUDE) -c $< 

example.out : example.f90
	$(FC) $(FFLAGS) $(INCLUDE) $< -o $@ $(OBJECTS) $(LIBS)

.PHONY: clean
clean:
	rm -f $(TARGET) mod_example_readmtx.mod
