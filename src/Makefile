TARGET = \
		kep.o \
		kep_lanczos.o \
		kep_silanczos.o \
		kep_utility.o

all: $(TARGET)

include ../Makefile.inc

%.o : %.f90
	$(FC) $(FFLAGS) $(INCLUDES) -c $< -o $@

.PHONY: clean
clean:
	rm -f $(TARGET) mod_kep.mod
