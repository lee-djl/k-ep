TARGET = \
		kep.o \
		kep_lanczos.o \
		kep_silanczos.o \
		kep_utility.o

all: $(TARGET)

include ../Makefile.inc

%.o : %.f
	$(FC) $(FFLAGS) $(INCLUDE) -c $< 
%.o : %.f90
	$(FC) $(FFLAGS) $(INCLUDE) -c $< 

.PHONY: clean
clean:
	rm -f $(TARGET)
