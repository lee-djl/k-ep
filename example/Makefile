include ../Makefile.inc

# k2l 
K2LPATHLIB = ../lib
K2LLIB = -L$(K2LPATHLIB) -lk2l
K2LPATHINC = ../include
K2LINC = -I$(K2LPATHINC)

TARGET_OBJ = \
		mmio.o \
		example_utility.o

TARGET = \
		example.out

.PHONY: all 
all: obj main

.PHONY: obj 
obj: $(TARGET_OBJ)

mmio.o : mmio.f
	$(FC) $(FFLAGS) -c $< -o $@

example_utility.o : example_utility.f90 
	$(FC) $(FFLAGS) -c $< -o $@

.PHONY: main
main: $(TARGET)

example.out : example.f90
	$(FC) $(FFLAGS) $(INCLUDES) $(K2LINC) $< -o $@ $(TARGET_OBJ) $(K2LLIB) $(LIBS)

.PHONY: clean
clean:
	rm -f $(TARGET_OBJ) example_utility.mod $(TARGET) 
