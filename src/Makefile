include ../Makefile.inc

TARGET_OBJ = \
		k2l_iostruct.o \
		k2l_intstruct.o \
		xorshift128plus.o \
		k2l_utility.o \
		k2l_solver.o \
		k2l_lanczos.o \
		k2l_stage.o \
		k2l.o 

TARGET_MOD = \
		k2l_iotype.mod \
		k2l_inttype.mod \
		xorshift128plus.mod \
		k2l_utility.mod \
		k2l_solver.mod \
		k2l_lanczos_mod.mod \
		k2l_stage.mod \
		k2l_mod.mod 

TARGET_LIB = \
		libk2l.a

.PHONY: all
all: obj mod lib

.PHONY: obj
obj: $(TARGET_OBJ)

%.o: %.f90
	$(FC) $(FFLAGS) $(INCLUDES) -c $< -o $@

.PHONY: mod
mod:
	mkdir -p ../include
	cp ./$(TARGET_MOD) ../include/

.PHONY: lib
lib: $(TARGET_OBJ)
	mkdir -p ../lib
	$(ARCHIVER) r ../lib/$(TARGET_LIB) $^
	
.PHONY: clean
clean:
	rm -f $(TARGET_OBJ) $(TARGET_MOD)
	rm -rf ../lib ../include
