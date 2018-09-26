# Fortran compiler and its options
FC = 
FFLAGS = 

# Archiver to generate a static library
ARCHIVER = 

# MUMPS
MUMPSPATH = 
MUMPSLIB = -L$(MUMPSPATH) -ldmumps -lmumps_common

# METIS
METISPATH = 
METISLIB = -L$(METISPATH) -lmetis

# PORD (MUMPS built-in library)
PORDPATH = 
PORDLIB = -L$(PORDPATH) -lpord

# MPI (dummy MPI for sequential MUMPS)
MPIPATH = 
MPILIB = -L$(MPIPATH) -lmpiseq 

# MPI (for parallel MUMPS)
# MPIPATH = 
# MPILIB = -L$(MPIPATH) -lmpi 

# ScaLAPACK and BLACS (for parallel MUMPS)
# SCALAPACKPATH = 
# SCALAPACKLIB = -L$(SCALAPACKPATH) -lscalapack
# BLACSPATH = 
# BLACSLIB = -L$(BLACSPATH) -lblacs

# LAPACK and BLAS
LAPACKPATH = 
LAPACKLIB = -L$(LAPACKPATH) -llapack
BLASPATH =
BLASLIB = -L$(BLASPATH) -lblas

# Others
OTHERLIB = -lpthread

LIBS =  $(MUMPSLIB) $(METISLIB) $(PORDLIB) \
        $(MPILIB) $(SCALAPACKLIB) $(BLACSLIB) \
        $(LAPACKLIB) $(BLASLIB) $(OTHERLIB)

# Search paths (for MUMPS)
MUMPSMPI = -I
MUMPSINC = -I
INCLUDES = $(MUMPSMPI) $(MUMPSINC)
