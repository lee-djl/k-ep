# Fortran compiler and its options
FC = mpif90
FFLAGS = -O3 -qopenmp -g -traceback

# Archiver to generate a static library
ARCHIVER = ar

# MUMPS
MUMPSPATH = /set_the_path_to/MUMPS_5.1.2/lib
MUMPSLIB = -L$(MUMPSPATH) -ldmumps -lmumps_common

# METIS
METISPATH = /set_the_path_to/metis-5.1.0/lib
METISLIB = -L$(METISPATH) -lmetis

# PORD (MUMPS built-in library)
PORDPATH = /set_the_path_to/MUMPS_5.1.2/PORD/lib
PORDLIB = -L$(PORDPATH) -lpord

# MPI (dummy MPI for sequential MUMPS)
# MPIPATH = /set_the_path_to/MUMPS_5.1.2/libseq
# MPILIB = -L$(MPIPATH) -lmpiseq 

# MPI (for parallel MUMPS)
MPIPATH = 
MPILIB = -lmpi 

# ScaLAPACK and BLACS (for parallel MUMPS)
SCALAPACKPATH = 
SCALAPACKLIB = -lmkl_scalapack_lp64 
BLACSPATH = 
BLACSLIB = -lmkl_blacs_intelmpi_lp64

# LAPACK and BLAS
LAPACKPATH = 
LAPACKLIB = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
BLASPATH =
BLASLIB = 

# Others
OTHERLIB = -lpthread

LIBS =  $(MUMPSLIB) $(METISLIB) $(PORDLIB) \
        $(MPILIB) $(SCALAPACKLIB) $(BLACSLIB) \
        $(LAPACKLIB) $(BLASLIB) $(OTHERLIB)

# Search paths (for MUMPS)
# MUMPSMPI = -I/set_the_path_to/MUMPS_5.1.2/libseq
MUMPSINC = -I/set_the_path_to/MUMPS_5.1.2/include
INCLUDES = $(MUMPSINC)
