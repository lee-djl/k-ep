# Fortran compiler and its options
FC = gfortran
FFLAGS = -fimplicit-none -O3

# MUMPS
MUMPSPATH = set_path_if_necessary_/MUMPS_5.0.1/lib
MUMPSFLAGS = -L$(MUMPSPATH) -ldmumps -lmumps_common 

# METIS
METISPATH = set_path_if_necessary_/metis-5.1.0/lib
METISFLAGS = -L$(METISPATH) -lmetis

# PORD (MUMPS built-in library)
PORDPATH = set_path_if_necessary_/MUMPS_5.0.1/PORD/lib
PORDFLAGS = -L$(PORDPATH) -lpord

# MPI for MUMPS (Sequential version of MUMPS requires MPI)
MPIPATH = set_path_if_necessary/MUMPS_5.0.1/libseq
MPIFLAGS = -L$(MPIPATH) -lmpiseq

# LAPACK
LAPACKPATH = set_path_if_necessary_/lib
LAPACKFLAGS = -L$(LAPACKPATH) -llapack

# BLAS
BLASPATH = set_path_if_necessary_/lib
BLASFLAGS = -L$(BLASPATH) -lblas

OTHERFLAGS = -lpthread -lm

LIBS = $(MUMPSFLAGS) $(METISFLAGS) $(PORDFLAGS) $(MPIFLAGS) $(LAPACKFLAGS) $(BLASFLAGS) $(OTHERFLAGS)

# For MUMPS
INCMPI = set_path_if_necessary_/MUMPS_5.0.1/libseq
INCMUMPS = set_path_if_necessary_/MUMPS_5.0.1/include
INCLUDES = -I. -I$(INCMPI) -I$(INCMUMPS)
