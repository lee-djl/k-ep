#!/bin/sh
#QSUB -queue set_the_queue
#QSUB -node 1
#QSUB -mpi 1
#QSUB -omp 1
#QSUB -place pack
#QSUB -over false
#PBS -l walltime=00:30:00
#PBS -N k2l
#PBS -j oe
#PBS -m abe

EXECPATH=/set_the_path_to/k-ep-k2l/example
EXECFILE=example.out

DATAPATH=/set_the_path_to/ELSES_MATRIX_APF4686_20170505
DATAFILEA=ELSES_MATRIX_APF4686_A.mtx
DATAFILEB=ELSES_MATRIX_APF4686_B.mtx

KNDXLOW=2323
KNDXUPP=2343

INTVLOW=-5.2D-1
INTVUPP=1.8D-1

LOGFILE=log.txt

cd ${PBS_O_WORKDIR}

${EXECPATH}/${EXECFILE} \
${DATAPATH}/${DATAFILEA} ${DATAPATH}/${DATAFILEB} \
${KNDXLOW} ${KNDXUPP} \
${INTVLOW} ${INTVUPP} \
&> ${LOGFILE}
