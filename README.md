# *k*2*l*

Fortran codes for computing the *k*-th to *l*-th eigenpair(s) of [symmetric definite generalized eigenproblems](http://www.netlib.org/lapack/lug/node54.html) of matrices *A* and *B*.

*k*2*l* is based on [k-ep](https://github.com/lee-djl/k-ep), a collection of fortran codes for computing the *k*-th eigenpair of the eigenproblem.

See [References](#references) for details about algorithms of k-ep.

## Getting started

### Prerequisites

* [LAPACK](http://www.netlib.org/lapack/) ver. 3.0 or later  
* [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) ver. 5.1.0 or later  
* [MUMPS](http://mumps.enseeiht.fr/) ver. 5.0.0 or later
    * Double-precision MUMPS (dMUMPS) is required.
    * Both sequential and parallel versions of MUMPS can be used. 
    * [ScaLAPACK](http://www.netlib.org/scalapack/), [BLACS](https://www.netlib.org/blacs/), and MPI are required to use parallel version of MUMPS.

### Installation

1. Download and uncompress the package.  
    ```
    wget --no-check-certificate https://codeload.github.com/lee-djl/k-ep/tar.gz/k2l -O k-ep-k2l.tar.gz
    tar zxf k-ep-*
    cd k-ep-*
    ```
2. Edit `Makefile.inc` to set compiler options and link the prerequisite libraries.
    * If sequential MUMPS is used, refer to [Makefile_seq.inc](/Makefile_seq.inc). 
    * If parallel MUMPS is used, refer to [Makefile_par.inc](/Makefile_par.inc). 

3. `make`.

### Run the example program

1. Example program uses matrix data in [ELSES Matrix Library](http://www.elses.jp/matrix/).

   Download matrix data [APF4686](http://www.damp.tottori-u.ac.jp/~hoshi/elses_matrix/ELSES_MATRIX_APF4686_20170505.tgz) to `example` directory and uncompress it.
    ```
    cd example
    wget http://www.damp.tottori-u.ac.jp/~hoshi/elses_matrix/ELSES_MATRIX_APF4686_20170505.tgz
    tar zxf ELSES_MATRIX_APF4686_*
    ```

2. Execute the example program `example.out` in `example` directory.  
    ```
    ./example.out \
    ./ELSES_MATRIX_APF4686_20170505/ELSES_MATRIX_APF4686_A.mtx \
    ./ELSES_MATRIX_APF4686_20170505/ELSES_MATRIX_APF4686_B.mtx \
    2323 \
    2343 \
    > log.txt
    ```

3. Check results `eval00002323to00002343.txt` and `ipr00002323to00002343.txt` in `example` directory.
   
    * `eval00002323to00002343.txt`: Contains the 2323-th to 2343-th eigenvalues.
    * `ipr00002323to00002343.txt`: Contains the 4-th power of the 4-norm of the 2323-th to 2343-th eigenvectors. (In quantum physics, these scalar values are referred to as the [inverse participation ratio](https://en.wikipedia.org/wiki/Purity_(quantum_mechanics)#Inverse_Participation_Ratio_(IPR)).)
   
    Compare the results with samples [s_eval00002323to00002343.txt](/example/s_eval00002323to00002343.txt) and [s_ipr00002323to00002343.txt](/example/s_ipr00002323to00002343.txt) in the same directory.  
    Note that the contents of the results may differ from those of the samples depending on computational environment.

### Deletion
1. Move to the top directory of the package.

2. `make clean`.

## How to use *k*2*l*

See [README.md](/example/README.md) in `example` directory.

## Acknowledgement

Matrix data were achieved from [ELSES Matrix Library](http://www.elses.jp/matrix/).

## References

D. Lee, T. Hoshi, T. Sogabe, Y. Miyatake, S.-L. Zhang, J. Comput. Phys., 371 (2018), 618-632. ([Journal](https://doi.org/10.1016/j.jcp.2018.06.002), [Preprint](https://arxiv.org/abs/1710.05134))
