# k-ep

Fortran codes for computing the k-th eigenpair of [generalized symmetric definite eigenproblems](http://www.netlib.org/lapack/lug/node54.html).

Input matrices in [Matrix Market exchange (MTX) format](http://math.nist.gov/MatrixMarket/formats.html) are currently supported.

See [References](#references) for details about the k-th eigenvalue problem and algorithms for it.

## Getting started

### Prerequisites

* [LAPACK](http://www.netlib.org/lapack/) ver. 3.0 or later  
* [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) ver. 5.1.0 or later  
* [MUMPS](http://mumps.enseeiht.fr/) ver. 5.0.0 or later

### Installing

1. Download and uncompress the package.  
```
wget --no-check-certificate https://codeload.github.com/lee-djl/k-ep/tar.gz/master -O k-ep-master.tar.gz
tar zxf k-ep-*
cd k-ep-*
```
2. Edit `Makefile.inc` to set compiler options and link the prerequisite libraries.

3. `make`.

### Running the example

1. Example program uses matrix data in [ELSES Matrix Library](http://www.elses.jp/matrix/).

   Download matrix data [APF4686](http://www.damp.tottori-u.ac.jp/~hoshi/elses_matrix/ELSES_MATRIX_APF4686_20170505.tgz) to `example` directory and uncompress it.
```
cd example
wget http://www.damp.tottori-u.ac.jp/~hoshi/elses_matrix/ELSES_MATRIX_APF4686_20170505.tgz
tar zxf ELSES_MATRIX_APF4686_*
```

2. Execute the example program `example.out` in `example` directory.  
```
./example.out ./ELSES_MATRIX_APF4686_20170505/ELSES_MATRIX_APF4686_A.mtx \
./ELSES_MATRIX_APF4686_20170505/ELSES_MATRIX_APF4686_B.mtx > output.txt
```

3. Check results `output.txt`, `eval00002343.txt`, and `evec00002343.txt` in `example` directory.

   * `output.txt` shows some implementation details of the program.
   * `eval00002343.txt` contains the index and value of the k-th eigenvalue (k = 2343).
   * `evec00002343.txt` contains elements of the k-th eigenvector normalized with respect to 2-norm.
   
   Compare the results with samples [s_output.txt](/example/s_output.txt), [s_eval00002343.txt](/example/s_eval00002343.txt), and [s_evec00002343.txt](/example/s_evec00002343.txt) in the same directory.  
   Note that the contents of the results may differ from those of the samples depending on computational environment.

## Acknowledgement

Matrix data used in the example program are achieved from [ELSES Matrix Library](http://www.elses.jp/matrix/).

## References

D. Lee, T. Hoshi, T. Sogabe, Y. Miyatake, S.-L. Zhang, [Solution of the k-th eigenvalue problem in large-scale electronic structure calculations](http://arxiv.org/abs/1710.05134), submitted.
