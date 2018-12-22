## How to use *k*2*l*

### Coding

See description below together with [example.f90](/example/example.f90).

1. Use the following modules.

    * `k2l_iotype`: Defines a derived type for using *k*2*l* solver.
    * `k2l_mod`: Includes *k*2*l* solver.

2. Declare a variable of the derived type. 
    ```
    TYPE(k2l_io_type) :: k2l_io
    ```

3. Provide information of input matrices *A* and *B*.

    * `k2l_io%n`: INTEGER. Size of *A* and *B*.
    * `k2l_io%nz_a`, `k2l_io%nz_b`: INTEGER. Number of nonzero elements of *A* and *B*.

4. Initialize *k*2*l* solver.
    ```
    k2l_io%job=0
    CALL k2l(k2l_io)
    ```

5. Assign values of input matrices *A* and *B* to `k2l_io`.

    Input matrices in [coordinate format](https://math.nist.gov/MatrixMarket/formats.html) (also referred to as [triplet format](https://people.sc.fsu.edu/~jburkardt/data/st/st.html) or assembled format) are supported.

    Assign values to the following variables.

    * `k2l_io%indx_a`, `k2l_io%indx_b`: INTEGER. One-dimensional array. Row indices of nonzero elements of *A* and *B*.
    * `k2l_io%jndx_a`, `k2l_io%jndx_b`: INTEGER. One-dimensional array. Column indices of nonzero elements of *A* and *B*.
    * `k2l_io%rval_a`, `k2l_io%rval_b`: DOUBLE PRECISION. One-dimensional array. Values of nonzero elements of *A* and *B*.

6. Set the target index range.

    * `k2l_io%k_lower`: INTEGER. Lower end of the index range.
    * `k2l_io%k_upper`: INTEGER. Upper end of the index range. Needs to be larger than or equal to `k2l_io%k_lower`.

    It is recommended that `k2l_io%k_upper`-`k2l_io%k_lower` be small (say, twenty) and be not too large (say, less than one hundred).

7. Optional: Set parameters.

    * `k2l_io%cprm(1)`: CHARACTER(16). Default value is null. To print details of computation to terminal, set the value as follows.
        ```
        k2l_io%cprm(1)='print'
        ```
    * `k2l_io%cprm(2)`: CHARACTER(16). Default value is null. To provide an initial interval containing the `k2l_io%k_lower`-th to `k2l_io%k_upper`-th eigenvalue(s), set the value as follows.
        ```
        k2l_io%cprm(2)='user'
        ```
    * `k2l_io%s_lower`: DOUBLE PRECISION. Lower endpoint of the user-defined initial interval. 
    * `k2l_io%s_upper`: DOUBLE PRECISION. Upper endpoint of the user-defined initial interval. Needs to be larger than `k2l_io%s_lower`.
    * `k2l_io%cprm(3)`: CHARACTER(16). Default value is null. To compute the 4-th power of the 4-norm of computed eigenvectors ([inverse participation ratio](https://en.wikipedia.org/wiki/Purity_(quantum_mechanics)#Inverse_Participation_Ratio_(IPR))s), set the value as follows.
        ```
        k2l_io%cprm(3)='ipr'
        ```
    * `k2l_io%cprm(6)`: CHARACTER(16). Default value is null. To compute only an interval containing the `k2l_io%k_lower`-th to `k2l_io%k_upper`-th eigenvalue(s), set the value as follows. Note that eigenpairs are not computed.
        ```
        k2l_io%cprm(6)='second'
        ```
    * `k2l_io%iprm(10)` and `k2l_io%dprm(1)`: INTEGER(8) and DOUBLE PRECISION. Stopping criterion for the second stage (bisection). Default values are 20 and 2.0, respectively. The initial interval is narrowed down until the number of eigenvalues in the interval becomes smaller than or equal to MAX(`k2l_io%iprm(10)`,CEILING(`k2l_io%dprm(1)`*(`k2l_io%k_upper`-`k2l_io%k_lower`+1)). If users want for the initial interval to be narrowed down to contain only the `k2l_io%k_lower`-th to `k2l_io%k_upper`-th eigenvalue(s), set the values as follows. Note that it is recommended to change these parameters only if users know how the three-stage algorithm works.
        ```
        k2l_io%iprm(10)=1
        k2l_io%dprm(1)=1.0D0
        ```

8. Solve.
    ```
    CALL mpi_init(ierr)

    k2l_io%job=1
    CALL k2l(k2l_io)

    CALL mpi_finalize(ierr)
    ```

    It is necessary to call MPI subroutines `mpi_init` and `mpi_finalize` to use MUMPS.

    * If sequential MUMPS is used, it is not required to link MPI library during compilation because MUMPS has its own dummy MPI subroutines. 
    * If parallel MUMPS is used, it is required link MPI library during compilation.
    * See [Compiling](#compiling) section below for details.

9. Extract output results.

    * `k2l_io%kndx`: INTEGER. One-dimensional array. Contains the target index (indices) in increasing order.
    * `k2l_io%kval`: DOUBLE PRECISION. One-dimensional array. Contains the `k2l_io%k_lower`-th to `k2l_io%k_upper`-th eigenvalue(s) in increasing order.
    * `k2l_io%kvec`: DOUBLE PRECISION. Two-dimensional array. Columns of the array are the `k2l_io%k_lower`-th to `k2l_io%k_upper`-th eigenvectors, normalized with respect to 2-norm.
    * `k2l_io%kipr`: DOUBLE PRECISION. One-dimensional array. Contains the 4-th power of the 4-norm of the `k2l_io%k_lower`-th to `k2l_io%k_upper`-th eigenvectors, normalized with respect to *B*-norm.

    Following values are accessible if `k2l_io%cprm(6)`='second'. After the second stage (bisection), interval [`k2l_io%s_lower2`,`k2l_io%s_upper2`) containing the [`k2l_io%k_lower2`,`k2l_io%k_upper2`]-th eigenvalue(s) is obtained.

    * `k2l_io%s_lower2`: DOUBLE PRECISION. Lower endpoint of the interval after the second stage. 
    * `k2l_io%s_upper2`: DOUBLE PRECISION. Upper endpoint of the interval after the second stage.
    * `k2l_io%k_lower2`: INTEGER. Lower end of the index range after the second stage.
    * `k2l_io%k_upper2`: INTEGER. Upper end of the index range after the second stage.

10. Finalize.
    ```
    k2l_io%job=-1
    CALL k2l(k2l_io)
    ```

### Compiling

See description below. Refer to [Makefile.inc](/Makefile.inc) and [Makefile](/example/Makefile) for details.

* Prerequisites

    * Link the prerequisite libraries, including MUMPS, METIS, and LAPACK, by setting `LIBS` and related values in [Makefile.inc](/Makefile.inc).
    * Add search paths for MUMPS by setting `MUMPSINC` in [Makefile.inc](/Makefile.inc).

* MPI    

    * If sequential MUMPS is used, setting `MPILIB`, `MUMPSMPI`, and their related values as [Makefile_seq.inc](/Makefile_seq.inc).
    * If parallel MUMPS is used, setting `MPILIB` as [Makefile_par.inc](/Makefile_par.inc). Setting `MUMPSMPI` is unnecessary.

* *k*2*l*

    * Link *k*2*l* by setting `K2LLIB` and related values in [Makefile](/example/Makefile).
    * Add search paths for *k*2*l* by setting `K2LINC` and related values in [Makefile](/example/Makefile).

### Execution

Example job scripts can be found in [script](/example/script) directory.
