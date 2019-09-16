# simplex

Implementation of basic LP solver using simplex algorithm. It supports parsing .lp and .MPS input file formats.

## Compilation

Run `make` command (macOS, Linux). 

**Note** that BLAS and LAPACK are required. 

### BLAS/LAPACK on macOS

It is enough to link to Accelerate Framework.

### BLAS/LAPACK on Linux

On Linux it is recommended to build OpenBLAS from sources.

`git clone https://github.com/xianyi/OpenBLAS.git`

Change into cloned `OpenBLAS`, and run

`make FC=gfortran NO_LAPACK=0`

`make PREFIX=path_to_simplex/extern install`

`export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:path_to_libopenblas`

## Use

`simplex [LP/MPS filename]`
