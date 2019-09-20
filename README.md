# simplex

Implementation of LP solver based on simplex algorithm. It supports parsing .lp and .MPS input file formats.

## Compilation

Run `make` command (macOS, Linux). 

## External dependencies

* BLAS/LAPACK
* Boost Spirit X3

### BLAS/LAPACK on macOS

It is enough to link to Accelerate Framework.

### BLAS/LAPACK on Linux

On Linux it is recommended to build OpenBLAS from sources.

`git clone https://github.com/xianyi/OpenBLAS.git`

Change into cloned `OpenBLAS`, and run

`make FC=gfortran NO_LAPACK=0`

`make PREFIX=path_to_simplex/extern install`

`export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:path_to_libopenblas`

### Boost Spirit X3

Unzip and put `boost_1_71_0` (or later) into `extern/` directory (no build/linking is required for Spirit).

## Use

`simplex [LP/MPS filename]`
