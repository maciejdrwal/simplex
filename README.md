# simplex

Implementation of Linear Programming solver based on Simplex algorithm. It supports parsing .lp file format.

## Building

Make sure [CMake](http://cmake.org) is installed and available in PATH. To build please use:

`mkdir build & cd build & cmake .. & cmake --build .`

Alternatively, you can also use the Makefile provided and `make` command (macOS, Linux).

CMake will also generate Visual Studio solution for Windows.

## External dependencies

- [Eigen](http://eigen.tuxfamily.org/) - C++ template library for linear algebra

## Use

`simplex [LP filename]`