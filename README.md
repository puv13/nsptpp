# A C++ Framework for NSPT Simulations

NSPTPP is a generic framework for Numerical Stochastic Perturbation Theory (NSPT)
simulations. The code is written in C++ and heavily relies on templates to make an
adaption for different physical models as simple as possible. 

## Dependencies 

###  Googletest 

The [googletest](https://github.com/google/googletest) library is used for unit
testing. ([BSD-3-Clause
License](https://github.com/google/googletest/blob/master/LICENSE)).


### HighFive

HDF5 file support is provided through the [HighFive
library](https://github.com/BlueBrain/HighFive). ([Boost
Software License 1.0](https://github.com/BlueBrain/HighFive/blob/master/LICENSE)).

### CMake 

CMake (in a version > 3.1) is needed to build the NSPTPP code. CMake is 
available from [cmake.org](https://cmake.org/) ([BSD 3-clause License](https://cmake.org/licensing/)).

### Doxygen 

The documentation for NSPTPP is generated with
[doxygen](http://www.doxygen.nl/index.html). ([GPL 2.0](http://www.gnu.org/licenses/old-licenses/gpl-2.0.html)). 

### Compiler 

The C++ compiler used to build NSPTPP has to support the C++11 standard. 


## Getting started 

### Obtain the NSPTPP code and the HighFive library  

```shell
# Get NSPTPP
git clone git@github.com:puv13/nsptpp.git

# Get HighFive
cd nsptpp  
git submodule update --init --recursive
```

### Build the code and run the tests 

```shell

mkdir build 
cd build 

# Run CMake 
# Replace <gcc> and <g++> with your C and C++ compiler, respectively.
CC=<gcc> CXX=<gpp> cmake ../ 

# Make tests 
make alltests 

# Run tests 
./tests/alltests 

```

### Build the documentation 

```shell

# Run this inside the nsptpp folder 
doxygen Doxyfile

```


## Contributors
Unless otherwise stated the code was written by Jakob Simeth and Matthias Puhr. The code
originated from research projects in the area of theoretical high energy physics. In the
[contributors list](Contributors.md) we therefore also include people who contributed
directly or indirectly to the projects and not just the authors of the source code.

## [Literature and References](Literature.md)
















