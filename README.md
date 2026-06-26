apoCHARMM: High-performance molecular dynamics simulations on GPUs for advanced
simulation methods
================================================================================

## About

apoCHARMM is a GPU-only molecular dynamics package. Features that are used in
the working unit tests (c.f. run_tests.sh) can be used for performing molecular
dynamics simulations. The Python API is currently being overhauled and is
deprecated.

Any work that uses apoCHARMM should cite the following publication:

[S. Prasad, F. Aviat, J. E. Gonzales II, and B. R. Brooks, "apoCHARMM:
High-performance molecular dynamics simulations on GPUs for advanced simulation
methods," J. Chem. Phys. 162, 182501
(2025).](https://pubs.aip.org/aip/jcp/article/162/18/182501/3346618)

## License

apoCHARMM is distributed under the
[BSD-3-clause](https://opensource.org/licenses/BSD-3-Clause) open source
license, as described in the LICENSE file in the top level of the repository.
The repository also ships a single-header copy of
[Catch2](https://github.com/catchorg/Catch2) for the C++ unit tests. Catch2 is
distributed under the
[Boost Software License 1.0](https://opensource.org/licenses/BSL-1.0).

## Dependencies

apoCHARMM does not require any git submodules.

The C++ unit tests use the vendored single-header Catch2 file at
`test/unittests/catch.hpp`. This file is shipped with the test sources and is
not fetched as an external library during configuration or build.

The source code was developed using the tool and compiler versions listed in the
installation section below. Other versions may work.

## Authors

Samarjeet Prasad (Nvidia)

James E. Gonzales II (NIH)

Bernard R. Brooks (NIH)

## Installation

The source code was developed using the following tool and compiler versions.
Other versions may work.

* GCC [12.2.0]
* CUDA [12.2.140]
* CMake [3.25.1]

### 0. Clone this repository

```
git clone git@github.com:jeg7/apocharmm
cd apocharmm/
```

### 1. Compile the source code

By default, the Python API and unit tests are built.

For a standard release build:

```
cmake -S . -B build
cmake --build build --parallel
```

For a debug build:

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel
```

To build without the Python C API:

```
cmake -S . -B build -DAPOCHARMM_ENABLE_PYTHON_API=OFF
cmake --build build --parallel
```

To build without unit tests:

```
cmake -S . -B build -DBUILD_TESTING=OFF
cmake --build build --parallel
```

To restrict CUDA architectures, set `APOCHARMM_CUDA_ARCHITECTURES`. For example,
to build only for NVIDIA Ampere and Hopper GPUs:

```
cmake -S . -B build -DAPOCHARMM_CUDA_ARCHITECTURES="80;90"
cmake --build build --parallel
```

### 2. Enable use of the Python API

In your shells's configuration script (typically `~/.bashrc`) add the following
lines to the end.

**Note that you should replace the `/absolute/path/to/` with the absolute path
that leads to the `apocharmm/` directory (i.e. Where you cloned it).**

```
# For apoCHARMM
export APOCHARMM_LIBRARY_PATH=/absolute/path/to/apocharmm/lib/libapocharmm_c.so
export PYTHONPATH=/absolute/path/to/apocharmm/python:$PYTHONPATH
```

## Testing and Examples

apoCHARMM uses CTest to run both the C++ unit tests and the Python API tests.
CTest is enabled when `BUILD_TESTING=ON`, which is the default CMake behavior.

After configuring and building apoCHARMM run all registered tests with:

```
ctest --test-dir build --output-on-failure
```

To run only the C++ API tests:

```
ctest --test-dir build -L cpp_api --output-on-failure
```

To run only the Python API tests:

```
ctest --test-dir build -L python_api --output-on-failure
```

To run a single test, use `-R` with the test name or a regular expression.

For example:

```
ctest --test-dir build -R python_api_charmm_psf --output-on-failure
```

or:

```
ctest --test-dir build -R unittest-cudaNoseHooverIntegrator --output-on-failure
```

The Python API tests are only registered when the Python C API is built. By
default, this is enabled. To disable the Python C API and Python API tests,
configure with:

```
cmake -S . -B build -DAPOCHARMM_ENABLE_PYTHON_API=OFF
```

To configure apoCHARMM without any tests:

```
cmake -S . -B build -DBUILD_TESTING=OFF
```

The convenience script [run_tests.sh](run_tests.sh) can also be used to run the
test suite:

```
./run_tests.sh
```

Example scripts of how to perform simple tasks in apoCHARMM can be found in the
[example directory](example).

## Getting Help, Feature Requests, and Contributing

Please contact <james.gonzales@nih.gov> for any of the following:

* Trouble installing
* Trouble using
* Bugs
* Feature requests
* Code/feature contributions
