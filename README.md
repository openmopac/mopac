# Molecular Orbital PACkage (MOPAC)

[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![DOI](https://zenodo.org/badge/177640376.svg)](https://zenodo.org/badge/latestdoi/177640376)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/mopac/badges/version.svg)](https://anaconda.org/conda-forge/mopac)
![build](https://github.com/openmopac/mopac/actions/workflows/CI.yaml/badge.svg)
[![codecov](https://codecov.io/gh/openmopac/mopac/branch/main/graph/badge.svg?token=qM2KeRvw06)](https://codecov.io/gh/openmopac/mopac)

This is the official repository of the modern open-source version of MOPAC, which is now released under an Apache license
(versions 22.0.0 through 23.0.3 are available under an LGPL license).
This is a direct continuation of the commercial development and distribution of MOPAC, which ended at MOPAC 2016.
Commercial versions of MOPAC are no longer supported, and all MOPAC users are encouraged to switch to the most recent open-source version.

[![mopac_at_molssi](.github/mopac_at_molssi.png)](https://molssi.org)

MOPAC is actively maintained and curated by the [Molecular Sciences Software Institute (MolSSI)](https://molssi.org).

For detailed information about MOPAC, see its [website](https://openmopac.github.io). A brief summary for developers is provided below.

## About

MOPAC is a popular semiempirical quantum chemistry package that was first released in 1983.
MOPAC can perform quantum mechanical calculations of molecules and materials, much like *ab initio* quantum chemistry packages such as [Gaussian](https://www.gaussian.com).
Semiempirical calculations are around 1000x faster than *ab initio* calculations, but they are usually less accurate because they rely on
simplified, minimal-basis model Hamiltonians that are fit to experimental data rather than more predictive theories.
The speed and ease-of-use of MOPAC makes it a good teaching tool, either to teach physical chemistry concepts to students or to prepare yourself
to use more expensive quantum chemistry software. It can also be useful for high-throughput calculations in screening and informatics applications.

## Building

Pre-built versions of MOPAC on Linux, Mac, and Windows are available with every [GitHub Release](https://github.com/openmopac/mopac/releases)
and also available on the Conda package manager using the command:
```
conda install -c conda-forge mopac
```

MOPAC can be built using its CMake build system. This repository has a [GitHub Actions workflow](https://github.com/openmopac/mopac/blob/main/.github/workflows/CI.yaml)
to build MOPAC, and you can build it locally using the standard CMake out-of-source invocation:
```
mkdir build
cd build
cmake ..
make
```
MOPAC will build without any other CMake options if a Fortran compiler, BLAS/LAPACK, Python 3, and NumPy are found in the software environment.
To build with optional [MolSSI Driver Interface)(MDI)](https://molssi-mdi.github.io/MDI_Library/) engine support, use the CMake command-line option `-DMDI=ON`.

## Usage

While MOPAC is often used through its integration with other software such as [WebMO](https://www.webmo.net), it is primarily a command-line program
that operates through input and output files. Recent versions of MOPAC also has an application programming interface (API) for a subset of functionality.

Examples of command-line and API usage are provided in the [`/examples` folder](https://github.com/openmopac/mopac/tree/main/examples) of this repository.

## Citation

The recommended citation of open-source MOPAC in scientific publications is its archival copy on Zenodo at
[DOI:10.5281/zenodo.6511958](https://doi.org/10.5281/zenodo.6511958).
