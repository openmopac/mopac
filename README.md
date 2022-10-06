# Molecular Orbital PACkage (MOPAC)

[![License: LGPL v3](https://img.shields.io/badge/License-LGPL_v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)
[![DOI](https://zenodo.org/badge/177640376.svg)](https://zenodo.org/badge/latestdoi/177640376)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/mopac/badges/version.svg)](https://anaconda.org/conda-forge/mopac)
![build](https://github.com/openmopac/mopac/actions/workflows/CI.yaml/badge.svg)
[![codecov](https://codecov.io/gh/openmopac/mopac/branch/main/graph/badge.svg?token=qM2KeRvw06)](https://codecov.io/gh/openmopac/mopac)

This is the official repository of the modern open-source version of MOPAC, which is now released under the LGPL license.
This is a direct continuation of the commercial development and distribution of MOPAC, ending at MOPAC 2016.
Commercial versions of MOPAC are no longer supported, and all MOPAC users are encouraged to switch to the most recent open-source version.

## Installation

Self-contained installers for Linux, Mac, and Windows are available on GitHub for each release,
which are constructed using the [Qt Installer Framework](https://doc.qt.io/qtinstallerframework/).

While the installers are meant to be run from a desktop environment by default, they can also be run from a command line without user input.
On Linux, the basic command-line installation syntax is:

`./mopac-x.y.z-linux.run install --accept-licenses --confirm-command --root type_installation_directory_here`

For more information on command-line installation, see the [Qt Installer Framework Documentation](https://doc.qt.io/qtinstallerframework/ifw-cli.html).

### Package managers

The officially supported package manager for MOPAC is the [conda-forge channel of Conda](https://anaconda.org/conda-forge/mopac).
MOPAC is also available on [Fedora](https://packages.fedoraproject.org/pkgs/mopac/mopac/)
and the [Google Play store](https://play.google.com/store/apps/details?id=cz.m).

We welcome the addition of MOPAC to other package managers and will add them to this list as they are made available and known.

### CMake

MOPAC is now built using a CMake 3.x build system with tests orchestrated using CTest.
The minimum required CMake version is presently 3.14.

CMake performs out-of-source builds, with the canonical sequence of commands:

```
mkdir build
cd build
cmake ..
make
```

starting from the root directory of the MOPAC repository. MOPAC should build without any additional options
if CMake successfully detects a Fortran compiler and BLAS/LAPACK libraries. Otherwise, the `cmake ..` command
may require additional command-line options to specify a Fortran compiler (`-DCMAKE_Fortran_COMPILER=...`)
or the path (`-DMOPAC_LINK_PATH=...`) and linker options (`-DMOPAC_LINK=...`) to link BLAS and LAPACK libraries to the MOPAC executable.

The CTest-based testing requires Python 3.x and Numpy.

## Documentation

The main source for MOPAC documentation is presently its old [online user manual](http://openmopac.net/manual/index.html).

There is a [new documentation website](https://openmopac.github.io) under development, but it is not yet ready for general use.
