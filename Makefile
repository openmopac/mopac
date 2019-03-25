##############################################################################
#
#   Compile and build MOPAC2016 for Mac OS-X
#
##############################################################################
#
#  The following set of definitions apply to this OS.
#
##############################################################################
MOPAC_SRC = ./src
OBJ       = ./obj
MKL_ROOT   = /opt/intel/compilers_and_libraries_2018.3.185/mac/mkl
COMPILER_ROOT = /opt/intel/compilers_and_libraries_2018.3.185/mac/compiler
FORTRAN_COMPILER = ifort -c -O3 -I$(OBJ) -module $(OBJ) -o $@
FORTRAN_LINKER   = ifort -lpthread -lstdc++ -O3
O                = o 
MATH_LINK        =   $(MKL_ROOT)/lib/libmkl_intel.a \
                     $(MKL_ROOT)/lib/libmkl_intel_thread.a \
                     $(MKL_ROOT)/lib/libmkl_core.a  \
                     $(COMPILER_ROOT)/lib/libiomp5.a \
                     -lpthread \
                     -lm
EXE              =   ./bin/mopac2016.exe
include MOPAC_Makefile_files.txt
