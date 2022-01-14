# Molecular Orbital PACkage (MOPAC)
# Copyright (C) 2021, Virginia Polytechnic Institute and State University
#
# MOPAC is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MOPAC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# This CMake file contains the list of source files, mostly for the library

# source files for the executable stubs
target_sources(mopac PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/mopac.F90)
target_sources(param PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/param.F90)

# source subdirectories
add_subdirectory(INDO)
add_subdirectory(MOZYME)
add_subdirectory(PARAM)
add_subdirectory(SCF)
add_subdirectory(chemistry)
add_subdirectory(corrections)
add_subdirectory(forces)
add_subdirectory(geometry)
add_subdirectory(input)
add_subdirectory(integrals)
add_subdirectory(matrix)
add_subdirectory(meci)
add_subdirectory(models)
add_subdirectory(optimization)
add_subdirectory(output)
add_subdirectory(potentials)
add_subdirectory(properties)
add_subdirectory(reactions)
add_subdirectory(solvation)
add_subdirectory(symmetry)
add_subdirectory(utilities)

# source files for the library
target_sources(mopac_core PRIVATE
# modules
        ${CMAKE_CURRENT_SOURCE_DIR}/Common_arrays_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/chanel_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/conref_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/elemts_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/funcon_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/molkst_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/molmec_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/overlaps_C.F90
# source
        ${CMAKE_CURRENT_SOURCE_DIR}/GetDateStamp.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/compfg.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/moldat.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/mopend.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/run_mopac.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/run_param.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/setup_mopac_arrays.F90
)
