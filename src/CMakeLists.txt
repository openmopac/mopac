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

file(GLOB main_src_list ${CMAKE_CURRENT_SOURCE_DIR} *.F90)

list(REMOVE_ITEM
     main_src_list
     ${CMAKE_CURRENT_SOURCE_DIR}/mopac.F90
     ${CMAKE_CURRENT_SOURCE_DIR}/param.F90
    )

multi_src_target(mopac_core PRIVATE main_src_list)

set(dir_list
    INDO MOZYME PARAM SCF chemistry corrections forces geometry
    input integrals matrix meci models optimization output
    potentials properties reactions solvation symmetry utilities
   )

foreach(idx IN LISTS dir_list)
	add_subdirectory(${idx})
endforeach()
