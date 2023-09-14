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

#===============================================
# source files for the executable stubs
target_sources(mopac PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/mopac.F90)
target_sources(mopac-param PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/param.F90)
#===============================================
# Take a list of source file names
set(main_src_list 
    conref_C run_mopac funcon_C
    molkst_C setup_mopac_arrays overlaps_C 
    molmec_C Common_arrays_C elemts_C run_param
    mopend moldat chanel_C compfg
   )
#-----------------------------------------------
# Add a list of source files to the target
foreach(idx IN LISTS main_src_list)
	target_sources(mopac-core PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/${idx}.F90)
endforeach()
#===============================================
# Make a list of subdirectories to be added to the build structure  
set(dir_list
    INDO MOZYME PARAM SCF bz chemistry corrections forces geometry
    input integrals interface makpol matrix meci models optimization
    output potentials properties reactions solvation symmetry utilities
   )
#-----------------------------------------------
# Add the subdirectories to the build directory
foreach(idx IN LISTS dir_list)
	add_subdirectory(${idx})
endforeach()
#===============================================
