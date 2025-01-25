# Molecular Orbital PACkage (MOPAC)
# Copyright 2021 Virginia Polytechnic Institute and State University
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

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
if(MDI)
  add_subdirectory(mdi)
endif()
#===============================================
