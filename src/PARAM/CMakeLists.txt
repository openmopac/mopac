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
# Create a list of source files (src_list) with the .F90 extension
set(src_list
    getpar partab extract_parameter pardip psort
    rapid0 rapid1 rapid3 resetp getmol datinp
    direct parkey rapid2 parfg Param_global_C
    optgeo pparsav finish depfn savgeo pargeo
    empire parips filusp chkpar
   )
#-----------------------------------------------
# Add a list of source files to the target
foreach(idx IN LISTS src_list)
        target_sources(mopac-core PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/${idx}.F90)
endforeach()
#===============================================

