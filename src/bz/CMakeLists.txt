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
    HEigensystem add_missing_cells blas bz calculate
    cntour common_common digit draw euler graphics
    kpoint lower make_line minv orient point reada rot rotcel rotsec rsp
    setup solir sort surfac symtrz upcase vastkind
   )
#-----------------------------------------------
# Add a list of source files to the target
if(BUILD_BZ)
  foreach(idx IN LISTS src_list)
          target_sources(mopac-bz PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/${idx}.F90)
  endforeach()
endif()
#===============================================
