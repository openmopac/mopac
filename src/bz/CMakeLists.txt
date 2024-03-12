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
