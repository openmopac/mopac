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
    deriv deri1 freqcy jcarin dcart deri23 delsta
    dfield dfock2 xyzcry forsav deri22 anavib
    deri0 intfc dhcore dhc deri2 force dijkl1
    dijkl2 dernvo deri21 fmat frame deritr
   )
#-----------------------------------------------
# Add a list of source files to the target
foreach(idx IN LISTS src_list)
        target_sources(mopac-core PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/${idx}.F90)
endforeach()
#===============================================
