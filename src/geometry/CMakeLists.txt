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

target_sources(mopac-all PRIVATE
# modules
# source
        ${CMAKE_CURRENT_SOURCE_DIR}/axis.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/bangle.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/dang.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/dihed.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/gmetry.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/symtry.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/volume.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/xyzint.F90
)
