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
        ${CMAKE_CURRENT_SOURCE_DIR}/symmetry_C.F90
# source
        ${CMAKE_CURRENT_SOURCE_DIR}/bldsym.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/charmo.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/charst.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/charvi.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/dtran2.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/dtrans.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/maksym.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/molsym.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/mult33.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/rotmol.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/symh.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/symoir.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/symopr.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/sympop.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/symr.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/symt.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/symtrz.F90
)
