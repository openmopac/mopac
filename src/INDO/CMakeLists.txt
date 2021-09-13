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
        ${CMAKE_CURRENT_SOURCE_DIR}/reimers_C.F90
# interfaces
# source
        ${CMAKE_CURRENT_SOURCE_DIR}/ci.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/cosci.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/excite.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/input.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/integral.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/ovlp.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/rci.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/scf.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/subs.F90
)
