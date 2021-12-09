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

target_sources(mopac_core PRIVATE
# modules
        ${CMAKE_CURRENT_SOURCE_DIR}/iter_C.F90
# source
        ${CMAKE_CURRENT_SOURCE_DIR}/capcor.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/cnvg.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/fock1.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/fock2.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/iter.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/local.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/pulay.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/resolv.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/swap.F90
)
