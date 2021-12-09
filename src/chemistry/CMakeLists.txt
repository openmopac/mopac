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

target_sources(mopac-common PRIVATE
# modules
        ${CMAKE_CURRENT_SOURCE_DIR}/atomradii_C.F90
# source
        ${CMAKE_CURRENT_SOURCE_DIR}/add_hydrogen_atoms.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/atomrs.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/chkion.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/chklew.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/findn1.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/geochk.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/geochk_bits_1.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/geochk_bits_2.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/greek.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/lewis.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/ligand.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/lyse.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/names.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/nxtmer.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/orient_water.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/reseq.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/txtype.F90
)
