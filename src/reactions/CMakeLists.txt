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
        ${CMAKE_CURRENT_SOURCE_DIR}/ef_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/maps_C.F90
# source
        ${CMAKE_CURRENT_SOURCE_DIR}/big_swap.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/drc.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/drcout.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/ef.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/grid.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/pathk.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/paths.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/picopt.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/prtdrc.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/quadr.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/react1.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/set_up_RAPID.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/write_trajectory.F90
)
