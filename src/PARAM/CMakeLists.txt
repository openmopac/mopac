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
        ${CMAKE_CURRENT_SOURCE_DIR}/Param_global_C.F90
# source
        ${CMAKE_CURRENT_SOURCE_DIR}/chkpar.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/datinp.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/depfn.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/direct.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/empire.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/extract_parameter.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/filusp.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/finish.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/getmol.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/getpar.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/optgeo.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/pardip.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/parfg.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/pargeo.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/parips.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/parkey.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/partab.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/pparsav.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/psort.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/rapid0.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/rapid1.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/rapid2.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/rapid3.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/resetp.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/savgeo.F90
)
