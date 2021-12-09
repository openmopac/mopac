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
        ${CMAKE_CURRENT_SOURCE_DIR}/meci_C.F90
# source
        ${CMAKE_CURRENT_SOURCE_DIR}/aababc.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/aabacd.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/aabbcd.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/babbbc.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/babbcd.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/ciosci.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/diagi.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/formxy.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/ijkl.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/meci.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/mecid.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/mecih.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/mecip.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/partxy.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/perm.F90
)
