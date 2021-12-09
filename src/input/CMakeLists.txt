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
# source
        ${CMAKE_CURRENT_SOURCE_DIR}/datin.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/digit.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/geo_ref.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/getdat.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/getgeg.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/getgeo.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/getpdb.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/getsym.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/gettxt.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/getval.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/init_filenames.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/myword.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/nuchar.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/readmo.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/wrt_diffs.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/wrtkey.F90
)
