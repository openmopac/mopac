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
# source
        ${CMAKE_CURRENT_SOURCE_DIR}/anavib.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/dcart.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/delsta.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/deri0.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/deri1.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/deri2.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/deri21.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/deri22.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/deri23.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/deritr.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/deriv.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/dernvo.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/dfield.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/dfock2.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/dhc.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/dhcore.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/dijkl1.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/dijkl2.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/fmat.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/force.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/forsav.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/frame.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/freqcy.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/intfc.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/jcarin.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/xyzcry.F90
)
