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
        ${CMAKE_CURRENT_SOURCE_DIR}/ddiag.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/densit.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/density_for_GPU.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/diag.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/diag_for_GPU.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/dot.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/eigenvectors_LAPACK.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/esp_utilities.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/interp.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/linpack.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/mamult.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/mat33.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/minv.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/molval.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/mtxm.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/mtxmc.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/mult.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/mult_symm_AB.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/mxm.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/mxmt.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/mxv.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/osinv.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/rsp.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/schmib.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/schmit.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/supdot.F90
)
