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
# source
        ${CMAKE_CURRENT_SOURCE_DIR}/H_bond_correction_EC_plus_ER.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/H_bond_correction_EH_plus.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/H_bond_correction_PM6_DH_Dispersion.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/H_bond_correction_PM6_DH_type.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/H_bond_correction_bits.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/H_bonds4.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/Hydrogen_bond_corrections.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/copyc6.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/dftd3.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/dftd3_bits.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/disp_DnX.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/gdisp.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/post_scf_corrections.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/set_up_dentate.F90
)
