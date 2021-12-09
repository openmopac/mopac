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
        ${CMAKE_CURRENT_SOURCE_DIR}/journal_references_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/parameters_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/parameters_for_AM1_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/parameters_for_AM1_Sparkles_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/parameters_for_INDO_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/parameters_for_PM3_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/parameters_for_PM3_Sparkles_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/parameters_for_PM6_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/parameters_for_PM6_Sparkles_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/parameters_for_PM7_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/parameters_for_PM7_Sparkles_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/parameters_for_PM7_TS_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/parameters_for_PM8_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/parameters_for_RM1_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/parameters_for_RM1_Sparkles_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/parameters_for_mndo_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/parameters_for_mndod_C.F90
# source
        ${CMAKE_CURRENT_SOURCE_DIR}/alpb_and_xfac_am1.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/alpb_and_xfac_mndo.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/alpb_and_xfac_mndod.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/alpb_and_xfac_pm3.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/calpar.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/refer.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/switch.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/update.F90
)
