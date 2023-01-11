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

#===============================================
# Create a list of source files (src_list) with the .F90 extension
set(src_list 
    parameters_for_RM1_C parameters_for_PM3_Sparkles_C
    parameters_for_PM7_C parameters_for_PM8_C
    parameters_for_PM6_ORG_C
    alpb_and_xfac_mndod refer alpb_and_xfac_mndo
    calpar journal_references_C alpb_and_xfac_pm3
    parameters_for_AM1_C parameters_for_PM7_Sparkles_C
    parameters_for_AM1_Sparkles_C parameters_C
    alpb_and_xfac_am1 parameters_for_mndo_C
    parameters_for_PM6_Sparkles_C update
    parameters_for_RM1_Sparkles_C
    parameters_for_PM6_C parameters_for_INDO_C
    parameters_for_mndod_C parameters_for_PM7_TS_C
    parameters_for_PM3_C switch
   )
#-----------------------------------------------
# Add a list of source files to the target
foreach(idx IN LISTS src_list)
        target_sources(mopac-core PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/${idx}.F90)
endforeach()
#===============================================
