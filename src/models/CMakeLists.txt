# Molecular Orbital PACkage (MOPAC)
# Copyright 2021 Virginia Polytechnic Institute and State University
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

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
