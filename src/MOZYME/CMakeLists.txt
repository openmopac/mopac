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
    mbonds fillij check hbonds buildf add_more_interactions
    isitsc setupk chrge_for_MOZYME reorth set_up_MOZYME_arrays
    eimp vecprt_for_MOZYME outer2 compct cnvgz diagg
    kab_for_MOZYME bonds_for_MOZYME prtlmo density_for_MOZYME
    addhb dipole_for_MOZYME local_for_MOZYME fock2z local2
    helecz scfcri MOZYME_C epseta iter_for_MOZYME pinout
    denrot_for_MOZYME fock1_for_MOZYME jab_for_MOZYME hybrid
    hcore_for_MOZYME tidy makvec selmos diagg2 mlmo
    convert_storage outer1 diagg1 adjvec MOZYME_eigs
   )
#-----------------------------------------------
# Add a list of source files to the target
foreach(idx IN LISTS src_list)
        target_sources(mopac-core PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/${idx}.F90)
endforeach()
#===============================================
