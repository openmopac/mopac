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

#Basic input/output comparison tests

set(HOF_ERROR 0.01) # Default numerical threshold to flag HOF errors

macro(add_mopac_test _name _files)
  # _files is not a list, it's a string... Transform it into a list
  set(files)
  string(REPLACE ";" " " _files "${_files}")
  foreach(_label "${_files}")
    list(APPEND files ${_label})
  endforeach()
  unset(_files)

  add_test(NAME "${_name}" COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/run_test.py ${CMAKE_CURRENT_SOURCE_DIR}
           ${CMAKE_BINARY_DIR}/mopac${CMAKE_EXECUTABLE_SUFFIX} ${HOF_ERROR} ${_files})
endmacro()

add_subdirectory(INDO-dev)
add_subdirectory(keywords)

add_mopac_test(port "port.mop;port.parameters;try.txt;aa.txt;mol.in;SETUP")
add_mopac_test(crambin "Crambin_1SCF.mop")
add_mopac_test(mozyme "test_Lewis_for_Proteins.mop")
add_mopac_test(xeno "Test_XENO.mop")
add_mopac_test(hetero_group "Test_Hetero_group.mop")
add_mopac_test(convert_edited_pdb "convert_an_edited_PDB_file.mop")
add_mopac_test(convert_unedited_pdb "convert_an_un-edited_PDB_file.mop")
add_mopac_test(ters_preserved "Check_that_TERs_are_preserved.mop;Check_TERs.pdb")
add_mopac_test(compare_proteins "Compare_two_protein_structures.mop;Structure-1.txt;Structure-2.txt")
add_mopac_test(pdb1cbn "PDB1CBN.ent")
add_mopac_test(pdb1ejg "PDB1EJG.ent")
add_mopac_test(add-h "ADD-H_Bacteriorhodopsin_1AP9.mop;Bacteriorhodopsin_1AP9.ent")
add_mopac_test(add-h2 "ADD-H_Bacteriorhodopsin_1AT9.mop;Bacteriorhodopsin_1AT9.ent")
add_mopac_test(add-h3 "Au_Add-H.mop")
add_mopac_test(compare_add-h "Compare_ADD-H_1EJG_and_1CBN.mop")
add_mopac_test(pdb_compare "Compare_PDB_1EJG_and_1CBN.mop")
set_tests_properties(compare_add-h pdb_compare PROPERTIES DEPENDS "pdb1cbn;pdb1ejg") # uses PDB1EJG.arc/ent & PDB1CBN.arc/ent
add_mopac_test(pdb_compare2 "Compare_PDB_Bacteriorhodopsin_1AP9_and_2AT9.mop;Bacteriorhodopsin_2AT9.ent")
add_mopac_test(pdb_compare3 "Compare_PDB_Bacteriorhodopsin_1AP9_and_1C3W.mop;Bacteriorhodopsin_1C3W.ent")
add_mopac_test(pdb_compare4 "Compare_PDB_Bacteriorhodopsin_1AP9_and_1BRX.mop;Bacteriorhodopsin_1BRX.ent")
add_mopac_test(pdb_compare5 "Compare_PDB_Bacteriorhodopsin_1AP9_and_1BRR.mop;Bacteriorhodopsin_1BRR.ent")
add_mopac_test(pdb_compare6 "Compare_PDB_Bacteriorhodopsin_1AP9_and_1AT9.mop")
set_tests_properties(pdb_compare2 pdb_compare3 pdb_compare4 pdb_compare5 pdb_compare6 PROPERTIES DEPENDS "add-h;add-h2") # uses ADD-H_Bacteriorhodopsin_1AP9.arc & ADD-H_Bacteriorhodopsin_1AT9.arc or Bacteriorhodopsin_1AP9.ent & Bacteriorhodopsin_1AT9.ent
add_mopac_test(locate_ts "Locate_TS.mop;Cyclobutadiene_1.arc;Cyclobutadiene_2.arc")
add_mopac_test(ionize "IONIZE.mop;Neutral.mop")
add_mopac_test(salt "SALT.mop")
set_tests_properties(salt PROPERTIES DEPENDS ionize) # uses Neutral.mop
add_mopac_test(salt2 "SALT_2.mop")
add_mopac_test(residues "SITE_Individual_residues.mop")
set_tests_properties(residues PROPERTIES DEPENDS ionize) # uses IONIZE.arc

set(HOF_ERROR 1.0) # Raise threshold, these tests are only looking for crashes
add_mopac_test(erratic-1 "erratic1.mop")
add_mopac_test(erratic-2 "erratic2.mop")
add_mopac_test(erratic-3 "erratic3.mop")
