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

macro(add_mopac_setup _name _files)
  # _files is not a list, it's a string... Transform it into a list
  set(files)
  string(REPLACE ";" " " _files "${_files}")
  foreach(_label "${_files}")
    list(APPEND files ${_label})
  endforeach()
  unset(_files)

  add_test(NAME "${_name}" COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX} ${_files} SETUP)
endmacro()

macro(add_mopac_test _name _files)
  # _files is not a list, it's a string... Transform it into a list
  set(files)
  string(REPLACE ";" " " _files "${_files}")
  foreach(_label "${_files}")
    list(APPEND files ${_label})
  endforeach()
  unset(_files)

  add_test(NAME "${_name}" COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX} ${_files})
endmacro()

add_subdirectory(keywords)

add_mopac_setup(port "port.mop;port.parameters;try.txt;aa.txt;mol.in")

add_mopac_test(crambin "Crambin_1SCF.mop;Crambin_1SCF.den")
add_mopac_test(mozyme "test_Lewis_for_Proteins.mop")
add_mopac_test(xeno "Test_XENO.mop")
add_mopac_test(hetero_group "Test_Hetero_group.mop")
add_mopac_test(convert_edited_pdb "convert_an_edited_PDB_file.mop")
add_mopac_test(convert_unedited_pdb "convert_an_un-edited_PDB_file.mop")
add_mopac_test(ters_preserved "Check_that_TERs_are_preserved.mop;Check_TERs.pdb")
add_mopac_test(compare_proteins "Compare_two_protein_structures.mop;Structure-1.txt;Structure-2.txt")
add_mopac_test(pdb1cbn "PDB1CBN.ent")
add_mopac_test(pdb1ejg "PDB1EJG.ent")
add_mopac_test(pdb_compare "Compare_PDB_1EJG_and_1CBN.mop;PDB1CBN.ent;PDB1EJG.ent")
add_mopac_test(pdb_compare2 "Compare_PDB_Bacteriorhodopsin_1AP9_and_2AT9.mop;Bacteriorhodopsin_1AP9.ent;Bacteriorhodopsin_2AT9.ent")
add_mopac_test(pdb_compare3 "Compare_PDB_Bacteriorhodopsin_1AP9_and_1C3W.mop;Bacteriorhodopsin_1AP9.ent;Bacteriorhodopsin_1C3W.ent")
add_mopac_test(pdb_compare4 "Compare_PDB_Bacteriorhodopsin_1AP9_and_1BRX.mop;Bacteriorhodopsin_1AP9.ent;Bacteriorhodopsin_1BRX.ent")
add_mopac_test(pdb_compare5 "Compare_PDB_Bacteriorhodopsin_1AP9_and_1BRR.mop;Bacteriorhodopsin_1AP9.ent;Bacteriorhodopsin_1BRR.ent")
add_mopac_test(pdb_compare6 "Compare_PDB_Bacteriorhodopsin_1AP9_and_1AT9.mop;ADD-H_Bacteriorhodopsin_1AP9.arc;ADD-H_Bacteriorhodopsin_1AT9.arc")
add_mopac_test(add-h "ADD-H_Bacteriorhodopsin_1AP9.mop;Bacteriorhodopsin_1AP9.ent")
add_mopac_test(add-h2 "ADD-H_Bacteriorhodopsin_1AT9.mop;Bacteriorhodopsin_1AT9.ent")
add_mopac_test(add-h3 "Au_Add-H.mop")
add_mopac_test(compare_add-h "Compare_ADD-H_1EJG_and_1CBN.mop;PDB1EJG.arc;PDB1CBN.arc")
add_mopac_test(locate_ts "Locate_TS.mop;Cyclobutadiene_1.arc;Cyclobutadiene_2.arc")
add_mopac_test(ionize "IONIZE.mop;Neutral.mop")
add_mopac_test(salt "SALT.mop;Neutral.mop")
add_mopac_test(salt2 "SALT_2.mop")
add_mopac_test(residues "SITE_Individual_residues.mop;IONIZE.arc")
add_mopac_test(ag20 "ag20.mop")
add_mopac_test(ag20_cis_active "ag20_cis_active.mop")
add_mopac_test(benzene "benzene.mop")
add_mopac_test(co_cas "co_cas.mop")
add_mopac_test(co_cis "co_cis.mop")
add_mopac_test(co_cis_2 "co_cis_2.mop")
add_mopac_test(co_cis_PM7 "co_cis_PM7.mop")
add_mopac_test(co_cis_active "co_cis_active.mop")
add_mopac_test(co_cis_solv "co_cis_solv.mop")
add_mopac_test(co_cis_tdip "co_cis_tdip.mop")
add_mopac_test(co_cis_trip "co_cis_trip.mop")
add_mopac_test(co_extern "co_extern.mop;param.txt")
add_mopac_test(co_fail "co_fail.mop")
add_mopac_test(co_field "co_field.mop")
add_mopac_test(co_mrci "co_mrci.mop")
add_mopac_test(co_mrci_maxci "co_mrci_maxci.mop")
add_mopac_test(co_mrci_solv "co_mrci_solv.mop")
add_mopac_test(co_pm7 "co_pm7.mop")
add_mopac_test(co_pm7_solv "co_pm7_solv.mop")
add_mopac_test(co_print "co_print.mop")
add_mopac_test(no_cis_quar "no_cis_quar.mop")
add_mopac_test(no_mrci "no_mrci.mop")
add_mopac_test(noparam "noparam.mop")
add_mopac_test(poly_c8 "poly_c8.mop")
add_mopac_test(poly_c8_mr "poly_c8_mr.mop")
