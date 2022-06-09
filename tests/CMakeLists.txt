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

add_subdirectory(keywords)

add_test(NAME port COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         port.mop port.parameters try.txt aa.txt mol.in SETUP)

add_test(NAME crambin COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         Crambin_1SCF.mop Crambin_1SCF.den)

add_test(NAME mozyme COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         test_Lewis_for_Proteins.mop)

add_test(NAME xeno COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         Test_XENO.mop)

add_test(NAME hetero_group COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         Test_Hetero_group.mop)

add_test(NAME convert_edited_pdb COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         convert_an_edited_PDB_file.mop)

add_test(NAME convert_unedited_pdb COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         convert_an_un-edited_PDB_file.mop)

add_test(NAME ters_preserved COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         Check_that_TERs_are_preserved.mop Check_TERs.pdb)

add_test(NAME compare_proteins COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         Compare_two_protein_structures.mop Structure-1.txt Structure-2.txt)

add_test(NAME pdb1cbn COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         PDB1CBN.ent)

add_test(NAME pdb1ejg COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         PDB1EJG.ent)

add_test(NAME pdb_compare COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         Compare_PDB_1EJG_and_1CBN.mop PDB1CBN.ent PDB1EJG.ent)

add_test(NAME pdb_compare2 COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         Compare_PDB_Bacteriorhodopsin_1AP9_and_2AT9.mop Bacteriorhodopsin_1AP9.ent Bacteriorhodopsin_2AT9.ent)

add_test(NAME pdb_compare3 COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         Compare_PDB_Bacteriorhodopsin_1AP9_and_1C3W.mop Bacteriorhodopsin_1AP9.ent Bacteriorhodopsin_1C3W.ent)

add_test(NAME pdb_compare4 COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         Compare_PDB_Bacteriorhodopsin_1AP9_and_1BRX.mop Bacteriorhodopsin_1AP9.ent Bacteriorhodopsin_1BRX.ent)

add_test(NAME pdb_compare5 COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         Compare_PDB_Bacteriorhodopsin_1AP9_and_1BRR.mop Bacteriorhodopsin_1AP9.ent Bacteriorhodopsin_1BRR.ent)

add_test(NAME pdb_compare6 COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         Compare_PDB_Bacteriorhodopsin_1AP9_and_1AT9.mop ADD-H_Bacteriorhodopsin_1AP9.arc ADD-H_Bacteriorhodopsin_1AT9.arc)

add_test(NAME add-h COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         ADD-H_Bacteriorhodopsin_1AP9.mop Bacteriorhodopsin_1AP9.ent)

add_test(NAME add-h2 COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         ADD-H_Bacteriorhodopsin_1AT9.mop Bacteriorhodopsin_1AT9.ent)

add_test(NAME add-h3 COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         Au_Add-H.mop)

add_test(NAME compare_add-h COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         Compare_ADD-H_1EJG_and_1CBN.mop PDB1EJG.arc PDB1CBN.arc)

add_test(NAME locate_ts COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         Locate_TS.mop Cyclobutadiene_1.arc Cyclobutadiene_2.arc)

add_test(NAME ionize COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         IONIZE.mop Neutral.mop)

add_test(NAME salt COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         SALT.mop Neutral.mop)

add_test(NAME salt2 COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         SALT_2.mop)

add_test(NAME residues COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         SITE_Individual_residues.mop IONIZE.arc)

add_test(NAME ag20 COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         ag20.mop)

add_test(NAME ag20_cis_active COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         ag20_cis_active.mop)

add_test(NAME benzene COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         benzene.mop)

add_test(NAME co_cas COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         co_cas.mop)

add_test(NAME co_cis COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         co_cis.mop)

add_test(NAME co_cis_2 COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         co_cis_2.mop)

add_test(NAME co_cis_PM7 COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         co_cis_PM7.mop)

add_test(NAME co_cis_active COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         co_cis_active.mop)

add_test(NAME co_cis_solv COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         co_cis_solv.mop)

add_test(NAME co_cis_tdip COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         co_cis_tdip.mop)

add_test(NAME co_cis_trip COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         co_cis_trip.mop)

add_test(NAME co_extern COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         co_extern.mop param.txt)

add_test(NAME co_fail COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         co_fail.mop)

add_test(NAME co_field COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         co_field.mop)

add_test(NAME co_mrci COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         co_mrci.mop)

add_test(NAME co_mrci_maxci COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         co_mrci_maxci.mop)

add_test(NAME co_mrci_solv COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         co_mrci_solv.mop)

add_test(NAME co_pm7 COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         co_pm7.mop)

add_test(NAME co_pm7_solv COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         co_pm7_solv.mop)

add_test(NAME co_print COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         co_print.mop)

add_test(NAME no_cis_quar COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         no_cis_quar.mop)

add_test(NAME no_mrci COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         no_mrci.mop)

add_test(NAME noparam COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         noparam.mop)

add_test(NAME poly_c8 COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         poly_c8.mop)

add_test(NAME poly_c8_mr COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../mopac${CMAKE_EXECUTABLE_SUFFIX}
         poly_c8_mr.mop)
