#Basic input/output comparison tests

add_test(NAME port COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../mopac
         port.mop port.parameters setup.iter try.txt aa.txt mol.in SETUP)

add_test(NAME crambin COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../mopac
         Crambin_1SCF.mop Crambin_1SCF.den)

add_test(NAME mozyme COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../mopac
         test_Lewis_for_Proteins.mop)

add_test(NAME xeno COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../mopac
         Test_XENO.mop)

add_test(NAME hetero_group COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../mopac
         Test_Hetero_group.mop)

add_test(NAME convert_edited_pdb COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../mopac
         convert_an_edited_PDB_file.mop)

add_test(NAME convert_unedited_pdb COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../mopac
         convert_an_un-edited_PDB_file.mop)

add_test(NAME ters_preserved COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../mopac
         Check_that_TERs_are_preserved.mop Check_TERs.pdb)

add_test(NAME compare_proteins COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../mopac
         Compare_two_protein_structures.mop Structure-1.txt Structure-2.txt)

add_test(NAME pdb1cbn COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../mopac
         PDB1CBN.ent)

add_test(NAME pdb1ejg COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../mopac
         PDB1EJG.ent)

add_test(NAME pdb_compare COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../mopac
         Compare_PDB_1EJG_and_1CBN.mop PDB1CBN.ent PDB1EJG.ent)

add_test(NAME pdb_compare2 COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../mopac
         Compare_PDB_Bacteriorhodopsin_1AP9_and_2AT9.mop Bacteriorhodopsin_1AP9.ent Bacteriorhodopsin_2AT9.ent)

add_test(NAME pdb_compare3 COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../mopac
         Compare_PDB_Bacteriorhodopsin_1AP9_and_1C3W.mop Bacteriorhodopsin_1AP9.ent Bacteriorhodopsin_1C3W.ent)

add_test(NAME pdb_compare4 COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../mopac
         Compare_PDB_Bacteriorhodopsin_1AP9_and_1BRX.mop Bacteriorhodopsin_1AP9.ent Bacteriorhodopsin_1BRX.ent)

add_test(NAME pdb_compare5 COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../mopac
         Compare_PDB_Bacteriorhodopsin_1AP9_and_1BRR.mop Bacteriorhodopsin_1AP9.ent Bacteriorhodopsin_1BRR.ent)

add_test(NAME pdb_compare6 COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../mopac
         Compare_PDB_Bacteriorhodopsin_1AP9_and_1AT9.mop ADD-H_Bacteriorhodopsin_1AP9.arc ADD-H_Bacteriorhodopsin_1AT9.arc)

add_test(NAME add-h COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../mopac
         ADD-H_Bacteriorhodopsin_1AP9.mop Bacteriorhodopsin_1AP9.ent)

add_test(NAME add-h2 COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../mopac
         ADD-H_Bacteriorhodopsin_1AT9.mop Bacteriorhodopsin_1AT9.ent)

add_test(NAME add-h3 COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../mopac
         Au_Add-H.mop)

add_test(NAME compare_add-h COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../mopac
         Compare_ADD-H_1EJG_and_1CBN.mop PDB1EJG.arc PDB1CBN.arc)

add_test(NAME locate_ts COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../mopac
         Locate_TS.mop Cyclobutadiene_1.arc Cyclobutadiene_2.arc)

add_test(NAME ionize COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../mopac
         IONIZE.mop Neutral.mop)

add_test(NAME salt COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../mopac
         SALT.mop Neutral.mop)

add_test(NAME salt2 COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../mopac
         SALT_2.mop)

add_test(NAME residues COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../mopac
         SITE_Individual_residues.mop IONIZE.arc)
