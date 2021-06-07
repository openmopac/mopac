# Basic tests for parsing/using each MOPAC keyword, including obsolete keywords & old syntax for backward compatibility

add_test(NAME key-ampersand COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         ampersand.mop)

add_test(NAME key-plus COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         plus.mop)

add_test(NAME key-plusplus COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         plusplus.mop)

add_test(NAME key-0SCF COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         0SCF.mop)

add_test(NAME key-1ELECTRON COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         1ELECTRON.mop)

add_test(NAME key-1SCF COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         1SCF.mop)

add_test(NAME key-A0 COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         A0.mop)

add_test(NAME key-AIGIN COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         AIGIN.mop)

add_test(NAME key-AIGOUT COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         AIGOUT.mop)

add_test(NAME key-ALLBONDS COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         ALLBONDS.mop)

add_test(NAME key-ALLVEC COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         ALLVEC.mop)

add_test(NAME key-AM1 COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         AM1.mop)

add_test(NAME key-ANGSTROMS COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         ANGSTROMS.mop)

add_test(NAME key-MNDO COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         MNDO.mop)

add_test(NAME key-MNDOD COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         MNDOD.mop)

add_test(NAME key-PM3 COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         PM3.mop)

add_test(NAME key-PM6 COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         PM6.mop)

add_test(NAME key-PM6-D3 COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         PM6-D3.mop)

add_test(NAME key-PM6-D3H4 COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         PM6-D3H4.mop)

add_test(NAME key-PM6-D3H4X COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         PM6-D3H4X.mop)

add_test(NAME key-PM6-DHplus COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         PM6-DHplus.mop)

add_test(NAME key-PM6-DH2 COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         PM6-DH2.mop)

add_test(NAME key-PM6-DH2X COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         PM6-DH2X.mop)

add_test(NAME key-PM7 COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         PM7.mop)

add_test(NAME key-PM7-TS COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         PM7-TS.mop)

add_test(NAME key-RM1 COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         RM1.mop)
