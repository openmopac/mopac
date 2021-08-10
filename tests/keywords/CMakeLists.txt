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

add_test(NAME key-AUTOSYM COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         AUTOSYM.mop)

add_test(NAME key-BIRADICAL COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         BIRADICAL.mop)

add_test(NAME key-CdotAdotSdot COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         CdotAdotSdot.mop)

add_test(NAME key-CdotIdot COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         CdotIdot.mop)

add_test(NAME key-CdotIdotDdot COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         CdotIdotDdot.mop)

add_test(NAME key-CHARGE COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         CHARGE.mop)

add_test(NAME key-CIS COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         CIS.mop)

add_test(NAME key-CISD COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         CISD.mop)

add_test(NAME key-CISDT COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         CISDT.mop)

add_test(NAME key-DISP COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         DISP.mop)

add_test(NAME key-DOUBLET COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         DOUBLET.mop)

add_test(NAME key-EXCITED COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         EXCITED.mop)

add_test(NAME key-INDO COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         INDO.mop)

add_test(NAME key-MAXCI COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         MAXCI.mop)

add_test(NAME key-MECI COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         MECI.mop)

add_test(NAME key-MNDO COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         MNDO.mop)

add_test(NAME key-MNDOD COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         MNDOD.mop)

add_test(NAME key-MRCI COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         MRCI.mop)

add_test(NAME key-MS COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         MS.mop)

add_test(NAME key-NOANCI COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         NOANCI.mop)

add_test(NAME key-NONET COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         NONET.mop)

add_test(NAME key-NOSYM COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         NOSYM.mop)

add_test(NAME key-OCTET COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         OCTET.mop)

add_test(NAME key-OPEN COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         OPEN.mop)

add_test(NAME key-PECI COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         PECI.mop)

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

add_test(NAME key-QUARTET COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         QUARTET.mop)

add_test(NAME key-QUINTET COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         QUINTET.mop)

add_test(NAME key-RHF COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         RHF.mop)

add_test(NAME key-RM1 COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         RM1.mop)

add_test(NAME key-ROOT COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         ROOT.mop)

add_test(NAME key-SEPTET COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         SEPTET.mop)

add_test(NAME key-SEXTET COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         SEXTET.mop)

add_test(NAME key-SINGLET COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         SINGLET.mop)

add_test(NAME key-SNAP COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         SNAP.mop)

add_test(NAME key-SYMMETRY COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         SYMMETRY.mop)

add_test(NAME key-TDIP COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         TDIP.mop)

add_test(NAME key-TRIPLET COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         TRIPLET.mop)

add_test(NAME key-UHF COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         UHF.mop)

add_test(NAME key-WRTCI COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         WRTCI.mop)

add_test(NAME key-WRTCONF COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ../../mopac
         WRTCONF.mop)
