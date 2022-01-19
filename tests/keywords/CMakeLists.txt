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

# Basic tests for parsing/using each MOPAC keyword, including obsolete keywords & old syntax for backward compatibility

add_test(NAME key-ampersand COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         ampersand.mop)

add_test(NAME key-plus COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         plus.mop)

add_test(NAME key-plusplus COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         plusplus.mop)

add_test(NAME key-0SCF COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         0SCF.mop)

add_test(NAME key-1ELECTRON COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         1ELECTRON.mop)

add_test(NAME key-1SCF COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         1SCF.mop)

add_test(NAME key-A0 COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         A0.mop)

add_test(NAME key-AIGIN COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         AIGIN.mop)

add_test(NAME key-AIGOUT COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         AIGOUT.mop)

add_test(NAME key-ALLBONDS COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         ALLBONDS.mop)

add_test(NAME key-ALLVEC COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         ALLVEC.mop)

add_test(NAME key-AM1 COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         AM1.mop)

add_test(NAME key-ANGSTROMS COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         ANGSTROMS.mop)

add_test(NAME key-AUTOSYM COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         AUTOSYM.mop)

#add_test(NAME key-BANANA COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         BANANA.mop)

add_test(NAME key-BIRADICAL COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         BIRADICAL.mop)

#add_test(NAME key-BONDS COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         BONDS.mop)

add_test(NAME key-CdotAdotSdot COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         CdotAdotSdot.mop)

add_test(NAME key-CdotIdot COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         CdotIdot.mop)

add_test(NAME key-CdotIdotDdot COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         CdotIdotDdot.mop)

#add_test(NAME key-CARTAB COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         CARTAB.mop)

add_test(NAME key-CHARGE COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         CHARGE.mop)

#add_test(NAME key-CHARGES COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         CHARGES.mop)

#add_test(NAME key-CHARST COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         CHARST.mop)

add_test(NAME key-CIS COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         CIS.mop)

add_test(NAME key-CISD COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         CISD.mop)

add_test(NAME key-CISDT COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         CISDT.mop)

#add_test(NAME key-COMPFG COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         COMPFG.mop)

#add_test(NAME key-DCART COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         DCART.mop)

#add_test(NAME key-DEBUG COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         DEBUG.mop)

#add_test(NAME key-DEBUGPULAY COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         DEBUGPULAY.mop)

#add_test(NAME key-DERI1 COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         DERI1.mop)

#add_test(NAME key-DERI2 COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         DERI2.mop)

#add_test(NAME key-DERITR COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         DERITR.mop)

#add_test(NAME key-DERIV COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         DERIV.mop)

#add_test(NAME key-DERNVO COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         DERNVO.mop)

#add_test(NAME key-DFORCE COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         DFORCE.mop)

add_test(NAME key-DISP COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         DISP.mop)

add_test(NAME key-DOUBLET COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         DOUBLET.mop)

#add_test(NAME key-ECHO COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         ECHO.mop)

#add_test(NAME key-EIGEN COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         EIGEN.mop)

#add_test(NAME key-EIGS COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         EIGS.mop)

add_test(NAME key-EXCITED COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         EXCITED.mop)

#add_test(NAME key-FLEPO COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         FLEPO.mop)

#add_test(NAME key-FMAT COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         FMAT.mop)

#add_test(NAME key-FOCK COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         FOCK.mop)

#add_test(NAME key-FREQCY COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         FREQCY.mop)

#add_test(NAME key-GRADIENTS COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         GRADIENTS.mop)

#add_test(NAME key-HCORE COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         HCORE.mop)

#add_test(NAME key-HESSIAN COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         HESSIAN.mop)

add_test(NAME key-INDO COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         INDO.mop)

#add_test(NAME key-ITER COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         ITER.mop)

#add_test(NAME key-LARGE COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         LARGE.mop)

#add_test(NAME key-LINMIN COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         LINMIN.mop)

#add_test(NAME key-LOCAL COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         LOCAL.mop)

add_test(NAME key-MAXCI COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         MAXCI.mop)

add_test(NAME key-MECI COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         MECI.mop)

#add_test(NAME key-MINI COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         MINI.mop)

add_test(NAME key-MNDO COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         MNDO.mop)

add_test(NAME key-MNDOD COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         MNDOD.mop)

#add_test(NAME key-MOLDAT COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         MOLDAT.mop)

#add_test(NAME key-MOLSYM COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         MOLSYM.mop)

add_test(NAME key-MRCI COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         MRCI.mop)

add_test(NAME key-MS COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         MS.mop)

#add_test(NAME key-MULLIK COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         MULLIK.mop)

add_test(NAME key-NOANCI COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         NOANCI.mop)

add_test(NAME key-NONET COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         NONET.mop)

add_test(NAME key-NOSYM COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         NOSYM.mop)

#add_test(NAME key-NOXYZ COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         NOXYZ.mop)

add_test(NAME key-OCTET COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         OCTET.mop)

add_test(NAME key-OPEN COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         OPEN.mop)

add_test(NAME key-PECI COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         PECI.mop)

#add_test(NAME key-PKA COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         PKA.mop)

add_test(NAME key-PM3 COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         PM3.mop)

add_test(NAME key-PM6 COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         PM6.mop)

add_test(NAME key-PM6-D3 COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         PM6-D3.mop)

add_test(NAME key-PM6-D3H4 COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         PM6-D3H4.mop)

add_test(NAME key-PM6-D3H4X COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         PM6-D3H4X.mop)

add_test(NAME key-PM6-DHplus COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         PM6-DHplus.mop)

add_test(NAME key-PM6-DH2 COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         PM6-DH2.mop)

add_test(NAME key-PM6-DH2X COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         PM6-DH2X.mop)

add_test(NAME key-PM7 COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         PM7.mop)

add_test(NAME key-PM7-TS COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         PM7-TS.mop)

#add_test(NAME key-POWSQ COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         POWSQ.mop)

#add_test(NAME key-PRNT COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         PRNT.mop)

#add_test(NAME key-PRTCHAR COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         PRTCHAR.mop)

#add_test(NAME key-PRTINT COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         PRTINT.mop)

#add_test(NAME key-PRTMEP COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         PRTMEP.mop)

#add_test(NAME key-PRTXYZ COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         PRTXYZ.mop)

add_test(NAME key-QUARTET COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         QUARTET.mop)

add_test(NAME key-QUINTET COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         QUINTET.mop)

#add_test(NAME key-RABBIT COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         RABBIT.mop)

#add_test(NAME key-RAMA COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         RAMA.mop)

#add_test(NAME key-RE-LOCAL COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         RE-LOCAL.mop)

add_test(NAME key-RHF COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         RHF.mop)

add_test(NAME key-RM1 COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         RM1.mop)

add_test(NAME key-ROOT COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         ROOT.mop)

add_test(NAME key-SEPTET COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         SEPTET.mop)

add_test(NAME key-SEXTET COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         SEXTET.mop)

add_test(NAME key-SINGLET COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         SINGLET.mop)

add_test(NAME key-SNAP COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         SNAP.mop)

#add_test(NAME key-SPIN COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         SPIN.mop)

#add_test(NAME key-SUPER COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         SUPER.mop)

add_test(NAME key-SYMMETRY COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         SYMMETRY.mop)

#add_test(NAME key-SYMOIR COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         SYMOIR.mop)

#add_test(NAME key-SYMTRZ COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         SYMTRZ.mop)

add_test(NAME key-TDIP COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         TDIP.mop)

#add_test(NAME key-TIMES COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         TIMES.mop)

add_test(NAME key-TRIPLET COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         TRIPLET.mop)

add_test(NAME key-UHF COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         UHF.mop)

#add_test(NAME key-VECTORS COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
#         VECTORS.mop)

add_test(NAME key-WRTCI COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         WRTCI.mop)

add_test(NAME key-WRTCONF COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/../output_compare.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../../mopac${CMAKE_EXECUTABLE_SUFFIX}
         WRTCONF.mop)
