# MAKPOL

This directory contains the MAKPOL program for generating supercells for use by MOPAC.
It is a useful tool for converging finite-size effects for  periodic calculations in MOPAC,
since there is no support in MOPAC for Brillouin zone sampling (all periodic calculations
are performed at the Gamma point). There is some overlap between the source code of MOPAC
and MAKPOL, but they are being kept separate for now (merge efforts are welcome).
