! Molecular Orbital PACkage (MOPAC)
! Copyright (C) 2021, Virginia Polytechnic Institute and State University
!
! MOPAC is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! MOPAC is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

      module journal_references_C
      character, dimension(107,10) :: allref*100
! 2nd index denotes the model:
! 1 = MNDO
! 2 = PM6
! 3 = AM1
! 4 = PM3
! 5 = MNDO/d
! 6 = RM1
! 7 = PM7
! 8 = RM1 sparkles
! 9 = PM6-ORG
! 10 = PM8

 data  allref( 1,6)/ '  H: (RM1): G.B. ROCHA, R.O. FREIRE, A.M.SIMAS, J.J.P. STEWART, J. COMP. CHEM. 27, 1101 (2006)'/
 data  allref( 6,6)/ '  C: (RM1): G.B. ROCHA, R.O. FREIRE, A.M.SIMAS, J.J.P. STEWART, J. COMP. CHEM. 27, 1101 (2006)'/
 data  allref( 7,6)/ '  N: (RM1): G.B. ROCHA, R.O. FREIRE, A.M.SIMAS, J.J.P. STEWART, J. COMP. CHEM. 27, 1101 (2006)'/
 data  allref( 8,6)/ '  O: (RM1): G.B. ROCHA, R.O. FREIRE, A.M.SIMAS, J.J.P. STEWART, J. COMP. CHEM. 27, 1101 (2006)'/
 data  allref( 9,6)/ '  F: (RM1): G.B. ROCHA, R.O. FREIRE, A.M.SIMAS, J.J.P. STEWART, J. COMP. CHEM. 27, 1101 (2006)'/
 data  allref(15,6)/ '  P: (RM1): G.B. ROCHA, R.O. FREIRE, A.M.SIMAS, J.J.P. STEWART, J. COMP. CHEM. 27, 1101 (2006)'/
 data  allref(16,6)/ '  S: (RM1): G.B. ROCHA, R.O. FREIRE, A.M.SIMAS, J.J.P. STEWART, J. COMP. CHEM. 27, 1101 (2006)'/
 data  allref(17,6)/ ' Cl: (RM1): G.B. ROCHA, R.O. FREIRE, A.M.SIMAS, J.J.P. STEWART, J. COMP. CHEM. 27, 1101 (2006)'/
 data  allref(35,6)/ ' Br: (RM1): G.B. ROCHA, R.O. FREIRE, A.M.SIMAS, J.J.P. STEWART, J. COMP. CHEM. 27, 1101 (2006)'/
 data  allref(53,6)/ '  I: (RM1): G.B. ROCHA, R.O. FREIRE, A.M.SIMAS, J.J.P. STEWART, J. COMP. CHEM. 27, 1101 (2006)'/

 data  allref(11,5)/ " Na: (MNDO/d): W.THIEL AND A.A.VOITYUK, J. PHYS. CHEM., 100, 616 (1996)                      "/
 data  allref(12,5)/ " Mg: (MNDO/d): W.THIEL AND A.A.VOITYUK, J. PHYS. CHEM., 100, 616 (1996)                      "/
 data  allref(13,5)/ " Al: (MNDO/d): W.THIEL AND A.A.VOITYUK, J. PHYS. CHEM., 100, 616 (1996)                      "/
 data  allref(14,5)/ " Si: (MNDO/d): W.THIEL AND A.A.VOITYUK, J. PHYS. CHEM., 100, 616 (1996)                      "/
 data  allref(15,5)/ "  P: (MNDO/d): W.THIEL AND A.A.VOITYUK, J. PHYS. CHEM., 100, 616 (1996)                      "/
 data  allref(16,5)/ "  S: (MNDO/d): W.THIEL AND A.A.VOITYUK, J. PHYS. CHEM., 100, 616 (1996)                      "/
 data  allref(17,5)/ " Cl: (MNDO/d): W.THIEL AND A.A.VOITYUK, J. PHYS. CHEM., 100, 616 (1996)                      "/
 data  allref(30,5)/ " Zn: (MNDO/d): W.THIEL AND A.A.VOITYUK, J. PHYS. CHEM., 100, 616 (1996)                      "/
 data  allref(35,5)/ " Br: (MNDO/d): W.THIEL AND A.A.VOITYUK, J. PHYS. CHEM., 100, 616 (1996)                      "/
 data  allref(48,5)/ " Cd: (MNDO/d): W.THIEL AND A.A.VOITYUK, J. PHYS. CHEM., 100, 616 (1996)                      "/
 data  allref(53,5)/ "  I: (MNDO/d): W.THIEL AND A.A.VOITYUK, J. PHYS. CHEM., 100, 616 (1996)                      "/
 data  allref(80,5)/ " Hg: (MNDO/d): W.THIEL AND A.A.VOITYUK, J. PHYS. CHEM., 100, 616 (1996)                      "/
!
!  References for the RM1 Lanthanides
!
 data  allref(57,6)/ " La: (RM1): Dutra, J.D.L., Filho, M.A.M., et al., PLOS ONE, 2015, DOI: 10.1371/journal.pone.0124372 "/
 data  allref(58,6)/ " Ce: (RM1): Dutra, J.D.L., Filho, M.A.M., et al., PLOS ONE, 2015, DOI: 10.1371/journal.pone.0124372 "/
 data  allref(59,6)/ " Pr: (RM1): Dutra, J.D.L., Filho, M.A.M., et al., PLOS ONE, 2015, DOI: 10.1371/journal.pone.0124372 "/
 data  allref(60,6)/ " Nd: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., RSC Adv., 2015, DOI: 10.1039/C4RA12682C "/
 data  allref(61,6)/ " Pm: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., RSC Adv., 2015, DOI: 10.1039/C4RA12682C "/
 data  allref(62,6)/ " Sm: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., RSC Adv., 2015, DOI: 10.1039/C4RA12682C "/
 data  allref(63,6)/ " Eu: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., JCTC, 2014, DOI: 10.1021/ct400909w "/
 data  allref(64,6)/ " Gd: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., JCTC, 2014, DOI: 10.1021/ct400909w "/
 data  allref(65,6)/ " Tb: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., JCTC, 2014, DOI: 10.1021/ct400909w "/
 data  allref(66,6)/ " Dy: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., PLOS ONE, 2014, DOI: 10.1371/journal.pone.0086376 "/
 data  allref(67,6)/ " Ho: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., PLOS ONE, 2014, DOI: 10.1371/journal.pone.0086376 "/
 data  allref(68,6)/ " Er: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., PLOS ONE, 2014, DOI: 10.1371/journal.pone.0086376 "/
 data  allref(69,6)/ " Tm: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., PLOS ONE, 2016, DOI: 10.1371/journal.pone.0154500 "/
 data  allref(70,6)/ " Yb: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., PLOS ONE, 2016, DOI: 10.1371/journal.pone.0154500 "/
 data  allref(71,6)/ " Lu: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., PLOS ONE, 2016, DOI: 10.1371/journal.pone.0154500 "/
!
!  References for the RM1 Lanthanides as sparkles
!
 data  allref(57,8)/ " La: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., RSC Adv., 2013, DOI: 10.1039/C3RA41406J "/
 data  allref(58,8)/ " Ce: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., RSC Adv., 2013, DOI: 10.1039/C3RA41406J "/
 data  allref(59,8)/ " Pr: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., RSC Adv., 2013, DOI: 10.1039/C3RA41406J "/
 data  allref(60,8)/ " Nd: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., RSC Adv., 2013, DOI: 10.1039/C3RA41406J "/
 data  allref(61,8)/ " Pm: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., RSC Adv., 2013, DOI: 10.1039/C3RA41406J "/
 data  allref(62,8)/ " Sm: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., RSC Adv., 2013, DOI: 10.1039/C3RA41406J "/
 data  allref(63,8)/ " Eu: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., RSC Adv., 2013, DOI: 10.1039/C3RA41406J "/
 data  allref(64,8)/ " Gd: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., RSC Adv., 2013, DOI: 10.1039/C3RA41406J "/
 data  allref(65,8)/ " Tb: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., RSC Adv., 2013, DOI: 10.1039/C3RA41406J "/
 data  allref(66,8)/ " Dy: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., RSC Adv., 2013, DOI: 10.1039/C3RA41406J "/
 data  allref(67,8)/ " Ho: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., RSC Adv., 2013, DOI: 10.1039/C3RA41406J "/
 data  allref(68,8)/ " Er: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., RSC Adv., 2013, DOI: 10.1039/C3RA41406J "/
 data  allref(69,8)/ " Tm: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., RSC Adv., 2013, DOI: 10.1039/C3RA41406J "/
 data  allref(70,8)/ " Yb: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., RSC Adv., 2013, DOI: 10.1039/C3RA41406J "/
 data  allref(71,8)/ " Lu: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., RSC Adv., 2013, DOI: 10.1039/C3RA41406J "/
!
!  References for the PM6 Lanthanides
!
 data  allref(57,2)/ " La: (PM6): Freire R.O. and Simas, A.M.  J. Chem. Theory and Comput., 2010, 6, 2019-2023     "/
 data  allref(58,2)/ " Ce: (PM6): Freire R.O. and Simas, A.M.  J. Chem. Theory and Comput., 2010, 6, 2019-2023     "/
 data  allref(59,2)/ " Pr: (PM6): Freire R.O. and Simas, A.M.  J. Chem. Theory and Comput., 2010, 6, 2019-2023     "/
 data  allref(60,2)/ " Nd: (PM6): Freire R.O. and Simas, A.M.  J. Chem. Theory and Comput., 2010, 6, 2019-2023     "/
 data  allref(61,2)/ " Pm: (PM6): Freire R.O. and Simas, A.M.  J. Chem. Theory and Comput., 2010, 6, 2019-2023     "/
 data  allref(62,2)/ " Sm: (PM6): Freire R.O. and Simas, A.M.  J. Chem. Theory and Comput., 2010, 6, 2019-2023     "/
 data  allref(63,2)/ " Eu: (PM6): Freire R.O. and Simas, A.M.  J. Chem. Theory and Comput., 2010, 6, 2019-2023     "/
 data  allref(64,2)/ " Gd: (PM6): Freire R.O. and Simas, A.M.  J. Chem. Theory and Comput., 2010, 6, 2019-2023     "/
 data  allref(65,2)/ " Tb: (PM6): Freire R.O. and Simas, A.M.  J. Chem. Theory and Comput., 2010, 6, 2019-2023     "/
 data  allref(66,2)/ " Dy: (PM6): Freire R.O. and Simas, A.M.  J. Chem. Theory and Comput., 2010, 6, 2019-2023     "/
 data  allref(67,2)/ " Ho: (PM6): Freire R.O. and Simas, A.M.  J. Chem. Theory and Comput., 2010, 6, 2019-2023     "/
 data  allref(68,2)/ " Er: (PM6): Freire R.O. and Simas, A.M.  J. Chem. Theory and Comput., 2010, 6, 2019-2023     "/
 data  allref(69,2)/ " Tm: (PM6): Freire R.O. and Simas, A.M.  J. Chem. Theory and Comput., 2010, 6, 2019-2023     "/
 data  allref(70,2)/ " Yb: (PM6): Freire R.O. and Simas, A.M.  J. Chem. Theory and Comput., 2010, 6, 2019-2023     "/
 data  allref(71,2)/ " Lu: (PM6): Freire R.O. and Simas, A.M.  J. Chem. Theory and Comput., 2010, 6, 2019-2023     "/
!
!
!  References for the PM7 Lanthanides
!
!
 data  allref(57,7)/ " La: (PM7): Dutra, J.D.L., Filho, M.A.M., et al J. Chem. Theory and Comput., 2013, 9, 3333-3341 "/
 data  allref(58,7)/ " Ce: (PM7): Dutra, J.D.L., Filho, M.A.M., et al J. Chem. Theory and Comput., 2013, 9, 3333-3341 "/
 data  allref(59,7)/ " Pr: (PM7): Dutra, J.D.L., Filho, M.A.M., et al J. Chem. Theory and Comput., 2013, 9, 3333-3341 "/
 data  allref(60,7)/ " Nd: (PM7): Dutra, J.D.L., Filho, M.A.M., et al J. Chem. Theory and Comput., 2013, 9, 3333-3341 "/
 data  allref(61,7)/ " Pm: (PM7): Dutra, J.D.L., Filho, M.A.M., et al J. Chem. Theory and Comput., 2013, 9, 3333-3341 "/
 data  allref(62,7)/ " Sm: (PM7): Dutra, J.D.L., Filho, M.A.M., et al J. Chem. Theory and Comput., 2013, 9, 3333-3341 "/
 data  allref(63,7)/ " Eu: (PM7): Dutra, J.D.L., Filho, M.A.M., et al J. Chem. Theory and Comput., 2013, 9, 3333-3341 "/
 data  allref(64,7)/ " Gd: (PM7): Dutra, J.D.L., Filho, M.A.M., et al J. Chem. Theory and Comput., 2013, 9, 3333-3341 "/
 data  allref(65,7)/ " Tb: (PM7): Dutra, J.D.L., Filho, M.A.M., et al J. Chem. Theory and Comput., 2013, 9, 3333-3341 "/
 data  allref(66,7)/ " Dy: (PM7): Dutra, J.D.L., Filho, M.A.M., et al J. Chem. Theory and Comput., 2013, 9, 3333-3341 "/
 data  allref(67,7)/ " Ho: (PM7): Dutra, J.D.L., Filho, M.A.M., et al J. Chem. Theory and Comput., 2013, 9, 3333-3341 "/
 data  allref(68,7)/ " Er: (PM7): Dutra, J.D.L., Filho, M.A.M., et al J. Chem. Theory and Comput., 2013, 9, 3333-3341 "/
 data  allref(69,7)/ " Tm: (PM7): Dutra, J.D.L., Filho, M.A.M., et al J. Chem. Theory and Comput., 2013, 9, 3333-3341 "/
 data  allref(70,7)/ " Yb: (PM7): Dutra, J.D.L., Filho, M.A.M., et al J. Chem. Theory and Comput., 2013, 9, 3333-3341 "/
 data  allref(71,7)/ " Lu: (PM7): Dutra, J.D.L., Filho, M.A.M., et al J. Chem. Theory and Comput., 2013, 9, 3333-3341 "/
!
! MNDO
!
 data  allref(  1,1)/ '  H: (MNDO): M.J.S. DEWAR, W. THIEL, J. AM. CHEM. SOC., 99, 4899, (1977)        '/
 data  allref(  2,1)/ ' He: (MNDO):                                                                    '/
 data  allref(  3,1)/ ' Li: (MNDO): TAKEN FROM MNDOC BY W.THIEL (QCPE NO.438, V. 2, P.63, (1982))      '/
 data  allref(  4,1)/ ' Be: (MNDO): M.J.S. DEWAR, H.S. RZEPA, J. AM. CHEM. SOC., 100, 777, (1978)      '/
 data  allref(  5,1)/ '  B: (MNDO): M.J.S. DEWAR, M.L. MCKEE, J. AM. CHEM. SOC., 99, 5231, (1977)      '/
 data  allref(  6,1)/ '  C: (MNDO): M.J.S. DEWAR, W. THIEL, J. AM. CHEM. SOC., 99, 4899, (1977)        '/
 data  allref(  7,1)/ '  N: (MNDO): M.J.S. DEWAR, W. THIEL, J. AM. CHEM. SOC., 99, 4899, (1977)        '/
 data  allref(  8,1)/ '  O: (MNDO): M.J.S. DEWAR, W. THIEL, J. AM. CHEM. SOC., 99, 4899, (1977)        '/
 data  allref(  9,1)/ '  F: (MNDO): M.J.S. DEWAR, H.S. RZEPA, J. AM. CHEM. SOC., 100, 777, (1978)      '/
 data  allref( 10,1)/ ' Ne: (MNDO):                                                                    '/
 data  allref( 11,1)/ ' Na: (MNDO): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                 '/
 data  allref( 12,1)/ ' Mg: (MNDO): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                 '/
 data  allref( 13,1)/ ' Al: (MNDO): L.P. DAVIS, ET.AL.  J. COMP. CHEM., 2, 433, (1981)                 '/
 data  allref( 14,1)/ ' Si: (MNDO): M.J.S.DEWAR, ET. AL. ORGANOMETALLICS  5, 375 (1986)                '/
 data  allref( 15,1)/ '  P: (MNDO): M.J.S.DEWAR, M.L.MCKEE, H.S. RZEPA,J. AM. CHEM. SOC., 100 3607 1978'/
 data  allref( 16,1)/ '  S: (MNDO): M.J.S.DEWAR, C.H. REYNOLDS, J. COMP. CHEM. 7, 140-143 (1986)       '/
 data  allref( 17,1)/ ' Cl: (MNDO): M.J.S.DEWAR, H.S.RZEPA, J. COMP. CHEM., 4, 158, (1983)             '/
 data  allref( 18,1)/ ' Ar: (MNDO):                                                                    '/
 data  allref( 19,1)/ '  K: (MNDO): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                 '/
 data  allref( 20,1)/ ' Ca: (MNDO): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                 '/
 data  allref( 30,1)/ ' Zn: (MNDO): M.J.S. DEWAR, K.M. MERZ, ORGANOMETALLICS, 5, 1494-1496 (1986)      '/
 data  allref( 31,1)/ ' Ga: (MNDO): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                 '/
 data  allref( 32,1)/ ' Ge: (MNDO): M.J.S.DEWAR, G.L.GRADY, E.F.HEALY,ORGANOMETALLICS 6 186-189 (1987) '/
 data  allref( 33,1)/ ' As: (MNDO): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                 '/
 data  allref( 34,1)/ ' Se: (MNDO): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                 '/
 data  allref( 35,1)/ ' Br: (MNDO): M.J.S.DEWAR, E.F. HEALY, J. COMP. CHEM., 4, 542, (1983)            '/
 data  allref( 36,1)/ ' Kr: (MNDO):                                                                    '/
 data  allref( 37,1)/ ' Rb: (MNDO): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                 '/
 data  allref( 38,1)/ ' Sr: (MNDO): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                 '/
 data  allref( 49,1)/ ' In: (MNDO): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                 '/
 data  allref( 50,1)/ ' Sn: (MNDO): M.J.S.DEWAR,G.L.GRADY,J.J.P.STEWART, J.AM.CHEM.SOC.,106 6771 (1984)'/
 data  allref( 51,1)/ ' Sb: (MNDO): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                 '/
 data  allref( 52,1)/ ' Te: (MNDO): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                 '/
 data  allref( 53,1)/ '  I: (MNDO): M.J.S.DEWAR, E.F. HEALY, J.J.P. STEWART, J.COMP.CHEM., 5,358,(1984)'/
 data  allref( 54,1)/ ' Xe: (MNDO):                                                                    '/
 data  allref( 55,1)/ ' Cs: (MNDO): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                 '/
 data  allref( 56,1)/ ' Ba: (MNDO): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                 '/
 data  allref( 80,1)/ ' Hg: (MNDO): M.J.S.DEWAR,  ET. AL. ORGANOMETALLICS 4, 1964, (1985)              '/
 data  allref( 81,1)/ ' Tl: (MNDO): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                 '/
 data  allref( 82,1)/ ' Pb: (MNDO): M.J.S.DEWAR, ET.AL ORGANOMETALLICS 4, 1973-1980 (1985)             '/
 data  allref( 83,1)/ ' Bi: (MNDO): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                 '/
 data  allref(102,1)/ ' Cb: (MNDO): Capped Bond  (Hydrogen-like, takes on a  zero charge.)             '/
 data  allref(104,1)/ '  +: (MNDO): Sparkle with charge of +1                                          '/
 data  allref(106,1)/ '  -: (MNDO): Sparkle with charge of -1                                          '/
!
!                                          AM1
!
 data  allref(  1,3)/ '  H: (AM1): M.J.S. DEWAR ET AL, J. AM. CHEM. SOC. 107 3902-3909 (1985)          '/
 data  allref(  2,3)/ ' He: (AM1):                                                                     '/
 data  allref(  3,3)/ ' Li: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                  '/
 data  allref(  4,3)/ ' Be: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                  '/
 data  allref(  5,3)/ '  B: (AM1): M.J.S. DEWAR ET AL, ORGANOMETALLICS, 7, 513-521 (1988)              '/
 data  allref(  6,3)/ '  C: (AM1): M.J.S. DEWAR ET AL, J. AM. CHEM. SOC. 107 3902-3909 (1985)          '/
 data  allref(  7,3)/ '  N: (AM1): M.J.S. DEWAR ET AL, J. AM. CHEM. SOC. 107 3902-3909 (1985)          '/
 data  allref(  8,3)/ '  O: (AM1): M.J.S. DEWAR ET AL, J. AM. CHEM. SOC. 107 3902-3909 (1985)          '/
 data  allref(  9,3)/ '  F: (AM1): M.J.S. DEWAR AND E. G. ZOEBISCH, THEOCHEM, 180, 1 (1988).           '/
 data  allref( 10,3)/ ' Ne: (AM1):                                                                     '/
 data  allref( 11,3)/ ' Na: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                  '/
 data  allref( 12,3)/ ' Mg: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                  '/
 data  allref( 13,3)/ ' Al: (AM1): M. J. S. Dewar, A. J. Holder, Organometallics, 9, 508-511 (1990).   '/
 data  allref( 14,3)/ ' Si: (AM1): M.J.S.DEWAR, C. JIE, ORGANOMETALLICS, 6, 1486-1490 (1987).          '/
 data  allref( 15,3)/ '  P: (AM1): M.J.S.DEWAR, JIE, C, THEOCHEM, 187, 1 (1989)                        '/
 data  allref( 16,3)/ '  S: (AM1): M.J.S. DEWAR, Y-C YUAN, INORGANIC CHEMISTRY, 29, 3881:3890,  (1990) '/
 data  allref( 17,3)/ ' Cl: (AM1): M.J.S. DEWAR AND E. G. ZOEBISCH, THEOCHEM, 180, 1 (1988).           '/
 data  allref( 18,3)/ ' Ar: (AM1):                                                                     '/
 data  allref( 19,3)/ '  K: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                  '/
 data  allref( 20,3)/ ' Ca: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                  '/
 data  allref( 30,3)/ ' Zn: (AM1): M.J.S. DEWAR, K.M. MERZ, ORGANOMETALLICS, 7, 522-524 (1988).        '/
 data  allref( 31,3)/ ' Ga: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                  '/
 data  allref( 32,3)/ ' Ge: (AM1): M.J.S.Dewar and C.Jie, Organometallics, 8, 1544, (1989)             '/
 data  allref( 33,3)/ ' As: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                  '/
 data  allref( 34,3)/ ' Se: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                  '/
 data  allref( 35,3)/ ' Br: (AM1): M.J.S. DEWAR AND E. G. ZOEBISCH, THEOCHEM, 180, 1 (1988).           '/
 data  allref( 36,3)/ ' Kr: (AM1):                                                                     '/
 data  allref( 37,3)/ ' Rb: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                  '/
 data  allref( 38,3)/ ' Sr: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                  '/
 data  allref( 42,3)/ ' Mo: (AM1): A.A. Voityuk, N. Roesch J. Phys. Chem. A 104, 4089 (2000).          '/
 data  allref( 49,3)/ ' In: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                  '/
 data  allref( 50,3)/ ' Sn: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                  '/
 data  allref( 51,3)/ ' Sb: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                  '/
 data  allref( 52,3)/ ' Te: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                  '/
 data  allref( 53,3)/ '  I: (AM1): M.J.S. DEWAR AND E. G. ZOEBISCH, THEOCHEM, 180, 1 (1988).           '/
 data  allref( 54,3)/ ' Xe: (AM1):                                                                     '/
 data  allref( 55,3)/ ' Cs: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                  '/
 data  allref( 56,3)/ ' Ba: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                  '/
 data  allref( 57,3)/ ' La: (AM1): R.O. FREIRE, ET AL, J. PHYS. CHEM. A 110 (17) (2006) 5897.          '/
 data  allref( 58,3)/ ' Ce: (AM1): R.O. FREIRE, ET AL, J. ORGANOMETALLIC CHEMISTRY 691 (2006) 2584.    '/
 data  allref( 59,3)/ ' Pr: (AM1): R.O. FREIRE, ET AL, J. ORGANOMETALLIC CHEMISTRY 690 (2005) 4099.    '/
 data  allref( 60,3)/ ' Nd: (AM1): C.C. BASTOS, ET AL, J. PHOTOCHEM. PHOTOBIOL. A 177 (2006) 225.      '/
 data  allref( 61,3)/ ' Pm: (AM1): R.O. FREIRE, ET AL, J. CHEM. THEORY COMPUT. 2 (2006) 64.            '/
 data  allref( 62,3)/ ' Sm: (AM1): R.O. FREIRE, ET AL, J. CHEM. THEORY COMPUT. 2 (2006) 64.            '/
 data  allref( 63,3)/ ' Eu: (AM1): R.O. FREIRE, G.B. ROCHA, A.M. SIMAS, INORG. CHEM.44 (2005) 3299.    '/
 data  allref( 64,3)/ ' Gd: (AM1): R.O. FREIRE, G.B. ROCHA, A.M. SIMAS, INORG. CHEM.44 (2005) 3299.    '/
 data  allref( 65,3)/ ' Tb: (AM1): R.O. FREIRE, G.B. ROCHA, A.M. SIMAS, INORG. CHEM.44 (2005) 3299.    '/
 data  allref( 66,3)/ ' Dy: (AM1): N.B. DA COSTA JR, ET AL, INORG. CHEM. COMM. 8 (2005) 831.           '/
 data  allref( 67,3)/ ' Ho: (AM1): N.B. DA COSTA JR, ET AL, POLYHEDRON 24 (2005) 3046.                 '/
 data  allref( 68,3)/ ' Er: (AM1): R.O. FREIRE, ET AL, J. ORGANOMETALLIC CHEMISTRY 691 (2006) 2584.    '/
 data  allref( 69,3)/ ' Tm: (AM1): R.O. FREIRE, ET AL, CHEM. PHYS. LETT., 411 (2005) 61.               '/
 data  allref( 70,3)/ ' Yb: (AM1): R.O. FREIRE, G.B. ROCHA, A.M. SIMAS, J. COMP. CHEM., 26 (2005) 1524.'/
 data  allref( 71,3)/ ' Lu: (AM1): R.O. FREIRE, ET AL, J. PHYS. CHEM. A 110 (17) (2006) 5897.          '/
 data  allref( 80,3)/ ' Hg: (AM1): M.J.S.Dewar and C.Jie, Organometallics 8, 1547, (1989)              '/
 data  allref( 81,3)/ ' Tl: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                  '/
 data  allref( 82,3)/ ' Pb: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                  '/
 data  allref( 83,3)/ ' Bi: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                  '/
 data  allref( 90,3)/ ' Th: (AM1): Rocha, G.B., D.Sc. Thesis, Dept. Quim. U. Pernambuco, Brazil  (2002)'/
 data  allref(102,3)/ ' Cb: (AM1): Capped Bond  (Hydrogen-like, takes on a zero charge.)               '/
 data  allref(104,3)/ '  +: (AM1): Sparkle with charge of +1                                           '/
 data  allref(106,3)/ '  -: (AM1): Sparkle with charge of -1                                           '/
!
!                                          PM3
!
 data  allref(  1,4)/ '  H: (PM3): J. J. P. STEWART, J. COMP. CHEM. 10, 209 (1989).                    '/
 data  allref(  2,4)/ ' He: (PM3):                                                                     '/
 data  allref(  3,4)/ ' Li: (PM3): E. ANDERS, R. KOCH, P. FREUNSCHT, J. COMP. CHEM 14 1301-1312  1993  '/
 data  allref(  4,4)/ ' Be: (PM3): J. J. P. STEWART, J. COMP. CHEM. 10, 209 (1989).                    '/
 data  allref(  5,4)/ '  B: (PM3): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                  '/
 data  allref(  6,4)/ '  C: (PM3): J. J. P. STEWART, J. COMP. CHEM. 10, 209 (1989).                    '/
 data  allref(  7,4)/ '  N: (PM3): J. J. P. STEWART, J. COMP. CHEM. 10, 209 (1989).                    '/
 data  allref(  8,4)/ '  O: (PM3): J. J. P. STEWART, J. COMP. CHEM. 10, 209 (1989).                    '/
 data  allref(  9,4)/ '  F: (PM3): J. J. P. STEWART, J. COMP. CHEM. 10, 209 (1989).                    '/
 data  allref( 10,4)/ ' Ne: (PM3):                                                                     '/
 data  allref( 11,4)/ ' Na: (PM3): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                  '/
 data  allref( 12,4)/ ' Mg: (PM3): J. J. P. STEWART, J. COMP. CHEM. 12, 320-341 (1991).                '/
 data  allref( 13,4)/ ' Al: (PM3): J. J. P. STEWART, J. COMP. CHEM. 10, 209 (1989).                    '/
 data  allref( 14,4)/ ' Si: (PM3): J. J. P. STEWART, J. COMP. CHEM. 10, 209 (1989).                    '/
 data  allref( 15,4)/ '  P: (PM3): J. J. P. STEWART, J. COMP. CHEM. 10, 209 (1989).                    '/
 data  allref( 16,4)/ '  S: (PM3): J. J. P. STEWART, J. COMP. CHEM. 10, 209 (1989).                    '/
 data  allref( 17,4)/ ' Cl: (PM3): J. J. P. STEWART, J. COMP. CHEM. 10, 209 (1989).                    '/
 data  allref( 18,4)/ ' Ar: (PM3):                                                                     '/
 data  allref( 19,4)/ '  K: (PM3): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                  '/
 data  allref( 20,4)/ ' Ca: (PM3): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                  '/
 data  allref( 30,4)/ ' Zn: (PM3): J. J. P. STEWART, J. COMP. CHEM. 12, 320-341 (1991).                '/
 data  allref( 31,4)/ ' Ga: (PM3): J. J. P. STEWART, J. COMP. CHEM. 12, 320-341 (1991).                '/
 data  allref( 32,4)/ ' Ge: (PM3): J. J. P. STEWART, J. COMP. CHEM. 12, 320-341 (1991).                '/
 data  allref( 33,4)/ ' As: (PM3): J. J. P. STEWART, J. COMP. CHEM. 12, 320-341 (1991).                '/
 data  allref( 34,4)/ ' Se: (PM3): J. J. P. STEWART, J. COMP. CHEM. 12, 320-341 (1991).                '/
 data  allref( 35,4)/ ' Br: (PM3): J. J. P. STEWART, J. COMP. CHEM. 10, 209 (1989).                    '/
 data  allref( 36,4)/ ' Kr: (PM3):                                                                     '/
 data  allref( 37,4)/ ' Rb: (PM3): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                  '/
 data  allref( 38,4)/ ' Sr: (PM3): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                  '/
 data  allref( 48,4)/ ' Cd: (PM3): J. J. P. STEWART, J. COMP. CHEM. 12, 320-341 (1991).                '/
 data  allref( 49,4)/ ' In: (PM3): J. J. P. STEWART, J. COMP. CHEM. 12, 320-341 (1991).                '/
 data  allref( 50,4)/ ' Sn: (PM3): J. J. P. STEWART, J. COMP. CHEM. 12, 320-341 (1991).                '/
 data  allref( 51,4)/ ' Sb: (PM3): J. J. P. STEWART, J. COMP. CHEM. 12, 320-341 (1991).                '/
 data  allref( 52,4)/ ' Te: (PM3): J. J. P. STEWART, J. COMP. CHEM. 12, 320-341 (1991).                '/
 data  allref( 53,4)/ '  I: (PM3): J. J. P. STEWART, J. COMP. CHEM. 10, 209 (1989).                    '/
 data  allref( 54,4)/ ' Xe: (PM3):                                                                     '/
 data  allref( 55,4)/ ' Cs: (PM3): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                  '/
 data  allref( 56,4)/ ' Ba: (PM3): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)                  '/
 data  allref( 57,4)/ ' La: (PM3): R.O. FREIRE, ET AL, J. PHYS. CHEM. A 110 (17) (2006) 5897.          '/
 data  allref( 58,4)/ ' Ce: (PM3): R.O. FREIRE, ET AL, J. ORGANOMETALLIC CHEMISTRY 691 (2006) 2584.    '/
 data  allref( 59,4)/ ' Pr: (PM3): R.O. FREIRE, ET AL, J. ORGANOMETALLIC CHEMISTRY 690 (2005) 4099.    '/
 data  allref( 60,4)/ ' Nd: (PM3): C.C. BASTOS, ET AL, J. PHOTOCHEM. PHOTOBIOL. A 177 (2006) 225.      '/
 data  allref( 61,4)/ ' Pm: (PM3): R.O. FREIRE, ET AL, J. CHEM. THEORY COMPUT. 2 (2006) 64.            '/
 data  allref( 62,4)/ ' Sm: (PM3): R.O. FREIRE, ET AL, J. CHEM. THEORY COMPUT. 2 (2006) 64.            '/
 data  allref( 63,4)/ ' Eu: (PM3): R.O. FREIRE, G.B. ROCHA, A.M. SIMAS, INORG. CHEM.44 (2005) 3299.    '/
 data  allref( 64,4)/ ' Gd: (PM3): R.O. FREIRE, G.B. ROCHA, A.M. SIMAS, INORG. CHEM.44 (2005) 3299.    '/
 data  allref( 65,4)/ ' Tb: (PM3): R.O. FREIRE, G.B. ROCHA, A.M. SIMAS, INORG. CHEM.44 (2005) 3299.    '/
 data  allref( 66,4)/ ' Dy: (PM3): N.B. DA COSTA JR, ET AL, INORG. CHEM. COMM. 8 (2005) 831.           '/
 data  allref( 67,4)/ ' Ho: (PM3): N.B. DA COSTA JR, ET AL, POLYHEDRON 24 (2005) 3046.                 '/
 data  allref( 68,4)/ ' Er: (PM3): R.O. FREIRE, ET AL, J. ORGANOMETALLIC CHEMISTRY 691 (2006) 2584.    '/
 data  allref( 69,4)/ ' Tm: (PM3): R.O. FREIRE, ET AL, CHEM. PHYS. LETT., 425 (2006) 138.              '/
 data  allref( 70,4)/ ' Yb: (PM3): R.O. FREIRE, G.B. ROCHA, A.M. SIMAS, J. COMP. CHEM., 26 (2005) 1524.'/
 data  allref( 71,4)/ ' Lu: (PM3): R.O. FREIRE, ET AL, J. PHYS. CHEM. A 110 (17) (2006) 5897.          '/
 data  allref( 80,4)/ ' Hg: (PM3): J. J. P. STEWART, J. COMP. CHEM. 12, 320-341 (1991).                '/
 data  allref( 81,4)/ ' Tl: (PM3): J. J. P. STEWART, J. COMP. CHEM. 12, 320-341 (1991).                '/
 data  allref( 82,4)/ ' Pb: (PM3): J. J. P. STEWART, J. COMP. CHEM. 12, 320-341 (1991).                '/
 data  allref( 83,4)/ ' Bi: (PM3): J. J. P. STEWART, J. COMP. CHEM. 12, 320-341 (1991).                '/
 data  allref(102,4)/ ' Cb: (PM3): Capped Bond  (Hydrogen-like, takes on a  zero charge.)              '/
 data  allref(104,4)/ '  +: (PM3): Sparkle with charge of +1                                           '/
 data  allref(106,4)/ '  -: (PM3): Sparkle with charge of -1                                           '/
!
! PM8
!
 data  allref(1:107,9)/ 107*" PM8 from PARAM (Development version)                                     "/
      end module journal_references_C
