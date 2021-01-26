      module journal_references_C 
      character (len=100), dimension(107) :: refmn, refpm6, refam, refpm3, refmd, refrm1, refpm7, refrm1_sparkles, refpm8
      character, dimension(107,9) :: allref*100
      equivalence (refmn(1),  allref(1,1)), &
                  (refpm6(1), allref(1,2)), &
                  (refam(1),  allref(1,3)), &
                  (refpm3(1), allref(1,4)), &
                  (refmd(1),  allref(1,5)), &
                  (refrm1(1), allref(1,6)), &
                  (refrm1_sparkles (1), allref(1,8)), &
                  (refpm8(1), allref(1,9))
      
      
 data refrm1(  1)/'  H: (RM1): G.B. ROCHA, R.O. FREIRE, A.M.SIMAS, J.J.P. STEWART, J. COMP. CHEM. 27, 1101 (2006)'/
 data refrm1(  6)/'  C: (RM1): G.B. ROCHA, R.O. FREIRE, A.M.SIMAS, J.J.P. STEWART, J. COMP. CHEM. 27, 1101 (2006)'/
 data refrm1(  7)/'  N: (RM1): G.B. ROCHA, R.O. FREIRE, A.M.SIMAS, J.J.P. STEWART, J. COMP. CHEM. 27, 1101 (2006)'/
 data refrm1(  8)/'  O: (RM1): G.B. ROCHA, R.O. FREIRE, A.M.SIMAS, J.J.P. STEWART, J. COMP. CHEM. 27, 1101 (2006)'/
 data refrm1(  9)/'  F: (RM1): G.B. ROCHA, R.O. FREIRE, A.M.SIMAS, J.J.P. STEWART, J. COMP. CHEM. 27, 1101 (2006)'/
 data refrm1( 15)/'  P: (RM1): G.B. ROCHA, R.O. FREIRE, A.M.SIMAS, J.J.P. STEWART, J. COMP. CHEM. 27, 1101 (2006)'/
 data refrm1( 16)/'  S: (RM1): G.B. ROCHA, R.O. FREIRE, A.M.SIMAS, J.J.P. STEWART, J. COMP. CHEM. 27, 1101 (2006)'/
 data refrm1( 17)/' Cl: (RM1): G.B. ROCHA, R.O. FREIRE, A.M.SIMAS, J.J.P. STEWART, J. COMP. CHEM. 27, 1101 (2006)'/ 
 data refrm1( 35)/' Br: (RM1): G.B. ROCHA, R.O. FREIRE, A.M.SIMAS, J.J.P. STEWART, J. COMP. CHEM. 27, 1101 (2006)'/
 data refrm1( 53)/'  I: (RM1): G.B. ROCHA, R.O. FREIRE, A.M.SIMAS, J.J.P. STEWART, J. COMP. CHEM. 27, 1101 (2006)'/
 data  refmd( 11)/ " Na: (MNDO/d): W.THIEL AND A.A.VOITYUK, J. PHYS. CHEM., 100, 616 (1996)                      "/
 data  refmd( 12)/ " Mg: (MNDO/d): W.THIEL AND A.A.VOITYUK, J. PHYS. CHEM., 100, 616 (1996)                      "/
 data  refmd( 13)/ " Al: (MNDO/d): W.THIEL AND A.A.VOITYUK, J. PHYS. CHEM., 100, 616 (1996)                      "/
 data  refmd( 14)/ " Si: (MNDO/d): W.THIEL AND A.A.VOITYUK, J. PHYS. CHEM., 100, 616 (1996)                      "/
 data  refmd( 15)/ "  P: (MNDO/d): W.THIEL AND A.A.VOITYUK, J. PHYS. CHEM., 100, 616 (1996)                      "/
 data  refmd( 16)/ "  S: (MNDO/d): W.THIEL AND A.A.VOITYUK, J. PHYS. CHEM., 100, 616 (1996)                      "/
 data  refmd( 17)/ " Cl: (MNDO/d): W.THIEL AND A.A.VOITYUK, J. PHYS. CHEM., 100, 616 (1996)                      "/
 data  refmd( 30)/ " Zn: (MNDO/d): W.THIEL AND A.A.VOITYUK, J. PHYS. CHEM., 100, 616 (1996)                      "/
 data  refmd( 35)/ " Br: (MNDO/d): W.THIEL AND A.A.VOITYUK, J. PHYS. CHEM., 100, 616 (1996)                      "/
 data  refmd( 48)/ " Cd: (MNDO/d): W.THIEL AND A.A.VOITYUK, J. PHYS. CHEM., 100, 616 (1996)                      "/
 data  refmd( 53)/ "  I: (MNDO/d): W.THIEL AND A.A.VOITYUK, J. PHYS. CHEM., 100, 616 (1996)                      "/
 data  refmd( 80)/ " Hg: (MNDO/d): W.THIEL AND A.A.VOITYUK, J. PHYS. CHEM., 100, 616 (1996)                      "/
!
!  References for the RM1 Lanthanides
!
 data  refrm1(57)/ " La: (RM1): Dutra, J.D.L., Filho, M.A.M., et al., PLOS ONE, 2015, DOI: 10.1371/journal.pone.0124372 "/
 data  refrm1(58)/ " Ce: (RM1): Dutra, J.D.L., Filho, M.A.M., et al., PLOS ONE, 2015, DOI: 10.1371/journal.pone.0124372 "/
 data  refrm1(59)/ " Pr: (RM1): Dutra, J.D.L., Filho, M.A.M., et al., PLOS ONE, 2015, DOI: 10.1371/journal.pone.0124372 "/
 data  refrm1(60)/ " Nd: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., RSC Adv., 2015, DOI: 10.1039/C4RA12682C "/
 data  refrm1(61)/ " Pm: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., RSC Adv., 2015, DOI: 10.1039/C4RA12682C "/
 data  refrm1(62)/ " Sm: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., RSC Adv., 2015, DOI: 10.1039/C4RA12682C "/
 data  refrm1(63)/ " Eu: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., JCTC, 2014, DOI: 10.1021/ct400909w "/
 data  refrm1(64)/ " Gd: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., JCTC, 2014, DOI: 10.1021/ct400909w "/
 data  refrm1(65)/ " Tb: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., JCTC, 2014, DOI: 10.1021/ct400909w "/
 data  refrm1(66)/ " Dy: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., PLOS ONE, 2014, DOI: 10.1371/journal.pone.0086376 "/
 data  refrm1(67)/ " Ho: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., PLOS ONE, 2014, DOI: 10.1371/journal.pone.0086376 "/
 data  refrm1(68)/ " Er: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., PLOS ONE, 2014, DOI: 10.1371/journal.pone.0086376 "/
 data  refrm1(69)/ " Tm: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., PLOS ONE, 2016, DOI: 10.1371/journal.pone.0154500 "/
 data  refrm1(70)/ " Yb: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., PLOS ONE, 2016, DOI: 10.1371/journal.pone.0154500 "/
 data  refrm1(71)/ " Lu: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., PLOS ONE, 2016, DOI: 10.1371/journal.pone.0154500 "/
!
!  References for the RM1 Lanthanides as sparkles
!
 data  refrm1_sparkles(57)/ " La: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., RSC Adv., 2013, DOI: 10.1039/C3RA41406J "/
 data  refrm1_sparkles(58)/ " Ce: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., RSC Adv., 2013, DOI: 10.1039/C3RA41406J "/
 data  refrm1_sparkles(59)/ " Pr: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., RSC Adv., 2013, DOI: 10.1039/C3RA41406J "/
 data  refrm1_sparkles(60)/ " Nd: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., RSC Adv., 2013, DOI: 10.1039/C3RA41406J "/
 data  refrm1_sparkles(61)/ " Pm: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., RSC Adv., 2013, DOI: 10.1039/C3RA41406J "/
 data  refrm1_sparkles(62)/ " Sm: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., RSC Adv., 2013, DOI: 10.1039/C3RA41406J "/
 data  refrm1_sparkles(63)/ " Eu: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., RSC Adv., 2013, DOI: 10.1039/C3RA41406J "/
 data  refrm1_sparkles(64)/ " Gd: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., RSC Adv., 2013, DOI: 10.1039/C3RA41406J "/
 data  refrm1_sparkles(65)/ " Tb: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., RSC Adv., 2013, DOI: 10.1039/C3RA41406J "/
 data  refrm1_sparkles(66)/ " Dy: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., RSC Adv., 2013, DOI: 10.1039/C3RA41406J "/
 data  refrm1_sparkles(67)/ " Ho: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., RSC Adv., 2013, DOI: 10.1039/C3RA41406J "/
 data  refrm1_sparkles(68)/ " Er: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., RSC Adv., 2013, DOI: 10.1039/C3RA41406J "/
 data  refrm1_sparkles(69)/ " Tm: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., RSC Adv., 2013, DOI: 10.1039/C3RA41406J "/
 data  refrm1_sparkles(70)/ " Yb: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., RSC Adv., 2013, DOI: 10.1039/C3RA41406J "/
 data  refrm1_sparkles(71)/ " Lu: (RM1): Filho, M.A.M., Dutra, J.D.L., et al., RSC Adv., 2013, DOI: 10.1039/C3RA41406J "/
!
!  References for the PM6 Lanthanides
!
 data  refpm6(57)/ " La: (PM6): Freire R.O. and Simas, A.M.  J. Chem. Theory and Comput., 2010, 6, 2019-2023     "/
 data  refpm6(58)/ " Ce: (PM6): Freire R.O. and Simas, A.M.  J. Chem. Theory and Comput., 2010, 6, 2019-2023     "/
 data  refpm6(59)/ " Pr: (PM6): Freire R.O. and Simas, A.M.  J. Chem. Theory and Comput., 2010, 6, 2019-2023     "/
 data  refpm6(60)/ " Nd: (PM6): Freire R.O. and Simas, A.M.  J. Chem. Theory and Comput., 2010, 6, 2019-2023     "/
 data  refpm6(61)/ " Pm: (PM6): Freire R.O. and Simas, A.M.  J. Chem. Theory and Comput., 2010, 6, 2019-2023     "/
 data  refpm6(62)/ " Sm: (PM6): Freire R.O. and Simas, A.M.  J. Chem. Theory and Comput., 2010, 6, 2019-2023     "/
 data  refpm6(63)/ " Eu: (PM6): Freire R.O. and Simas, A.M.  J. Chem. Theory and Comput., 2010, 6, 2019-2023     "/
 data  refpm6(64)/ " Gd: (PM6): Freire R.O. and Simas, A.M.  J. Chem. Theory and Comput., 2010, 6, 2019-2023     "/
 data  refpm6(65)/ " Tb: (PM6): Freire R.O. and Simas, A.M.  J. Chem. Theory and Comput., 2010, 6, 2019-2023     "/
 data  refpm6(66)/ " Dy: (PM6): Freire R.O. and Simas, A.M.  J. Chem. Theory and Comput., 2010, 6, 2019-2023     "/
 data  refpm6(67)/ " Ho: (PM6): Freire R.O. and Simas, A.M.  J. Chem. Theory and Comput., 2010, 6, 2019-2023     "/
 data  refpm6(68)/ " Er: (PM6): Freire R.O. and Simas, A.M.  J. Chem. Theory and Comput., 2010, 6, 2019-2023     "/
 data  refpm6(69)/ " Tm: (PM6): Freire R.O. and Simas, A.M.  J. Chem. Theory and Comput., 2010, 6, 2019-2023     "/
 data  refpm6(70)/ " Yb: (PM6): Freire R.O. and Simas, A.M.  J. Chem. Theory and Comput., 2010, 6, 2019-2023     "/
 data  refpm6(71)/ " Lu: (PM6): Freire R.O. and Simas, A.M.  J. Chem. Theory and Comput., 2010, 6, 2019-2023     "/
!
!
!  References for the PM7 Lanthanides
!
!
 data  refpm7(57)/ " La: (PM7): Dutra, J.D.L., Filho, M.A.M., et al J. Chem. Theory and Comput., 2013, 9, 3333-3341 "/
 data  refpm7(58)/ " Ce: (PM7): Dutra, J.D.L., Filho, M.A.M., et al J. Chem. Theory and Comput., 2013, 9, 3333-3341 "/
 data  refpm7(59)/ " Pr: (PM7): Dutra, J.D.L., Filho, M.A.M., et al J. Chem. Theory and Comput., 2013, 9, 3333-3341 "/
 data  refpm7(60)/ " Nd: (PM7): Dutra, J.D.L., Filho, M.A.M., et al J. Chem. Theory and Comput., 2013, 9, 3333-3341 "/
 data  refpm7(61)/ " Pm: (PM7): Dutra, J.D.L., Filho, M.A.M., et al J. Chem. Theory and Comput., 2013, 9, 3333-3341 "/
 data  refpm7(62)/ " Sm: (PM7): Dutra, J.D.L., Filho, M.A.M., et al J. Chem. Theory and Comput., 2013, 9, 3333-3341 "/
 data  refpm7(63)/ " Eu: (PM7): Dutra, J.D.L., Filho, M.A.M., et al J. Chem. Theory and Comput., 2013, 9, 3333-3341 "/
 data  refpm7(64)/ " Gd: (PM7): Dutra, J.D.L., Filho, M.A.M., et al J. Chem. Theory and Comput., 2013, 9, 3333-3341 "/
 data  refpm7(65)/ " Tb: (PM7): Dutra, J.D.L., Filho, M.A.M., et al J. Chem. Theory and Comput., 2013, 9, 3333-3341 "/
 data  refpm7(66)/ " Dy: (PM7): Dutra, J.D.L., Filho, M.A.M., et al J. Chem. Theory and Comput., 2013, 9, 3333-3341 "/
 data  refpm7(67)/ " Ho: (PM7): Dutra, J.D.L., Filho, M.A.M., et al J. Chem. Theory and Comput., 2013, 9, 3333-3341 "/
 data  refpm7(68)/ " Er: (PM7): Dutra, J.D.L., Filho, M.A.M., et al J. Chem. Theory and Comput., 2013, 9, 3333-3341 "/
 data  refpm7(69)/ " Tm: (PM7): Dutra, J.D.L., Filho, M.A.M., et al J. Chem. Theory and Comput., 2013, 9, 3333-3341 "/
 data  refpm7(70)/ " Yb: (PM7): Dutra, J.D.L., Filho, M.A.M., et al J. Chem. Theory and Comput., 2013, 9, 3333-3341 "/
 data  refpm7(71)/ " Lu: (PM7): Dutra, J.D.L., Filho, M.A.M., et al J. Chem. Theory and Comput., 2013, 9, 3333-3341 "/

      data  refmn(  1)  /  &
       '  H: (MNDO): M.J.S. DEWAR, W. THIEL, J. AM. CHEM. SOC., 99, 4899, (1977)&
      &       '/
      data  refmn(  2)  /  &
       ' He: (MNDO):                                                            &
      &       '/
      data  refmn(  3)  /  &
       ' Li: (MNDO): TAKEN FROM MNDOC BY W.THIEL (QCPE NO.438, V. 2, P.63, (1982&
      & ).)   '/
      data  refmn(  4)  /  &
       ' Be: (MNDO): M.J.S. DEWAR, H.S. RZEPA, J. AM. CHEM. SOC., 100, 777, (197&
      & 8)    '/
      data  refmn(  5)  /  &
       '  B: (MNDO): M.J.S. DEWAR, M.L. MCKEE, J. AM. CHEM. SOC., 99, 5231, (197&
      & 7)     '/
      data  refmn(  6)  /  &
       '  C: (MNDO): M.J.S. DEWAR, W. THIEL, J. AM. CHEM. SOC., 99, 4899, (1977)&
      &        '/
      data  refmn(  7)  /  &
       '  N: (MNDO): M.J.S. DEWAR, W. THIEL, J. AM. CHEM. SOC., 99, 4899, (1977)&
      &        '/
      data  refmn(  8)  /  &
       '  O: (MNDO): M.J.S. DEWAR, W. THIEL, J. AM. CHEM. SOC., 99, 4899, (1977)&
      &        '/
      data  refmn(  9)  /  &
       '  F: (MNDO): M.J.S. DEWAR, H.S. RZEPA, J. AM. CHEM. SOC., 100, 777, (197&
      & 8)     '/
      data  refmn( 10)  /  &
       ' Ne: (MNDO):                                                            &
      &        '/
      data  refmn( 11)  /  &
       ' Na: (MNDO): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)         &
      &        '/
      data  refmn( 12)  /  &
       ' Mg: (MNDO): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)         &
      &        '/
      data  refmn( 13)  /  &
       ' Al: (MNDO): L.P. DAVIS, ET.AL.  J. COMP. CHEM., 2, 433, (1981) SEE MANU&
      & AL.    '/
      data  refmn( 14)  /  &
       ' Si: (MNDO): M.J.S.DEWAR, ET. AL. ORGANOMETALLICS  5, 375 (1986)        &
      &        '/
      data  refmn( 15)  /  &
       '  P: (MNDO): M.J.S.DEWAR, M.L.MCKEE, H.S. RZEPA,J. AM. CHEM. SOC., 100 3&
      &607 1978'/
      data  refmn( 16)  /  &
       '  S: (MNDO): M.J.S.DEWAR, C.H. REYNOLDS, J. COMP. CHEM. 7, 140-143 (1986&
      & )      '/
      data  refmn( 17)  /  &
       ' Cl: (MNDO): M.J.S.DEWAR, H.S.RZEPA, J. COMP. CHEM., 4, 158, (1983)     &
      &        '/
      data  refmn( 18)  /  &
       ' Ar: (MNDO):                                                            &
      &        '/
      data  refmn( 19)  /  &
       '  K: (MNDO): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)         &
      &        '/
      data  refmn( 20)  /  &
       ' Ca: (MNDO): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)         &
      &        '/
      data  refmn( 30)  /  &
       ' Zn: (MNDO): M.J.S. DEWAR, K.M. MERZ, ORGANOMETALLICS, 5, 1494-1496 (198&
      & 6)     '/
      data  refmn( 31)  /  &
       ' Ga: (MNDO): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)         &
      &        '/
      data  refmn( 32)  /  &
       ' Ge: (MNDO): M.J.S.DEWAR, G.L.GRADY, E.F.HEALY,ORGANOMETALLICS 6 186-189&
      & (1987)'/
      data  refmn( 33)  /  &
       ' As: (MNDO): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)         &
      &        '/
      data  refmn( 34)  /  &
       ' Se: (MNDO): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)         &
      &        '/
      data  refmn( 35)  /  &
       ' Br: (MNDO): M.J.S.DEWAR, E.F. HEALY, J. COMP. CHEM., 4, 542, (1983)    &
      &        '/
      data  refmn( 36)  /  &
       ' Kr: (MNDO):                                                            &
      &        '/
      data  refmn( 37)  /  &
       ' Rb: (MNDO): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)         &
      &        '/
      data  refmn( 38)  /  &
       ' Sr: (MNDO): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)         &
      &        '/
      data  refmn( 49)  /  &
       ' In: (MNDO): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)         &
      &        '/
      data  refmn( 50)  /  &
       ' Sn: (MNDO): M.J.S.DEWAR,G.L.GRADY,J.J.P.STEWART, J.AM.CHEM.SOC.,106 677&
      &1 (1984)'/
      data  refmn( 51)  /  &
       ' Sb: (MNDO): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)         &
      &        '/
      data  refmn( 52)  /  &
       ' Te: (MNDO): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)         &
      &        '/
      data  refmn( 53)  /  &
       '  I: (MNDO): M.J.S.DEWAR, E.F. HEALY, J.J.P. STEWART, J.COMP.CHEM., 5,35&
      &8,(1984)'/
      data  refmn( 54)  /  &
       ' Xe: (MNDO):                                                            &
      &        '/
      data  refmn( 55)  /  &
       ' Cs: (MNDO): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)         &
      &        '/
      data  refmn( 56)  /  &
       ' Ba: (MNDO): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)         &
      &        '/
      data  refmn( 80)  /  &
       ' Hg: (MNDO): M.J.S.DEWAR,  ET. AL. ORGANOMETALLICS 4, 1964, (1985) SEE M&
      & ANUAL  '/
      data  refmn( 81)  /  &
       ' Tl: (MNDO): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)         &
      &        '/
      data  refmn( 82)  /  &
       ' Pb: (MNDO): M.J.S.DEWAR, ET.AL ORGANOMETALLICS 4, 1973-1980 (1985)     &
      &        '/
      data  refmn( 83)  /  &
       ' Bi: (MNDO): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)         &
      &        '/
      data  refmn(102)  /  &
       ' Cb: (MNDO): Capped Bond  (Hydrogen-like, takes on a  zero charge.)     &
      &        '/
      data  refmn(104)  /  &
       '  +: (MNDO): Sparkle with charge of +1                                  &
      &        '/
      data  refmn(106)  /  &
       '  -: (MNDO): Sparkle with charge of -1                                  &
      &        '/
!
!                                          AM1
!
      data  refam(  1)  /  &
       '  H: (AM1): M.J.S. DEWAR ET AL, J. AM. CHEM. SOC. 107 3902-3909 (1985)  &
      &        '/
      data  refam(  2)  /  &
       ' He: (AM1):                                                             &
      &        '/
      data  refam(  3)  /  &
       ' Li: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)          &
      &        '/
      data  refam(  4)  /  &
       ' Be: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)          &
      &        '/
      data  refam(  6)  /  &
       '  C: (AM1): M.J.S. DEWAR ET AL, J. AM. CHEM. SOC. 107 3902-3909 (1985)  &
      &        '/
      data  refam(  7)  /  &
       '  N: (AM1): M.J.S. DEWAR ET AL, J. AM. CHEM. SOC. 107 3902-3909 (1985)  &
      &        '/
      data  refam(  8)  /  &
       '  O: (AM1): M.J.S. DEWAR ET AL, J. AM. CHEM. SOC. 107 3902-3909 (1985)  &
      &        '/
      data  refam(  9)  /  &
       '  F: (AM1): M.J.S. DEWAR AND E. G. ZOEBISCH, THEOCHEM, 180, 1 (1988).   &
      &        '/
      data  refam( 10)  /  &
       ' Ne: (AM1):                                                             &
      &        '/
      data  refam( 11)  /  &
       ' Na: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)          &
      &        '/
      data  refam( 12)  /  &
       ' Mg: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)          &
      &        '/
      data  refam( 13)  /  &
       ' Al: (AM1): M. J. S. Dewar, A. J. Holder, Organometallics, 9, 508-511 (1&
      & 990).  '/
      data  refam( 14)  /  &
       ' Si: (AM1): M.J.S.DEWAR, C. JIE, ORGANOMETALLICS, 6, 1486-1490 (1987).  &
      &        '/
      data  refam( 15)  /  &
       '  P: (AM1): M.J.S.DEWAR, JIE, C, THEOCHEM, 187, 1 (1989)                &
      &        '/
      data  refam( 16)  /  &
       '  S: (AM1): M.J.S. DEWAR, Y-C YUAN, INORGANIC CHEMISTRY, 29, 3881:3890, &
      & (1990) '/
      data  refam( 17)  /  &
       ' Cl: (AM1): M.J.S. DEWAR AND E. G. ZOEBISCH, THEOCHEM, 180, 1 (1988).   &
      &        '/
      data  refam( 18)  /  &
       ' Ar: (AM1):                                                             &
      &        '/
      data  refam( 19)  /  &
       '  K: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)          &
      &        '/
      data  refam( 20)  /  &
       ' Ca: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)          &
      &        '/
      data  refam( 30)  /  &
       ' Zn: (AM1): M.J.S. DEWAR, K.M. MERZ, ORGANOMETALLICS, 7, 522-524 (1988).&
      &        '/
      data  refam( 31)  /  &
       ' Ga: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)          &
      &        '/
      data  refam( 32)  /  &
       ' Ge: (AM1): M.J.S.Dewar and C.Jie, Organometallics, 8, 1544, (1989)     &
      &        '/
      data  refam( 33)  /  &
       ' As: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)          &
      &        '/
      data  refam( 34)  /  &
       ' Se: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)          &
      &        '/
      data  refam( 35)  /  &
       ' Br: (AM1): M.J.S. DEWAR AND E. G. ZOEBISCH, THEOCHEM, 180, 1 (1988).   &
      &        '/
      data  refam( 36)  /  &
       ' Kr: (AM1):                                                             &
      &        '/
      data  refam( 37)  /  &
       ' Rb: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)          &
      &        '/
      data  refam( 38)  /  &
       ' Sr: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)          &
      &        '/
      data  refam( 42)  /  &
       ' Mo: (AM1): A.A. Voityuk, N. Roesch J. Phys. Chem. A 104, 4089 (2000).  &
      &        '/
      data  refam( 49)  /  &
       ' In: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)          &
      &        '/
      data  refam( 50)  /  &
       ' Sn: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)          &
      &        '/
      data  refam( 51)  /  &
       ' Sb: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)          &
      &        '/
      data  refam( 52)  /  &
       ' Te: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)          &
      &        '/
      data  refam( 53)  /  &
       '  I: (AM1): M.J.S. DEWAR AND E. G. ZOEBISCH, THEOCHEM, 180, 1 (1988).   &
      &        '/
      data  refam( 54)  /  &
       ' Xe: (AM1):                                                             &
      &        '/
      data  refam( 55)  /  &
       ' Cs: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)          &
      &        '/
      data  refam( 56)  /  &
       ' Ba: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)          &
      &        '/
      data  refam( 57) / &
     ' La: (AM1): R.O. FREIRE, ET AL, J. PHYS. CHEM. A 110 (17) (2006) 5897. '/
      data  refam( 58) / & 
     ' Ce: (AM1): R.O. FREIRE, ET AL, J. ORGANOMETALLIC CHEMISTRY 691 (2006) 2&
      &584.          '/   
      data  refam( 59) / & 
     ' Pr: (AM1): R.O. FREIRE, ET AL, J. ORGANOMETALLIC CHEMISTRY 690 (2005) 4&
      &099.    '/  
      data  refam( 60) / & 
     ' Nd: (AM1): C.C. BASTOS, ET AL, J. PHOTOCHEM. PHOTOBIOL. A 177 (2006) 22&
	  &5.      '/ 
      data  refam( 61) / & 
     ' Pm: (AM1): R.O. FREIRE, ET AL, J. CHEM. THEORY COMPUT. 2 (2006) 64.  '/
      data  refam( 62) / &
     ' Sm: (AM1): R.O. FREIRE, ET AL, J. CHEM. THEORY COMPUT. 2 (2006) 64.  '/  
      data  refam( 63) / &
     ' Eu: (AM1): R.O. FREIRE, G.B. ROCHA, A.M. SIMAS, INORG. CHEM.44 (2005) 3&
      &299.    '/ 
      data  refam( 64) / &
     ' Gd: (AM1): R.O. FREIRE, G.B. ROCHA, A.M. SIMAS, INORG. CHEM.44 (2005) 3&
      &299.    '/ 
      data  refam( 65) / &
     ' Tb: (AM1): R.O. FREIRE, G.B. ROCHA, A.M. SIMAS, INORG. CHEM.44 (2005) 3&
      &299.    '/ 
      data  refam( 66) / &
     ' Dy: (AM1): N.B. DA COSTA JR, ET AL, INORG. CHEM. COMM. 8 (2005) 831.   &
      &        '/ 
      data  refam( 67) / &
     ' Ho: (AM1): N.B. DA COSTA JR, ET AL, POLYHEDRON 24 (2005) 3046. '/  
      data  refam( 68) / &
     ' Er: (AM1): R.O. FREIRE, ET AL, J. ORGANOMETALLIC CHEMISTRY 691 (2006) 2&
	  &584.          '/
      data  refam( 69) / &
     ' Tm: (AM1): R.O. FREIRE, ET AL, CHEM. PHYS. LETT., 411 (2005) 61.     '/ 
      data  refam( 70) / &
     ' Yb: (AM1): R.O. FREIRE, G.B. ROCHA, A.M. SIMAS, J. COMP. CHEM., 26 (20&
      &05) 1524.      '/        
      data  refam( 71) / & 
     ' Lu: (AM1): R.O. FREIRE, ET AL, J. PHYS. CHEM. A 110 (17) (2006) 5897. '/
      data  refam( 80)  /  &
       ' Hg: (AM1): M.J.S.Dewar and C.Jie, Organometallics 8, 1547, (1989)      &
      &        '/
      data  refam( 81)  /  &
       ' Tl: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)          &
      &        '/
      data  refam( 82)  /  &
       ' Pb: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)          &
      &        '/
      data  refam( 83)  /  &
       ' Bi: (AM1): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)          &
      &        '/
      data  refam( 90)  /  &
       ' Th: (AM1): Rocha, G.B., D.Sc. Thesis, Dept. Quim. U. Pernambuco, Brazil&
      &  (2002)'/
      data  refam(102)  /  &
       ' Cb: (AM1): Capped Bond  (Hydrogen-like, takes on a  zero charge.)      &
      &        '/
      data  refam(104)  /  &
       '  +: (AM1): Sparkle with charge of +1                                   &
      &        '/
      data  refam(106)  /  &
       '  -: (AM1): Sparkle with charge of -1                                   &
      &        '/
!
!                                          PM3
!
      data  refpm3(  1)  /  &
       '  H: (PM3): J. J. P. STEWART, J. COMP. CHEM. 10, 209 (1989).            &
      &        '/
      data  refpm3(  2)  /  &
       ' He: (PM3):                                                             &
      &        '/
      data  refpm3(  3)  /  &
       ' Li: (PM3): E. ANDERS, R. KOCH, P. FREUNSCHT, J. COMP. CHEM 14 1301-1312&
      &  1993  '/
      data  refpm3(  4)  /  &
       ' Be: (PM3): J. J. P. STEWART, J. COMP. CHEM. 10, 209 (1989).            &
      &        '/
      data  refpm3(  5)  /  &
       '  B: (PM3): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)          &
      &        '/
      data  refpm3(  6)  /  &
       '  C: (PM3): J. J. P. STEWART, J. COMP. CHEM. 10, 209 (1989).            &
      &        '/
      data  refpm3(  7)  /  &
       '  N: (PM3): J. J. P. STEWART, J. COMP. CHEM. 10, 209 (1989).            &
      &        '/
      data  refpm3(  8)  /  &
       '  O: (PM3): J. J. P. STEWART, J. COMP. CHEM. 10, 209 (1989).            &
      &        '/
      data  refpm3(  9)  /  &
       '  F: (PM3): J. J. P. STEWART, J. COMP. CHEM. 10, 209 (1989).            &
      &        '/
      data  refpm3( 10)  /  &
       ' Ne: (PM3):                                                             &
      &        '/
      data  refpm3( 11)  /  &
       ' Na: (PM3): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)          &
      &        '/
      data  refpm3( 12)  /  &
       ' Mg: (PM3): J. J. P. STEWART, J. COMP. CHEM. 12, 320-341 (1991).        &
      &        '/
      data  refpm3( 13)  /  &
       ' Al: (PM3): J. J. P. STEWART, J. COMP. CHEM. 10, 209 (1989).            &
      &        '/
      data  refpm3( 14)  /  &
       ' Si: (PM3): J. J. P. STEWART, J. COMP. CHEM. 10, 209 (1989).            &
      &        '/
      data  refpm3( 15)  /  &
       '  P: (PM3): J. J. P. STEWART, J. COMP. CHEM. 10, 209 (1989).            &
      &        '/
      data  refpm3( 16)  /  &
       '  S: (PM3): J. J. P. STEWART, J. COMP. CHEM. 10, 209 (1989).            &
      &        '/
      data  refpm3( 17)  /  &
       ' Cl: (PM3): J. J. P. STEWART, J. COMP. CHEM. 10, 209 (1989).            &
      &        '/
      data  refpm3( 18)  /  &
       ' Ar: (PM3):                                                             &
      &        '/
      data  refpm3( 19)  /  &
       '  K: (PM3): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)          &
      &        '/
      data  refpm3( 20)  /  &
       ' Ca: (PM3): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)          &
      &        '/
      data  refpm3( 30)  /  &
       ' Zn: (PM3): J. J. P. STEWART, J. COMP. CHEM. 12, 320-341 (1991).        &
      &        '/
      data  refpm3( 31)  /  &
       ' Ga: (PM3): J. J. P. STEWART, J. COMP. CHEM. 12, 320-341 (1991).        &
      &        '/
      data  refpm3( 32)  /  &
       ' Ge: (PM3): J. J. P. STEWART, J. COMP. CHEM. 12, 320-341 (1991).        &
      &        '/
      data  refpm3( 33)  /  &
       ' As: (PM3): J. J. P. STEWART, J. COMP. CHEM. 12, 320-341 (1991).        &
      &        '/
      data  refpm3( 34)  /  &
       ' Se: (PM3): J. J. P. STEWART, J. COMP. CHEM. 12, 320-341 (1991).        &
      &        '/
      data  refpm3( 35)  /  &
       ' Br: (PM3): J. J. P. STEWART, J. COMP. CHEM. 10, 209 (1989).            &
      &        '/
      data  refpm3( 36)  /  &
       ' Kr: (PM3):                                                             &
      &        '/
      data  refpm3( 37)  /  &
       ' Rb: (PM3): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)          &
      &        '/
      data  refpm3( 38)  /  &
       ' Sr: (PM3): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)          &
      &        '/
      data  refpm3( 48)  /  &
       ' Cd: (PM3): J. J. P. STEWART, J. COMP. CHEM. 12, 320-341 (1991).        &
      &        '/
      data  refpm3( 49)  /  &
       ' In: (PM3): J. J. P. STEWART, J. COMP. CHEM. 12, 320-341 (1991).        &
      &        '/
      data  refpm3( 50)  /  &
       ' Sn: (PM3): J. J. P. STEWART, J. COMP. CHEM. 12, 320-341 (1991).        &
      &        '/
      data  refpm3( 51)  /  &
       ' Sb: (PM3): J. J. P. STEWART, J. COMP. CHEM. 12, 320-341 (1991).        &
      &        '/
      data  refpm3( 52)  /  &
       ' Te: (PM3): J. J. P. STEWART, J. COMP. CHEM. 12, 320-341 (1991).        &
      &        '/
      data  refpm3( 53)  /  &
       '  I: (PM3): J. J. P. STEWART, J. COMP. CHEM. 10, 209 (1989).            &
      &        '/
      data  refpm3( 54)  /  &
       ' Xe: (PM3):                                                             &
      &        '/
      data  refpm3( 55)  /  &
       ' Cs: (PM3): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)          &
      &        '/
      data  refpm3( 56)  /  &
       ' Ba: (PM3): J. J. P. STEWART, J. MOL. MOD., 10, 155-164 (2004)          &
      &        '/
      data  refpm3( 57) / &
     ' La: (PM3): R.O. FREIRE, ET AL, J. PHYS. CHEM. A 110 (17) (2006) 5897. '/
      data  refpm3( 58) / & 
     ' Ce: (PM3): R.O. FREIRE, ET AL, J. ORGANOMETALLIC CHEMISTRY 691 (2006) 2&
      &584.          '/   
      data  refpm3( 59) / & 
     ' Pr: (PM3): R.O. FREIRE, ET AL, J. ORGANOMETALLIC CHEMISTRY 690 (2005) 4&
      &099.    '/  
      data  refpm3( 60) / & 
     ' Nd: (PM3): C.C. BASTOS, ET AL, J. PHOTOCHEM. PHOTOBIOL. A 177 (2006) 22&
	  &5.      '/ 
      data  refpm3( 61) / & 
     ' Pm: (PM3): R.O. FREIRE, ET AL, J. CHEM. THEORY COMPUT. 2 (2006) 64.  '/
      data  refpm3( 62) / &
     ' Sm: (PM3): R.O. FREIRE, ET AL, J. CHEM. THEORY COMPUT. 2 (2006) 64.  '/  
      data  refpm3( 63) / &
     ' Eu: (PM3): R.O. FREIRE, G.B. ROCHA, A.M. SIMAS, INORG. CHEM.44 (2005) 3&
      &299.    '/ 
      data  refpm3( 64) / &
     ' Gd: (PM3): R.O. FREIRE, G.B. ROCHA, A.M. SIMAS, INORG. CHEM.44 (2005) 3&
      &299.    '/ 
      data  refpm3( 65) / &
     ' Tb: (PM3): R.O. FREIRE, G.B. ROCHA, A.M. SIMAS, INORG. CHEM.44 (2005) 3&
      &299.    '/ 
      data  refpm3( 66) / &
     ' Dy: (PM3): N.B. DA COSTA JR, ET AL, INORG. CHEM. COMM. 8 (2005) 831.   &
      &        '/ 
      data  refpm3( 67) / &
     ' Ho: (PM3): N.B. DA COSTA JR, ET AL, POLYHEDRON 24 (2005) 3046. '/  
      data  refpm3( 68) / &
     ' Er: (PM3): R.O. FREIRE, ET AL, J. ORGANOMETALLIC CHEMISTRY 691 (2006) 2&
	  &584.          '/
      data  refpm3( 69) / &
     ' Tm: (PM3): R.O. FREIRE, ET AL, CHEM. PHYS. LETT., 425 (2006) 138.    '/ 
      data  refpm3( 70) / &
     ' Yb: (PM3): R.O. FREIRE, G.B. ROCHA, A.M. SIMAS, J. COMP. CHEM., 26 (20&
      &05) 1524.      '/        
      data  refpm3( 71) / & 
     ' Lu: (PM3): R.O. FREIRE, ET AL, J. PHYS. CHEM. A 110 (17) (2006) 5897. '/
      data  refpm3( 80)  /  &
       ' Hg: (PM3): J. J. P. STEWART, J. COMP. CHEM. 12, 320-341 (1991).        &
      &        '/
      data  refpm3( 81)  /  &
       ' Tl: (PM3): J. J. P. STEWART, J. COMP. CHEM. 12, 320-341 (1991).        &
      &        '/
      data  refpm3( 82)  /  &
       ' Pb: (PM3): J. J. P. STEWART, J. COMP. CHEM. 12, 320-341 (1991).        &
      &        '/
      data  refpm3( 83)  /  &
       ' Bi: (PM3): J. J. P. STEWART, J. COMP. CHEM. 12, 320-341 (1991).        &
      &        '/
      data  refpm3(102)  /  &
       ' Cb: (PM3): Capped Bond  (Hydrogen-like, takes on a  zero charge.)      &
      &        '/
      data  refpm3(104)  /  &
       '  +: (PM3): Sparkle with charge of +1                                   &
      &        '/
      data  refpm3(106)  /  &
       '  -: (PM3): Sparkle with charge of -1                                   &
      &        '/
      data refpm8 /107*" PM8 from PARAM (Development version)                                                           "/
      end module journal_references_C
