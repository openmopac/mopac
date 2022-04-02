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

  module Parameters_for_INDO_C
    integer      :: isoki(80),nbfai(80)
    double precision :: zetai(80),zetadi(2,80),zetawti(2,80), &
                    zcoreai(80),betaai(3,80), fgi(24,80)

! INDO parameter structure:
! isok      = 1 (full parameter set exists for element) or 0 (otherwise)
! nfba      = # basis functions on atom (corresponds to natorb)
! zcorea    = nuclear charge (corresponds to tore)
!
! zeta      = sp Slater orbital exponent
! zetad(1)  =  d Slater orbital exponent 1 (d orbitals are a linear combination
!              of two Slater functions)
! zetad(2)  =  d Slater orbital exponent 2
! zetawt(1) =  d Slater orbital coefficient 1
! zetawt(2) =  d Slater orbital coefficient 2
!
! betaa(1)  = Resonance integral for s orbitals
! betaa(2)  = Resonance integral for p orbitals
! betaa(3)  = Resonance integral for d orbitals
!
! fg(1)     = s ioniz pot from d(n-2)s2 or from sp
! fg(2)     = s ioniz pot from d(n-1)s1
! fg(3)     = p ioniz pot from d(n-2)s2 or sp
! fg(4)     = p ioniz pot from d(n-1)s1
! fg(5)     = d ioniz pot from d(n-2)s2 or sp
! fg(6)     = d ioniz pot from d(n-1)s1
! fg(7)     = d ioniz pot from dns0
! fg(8)     = Probability of d(n-2)s2 or sp
! fg(9)     = Probability of d(n-1)s1
! fg(10)    = Probability of dns0
!
! fg(11)    = F0(s,s) (=F0sp,F0pp) Slater-Condon one-center-two-electron integral
! fg(12)    = F0(s,d)
! fg(13)    = F0(d,d)
! fg(14)    = G1(s,p) [Higher-order Slater-Condon factors]
! fg(15)    = F2(p,p)
! fg(16)    = G2(s,d)
! fg(17)    = G1(p,d)
! fg(18)    = F2(p,d)
! fg(19)    = G3(p,d)
! fg(20)    = F2(d,d)
! fg(21)    = F4(d,d)
! fg(22)    = R1(sppd)
! fg(23)    = R3(sddd)
! fg(24)    = R2(sdpp)

!
!                    Data for Element   1
!uu
      data        isoki(  1)/    1         /
      data        nbfai(  1)/    1         /
      data      zcoreai(  1)/    1.000000d0/
      data        zetai(  1)/    1.200000d0/
      data   betaai(  1,  1)/  -12.000000d0/
      data      fgi(  1,  1)/  -13.060000d0/
      data      fgi(  8,  1)/    1.000000d0/
      data      fgi( 11,  1)/   12.850000d0/
      data      fgi( 12,  1)/   12.850000d0/
!
!                    Data for Element   2
!
      data        isoki(  2)/    0         /
      data        nbfai(  2)/    1         /
      data      zcoreai(  2)/    2.000000d0/
      data        zetai(  2)/    1.700000d0/
      data   betaai(  1,  2)/ -100.000000d0/
      data      fgi(  8,  2)/    1.000000d0/
!
!                    Data for Element   3
!
      data        isoki(  3)/    1         /
      data        nbfai(  3)/    4         /
      data      zcoreai(  3)/    1.000000d0/
      data        zetai(  3)/    0.650000d0/
      data   betaai(  1,  3)/   -1.000000d0/
      data   betaai(  2,  3)/   -1.000000d0/
      data      fgi(  1,  3)/   -5.410000d0/
      data      fgi(  3,  3)/   -3.610000d0/
      data      fgi(  8,  3)/    1.000000d0/
      data      fgi( 11,  3)/    4.570000d0/
      data      fgi( 12,  3)/    4.570000d0/
      data      fgi( 14,  3)/    2.503730d0/
      data      fgi( 15,  3)/    1.356879d0/
!
!                    Data for Element   4
!
      data        isoki(  4)/    1         /
      data        nbfai(  4)/    4         /
      data      zcoreai(  4)/    2.000000d0/
      data        zetai(  4)/    0.975000d0/
      data   betaai(  1,  4)/  -13.000000d0/
      data   betaai(  2,  4)/  -13.000000d0/
      data      fgi(  1,  4)/   -9.330000d0/
      data      fgi(  3,  4)/   -5.880000d0/
      data      fgi(  8,  4)/    1.000000d0/
      data      fgi( 11,  4)/    6.780000d0/
      data      fgi( 12,  4)/    6.780000d0/
      data      fgi( 14,  4)/    3.828126d0/
      data      fgi( 15,  4)/    2.656354d0/
!
!                    Data for Element   5
!
      data        isoki(  5)/    1         /
      data        nbfai(  5)/    4         /
      data      zcoreai(  5)/    3.000000d0/
      data        zetai(  5)/    1.300000d0/
      data   betaai(  1,  5)/   -8.000000d0/
      data   betaai(  2,  5)/   -8.000000d0/
      data      fgi(  1,  5)/  -14.000000d0/
      data      fgi(  3,  5)/   -8.240000d0/
      data      fgi(  8,  5)/    1.000000d0/
      data      fgi( 11,  5)/    8.680000d0/
      data      fgi( 12,  5)/    8.680000d0/
      data      fgi( 14,  5)/    5.401481d0/
      data      fgi( 15,  5)/    3.480847d0/
!
!                    Data for Element   6
!
      data        isoki(  6)/    1         /
      data        nbfai(  6)/    4         /
      data      zcoreai(  6)/    4.000000d0/
      data        zetai(  6)/    1.625000d0/
      data   betaai(  1,  6)/  -17.000000d0/
      data   betaai(  2,  6)/  -17.000000d0/
      data      fgi(  1,  6)/  -19.420000d0/
      data      fgi(  3,  6)/  -10.700000d0/
      data      fgi(  8,  6)/    1.000000d0/
      data      fgi( 11,  6)/   11.110000d0/
      data      fgi( 12,  6)/   11.110000d0/
      data      fgi( 14,  6)/    6.897842d0/
      data      fgi( 15,  6)/    4.509913d0/
!
!                    Data for Element   7
!
      data        isoki(  7)/    1         /
      data        nbfai(  7)/    4         /
      data      zcoreai(  7)/    5.000000d0/
      data        zetai(  7)/    1.950000d0/
      data   betaai(  1,  7)/  -26.000000d0/
      data   betaai(  2,  7)/  -26.000000d0/
      data      fgi(  1,  7)/  -25.580000d0/
      data      fgi(  3,  7)/  -13.250000d0/
      data      fgi(  8,  7)/    1.000000d0/
      data      fgi( 11,  7)/   12.010000d0/
      data      fgi( 12,  7)/   12.010000d0/
      data      fgi( 14,  7)/    8.958454d0/
      data      fgi( 15,  7)/    6.459559d0/
!
!                    Data for Element   8
!
      data        isoki(  8)/    1         /
      data        nbfai(  8)/    4         /
      data      zcoreai(  8)/    6.000000d0/
      data        zetai(  8)/    2.275000d0/
      data   betaai(  1,  8)/  -34.000000d0/
      data   betaai(  2,  8)/  -34.000000d0/
      data      fgi(  1,  8)/  -32.490000d0/
      data      fgi(  3,  8)/  -15.880000d0/
      data      fgi(  8,  8)/    1.000000d0/
      data      fgi( 11,  8)/   13.000000d0/
      data      fgi( 12,  8)/   13.000000d0/
      data      fgi( 14,  8)/   11.815414d0/
      data      fgi( 15,  8)/    6.902802d0/
!
!                    Data for Element   9
!
      data        isoki(  9)/    1         /
      data        nbfai(  9)/    4         /
      data      zcoreai(  9)/    7.000000d0/
      data        zetai(  9)/    2.600000d0/
      data   betaai(  1,  9)/  -44.000000d0/
      data   betaai(  2,  9)/  -44.000000d0/
      data      fgi(  1,  9)/  -40.140000d0/
      data      fgi(  3,  9)/  -18.610000d0/
      data      fgi(  8,  9)/    1.000000d0/
      data      fgi( 11,  9)/   14.000000d0/
      data      fgi( 12,  9)/   14.000000d0/
      data      fgi( 14,  9)/   14.484415d0/
      data      fgi( 15,  9)/    8.593198d0/
!
!                    Data for Element  10
!
      data        isoki( 10)/    0         /
      data        nbfai( 10)/    4         /
      data      zcoreai( 10)/    8.000000d0/
      data        zetai( 10)/    2.925000d0/
      data   betaai(  1, 10)/ -100.000000d0/
      data   betaai(  2, 10)/ -100.000000d0/
      data      fgi(  8, 10)/    1.000000d0/
      data      fgi( 14, 10)/   12.781917d0/
      data      fgi( 15, 10)/    9.327345d0/
!
!                    Data for Element  11
!
      data        isoki( 11)/    1         /
      data        nbfai( 11)/    4         /
      data      zcoreai( 11)/    1.000000d0/
      data        zetai( 11)/    0.836000d0/
      data   betaai(  1, 11)/   -5.000000d0/
      data   betaai(  2, 11)/   -5.000000d0/
      data      fgi(  1, 11)/   -4.860000d0/
      data      fgi(  3, 11)/   -2.860000d0/
      data      fgi(  8, 11)/    1.000000d0/
      data      fgi( 11, 11)/    3.310000d0/
      data      fgi( 12, 11)/    3.310000d0/
      data      fgi( 13, 11)/    1.670000d0/
      data      fgi( 14, 11)/    1.667583d0/
      data      fgi( 15, 11)/    0.743903d0/
!
!                    Data for Element  12
!
      data        isoki( 12)/    1         /
      data        nbfai( 12)/    4         /
      data      zcoreai( 12)/    2.000000d0/
      data        zetai( 12)/    1.103000d0/
      data   betaai(  1, 12)/   -6.000000d0/
      data   betaai(  2, 12)/   -6.000000d0/
      data      fgi(  1, 12)/   -8.110000d0/
      data      fgi(  3, 12)/   -4.550000d0/
      data      fgi(  8, 12)/    1.000000d0/
      data      fgi( 11, 12)/    4.790000d0/
      data      fgi( 12, 12)/    4.790000d0/
      data      fgi( 13, 12)/    2.430000d0/
      data      fgi( 14, 12)/    2.476454d0/
      data      fgi( 15, 12)/    3.273174d0/
!
!                    Data for Element  13
!
      data        isoki( 13)/    1         /
      data        nbfai( 13)/    4         /
      data      zcoreai( 13)/    3.000000d0/
      data        zetai( 13)/    1.370000d0/
      data   betaai(  1, 13)/   -7.000000d0/
      data   betaai(  2, 13)/   -7.000000d0/
      data      fgi(  1, 13)/  -11.420000d0/
      data      fgi(  3, 13)/   -6.290000d0/
      data      fgi(  8, 13)/    1.000000d0/
      data      fgi( 11, 13)/    6.210000d0/
      data      fgi( 12, 13)/    6.210000d0/
      data      fgi( 13, 13)/    3.420000d0/
      data      fgi( 14, 13)/    3.359095d0/
      data      fgi( 15, 13)/    1.602491d0/
!
!                    Data for Element  14
!
      data        isoki( 14)/    1         /
      data        nbfai( 14)/    4         /
      data      zcoreai( 14)/    4.000000d0/
      data        zetai( 14)/    1.520000d0/
      data   betaai(  1, 14)/   -9.000000d0/
      data   betaai(  2, 14)/   -9.000000d0/
      data      fgi(  1, 14)/  -14.790000d0/
      data      fgi(  3, 14)/   -8.100000d0/
      data      fgi(  8, 14)/    1.000000d0/
      data      fgi( 11, 14)/    7.570000d0/
      data      fgi( 12, 14)/    7.570000d0/
      data      fgi( 13, 14)/    4.630000d0/
      data      fgi( 14, 14)/    4.812310d0/
      data      fgi( 15, 14)/    2.262706d0/
!
!                    Data for Element  15
!
      data        isoki( 15)/    1         /
      data        nbfai( 15)/    4         /
      data      zcoreai( 15)/    5.000000d0/
      data        zetai( 15)/    1.730000d0/
      data   betaai(  1, 15)/  -15.000000d0/
      data   betaai(  2, 15)/  -15.000000d0/
      data      fgi(  1, 15)/  -18.230000d0/
      data      fgi(  3, 15)/   -9.980000d0/
      data      fgi(  8, 15)/    1.000000d0/
      data      fgi( 11, 15)/    8.860000d0/
      data      fgi( 12, 15)/    8.860000d0/
      data      fgi( 13, 15)/    6.090000d0/
      data      fgi( 14, 15)/    1.047788d0/
      data      fgi( 15, 15)/    2.947716d0/
!
!                    Data for Element  16
!
      data        isoki( 16)/    1         /
      data        nbfai( 16)/    4         /
      data      zcoreai( 16)/    6.000000d0/
      data        zetai( 16)/    1.925000d0/
      data   betaai(  1, 16)/  -15.000000d0/
      data   betaai(  2, 16)/  -15.000000d0/
      data      fgi(  1, 16)/  -21.730000d0/
      data      fgi(  3, 16)/  -11.920000d0/
      data      fgi(  8, 16)/    1.000000d0/
      data      fgi( 11, 16)/   10.090000d0/
      data      fgi( 12, 16)/   10.090000d0/
      data      fgi( 13, 16)/    7.770000d0/
      data      fgi( 14, 16)/    3.075668d0/
      data      fgi( 15, 16)/    4.537809d0/
!
!                    Data for Element  17
!
      data        isoki( 17)/    1         /
      data        nbfai( 17)/    4         /
      data      zcoreai( 17)/    7.000000d0/
      data        zetai( 17)/    2.130000d0/
      data   betaai(  1, 17)/  -11.000000d0/
      data   betaai(  2, 17)/  -11.000000d0/
      data      fgi(  1, 17)/  -25.290000d0/
      data      fgi(  3, 17)/  -13.930000d0/
      data      fgi(  8, 17)/    1.000000d0/
      data      fgi( 11, 17)/   11.250000d0/
      data      fgi( 12, 17)/   11.250000d0/
      data      fgi( 13, 17)/    9.680000d0/
      data      fgi( 14, 17)/    8.802854d0/
      data      fgi( 15, 17)/    6.447161d0/
!
!                    Data for Element  18
!
      data        isoki( 18)/    0         /
      data        nbfai( 18)/    4         /
      data      zcoreai( 18)/    8.000000d0/
      data        zetai( 18)/    2.365000d0/
      data   betaai(  1, 18)/ -100.000000d0/
      data   betaai(  2, 18)/ -100.000000d0/
      data      fgi(  8, 18)/    1.000000d0/
      data      fgi( 14, 18)/    7.762259d0/
      data      fgi( 15, 18)/    5.846134d0/
!
!                    Data for Element  19
!
      data        isoki( 19)/    1         /
      data        nbfai( 19)/    4         /
      data      zcoreai( 19)/    1.000000d0/
      data        zetai( 19)/    1.180000d0/
      data   betaai(  1, 19)/   -1.000000d0/
      data   betaai(  2, 19)/   -1.000000d0/
      data      fgi(  1, 19)/   -4.340000d0/
      data      fgi(  3, 19)/   -2.730000d0/
      data      fgi(  8, 19)/    1.000000d0/
      data      fgi( 11, 19)/    3.180000d0/
      data      fgi( 12, 19)/    3.180000d0/
      data      fgi( 13, 19)/    5.030000d0/
      data      fgi( 14, 19)/    1.110895d0/
      data      fgi( 15, 19)/    0.495935d0/
!
!                    Data for Element  20
!
      data        isoki( 20)/    1         /
      data        nbfai( 20)/    9         /
      data      zcoreai( 20)/    2.000000d0/
      data        zetai( 20)/    1.210000d0/
      data   zetadi(  1, 20)/    1.850000d0/
      data  zetawti(  1, 20)/    1.000000d0/
      data   betaai(  1, 20)/    2.000000d0/
      data   betaai(  2, 20)/    2.000000d0/
      data   betaai(  3, 20)/  -11.400000d0/
      data      fgi(  1, 20)/   -6.030000d0/
      data      fgi(  2, 20)/   -5.130000d0/
      data      fgi(  3, 20)/   -3.960000d0/
      data      fgi(  4, 20)/   -2.990000d0/
      data      fgi(  5, 20)/   -3.440000d0/
      data      fgi(  6, 20)/   -3.440000d0/
      data      fgi(  7, 20)/ -100.000000d0/
      data      fgi(  8, 20)/    0.960800d0/
      data      fgi(  9, 20)/    0.039200d0/
      data      fgi( 11, 20)/    3.250000d0/
      data      fgi( 12, 20)/    4.000000d0/
      data      fgi( 13, 20)/    6.030000d0/
      data      fgi( 14, 20)/    1.562197d0/
      data      fgi( 15, 20)/    0.288262d0/
      data      fgi( 16, 20)/    0.462460d0/
      data      fgi( 17, 20)/    0.730265d0/
      data      fgi( 18, 20)/    0.555448d0/
      data      fgi( 19, 20)/   -0.029508d0/
      data      fgi( 20, 20)/    2.343295d0/
      data      fgi( 21, 20)/    1.177847d0/
      data      fgi( 22, 20)/    2.146212d0/
      data      fgi( 23, 20)/    2.216599d0/
      data      fgi( 24, 20)/    1.575989d0/
!
!                    Data for Element  21
!
      data        isoki( 21)/    1         /
      data        nbfai( 21)/    9         /
      data      zcoreai( 21)/    3.000000d0/
      data        zetai( 21)/    1.230000d0/
      data   zetadi(  1, 21)/    4.222400d0/
      data   zetadi(  2, 21)/    1.746500d0/
      data  zetawti(  1, 21)/    0.359220d0/
      data  zetawti(  2, 21)/    0.766010d0/
      data   betaai(  1, 21)/   -1.000000d0/
      data   betaai(  2, 21)/   -1.000000d0/
      data   betaai(  3, 21)/  -18.000000d0/
      data      fgi(  1, 21)/   -6.720000d0/
      data      fgi(  2, 21)/   -5.830000d0/
      data      fgi(  3, 21)/   -4.200000d0/
      data      fgi(  4, 21)/   -3.430000d0/
      data      fgi(  5, 21)/   -8.160000d0/
      data      fgi(  6, 21)/   -4.850000d0/
      data      fgi(  7, 21)/ -100.000000d0/
      data      fgi(  8, 21)/    0.939900d0/
      data      fgi(  9, 21)/    0.060100d0/
      data      fgi( 11, 21)/    3.890000d0/
      data      fgi( 12, 21)/    4.710000d0/
      data      fgi( 13, 21)/    7.020000d0/
      data      fgi( 14, 21)/    1.500205d0/
      data      fgi( 15, 21)/    0.619919d0/
      data      fgi( 16, 21)/    0.727785d0/
      data      fgi( 17, 21)/    0.700509d0/
      data      fgi( 18, 21)/    1.363822d0/
      data      fgi( 19, 21)/    0.274004d0/
      data      fgi( 20, 21)/    3.657524d0/
      data      fgi( 21, 21)/    1.810164d0/
      data      fgi( 22, 21)/    1.961795d0/
      data      fgi( 23, 21)/    2.115305d0/
      data      fgi( 24, 21)/    1.424637d0/
!
!                    Data for Element  22
!
      data        isoki( 22)/    1         /
      data        nbfai( 22)/    9         /
      data      zcoreai( 22)/    4.000000d0/
      data        zetai( 22)/    1.300000d0/
      data   zetadi(  1, 22)/    4.670000d0/
      data   zetadi(  2, 22)/    1.986000d0/
      data  zetawti(  1, 22)/    0.364610d0/
      data  zetawti(  2, 22)/    0.755610d0/
      data   betaai(  1, 22)/   -1.000000d0/
      data   betaai(  2, 22)/   -1.000000d0/
      data   betaai(  3, 22)/  -19.000000d0/
      data      fgi(  1, 22)/   -7.280000d0/
      data      fgi(  2, 22)/   -6.340000d0/
      data      fgi(  3, 22)/   -4.480000d0/
      data      fgi(  4, 22)/   -3.750000d0/
      data      fgi(  5, 22)/   -9.070000d0/
      data      fgi(  6, 22)/   -5.930000d0/
      data      fgi(  7, 22)/ -100.000000d0/
      data      fgi(  8, 22)/    0.906900d0/
      data      fgi(  9, 22)/    0.093100d0/
      data      fgi( 11, 22)/    4.500000d0/
      data      fgi( 12, 22)/    5.380000d0/
      data      fgi( 13, 22)/    7.980000d0/
      data      fgi( 14, 22)/    1.624189d0/
      data      fgi( 15, 22)/    0.681911d0/
      data      fgi( 16, 22)/    0.768700d0/
      data      fgi( 17, 22)/    0.907562d0/
      data      fgi( 18, 22)/    1.698579d0/
      data      fgi( 19, 22)/    1.277034d0/
      data      fgi( 20, 22)/    5.566875d0/
      data      fgi( 21, 22)/    3.682321d0/
      data      fgi( 22, 22)/    1.571618d0/
      data      fgi( 23, 22)/    1.883618d0/
      data      fgi( 24, 22)/    1.107171d0/
!
!                    Data for Element  23
!
      data        isoki( 23)/    1         /
      data        nbfai( 23)/    9         /
      data      zcoreai( 23)/    5.000000d0/
      data        zetai( 23)/    1.300000d0/
      data   zetadi(  1, 23)/    5.052000d0/
      data   zetadi(  2, 23)/    2.173000d0/
      data  zetawti(  1, 23)/    0.373780d0/
      data  zetawti(  2, 23)/    0.745640d0/
      data   betaai(  1, 23)/   -1.000000d0/
      data   betaai(  2, 23)/   -1.000000d0/
      data   betaai(  3, 23)/  -20.000000d0/
      data      fgi(  1, 23)/   -7.730000d0/
      data      fgi(  2, 23)/   -6.710000d0/
      data      fgi(  3, 23)/   -4.770000d0/
      data      fgi(  4, 23)/   -3.950000d0/
      data      fgi(  5, 23)/   -9.890000d0/
      data      fgi(  6, 23)/   -6.770000d0/
      data      fgi(  7, 23)/ -100.000000d0/
      data      fgi(  8, 23)/    0.839500d0/
      data      fgi(  9, 23)/    0.160500d0/
      data      fgi( 11, 23)/    5.070000d0/
      data      fgi( 12, 23)/    6.010000d0/
      data      fgi( 13, 23)/    8.910000d0/
      data      fgi( 14, 23)/    1.872156d0/
      data      fgi( 15, 23)/    0.743903d0/
      data      fgi( 16, 23)/    0.773659d0/
      data      fgi( 17, 23)/    0.642236d0/
      data      fgi( 18, 23)/    1.388619d0/
      data      fgi( 19, 23)/    0.212012d0/
      data      fgi( 20, 23)/    6.298380d0/
      data      fgi( 21, 23)/    4.389029d0/
      data      fgi( 22, 23)/    1.298119d0/
      data      fgi( 23, 23)/    1.669320d0/
      data      fgi( 24, 23)/    0.894537d0/
!
!                    Data for Element  24
!
      data        isoki( 24)/    1         /
      data        nbfai( 24)/    9         /
      data      zcoreai( 24)/    6.000000d0/
      data        zetai( 24)/    1.320000d0/
      data   zetadi(  1, 24)/    5.138000d0/
      data   zetadi(  2, 24)/    2.077000d0/
      data  zetawti(  1, 24)/    0.407140d0/
      data  zetawti(  2, 24)/    0.732420d0/
      data   betaai(  1, 24)/   -1.000000d0/
      data   betaai(  2, 24)/   -1.000000d0/
      data   betaai(  3, 24)/  -21.000000d0/
      data      fgi(  1, 24)/   -8.070000d0/
      data      fgi(  2, 24)/   -6.970000d0/
      data      fgi(  3, 24)/   -5.040000d0/
      data      fgi(  4, 24)/   -4.060000d0/
      data      fgi(  5, 24)/  -10.660000d0/
      data      fgi(  6, 24)/   -7.430000d0/
      data      fgi(  7, 24)/ -100.000000d0/
      data      fgi(  8, 24)/    0.705200d0/
      data      fgi(  9, 24)/    0.294800d0/
      data      fgi( 11, 24)/    5.600000d0/
      data      fgi( 12, 24)/    6.600000d0/
      data      fgi( 13, 24)/    9.810000d0/
      data      fgi( 14, 24)/    1.785368d0/
      data      fgi( 15, 24)/    0.805895d0/
      data      fgi( 16, 24)/    0.647196d0/
      data      fgi( 17, 24)/    0.691830d0/
      data      fgi( 18, 24)/    1.413416d0/
      data      fgi( 19, 24)/    0.036823d0/
      data      fgi( 20, 24)/    7.872975d0/
      data      fgi( 21, 24)/    4.562606d0/
      data      fgi( 22, 24)/    1.138395d0/
      data      fgi( 23, 24)/    1.544197d0/
      data      fgi( 24, 24)/    0.770846d0/
!
!                    Data for Element  25
!
      data        isoki( 25)/    1         /
      data        nbfai( 25)/    9         /
      data      zcoreai( 25)/    7.000000d0/
      data        zetai( 25)/    1.360000d0/
      data   zetadi(  1, 25)/    5.767000d0/
      data   zetadi(  2, 25)/    2.510000d0/
      data  zetawti(  1, 25)/    0.389840d0/
      data  zetawti(  2, 25)/    0.729650d0/
      data   betaai(  1, 25)/   -1.000000d0/
      data   betaai(  2, 25)/   -1.000000d0/
      data   betaai(  3, 25)/  -22.000000d0/
      data      fgi(  1, 25)/   -8.350000d0/
      data      fgi(  2, 25)/   -7.150000d0/
      data      fgi(  3, 25)/   -5.270000d0/
      data      fgi(  4, 25)/   -4.100000d0/
      data      fgi(  5, 25)/  -11.450000d0/
      data      fgi(  6, 25)/   -7.990000d0/
      data      fgi(  7, 25)/ -100.000000d0/
      data      fgi(  8, 25)/    0.665200d0/
      data      fgi(  9, 25)/    0.334800d0/
      data      fgi( 11, 25)/    6.090000d0/
      data      fgi( 12, 25)/    7.160000d0/
      data      fgi( 13, 25)/   10.680000d0/
      data      fgi( 14, 25)/    2.343295d0/
      data      fgi( 15, 25)/    0.867887d0/
      data      fgi( 16, 25)/    0.757541d0/
      data      fgi( 17, 25)/    0.153740d0/
      data      fgi( 18, 25)/    0.993111d0/
      data      fgi( 19, 25)/    0.616200d0/
      data      fgi( 20, 25)/    8.182935d0/
      data      fgi( 21, 25)/    4.698988d0/
      data      fgi( 22, 25)/    1.054290d0/
      data      fgi( 23, 25)/    1.486488d0/
      data      fgi( 24, 25)/    0.704629d0/
!
!                    Data for Element  26
!
      data        isoki( 26)/    1         /
      data        nbfai( 26)/    9         /
      data      zcoreai( 26)/    8.000000d0/
      data        zetai( 26)/    1.370000d0/
      data   zetadi(  1, 26)/    6.068000d0/
      data   zetadi(  2, 26)/    2.618000d0/
      data  zetawti(  1, 26)/    0.403790d0/
      data  zetawti(  2, 26)/    0.719840d0/
      data   betaai(  1, 26)/   -1.000000d0/
      data   betaai(  2, 26)/   -1.000000d0/
      data   betaai(  3, 26)/  -23.000000d0/
      data      fgi(  1, 26)/   -8.570000d0/
      data      fgi(  2, 26)/   -7.270000d0/
      data      fgi(  3, 26)/   -5.420000d0/
      data      fgi(  4, 26)/   -4.080000d0/
      data      fgi(  5, 26)/  -12.310000d0/
      data      fgi(  6, 26)/   -8.530000d0/
      data      fgi(  7, 26)/ -100.000000d0/
      data      fgi(  8, 26)/    0.314300d0/
      data      fgi(  9, 26)/    0.685700d0/
      data      fgi( 11, 26)/    6.540000d0/
      data      fgi( 12, 26)/    7.680000d0/
      data      fgi( 13, 26)/   11.520000d0/
      data      fgi( 14, 26)/    2.020937d0/
      data      fgi( 15, 26)/    0.929879d0/
      data      fgi( 16, 26)/    0.823253d0/
      data      fgi( 17, 26)/    0.303760d0/
      data      fgi( 18, 26)/    0.622399d0/
      data      fgi( 19, 26)/    0.436423d0/
      data      fgi( 20, 26)/    7.563016d0/
      data      fgi( 21, 26)/    4.760980d0/
      data      fgi( 22, 26)/    0.937362d0/
      data      fgi( 23, 26)/    1.382743d0/
      data      fgi( 24, 26)/    0.616758d0/
!
!                    Data for Element  27
!
      data        isoki( 27)/    1         /
      data        nbfai( 27)/    9         /
      data      zcoreai( 27)/    9.000000d0/
      data        zetai( 27)/    1.423000d0/
      data   zetadi(  1, 27)/    6.386000d0/
      data   zetadi(  2, 27)/    2.745000d0/
      data  zetawti(  1, 27)/    0.413330d0/
      data  zetawti(  2, 27)/    0.712620d0/
      data   betaai(  1, 27)/   -1.000000d0/
      data   betaai(  2, 27)/   -1.000000d0/
      data   betaai(  3, 27)/  -31.000000d0/
      data      fgi(  1, 27)/   -8.760000d0/
      data      fgi(  2, 27)/   -7.380000d0/
      data      fgi(  3, 27)/   -5.480000d0/
      data      fgi(  4, 27)/   -4.020000d0/
      data      fgi(  5, 27)/  -13.300000d0/
      data      fgi(  6, 27)/   -9.100000d0/
      data      fgi(  7, 27)/ -100.000000d0/
      data      fgi(  8, 27)/    0.206500d0/
      data      fgi(  9, 27)/    0.793500d0/
      data      fgi( 11, 27)/    6.960000d0/
      data      fgi( 12, 27)/    8.160000d0/
      data      fgi( 13, 27)/   12.320000d0/
      data      fgi( 14, 27)/    2.814434d0/
      data      fgi( 15, 27)/    0.991871d0/
      data      fgi( 16, 27)/    0.786058d0/
      data      fgi( 17, 27)/    0.393029d0/
      data      fgi( 18, 27)/    0.779859d0/
      data      fgi( 19, 27)/    0.280204d0/
      data      fgi( 20, 27)/    7.996959d0/
      data      fgi( 21, 27)/    5.963624d0/
      data      fgi( 22, 27)/    0.927827d0/
      data      fgi( 23, 27)/    1.392741d0/
      data      fgi( 24, 27)/    0.606753d0/
!
!                    Data for Element  28
!
      data        isoki( 28)/    1         /
      data        nbfai( 28)/    9         /
      data      zcoreai( 28)/   10.000000d0/
      data        zetai( 28)/    1.473000d0/
      data   zetadi(  1, 28)/    6.706000d0/
      data   zetadi(  2, 28)/    2.874000d0/
      data  zetawti(  1, 28)/    0.421200d0/
      data  zetawti(  2, 28)/    0.706580d0/
      data   betaai(  1, 28)/   -1.000000d0/
      data   betaai(  2, 28)/   -1.000000d0/
      data   betaai(  3, 28)/  -35.000000d0/
      data      fgi(  1, 28)/   -8.940000d0/
      data      fgi(  2, 28)/   -7.510000d0/
      data      fgi(  3, 28)/   -5.410000d0/
      data      fgi(  4, 28)/   -3.930000d0/
      data      fgi(  5, 28)/  -14.460000d0/
      data      fgi(  6, 28)/   -9.790000d0/
      data      fgi(  7, 28)/ -100.000000d0/
      data      fgi(  8, 28)/    0.142100d0/
      data      fgi(  9, 28)/    0.857900d0/
      data      fgi( 11, 28)/    7.340000d0/
      data      fgi( 12, 28)/    8.610000d0/
      data      fgi( 13, 28)/   13.100000d0/
      data      fgi( 14, 28)/    2.405287d0/
      data      fgi( 15, 28)/    1.053863d0/
      data      fgi( 16, 28)/    0.830692d0/
      data      fgi( 17, 28)/    0.373191d0/
      data      fgi( 18, 28)/    0.750102d0/
      data      fgi( 19, 28)/    0.402948d0/
      data      fgi( 20, 28)/    9.893912d0/
      data      fgi( 21, 28)/    6.608340d0/
      data      fgi( 22, 28)/    0.914218d0/
      data      fgi( 23, 28)/    1.397009d0/
      data      fgi( 24, 28)/    0.594082d0/
!
!                    Data for Element  29
!
      data        isoki( 29)/    1         /
      data        nbfai( 29)/    9         /
      data      zcoreai( 29)/   11.000000d0/
      data        zetai( 29)/    1.482000d0/
      data   zetadi(  1, 29)/    6.795000d0/
      data   zetadi(  2, 29)/    2.765000d0/
      data  zetawti(  1, 29)/    0.447290d0/
      data  zetawti(  2, 29)/    0.696830d0/
      data   betaai(  1, 29)/   -1.000000d0/
      data   betaai(  2, 29)/   -1.000000d0/
      data   betaai(  3, 29)/  -40.000000d0/
      data      fgi(  1, 29)/   -9.130000d0/
      data      fgi(  2, 29)/   -7.690000d0/
      data      fgi(  3, 29)/   -5.180000d0/
      data      fgi(  4, 29)/   -3.840000d0/
      data      fgi(  5, 29)/  -15.870000d0/
      data      fgi(  6, 29)/  -10.670000d0/
      data      fgi(  7, 29)/ -100.000000d0/
      data      fgi(  8, 29)/    0.095600d0/
      data      fgi(  9, 29)/    0.904400d0/
      data      fgi( 11, 29)/    7.680000d0/
      data      fgi( 12, 29)/    9.010000d0/
      data      fgi( 13, 29)/   13.840000d0/
      data      fgi( 14, 29)/    2.566466d0/
      data      fgi( 15, 29)/    1.115855d0/
      data      fgi( 16, 29)/    0.552968d0/
      data      fgi( 17, 29)/    0.696789d0/
      data      fgi( 18, 29)/    1.326627d0/
      data      fgi( 19, 29)/    0.859208d0/
      data      fgi( 20, 29)/   10.660133d0/
      data      fgi( 21, 29)/    7.186725d0/
      data      fgi( 22, 29)/    0.815680d0/
      data      fgi( 23, 29)/    1.301797d0/
      data      fgi( 24, 29)/    0.521822d0/
!
!                    Data for Element  30
!
      data        isoki( 30)/    1         /
      data        nbfai( 30)/    4         /
      data      zcoreai( 30)/    2.000000d0/
      data        zetai( 30)/    1.509000d0/
      data   betaai(  1, 30)/  -10.000000d0/
      data   betaai(  2, 30)/  -10.000000d0/
      data      fgi(  1, 30)/   -9.360000d0/
      data      fgi(  2, 30)/   -9.360000d0/
      data      fgi(  3, 30)/   -4.770000d0/
      data      fgi(  4, 30)/   -4.770000d0/
      data      fgi(  8, 30)/    1.000000d0/
      data      fgi( 11, 30)/    7.980000d0/
      data      fgi( 12, 30)/    7.980000d0/
      data      fgi( 13, 30)/   14.550000d0/
      data      fgi( 14, 30)/    2.529271d0/
      data      fgi( 15, 30)/    1.177847d0/
!
!                    Data for Element  31
!
      data        isoki( 31)/    0         /
      data        nbfai( 31)/    4         /
      data      zcoreai( 31)/    3.000000d0/
      data   betaai(  1, 31)/ -100.000000d0/
      data   betaai(  2, 31)/ -100.000000d0/
      data      fgi(  1, 31)/ -100.000000d0/
      data      fgi(  3, 31)/ -100.000000d0/
      data      fgi(  8, 31)/    1.000000d0/
!
!                    Data for Element  32
!
      data        isoki( 32)/    0         /
      data        nbfai( 32)/    4         /
      data      zcoreai( 32)/    4.000000d0/
      data   betaai(  1, 32)/ -100.000000d0/
      data   betaai(  2, 32)/ -100.000000d0/
      data      fgi(  1, 32)/ -100.000000d0/
      data      fgi(  3, 32)/ -100.000000d0/
      data      fgi(  8, 32)/    1.000000d0/
!
!                    Data for Element  33
!
      data        isoki( 33)/    0         /
      data        nbfai( 33)/    4         /
      data      zcoreai( 33)/    5.000000d0/
      data        zetai( 33)/    2.921000d0/
      data   betaai(  1, 33)/  -10.000000d0/
      data   betaai(  2, 33)/  -10.000000d0/
      data      fgi(  1, 33)/  -18.180000d0/
      data      fgi(  3, 33)/   -9.190000d0/
      data      fgi(  8, 33)/    1.000000d0/
      data      fgi( 14, 33)/    5.036471d0/
      data      fgi( 15, 33)/    3.802126d0/
!
!                    Data for Element  34
!
      data        isoki( 34)/    0         /
      data        nbfai( 34)/    4         /
      data      zcoreai( 34)/    6.000000d0/
      data        zetai( 34)/    2.439000d0/
      data   betaai(  1, 34)/  -11.670000d0/
      data   betaai(  2, 34)/  -11.670000d0/
      data      fgi(  1, 34)/  -20.950000d0/
      data      fgi(  3, 34)/  -10.370000d0/
      data      fgi(  8, 34)/    1.000000d0/
      data      fgi( 14, 34)/    5.630697d0/
      data      fgi( 15, 34)/    4.230938d0/
!
!                    Data for Element  35
!
      data        isoki( 35)/    1         /
      data        nbfai( 35)/    4         /
      data      zcoreai( 35)/    7.000000d0/
      data        zetai( 35)/    2.638000d0/
      data   betaai(  1, 35)/   -8.000000d0/
      data   betaai(  2, 35)/   -8.000000d0/
      data      fgi(  1, 35)/  -23.940000d0/
      data      fgi(  3, 35)/  -12.440000d0/
      data      fgi(  8, 35)/    1.000000d0/
      data      fgi( 11, 35)/    9.080000d0/
      data      fgi( 12, 35)/    9.080000d0/
      data      fgi( 14, 35)/    6.141120d0/
      data      fgi( 15, 35)/    4.608700d0/
!
!                    Data for Element  36
!
      data        isoki( 36)/    0         /
      data        nbfai( 36)/    4         /
      data      zcoreai( 36)/    8.000000d0/
      data   betaai(  1, 36)/ -100.000000d0/
      data   betaai(  2, 36)/ -100.000000d0/
      data      fgi(  1, 36)/ -100.000000d0/
      data      fgi(  3, 36)/ -100.000000d0/
      data      fgi(  8, 36)/    1.000000d0/
      data      fgi( 13, 36)/    2.480000d0/
!
!                    Data for Element  37
!
      data        isoki( 37)/    0         /
      data        nbfai( 37)/    4         /
      data      zcoreai( 37)/    1.000000d0/
      data   betaai(  1, 37)/ -100.000000d0/
      data   betaai(  2, 37)/ -100.000000d0/
      data      fgi(  1, 37)/   -4.430000d0/
      data      fgi(  3, 37)/   -2.650000d0/
      data      fgi(  8, 37)/    1.000000d0/
      data      fgi( 11, 37)/    1.500000d0/
      data      fgi( 12, 37)/    1.500000d0/
      data      fgi( 13, 37)/    3.380000d0/
!
!                    Data for Element  38
!
      data        isoki( 38)/    0         /
      data        nbfai( 38)/    9         /
      data      zcoreai( 38)/    2.000000d0/
      data  zetawti(  1, 38)/    1.000000d0/
      data   betaai(  1, 38)/ -100.000000d0/
      data   betaai(  2, 38)/ -100.000000d0/
      data   betaai(  3, 38)/ -100.000000d0/
      data      fgi(  1, 38)/   -5.840000d0/
      data      fgi(  2, 38)/   -5.190000d0/
      data      fgi(  3, 38)/   -3.760000d0/
      data      fgi(  4, 38)/   -3.160000d0/
      data      fgi(  5, 38)/   -3.660000d0/
      data      fgi(  6, 38)/   -3.660000d0/
      data      fgi(  7, 38)/   -2.490000d0/
      data      fgi(  8, 38)/    0.939000d0/
      data      fgi(  9, 38)/    0.055000d0/
      data      fgi( 10, 38)/    0.007000d0/
      data      fgi( 11, 38)/    2.050000d0/
      data      fgi( 12, 38)/    1.970000d0/
      data      fgi( 13, 38)/    4.280000d0/
!
!                    Data for Element  39
!
      data        isoki( 39)/    1         /
      data        nbfai( 39)/    9         /
      data      zcoreai( 39)/    3.000000d0/
      data        zetai( 39)/    1.275000d0/
      data   zetadi(  1, 39)/    3.837000d0/
      data   zetadi(  2, 39)/    1.739000d0/
      data  zetawti(  1, 39)/    0.286550d0/
      data  zetawti(  2, 39)/    0.824680d0/
      data   betaai(  1, 39)/   -1.000000d0/
      data   betaai(  2, 39)/   -1.000000d0/
      data   betaai(  3, 39)/   -7.000000d0/
      data      fgi(  1, 39)/   -6.550000d0/
      data      fgi(  2, 39)/   -5.850000d0/
      data      fgi(  3, 39)/   -4.130000d0/
      data      fgi(  4, 39)/   -3.580000d0/
      data      fgi(  5, 39)/   -6.610000d0/
      data      fgi(  6, 39)/   -4.740000d0/
      data      fgi(  7, 39)/   -3.220000d0/
      data      fgi(  8, 39)/    0.920000d0/
      data      fgi(  9, 39)/    0.073000d0/
      data      fgi( 10, 39)/    0.007000d0/
      data      fgi( 11, 39)/    2.570000d0/
      data      fgi( 12, 39)/    2.410000d0/
      data      fgi( 13, 39)/    5.220000d0/
      data      fgi( 14, 39)/    2.846438d0/
      data      fgi( 15, 39)/    2.231965d0/
      data      fgi( 16, 39)/    1.250107d0/
      data      fgi( 17, 39)/    1.607476d0/
      data      fgi( 18, 39)/    1.833888d0/
      data      fgi( 19, 39)/    1.011633d0/
      data      fgi( 20, 39)/    3.834798d0/
      data      fgi( 21, 39)/    2.562213d0/
      data      fgi( 22, 39)/    1.942611d0/
      data      fgi( 23, 39)/    1.968347d0/
      data      fgi( 24, 39)/    1.458400d0/
!
!                    Data for Element  40
!
      data        isoki( 40)/    1         /
      data        nbfai( 40)/    9         /
      data      zcoreai( 40)/    4.000000d0/
      data        zetai( 40)/    1.337000d0/
      data   zetadi(  1, 40)/    3.639000d0/
      data   zetadi(  2, 40)/    1.804000d0/
      data  zetawti(  1, 40)/    0.399240d0/
      data  zetawti(  2, 40)/    0.713810d0/
      data   betaai(  1, 40)/   -1.000000d0/
      data   betaai(  2, 40)/   -1.000000d0/
      data   betaai(  3, 40)/  -10.000000d0/
      data      fgi(  1, 40)/   -7.160000d0/
      data      fgi(  2, 40)/   -6.400000d0/
      data      fgi(  3, 40)/   -4.430000d0/
      data      fgi(  4, 40)/   -3.920000d0/
      data      fgi(  5, 40)/   -8.080000d0/
      data      fgi(  6, 40)/   -5.790000d0/
      data      fgi(  7, 40)/   -3.950000d0/
      data      fgi(  8, 40)/    0.884000d0/
      data      fgi(  9, 40)/    0.108000d0/
      data      fgi( 10, 40)/    0.008000d0/
      data      fgi( 11, 40)/    3.060000d0/
      data      fgi( 12, 40)/    2.930000d0/
      data      fgi( 13, 40)/    6.120000d0/
      data      fgi( 14, 40)/    2.984853d0/
      data      fgi( 15, 40)/    2.340500d0/
      data      fgi( 16, 40)/    1.064414d0/
      data      fgi( 17, 40)/    1.368699d0/
      data      fgi( 18, 40)/    1.759978d0/
      data      fgi( 19, 40)/    0.861363d0/
      data      fgi( 20, 40)/    4.369791d0/
      data      fgi( 21, 40)/    2.919668d0/
      data      fgi( 22, 40)/    1.777103d0/
      data      fgi( 23, 40)/    1.890505d0/
      data      fgi( 24, 40)/    1.315937d0/
!
!                    Data for Element  41
!
      data        isoki( 41)/    1         /
      data        nbfai( 41)/    9         /
      data      zcoreai( 41)/    5.000000d0/
      data        zetai( 41)/    1.393000d0/
      data   zetadi(  1, 41)/    3.774100d0/
      data   zetadi(  2, 41)/    1.925100d0/
      data  zetawti(  1, 41)/    0.447750d0/
      data  zetawti(  2, 41)/    0.663010d0/
      data   betaai(  1, 41)/   -1.000000d0/
      data   betaai(  2, 41)/   -1.000000d0/
      data   betaai(  3, 41)/  -13.000000d0/
      data      fgi(  1, 41)/   -7.690000d0/
      data      fgi(  2, 41)/   -6.840000d0/
      data      fgi(  3, 41)/   -4.670000d0/
      data      fgi(  4, 41)/   -4.150000d0/
      data      fgi(  5, 41)/   -9.500000d0/
      data      fgi(  6, 41)/   -6.820000d0/
      data      fgi(  7, 41)/   -4.700000d0/
      data      fgi(  8, 41)/    0.806000d0/
      data      fgi(  9, 41)/    0.187000d0/
      data      fgi( 10, 41)/    0.007000d0/
      data      fgi( 11, 41)/    3.520000d0/
      data      fgi( 12, 41)/    3.450000d0/
      data      fgi( 13, 41)/    7.010000d0/
      data      fgi( 14, 41)/    3.109873d0/
      data      fgi( 15, 41)/    2.438532d0/
      data      fgi( 16, 41)/    0.950646d0/
      data      fgi( 17, 41)/    1.222408d0/
      data      fgi( 18, 41)/    1.717580d0/
      data      fgi( 19, 41)/    0.769298d0/
      data      fgi( 20, 41)/    4.808812d0/
      data      fgi( 21, 41)/    3.212999d0/
      data      fgi( 22, 41)/    1.673072d0/
      data      fgi( 23, 41)/    1.842142d0/
      data      fgi( 24, 41)/    1.226002d0/
!
!                    Data for Element  42
!
      data        isoki( 42)/    1         /
      data        nbfai( 42)/    9         /
      data      zcoreai( 42)/    6.000000d0/
      data        zetai( 42)/    1.440000d0/
      data   zetadi(  1, 42)/    3.954000d0/
      data   zetadi(  2, 42)/    2.047000d0/
      data  zetawti(  1, 42)/    0.480130d0/
      data  zetawti(  2, 42)/    0.628920d0/
      data   betaai(  1, 42)/   -1.000000d0/
      data   betaai(  2, 42)/   -1.000000d0/
      data   betaai(  3, 42)/  -15.000000d0/
      data      fgi(  1, 42)/   -8.130000d0/
      data      fgi(  2, 42)/   -7.180000d0/
      data      fgi(  3, 42)/   -4.870000d0/
      data      fgi(  4, 42)/   -4.290000d0/
      data      fgi(  5, 42)/  -10.850000d0/
      data      fgi(  6, 42)/   -7.820000d0/
      data      fgi(  7, 42)/   -5.460000d0/
      data      fgi(  8, 42)/    0.616000d0/
      data      fgi(  9, 42)/    0.382000d0/
      data      fgi( 10, 42)/    0.002000d0/
      data      fgi( 11, 42)/    3.980000d0/
      data      fgi( 12, 42)/    3.990000d0/
      data      fgi( 13, 42)/    7.880000d0/
      data      fgi( 14, 42)/    3.214801d0/
      data      fgi( 15, 42)/    2.520808d0/
      data      fgi( 16, 42)/    0.872601d0/
      data      fgi( 17, 42)/    1.122052d0/
      data      fgi( 18, 42)/    1.688591d0/
      data      fgi( 19, 42)/    0.706140d0/
      data      fgi( 20, 42)/    5.170238d0/
      data      fgi( 21, 42)/    3.454485d0/
      data      fgi( 22, 42)/    1.599116d0/
      data      fgi( 23, 42)/    1.807038d0/
      data      fgi( 24, 42)/    1.162179d0/
!
!                    Data for Element  43
!
      data        isoki( 43)/    1         /
      data        nbfai( 43)/    9         /
      data      zcoreai( 43)/    7.000000d0/
      data        zetai( 43)/    1.461000d0/
      data   zetadi(  1, 43)/    4.124000d0/
      data   zetadi(  2, 43)/    2.155000d0/
      data  zetawti(  1, 43)/    0.512410d0/
      data  zetawti(  2, 43)/    0.595420d0/
      data   betaai(  1, 43)/   -1.000000d0/
      data   betaai(  2, 43)/   -1.000000d0/
      data   betaai(  3, 43)/  -17.000000d0/
      data      fgi(  1, 43)/   -8.490000d0/
      data      fgi(  2, 43)/   -7.420000d0/
      data      fgi(  3, 43)/   -5.000000d0/
      data      fgi(  4, 43)/   -4.340000d0/
      data      fgi(  5, 43)/  -12.130000d0/
      data      fgi(  6, 43)/   -8.800000d0/
      data      fgi(  7, 43)/   -6.210000d0/
      data      fgi(  8, 43)/    0.242000d0/
      data      fgi(  9, 43)/    0.736000d0/
      data      fgi( 10, 43)/    0.022000d0/
      data      fgi( 11, 43)/    4.400000d0/
      data      fgi( 12, 43)/    4.540000d0/
      data      fgi( 13, 43)/    8.730000d0/
      data      fgi( 14, 43)/    3.261683d0/
      data      fgi( 15, 43)/    2.557570d0/
      data      fgi( 16, 43)/    0.792260d0/
      data      fgi( 17, 43)/    1.018744d0/
      data      fgi( 18, 43)/    1.635130d0/
      data      fgi( 19, 43)/    0.641126d0/
      data      fgi( 20, 43)/    5.431609d0/
      data      fgi( 21, 43)/    3.629120d0/
      data      fgi( 22, 43)/    1.507582d0/
      data      fgi( 23, 43)/    1.744914d0/
      data      fgi( 24, 43)/    1.087087d0/
!
!                    Data for Element  44
!
      data        isoki( 44)/    1         /
      data        nbfai( 44)/    9         /
      data      zcoreai( 44)/    8.000000d0/
      data        zetai( 44)/    1.470000d0/
      data   zetadi(  1, 44)/    4.259000d0/
      data   zetadi(  2, 44)/    2.094000d0/
      data  zetawti(  1, 44)/    0.534230d0/
      data  zetawti(  2, 44)/    0.592710d0/
      data   betaai(  1, 44)/   -2.000000d0/
      data   betaai(  2, 44)/   -2.000000d0/
      data   betaai(  3, 44)/  -20.000000d0/
      data      fgi(  1, 44)/   -8.750000d0/
      data      fgi(  2, 44)/   -7.550000d0/
      data      fgi(  3, 44)/   -5.070000d0/
      data      fgi(  4, 44)/   -4.300000d0/
      data      fgi(  5, 44)/  -13.380000d0/
      data      fgi(  6, 44)/   -9.770000d0/
      data      fgi(  7, 44)/   -6.980000d0/
      data      fgi(  8, 44)/    0.016000d0/
      data      fgi(  9, 44)/    0.717000d0/
      data      fgi( 10, 44)/    0.266000d0/
      data      fgi( 11, 44)/    4.810000d0/
      data      fgi( 12, 44)/    5.120000d0/
      data      fgi( 13, 44)/    9.550000d0/
      data      fgi( 14, 44)/    2.479677d0/
      data      fgi( 15, 44)/    1.120070d0/
      data      fgi( 16, 44)/    0.990011d0/
      data      fgi( 17, 44)/    0.330045d0/
      data      fgi( 18, 44)/    1.440073d0/
      data      fgi( 19, 44)/    0.180025d0/
      data      fgi( 20, 44)/    6.170305d0/
      data      fgi( 21, 44)/    4.630177d0/
      data      fgi( 22, 44)/    1.431276d0/
      data      fgi( 23, 44)/    1.687971d0/
      data      fgi( 24, 44)/    1.025590d0/
!
!                    Data for Element  45
!
      data        isoki( 45)/    1         /
      data        nbfai( 45)/    9         /
      data      zcoreai( 45)/    9.000000d0/
      data        zetai( 45)/    1.482000d0/
      data   zetadi(  1, 45)/    4.485000d0/
      data   zetadi(  2, 45)/    2.217000d0/
      data  zetawti(  1, 45)/    0.544090d0/
      data  zetawti(  2, 45)/    0.581370d0/
      data   betaai(  1, 45)/   -1.000000d0/
      data   betaai(  2, 45)/   -1.000000d0/
      data   betaai(  3, 45)/  -21.000000d0/
      data      fgi(  1, 45)/   -8.930000d0/
      data      fgi(  2, 45)/   -7.580000d0/
      data      fgi(  3, 45)/   -5.070000d0/
      data      fgi(  4, 45)/   -4.170000d0/
      data      fgi(  5, 45)/  -14.540000d0/
      data      fgi(  6, 45)/  -10.700000d0/
      data      fgi(  7, 45)/   -7.750000d0/
      data      fgi(  8, 45)/    0.003000d0/
      data      fgi(  9, 45)/    0.319000d0/
      data      fgi( 10, 45)/    0.678000d0/
      data      fgi( 11, 45)/    5.190000d0/
      data      fgi( 12, 45)/    5.710000d0/
      data      fgi( 13, 45)/   10.360000d0/
      data      fgi( 14, 45)/    3.308566d0/
      data      fgi( 15, 45)/    2.594332d0/
      data      fgi( 16, 45)/    0.650041d0/
      data      fgi( 17, 45)/    0.835868d0/
      data      fgi( 18, 45)/    1.518541d0/
      data      fgi( 19, 45)/    0.526037d0/
      data      fgi( 20, 45)/    5.864504d0/
      data      fgi( 21, 45)/    3.918358d0/
      data      fgi( 22, 45)/    1.328619d0/
      data      fgi( 23, 45)/    1.608894d0/
      data      fgi( 24, 45)/    0.943435d0/
!
!                    Data for Element  46
!
      data        isoki( 46)/    0         /
      data        nbfai( 46)/    9         /
      data      zcoreai( 46)/   10.000000d0/
      data        zetai( 46)/    1.497000d0/
      data   zetadi(  1, 46)/    4.723000d0/
      data   zetadi(  2, 46)/    2.342000d0/
      data  zetawti(  1, 46)/    0.549020d0/
      data  zetawti(  2, 46)/    0.575610d0/
      data   betaai(  1, 46)/   -1.000000d0/
      data   betaai(  2, 46)/   -1.000000d0/
      data   betaai(  3, 46)/  -21.500000d0/
      data      fgi(  1, 46)/   -9.020000d0/
      data      fgi(  2, 46)/   -7.500000d0/
      data      fgi(  3, 46)/   -5.020000d0/
      data      fgi(  4, 46)/   -3.940000d0/
      data      fgi(  5, 46)/  -15.650000d0/
      data      fgi(  6, 46)/  -11.620000d0/
      data      fgi(  7, 46)/   -8.530000d0/
      data      fgi(  8, 46)/    0.007000d0/
      data      fgi(  9, 46)/    0.121000d0/
      data      fgi( 10, 46)/    0.872000d0/
      data      fgi( 12, 46)/    6.330000d0/
      data      fgi( 13, 46)/   11.150000d0/
      data      fgi( 14, 46)/    3.342053d0/
      data      fgi( 15, 46)/    2.620590d0/
      data      fgi( 16, 46)/    0.583298d0/
      data      fgi( 17, 46)/    0.750046d0/
      data      fgi( 18, 46)/    1.460972d0/
      data      fgi( 19, 46)/    0.472026d0/
      data      fgi( 20, 46)/    6.121791d0/
      data      fgi( 21, 46)/    4.090263d0/
      data      fgi( 22, 46)/    1.240469d0/
      data      fgi( 23, 46)/    1.540221d0/
      data      fgi( 24, 46)/    0.873138d0/
!
!                    Data for Element  47
!
      data        isoki( 47)/    1         /
      data        nbfai( 47)/    9         /
      data      zcoreai( 47)/   11.000000d0/
      data        zetai( 47)/    1.507000d0/
      data   zetadi(  1, 47)/    4.988960d0/
      data   zetadi(  2, 47)/    2.583740d0/
      data  zetawti(  1, 47)/    0.557630d0/
      data  zetawti(  2, 47)/    0.553590d0/
      data   betaai(  1, 47)/   -2.500000d0/
      data   betaai(  2, 47)/   -2.500000d0/
      data   betaai(  3, 47)/  -30.000000d0/
      data      fgi(  2, 47)/   -7.310000d0/
      data      fgi(  4, 47)/   -3.620000d0/
      data      fgi(  6, 47)/  -12.510000d0/
      data      fgi(  9, 47)/    1.000000d0/
      data      fgi( 11, 47)/    6.000000d0/
      data      fgi( 12, 47)/    6.960000d0/
      data      fgi( 13, 47)/    9.550000d0/
      data      fgi( 14, 47)/    3.027900d0/
      data      fgi( 15, 47)/    2.374200d0/
      data      fgi( 16, 47)/    0.386800d0/
      data      fgi( 17, 47)/    0.508300d0/
      data      fgi( 18, 47)/    1.024000d0/
      data      fgi( 19, 47)/    0.308600d0/
      data      fgi( 20, 47)/    5.965000d0/
      data      fgi( 21, 47)/    3.828500d0/
!
!                    Data for Element  48
!
      data        isoki( 48)/    0         /
      data        nbfai( 48)/    4         /
      data      zcoreai( 48)/    2.000000d0/
      data        zetai( 48)/    1.699000d0/
      data   betaai(  1, 48)/   -1.000000d0/
      data   betaai(  2, 48)/   -1.000000d0/
      data      fgi(  1, 48)/   -8.940000d0/
      data      fgi(  2, 48)/   -8.940000d0/
      data      fgi(  3, 48)/   -4.750000d0/
      data      fgi(  4, 48)/   -4.750000d0/
      data      fgi(  8, 48)/    1.000000d0/
      data      fgi( 14, 48)/    3.793018d0/
      data      fgi( 15, 48)/    2.974203d0/
!
!                    Data for Element  49
!
      data        isoki( 49)/    0         /
      data        nbfai( 49)/    4         /
      data      zcoreai( 49)/    3.000000d0/
      data   betaai(  1, 49)/ -100.000000d0/
      data   betaai(  2, 49)/ -100.000000d0/
      data      fgi(  1, 49)/ -100.000000d0/
      data      fgi(  3, 49)/ -100.000000d0/
      data      fgi(  8, 49)/    1.000000d0/
!
!                    Data for Element  50
!
      data        isoki( 50)/    0         /
      data        nbfai( 50)/    4         /
      data      zcoreai( 50)/    4.000000d0/
      data   betaai(  1, 50)/ -100.000000d0/
      data   betaai(  2, 50)/ -100.000000d0/
      data      fgi(  1, 50)/ -100.000000d0/
      data      fgi(  3, 50)/ -100.000000d0/
      data      fgi(  8, 50)/    1.000000d0/
!
!                    Data for Element  51
!
      data        isoki( 51)/    0         /
      data        nbfai( 51)/    4         /
      data      zcoreai( 51)/    5.000000d0/
      data   betaai(  1, 51)/ -100.000000d0/
      data   betaai(  2, 51)/ -100.000000d0/
      data      fgi(  1, 51)/ -100.000000d0/
      data      fgi(  3, 51)/ -100.000000d0/
      data      fgi(  8, 51)/    1.000000d0/
!
!                    Data for Element  52
!
      data        isoki( 52)/    0         /
      data        nbfai( 52)/    4         /
      data      zcoreai( 52)/    6.000000d0/
      data   betaai(  1, 52)/ -100.000000d0/
      data   betaai(  2, 52)/ -100.000000d0/
      data      fgi(  1, 52)/ -100.000000d0/
      data      fgi(  3, 52)/ -100.000000d0/
      data      fgi(  8, 52)/    1.000000d0/
!
!                    Data for Element  53
!
      data        isoki( 53)/    1         /
      data        nbfai( 53)/    4         /
      data      zcoreai( 53)/    7.000000d0/
      data        zetai( 53)/    3.341000d0/
      data   betaai(  1, 53)/   -8.000000d0/
      data   betaai(  2, 53)/   -8.000000d0/
      data      fgi(  1, 53)/  -20.840000d0/
      data      fgi(  3, 53)/  -11.210000d0/
      data      fgi(  8, 53)/    1.000000d0/
      data      fgi( 11, 53)/    8.150000d0/
      data      fgi( 12, 53)/    8.150000d0/
      data      fgi( 14, 53)/    4.886521d0/
      data      fgi( 15, 53)/    3.842482d0/
!
!                    Data for Element  54
!
      data        isoki( 54)/    0         /
      data        nbfai( 54)/    4         /
      data      zcoreai( 54)/    8.000000d0/
      data   betaai(  1, 54)/ -100.000000d0/
      data   betaai(  2, 54)/ -100.000000d0/
      data      fgi(  1, 54)/ -100.000000d0/
      data      fgi(  3, 54)/ -100.000000d0/
      data      fgi(  8, 54)/    1.000000d0/
!
!                    Data for Element  55
!
      data        isoki( 55)/    0         /
      data        nbfai( 55)/    4         /
      data      zcoreai( 55)/    1.000000d0/
      data   betaai(  1, 55)/ -100.000000d0/
      data   betaai(  2, 55)/ -100.000000d0/
      data      fgi(  1, 55)/ -100.000000d0/
      data      fgi(  3, 55)/ -100.000000d0/
      data      fgi(  8, 55)/    1.000000d0/
!
!                    Data for Element  56
!
      data        isoki( 56)/    0         /
      data        nbfai( 56)/    9         /
      data      zcoreai( 56)/    2.000000d0/
      data   betaai(  1, 56)/ -100.000000d0/
      data   betaai(  2, 56)/ -100.000000d0/
      data      fgi(  1, 56)/ -100.000000d0/
      data      fgi(  3, 56)/ -100.000000d0/
      data      fgi(  8, 56)/    1.000000d0/
!
!                    Data for Element  57
!
      data        isoki( 57)/    0         /
      data        nbfai( 57)/    9         /
      data      zcoreai( 57)/    3.000000d0/
      data   betaai(  1, 57)/ -100.000000d0/
      data   betaai(  2, 57)/ -100.000000d0/
      data      fgi(  1, 57)/ -100.000000d0/
      data      fgi(  3, 57)/ -100.000000d0/
      data      fgi(  8, 57)/    1.000000d0/
!
!                    Data for Element  58
!
      data        isoki( 58)/    0         /
      data        nbfai( 58)/    9         /
      data      zcoreai( 58)/    4.000000d0/
      data   betaai(  1, 58)/ -100.000000d0/
      data   betaai(  2, 58)/ -100.000000d0/
      data      fgi(  1, 58)/ -100.000000d0/
      data      fgi(  3, 58)/ -100.000000d0/
      data      fgi(  8, 58)/    1.000000d0/
!
!                    Data for Element  59
!
      data        isoki( 59)/    0         /
      data        nbfai( 59)/    9         /
      data      zcoreai( 59)/    5.000000d0/
      data   betaai(  1, 59)/ -100.000000d0/
      data   betaai(  2, 59)/ -100.000000d0/
      data      fgi(  1, 59)/ -100.000000d0/
      data      fgi(  3, 59)/ -100.000000d0/
      data      fgi(  8, 59)/    1.000000d0/
!
!                    Data for Element  60
!
      data        isoki( 60)/    0         /
      data        nbfai( 60)/    9         /
      data      zcoreai( 60)/    6.000000d0/
      data   betaai(  1, 60)/ -100.000000d0/
      data   betaai(  2, 60)/ -100.000000d0/
      data      fgi(  1, 60)/ -100.000000d0/
      data      fgi(  3, 60)/ -100.000000d0/
      data      fgi(  8, 60)/    1.000000d0/
!
!                    Data for Element  61
!
      data        isoki( 61)/    0         /
      data        nbfai( 61)/    9         /
      data      zcoreai( 61)/    7.000000d0/
      data   betaai(  1, 61)/ -100.000000d0/
      data   betaai(  2, 61)/ -100.000000d0/
      data      fgi(  1, 61)/ -100.000000d0/
      data      fgi(  3, 61)/ -100.000000d0/
      data      fgi(  8, 61)/    1.000000d0/
!
!                    Data for Element  62
!
      data        isoki( 62)/    0         /
      data        nbfai( 62)/    9         /
      data      zcoreai( 62)/    8.000000d0/
      data   betaai(  1, 62)/ -100.000000d0/
      data   betaai(  2, 62)/ -100.000000d0/
      data      fgi(  1, 62)/ -100.000000d0/
      data      fgi(  3, 62)/ -100.000000d0/
      data      fgi(  8, 62)/    1.000000d0/
!
!                    Data for Element  63
!
      data        isoki( 63)/    0         /
      data        nbfai( 63)/    9         /
      data      zcoreai( 63)/    9.000000d0/
      data   betaai(  1, 63)/ -100.000000d0/
      data   betaai(  2, 63)/ -100.000000d0/
      data      fgi(  1, 63)/ -100.000000d0/
      data      fgi(  3, 63)/ -100.000000d0/
      data      fgi(  8, 63)/    1.000000d0/
!
!                    Data for Element  64
!
      data        isoki( 64)/    0         /
      data        nbfai( 64)/    9         /
      data      zcoreai( 64)/   10.000000d0/
      data   betaai(  1, 64)/ -100.000000d0/
      data   betaai(  2, 64)/ -100.000000d0/
      data      fgi(  1, 64)/ -100.000000d0/
      data      fgi(  3, 64)/ -100.000000d0/
      data      fgi(  8, 64)/    1.000000d0/
!
!                    Data for Element  65
!
      data        isoki( 65)/    0         /
      data        nbfai( 65)/    9         /
      data      zcoreai( 65)/   11.000000d0/
      data   betaai(  1, 65)/ -100.000000d0/
      data   betaai(  2, 65)/ -100.000000d0/
      data      fgi(  1, 65)/ -100.000000d0/
      data      fgi(  3, 65)/ -100.000000d0/
      data      fgi(  8, 65)/    1.000000d0/
!
!                    Data for Element  66
!
      data        isoki( 66)/    0         /
      data        nbfai( 66)/    9         /
      data      zcoreai( 66)/   12.000000d0/
      data   betaai(  1, 66)/ -100.000000d0/
      data   betaai(  2, 66)/ -100.000000d0/
      data      fgi(  1, 66)/ -100.000000d0/
      data      fgi(  3, 66)/ -100.000000d0/
      data      fgi(  8, 66)/    1.000000d0/
!
!                    Data for Element  67
!
      data        isoki( 67)/    0         /
      data        nbfai( 67)/    9         /
      data      zcoreai( 67)/   13.000000d0/
      data   betaai(  1, 67)/ -100.000000d0/
      data   betaai(  2, 67)/ -100.000000d0/
      data      fgi(  1, 67)/ -100.000000d0/
      data      fgi(  3, 67)/ -100.000000d0/
      data      fgi(  8, 67)/    1.000000d0/
!
!                    Data for Element  68
!
      data        isoki( 68)/    0         /
      data        nbfai( 68)/    9         /
      data      zcoreai( 68)/   14.000000d0/
      data   betaai(  1, 68)/ -100.000000d0/
      data   betaai(  2, 68)/ -100.000000d0/
      data      fgi(  1, 68)/ -100.000000d0/
      data      fgi(  3, 68)/ -100.000000d0/
      data      fgi(  8, 68)/    1.000000d0/
!
!                    Data for Element  69
!
      data        isoki( 69)/    0         /
      data        nbfai( 69)/    9         /
      data      zcoreai( 69)/   15.000000d0/
      data   betaai(  1, 69)/ -100.000000d0/
      data   betaai(  2, 69)/ -100.000000d0/
      data      fgi(  1, 69)/ -100.000000d0/
      data      fgi(  3, 69)/ -100.000000d0/
      data      fgi(  8, 69)/    1.000000d0/
!
!                    Data for Element  70
!
      data        isoki( 70)/    0         /
      data        nbfai( 70)/    9         /
      data      zcoreai( 70)/   16.000000d0/
      data   betaai(  1, 70)/ -100.000000d0/
      data   betaai(  2, 70)/ -100.000000d0/
      data      fgi(  1, 70)/ -100.000000d0/
      data      fgi(  3, 70)/ -100.000000d0/
      data      fgi(  8, 70)/    1.000000d0/
!
!                    Data for Element  71
!
      data        isoki( 71)/    0         /
      data        nbfai( 71)/    9         /
      data      zcoreai( 71)/   17.000000d0/
      data   betaai(  1, 71)/ -100.000000d0/
      data   betaai(  2, 71)/ -100.000000d0/
      data      fgi(  1, 71)/ -100.000000d0/
      data      fgi(  3, 71)/ -100.000000d0/
      data      fgi(  8, 71)/    1.000000d0/
!
!                    Data for Element  72
!
      data        isoki( 72)/    0         /
      data        nbfai( 72)/    9         /
      data      zcoreai( 72)/    4.000000d0/
      data   betaai(  1, 72)/ -100.000000d0/
      data   betaai(  2, 72)/ -100.000000d0/
      data      fgi(  1, 72)/ -100.000000d0/
      data      fgi(  3, 72)/ -100.000000d0/
      data      fgi(  8, 72)/    1.000000d0/
!
!                    Data for Element  73
!
      data        isoki( 73)/    0         /
      data        nbfai( 73)/    9         /
      data      zcoreai( 73)/    5.000000d0/
      data   betaai(  1, 73)/ -100.000000d0/
      data   betaai(  2, 73)/ -100.000000d0/
      data      fgi(  1, 73)/ -100.000000d0/
      data      fgi(  3, 73)/ -100.000000d0/
      data      fgi(  8, 73)/    1.000000d0/
!
!                    Data for Element  74
!
      data        isoki( 74)/    0         /
      data        nbfai( 74)/    9         /
      data      zcoreai( 74)/    6.000000d0/
      data   betaai(  1, 74)/ -100.000000d0/
      data   betaai(  2, 74)/ -100.000000d0/
      data      fgi(  1, 74)/ -100.000000d0/
      data      fgi(  3, 74)/ -100.000000d0/
      data      fgi(  8, 74)/    1.000000d0/
!
!                    Data for Element  75
!
      data        isoki( 75)/    0         /
      data        nbfai( 75)/    9         /
      data      zcoreai( 75)/    7.000000d0/
      data   betaai(  1, 75)/ -100.000000d0/
      data   betaai(  2, 75)/ -100.000000d0/
      data      fgi(  1, 75)/ -100.000000d0/
      data      fgi(  3, 75)/ -100.000000d0/
      data      fgi(  8, 75)/    1.000000d0/
!
!                    Data for Element  76
!
      data        isoki( 76)/    0         /
      data        nbfai( 76)/    9         /
      data      zcoreai( 76)/    8.000000d0/
      data   betaai(  1, 76)/ -100.000000d0/
      data   betaai(  2, 76)/ -100.000000d0/
      data      fgi(  1, 76)/ -100.000000d0/
      data      fgi(  3, 76)/ -100.000000d0/
      data      fgi(  8, 76)/    1.000000d0/
!
!                    Data for Element  77
!
      data        isoki( 77)/    0         /
      data        nbfai( 77)/    9         /
      data      zcoreai( 77)/    9.000000d0/
      data   betaai(  1, 77)/ -100.000000d0/
      data   betaai(  2, 77)/ -100.000000d0/
      data      fgi(  1, 77)/ -100.000000d0/
      data      fgi(  3, 77)/ -100.000000d0/
      data      fgi(  8, 77)/    1.000000d0/
!
!                    Data for Element  78
!
      data        isoki( 78)/    0         /
      data        nbfai( 78)/    9         /
      data      zcoreai( 78)/   10.000000d0/
      data   betaai(  1, 78)/ -100.000000d0/
      data   betaai(  2, 78)/ -100.000000d0/
      data      fgi(  1, 78)/ -100.000000d0/
      data      fgi(  3, 78)/ -100.000000d0/
      data      fgi(  8, 78)/    1.000000d0/
!
!                    Data for Element  79
!
      data        isoki( 79)/    1         /
      data        nbfai( 79)/    9         /
      data      zcoreai( 79)/   11.000000d0/
      data        zetai( 79)/    1.986000d0/
      data   zetadi(  1, 79)/    3.288000d0/
      data  zetawti(  1, 79)/    1.000000d0/
      data   betaai(  1, 79)/   -4.000000d0/
      data   betaai(  2, 79)/   -4.000000d0/
      data   betaai(  3, 79)/  -50.000000d0/
      data      fgi(  1, 79)/  -10.120000d0/
      data      fgi(  2, 79)/   -9.226000d0/
      data      fgi(  3, 79)/   -5.432000d0/
      data      fgi(  4, 79)/   -4.358000d0/
      data      fgi(  5, 79)/  -14.692000d0/
      data      fgi(  6, 79)/  -12.017000d0/
      data      fgi(  9, 79)/    1.000000d0/
      data      fgi( 11, 79)/    7.364000d0/
      data      fgi( 12, 79)/    7.366000d0/
      data      fgi( 13, 79)/    9.149000d0/
      data      fgi( 14, 79)/    3.456875d0/
      data      fgi( 15, 79)/    2.749111d0/
      data      fgi( 16, 79)/    1.011275d0/
      data      fgi( 17, 79)/    1.280173d0/
      data      fgi( 18, 79)/    1.942687d0/
      data      fgi( 19, 79)/    0.827134d0/
      data      fgi( 20, 79)/    3.569660d0/
      data      fgi( 21, 79)/    3.530346d0/
!
!                    Data for Element  80
!
      data        isoki( 80)/    0         /
      data        nbfai( 80)/    4         /
      data      zcoreai( 80)/    2.000000d0/
      data   betaai(  1, 80)/ -100.000000d0/
      data   betaai(  2, 80)/ -100.000000d0/
      data      fgi(  1, 80)/ -100.000000d0/
      data      fgi(  3, 80)/ -100.000000d0/
      data      fgi(  8, 80)/    1.000000d0/

  end module Parameters_for_INDO_C
