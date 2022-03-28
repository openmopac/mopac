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

      module parameters_C
!
! This module holds all the parameters for the elements.  It is filled
! from the reference parameter sets
!
    logical, dimension(107) :: dorbs
    double precision, dimension(107) :: alp
    double precision, dimension(107,4) :: guess1, guess2, guess3
    double precision, dimension(60) :: v_par = 0.d0
    double precision, dimension(6,107) :: ddp
    double precision, dimension(9,107) :: po = 0.d0
    double precision, dimension(107) :: betas = 0.d0, betap = 0.d0, betad = 0.d0
    double precision, dimension(107) :: uss = 0.d0, upp = 0.d0, udd = 0.d0
    double precision, dimension(107) :: gpp = 0.d0, gp2 = 0.d0, hsp = 0.d0, gss = 0.d0, gsp = 0.d0
    double precision, dimension(107) :: am, ad, aq, dd, qq
    double precision, dimension(107) :: dsd = 0.d0, dpd = 0.d0, ddd = 0.d0
    double precision, dimension(107) :: zs, zp, zd, zsn = 0.d0, zpn = 0.d0, zdn = 0.d0
    double precision, dimension(107) :: eisol, eheat
    double precision, dimension(57:71) ::  eheat_sparkles
    double precision, dimension(107) :: ams
    double precision, dimension (107) :: tore, polvol, pocord, f0dd = 0.d0, f2dd = 0.d0, f4dd = 0.d0, &
    f0sd = 0.d0, g2sd = 0.d0, f0pd = 0.d0, f2pd = 0.d0, g1pd = 0.d0, g3pd = 0.d0, CPE_Zeta = 0.d0, &
      CPE_Z0 = 0.d0, CPE_B = 0.d0, CPE_Xlo = 0.d0, CPE_Xhi = 0.d0
    double precision, dimension (107) :: f0sd_store(107) ! Used by PARAM, not used in MOPAC
    double precision, dimension (107) :: g2sd_store(107) ! Used by PARAM, not used in MOPAC
    double precision, dimension (107) :: atom_radius_vdw, atom_radius_cosmo, &
    &  atom_radius_cosmo_oldcav
    double precision, dimension(100,100) :: xfac, alpb
    double precision :: dh2_a_parameters(6)
    integer, dimension (107,3) :: npq
    integer, dimension (107) :: iod, iop, ios, natorb, ndelec
    integer, parameter :: n_partyp = 72, n_partyp_alpb = 39, n_partyp_fn = 27
    logical, dimension (107) :: main_group
    double precision :: par1,  par2,  par3,  par4,  par5,  par6,  par7,  par8,  par9,  par10, &
                        par11, par12, par13, par14, par15, par16, par17, par18, par19, par20, &
                        par21, par22, par23, par24, par25, par26, par27, par28, par29, par30, &
                        par31, par32, par33, par34, par35, par36, par37, par38, par39, par40, &
                        par41, par42, par43, par44, par45, par46, par47, par48, par49, par50, &
                        par51, par52, par53, par54, par55, par56, par57, par58, par59, par60

    character :: partyp(n_partyp)*5
    character (len = 70) :: t_par(60), &
                        tpar1,  tpar2,  tpar3,  tpar4,  tpar5,  tpar6,  tpar7,  tpar8,  tpar9,  tpar10, &
                        tpar11, tpar12, tpar13, tpar14, tpar15, tpar16, tpar17, tpar18, tpar19, tpar20, &
                        tpar21, tpar22, tpar23, tpar24, tpar25, tpar26, tpar27, tpar28, tpar29, tpar30, &
                        tpar31, tpar32, tpar33, tpar34, tpar35, tpar36, tpar37, tpar38, tpar39, tpar40, &
                        tpar41, tpar42, tpar43, tpar44, tpar45, tpar46, tpar47, tpar48, tpar49, tpar50, &
                        tpar51, tpar52, tpar53, tpar54, tpar55, tpar56, tpar57, tpar58, tpar59, tpar60

    equivalence &
      (v_par(1),  par1),  (v_par(2), par2),   (v_par(3), par3),   (v_par(4), par4),   (v_par(5), par5),   &
      (v_par(6),  par6),  (v_par(7), par7),   (v_par(8), par8),   (v_par(9), par9),   (v_par(10), par10), &
      (v_par(11), par11), (v_par(12), par12), (v_par(13), par13), (v_par(14), par14), (v_par(15), par15), &
      (v_par(16), par16), (v_par(17), par17), (v_par(18), par18), (v_par(19), par19), (v_par(20), par20), &
      (v_par(21), par21), (v_par(22), par22), (v_par(23), par23), (v_par(24), par24), (v_par(25), par25), &
      (v_par(26), par26), (v_par(27), par27), (v_par(28), par28), (v_par(29), par29), (v_par(30), par30), &
      (v_par(31), par31), (v_par(32), par32), (v_par(33), par33), (v_par(34), par34), (v_par(35), par35), &
      (v_par(36), par36), (v_par(37), par37), (v_par(38), par38), (v_par(39), par39), (v_par(40), par40), &
      (v_par(41), par41), (v_par(42), par42), (v_par(43), par43), (v_par(44), par44), (v_par(45), par45), &
      (v_par(46), par46), (v_par(47), par47), (v_par(48), par48), (v_par(49), par49), (v_par(50), par50), &
      (v_par(51), par51), (v_par(52), par52), (v_par(53), par53), (v_par(54), par54), (v_par(55), par55), &
      (v_par(56), par56), (v_par(57), par57), (v_par(58), par58), (v_par(59), par59), (v_par(60), par60)
     equivalence &
      (t_par(1), tpar1),  (t_par(2),tpar2),   (t_par(3),tpar3),   (t_par(4),tpar4),   (t_par(5),tpar5),   &
      (t_par(6), tpar6),  (t_par(7),tpar7),   (t_par(8),tpar8),   (t_par(9),tpar9),   (t_par(10),tpar10), &
      (t_par(11),tpar11), (t_par(12),tpar12), (t_par(13),tpar13), (t_par(14),tpar14), (t_par(15),tpar15), &
      (t_par(16),tpar16), (t_par(17),tpar17), (t_par(18),tpar18), (t_par(19),tpar19), (t_par(20),tpar20), &
      (t_par(21),tpar21), (t_par(22),tpar22), (t_par(23),tpar23), (t_par(24),tpar24), (t_par(25),tpar25), &
      (t_par(26),tpar26), (t_par(27),tpar27), (t_par(28),tpar28), (t_par(29),tpar29), (t_par(30),tpar30), &
      (t_par(31),tpar31), (t_par(32),tpar32), (t_par(33),tpar33), (t_par(34),tpar34), (t_par(35),tpar35), &
      (t_par(36),tpar36), (t_par(37),tpar37), (t_par(38),tpar38), (t_par(39),tpar39), (t_par(40),tpar40), &
      (t_par(41),tpar41), (t_par(42),tpar42), (t_par(43),tpar43), (t_par(44),tpar44), (t_par(45),tpar45), &
      (t_par(46),tpar46), (t_par(47),tpar47), (t_par(48),tpar48), (t_par(49),tpar49), (t_par(50),tpar50), &
      (t_par(51),tpar51), (t_par(52),tpar52), (t_par(53),tpar53), (t_par(54),tpar54), (t_par(55),tpar55), &
      (t_par(56),tpar56), (t_par(57),tpar57), (t_par(58),tpar58), (t_par(59),tpar59), (t_par(60),tpar60)

    data ndelec / 20 * 0, &
!       Sc Ti V  Cr Mn Fe Co Ni Cu  Zn
      & 0, 0, 2, 2, 4, 4, 6, 8, 10, 10, 8 * 0, &   !  First transition series
      & 0, 0, 2, 2, 4, 4, 6, 8, 10, 10, 22* 0, &   !  Second transition series
      & 0, 0, 2, 2, 4, 4, 6, 8, 10, 10, 27 * 0 /   !  Third transition series
    double precision, dimension (n_partyp) :: defmax, defmin  ! Used by PARAM, not used in MOPAC

    save
!              H           Initial "s" Orbital Occupancies                     He
!              Li Be                                            B  C  N  O  F  Ne
!              Na Mg                                            Al Si P  S  Cl Ar
!              K  Ca Sc            Ti V  Cr Mn Fe Co Ni Cu Zn   Ga Ge As Se Br Kr
!              Rb Sr Y             Zr Nb Mo Tc Ru Rh Pd Ag Cd   In Sn Sb Te I  Xe
!              Cs Ba La Ce-Lu      Hf Ta W  Re Os Ir Pt Au Hg   Tl Pb Bi Po At Rn
!              Fr Ra Ac Th Pa U    Np Pu Am Cm Bk Cf            Cb ++ +  -- -  Tv
!                                      "s" shell
    data ios &
        &/ 1,                                                                2, &!    2
        &  1, 2,                                              2, 2, 2, 2, 2, 0, &!   10
        &  1, 2,                                              2, 2, 2, 2, 2, 0, &!   18
        &  1, 2, 2,              2, 2, 1, 2, 2, 2, 2, 1, 2,   2, 2, 2, 2, 2, 0, &!   36
        &  1, 2, 2,              2, 1, 1, 2, 1, 1, 0, 1, 2,   2, 2, 2, 2, 2, 0, &!   54
        &  1, 2, 2, 5*2,3*2,6*2, 2, 2, 1, 2, 2, 2, 1, 1, 2,   2, 2, 2, 2, 2, 0, &!   86
        &  1, 1, 2, 4, 2, 2,     2, 2, 2, 2, 2, 1, 0, 3,-3,   1, 2, 1,-2,-1, 0 /
!                                  /
!
!              H           Initial "p" Orbital Occupancies                   He
!              Li Be                                          B  C  N  O  F  Ne
!              Na Mg                                          Al Si P  S  Cl Ar
!              K  Ca Sc          Ti V  Cr Mn Fe Co Ni Cu Zn   Ga Ge As Se Br Kr
!              Rb Sr Y           Zr Nb Mo Tc Ru Rh Pd Ag Cd   In Sn Sb Te I  Xe
!              Cs Ba La Ce-Lu    Hf Ta W  Re Os Ir Pt Au Hg   Tl Pb Bi Po At Rn
!              Fr Ra Ac Th Pa U  Np Pu Am Cm Bk Cf (The rest are reserved for MOPAC)
!                                      "p" shell
    data iop / 0 ,                                                           0, &!    2
            &  0, 0,                                          1, 2, 3, 4, 5, 6, &!   10
            &  0, 0,                                          1, 2, 3, 4, 5, 6, &!   18
            &  0, 0, 0,          0, 0, 0, 0, 0, 0, 0, 0, 0,   1, 2, 3, 4, 5, 6, &!   36
            &  0, 0, 0,          0, 0, 0, 0, 0, 0, 0, 0, 0,   1, 2, 3, 4, 5, 6, &!   54
            &  0, 0, 0,  14*0,   0, 0, 0, 0, 0, 0, 0, 0, 0,   1, 2, 3, 4, 5, 6, &!   86
            &  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9*0                        /
!
!              H           Initial "d" Orbital Occupancies                   He
!              Li Be                                          B  C  N  O  F  Ne
!              Na Mg                                          Al Si P  S  Cl Ar
!              K  Ca Sc          Ti V  Cr Mn Fe Co Ni Cu Zn   Ga Ge As Se Br Kr
!              Rb Sr Y           Zr Nb Mo Tc Ru Rh Pd Ag Cd   In Sn Sb Te I  Xe
!              Cs Ba La Ce-Lu    Hf Ta W  Re Os Ir Pt Au Hg   Tl Pb Bi Po At Rn
!              Fr Ra Ac Th Pa U  Np Pu Am Cm Bk Cf (The rest are reserved for MOPAC)
!                                      "d" shell
    data iod / 0,                                                           0, &!    2
             & 0, 0,                                         0, 0, 0, 0, 0, 0, &!   10
             & 0, 0,                                         0, 0, 0, 0, 0, 0, &!   18
             & 0, 0, 1,          2, 3, 5, 5, 6, 7, 8, 10, 0, 0, 0, 0, 0, 0, 0, &!   36
             & 0, 0, 1,          2, 4, 5, 5, 7, 8,10, 10, 0, 0, 0, 0, 0, 0, 0, &!   54
             & 0, 0, 1,13*1,  1, 2, 3, 5, 5, 6, 7, 9, 10, 0, 0, 0, 0, 0, 0, 0, &!   86
             & 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9*0                     /
!
!                     Principal Quantum Numbers for all shells.
!
!              H                 "s"  shell                                  He
!              Li Be                                          B  C  N  O  F  Ne
!              Na Mg                                          Al Si P  S  Cl Ar
!              K  Ca Sc          Ti V  Cr Mn Fe Co Ni Cu Zn   Ga Ge As Se Br Kr
!              Rb Sr Y           Zr Nb Mo Tc Ru Rh Pd Ag Cd   In Sn Sb Te I  Xe
!              Cs Ba La Ce-Lu    Hf Ta W  Re Os Ir Pt Au Hg   Tl Pb Bi Po At Rn
!              Fr Ra Ac Th-Lr    ?? ?? ?? ??
!
data npq(1:107,1) / &
             & 1,                                                             1, &!  2
             & 2, 2,                                           2, 2, 2, 2, 2, 3, &! 10
             & 3, 3,                                           3, 3, 3, 3, 3, 4, &! 18
             & 4, 4,             4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, &! 36
             & 5, 5,             5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, &! 54
             & 6, 6, 14 * 6,     6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, &! 86
             & 11 * 0, 1, 0, 0, 0,3, 5 * 0 /
!
!              H                "p"  shell                                   He
!              Li Be                                          B  C  N  O  F  Ne
!              Na Mg                                          Al Si P  S  Cl Ar
!              K  Ca Sc          Ti V  Cr Mn Fe Co Ni Cu Zn   Ga Ge As Se Br Kr
!              Rb Sr Y           Zr Nb Mo Tc Ru Rh Pd Ag Cd   In Sn Sb Te I  Xe
!              Cs Ba La Ce-Lu    Hf Ta W  Re Os Ir Pt Au Hg   Tl Pb Bi Po At Rn
!              Fr Ra Ac Th-Lr    ?? ?? ?? ??
!
data npq(1:107,2) / &
             & 1,                                                             2, &!  2
             & 2, 2,                                           2, 2, 2, 2, 2, 2, &! 10
             & 3, 3,                                           3, 3, 3, 3, 3, 3, &! 18
             & 4, 4,             4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &! 36
             & 5, 5,             5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, &! 54
             & 6, 6, 14 * 6,     6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &! 86
             & 21 * 0 /
!
!              H                 "d"  shell                                  He
!              Li Be                                          B  C  N  O  F  Ne
!              Na Mg                                          Al Si P  S  Cl Ar
!              K  Ca Sc          Ti V  Cr Mn Fe Co Ni Cu Zn   Ga Ge As Se Br Kr
!              Rb Sr Y           Zr Nb Mo Tc Ru Rh Pd Ag Cd   In Sn Sb Te I  Xe
!              Cs Ba La Ce-Lu    Hf Ta W  Re Os Ir Pt Au Hg   Tl Pb Bi Po At Rn
!              Fr Ra Ac Th-Lr    ?? ?? ?? ??
!
data npq(1:107,3) / &
             & 0,                                                             0, &!  2
             & 0, 0,                                           0, 0, 0, 0, 0, 0, &! 10
             & 3, 3,                                           3, 3, 3, 3, 3, 4, &! 18
             & 3, 3,             3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, &! 36
             & 4, 4,             4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, &! 54
             & 5, 5, 14 * 5,     5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, &! 86
             & 21 * 0 /
!  Main_group is .true. if the Gss, Gpp, etc., are independent of zsn, zpn, and zdn,
!  .false. otherwise.
!
!  For main-group elements, the Gss value is important, for the transition metals,
!  Gss is less important.
!
    data main_group /&
   &  2*.true.,                              & ! H  - He
   &  8*.true.,                              & ! Li - Ne
   &  8*.true.,                              & ! Na - Ar
   &  2*.true., 9*.false., 7*.true.,         & ! K  - Kr
   &  2*.true., 9*.false., 7*.true.,         & ! Rb - Xe
   &  2*.true.,23*.false., 7*.true.,         & ! Cs - Rn
   & 21*.true.                               / ! Fr - Tv
!
!     ENTHALPIES OF FORMATION OF GASEOUS ATOMS ARE TAKEN FROM \ANNUAL
!     REPORTS,1974,71B,P 117\  THERE ARE SOME SIGNIFICANT DIFFERENCES
!     BETWEEN THE VALUES REPORTED THERE AND THE VALUES PREVIOUSLY IN
!     THE BLOCK DATA OF THIS PROGRAM.  ONLY THE THIRD  ROW ELEMENTS
!     HAVE BEEN UPDATED.
!
! ALL THE OTHER ELEMENTS ARE TAKEN FROM CRC HANDBOOK 1981-1982
    data eheat(1) / 52.102d0 /
    data eheat(2) / 0.000d0 /
    data eheat(3) / 38.410d0 /
    data eheat(4) / 76.960d0 /
    data eheat(5) / 135.700d0 /
    data eheat(6) / 170.890d0 /
    data eheat(7) / 113.000d0 /
    data eheat(8) / 59.559d0 /
    data eheat(9) / 18.890d0 /
    data eheat(10) / 0.000d0 /
    data eheat(11) / 25.650d0 /
    data eheat(12) / 35.000d0 /
    data eheat(13) / 79.490d0 /
    data eheat(14) / 108.390d0 /
    data eheat(15) / 75.570d0 /
    data eheat(16) / 66.400d0 /
    data eheat(17) / 28.990d0 /
    data eheat(18) / 0.000d0 /
    data eheat(19) / 21.420d0 /
    data eheat(20) / 42.600d0 /
    data eheat(21) / 90.300d0 /
    data eheat(22) / 112.300d0 /
    data eheat(23) / 122.900d0 /
    data eheat(24) / 95.000d0 /
    data eheat(25) / 67.700d0 /
    data eheat(26) / 99.300d0 /
    data eheat(27) / 102.400d0 /
    data eheat(28) / 102.800d0 /
    data eheat(29) / 80.700d0 /
    data eheat(30) / 31.170d0 /
    data eheat(31) / 65.400d0 /
    data eheat(32) / 89.500d0 /
    data eheat(33) / 72.300d0 /
    data eheat(34) / 54.300d0 /
    data eheat(35) / 26.740d0 /
    data eheat(36) / 0.000d0 /
    data eheat(37) / 19.600d0 /
    data eheat(38) / 39.100d0 /
    data eheat(39) / 101.500d0 /
    data eheat(40) / 145.500d0 /
    data eheat(41) / 172.400d0 /
    data eheat(42) / 157.300d0 /
    data eheat(43) / 162.000d0 /
    data eheat(44) / 155.500d0 /
    data eheat(45) / 133.000d0 /
    data eheat(46) / 90.000d0 /
    data eheat(47) / 68.100d0 /
    data eheat(48) / 26.720d0 /
    data eheat(49) / 58.000d0 /
    data eheat(50) / 72.200d0 /
    data eheat(51) / 63.200d0 /
    data eheat(52) / 47.000d0 /
    data eheat(53) / 25.517d0 /
    data eheat(54) / 0.000d0 /
    data eheat(55) / 18.700d0 /
    data eheat(56) / 42.500d0 /
    data eheat(57) /103.011d0/
    data eheat(58) /101.004d0/
    data eheat(59) /84.990d0/
    data eheat(60) /78.298d0/
    data eheat(61) /83.174d0/
    data eheat(62) /49.402d0/
    data eheat(63) /41.898d0/
    data eheat(64) /95.007d0/
    data eheat(65) /92.902d0/
    data eheat(66) /69.407d0/
    data eheat(67) /71.893d0/
    data eheat(68) /75.791d0/
    data eheat(69) /55.500d0/
    data eheat(70) /36.358d0/
    data eheat(71) /102.199d0/
    data eheat(72) / 148.000d0 /
    data eheat(73) / 186.900d0 /
    data eheat(74) / 203.100d0 /
    data eheat(75) / 185.000d0 /
    data eheat(76) / 188.000d0 /
    data eheat(77) / 160.000d0 /
    data eheat(78) / 135.200d0 /
    data eheat(79) / 88.000d0 /
    data eheat(80) / 14.690d0 /
    data eheat(81) / 43.550d0 /
    data eheat(82) / 46.620d0 /
    data eheat(83) / 50.100d0 /
    data eheat(86) / 0.000d0 /
    data eheat(90) / 1674.64d0/
    data eheat(102) / 207.0d0 /
    data eheat(103) / 0.0d0 /
    data eheat(104) / 0.0d0 /
    data eheat(105) / 0.0d0 /
    data eheat(106) / 0.0d0 /
    data eheat(107) / 0.0d0 /
    data eheat_sparkles(57) / 928.9D0/   !  Represents La(+++)
    data eheat_sparkles(58) / 944.7D0 /  !  Represents Ce(+++)
    data eheat_sparkles(59) / 952.9D0 /  !  Represents Pr(+++)
    data eheat_sparkles(60) / 962.8D0 /  !  Represents Nd(+++)
    data eheat_sparkles(61) / 976.9D0 /  !  Represents Pm(+++)
    data eheat_sparkles(62) / 974.4d0 /  !  Represents Sm(+++)
    data eheat_sparkles(63) / 1006.6d0/  !  Represents Eu(+++)
    data eheat_sparkles(64) / 991.37d0/  !  Represents Gd(+++)
    data eheat_sparkles(65) / 999.0d0/   !  Represents Tb(+++)
    data eheat_sparkles(66) / 1001.3d0 / !  Represents Dy(+++)
    data eheat_sparkles(67) / 1009.6d0 / !  Represents Ho(+++)
    data eheat_sparkles(68) / 1016.15d0 /!  Represents Er(+++)
    data eheat_sparkles(69) / 1022.06d0/ !  Represents Tm(+++)
    data eheat_sparkles(70) / 1039.03d0/ !  Represents Ho(+++)
    data eheat_sparkles(71) / 1031.2d0 / !  Represents Lu(+++)
    data eisol / 107 * 0.0d0 /
    data ams / 1.00790d0, 4.00260d0, 6.94000d0, 9.01218d0, 10.81000d0, &
   & 12.01100d0, 14.00670d0, 15.99940d0, 18.99840d0, 20.17900d0, 22.98977d0, &
   & 24.30500d0, 26.98154d0, 28.08550d0, 30.97376d0, 32.06000d0, 35.45300d0, &
   & 39.94800d0, 39.09830d0, 40.08000d0, 44.95590d0, 47.90000d0, 50.94150d0, &
   & 51.99600d0, 54.93800d0, 55.84700d0, 58.93320d0, 58.71000d0, 63.54600d0, &
   & 65.38000d0, 69.73500d0, 72.59000d0, 74.92160d0, 78.96000d0, 79.90400d0, &
   & 83.80000d0, 85.46780d0, 87.62000d0, 88.90590d0, 91.22000d0, 92.90640d0, &
   & 95.94000d0, 98.90620d0, 101.0700d0, 102.9055d0, 106.4000d0, 107.8680d0, &
   & 112.4100d0, 114.8200d0, 118.6900d0, 121.7500d0, 127.6000d0, 126.9045d0, &
   & 131.3000d0, 132.9054d0, 137.3300d0, 138.9060d0, 140.1160d0, 140.9077d0, &
   & 144.2400d0, 145.0000d0, 150.3600d0, 151.9640d0, 157.2500d0, 158.9253d0, &
   & 162.5000d0, 164.9303d0, 167.2600d0, 168.9342d0, 173.0400d0, 174.9670d0, &
   & 178.4900d0, 180.9479d0, &
   & 183.8500d0, 186.2070d0, 190.2000d0, 192.2200d0, 195.0900d0, 196.9665d0, &
   & 200.5900d0, 204.3700d0, 207.2000d0, 208.9804d0, 209.0000d0, 210.0000d0, &
   & 222.0000d0, 223.0000d0, 226.0000d0, 227.0000d0, 232.0381d0, 231.0359d0, &
   & 238.0289d0, 5 * 0.000d0, 0.00005d0, 3*0.d0, 1.0079d0, &
   & 5 * 0.000d0 /
    data partyp (1),  defmin (1),  defmax (1) / "USS  ", - 200.d0, 60.d0 /
    data partyp (2),  defmin (2),  defmax (2) / "UPP  ", - 200.d0, 60.d0 /
    data partyp (3),  defmin (3),  defmax (3) / "UDD  ", - 200.d0, 60.d0 /
    data partyp (4),  defmin (4),  defmax (4) / "ZS   ", 0.6d0, 6.d0 /
    data partyp (5),  defmin (5),  defmax (5) / "ZP   ", 0.6d0, 6.d0 /
    data partyp (6),  defmin (6),  defmax (6) / "ZD   ", 0.6d0, 6.d0 /
    data partyp (7),  defmin (7),  defmax (7) / "BETAS", - 70.d0, 10.d0 /
    data partyp (8),  defmin (8),  defmax (8) / "BETAP", - 70.d0, 10.d0 /
    data partyp (9),  defmin (9),  defmax (9) / "BETAD", - 70.d0, 10.d0 /
    data partyp (10), defmin (10), defmax (10) / "GSS  ", 1.d0, 20.d0 /
    data partyp (11), defmin (11), defmax (11) / "GSP  ", 1.d0, 20.d0 /
    data partyp (12), defmin (12), defmax (12) / "GPP  ", 1.d0, 20.d0 /
    data partyp (13), defmin (13), defmax (13) / "GP2  ", 1.0d0, 19.d0 /
    data partyp (14), defmin (14), defmax (14) / "HSP  ", 1.d0, 5.d0 /
    data partyp (15), defmin (15), defmax (15) / "F0SD ", 1.d0,10.d0 /
    data partyp (16), defmin (16), defmax (16) / "G2SD ", 1.d0,10.d0 /
    data partyp (17), defmin (17), defmax (17) / "POC  ", 1.d0,10.d0 /
    data partyp (18), defmin (18), defmax (18) / "ALP  ", 0.5d0, 6.d0 /
    data partyp (19), defmin (19), defmax (19) / "ZSN  ", 0.5d0,25.d0 /
    data partyp (20), defmin (20), defmax (20) / "ZPN  ", 0.5d0,25.d0 /
    data partyp (21), defmin (21), defmax (21) / "ZDN  ", 0.5d0,25.d0 /
    data partyp (22), defmin (22), defmax (22) / "C_ZET", -0.5d0,25.d0 /
    data partyp (23), defmin (23), defmax (23) / "C_Z0 ", -0.5d0,25.d0 /
    data partyp (24), defmin (24), defmax (24) / "C_B  ", -0.5d0,25.d0 /
    data partyp (25), defmin (25), defmax (25) / "C_XLO", -0.5d0,25.d0 /
    data partyp (26), defmin (26), defmax (26) / "C_XHI", -0.5d0,25.d0 /
    data partyp (27), defmin (27), defmax (27) / "FN11 ", -1.d0, 1.d0 /
    data partyp (28), defmin (28), defmax (28) / "FN21 ", 0.5d0, 3.d0 /
    data partyp (29), defmin (29), defmax (29) / "FN31 ", 0.5d0, 3.d0 /
    data partyp (30), defmin (30), defmax (30) / "FN12 ", - 1.d0, 1.d0 /
    data partyp (31), defmin (31), defmax (31) / "FN22 ", 0.15d0, 3.d0 /
    data partyp (32), defmin (32), defmax (32) / "FN32 ", 0.5d0, 3.5d0 /
    data partyp (33), defmin (33), defmax (33) / "FN13 ", - 1.d0, 1.d0 /
    data partyp (34), defmin (34), defmax (34) / "FN23 ", 1.d0, 3.d0 /
    data partyp (35), defmin (35), defmax (35) / "FN33 ", 0.5d0, 4.d0 /
    data partyp (36), defmin (36), defmax (36) / "FN14 ", - 1.d0, 1.d0 /
    data partyp (37), defmin (37), defmax (37) / "FN24 ", 1.d0, 3.d0 /
    data partyp (38), defmin (38), defmax (38) / "FN34 ", 0.5d0, 3.d0 /
    data partyp (39), defmin (39), defmax (39) / "ALPB_", 0.9d0, 3.0d0 /
    data partyp (40), defmin (40), defmax (40) / "XFAC_", 0.5d0, 30.d0 /
    data partyp (41), defmin (41), defmax (41) / "PAR  ", 0.5d0, 6.d0 /
    data partyp (42), defmin (42), defmax (42) / "NORBS", 0.9d0, 9.1d0/
    data partyp (43)                           / "ZCORE"/
    data partyp (44)                           / "ZSP  "/
    data partyp (45)                           / "ZD1  "/
    data partyp (46)                           / "ZD2  "/
    data partyp (47)                           / "ZWT1 "/
    data partyp (48)                           / "ZWT2 "/
    data partyp (49)                           / "FG1  "/
    data partyp (50)                           / "FG2  "/
    data partyp (51)                           / "FG3  "/
    data partyp (52)                           / "FG4  "/
    data partyp (53)                           / "FG5  "/
    data partyp (54)                           / "FG6  "/
    data partyp (55)                           / "FG7  "/
    data partyp (56)                           / "FG8  "/
    data partyp (57)                           / "FG9  "/
    data partyp (58)                           / "FG10 "/
    data partyp (59)                           / "FG11 "/
    data partyp (60)                           / "FG12 "/
    data partyp (61)                           / "FG13 "/
    data partyp (62)                           / "FG14 "/
    data partyp (63)                           / "FG15 "/
    data partyp (64)                           / "FG16 "/
    data partyp (65)                           / "FG17 "/
    data partyp (66)                           / "FG18 "/
    data partyp (67)                           / "FG19 "/
    data partyp (68)                           / "FG20 "/
    data partyp (69)                           / "FG21 "/
    data partyp (70)                           / "FG22 "/
    data partyp (71)                           / "FG23 "/
    data partyp (72)                           / "FG24 "/
      end module parameters_C
