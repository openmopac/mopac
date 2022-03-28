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

module symmetry_C
    integer, parameter :: ntbs = 38

    save
  integer :: &
  &  igroup, & !  Index of point-group
  &  nclass, & !  Number of classes used in the point-group
  &  nirred, & !  Number of irreducible representations
  &  nsym,   &
  &  nent,   &
  &  dummy

  integer, dimension(20)               :: ielem
  integer, dimension(6)                :: jy
  integer, dimension (348)             :: nallop
  integer, dimension (ntbs)            :: ntab
  integer, dimension (764)             :: nallg

  character, dimension(20)             :: jx*4      !
  character                            :: name*4
  character, dimension (406)           :: allrep*4
  character, dimension(:), allocatable :: namo*4             !  Names of irreducible representations
  character                            :: state_spin*8       !  Spin-state, e.g., "TRIPLET"
  character                            :: state_Irred_Rep*4  !  Irreducible representation, e.g. "T2g"
  integer                              :: state_QN           !  Quantum number, e.g., 2, 3, 4.

  integer, dimension(:), allocatable :: &
  &  jndex,  & !  Principal quantum number of irreducible representation
  &  dummys


  integer, dimension(:,:), allocatable :: &
  &  ipo,    & !
  &  jelem
  double precision, dimension(3,3)              :: cub
  double precision, dimension(9,120)            :: r
  double precision, dimension(3,3,20)           :: elem
  double precision, dimension(20,5)             :: group
  double precision, dimension(:), allocatable   :: depmul
  integer, dimension(:), allocatable        :: locpar, idepfn, locdep



    !*********************************************************************
   !
   !  BLOCK DATA FOR ALL THE POINT-GROUPS USED IN MOPAC
   !
   !  Point-Groups are stored in order of increasing symmetry.  Groups
   !   supported are, in order
   !
   ! C1, Cs, Ci, C2, D2, C2v, C2h, D2h, C3, C4, S4, D3, C3v, C3h, C5, C6,
   ! S6, C7, C8, S8, C4v, D4, D2d, C5v, C6v, D6, C4h, D3h, D3d, C7v, C7h,
   ! C8v, D8, D4d, D4h, C5h, D5, D5h, D5d, C6h, D6h, D6d, D7, D7h, D7d,
   ! C8h, D8h, T, Td, O, Th, Oh, I, Ih, Cv, Dh, R3
   !
   !  Character tables are stored separately from the point groups, so
   !  that more than one group can use the same character table.
   !  The totally symmetric representation is assumed to exist, and
   !  has been suppressed from the character tables.
   !
   !  Each point group is represented by two arrays, a I<group symbol>
   !  and a R<group symbol> array.  The structure of these tables is as
   !  follows:
   !
   !      I<group symbol>
   !
   ! (1):    Number of classes of operation used.
   ! (2):    Number of Irreducible Representations.
   ! (3):    The number of the associated character table.
   ! (4):    A "Magic number" identifying the group.
   ! (5 on): Numbers indicating the operations of the group
   !
   !      R<group symbol>
   !
   ! (1):    The name of the group
   ! (2 on): The names of the irreducible representations.
   !
   !  Character tables are named arbitarily.  Each table has an associated
   !  entry in the NTAB array.  NTAB(i) holds the number of elements in
   !  table i.
   !
   !*********************************************************************
   !   Storage of point-groups:  Point-Group Name,Number of Classes,
   !                             Names of Classes,Number of Irred. Reps.,
   !                             Names of Irred. Reps,Number of data in
   !                             point-group,Name of Representation
   !.. Implicit Declarations ..
 !
    ! IG11
    data ntab(1) / 1 /
    data nallg(1:1) / 1 /
    ! IG 21
    data ntab(2) / 2 /
    data nallg(2:3) / 1, -1 /
    ! IG42
    data ntab(3) / 9 /
    data nallg(4:12) / 1, 1, -1, 1, -1, 1, 1, -1, -1 /
    ! IG83
    data ntab(4) / 28 /
    data nallg(13:40) / 1, 1, -1, 1, 1, -1, 1, 1, 1, -1, -1, 1, 1, 1, 1, &
         & - 1, 1, 1, -1, -1, 1, -1, 1, -1, 1, -1, -1, -1 /
    ! IG31
    data ntab(5) / 2 /
    data nallg(41:42) / 2, -1 /
    ! IG41
    data ntab(6) / 4 /
    data nallg(43:46) / 1, -1, 2, 0 /
    ! IG51
    data ntab(7) / 4 /
    data nallg(47:50) / 2, 51, 2, 52 /
    ! IG61
    data ntab(8) / 6 /
    data nallg(51:56) / 1, -1, 2, 1, 2, -1 /
    ! IG71
    data ntab(9) / 6 /
    data nallg(57:62) / 2, 71, 2, 72, 2, 73 /
    ! IG81
    data ntab(10) / 8 /
    data nallg(63:70) / 1, -1, 2, 81, 2, 0, 2, 83 /
    ! IG84
    data ntab(11) / 12 /
    data nallg(71:82) / 1, 1, -1, 1, -1, 1, 1, -1, -1, 2, 0, 0 /
    ! IG52
    data ntab(12) / 9 /
    data nallg(83:91) / 1, 1, -1, 2, 51, 0, 2, 52, 0 /
    ! IG62
    data ntab(13) / 20 /
    data nallg(92:111) / 1, 1, 1, -1, 1, 1, -1, 1, 1, 1, -1, -1, 2, -1, &
         & 0, -2, 2, -1, 0, 2 /
    ! IG72
    data ntab(14) / 12 /
    data nallg(112:123) / 1, 1, -1, 2, 71, 0, 2, 72, 0, 2, 73, 0 /
    ! IG82
    data ntab(15) / 18 /
    data nallg(124:141) / 1, 1, -1, 1, -1, -1, 1, -1, 1, 2, 81, 0, 2, 0, &
         & 0, 2, 83, 0 /
    ! IG4H
    data ntab(16) / 36 /
    data nallg(142:177) / 1, 1, 1, -1, 1, 1, -1, 1, 1, 1, -1, -1, 1, -1, &
         & 1, 1, 1, -1, 1, -1, 1, -1, -1, 1, 1, -1, -1, -1, 2, 0, 0, -2, 2, &
         & 0, 0, 2 /
    ! IGC5H
    data ntab(17) / 15 /
    data nallg(178:192) / 1, 1, -1, 2, 51, 2, 2, 51, -2, 2, 52, 2, 2, 52, -2 /
    ! IGD5G
    data ntab(18) / 28 /
    data nallg(193:220) / 1, 1, -1, -1, 1, 1, 1, -1, 1, 1, -1, 1, 2, 51, &
         & 2, 0, 2, 51, -2, 0, 2, 52, 2, 0, 2, 52, -2, 0 /
    ! IGC6H
    data ntab(19) / 21 /
    data nallg(221:241) / 1, 1, -1, 1, -1, 1, 1, -1, -1, 2, 1, 2, 2, 1, &
         & - 2, 2, -1, 2, 2, -1, -2 /
    ! IGD6H
    data ntab(20) / 44 /
    data nallg(242:285) / 1, 1, 1, -1, 1, 1, -1, 1, 1, 1, -1, -1, 1, -1, &
         & 1, 1, 1, -1, 1, -1, 1, -1, -1, 1, 1, -1, -1, -1, 2, 1, 0, 2, 2, 1, &
         & 0, -2, 2, -1, 0, 2, 2, -1, 0, -2 /
    ! IGD6D
    data ntab(21) / 24 /
    data nallg(286:309) / 1, 1, -1, 1, -1, -1, 1, -1, 1, 2, 121, 0, 2, 1, &
         & 0, 2, 0, 0, 2, -1, 0, 2, 125, 0 /
    ! IGC7H
    data ntab(22) / 21 /
    data nallg(310:330) / 1, 1, -1, 2, 71, 2, 2, 71, -2, 2, 72, 2, 2, 72, &
         & - 2, 2, 73, 2, 2, 73, -2 /
    ! IGD7H
    data ntab(23) / 36 /
    data nallg(331:366) / 1, 1, -1, -1, 1, 1, 1, -1, 1, 1, -1, 1, 2, 71, &
         & 2, 0, 2, 71, -2, 0, 2, 72, 2, 0, 2, 72, -2, 0, 2, 73, 2, 0, 2, &
         & 73, -2, 0 /
    ! IGD8H
    data ntab(24) / 52 /
    data nallg(367:418) / 1, 1, -1, -1, 1, 1, 1, -1, 1, 1, -1, 1, 1, -1, &
         & - 1, -1, 1, -1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 2, 81, -2, 0, 2, &
         & 81, 2, 0, 2, 0, 2, 0, 2, 0, -2, 0, 2, 83, -2, 0, 2, 83, 2, 0 /
    ! IGC8H
    data ntab(25) / 27 /
    data nallg(419:445) / 1, 1, -1, 1, -1, 1, 1, -1, -1, 2, 81, 2, 2, 81, &
         & - 2, 2, 0, 2, 2, 0, -2, 2, 83, 2, 2, 83, -2 /
    ! IGT
    data ntab(26) / 6 /
    data nallg(446:451) / 2, -1, 2, 3, 0, -1 /
    ! IGTH
    data ntab(27) / 20 /
    data nallg(452:471) / 1, 1, -1, 1, 2, -1, 2, 2, 2, -1, -2, 2, 3, 0, &
         & 3, -1, 3, 0, -3, -1 /
    ! IGTD
    data ntab(28) / 16 /
    data nallg(472:487) / 1, -1, 1, 1, 2, 0, -1, 2, 3, -1, 0, -1, 3, 1, &
         & 0, -1 /
    ! IGOH
    data ntab(29) / 36 /
    data nallg(488:523) / 1, 1, -1, 1, 1, -1, 1, 1, 1, -1, -1, 1, 2, 0, &
         & 2, -1, 2, 0, -2, -1, 3, -1, 3, 0, 3, -1, -3, 0, 3, 1, 3, 0, 3, 1, -3, 0 /
    ! IGI
    data ntab(30) / 12 /
    data nallg(524:535) / 3, -1, 101, 3, -1, 103, 4, 0, -1, 5, 1, 0 /
    ! IGIH
    data ntab(31) / 36 /
    data nallg(536:571) / 1, 1, -1, 1, 3, 101, 3, 0, 3, 101, -3, 0, 3, 103, &
         & 3, 0, 3, 103, -3, 0, 4, -1, 4, 1, 4, -1, -4, 1, 5, 0, 5, -1, 5, 0, -5, -1 /
    ! IGCV
    data ntab(32) / 20 /
    data nallg(572:591) / 1, 1, 1, -1, 2, 71, 61, 0, 2, 72, 62, 0, 2, 73, &
         & 63, 0, 2, 74, 64, 0 /
    ! IGDH
    data ntab(33) / 55 /
    data nallg(592:646) / 1, 1, 1, -1, 1, 1, 1, 1, 1, -1, 1, 1, 1, -1, -1, &
         & 2, 71, 61, 2, 0, 2, 71, 61, -2, 0, 2, 72, 62, 2, 0, 2, 72, 62, -2, &
         & 0, 2, 73, 63, 2, 0, 2, 73, 63, -2, 0, 2, 74, 64, 2, 0, 2, 74, 64, -2, 0 /
    ! IGC4H
    data ntab(34) / 15 /
    data nallg(647:661) / 1, 1, -1, 1, -1, 1, 1, -1, -1, 2, 0, 2, 2, 0, -2 /
    ! IGD3
    data ntab(35) / 6 /
    data nallg(662:667) / 1, 1, -1, 2, -1, 0 /
    ! IGC3H
    data ntab(36) / 9 /
    data nallg(668:676) / 2, -1, 2, 1, 1, -1, 2, -1, -2 /
    ! IGO
    data ntab(37) / 12 /
    data nallg(677:688) / 1, -1, 1, 2, 0, -1, 3, -1, 0, 3, 1, 0 /
    ! IGR3
    data ntab(38) / 60 /
    data nallg(689:764) / 1, 1, 1, -1, 3, 101, 0, 3, 3, 101, 0, -3, 5, 0, -1, &
         & 5, 5, 0, -1, -5, 7, 106, 1, 7, 7, 106, 1, -7, 9, -1, 0, 9, 9, -1, &
         & 0, -9, 11, 1, -1, 11, 11, 1, -1, -11, 13, 101, 1, 13, 13, 101, 1, -13, &
         & 15, 0, 0, 15, 15, 0, 0, -15, 17, 106, -1, 17, 17, 106, -1, -17, &
         & 19, -1, 1, 19, 19, -1, 1, -19 /
    !
    ! C1
    data allrep (1:2) / "C1", "A" /
    data nallop(1:4) / 0, 1, 1, 0 /
    ! CS
    data allrep (3:5) / "Cs", "A'", "A""" /
    data nallop(5:9) / 1, 2, 2, 8, 4 /
    ! CI
    data allrep (6:8) / "Ci", "Ag", "Au" /
    data nallop(10:14) / 1, 2, 2, 64, 7 /
    ! C2
    data allrep (9:11) / "C2", "A", "B" /
    data nallop(15:19) / 1, 2, 2, 4, 3 /
    ! D2
    data allrep (12:16) / "D2", "A", "B1", "B2", "B3" /
    data nallop(20:25) / 2, 4, 3, 7, 3, 2 /
    ! C2V
    data allrep (17:21) / "C2v", "A1", "A2", "B1", "B2" /
    data nallop(26:31) / 2, 4, 3, 52, 3, 5 /
    ! C2H
    data allrep (22:26) / "C2h", "Ag", "Au", "Bg", "Bu" /
    data nallop(32:37) / 2, 4, 3, 76, 3, 7 /
    ! D2H
    data allrep (27:35) / "D2h", "Ag", "B1g", "B2g", "B3g", "Au", "B1u", &
         & "B2u", "B3u" /
    data nallop(38:44) / 3, 8, 4, 127, 3, 2, 7 /
    ! C3
    data allrep (36:38) / "C3", "A", "E" /
    data nallop(45:49) / 1, 2, 5, 128, 8 /
    ! C4
    data allrep (39:42) / "C4", "A", "B", "E" /
    data nallop(50:54) / 1, 3, 6, 260, 9 /
    ! S4
    data allrep (43:46) / "S4", "A", "B", "E" /
    data nallop(55:59) / 1, 3, 6, 8196, 14 /
    ! D3
    data allrep (47:50) / "D3", "A1", "A2", "E" /
    data nallop(60:65) / 2, 3, 35, 129, 8, 1 /
    ! C3V
    data allrep (51:54) / "C3v", "A1", "A2", "E" /
    data nallop(66:71) / 2, 3, 35, 144, 8, 5 /
    ! C3H
    data allrep (55:59) / "C3h", "A'", "E'", "A""", "E""" /
    data nallop(72:77) / 2, 4, 36, 136, 8, 4 /
    ! C5
    data allrep (60:63) / "C5", "A", "E1", "E2" /
    data nallop(78:82) / 1, 3, 7, 512, 10 /
    ! C6
    data allrep (64:68) / "C6", "A", "B", "E1", "E2" /
    data nallop(83:87) / 1, 4, 8, 1156, 11 /
    ! S6
    data allrep (69:73) / "S6", "Ag", "Au", "Eu", "Eg" /
    data nallop(88:92) / 1, 4, 8, 16576, 15 /
    ! C7
    data allrep (74:78) / "C7", "A", "E1", "E2", "E3" /
    data nallop(93:97) / 1, 4, 9, 2048, 12 /
    ! C8
    data allrep (79:84) / "C8", "A", "B", "E1", "E2", "E3" /
    data nallop(98:102) / 1, 5, 10, 4356, 13 /
    ! S8
    data allrep (85:90) / "S8", "A", "B", "E1", "E2", "E3" /
    data nallop(103:107) / 1, 5, 10, 33028, 16 /
    ! C4V
    data allrep (91:96) / "C4v", "A1", "A2", "B1", "B2", "E" /
    data nallop(108:113) / 2, 5, 11, 308, 9, 5 /
    ! D4
    data allrep (97:102) / "D4", "A1", "A2", "B1", "B2", "E" /
    data nallop(114:119) / 2, 5, 11, 263, 9, 1 /
    ! D2D
    data allrep (103:108) / "D2d", "A1", "A2", "B2", "B1", "E" /
    data nallop(120:125) / 2, 5, 11, 8244, 14, 5 /
    ! C5V
    data allrep (109:113) / "C5v", "A1", "A2", "E1", "E2" /
    data nallop(126:131) / 2, 4, 12, 528, 10, 5 /
    ! C6V
    data allrep (114:120) / "C6v", "A1", "B1", "A2", "B2", "E1", "E2" /
    data nallop(132:138) / 3, 6, 13, 1204, 8, 5, 3 /
    ! D6
    data allrep (121:127) / "D6", "A1", "B1", "A2", "B2", "E1", "E2" /
    data nallop(139:145) / 3, 6, 13, 1159, 8, 1, 3 /
    ! C4H
    data allrep (128:134) / "C4h", "Ag", "Au", "Bg", "Bu", "Eg", "Eu" /
    data nallop(146:151) / 2, 6, 34, 8524, 9, 7 /
    ! D3H
    data allrep (135:141) / "D3h", "A1'", "A1""", "A2'", "A2""", "E""", "E'" /
    data nallop(152:158) / 3, 6, 13, 153, 8, 1, 4 /
    ! D3D
    data allrep (142:148) / "D3d", "A1g", "A1u", "A2g", "A2u", "Eu", "Eg" /
    data nallop(159:165) / 3, 6, 13, 16594, 8, 2, 7 /
    ! C7V
    data allrep (149:154) / "C7v", "A1", "A2", "E1", "E2", "E3" /
    data nallop(166:171) / 2, 5, 14, 2064, 12, 5 /
    ! C7H
    data allrep (155:163) / "C7h", "A'", "A""", "E1'", "E1""", "E2'", "E2""", &
         & "E3'", "E3""" /
    data nallop(172:177) / 2, 8, 22, 2056, 12, 4 /
    ! C8V
    data allrep (164:171) / "C8v", "A1", "A2", "B2", "B1", "E1", "E2", "E3" /
    data nallop(178:183) / 2, 7, 15, 4404, 13, 5 /
    ! D8
    data allrep (172:179) / "D8", "A1", "A2", "B2", "B1", "E1", "E2", "E3" /
    data nallop(184:189) / 2, 7, 15, 4359, 13, 1 /
    ! D4D
    data allrep (180:187) / "D4d", "A1", "A2", "B1", "B2", "E1", "E2", "E3" /
    data nallop(190:195) / 2, 7, 15, 33076, 16, 5 /
    ! D4H
    data allrep (188:198) / "D4h", "A1g", "A1u", "A2g", "A2u", "B1g", "B1u", &
         & "B2g", "B2u", "Eg", "Eu" /
    data nallop(196:202) / 3, 10, 16, 8575, 9, 1, 4 /
    ! C5H
    data allrep (199:205) / "C5h", "A'", "A""", "E1'", "E1""", "E2'", "E2""" /
    data nallop(203:208) / 2, 6, 17, 520, 10, 4 /
    ! D5
    data allrep (206:210) / "D5", "A1", "A2", "E1", "E2" /
    data nallop(209:214) / 2, 4, 12, 513, 10, 1 /
    ! D5H
    data allrep (211:219) / "D5h", "A1'", "A1""", "A2'", "A2""", "E1'", "E1""",&
         &  "E2'", "E2""" /
    data nallop(215:221) / 3, 8, 18, 537, 10, 4, 5 /
    ! D5D
    data allrep (220:228) / "D5d", "A1g", "A1u", "A2g", "A2u", "E1g", "E1u", &
         & "E2g", "E2u" /
    data nallop(222:228) / 3, 8, 18, 66130, 10, 7, 5 /
    ! C6H
    data allrep (229:237) / "C6h", "Ag", "Au", "Bg", "Bu", "E1g", "E1u", "E2g",&
         &  "E2u" /
    data nallop(229:234) / 2, 8, 19, 17612, 11, 7 /
    ! D6H
    data allrep (238:250) / "D6h", "A1g", "A1u", "A2g", "A2u", "B1g", "B1u", &
         & "B2g", "B2u", "E1g", "E1u", "E2g", "E2u" /
    data nallop(235:241) / 3, 12, 20, 17663, 11, 1, 7 /
    ! D6D
    data allrep (251:260) / "D6d", "A1", "A2", "B1", "B2", "E1", "E2", "E3", &
         & "E4", "E5" /
    data nallop(242:247) / 2, 9, 21, 140468, 18, 5 /
    ! D7
    data allrep (261:266) / "D7", "A1", "A2", "E1", "E2", "E3" /
    data nallop(248:253) / 2, 5, 14, 2049, 12, 1 /
    ! D7H
    data allrep (267:277) / "D7h", "A1'", "A1""", "A2'", "A2""", "E1'", &
         & "E1""", "E2'", "E2""", "E3'", "E3""" /
    data nallop(254:260) / 3, 10, 23, 2073, 12, 4, 5 /
    ! D7D
    data allrep (278:288) / "D7d", "A1g", "A1u", "A2g", "A2u", "E1g", "E1u", &
         & "E2g", "E2u", "E3g", "E3u" /
    data nallop(261:267) / 3, 10, 23, 2130, 12, 7, 5 /
    ! C8H
    data allrep (289:299) / "C8h", "Ag", "Au", "Bg", "Bu", "E1g", "E1u", "E2g",&
         &  "E2u", "E3g", "E3u" /
    data nallop(268:273) / 2, 10, 25, 45388, 13, 7 /
    ! D8H
    data allrep (300:314) / "D8h", "A1g", "A1u", "A2g", "A2u", "B2u", "B2g", &
         & "B1u", "B1g", "E1g", "E1u", "E2g", "E2u", "E3g", "E3u" /
    data nallop(274:280) / 3, 14, 24, 45439, 13, 4, 5 /
    ! T
    data allrep (315:318) / "T", "A", "E", "T" /
    data nallop(281:286) / 2, 3, 26, 262276, 8, 3 /
    ! TD
    data allrep (319:324) / "Td", "A1", "A2", "E", "T1", "T2" /
    data nallop(287:293) / 3, 5, 28, 270516, 5, 8, 3 /
    ! O
    data allrep (325:330) / "O", "A1", "A2", "E", "T2", "T1" /
    data nallop(294:299) / 2, 5, 37, 262273, 1, 8 /
    ! TH
    data allrep (331:337) / "Th", "Ag", "Au", "Eg", "Eu", "Tg", "Tu" /
    data nallop(300:306) / 3, 6, 27, 278732, 8, 7, 3 /
    ! OH
    data allrep (338:348) / "Oh", "A1g", "A1u", "A2g", "A2u", "Eg", "Eu", &
         & "T1g", "T1u", "T2g", "T2u" /
    data nallop(307:313) / 3, 10, 29, 287231, 1, 7, 8 /
    ! I
    data allrep (349:354) / "I", "A", "T1", "T2", "G", "H" /
    data nallop(314:319) / 2, 5, 30, 262657, 1, 10 /
    ! IH
    data allrep (355:365) / "Ih", "Ag", "Au", "T1g", "T1u", "T2g", "T2u", "Gg",&
         &  "Gu", "Hg", "Hu" /
    data nallop(320:326) / 3, 10, 31, 344786, 10, 7, 8 /
    ! CV
    data allrep (366:372) / "C*v", "Sig", "Sig-", "Pi", "Del", "Phi", "Gam" /
    data nallop(327:333) / 3, 6, 32, 524340, 12, 11, 5 /
    ! DH
    data allrep (373:385) / "D*h", "Sig", "Siu", "Sig-", "Siu-", "Pig", "Piu", &
         & "Delg", "Delu", "Phig", "Phiu", "Gamg", "Gamu" /
    data nallop(334:341) / 4, 12, 33, 524415, 12, 11, 7, 5 /
    ! R3
    data allrep (386:406) / "R3", "S(g)", "S(u)", "P(g)", "P(u)", "D(g)", &
         & "D(u)", "F(g)", "F(u)", "G(g)", "G(u)", "H(g)", "H(u)", "I(g)", &
         & "I(u)", "K(g)", "K(u)", "L(g)", "L(u)", "M(g)", "M(u)" /
    data nallop(342:348) / 3, 20, 38, 524992, 10, 8, 7 /


end module symmetry_C
