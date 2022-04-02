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

  module Parameters_for_PM7_Sparkles_C
    double precision, dimension(107) ::  gss7sp, alp7sp
    double precision, dimension(107,2) :: gues7sp1, gues7sp2, gues7sp3
!
!       Data for Element  57:                       Lanthanum
!
      data     alp7sp( 57)/        2.72878366D0/
      data     gss7sp( 57)/       51.00521418D0/
      data gues7sp1( 57,1)/        0.27941646D0/
      data gues7sp2( 57,1)/       13.44993801D0/
      data gues7sp3( 57,1)/        3.03636844D0/
      data gues7sp1( 57,2)/        0.22388397D0/
      data gues7sp2( 57,2)/       13.33085266D0/
      data gues7sp3( 57,2)/        3.44240129D0/

!                    Data for Element  58           Cerium
!
      data     alp7sp( 58)/        3.67207094D0/
      data     gss7sp( 58)/       56.61996758D0/
      data gues7sp1( 58,1)/        0.23721649D0/
      data gues7sp2( 58,1)/        9.05598877D0/
      data gues7sp3( 58,1)/        3.05700455D0/
      data gues7sp1( 58,2)/        2.16960030D0/
      data gues7sp2( 58,2)/        8.57474567D0/
      data gues7sp3( 58,2)/        2.97137088D0/

!
!       Data for Element  59:                       Praseodymium
!
      data     alp7sp( 59)/        4.56658977D0/
      data     gss7sp( 59)/       56.12411820D0/
      data gues7sp1( 59,1)/        0.27951368D0/
      data gues7sp2( 59,1)/       14.17545394D0/
      data gues7sp3( 59,1)/        2.96980969D0/
      data gues7sp1( 59,2)/        2.01816742D0/
      data gues7sp2( 59,2)/        5.58007643D0/
      data gues7sp3( 59,2)/        3.13582779D0/

!
!       Data for Element  60:                       Neodymium
!
      data     alp7sp( 60)/        4.70304554D0/
      data     gss7sp( 60)/       51.86812242D0/
      data gues7sp1( 60,1)/        0.37179588D0/
      data gues7sp2( 60,1)/       11.01756863D0/
      data gues7sp3( 60,1)/        2.94211575D0/
      data gues7sp1( 60,2)/        2.01934079D0/
      data gues7sp2( 60,2)/       11.90273308D0/
      data gues7sp3( 60,2)/        3.31388588D0/

!
!       Data for Element  61:                       Promethium
!
      data     alp7sp( 61)/        2.62831297D0/
      data     gss7sp( 61)/       56.60105705D0/
      data gues7sp1( 61,1)/        0.25254367D0/
      data gues7sp2( 61,1)/       11.32762969D0/
      data gues7sp3( 61,1)/        2.91051149D0/
      data gues7sp1( 61,2)/        0.42096080D0/
      data gues7sp2( 61,2)/        9.59944282D0/
      data gues7sp3( 61,2)/        3.51758897D0/

!
!       Data for Element  62:                       Samarium
!
      data     alp7sp( 62)/        2.48035044D0/
      data     gss7sp( 62)/       56.91979708D0/
      data gues7sp1( 62,1)/        0.27745438D0/
      data gues7sp2( 62,1)/       11.39510076D0/
      data gues7sp3( 62,1)/        2.89944391D0/
      data gues7sp1( 62,2)/        0.38744285D0/
      data gues7sp2( 62,2)/        9.47981871D0/
      data gues7sp3( 62,2)/        3.73339418D0/
!
!       Data for Element  63:                       Europium
!
      data     alp7sp( 63)/        2.64501092D0/
      data     gss7sp( 63)/       57.16819322D0/
      data gues7sp1( 63,1)/        0.25290496D0/
      data gues7sp2( 63,1)/       12.42815090D0/
      data gues7sp3( 63,1)/        2.88466231D0/
      data gues7sp1( 63,2)/        0.56909102D0/
      data gues7sp2( 63,2)/        9.75538134D0/
      data gues7sp3( 63,2)/        3.00378806D0/

!
!       Data for Element  64:                       Gadolinium
!
      data     alp7sp( 64)/       4.11940733D0/
      data     gss7sp( 64)/      55.80448082D0/
      data gues7sp1( 64,1)/       0.33731573D0/
      data gues7sp2( 64,1)/       9.55199105D0/
      data gues7sp3( 64,1)/       2.90273188D0/
      data gues7sp1( 64,2)/       2.93935956D0/
      data gues7sp2( 64,2)/       8.93461822D0/
      data gues7sp3( 64,2)/       3.32168472D0/
!
!       Data for Element  65:                       Terbium
!
      data     alp7sp( 65)/       2.71687807D0/
      data     gss7sp( 65)/      57.37215107D0/
      data gues7sp1( 65,1)/       0.30634054D0/
      data gues7sp2( 65,1)/      10.88125130D0/
      data gues7sp3( 65,1)/       2.85828879D0/
      data gues7sp1( 65,2)/       0.66155895D0/
      data gues7sp2( 65,2)/       9.13341532D0/
      data gues7sp3( 65,2)/       3.41583594D0/

!
!       Data for Element  66:                       Dysprosium
!
      data     alp7sp( 66)/       2.82756925D0/
      data     gss7sp( 66)/      55.61069514D0/
      data gues7sp1( 66,1)/       0.35812416D0/
      data gues7sp2( 66,1)/       8.98941131D0/
      data gues7sp3( 66,1)/       2.86815735D0/
      data gues7sp1( 66,2)/       0.40611392D0/
      data gues7sp2( 66,2)/       8.95555466D0/
      data gues7sp3( 66,2)/       3.69733260D0/

!
!       Data for Element  67:                       Holmium
!
      data     alp7sp( 67)/       2.44068845D0/
      data     gss7sp( 67)/      56.68096021D0/
      data gues7sp1( 67,1)/       0.27966404D0/
      data gues7sp2( 67,1)/      11.42737486D0/
      data gues7sp3( 67,1)/       2.83804623D0/
      data gues7sp1( 67,2)/       2.52779706D0/
      data gues7sp2( 67,2)/       9.75311832D0/
      data gues7sp3( 67,2)/       3.41599480D0/
!
!       Data for Element  68:                       Erbium
!
      data     alp7sp( 68)/        3.79457816D0/
      data     gss7sp( 68)/       58.52071878D0/
      data gues7sp1( 68,1)/        0.35009817D0/
      data gues7sp2( 68,1)/        9.49506162D0/
      data gues7sp3( 68,1)/        2.81837640D0/
      data gues7sp1( 68,2)/        0.78014002D0/
      data gues7sp2( 68,2)/        8.47946317D0/
      data gues7sp3( 68,2)/        3.59076847D0/

!
!       Data for Element  69:                       Thulium
!
      data     alp7sp( 69)/       1.95626166D0/
      data     gss7sp( 69)/      55.97220279D0/
      data gues7sp1( 69,1)/       0.36171208D0/
      data gues7sp2( 69,1)/      10.49976955D0/
      data gues7sp3( 69,1)/       2.81614888D0/
      data gues7sp1( 69,2)/       0.63407320D0/
      data gues7sp2( 69,2)/       8.17764346D0/
      data gues7sp3( 69,2)/       3.52086186D0/

!
!       Data for Element  70:                       Ytterbium
!
      data     alp7sp( 70)/        3.82324875D0/
      data     gss7sp( 70)/       56.83760040D0/
      data gues7sp1( 70,1)/        0.40711494D0/
      data gues7sp2( 70,1)/        9.29829027D0/
      data gues7sp3( 70,1)/        2.82195825D0/
      data gues7sp1( 70,2)/        0.52341070D0/
      data gues7sp2( 70,2)/       10.33065656D0/
      data gues7sp3( 70,2)/        3.73059950D0/

!
!       Data for Element  71:                       Lutetium
!
      data     alp7sp( 71)/       5.33103337D0/
      data     gss7sp( 71)/      55.14010416D0/
      data gues7sp1( 71,1)/       0.34482438D0/
      data gues7sp2( 71,1)/      10.07236769D0/
      data gues7sp3( 71,1)/       2.80387860D0/
      data gues7sp1( 71,2)/       1.30538617D0/
      data gues7sp2( 71,2)/       7.21222726D0/
      data gues7sp3( 71,2)/       2.96294393D0/
    end module Parameters_for_PM7_Sparkles_C
