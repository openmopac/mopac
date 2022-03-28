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

      module Parameters_for_PM3_C
      double precision, dimension(107) :: usspm3, upppm3, zspm3, zppm3, &
        betasp, betapp, betadp, alppm3, gsspm3, gsppm3, gpppm3, gp2pm3, &
        hsppm3, polvolpm3, zsnpm3, zpnpm3, zdnpm3, f0sdpm3, g2sdpm3
      double precision, dimension(107,4) :: guesp1, guesp2, guesp3
!
!       Data for Element   1:     Hydrogen
!
      data      usspm3(  1)/     -13.0733210d0/
      data      betasp(  1)/      -5.6265120d0/
      data       zspm3(  1)/       0.9678070d0/
      data      alppm3(  1)/       3.3563860d0/
      data      gsspm3(  1)/      14.7942080d0/
      data   polvolpm3(  1)/       0.1785810d0/
      data    guesp1(  1,1)/       1.1287500d0/
      data    guesp2(  1,1)/       5.0962820d0/
      data    guesp3(  1,1)/       1.5374650d0/
      data    guesp1(  1,2)/      -1.0603290d0/
      data    guesp2(  1,2)/       6.0037880d0/
      data    guesp3(  1,2)/       1.5701890d0/
!
!       Data for Element   2:       Helium
!
      data      usspm3(  2)/     -35.4931170d0/
      data      upppm3(  2)/       9.9998433d0/
      data      betasp(  2)/     -17.5545020d0/
      data      betapp(  2)/     -26.0526064d0/
      data       zspm3(  2)/       1.7710761d0/
      data       zppm3(  2)/       6.9018258d0/
      data      alppm3(  2)/       5.7813931d0/
      data      gsspm3(  2)/       8.0590135d0/
      data      gsppm3(  2)/      11.0152229d0/
      data      gpppm3(  2)/      19.3843339d0/
      data      gp2pm3(  2)/      11.3517224d0/
      data      hsppm3(  2)/       0.5321205d0/
      data    guesp1(  2,1)/       0.1985976d0/
      data    guesp2(  2,1)/       1.0009593d0/
      data    guesp3(  2,1)/       0.5031509d0/
!
!       Data for Element   3:      Lithium
!
      data      usspm3(  3)/      -5.3000000d0/
      data      upppm3(  3)/      -3.4000000d0/
      data      betasp(  3)/      -0.5500000d0/
      data      betapp(  3)/      -1.5000000d0/
      data       zspm3(  3)/       0.6500000d0/
      data       zppm3(  3)/       0.7500000d0/
      data      alppm3(  3)/       1.2550000d0/
      data      gsspm3(  3)/       4.5000000d0/
      data      gsppm3(  3)/       3.0000000d0/
      data      gpppm3(  3)/       5.2500000d0/
      data      gp2pm3(  3)/       4.5000000d0/
      data      hsppm3(  3)/       0.1500000d0/
      data    guesp1(  3,1)/      -0.4500000d0/
      data    guesp2(  3,1)/       5.0000000d0/
      data    guesp3(  3,1)/       1.0000000d0/
      data    guesp1(  3,2)/       0.8000000d0/
      data    guesp2(  3,2)/       6.5000000d0/
      data    guesp3(  3,2)/       1.0000000d0/
!
!       Data for Element   4:    Beryllium
!
      data      usspm3(  4)/     -17.2647520d0/
      data      upppm3(  4)/     -11.3042430d0/
      data      betasp(  4)/      -3.9620530d0/
      data      betapp(  4)/      -2.7806840d0/
      data       zspm3(  4)/       0.8774390d0/
      data       zppm3(  4)/       1.5087550d0/
      data      alppm3(  4)/       1.5935360d0/
      data      gsspm3(  4)/       9.0128510d0/
      data      gsppm3(  4)/       6.5761990d0/
      data      gpppm3(  4)/       6.0571820d0/
      data      gp2pm3(  4)/       9.0052190d0/
      data      hsppm3(  4)/       0.5446790d0/
      data    guesp1(  4,1)/       1.6315720d0/
      data    guesp2(  4,1)/       2.6729620d0/
      data    guesp3(  4,1)/       1.7916860d0/
      data    guesp1(  4,2)/      -2.1109590d0/
      data    guesp2(  4,2)/       1.9685940d0/
      data    guesp3(  4,2)/       1.7558710d0/
!
!       Data for Element   5:        Boron
!
      data      usspm3(  5)/     -50.4776829d0/
      data      upppm3(  5)/     -37.4119835d0/
      data      betasp(  5)/     -10.5497263d0/
      data      betapp(  5)/      -3.9995953d0/
      data       zspm3(  5)/       1.5312597d0/
      data       zppm3(  5)/       1.1434597d0/
      data      alppm3(  5)/       2.2104163d0/
      data      gsspm3(  5)/      18.2782796d0/
      data      gsppm3(  5)/      15.3330673d0/
      data      gpppm3(  5)/      12.3158582d0/
      data      gp2pm3(  5)/      11.1785351d0/
      data      hsppm3(  5)/       0.5997885d0/
      data    guesp1(  5,1)/      -0.3518407d0/
      data    guesp2(  5,1)/       3.0008621d0/
      data    guesp3(  5,1)/       0.8241176d0/
!
!       Data for Element   6:       Carbon
!
      data      usspm3(  6)/     -47.2703200d0/
      data      upppm3(  6)/     -36.2669180d0/
      data      betasp(  6)/     -11.9100150d0/
      data      betapp(  6)/      -9.8027550d0/
      data       zspm3(  6)/       1.5650850d0/
      data       zppm3(  6)/       1.8423450d0/
      data      alppm3(  6)/       2.7078070d0/
      data      gsspm3(  6)/      11.2007080d0/
      data      gsppm3(  6)/      10.2650270d0/
      data      gpppm3(  6)/      10.7962920d0/
      data      gp2pm3(  6)/       9.0425660d0/
      data      hsppm3(  6)/       2.2909800d0/
      data   polvolpm3(  6)/       1.2427300d0/
      data    guesp1(  6,1)/       0.0501070d0/
      data    guesp2(  6,1)/       6.0031650d0/
      data    guesp3(  6,1)/       1.6422140d0/
      data    guesp1(  6,2)/       0.0507330d0/
      data    guesp2(  6,2)/       6.0029790d0/
      data    guesp3(  6,2)/       0.8924880d0/
!
!       Data for Element   7:     Nitrogen
!
      data      usspm3(  7)/     -49.3356720d0/
      data      upppm3(  7)/     -47.5097360d0/
      data      betasp(  7)/     -14.0625210d0/
      data      betapp(  7)/     -20.0438480d0/
      data       zspm3(  7)/       2.0280940d0/
      data       zppm3(  7)/       2.3137280d0/
      data      alppm3(  7)/       2.8305450d0/
      data      gsspm3(  7)/      11.9047870d0/
      data      gsppm3(  7)/       7.3485650d0/
      data      gpppm3(  7)/      11.7546720d0/
      data      gp2pm3(  7)/      10.8072770d0/
      data      hsppm3(  7)/       1.1367130d0/
      data   polvolpm3(  7)/       0.9201900d0/
      data    guesp1(  7,1)/       1.5016740d0/
      data    guesp2(  7,1)/       5.9011480d0/
      data    guesp3(  7,1)/       1.7107400d0/
      data    guesp1(  7,2)/      -1.5057720d0/
      data    guesp2(  7,2)/       6.0046580d0/
      data    guesp3(  7,2)/       1.7161490d0/
!
!       Data for Element   8:       Oxygen
!
      data      usspm3(  8)/     -86.9930020d0/
      data      upppm3(  8)/     -71.8795800d0/
      data      betasp(  8)/     -45.2026510d0/
      data      betapp(  8)/     -24.7525150d0/
      data       zspm3(  8)/       3.7965440d0/
      data       zppm3(  8)/       2.3894020d0/
      data      alppm3(  8)/       3.2171020d0/
      data      gsspm3(  8)/      15.7557600d0/
      data      gsppm3(  8)/      10.6211600d0/
      data      gpppm3(  8)/      13.6540160d0/
      data      gp2pm3(  8)/      12.4060950d0/
      data      hsppm3(  8)/       0.5938830d0/
      data   polvolpm3(  8)/       0.4346050d0/
      data    guesp1(  8,1)/      -1.1311280d0/
      data    guesp2(  8,1)/       6.0024770d0/
      data    guesp3(  8,1)/       1.6073110d0/
      data    guesp1(  8,2)/       1.1378910d0/
      data    guesp2(  8,2)/       5.9505120d0/
      data    guesp3(  8,2)/       1.5983950d0/
!
!       Data for Element   9:     Fluorine
!
      data      usspm3(  9)/    -110.4353030d0/
      data      upppm3(  9)/    -105.6850470d0/
      data      betasp(  9)/     -48.4059390d0/
      data      betapp(  9)/     -27.7446600d0/
      data       zspm3(  9)/       4.7085550d0/
      data       zppm3(  9)/       2.4911780d0/
      data      alppm3(  9)/       3.3589210d0/
      data      gsspm3(  9)/      10.4966670d0/
      data      gsppm3(  9)/      16.0736890d0/
      data      gpppm3(  9)/      14.8172560d0/
      data      gp2pm3(  9)/      14.4183930d0/
      data      hsppm3(  9)/       0.7277630d0/
      data   polvolpm3(  9)/       0.2358930d0/
      data    guesp1(  9,1)/      -0.0121660d0/

data    guesp2(  9,1)/       6.0235740d0/
      data    guesp3(  9,1)/       1.8568590d0/
      data    guesp1(  9,2)/      -0.0028520d0/
      data    guesp2(  9,2)/       6.0037170d0/
      data    guesp3(  9,2)/       2.6361580d0/
!
!       Data for Element  10:         Neon
!
      data      usspm3( 10)/       9.6554736d0/
      data      upppm3( 10)/     -71.1466474d0/
      data      betasp( 10)/      -0.1512657d0/
      data      betapp( 10)/     -23.8470187d0/
      data       zspm3( 10)/       5.9998745d0/
      data       zppm3( 10)/       4.1752600d0/
      data      alppm3( 10)/       2.5377255d0/
      data      gsspm3( 10)/       0.4989374d0/
      data      gsppm3( 10)/      10.1592044d0/
      data      gpppm3( 10)/      18.9449976d0/
      data      gp2pm3( 10)/       8.5648589d0/
      data      hsppm3( 10)/       0.2999991d0/
      data    guesp1( 10,1)/       0.2352109d0/
      data    guesp2( 10,1)/       8.8116607d0/
      data    guesp3( 10,1)/       1.0927884d0/
!
!       Data for Element  11:       Sodium
!
      data      usspm3( 11)/      -5.0140594d0/
      data      upppm3( 11)/      -2.9252301d0/
      data      betasp( 11)/      -3.9486341d0/
      data      betapp( 11)/      -4.2401809d0/
      data       zspm3( 11)/       2.6618938d0/
      data       zppm3( 11)/       0.8837425d0/
      data      alppm3( 11)/       0.9109177d0/
      data      gsspm3( 11)/       5.7975247d0/
      data      gsppm3( 11)/      10.0000062d0/
      data      gpppm3( 11)/       1.2247903d0/
      data      gp2pm3( 11)/       0.9996497d0/
      data      hsppm3( 11)/       0.3999363d0/
!
!       Data for Element  12:    Magnesium
!
      data      usspm3( 12)/     -14.6236880d0/
      data      upppm3( 12)/     -14.1734600d0/
      data      betasp( 12)/      -2.0716910d0/
      data      betapp( 12)/      -0.5695810d0/
      data       zspm3( 12)/       0.6985520d0/
      data       zppm3( 12)/       1.4834530d0/
      data      alppm3( 12)/       1.3291470d0/
      data      gsspm3( 12)/       6.6943000d0/
      data      gsppm3( 12)/       6.7939950d0/
      data      gpppm3( 12)/       6.9104460d0/
      data      gp2pm3( 12)/       7.0908230d0/
      data      hsppm3( 12)/       0.5433000d0/
      data    guesp1( 12,1)/       2.1170500d0/
      data    guesp2( 12,1)/       6.0094770d0/
      data    guesp3( 12,1)/       2.0844060d0/
      data    guesp1( 12,2)/      -2.5477670d0/
      data    guesp2( 12,2)/       4.3953700d0/
      data    guesp3( 12,2)/       2.0636740d0/
!
!       Data for Element  13:    Aluminium
!
      data      usspm3( 13)/     -24.8454040d0/
      data      upppm3( 13)/     -22.2641590d0/
      data      betasp( 13)/      -0.5943010d0/
      data      betapp( 13)/      -0.9565500d0/
      data       zspm3( 13)/       1.7028880d0/
      data       zppm3( 13)/       1.0736290d0/
      data      alppm3( 13)/       1.5217030d0/
      data      gsspm3( 13)/       5.7767370d0/
      data      gsppm3( 13)/      11.6598560d0/
      data      gpppm3( 13)/       6.3477900d0/
      data      gp2pm3( 13)/       6.1210770d0/
      data      hsppm3( 13)/       4.0062450d0/
      data    guesp1( 13,1)/      -0.4730900d0/
      data    guesp2( 13,1)/       1.9158250d0/
      data    guesp3( 13,1)/       1.4517280d0/
      data    guesp1( 13,2)/      -0.1540510d0/
      data    guesp2( 13,2)/       6.0050860d0/
      data    guesp3( 13,2)/       2.5199970d0/
!
!       Data for Element  14:      Silicon
!
      data      usspm3( 14)/     -26.7634830d0/
      data      upppm3( 14)/     -22.8136350d0/
      data      betasp( 14)/      -2.8621450d0/
      data      betapp( 14)/      -3.9331480d0/
      data       zspm3( 14)/       1.6350750d0/
      data       zppm3( 14)/       1.3130880d0/
      data      alppm3( 14)/       2.1358090d0/
      data      gsspm3( 14)/       5.0471960d0/
      data      gsppm3( 14)/       5.9490570d0/
      data      gpppm3( 14)/       6.7593670d0/
      data      gp2pm3( 14)/       5.1612970d0/
      data      hsppm3( 14)/       0.9198320d0/
      data    guesp1( 14,1)/      -0.3906000d0/
      data    guesp2( 14,1)/       6.0000540d0/
      data    guesp3( 14,1)/       0.6322620d0/
      data    guesp1( 14,2)/       0.0572590d0/
      data    guesp2( 14,2)/       6.0071830d0/
      data    guesp3( 14,2)/       2.0199870d0/
!
!       Data for Element  15:   Phosphorus
!
      data      usspm3( 15)/     -40.4130960d0/
      data      upppm3( 15)/     -29.5930520d0/
      data      betasp( 15)/     -12.6158790d0/
      data      betapp( 15)/      -4.1600400d0/
      data       zspm3( 15)/       2.0175630d0/
      data       zppm3( 15)/       1.5047320d0/
      data      alppm3( 15)/       1.9405340d0/
      data      gsspm3( 15)/       7.8016150d0/
      data      gsppm3( 15)/       5.1869490d0/
      data      gpppm3( 15)/       6.6184780d0/
      data      gp2pm3( 15)/       6.0620020d0/
      data      hsppm3( 15)/       1.5428090d0/
      data    guesp1( 15,1)/      -0.6114210d0/
      data    guesp2( 15,1)/       1.9972720d0/
      data    guesp3( 15,1)/       0.7946240d0/
      data    guesp1( 15,2)/      -0.0939350d0/
      data    guesp2( 15,2)/       1.9983600d0/
      data    guesp3( 15,2)/       1.9106770d0/
!
!       Data for Element  16:       Sulfur
!
      data      usspm3( 16)/     -49.8953710d0/
      data      upppm3( 16)/     -44.3925830d0/
      data      betasp( 16)/      -8.8274650d0/
      data      betapp( 16)/      -8.0914150d0/
      data       zspm3( 16)/       1.8911850d0/
      data       zppm3( 16)/       1.6589720d0/
      data      alppm3( 16)/       2.2697060d0/
      data      gsspm3( 16)/       8.9646670d0/
      data      gsppm3( 16)/       6.7859360d0/
      data      gpppm3( 16)/       9.9681640d0/
      data      gp2pm3( 16)/       7.9702470d0/
      data      hsppm3( 16)/       4.0418360d0/
      data    guesp1( 16,1)/      -0.3991910d0/
      data    guesp2( 16,1)/       6.0006690d0/
      data    guesp3( 16,1)/       0.9621230d0/
      data    guesp1( 16,2)/      -0.0548990d0/
      data    guesp2( 16,2)/       6.0018450d0/
      data    guesp3( 16,2)/       1.5799440d0/
!
!       Data for Element  17:     Chlorine
!
      data      usspm3( 17)/    -100.6267470d0/
      data      upppm3( 17)/     -53.6143960d0/
      data      betasp( 17)/     -27.5285600d0/
      data      betapp( 17)/     -11.5939220d0/
      data       zspm3( 17)/       2.2462100d0/
      data       zppm3( 17)/       2.1510100d0/
      data      alppm3( 17)/       2.5172960d0/
      data      gsspm3( 17)/      16.0136010d0/
      data      gsppm3( 17)/       8.0481150d0/
      data      gpppm3( 17)/       7.5222150d0/
      data      gp2pm3( 17)/       7.5041540d0/
      data      hsppm3( 17)/       3.4811530d0/
      data   polvolpm3( 17)/       1.8556000d0/
      data    guesp1( 17,1)/      -0.1715910d0/
      data    guesp2( 17,1)/       6.0008020d0/
      data    guesp3( 17,1)/       1.0875020d0/
      data    guesp1( 17,2)/      -0.0134580d0/
      data    guesp2( 17,2)/       1.9666180d0/
      data    guesp3( 17,2)/       2.2928910d0/
!
!       Data for Element  18:        Argon
!
      data      usspm3( 18)/       4.2759913d0/
      data      upppm3( 18)/     -75.1004104d0/
      data      betasp( 18)/      -2.4320045d0/
      data      betapp( 18)/     -22.9856937d0/
      data       zspm3( 18)/       0.9821697d0/
      data       zppm3( 18)/       5.9997150d0/
      data      alppm3( 18)/       2.2113324d0/
      data      gsspm3( 18)/       3.2377414d0/
      data      gsppm3( 18)/      16.3573969d0/
      data      gpppm3( 18)/      19.9984646d0/
      data      gp2pm3( 18)/      10.9343238d0/
      data      hsppm3( 18)/       0.8169151d0/
      data    guesp1( 18,1)/      -0.7376503d0/
      data    guesp2( 18,1)/       3.8882543d0/
      data    guesp3( 18,1)/       0.7150665d0/
!
!       Data for Element  19:    Potassium
!
      data      usspm3( 19)/      -4.2596475d0/
      data      upppm3( 19)/      -2.6425087d0/
      data      betasp( 19)/      -0.4250386d0/
      data      betapp( 19)/      -3.1998137d0/
      data       zspm3( 19)/       0.8101687d0/
      data       zppm3( 19)/       0.9578342d0/
      data      alppm3( 19)/       0.7252158d0/
      data      gsspm3( 19)/       6.7788311d0/
      data      gsppm3( 19)/       9.3472796d0/
      data      gpppm3( 19)/       3.4963514d0/
      data      gp2pm3( 19)/       2.7416795d0/
      data      hsppm3( 19)/       1.6592458d0/
!
!       Data for Element  20:      Calcium
!
      data      usspm3( 20)/     -11.9608156d0/
      data      upppm3( 20)/      -8.6350859d0/
      data      betasp( 20)/      -0.9868271d0/
      data      betapp( 20)/      -3.9684267d0/
      data       zspm3( 20)/       1.2087415d0/
      data       zppm3( 20)/       0.9409370d0/
      data      alppm3( 20)/       0.4999979d0/
      data      gsspm3( 20)/       5.9719910d0/
      data      gsppm3( 20)/       4.9607780d0/
      data      gpppm3( 20)/       3.7214897d0/
      data      gp2pm3( 20)/       3.7116499d0/
      data      hsppm3( 20)/       0.7928150d0/
!
!       Data for Element  30:         Zinc
!
      data      usspm3( 30)/     -18.5321980d0/
      data      upppm3( 30)/     -11.0474090d0/
      data      betasp( 30)/      -0.7155780d0/
      data      betapp( 30)/      -6.3518640d0/
      data       zspm3( 30)/       1.8199890d0/
      data       zppm3( 30)/       1.5069220d0/
      data      alppm3( 30)/       1.3501260d0/
      data      gsspm3( 30)/       9.6771960d0/
      data      gsppm3( 30)/       7.7362040d0/
      data      gpppm3( 30)/       4.9801740d0/
      data      gp2pm3( 30)/       4.6696560d0/
      data      hsppm3( 30)/       0.6004130d0/
      data    guesp1( 30,1)/      -0.1112340d0/
      data    guesp2( 30,1)/       6.0014780d0/
      data    guesp3( 30,1)/       1.5160320d0/
      data    guesp1( 30,2)/      -0.1323700d0/
      data    guesp2( 30,2)/       1.9958390d0/
      data    guesp3( 30,2)/       2.5196420d0/
!
!       Data for Element  31:      Gallium
!
      data      usspm3( 31)/     -29.8555930d0/
      data      upppm3( 31)/     -21.8753710d0/
      data      betasp( 31)/      -4.9456180d0/
      data      betapp( 31)/      -0.4070530d0/
      data       zspm3( 31)/       1.8470400d0/
      data       zppm3( 31)/       0.8394110d0/
      data      alppm3( 31)/       1.6051150d0/
      data      gsspm3( 31)/       8.4585540d0/
      data      gsppm3( 31)/       8.9256190d0/
      data      gpppm3( 31)/       5.0868550d0/
      data      gp2pm3( 31)/       4.9830450d0/
      data      hsppm3( 31)/       2.0512600d0/
      data    guesp1( 31,1)/      -0.5601790d0/
      data    guesp2( 31,1)/       5.6232730d0/
      data    guesp3( 31,1)/       1.5317800d0/
      data    guesp1( 31,2)/      -0.2727310d0/
      data    guesp2( 31,2)/       1.9918430d0/
      data    guesp3( 31,2)/       2.1838640d0/
!
!       Data for Element  32:    Germanium
!
      data      usspm3( 32)/     -35.4671955d0/
      data      upppm3( 32)/     -31.5863583d0/
      data      betasp( 32)/      -5.3250024d0/
      data      betapp( 32)/      -2.2501567d0/
      data       zspm3( 32)/       2.2373526d0/
      data       zppm3( 32)/       1.5924319d0/
      data      alppm3( 32)/       1.9723370d0/
      data      gsspm3( 32)/       5.3769635d0/
      data      gsppm3( 32)/      10.2095293d0/
      data      gpppm3( 32)/       7.6718647d0/
      data      gp2pm3( 32)/       6.9242663d0/
      data      hsppm3( 32)/       1.3370204d0/
      data    guesp1( 32,1)/       0.9631726d0/
      data    guesp2( 32,1)/       6.0120134d0/
      data    guesp3( 32,1)/       2.1633655d0/
      data    guesp1( 32,2)/      -0.9593891d0/
      data    guesp2( 32,2)/       5.7491802d0/
      data    guesp3( 32,2)/       2.1693724d0/
!
!       Data for Element  33:      Arsenic
!
      data      usspm3( 33)/     -38.5074240d0/
      data      upppm3( 33)/     -35.1524150d0/
      data      betasp( 33)/      -8.2321650d0/
      data      betapp( 33)/      -5.0173860d0/
      data       zspm3( 33)/       2.6361770d0/
      data       zppm3( 33)/       1.7038890d0/
      data      alppm3( 33)/       1.7944770d0/
      data      gsspm3( 33)/       8.7890010d0/
      data      gsppm3( 33)/       5.3979830d0/
      data      gpppm3( 33)/       8.2872500d0/
      data      gp2pm3( 33)/       8.2103460d0/
      data      hsppm3( 33)/       1.9510340d0/
      data    guesp1( 33,1)/      -0.4600950d0/
      data    guesp2( 33,1)/       1.9831150d0/
      data    guesp3( 33,1)/       1.0867930d0/
      data    guesp1( 33,2)/      -0.0889960d0/
      data    guesp2( 33,2)/       1.9929440d0/
      data    guesp3( 33,2)/       2.1400580d0/
!
!       Data for Element  34:     Selenium
!
      data      usspm3( 34)/     -55.3781350d0/
      data      upppm3( 34)/     -49.8230760d0/
      data      betasp( 34)/      -6.1578220d0/
      data      betapp( 34)/      -5.4930390d0/
      data       zspm3( 34)/       2.8280510d0/
      data       zppm3( 34)/       1.7325360d0/
      data      alppm3( 34)/       3.0439570d0/
      data      gsspm3( 34)/       7.4325910d0/
      data      gsppm3( 34)/      10.0604610d0/
      data      gpppm3( 34)/       9.5683260d0/
      data      gp2pm3( 34)/       7.7242890d0/
      data      hsppm3( 34)/       4.0165580d0/
      data    guesp1( 34,1)/       0.0478730d0/
      data    guesp2( 34,1)/       6.0074000d0/
      data    guesp3( 34,1)/       2.0817170d0/
      data    guesp1( 34,2)/       0.1147200d0/
      data    guesp2( 34,2)/       6.0086720d0/
      data    guesp3( 34,2)/       1.5164230d0/
!
!       Data for Element  35:      Bromine
!
      data      usspm3( 35)/    -116.6193110d0/
      data      upppm3( 35)/     -74.2271290d0/
      data      betasp( 35)/     -31.1713420d0/
      data      betapp( 35)/      -6.8140130d0/
      data       zspm3( 35)/       5.3484570d0/
      data       zppm3( 35)/       2.1275900d0/
      data      alppm3( 35)/       2.5118420d0/
      data      gsspm3( 35)/      15.9434250d0/
      data      gsppm3( 35)/      16.0616800d0/
      data      gpppm3( 35)/       8.2827630d0/
      data      gp2pm3( 35)/       7.8168490d0/
      data      hsppm3( 35)/       0.5788690d0/
      data   polvolpm3( 35)/       2.8454200d0/
      data    guesp1( 35,1)/       0.9604580d0/
      data    guesp2( 35,1)/       5.9765080d0/
      data    guesp3( 35,1)/       2.3216540d0/
      data    guesp1( 35,2)/      -0.9549160d0/
      data    guesp2( 35,2)/       5.9447030d0/
      data    guesp3( 35,2)/       2.3281420d0/
!
!       Data for Element  36:      Krypton
!
      data      usspm3( 36)/       9.7986937d0/
      data      upppm3( 36)/     -73.8595661d0/
      data      betasp( 36)/      -1.2094707d0/
      data      betapp( 36)/      -8.4328271d0/
      data       zspm3( 36)/       3.5608622d0/
      data       zppm3( 36)/       1.9832062d0/
      data      alppm3( 36)/       1.7025435d0/
      data      gsspm3( 36)/       0.4994169d0/
      data      gsppm3( 36)/       9.5261258d0/
      data      gpppm3( 36)/      19.9987286d0/
      data      gp2pm3( 36)/      11.1986052d0/
      data      hsppm3( 36)/       2.1818338d0/
      data    guesp1( 36,1)/       0.7804390d0/
      data    guesp2( 36,1)/       8.8371397d0/
      data    guesp3( 36,1)/       0.5000419d0/
!
!       Data for Element  37:     Rubidium
!
      data      usspm3( 37)/      -4.5920539d0/
      data      upppm3( 37)/      -3.0119211d0/
      data      betasp( 37)/     -12.0894660d0/
      data      betapp( 37)/      -1.9999268d0/
      data       zspm3( 37)/       4.0000415d0/
      data       zppm3( 37)/       1.0134590d0/
      data      alppm3( 37)/       0.9985324d0/
      data      gsspm3( 37)/       9.2757213d0/
      data      gsppm3( 37)/      20.0000037d0/
      data      gpppm3( 37)/      13.3694136d0/
      data      gp2pm3( 37)/      19.0000035d0/
      data      hsppm3( 37)/       4.9987146d0/
!
!       Data for Element  38:    Strontium
!
      data      usspm3( 38)/     -10.9183439d0/
      data      upppm3( 38)/      -7.9835223d0/
      data      betasp( 38)/     -10.0000071d0/
      data      betapp( 38)/      -5.6755392d0/
      data       zspm3( 38)/       1.2794532d0/
      data       zppm3( 38)/       1.3912500d0/
      data      alppm3( 38)/       1.6467554d0/
      data      gsspm3( 38)/       5.0361465d0/
      data      gsppm3( 38)/       3.9284011d0/
      data      gpppm3( 38)/       3.1533139d0/
      data      gp2pm3( 38)/       3.2457312d0/
      data      hsppm3( 38)/       0.7596212d0/
!
!       Data for Element  48:      Cadmium
!
      data      usspm3( 48)/     -15.8285840d0/
      data      upppm3( 48)/       8.7497950d0/
      data      betasp( 48)/      -8.5819440d0/
      data      betapp( 48)/      -0.6010340d0/
      data       zspm3( 48)/       1.6793510d0/
      data       zppm3( 48)/       2.0664120d0/
      data      alppm3( 48)/       1.5253820d0/
      data      gsspm3( 48)/       9.2069600d0/
      data      gsppm3( 48)/       8.2315390d0/
      data      gpppm3( 48)/       4.9481040d0/
      data      gp2pm3( 48)/       4.6696560d0/
      data      hsppm3( 48)/       1.6562340d0/
!
!       Data for Element  49:       Indium
!
      data      usspm3( 49)/     -26.1762050d0/
      data      upppm3( 49)/     -20.0058220d0/
      data      betasp( 49)/      -2.9933190d0/
      data      betapp( 49)/      -1.8289080d0/
      data       zspm3( 49)/       2.0161160d0/
      data       zppm3( 49)/       1.4453500d0/
      data      alppm3( 49)/       1.4183850d0/
      data      gsspm3( 49)/       6.5549000d0/
      data      gsppm3( 49)/       8.2298730d0/
      data      gpppm3( 49)/       6.2992690d0/
      data      gp2pm3( 49)/       4.9842110d0/
      data      hsppm3( 49)/       2.6314610d0/
      data    guesp1( 49,1)/      -0.3431380d0/
      data    guesp2( 49,1)/       1.9940340d0/
      data    guesp3( 49,1)/       1.6255160d0/
      data    guesp1( 49,2)/      -0.1095320d0/
      data    guesp2( 49,2)/       5.6832170d0/
      data    guesp3( 49,2)/       2.8670090d0/
!
!       Data for Element  50:          Tin
!
      data      usspm3( 50)/     -34.5501920d0/
      data      upppm3( 50)/     -25.8944190d0/
      data      betasp( 50)/      -2.7858020d0/
      data      betapp( 50)/      -2.0059990d0/
      data       zspm3( 50)/       2.3733280d0/
      data       zppm3( 50)/       1.6382330d0/
      data      alppm3( 50)/       1.6996500d0/
      data      gsspm3( 50)/      10.1900330d0/
      data      gsppm3( 50)/       7.2353270d0/
      data      gpppm3( 50)/       5.6738100d0/
      data      gp2pm3( 50)/       5.1822140d0/
      data      hsppm3( 50)/       1.0331570d0/
      data    guesp1( 50,1)/      -0.1503530d0/
      data    guesp2( 50,1)/       6.0056940d0/
      data    guesp3( 50,1)/       1.7046420d0/
      data    guesp1( 50,2)/      -0.0444170d0/
      data    guesp2( 50,2)/       2.2573810d0/
      data    guesp3( 50,2)/       2.4698690d0/
!
!       Data for Element  51:     Antimony
!
      data      usspm3( 51)/     -56.4321960d0/
      data      upppm3( 51)/     -29.4349540d0/
      data      betasp( 51)/     -14.7942170d0/
      data      betapp( 51)/      -2.8179480d0/
      data       zspm3( 51)/       2.3430390d0/
      data       zppm3( 51)/       1.8999920d0/
      data      alppm3( 51)/       2.0343010d0/
      data      gsspm3( 51)/       9.2382770d0/
      data      gsppm3( 51)/       5.2776800d0/
      data      gpppm3( 51)/       6.3500000d0/
      data      gp2pm3( 51)/       6.2500000d0/
      data      hsppm3( 51)/       2.4244640d0/
      data    guesp1( 51,1)/       3.0020280d0/
      data    guesp2( 51,1)/       6.0053420d0/
      data    guesp3( 51,1)/       0.8530600d0/
      data    guesp1( 51,2)/      -0.0188920d0/
      data    guesp2( 51,2)/       6.0114780d0/
      data    guesp3( 51,2)/       2.7933110d0/
!
!       Data for Element  52:    Tellurium
!
      data      usspm3( 52)/     -44.9380360d0/
      data      upppm3( 52)/     -46.3140990d0/
      data      betasp( 52)/      -2.6651460d0/
      data      betapp( 52)/      -3.8954300d0/
      data       zspm3( 52)/       4.1654920d0/
      data       zppm3( 52)/       1.6475550d0/
      data      alppm3( 52)/       2.4850190d0/
      data      gsspm3( 52)/      10.2550730d0/
      data      gsppm3( 52)/       8.1691450d0/
      data      gpppm3( 52)/       7.7775920d0/
      data      gp2pm3( 52)/       7.7551210d0/
      data      hsppm3( 52)/       3.7724620d0/
      data    guesp1( 52,1)/       0.0333910d0/
      data    guesp2( 52,1)/       5.9563790d0/
      data    guesp3( 52,1)/       2.2775750d0/
      data    guesp1( 52,2)/      -1.9218670d0/
      data    guesp2( 52,2)/       4.9732190d0/
      data    guesp3( 52,2)/       0.5242430d0/
!
!       Data for Element  53:       Iodine
!
      data      usspm3( 53)/     -96.4540370d0/
      data      upppm3( 53)/     -61.0915820d0/
      data      betasp( 53)/     -14.4942340d0/
      data      betapp( 53)/      -5.8947030d0/
      data       zspm3( 53)/       7.0010130d0/
      data       zppm3( 53)/       2.4543540d0/
      data      alppm3( 53)/       1.9901850d0/
      data      gsspm3( 53)/      13.6319430d0/
      data      gsppm3( 53)/      14.9904060d0/
      data      gpppm3( 53)/       7.2883300d0/
      data      gp2pm3( 53)/       5.9664070d0/
      data      hsppm3( 53)/       2.6300350d0/
      data   polvolpm3( 53)/       4.7997300d0/
      data    guesp1( 53,1)/      -0.1314810d0/
      data    guesp2( 53,1)/       5.2064170d0/
      data    guesp3( 53,1)/       1.7488240d0/
      data    guesp1( 53,2)/      -0.0368970d0/
      data    guesp2( 53,2)/       6.0101170d0/
      data    guesp3( 53,2)/       2.7103730d0/
!
!       Data for Element  54:        Xenon
!
      data      usspm3( 54)/       5.9205889d0/
      data      upppm3( 54)/     -86.9828315d0/
      data      betasp( 54)/      -3.6048445d0/
      data      betapp( 54)/      -4.9673738d0/
      data       zspm3( 54)/       4.9900791d0/
      data       zppm3( 54)/       2.6929255d0/
      data      alppm3( 54)/       1.7948555d0/
      data      gsspm3( 54)/       2.1874755d0/
      data      gsppm3( 54)/       4.9038680d0/
      data      gpppm3( 54)/      12.4939185d0/
      data      gp2pm3( 54)/      14.9527902d0/
      data      hsppm3( 54)/       2.4768528d0/
      data    guesp1( 54,1)/      -0.2455348d0/
      data    guesp2( 54,1)/       2.0580083d0/
      data    guesp3( 54,1)/       1.7173301d0/
!
!       Data for Element  55:       Cesium
!
      data      usspm3( 55)/      -3.2036564d0/
      data      upppm3( 55)/      -1.7451970d0/
      data      betasp( 55)/      -0.6028761d0/
      data      betapp( 55)/      -5.9386091d0/
      data       zspm3( 55)/       3.5960298d0/
      data       zppm3( 55)/       0.9255168d0/
      data      alppm3( 55)/       0.5238336d0/
      data      gsspm3( 55)/       2.1605658d0/
      data      gsppm3( 55)/       4.1667498d0/
      data      gpppm3( 55)/       5.4140484d0/
      data      gp2pm3( 55)/       6.2904939d0/
      data      hsppm3( 55)/       0.3995561d0/
!
!       Data for Element  56:       Barium
!
      data      usspm3( 56)/     -10.1028101d0/
      data      upppm3( 56)/      -6.6743838d0/
      data      betasp( 56)/     -10.0000046d0/
      data      betapp( 56)/     -10.0000103d0/
      data       zspm3( 56)/       1.9258219d0/
      data       zppm3( 56)/       1.4519912d0/
      data      alppm3( 56)/       0.4999966d0/
      data      gsspm3( 56)/       4.8372407d0/
      data      gsppm3( 56)/       3.1942225d0/
      data      gpppm3( 56)/       2.1243276d0/
      data      gp2pm3( 56)/       2.2153940d0/
      data      hsppm3( 56)/       0.3999242d0/
!
!       Data for Element  80:      Mercury
!
      data      usspm3( 80)/     -17.7622290d0/
      data      upppm3( 80)/     -18.3307510d0/
      data      betasp( 80)/      -3.1013650d0/
      data      betapp( 80)/      -3.4640310d0/
      data       zspm3( 80)/       1.4768850d0/
      data       zppm3( 80)/       2.4799510d0/
      data      alppm3( 80)/       1.5293770d0/
      data      gsspm3( 80)/       6.6247200d0/
      data      gsppm3( 80)/      10.6392970d0/
      data      gpppm3( 80)/      14.7092830d0/
      data      gp2pm3( 80)/      16.0007400d0/
      data      hsppm3( 80)/       2.0363110d0/
      data    guesp1( 80,1)/       1.0827200d0/
      data    guesp2( 80,1)/       6.4965980d0/
      data    guesp3( 80,1)/       1.1951460d0/
      data    guesp1( 80,2)/      -0.0965530d0/
      data    guesp2( 80,2)/       3.9262810d0/
      data    guesp3( 80,2)/       2.6271600d0/
!
!       Data for Element  81:     Thallium
!
      data      usspm3( 81)/     -30.0531700d0/
      data      upppm3( 81)/     -26.9206370d0/
      data      betasp( 81)/      -1.0844950d0/
      data      betapp( 81)/      -7.9467990d0/
      data       zspm3( 81)/       6.8679210d0/
      data       zppm3( 81)/       1.9694450d0/
      data      alppm3( 81)/       1.3409510d0/
      data      gsspm3( 81)/      10.4604120d0/
      data      gsppm3( 81)/      11.2238830d0/
      data      gpppm3( 81)/       4.9927850d0/
      data      gp2pm3( 81)/       8.9627270d0/
      data      hsppm3( 81)/       2.5304060d0/
      data    guesp1( 81,1)/      -1.3613990d0/
      data    guesp2( 81,1)/       3.5572260d0/
      data    guesp3( 81,1)/       1.0928020d0/
      data    guesp1( 81,2)/      -0.0454010d0/
      data    guesp2( 81,2)/       2.3069950d0/
      data    guesp3( 81,2)/       2.9650290d0/
!
!       Data for Element  82:         Lead
!
      data      usspm3( 82)/     -30.3227560d0/
      data      upppm3( 82)/     -24.4258340d0/
      data      betasp( 82)/      -6.1260240d0/
      data      betapp( 82)/      -1.3954300d0/
      data       zspm3( 82)/       3.1412890d0/
      data       zppm3( 82)/       1.8924180d0/
      data      alppm3( 82)/       1.6200450d0/
      data      gsspm3( 82)/       7.0119920d0/
      data      gsppm3( 82)/       6.7937820d0/
      data      gpppm3( 82)/       5.1837800d0/
      data      gp2pm3( 82)/       5.0456510d0/
      data      hsppm3( 82)/       1.5663020d0/
      data    guesp1( 82,1)/      -0.1225760d0/
      data    guesp2( 82,1)/       6.0030620d0/
      data    guesp3( 82,1)/       1.9015970d0/
      data    guesp1( 82,2)/      -0.0566480d0/
      data    guesp2( 82,2)/       4.7437050d0/
      data    guesp3( 82,2)/       2.8618790d0/
!
!       Data for Element  83:      Bismuth
!
      data      usspm3( 83)/     -33.4959380d0/
      data      upppm3( 83)/     -35.5210260d0/
      data      betasp( 83)/      -5.6072830d0/
      data      betapp( 83)/      -5.8001520d0/
      data       zspm3( 83)/       4.9164510d0/
      data       zppm3( 83)/       1.9349350d0/
      data      alppm3( 83)/       1.8574310d0/
      data      gsspm3( 83)/       4.9894800d0/
      data      gsppm3( 83)/       6.1033080d0/
      data      gpppm3( 83)/       8.6960070d0/
      data      gp2pm3( 83)/       8.3354470d0/
      data      hsppm3( 83)/       0.5991220d0/
      data    guesp1( 83,1)/       2.5816930d0/
      data    guesp2( 83,1)/       5.0940220d0/
      data    guesp3( 83,1)/       0.4997870d0/
      data    guesp1( 83,2)/       0.0603200d0/
      data    guesp2( 83,2)/       6.0015380d0/
      data    guesp3( 83,2)/       2.4279700d0/
!
!       Data for Element  87:     Francium
!
      data      alppm3( 87)/       3.0000000d0/
      data      gsspm3( 87)/      10.0000000d0/
!
!       Data for Element 102:  Capped bond
!
      data      usspm3(102)/     -11.9062760d0/
      data      betasp(102)/-9999999.0000000d0/
      data       zspm3(102)/       4.0000000d0/
      data       zppm3(102)/       0.3000000d0/
      data      alppm3(102)/       2.5441341d0/
      data      gsspm3(102)/      12.8480000d0/
      data      hsppm3(102)/       0.1000000d0/
!
!       Data for Element 104:    + Sparkle
!
      data      alppm3(104)/       1.5000000d0/
!
!       Data for Element 106:    - Sparkle
!
      data      alppm3(106)/       1.5000000d0/
      end module Parameters_for_PM3_C
