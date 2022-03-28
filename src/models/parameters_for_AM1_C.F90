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

      module Parameters_for_AM1_C
      double precision, dimension(107) :: ussam1, uppam1, uddam1, zsam1, zpam1, zdam1, &
        betasa, betapa, betada, alpam1, gssam1, gspam1, gppam1, gp2am1, hspam1, &
        polvolam1, zsnam1, zpnam1, zdnam1, f0sdam1, g2sdam1
      double precision, dimension(107,4) :: guesa1, guesa2, guesa3
!
!       Data for Element   1:     Hydrogen
!
      data      ussam1(  1)/     -11.3964270d0/
      data      betasa(  1)/      -6.1737870d0/
      data       zsam1(  1)/       1.1880780d0/
      data      alpam1(  1)/       2.8823240d0/
      data      gssam1(  1)/      12.8480000d0/
      data   polvolam1(  1)/       0.1719890d0/
      data    guesa1(  1,1)/       0.1227960d0/
      data    guesa2(  1,1)/       5.0000000d0/
      data    guesa3(  1,1)/       1.2000000d0/
      data    guesa1(  1,2)/       0.0050900d0/
      data    guesa2(  1,2)/       5.0000000d0/
      data    guesa3(  1,2)/       1.8000000d0/
      data    guesa1(  1,3)/      -0.0183360d0/
      data    guesa2(  1,3)/       2.0000000d0/
      data    guesa3(  1,3)/       2.1000000d0/
!
!       Data for Element   2:       Helium
!
      data      ussam1(  2)/     -35.2271054d0/
      data      uppam1(  2)/       9.9998070d0/
      data      betasa(  2)/     -17.5188861d0/
      data      betapa(  2)/     -26.0525879d0/
      data       zsam1(  2)/       2.1956103d0/
      data       zpam1(  2)/       6.9012486d0/
      data      alpam1(  2)/       5.7797256d0/
      data      gssam1(  2)/       8.2914834d0/
      data      gspam1(  2)/      11.0151986d0/
      data      gppam1(  2)/      19.3843171d0/
      data      gp2am1(  2)/      11.3517050d0/
      data      hspam1(  2)/       0.5320257d0/
      data    guesa1(  2,1)/       0.2060365d0/
      data    guesa2(  2,1)/       1.0637861d0/
      data    guesa3(  2,1)/       0.4997820d0/
!
!       Data for Element   3:      Lithium
!
      data      ussam1(  3)/      -4.9384384d0/
      data      uppam1(  3)/      -3.0959064d0/
      data      betasa(  3)/      -1.4598822d0/
      data      betapa(  3)/      -1.5278541d0/
      data       zsam1(  3)/       0.7973487d0/
      data       zpam1(  3)/       0.9045583d0/
      data      alpam1(  3)/       1.5522111d0/
      data      gssam1(  3)/       5.3999239d0/
      data      gspam1(  3)/       8.9521838d0/
      data      gppam1(  3)/       4.4594975d0/
      data      gp2am1(  3)/      10.0000285d0/
      data      hspam1(  3)/       0.3999650d0/
!
!       Data for Element   4:    Beryllium
!
      data      ussam1(  4)/     -17.1528584d0/
      data      uppam1(  4)/     -14.6296419d0/
      data      betasa(  4)/      -4.4963564d0/
      data      betapa(  4)/      -2.6466323d0/
      data       zsam1(  4)/       0.7425237d0/
      data       zpam1(  4)/       0.8080499d0/
      data      alpam1(  4)/       0.4979614d0/
      data      gssam1(  4)/       7.5260764d0/
      data      gspam1(  4)/       8.5518975d0/
      data      gppam1(  4)/       1.0247326d0/
      data      gp2am1(  4)/       9.3833365d0/
      data      hspam1(  4)/       0.3996483d0/
!
!       Data for Element   5:    Boron
!
      data      ussam1(  5)/     -34.492870d0/
      data      uppam1(  5)/     -22.631525d0/
      data      betasa(  5)/      -9.599114d0/
      data      betapa(  5)/      -6.273757d0/
      data       zsam1(  5)/       1.611709d0/
      data       zpam1(  5)/       1.555385d0/
      data      alpam1(  5)/       2.446909d0/
      data      gssam1(  5)/      10.590000d0/
      data      gspam1(  5)/       9.560000d0/
      data      gppam1(  5)/       8.860000d0/
      data      gp2am1(  5)/       7.860000d0/
      data      hspam1(  5)/       1.810000d0/
      data    guesa1(  5,1)/       0.182613d0/
      data    guesa2(  5,1)/       6.000000d0/
      data    guesa3(  5,1)/       0.727592d0/
      data    guesa1(  5,2)/       0.118587d0/
      data    guesa2(  5,2)/       6.000000d0/
      data    guesa3(  5,2)/       1.466639d0/
      data    guesa1(  5,3)/      -0.073280d0/
      data    guesa2(  5,3)/       5.000000d0/
      data    guesa3(  5,3)/       1.570975d0/
!
!       Data for Element   6:       Carbon
!
      data      ussam1(  6)/     -52.0286580d0/
      data      uppam1(  6)/     -39.6142390d0/
      data      betasa(  6)/     -15.7157830d0/
      data      betapa(  6)/      -7.7192830d0/
      data       zsam1(  6)/       1.8086650d0/
      data       zpam1(  6)/       1.6851160d0/
      data      alpam1(  6)/       2.6482740d0/
      data      gssam1(  6)/      12.2300000d0/
      data      gspam1(  6)/      11.4700000d0/
      data      gppam1(  6)/      11.0800000d0/
      data      gp2am1(  6)/       9.8400000d0/
      data      hspam1(  6)/       2.4300000d0/
      data   polvolam1(  6)/       0.8472950d0/
      data    guesa1(  6,1)/       0.0113550d0/
      data    guesa2(  6,1)/       5.0000000d0/
      data    guesa3(  6,1)/       1.6000000d0/
      data    guesa1(  6,2)/       0.0459240d0/
      data    guesa2(  6,2)/       5.0000000d0/
      data    guesa3(  6,2)/       1.8500000d0/
      data    guesa1(  6,3)/      -0.0200610d0/
      data    guesa2(  6,3)/       5.0000000d0/
      data    guesa3(  6,3)/       2.0500000d0/
      data    guesa1(  6,4)/      -0.0012600d0/
      data    guesa2(  6,4)/       5.0000000d0/
      data    guesa3(  6,4)/       2.6500000d0/
!
!       Data for Element   7:     Nitrogen
!
      data      ussam1(  7)/     -71.8600000d0/
      data      uppam1(  7)/     -57.1675810d0/
      data      betasa(  7)/     -20.2991100d0/
      data      betapa(  7)/     -18.2386660d0/
      data       zsam1(  7)/       2.3154100d0/
      data       zpam1(  7)/       2.1579400d0/
      data      alpam1(  7)/       2.9472860d0/
      data      gssam1(  7)/      13.5900000d0/
      data      gspam1(  7)/      12.6600000d0/
      data      gppam1(  7)/      12.9800000d0/
      data      gp2am1(  7)/      11.5900000d0/
      data      hspam1(  7)/       3.1400000d0/
      data   polvolam1(  7)/       0.6838130d0/
      data    guesa1(  7,1)/       0.0252510d0/
      data    guesa2(  7,1)/       5.0000000d0/
      data    guesa3(  7,1)/       1.5000000d0/
      data    guesa1(  7,2)/       0.0289530d0/
      data    guesa2(  7,2)/       5.0000000d0/
      data    guesa3(  7,2)/       2.1000000d0/
      data    guesa1(  7,3)/      -0.0058060d0/
      data    guesa2(  7,3)/       2.0000000d0/
      data    guesa3(  7,3)/       2.4000000d0/
!
!       Data for Element   8:       Oxygen
!
      data      ussam1(  8)/     -97.8300000d0/
      data      uppam1(  8)/     -78.2623800d0/
      data      betasa(  8)/     -29.2727730d0/
      data      betapa(  8)/     -29.2727730d0/
      data       zsam1(  8)/       3.1080320d0/
      data       zpam1(  8)/       2.5240390d0/
      data      alpam1(  8)/       4.4553710d0/
      data      gssam1(  8)/      15.4200000d0/
      data      gspam1(  8)/      14.4800000d0/
      data      gppam1(  8)/      14.5200000d0/
      data      gp2am1(  8)/      12.9800000d0/
      data      hspam1(  8)/       3.9400000d0/
      data   polvolam1(  8)/       0.3245890d0/
      data    guesa1(  8,1)/       0.2809620d0/
      data    guesa2(  8,1)/       5.0000000d0/
      data    guesa3(  8,1)/       0.8479180d0/
      data    guesa1(  8,2)/       0.0814300d0/
      data    guesa2(  8,2)/       7.0000000d0/
      data    guesa3(  8,2)/       1.4450710d0/
!
!       Data for Element   9:     Fluorine
!
      data      ussam1(  9)/    -136.1055790d0/
      data      uppam1(  9)/    -104.8898850d0/
      data      betasa(  9)/     -69.5902770d0/
      data      betapa(  9)/     -27.9223600d0/
      data       zsam1(  9)/       3.7700820d0/
      data       zpam1(  9)/       2.4946700d0/
      data      alpam1(  9)/       5.5178000d0/
      data      gssam1(  9)/      16.9200000d0/
      data      gspam1(  9)/      17.2500000d0/
      data      gppam1(  9)/      16.7100000d0/
      data      gp2am1(  9)/      14.9100000d0/
      data      hspam1(  9)/       4.8300000d0/
      data   polvolam1(  9)/       0.1286740d0/
      data    guesa1(  9,1)/       0.2420790d0/
      data    guesa2(  9,1)/       4.8000000d0/
      data    guesa3(  9,1)/       0.9300000d0/
      data    guesa1(  9,2)/       0.0036070d0/
      data    guesa2(  9,2)/       4.6000000d0/
      data    guesa3(  9,2)/       1.6600000d0/
!
!       Data for Element  10:         Neon
!
      data      ussam1( 10)/       9.6554213d0/
      data      uppam1( 10)/     -71.1385745d0/
      data      betasa( 10)/      -0.1512354d0/
      data      betapa( 10)/     -23.8476491d0/
      data       zsam1( 10)/       5.9983770d0/
      data       zpam1( 10)/       4.1699304d0/
      data      alpam1( 10)/       2.5799642d0/
      data      gssam1( 10)/       0.4984870d0/
      data      gspam1( 10)/      10.1589283d0/
      data      gppam1( 10)/      18.9450935d0/
      data      gp2am1( 10)/       8.6008669d0/
      data      hspam1( 10)/       0.3017990d0/
      data    guesa1( 10,1)/       0.2388679d0/
      data    guesa2( 10,1)/       8.8109263d0/
      data    guesa3( 10,1)/       1.1095395d0/
!
!       Data for Element  11:       Sodium
!
      data      ussam1( 11)/      -5.0711164d0/
      data      uppam1( 11)/      -2.9704119d0/
      data      betasa( 11)/      -1.1375097d0/
      data      betapa( 11)/      -2.1005594d0/
      data       zsam1( 11)/       0.7890090d0/
      data       zpam1( 11)/       1.1399864d0/
      data      alpam1( 11)/       6.0000025d0/
      data      gssam1( 11)/       6.4751360d0/
      data      gspam1( 11)/       5.4272286d0/
      data      gppam1( 11)/       9.5560913d0/
      data      gp2am1( 11)/       5.4192229d0/
      data      hspam1( 11)/       2.8677294d0/
      data    guesa1( 11,1)/       0.8576729d0/
      data    guesa2( 11,1)/       1.3171032d0/
      data    guesa3( 11,1)/       2.0000497d0/
!
!       Data for Element  12:    Magnesium
!
      data      ussam1( 12)/     -14.6688806d0/
      data      uppam1( 12)/     -11.8763861d0/
      data      betasa( 12)/      -1.1883090d0/
      data      betapa( 12)/      -5.2849791d0/
      data       zsam1( 12)/       1.0128928d0/
      data       zpam1( 12)/       1.1798191d0/
      data      alpam1( 12)/       5.8667260d0/
      data      gssam1( 12)/       6.6824002d0/
      data      gspam1( 12)/       7.1060848d0/
      data      gppam1( 12)/       9.3035830d0/
      data      gp2am1( 12)/       9.4772262d0/
      data      hspam1( 12)/       0.7866442d0/
      data    guesa1( 12,1)/       1.0000046d0/
      data    guesa2( 12,1)/       3.0000044d0/
      data    guesa3( 12,1)/       2.0000043d0/
!
!       Data for Element  13:    Aluminium
!
      data      ussam1( 13)/     -24.3535850d0/
      data      uppam1( 13)/     -18.3636450d0/
      data      betasa( 13)/      -3.8668220d0/
      data      betapa( 13)/      -2.3171460d0/
      data       zsam1( 13)/       1.5165930d0/
      data       zpam1( 13)/       1.3063470d0/
      data      alpam1( 13)/       1.9765860d0/
      data      gssam1( 13)/       8.0900000d0/
      data      gspam1( 13)/       6.6300000d0/
      data      gppam1( 13)/       5.9800000d0/
      data      gp2am1( 13)/       5.4000000d0/
      data      hspam1( 13)/       0.7000000d0/
      data    guesa1( 13,1)/       0.0900000d0/
      data    guesa2( 13,1)/      12.3924430d0/
      data    guesa3( 13,1)/       2.0503940d0/
!
!       Data for Element  14:      Silicon
!
      data      ussam1( 14)/     -33.9536220d0/
      data      uppam1( 14)/     -28.9347490d0/
      data      betasa( 14)/      -3.7848520d0/
      data      betapa( 14)/      -1.9681230d0/
      data       zsam1( 14)/       1.8306970d0/
      data       zpam1( 14)/       1.2849530d0/
      data      alpam1( 14)/       2.2578160d0/
      data      gssam1( 14)/       9.8200000d0/
      data      gspam1( 14)/       8.3600000d0/
      data      gppam1( 14)/       7.3100000d0/
      data      gp2am1( 14)/       6.5400000d0/
      data      hspam1( 14)/       1.3200000d0/
      data    guesa1( 14,1)/       0.2500000d0/
      data    guesa2( 14,1)/       9.0000000d0/
      data    guesa3( 14,1)/       0.9114530d0/
      data    guesa1( 14,2)/       0.0615130d0/
      data    guesa2( 14,2)/       5.0000000d0/
      data    guesa3( 14,2)/       1.9955690d0/
      data    guesa1( 14,3)/       0.0207890d0/
      data    guesa2( 14,3)/       5.0000000d0/
      data    guesa3( 14,3)/       2.9906100d0/
!
!       Data for Element  15:   Phosphorus
!
      data      ussam1( 15)/     -42.0298630d0/
      data      uppam1( 15)/     -34.0307090d0/
      data      betasa( 15)/      -6.3537640d0/
      data      betapa( 15)/      -6.5907090d0/
      data       zsam1( 15)/       1.9812800d0/
      data       zpam1( 15)/       1.8751500d0/
      data      alpam1( 15)/       2.4553220d0/
      data      gssam1( 15)/      11.5600050d0/
      data      gspam1( 15)/       5.2374490d0/
      data      gppam1( 15)/       7.8775890d0/
      data      gp2am1( 15)/       7.3076480d0/
      data      hspam1( 15)/       0.7792380d0/
      data    guesa1( 15,1)/      -0.0318270d0/
      data    guesa2( 15,1)/       6.0000000d0/
      data    guesa3( 15,1)/       1.4743230d0/
      data    guesa1( 15,2)/       0.0184700d0/
      data    guesa2( 15,2)/       7.0000000d0/
      data    guesa3( 15,2)/       1.7793540d0/
      data    guesa1( 15,3)/       0.0332900d0/
      data    guesa2( 15,3)/       9.0000000d0/
      data    guesa3( 15,3)/       3.0065760d0/
!
!       Data for Element  16:       Sulfur
!
      data      ussam1( 16)/     -56.6940560d0/
      data      uppam1( 16)/     -48.7170490d0/
      data      betasa( 16)/      -3.9205660d0/
      data      betapa( 16)/      -7.9052780d0/
      data       zsam1( 16)/       2.3665150d0/
      data       zpam1( 16)/       1.6672630d0/
      data      alpam1( 16)/       2.4616480d0/
      data      gssam1( 16)/      11.7863290d0/
      data      gspam1( 16)/       8.6631270d0/
      data      gppam1( 16)/      10.0393080d0/
      data      gp2am1( 16)/       7.7816880d0/
      data      hspam1( 16)/       2.5321370d0/
      data    guesa1( 16,1)/      -0.5091950d0/
      data    guesa2( 16,1)/       4.5936910d0/
      data    guesa3( 16,1)/       0.7706650d0/
      data    guesa1( 16,2)/      -0.0118630d0/
      data    guesa2( 16,2)/       5.8657310d0/
      data    guesa3( 16,2)/       1.5033130d0/
      data    guesa1( 16,3)/       0.0123340d0/
      data    guesa2( 16,3)/      13.5573360d0/
      data    guesa3( 16,3)/       2.0091730d0/
!
!       Data for Element  17:     Chlorine
!
      data      ussam1( 17)/    -111.6139480d0/
      data      uppam1( 17)/     -76.6401070d0/
      data      betasa( 17)/     -24.5946700d0/
      data      betapa( 17)/     -14.6372160d0/
      data       zsam1( 17)/       3.6313760d0/
      data       zpam1( 17)/       2.0767990d0/
      data      alpam1( 17)/       2.9193680d0/
      data      gssam1( 17)/      15.0300000d0/
      data      gspam1( 17)/      13.1600000d0/
      data      gppam1( 17)/      11.3000000d0/
      data      gp2am1( 17)/       9.9700000d0/
      data      hspam1( 17)/       2.4200000d0/
      data   polvolam1( 17)/       1.6635500d0/
      data    guesa1( 17,1)/       0.0942430d0/
      data    guesa2( 17,1)/       4.0000000d0/
      data    guesa3( 17,1)/       1.3000000d0/
      data    guesa1( 17,2)/       0.0271680d0/
      data    guesa2( 17,2)/       4.0000000d0/
      data    guesa3( 17,2)/       2.1000000d0/
!
!       Data for Element  18:        Argon
!
      data      ussam1( 18)/      14.2755779d0/
      data      uppam1( 18)/     -84.6763032d0/
      data      betasa( 18)/      -2.4380306d0/
      data      betapa( 18)/     -22.9801169d0/
      data       zsam1( 18)/       0.9714216d0/
      data       zpam1( 18)/       5.9236231d0/
      data      alpam1( 18)/       1.8714620d0/
      data      gssam1( 18)/       0.9106879d0/
      data      gspam1( 18)/      16.3549177d0/
      data      gppam1( 18)/      19.9993788d0/
      data      gp2am1( 18)/      12.9297003d0/
      data      hspam1( 18)/       0.8181126d0/
      data    guesa1( 18,1)/       0.3969461d0/
      data    guesa2( 18,1)/       3.8771587d0/
      data    guesa3( 18,1)/       1.0958585d0/
!
!       Data for Element  19:    Potassium
!
      data      ussam1( 19)/      -4.2628511d0/
      data      uppam1( 19)/      -2.6669543d0/
      data      betasa( 19)/      -0.2601130d0/
      data      betapa( 19)/      -1.6603661d0/
      data       zsam1( 19)/       1.2660244d0/
      data       zpam1( 19)/       0.9555939d0/
      data      alpam1( 19)/       5.8806897d0/
      data      gssam1( 19)/      10.0000250d0/
      data      gspam1( 19)/       4.8402183d0/
      data      gppam1( 19)/       2.9327889d0/
      data      gp2am1( 19)/       2.7340015d0/
      data      hspam1( 19)/       3.0014514d0/
!
!       Data for Element  20:      Calcium
!
      data      ussam1( 20)/     -12.3085333d0/
      data      uppam1( 20)/      -9.4760505d0/
      data      betasa( 20)/      -4.2657396d0/
      data      betapa( 20)/      -6.2934710d0/
      data       zsam1( 20)/       1.1767754d0/
      data       zpam1( 20)/       1.2738520d0/
      data      alpam1( 20)/       1.2137738d0/
      data      gssam1( 20)/       6.4320360d0/
      data      gspam1( 20)/       6.0410623d0/
      data      gppam1( 20)/       5.4337615d0/
      data      gp2am1( 20)/       5.3548430d0/
      data      hspam1( 20)/       0.6525142d0/
!
!       Data for Element  30:         Zinc
!
      data      ussam1( 30)/     -21.0400080d0/
      data      uppam1( 30)/     -17.6555740d0/
      data      betasa( 30)/      -1.9974290d0/
      data      betapa( 30)/      -4.7581190d0/
      data       zsam1( 30)/       1.9542990d0/
      data       zpam1( 30)/       1.3723650d0/
      data      alpam1( 30)/       1.4845630d0/
      data      gssam1( 30)/      11.8000000d0/
      data      gspam1( 30)/      11.1820180d0/
      data      gppam1( 30)/      13.3000000d0/
      data      gp2am1( 30)/      12.9305200d0/
      data      hspam1( 30)/       0.4846060d0/
!
!       Data for Element  31:      Gallium
!
      data      ussam1( 31)/     -29.7425311d0/
      data      uppam1( 31)/     -27.8983288d0/
      data      betasa( 31)/      -3.9999888d0/
      data      betapa( 31)/      -3.9993727d0/
      data       zsam1( 31)/       4.0002160d0/
      data       zpam1( 31)/       1.3540466d0/
      data      alpam1( 31)/       2.0566626d0/
      data      gssam1( 31)/       8.9143011d0/
      data      gspam1( 31)/      10.9447570d0/
      data      gppam1( 31)/       6.8859524d0/
      data      gp2am1( 31)/       8.4778025d0/
      data      hspam1( 31)/       0.5996245d0/
      data    guesa1( 31,1)/       0.0541706d0/
      data    guesa2( 31,1)/       0.9999031d0/
      data    guesa3( 31,1)/       2.0000642d0/
!
!       Data for Element  32:    Germanium
!
      data      ussam1( 32)/     -34.1838890d0/
      data      uppam1( 32)/     -28.6408110d0/
      data      betasa( 32)/      -4.3566070d0/
      data      betapa( 32)/      -0.9910910d0/
      data       zsam1( 32)/       1.2196310d0/
      data       zpam1( 32)/       1.9827940d0/
      data      alpam1( 32)/       2.1364050d0/
      data      gssam1( 32)/      10.1686050d0/
      data      gspam1( 32)/       8.1444730d0/
      data      gppam1( 32)/       6.6719020d0/
      data      gp2am1( 32)/       6.2697060d0/
      data      hspam1( 32)/       0.9370930d0/
!
!       Data for Element  33:      Arsenic
!
      data      ussam1( 33)/     -41.6817510d0/
      data      uppam1( 33)/     -33.4506152d0/
      data      betasa( 33)/      -5.6481504d0/
      data      betapa( 33)/      -4.9979109d0/
      data       zsam1( 33)/       2.2576897d0/
      data       zpam1( 33)/       1.7249710d0/
      data      alpam1( 33)/       2.2405380d0/
      data      gssam1( 33)/      11.0962258d0/
      data      gspam1( 33)/       4.9259328d0/
      data      gppam1( 33)/       7.8781648d0/
      data      gp2am1( 33)/       7.5961088d0/
      data      hspam1( 33)/       0.6246173d0/
      data    guesa1( 33,1)/      -0.0073614d0/
      data    guesa2( 33,1)/       4.9433993d0/
      data    guesa3( 33,1)/       1.4544264d0/
      data    guesa1( 33,2)/       0.0437629d0/
      data    guesa2( 33,2)/       3.1944613d0/
      data    guesa3( 33,2)/       2.0144939d0/
!
!       Data for Element  34:     Selenium
!
      data      ussam1( 34)/     -41.9984056d0/
      data      uppam1( 34)/     -32.8575485d0/
      data      betasa( 34)/      -3.1470826d0/
      data      betapa( 34)/      -6.1468406d0/
      data       zsam1( 34)/       2.6841570d0/
      data       zpam1( 34)/       2.0506164d0/
      data      alpam1( 34)/       2.6375694d0/
      data      gssam1( 34)/       6.7908891d0/
      data      gspam1( 34)/       6.4812786d0/
      data      gppam1( 34)/       6.4769273d0/
      data      gp2am1( 34)/       5.2796993d0/
      data      hspam1( 34)/       4.4548356d0/
      data    guesa1( 34,1)/       0.1116681d0/
      data    guesa2( 34,1)/       6.5086644d0/
      data    guesa3( 34,1)/       1.4981077d0/
      data    guesa1( 34,2)/       0.0396143d0/
      data    guesa2( 34,2)/       6.5241228d0/
      data    guesa3( 34,2)/       2.0751916d0/
!
!       Data for Element  35:      Bromine
!
      data      ussam1( 35)/    -104.6560630d0/
      data      uppam1( 35)/     -74.9300520d0/
      data      betasa( 35)/     -19.3998800d0/
      data      betapa( 35)/      -8.9571950d0/
      data       zsam1( 35)/       3.0641330d0/
      data       zpam1( 35)/       2.0383330d0/
      data      alpam1( 35)/       2.5765460d0/
      data      gssam1( 35)/      15.0364395d0/
      data      gspam1( 35)/      13.0346824d0/
      data      gppam1( 35)/      11.2763254d0/
      data      gp2am1( 35)/       9.8544255d0/
      data      hspam1( 35)/       2.4558683d0/
      data   polvolam1( 35)/       2.5941900d0/
      data    guesa1( 35,1)/       0.0666850d0/
      data    guesa2( 35,1)/       4.0000000d0/
      data    guesa3( 35,1)/       1.5000000d0/
      data    guesa1( 35,2)/       0.0255680d0/
      data    guesa2( 35,2)/       4.0000000d0/
      data    guesa3( 35,2)/       2.3000000d0/
!
!       Data for Element  36:      Krypton
!
      data      ussam1( 36)/       9.7991938d0/
      data      uppam1( 36)/     -73.8601784d0/
      data      betasa( 36)/      -1.2088255d0/
      data      betapa( 36)/      -8.1956064d0/
      data       zsam1( 36)/       3.5931632d0/
      data       zpam1( 36)/       2.0944633d0/
      data      alpam1( 36)/       1.7973083d0/
      data      gssam1( 36)/       0.7171517d0/
      data      gspam1( 36)/       9.5290790d0/
      data      gppam1( 36)/      19.9666844d0/
      data      gp2am1( 36)/      11.2405905d0/
      data      hspam1( 36)/       2.1732260d0/
      data    guesa1( 36,1)/       0.8882237d0/
      data    guesa2( 36,1)/       8.7917499d0/
      data    guesa3( 36,1)/       1.6319435d0/
!
!       Data for Element  37:     Rubidium
!
      data      ussam1( 37)/      -4.4990147d0/
      data      uppam1( 37)/      -2.9263643d0/
      data      betasa( 37)/      -1.9999892d0/
      data      betapa( 37)/      -4.4131246d0/
      data       zsam1( 37)/       4.0000187d0/
      data       zpam1( 37)/       1.0140619d0/
      data      alpam1( 37)/       1.1550020d0/
      data      gssam1( 37)/      18.7604025d0/
      data      gspam1( 37)/      18.0931959d0/
      data      gppam1( 37)/      10.8002500d0/
      data      gp2am1( 37)/       9.5613216d0/
      data      hspam1( 37)/       0.7084525d0/
      data    guesa1( 37,1)/       0.6444472d0/
      data    guesa2( 37,1)/       0.9994819d0/
      data    guesa3( 37,1)/       2.0004780d0/
!
!       Data for Element  38:    Strontium
!
      data      ussam1( 38)/     -10.9278146d0/
      data      uppam1( 38)/      -8.5185910d0/
      data      betasa( 38)/      -9.6008645d0/
      data      betapa( 38)/      -3.0661804d0/
      data       zsam1( 38)/       1.5236848d0/
      data       zpam1( 38)/       1.5723524d0/
      data      alpam1( 38)/       4.6716058d0/
      data      gssam1( 38)/       5.1033321d0/
      data      gspam1( 38)/       4.4927652d0/
      data      gppam1( 38)/       4.2101543d0/
      data      gp2am1( 38)/       4.3004995d0/
      data      hspam1( 38)/       0.7724969d0/
!
!       Data for Element  42:   Molybdenum
!
      data     ussam1( 42)/     -44.4880000d0/
      data     uppam1( 42)/     -20.2950000d0/
      data     uddam1( 42)/     -55.9520000d0/
      data     betasa( 42)/      -9.4140000d0/
      data     betapa( 42)/      -6.1800000d0/
      data     betada( 42)/     -15.4890000d0/
      data      zsam1( 42)/       1.9450000d0/
      data      zpam1( 42)/       1.4770000d0/
      data      zdam1( 42)/       2.4680000d0/
      data     zsnam1( 42)/       1.4240000d0/
      data     zpnam1( 42)/       1.2500000d0/
      data     zdnam1( 42)/       1.9470000d0/
      data    f0sdam1( 42)/       7.7050000d0/
      data    g2sdam1( 42)/       1.2000000d0/
      data     gssam1( 42)/       6.384309727d0/
!
!       Data for Element  49:       Indium
!
      data      ussam1( 49)/     -28.2223064d0/
      data      uppam1( 49)/     -18.3287837d0/
      data      betasa( 49)/      -6.1333658d0/
      data      betapa( 49)/      -0.9999602d0/
      data       zsam1( 49)/       1.8281576d0/
      data       zpam1( 49)/       1.4847500d0/
      data      alpam1( 49)/       1.8590637d0/
      data      gssam1( 49)/       9.3685202d0/
      data      gspam1( 49)/       6.6873024d0/
      data      gppam1( 49)/       5.9406805d0/
      data      gp2am1( 49)/       4.9356943d0/
      data      hspam1( 49)/       0.5998997d0/
      data    guesa1( 49,1)/       0.1182997d0/
      data    guesa2( 49,1)/       1.0033833d0/
      data    guesa3( 49,1)/       1.8646418d0/
!
!       Data for Element  50:          Tin
!
      data      ussam1( 50)/     -26.6529104d0/
      data      uppam1( 50)/     -12.7840857d0/
      data      betasa( 50)/      -1.9999126d0/
      data      betapa( 50)/      -2.1702085d0/
      data       zsam1( 50)/       1.6182807d0/
      data       zpam1( 50)/       1.5084984d0/
      data      alpam1( 50)/       1.6753624d0/
      data      gssam1( 50)/       7.0918140d0/
      data      gspam1( 50)/       2.9999326d0/
      data      gppam1( 50)/       5.3314764d0/
      data      gp2am1( 50)/       3.5204737d0/
      data      hspam1( 50)/       2.9523812d0/
      data    guesa1( 50,1)/      -0.3446636d0/
      data    guesa2( 50,1)/       1.9822018d0/
      data    guesa3( 50,1)/       1.3433163d0/
!
!       Data for Element  51:     Antimony
!
      data      ussam1( 51)/     -44.4381620d0/
      data      uppam1( 51)/     -32.3895140d0/
      data      betasa( 51)/      -7.3823300d0/
      data      betapa( 51)/      -3.6331190d0/
      data       zsam1( 51)/       2.2548230d0/
      data       zpam1( 51)/       2.2185920d0/
      data      alpam1( 51)/       2.2763310d0/
      data      gssam1( 51)/      11.4302510d0/
      data      gspam1( 51)/       5.7879220d0/
      data      gppam1( 51)/       6.4240940d0/
      data      gp2am1( 51)/       6.8491810d0/
      data      hspam1( 51)/       0.5883400d0/
      data    guesa1( 51,1)/      -0.5964470d0/
      data    guesa2( 51,1)/       6.0279500d0/
      data    guesa3( 51,1)/       1.7103670d0/
      data    guesa1( 51,2)/       0.8955130d0/
      data    guesa2( 51,2)/       3.0281090d0/
      data    guesa3( 51,2)/       1.5383180d0/
!
!       Data for Element  52:    Tellurium
!
      data      ussam1( 52)/     -39.2454230d0/
      data      uppam1( 52)/     -30.8515845d0/
      data      betasa( 52)/      -8.3897294d0/
      data      betapa( 52)/      -5.1065429d0/
      data       zsam1( 52)/       2.1321165d0/
      data       zpam1( 52)/       1.9712680d0/
      data      alpam1( 52)/       6.0171167d0/
      data      gssam1( 52)/       4.9925231d0/
      data      gspam1( 52)/       4.9721484d0/
      data      gppam1( 52)/       7.2097852d0/
      data      gp2am1( 52)/       5.6211521d0/
      data      hspam1( 52)/       4.0071821d0/
      data    guesa1( 52,1)/       0.4873378d0/
      data    guesa2( 52,1)/       6.0519413d0/
      data    guesa3( 52,1)/       1.3079857d0/
      data    guesa1( 52,2)/       0.1520464d0/
      data    guesa2( 52,2)/       3.8304067d0/
      data    guesa3( 52,2)/       2.0899707d0/
!
!       Data for Element  53:       Iodine
!
      data      ussam1( 53)/    -103.5896630d0/
      data      uppam1( 53)/     -74.4299970d0/
      data      betasa( 53)/      -8.4433270d0/
      data      betapa( 53)/      -6.3234050d0/
      data       zsam1( 53)/       2.1028580d0/
      data       zpam1( 53)/       2.1611530d0/
      data      alpam1( 53)/       2.2994240d0/
      data      gssam1( 53)/      15.0404486d0/
      data      gspam1( 53)/      13.0565580d0/
      data      gppam1( 53)/      11.1477837d0/
      data      gp2am1( 53)/       9.9140907d0/
      data      hspam1( 53)/       2.4563820d0/
      data   polvolam1( 53)/       4.5523800d0/
      data    guesa1( 53,1)/       0.0043610d0/
      data    guesa2( 53,1)/       2.3000000d0/
      data    guesa3( 53,1)/       1.8000000d0/
      data    guesa1( 53,2)/       0.0157060d0/
      data    guesa2( 53,2)/       3.0000000d0/
      data    guesa3( 53,2)/       2.2400000d0/
!
!       Data for Element  54:        Xenon
!
      data      ussam1( 54)/      29.9974173d0/
      data      uppam1( 54)/    -109.1315491d0/
      data      betasa( 54)/      -3.6467952d0/
      data      betapa( 54)/      -4.9334990d0/
      data       zsam1( 54)/       4.9675243d0/
      data       zpam1( 54)/       3.1432142d0/
      data      alpam1( 54)/       3.8260675d0/
      data      gssam1( 54)/       5.3202333d0/
      data      gspam1( 54)/       8.8602024d0/
      data      gppam1( 54)/      12.4851787d0/
      data      gp2am1( 54)/      19.1130026d0/
      data      hspam1( 54)/       2.4600467d0/
      data    guesa1( 54,1)/       1.0074929d0/
      data    guesa2( 54,1)/       2.7775892d0/
      data    guesa3( 54,1)/       0.5351145d0/
!
!       Data for Element  55:       Cesium
!
      data      ussam1( 55)/      -3.1358230d0/
      data      uppam1( 55)/      -1.6791847d0/
      data      betasa( 55)/      -4.4412054d0/
      data      betapa( 55)/      -4.3246899d0/
      data       zsam1( 55)/       5.7873708d0/
      data       zpam1( 55)/       1.0311693d0/
      data      alpam1( 55)/       0.5267821d0/
      data      gssam1( 55)/       3.8928349d0/
      data      gspam1( 55)/       2.9638098d0/
      data      gppam1( 55)/       5.6069289d0/
      data      gp2am1( 55)/       3.5192887d0/
      data      hspam1( 55)/       0.3994827d0/
      data    guesa1( 55,1)/      -1.0009018d0/
      data    guesa2( 55,1)/       2.4474604d0/
      data    guesa3( 55,1)/       0.6728225d0/
!
!       Data for Element  56:       Barium
!
      data      ussam1( 56)/     -10.1164434d0/
      data      uppam1( 56)/      -8.0393806d0/
      data      betasa( 56)/      -9.9997673d0/
      data      betapa( 56)/      -9.7724365d0/
      data       zsam1( 56)/       1.9136517d0/
      data       zpam1( 56)/       1.3948894d0/
      data      alpam1( 56)/       0.9963852d0/
      data      gssam1( 56)/       4.8572599d0/
      data      gspam1( 56)/       4.4042932d0/
      data      gppam1( 56)/       4.7218273d0/
      data      gp2am1( 56)/       4.8406105d0/
      data      hspam1( 56)/       0.5159824d0/
!
!       Data for Element  80:      Mercury
!
      data      ussam1( 80)/     -19.9415780d0/
      data      uppam1( 80)/     -11.1108700d0/
      data      betasa( 80)/      -0.9086570d0/
      data      betapa( 80)/      -4.9093840d0/
      data       zsam1( 80)/       2.0364130d0/
      data       zpam1( 80)/       1.9557660d0/
      data      alpam1( 80)/       1.4847340d0/
      data      gssam1( 80)/      10.8000000d0/
      data      gspam1( 80)/       9.3000000d0/
      data      gppam1( 80)/      14.3000000d0/
      data      gp2am1( 80)/      13.5000000d0/
      data      hspam1( 80)/       1.3000000d0/
!
!       Data for Element  81:     Thallium
!
      data      ussam1( 81)/     -29.8282621d0/
      data      uppam1( 81)/     -30.5358091d0/
      data      betasa( 81)/      -6.6096803d0/
      data      betapa( 81)/      -6.5157709d0/
      data       zsam1( 81)/       3.8077333d0/
      data       zpam1( 81)/       1.5511578d0/
      data      alpam1( 81)/       1.2571916d0/
      data      gssam1( 81)/       9.0641669d0/
      data      gspam1( 81)/      12.5941972d0/
      data      gppam1( 81)/      10.2189635d0/
      data      gp2am1( 81)/      13.0987769d0/
      data      hspam1( 81)/       0.5997632d0/
      data    guesa1( 81,1)/      -0.5293156d0/
      data    guesa2( 81,1)/       1.2083491d0/
      data    guesa3( 81,1)/       1.4794195d0/
!
!       Data for Element  82:         Lead
!
      data      ussam1( 82)/     -38.6798569d0/
      data      uppam1( 82)/     -26.4559953d0/
      data      betasa( 82)/      -6.5924919d0/
      data      betapa( 82)/      -1.3368867d0/
      data       zsam1( 82)/       2.4432161d0/
      data       zpam1( 82)/       1.5506706d0/
      data      alpam1( 82)/       1.6534073d0/
      data      gssam1( 82)/       8.6199280d0/
      data      gspam1( 82)/       7.6465534d0/
      data      gppam1( 82)/       6.7366709d0/
      data      gp2am1( 82)/       5.4967156d0/
      data      hspam1( 82)/       1.2176598d0/
      data    guesa1( 82,1)/      -0.3085992d0/
      data    guesa2( 82,1)/       3.0001372d0/
      data    guesa3( 82,1)/       1.6877190d0/
!
!       Data for Element  83:      Bismuth
!
      data      ussam1( 83)/     -42.0556490d0/
      data      uppam1( 83)/     -34.9221058d0/
      data      betasa( 83)/      -0.9993474d0/
      data      betapa( 83)/      -1.8948197d0/
      data       zsam1( 83)/       4.0007862d0/
      data       zpam1( 83)/       0.9547714d0/
      data      alpam1( 83)/       1.9060635d0/
      data      gssam1( 83)/      10.3839608d0/
      data      gspam1( 83)/       5.7403240d0/
      data      gppam1( 83)/      12.2196363d0/
      data      gp2am1( 83)/      11.2050063d0/
      data      hspam1( 83)/       5.0004083d0/
      data    guesa1( 83,1)/      -1.0004931d0/
      data    guesa2( 83,1)/       1.5860780d0/
      data    guesa3( 83,1)/       1.1085026d0/
!
!       Data for Element 102:  Capped bond
!
      data      ussam1(102)/     -11.9062760d0/
      data      betasa(102)/-9999999.0000000d0/
      data       zsam1(102)/       4.0000000d0/
      data       zpam1(102)/       0.3000000d0/
      data       zdam1(102)/       0.3000000d0/
      data      alpam1(102)/       2.5441341d0/
      data      gssam1(102)/      12.8480000d0/
      data      hspam1(102)/       0.1000000d0/
!
!       Data for Element 104:    + Sparkle
!
      data      alpam1(104)/       1.5000000d0/
!
!       Data for Element 106:    - Sparkle
!
      data      alpam1(106)/       1.5000000d0/
   end module Parameters_for_AM1_C
