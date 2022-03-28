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

      module Parameters_for_MNDO_C
      double precision, dimension(107) :: ussm, uppm, uddm, zsm, zpm, zdm, betasm, &
        betapm, betadm, alpm, gssm, gspm, gppm, gp2m, hspm, polvom, &
        zsnm, zpnm, zdnm, polvolm, f0sdm, g2sdm, pocm
      double precision, dimension(107,4) :: guesm1, guesm2, guesm3
!
!       Data for Element   1:     Hydrogen
!
      data      ussm(  1)/     -11.9062760d0/
      data    betasm(  1)/      -6.9890640d0/
      data       zsm(  1)/       1.3319670d0/
      data      alpm(  1)/       2.5441341d0/
      data      gssm(  1)/      12.8480000d0/
      data   polvolm(  1)/       0.1836440d0/
!
!       Data for Element   2:       Helium
!
      data      ussm(  2)/     -35.4931170d0/
      data      uppm(  2)/       9.9998433d0/
      data    betasm(  2)/     -17.5545020d0/
      data    betapm(  2)/     -26.0526064d0/
      data       zsm(  2)/       1.7710761d0/
      data       zpm(  2)/       6.9018258d0/
      data      alpm(  2)/       5.7813931d0/
      data      gssm(  2)/       8.0590135d0/
      data      gspm(  2)/      11.0152229d0/
      data      gppm(  2)/      19.3843339d0/
      data      gp2m(  2)/      11.3517224d0/
      data      hspm(  2)/       0.5321205d0/
      data guesm1(  2,1)/       0.1985976d0/
      data guesm2(  2,1)/       1.0009593d0/
      data guesm3(  2,1)/       0.5031509d0/
!
!       Data for Element   3:      Lithium
!
      data      ussm(  3)/      -4.8578570d0/
      data      uppm(  3)/      -2.0266084d0/
      data    betasm(  3)/      -0.1904126d0/
      data    betapm(  3)/      -1.6081355d0/
      data       zsm(  3)/       0.4296141d0/
      data       zpm(  3)/       0.7554884d0/
      data      alpm(  3)/       1.2083987d0/
      data      gssm(  3)/       7.5947069d0/
      data      gspm(  3)/       6.7259856d0/
      data      gppm(  3)/       8.6596829d0/
      data      gp2m(  3)/       3.8714751d0/
      data      hspm(  3)/       5.0003381d0/
!
!       Data for Element   4:    Beryllium
!
      data      ussm(  4)/     -16.6023780d0/
      data      uppm(  4)/     -10.7037710d0/
      data    betasm(  4)/      -4.0170960d0/
      data    betapm(  4)/      -4.0170960d0/
      data       zsm(  4)/       1.0042100d0/
      data       zpm(  4)/       1.0042100d0/
      data      alpm(  4)/       1.6694340d0/
      data      gssm(  4)/       9.0000000d0/
      data      gspm(  4)/       7.4300000d0/
      data      gppm(  4)/       6.9700000d0/
      data      gp2m(  4)/       6.2200000d0/
      data      hspm(  4)/       1.2800000d0/
!
!       Data for Element   5:        Boron
!
      data      ussm(  5)/     -34.5471300d0/
      data      uppm(  5)/     -23.1216900d0/
      data    betasm(  5)/      -8.2520540d0/
      data    betapm(  5)/      -8.2520540d0/
      data       zsm(  5)/       1.5068010d0/
      data       zpm(  5)/       1.5068010d0/
      data      alpm(  5)/       2.1349930d0/
      data      gssm(  5)/      10.5900000d0/
      data      gspm(  5)/       9.5600000d0/
      data      gppm(  5)/       8.8600000d0/
      data      gp2m(  5)/       7.8600000d0/
      data      hspm(  5)/       1.8100000d0/
!
!       Data for Element   6:       Carbon
!
      data      ussm(  6)/     -52.2797450d0/
      data      uppm(  6)/     -39.2055580d0/
      data    betasm(  6)/     -18.9850440d0/
      data    betapm(  6)/      -7.9341220d0/
      data       zsm(  6)/       1.7875370d0/
      data       zpm(  6)/       1.7875370d0/
      data      alpm(  6)/       2.5463800d0/
      data      gssm(  6)/      12.2300000d0/
      data      gspm(  6)/      11.4700000d0/
      data      gppm(  6)/      11.0800000d0/
      data      gp2m(  6)/       9.8400000d0/
      data      hspm(  6)/       2.4300000d0/
      data   polvolm(  6)/       0.8215720d0/
!
!       Data for Element   7:     Nitrogen
!
      data      ussm(  7)/     -71.9321220d0/
      data      uppm(  7)/     -57.1723190d0/
      data    betasm(  7)/     -20.4957580d0/
      data    betapm(  7)/     -20.4957580d0/
      data       zsm(  7)/       2.2556140d0/
      data       zpm(  7)/       2.2556140d0/
      data      alpm(  7)/       2.8613420d0/
      data      gssm(  7)/      13.5900000d0/
      data      gspm(  7)/      12.6600000d0/
      data      gppm(  7)/      12.9800000d0/
      data      gp2m(  7)/      11.5900000d0/
      data      hspm(  7)/       3.1400000d0/
      data   polvolm(  7)/       0.7031620d0/
!
!       Data for Element   8:       Oxygen
!
      data      ussm(  8)/     -99.6443090d0/
      data      uppm(  8)/     -77.7974720d0/
      data    betasm(  8)/     -32.6880820d0/
      data    betapm(  8)/     -32.6880820d0/
      data       zsm(  8)/       2.6999050d0/
      data       zpm(  8)/       2.6999050d0/
      data      alpm(  8)/       3.1606040d0/
      data      gssm(  8)/      15.4200000d0/
      data      gspm(  8)/      14.4800000d0/
      data      gppm(  8)/      14.5200000d0/
      data      gp2m(  8)/      12.9800000d0/
      data      hspm(  8)/       3.9400000d0/
      data   polvolm(  8)/       0.3615210d0/
!
!       Data for Element   9:     Fluorine
!
      data      ussm(  9)/    -131.0715480d0/
      data      uppm(  9)/    -105.7821370d0/
      data    betasm(  9)/     -48.2904660d0/
      data    betapm(  9)/     -36.5085400d0/
      data       zsm(  9)/       2.8484870d0/
      data       zpm(  9)/       2.8484870d0/
      data      alpm(  9)/       3.4196606d0/
      data      gssm(  9)/      16.9200000d0/
      data      gspm(  9)/      17.2500000d0/
      data      gppm(  9)/      16.7100000d0/
      data      gp2m(  9)/      14.9100000d0/
      data      hspm(  9)/       4.8300000d0/
      data   polvolm(  9)/       0.2000100d0/
!
!       Data for Element  10:         Neon
!
      data      ussm( 10)/       9.6554736d0/
      data      uppm( 10)/     -71.1466474d0/
      data    betasm( 10)/      -0.1512657d0/
      data    betapm( 10)/     -23.8470187d0/
      data       zsm( 10)/       5.9998745d0/
      data       zpm( 10)/       4.1752600d0/
      data      alpm( 10)/       2.5377255d0/
      data      gssm( 10)/       0.4989374d0/
      data      gspm( 10)/      10.1592044d0/
      data      gppm( 10)/      18.9449976d0/
      data      gp2m( 10)/       8.5648589d0/
      data      hspm( 10)/       0.2999991d0/
      data guesm1( 10,1)/       0.2352109d0/
      data guesm2( 10,1)/       8.8116607d0/
      data guesm3( 10,1)/       1.0927884d0/
!
!       Data for Element  11:       Sodium
!
      data      ussm( 11)/      -5.1235942d0/
      data      uppm( 11)/      -3.0124713d0/
      data    betasm( 11)/      -1.4916657d0/
      data    betapm( 11)/      -0.2208233d0/
      data       zsm( 11)/       0.8213124d0/
      data       zpm( 11)/       1.0303270d0/
      data      alpm( 11)/       5.9940638d0/
      data      gssm( 11)/       6.9934520d0/
      data      gspm( 11)/       5.4380530d0/
      data      gppm( 11)/       6.9285408d0/
      data      gp2m( 11)/       2.4299952d0/
      data      hspm( 11)/       3.1374231d0/
!
!       Data for Element  12:    Magnesium
!
      data      ussm( 12)/     -14.8057921d0/
      data      uppm( 12)/     -12.7545494d0/
      data    betasm( 12)/      -0.0999626d0/
      data    betapm( 12)/      -7.6952692d0/
      data       zsm( 12)/       0.9394811d0/
      data       zpm( 12)/       1.3103428d0/
      data      alpm( 12)/       0.5354201d0/
      data      gssm( 12)/       6.9632774d0/
      data      gspm( 12)/       8.2410934d0/
      data      gppm( 12)/      10.0000349d0/
      data      gp2m( 12)/       9.0593304d0/
      data      hspm( 12)/       0.7165063d0/
!
!       Data for Element  13:    Aluminium
!
      data      ussm( 13)/     -23.8070970d0/
      data      uppm( 13)/     -17.5198780d0/
      data    betasm( 13)/      -2.6702840d0/
      data    betapm( 13)/      -2.6702840d0/
      data       zsm( 13)/       1.4441610d0/
      data       zpm( 13)/       1.4441610d0/
      data      alpm( 13)/       1.8688394d0/
      data      gssm( 13)/       8.0900000d0/
      data      gspm( 13)/       6.6300000d0/
      data      gppm( 13)/       5.9800000d0/
      data      gp2m( 13)/       5.4000000d0/
      data      hspm( 13)/       0.7000000d0/
!
!       Data for Element  14:      Silicon
!
      data      ussm( 14)/     -37.0375330d0/
      data      uppm( 14)/     -27.7696780d0/
      data    betasm( 14)/      -9.0868040d0/
      data    betapm( 14)/      -1.0758270d0/
      data       zsm( 14)/       1.3159860d0/
      data       zpm( 14)/       1.7099430d0/
      data      alpm( 14)/       2.2053160d0/
      data      gssm( 14)/       9.8200000d0/
      data      gspm( 14)/       8.3600000d0/
      data      gppm( 14)/       7.3100000d0/
      data      gp2m( 14)/       6.5400000d0/
      data      hspm( 14)/       1.3200000d0/
!
!       Data for Element  15:   Phosphorus
!
      data      ussm( 15)/     -56.1433600d0/
      data      uppm( 15)/     -42.8510800d0/
      data    betasm( 15)/      -6.7916000d0/
      data    betapm( 15)/      -6.7916000d0/
      data       zsm( 15)/       2.1087200d0/
      data       zpm( 15)/       1.7858100d0/
      data      alpm( 15)/       2.4152800d0/
      data      gssm( 15)/      11.5600000d0/
      data      gspm( 15)/      10.0800000d0/
      data      gppm( 15)/       8.6400000d0/
      data      gp2m( 15)/       7.6800000d0/
      data      hspm( 15)/       1.9200000d0/
!
!       Data for Element  16:       Sulfur
!
      data      ussm( 16)/     -72.2422810d0/
      data      uppm( 16)/     -56.9732070d0/
      data    betasm( 16)/     -10.7616700d0/
      data    betapm( 16)/     -10.1084330d0/
      data       zsm( 16)/       2.3129620d0/
      data       zpm( 16)/       2.0091460d0/
      data      alpm( 16)/       2.4780260d0/
      data      gssm( 16)/      12.8800000d0/
      data      gspm( 16)/      11.2600000d0/
      data      gppm( 16)/       9.9000000d0/
      data      gp2m( 16)/       8.8300000d0/
      data      hspm( 16)/       2.2600000d0/
!
!       Data for Element  17:     Chlorine
!
      data      ussm( 17)/    -100.2271660d0/
      data      uppm( 17)/     -77.3786670d0/
      data    betasm( 17)/     -14.2623200d0/
      data    betapm( 17)/     -14.2623200d0/
      data       zsm( 17)/       3.7846450d0/
      data       zpm( 17)/       2.0362630d0/
      data      alpm( 17)/       2.5422010d0/
      data      gssm( 17)/      15.0300000d0/
      data      gspm( 17)/      13.1600000d0/
      data      gppm( 17)/      11.3000000d0/
      data      gp2m( 17)/       9.9700000d0/
      data      hspm( 17)/       2.4200000d0/
      data   polvolm( 17)/       1.6618700d0/
!
!       Data for Element  18:        Argon
!
      data      ussm( 18)/       4.2759913d0/
      data      uppm( 18)/     -75.1004104d0/
      data    betasm( 18)/      -2.4320045d0/
      data    betapm( 18)/     -22.9856937d0/
      data       zsm( 18)/       0.9821697d0/
      data       zpm( 18)/       5.9997150d0/
      data      alpm( 18)/       2.2113324d0/
      data      gssm( 18)/       3.2377414d0/
      data      gspm( 18)/      16.3573969d0/
      data      gppm( 18)/      19.9984646d0/
      data      gp2m( 18)/      10.9343238d0/
      data      hspm( 18)/       0.8169151d0/
      data guesm1( 18,1)/      -0.7376503d0/
      data guesm2( 18,1)/       3.8882543d0/
      data guesm3( 18,1)/       0.7150665d0/
!
!       Data for Element  19:    Potassium
!
      data      ussm( 19)/      -3.6401731d0/
      data      uppm( 19)/      -2.0249253d0/
      data    betasm( 19)/      -0.1361851d0/
      data    betapm( 19)/      -2.8142350d0/
      data       zsm( 19)/       0.7276039d0/
      data       zpm( 19)/       0.9871174d0/
      data      alpm( 19)/       0.5616422d0/
      data      gssm( 19)/       3.7939792d0/
      data      gspm( 19)/       6.4170233d0/
      data      gppm( 19)/       5.0972823d0/
      data      gp2m( 19)/       2.1945567d0/
      data      hspm( 19)/       1.5788130d0/
!
!       Data for Element  20:      Calcium
!
      data      ussm( 20)/     -12.3919094d0/
      data      uppm( 20)/      -9.9348289d0/
      data    betasm( 20)/      -8.6404687d0/
      data    betapm( 20)/      -9.9515712d0/
      data       zsm( 20)/       1.0034161d0/
      data       zpm( 20)/       1.3102564d0/
      data      alpm( 20)/       0.4999997d0/
      data      gssm( 20)/       6.5321649d0/
      data      gspm( 20)/       6.5424442d0/
      data      gppm( 20)/       6.4627059d0/
      data      gp2m( 20)/       6.3842472d0/
      data      hspm( 20)/       0.5789676d0/
!
!       Data for Element  21:     Scandium
!
      data      ussm( 21)/     -15.0449253d0/
      data      uppm( 21)/     -13.3798697d0/
      data      uddm( 21)/     -22.0180116d0/
      data    betasm( 21)/      -0.1667211d0/
      data    betapm( 21)/      -1.1476882d0/
      data    betadm( 21)/      -0.9999564d0/
      data       zsm( 21)/       1.3951231d0/
      data       zpm( 21)/       5.0160943d0/
      data       zdm( 21)/       0.9264186d0/
      data      zsnm( 21)/       0.5018770d0/
      data      zpnm( 21)/       1.1053051d0/
      data      zdnm( 21)/       1.6874974d0/
      data      alpm( 21)/       0.5000735d0/
      data     f0sdm( 21)/       6.9604525d0/
      data     g2sdm( 21)/       1.1643214d0/
      data      pocm( 21)/       1.6013033d0/
!
!       Data for Element  22:     Titanium
!
      data      ussm( 22)/     -38.4842092d0/
      data      uppm( 22)/     -29.7533675d0/
      data      uddm( 22)/     -43.1364568d0/
      data    betasm( 22)/      -4.6294107d0/
      data    betapm( 22)/     -12.1879220d0/
      data    betadm( 22)/     -12.0696917d0/
      data       zsm( 22)/       0.8961552d0/
      data       zpm( 22)/       0.9676159d0/
      data       zdm( 22)/       1.8698884d0/
      data      zsnm( 22)/       1.7784012d0/
      data      zpnm( 22)/       1.5647427d0/
      data      zdnm( 22)/       1.9578396d0/
      data      alpm( 22)/       0.4999658d0/
!
!       Data for Element  23:     Vanadium
!
      data      ussm( 23)/     -36.6155632d0/
      data      uppm( 23)/     -13.7450846d0/
      data      uddm( 23)/     -49.3642151d0/
      data    betasm( 23)/      -2.7525182d0/
      data    betapm( 23)/      -3.8379043d0/
      data    betadm( 23)/      -4.8620752d0/
      data       zsm( 23)/       1.2873544d0/
      data       zpm( 23)/       1.1744379d0/
      data       zdm( 23)/       2.0150220d0/
      data      zsnm( 23)/       1.1517166d0/
      data      zpnm( 23)/       0.6496939d0/
      data      zdnm( 23)/       1.7060710d0/
      data      alpm( 23)/       2.6556548d0/
      data     f0sdm( 23)/       7.6311252d0/
      data     g2sdm( 23)/       1.3877486d0/
      data      pocm( 23)/       1.7704111d0/
!
!       Data for Element  24:     Chromium
!
      data      ussm( 24)/     -43.9655238d0/
      data      uppm( 24)/     -20.2837004d0/
      data      uddm( 24)/     -55.6262326d0/
      data    betasm( 24)/      -9.1351851d0/
      data    betapm( 24)/      -5.3477764d0/
      data    betadm( 24)/     -15.6910401d0/
      data       zsm( 24)/       2.1495003d0/
      data       zpm( 24)/       1.3131074d0/
      data       zdm( 24)/       2.3289346d0/
      data      zsnm( 24)/       0.9189804d0/
      data      zpnm( 24)/       0.6655499d0/
      data      zdnm( 24)/       1.5044456d0/
      data      alpm( 24)/       2.4305854d0/
      data     f0sdm( 24)/       7.9331069d0/
      data     g2sdm( 24)/       1.0113010d0/
      data      pocm( 24)/       1.6912431d0/
!
!       Data for Element  26:         Iron
!
      data      ussm( 26)/     -67.9422758d0/
      data      uppm( 26)/     -58.5634874d0/
      data      uddm( 26)/     -95.6743362d0/
      data    betasm( 26)/      -0.3135717d0/
      data    betapm( 26)/      -0.5122152d0/
      data    betadm( 26)/     -11.9511998d0/
      data       zsm( 26)/       1.4536275d0/
      data       zpm( 26)/       0.8933716d0/
      data       zdm( 26)/       1.8691105d0/
      data      zsnm( 26)/       1.2456921d0/
      data      zpnm( 26)/       1.8214074d0/
      data      zdnm( 26)/       1.9356440d0/
      data      alpm( 26)/       2.7722184d0/
      data     f0sdm( 26)/       8.9179096d0/
      data     g2sdm( 26)/       1.1871586d0/
      data      pocm( 26)/       1.0444762d0/
!
!       Data for Element  27:       Cobalt
!
      data      ussm( 27)/     -55.6182570d0/
      data      uppm( 27)/     -30.8068730d0/
      data      uddm( 27)/     -88.9567800d0/
      data    betasm( 27)/      -5.4461890d0/
      data    betapm( 27)/      -4.1051470d0/
      data    betadm( 27)/     -17.2059100d0/
      data       zsm( 27)/       0.5997500d0/
      data       zpm( 27)/       0.6073140d0/
      data       zdm( 27)/       1.8567970d0/
      data      zsnm( 27)/       0.8172170d0/
      data      zpnm( 27)/       2.0251460d0/
      data      zdnm( 27)/       1.6008320d0/
      data      alpm( 27)/       5.7540470d0/
      data     f0sdm( 27)/       6.7895950d0/
      data     g2sdm( 27)/       1.3490200d0/
!
!       Data for Element  28:       Nickel
!
      data      ussm( 28)/     -82.4828721d0/
      data      uppm( 28)/      10.0000008d0/
      data      uddm( 28)/     -68.0786730d0/
      data    betasm( 28)/     -12.1996627d0/
      data    betapm( 28)/      -6.8601153d0/
      data    betadm( 28)/     -26.8807821d0/
      data       zsm( 28)/       0.7735888d0/
      data       zpm( 28)/       6.0000132d0/
      data       zdm( 28)/       2.7857108d0/
      data      zsnm( 28)/       2.0155490d0/
      data      zpnm( 28)/       4.0000022d0/
      data      zdnm( 28)/       0.9297882d0/
      data      alpm( 28)/       3.8709519d0/
      data     f0sdm( 28)/       8.0280725d0/
      data     g2sdm( 28)/       2.4681665d0/
!
!       Data for Element  29:       Copper
!
      data      ussm( 29)/    -188.3951037d0/
      data      uppm( 29)/     -75.7592933d0/
      data      uddm( 29)/    -173.9068111d0/
      data    betasm( 29)/      -9.8202653d0/
      data    betapm( 29)/      -8.3686593d0/
      data    betadm( 29)/     -25.6520270d0/
      data       zsm( 29)/       3.3957872d0/
      data       zpm( 29)/       1.7861780d0/
      data       zdm( 29)/       3.3573266d0/
      data      zsnm( 29)/       3.9457339d0/
      data      zpnm( 29)/       1.1050853d0/
      data      zdnm( 29)/       2.4262504d0/
      data      alpm( 29)/       2.7857955d0/
!
!       Data for Element  30:         Zinc
!
      data      ussm( 30)/     -20.8397160d0/
      data      uppm( 30)/     -19.6252240d0/
      data    betasm( 30)/      -1.0000000d0/
      data    betapm( 30)/      -2.0000000d0/
      data       zsm( 30)/       2.0473590d0/
      data       zpm( 30)/       1.4609460d0/
      data      alpm( 30)/       1.5064570d0/
      data      gssm( 30)/      11.8000000d0/
      data      gspm( 30)/      11.1820180d0/
      data      gppm( 30)/      13.3000000d0/
      data      gp2m( 30)/      12.9305200d0/
      data      hspm( 30)/       0.4846060d0/
!
!       Data for Element  31:      Gallium
!
      data      ussm( 31)/     -28.3044924d0/
      data      uppm( 31)/     -27.2063910d0/
      data    betasm( 31)/      -3.9987435d0/
      data    betapm( 31)/      -4.3146711d0/
      data       zsm( 31)/       0.6986316d0/
      data       zpm( 31)/       1.8386933d0/
      data      alpm( 31)/       2.7577991d0/
      data      gssm( 31)/       7.5468114d0/
      data      gspm( 31)/      10.4697612d0/
      data      gppm( 31)/       8.4599454d0/
      data      gp2m( 31)/      10.4251148d0/
      data      hspm( 31)/       1.0628013d0/
      data guesm1( 31,1)/       0.6265886d0/
      data guesm2( 31,1)/       3.0001279d0/
      data guesm3( 31,1)/       1.2564374d0/
!
!       Data for Element  32:    Germanium
!
      data      ussm( 32)/     -33.9493670d0/
      data      uppm( 32)/     -27.4251050d0/
      data    betasm( 32)/      -4.5164790d0/
      data    betapm( 32)/      -1.7555170d0/
      data       zsm( 32)/       1.2931800d0/
      data       zpm( 32)/       2.0205640d0/
      data      alpm( 32)/       1.9784980d0/
      data      gssm( 32)/       9.8000000d0/
      data      gspm( 32)/       8.3000000d0/
      data      gppm( 32)/       7.3000000d0/
      data      gp2m( 32)/       6.5000000d0/
      data      hspm( 32)/       1.3000000d0/
!
!       Data for Element  33:      Arsenic
!
      data      ussm( 33)/     -38.6240790d0/
      data      uppm( 33)/     -33.4995395d0/
      data    betasm( 33)/      -3.9998231d0/
      data    betapm( 33)/      -4.9056176d0/
      data       zsm( 33)/       2.5614338d0/
      data       zpm( 33)/       1.6117315d0/
      data      alpm( 33)/       1.9381219d0/
      data      gssm( 33)/       6.7464011d0/
      data      gspm( 33)/       5.6174985d0/
      data      gppm( 33)/       6.9333925d0/
      data      gp2m( 33)/       6.3054798d0/
      data      hspm( 33)/       0.5994346d0/
      data guesm1( 33,1)/      -0.2025391d0/
      data guesm2( 33,1)/       3.0002200d0/
      data guesm3( 33,1)/       1.3902090d0/
!
!       Data for Element  34:     Selenium
!
      data      ussm( 34)/     -49.8117347d0/
      data      uppm( 34)/     -38.0475911d0/
      data    betasm( 34)/     -12.4685178d0/
      data    betapm( 34)/      -5.1744376d0/
      data       zsm( 34)/       0.7242956d0/
      data       zpm( 34)/       1.9267288d0/
      data      alpm( 34)/       2.3513155d0/
      data      gssm( 34)/      10.3549483d0/
      data      gspm( 34)/       5.2801360d0/
      data      gppm( 34)/       7.3611317d0/
      data      gp2m( 34)/       6.1897284d0/
      data      hspm( 34)/       0.5996560d0/
      data guesm1( 34,1)/      -1.0001121d0/
      data guesm2( 34,1)/       2.7043270d0/
      data guesm3( 34,1)/       0.4993772d0/
!
!       Data for Element  35:      Bromine
!
      data      ussm( 35)/     -99.9864405d0/
      data      uppm( 35)/     -75.6713075d0/
      data    betasm( 35)/      -8.9171070d0/
      data    betapm( 35)/      -9.9437400d0/
      data       zsm( 35)/       3.8543019d0/
      data       zpm( 35)/       2.1992091d0/
      data      alpm( 35)/       2.4457051d0/
      data      gssm( 35)/      15.0364395d0/
      data      gspm( 35)/      13.0346824d0/
      data      gppm( 35)/      11.2763254d0/
      data      gp2m( 35)/       9.8544255d0/
      data      hspm( 35)/       2.4558683d0/
      data   polvolm( 35)/       2.6532400d0/
!
!       Data for Element  36:      Krypton
!
      data      ussm( 36)/       9.7986937d0/
      data      uppm( 36)/     -73.8595661d0/
      data    betasm( 36)/      -1.2094707d0/
      data    betapm( 36)/      -8.4328271d0/
      data       zsm( 36)/       3.5608622d0/
      data       zpm( 36)/       1.9832062d0/
      data      alpm( 36)/       1.7025435d0/
      data      gssm( 36)/       0.4994169d0/
      data      gspm( 36)/       9.5261258d0/
      data      gppm( 36)/      19.9987286d0/
      data      gp2m( 36)/      11.1986052d0/
      data      hspm( 36)/       2.1818338d0/
      data guesm1( 36,1)/       0.7804390d0/
      data guesm2( 36,1)/       8.8371397d0/
      data guesm3( 36,1)/       0.5000419d0/
!
!       Data for Element  37:     Rubidium
!
      data      ussm( 37)/      -4.3098071d0/
      data      uppm( 37)/      -2.7381921d0/
      data    betasm( 37)/      -2.2062173d0/
      data    betapm( 37)/      -6.2176392d0/
      data       zsm( 37)/       4.0001632d0/
      data       zpm( 37)/       0.9187408d0/
      data      alpm( 37)/       0.9976197d0/
      data      gssm( 37)/      10.7409462d0/
      data      gspm( 37)/      11.4853623d0/
      data      gppm( 37)/       8.9878480d0/
      data      gp2m( 37)/       7.7258368d0/
      data      hspm( 37)/       0.1999294d0/
!
!       Data for Element  38:    Strontium
!
      data      ussm( 38)/     -10.8451287d0/
      data      uppm( 38)/      -8.3129821d0/
      data    betasm( 38)/      -9.9683427d0/
      data    betapm( 38)/      -9.9946390d0/
      data       zsm( 38)/       1.3729266d0/
      data       zpm( 38)/       1.1118128d0/
      data      alpm( 38)/       0.5082703d0/
      data      gssm( 38)/       4.9305520d0/
      data      gspm( 38)/       4.4249843d0/
      data      gppm( 38)/       4.0178045d0/
      data      gp2m( 38)/       4.0335640d0/
      data      hspm( 38)/       0.6272993d0/
!
!       Data for Element  40:    Zirconium
!
      data      ussm( 40)/     -37.7166490d0/
      data      uppm( 40)/     -29.3576225d0/
      data      uddm( 40)/     -42.6294048d0/
      data    betasm( 40)/      -3.6662593d0/
      data    betapm( 40)/     -10.8653809d0/
      data    betadm( 40)/     -11.1948555d0/
      data       zsm( 40)/       1.5386288d0/
      data       zpm( 40)/       1.1472515d0/
      data       zdm( 40)/       1.8744783d0/
      data      zsnm( 40)/       1.8931992d0/
      data      zpnm( 40)/       1.4185044d0/
      data      zdnm( 40)/       2.6877825d0/
      data      alpm( 40)/       1.0006438d0/
!
!       Data for Element  42:   Molybdenum
!
      data      ussm( 42)/     -44.4974811d0/
      data      uppm( 42)/     -19.4687724d0/
      data      uddm( 42)/     -55.4804886d0/
      data    betasm( 42)/      -8.9510879d0/
      data    betapm( 42)/      -5.5789800d0/
      data    betadm( 42)/     -15.2348605d0/
      data       zsm( 42)/       2.0001083d0/
      data       zpm( 42)/       1.4112837d0/
      data       zdm( 42)/       2.1944707d0/
      data      zsnm( 42)/       1.7997044d0/
      data      zpnm( 42)/       1.7801669d0/
      data      zdnm( 42)/       1.9231611d0/
      data      alpm( 42)/       0.9450106d0/
      data     f0sdm( 42)/       7.6696933d0/
      data     g2sdm( 42)/       0.9999819d0/
      data      pocm( 42)/       1.3526363d0/
!
!       Data for Element  46:    Palladium
!
      data      ussm( 46)/     -83.7375043d0/
      data      uppm( 46)/     -37.5463239d0/
      data      uddm( 46)/    -114.7045474d0/
      data    betasm( 46)/      -1.6390668d0/
      data    betapm( 46)/      -0.4999895d0/
      data    betadm( 46)/     -14.1763613d0/
      data       zsm( 46)/       1.6942397d0/
      data       zpm( 46)/       6.0000131d0/
      data       zdm( 46)/       2.2314824d0/
      data      zsnm( 46)/       2.3416682d0/
      data      zpnm( 46)/       3.9910252d0/
      data      zdnm( 46)/       2.2412635d0/
      data      alpm( 46)/       5.1564164d0/
      data     f0sdm( 46)/       8.9626297d0/
      data     g2sdm( 46)/       2.8184354d0/
      data      pocm( 46)/       1.3059919d0/
!
!       Data for Element  47:       Silver
!
      data      ussm( 47)/     -85.4105702d0/
      data      uppm( 47)/     -66.6714256d0/
      data      uddm( 47)/    -137.0858911d0/
      data    betasm( 47)/      -3.7551885d0/
      data    betapm( 47)/      -1.4435948d0/
      data    betadm( 47)/     -15.0893933d0/
      data       zsm( 47)/       2.6156672d0/
      data       zpm( 47)/       1.5209942d0/
      data       zdm( 47)/       3.1178537d0/
      data      zsnm( 47)/       1.4880885d0/
      data      zpnm( 47)/       1.1054787d0/
      data      zdnm( 47)/       2.4424624d0/
      data      alpm( 47)/       2.1347344d0/
      data      pocm( 47)/       1.5299649d0/
!
!       Data for Element  48:      Cadmium
!
      data      ussm( 48)/     -26.1908325d0/
      data      uppm( 48)/     -22.4523736d0/
      data    betasm( 48)/     -11.9610608d0/
      data    betapm( 48)/      -3.9999848d0/
      data       zsm( 48)/       1.4192491d0/
      data       zpm( 48)/       1.0480637d0/
      data      alpm( 48)/       1.1507745d0/
      data      gssm( 48)/      17.2196544d0/
      data      gspm( 48)/      17.9900180d0/
      data      gppm( 48)/      19.2007795d0/
      data      gp2m( 48)/      19.0000212d0/
      data      hspm( 48)/       0.7384511d0/
      data guesm1( 48,1)/      -1.0001225d0/
      data guesm2( 48,1)/       2.3120365d0/
      data guesm3( 48,1)/       1.6572325d0/
!
!       Data for Element  49:       Indium
!
      data      ussm( 49)/     -28.0980892d0/
      data      uppm( 49)/     -19.2780588d0/
      data    betasm( 49)/      -6.3107479d0/
      data    betapm( 49)/      -2.7025837d0/
      data       zsm( 49)/       1.7625740d0/
      data       zpm( 49)/       1.8648962d0/
      data      alpm( 49)/       2.3438756d0/
      data      gssm( 49)/       9.4928794d0/
      data      gspm( 49)/       7.0094241d0/
      data      gppm( 49)/       9.6640986d0/
      data      gp2m( 49)/       7.0100315d0/
      data      hspm( 49)/       0.5995894d0/
      data  guesm1( 49,1)/       0.9277024d0/
      data  guesm2( 49,1)/       0.9999001d0/
      data  guesm3( 49,1)/       1.1829906d0/

!
!       Data for Element  50:          Tin
!
      data      ussm( 50)/     -40.8518020d0/
      data      uppm( 50)/     -28.5602490d0/
      data    betasm( 50)/      -3.2351470d0/
      data    betapm( 50)/      -4.2904160d0/
      data       zsm( 50)/       2.0803800d0/
      data       zpm( 50)/       1.9371060d0/
      data      alpm( 50)/       1.8008140d0/
      data      gssm( 50)/       9.8000000d0/
      data      gspm( 50)/       8.3000000d0/
      data      gppm( 50)/       7.3000000d0/
      data      gp2m( 50)/       6.5000000d0/
      data      hspm( 50)/       1.3000000d0/
!
!       Data for Element  51:     Antimony
!
      data      ussm( 51)/     -42.0643435d0/
      data      uppm( 51)/     -35.0626031d0/
      data    betasm( 51)/      -0.9999715d0/
      data    betapm( 51)/      -4.0920176d0/
      data       zsm( 51)/       3.6458835d0/
      data       zpm( 51)/       1.9733156d0/
      data      alpm( 51)/       1.9763403d0/
      data      gssm( 51)/      10.6739308d0/
      data      gspm( 51)/       7.0477648d0/
      data      gppm( 51)/       6.7446162d0/
      data      gp2m( 51)/       6.3408531d0/
      data      hspm( 51)/       0.5997512d0/
      data guesm1( 51,1)/      -1.0003602d0/
      data guesm2( 51,1)/       0.9992881d0/
      data guesm3( 51,1)/       0.0000000d0/
!
!       Data for Element  52:    Tellurium
!
      data      ussm( 52)/     -84.2274722d0/
      data      uppm( 52)/     -46.5332871d0/
      data    betasm( 52)/      -8.5622652d0/
      data    betapm( 52)/      -2.6942963d0/
      data       zsm( 52)/       2.7461609d0/
      data       zpm( 52)/       1.6160376d0/
      data      alpm( 52)/       2.2924145d0/
      data      gssm( 52)/       5.1367706d0/
      data      gspm( 52)/      11.0720752d0/
      data      gppm( 52)/       5.8447934d0/
      data      gp2m( 52)/       5.0720495d0/
      data      hspm( 52)/       0.5997994d0/
      data guesm1( 52,1)/      -0.6033681d0/
      data guesm2( 52,1)/       1.4127317d0/
      data guesm3( 52,1)/       0.4996755d0/
!
!       Data for Element  53:       Iodine
!
      data      ussm( 53)/    -100.0030538d0/
      data      uppm( 53)/     -74.6114692d0/
      data    betasm( 53)/      -7.4144510d0/
      data    betapm( 53)/      -6.1967810d0/
      data       zsm( 53)/       2.2729610d0/
      data       zpm( 53)/       2.1694980d0/
      data      alpm( 53)/       2.2073200d0/
      data      gssm( 53)/      15.0404486d0/
      data      gspm( 53)/      13.0565580d0/
      data      gppm( 53)/      11.1477837d0/
      data      gp2m( 53)/       9.9140907d0/
      data      hspm( 53)/       2.4563820d0/
      data   polvolm( 53)/       4.5973300d0/
!
!       Data for Element  54:        Xenon
!
      data      ussm( 54)/       5.9205889d0/
      data      uppm( 54)/     -86.9828315d0/
      data    betasm( 54)/      -3.6048445d0/
      data    betapm( 54)/      -4.9673738d0/
      data       zsm( 54)/       4.9900791d0/
      data       zpm( 54)/       2.6929255d0/
      data      alpm( 54)/       1.7948555d0/
      data      gssm( 54)/       2.1874755d0/
      data      gspm( 54)/       4.9038680d0/
      data      gppm( 54)/      12.4939185d0/
      data      gp2m( 54)/      14.9527902d0/
      data      hspm( 54)/       2.4768528d0/
      data guesm1( 54,1)/      -0.2455348d0/
      data guesm2( 54,1)/       2.0580083d0/
      data guesm3( 54,1)/       1.7173301d0/
!
!       Data for Element  55:       Cesium
!
      data      ussm( 55)/      -3.2184078d0/
      data      uppm( 55)/      -1.7699113d0/
      data    betasm( 55)/      -1.6043600d0/
      data    betapm( 55)/      -4.2698040d0/
      data       zsm( 55)/       6.0004170d0/
      data       zpm( 55)/       0.8986916d0/
      data      alpm( 55)/       0.4981646d0/
      data      gssm( 55)/       7.6447851d0/
      data      gspm( 55)/       3.0454989d0/
      data      gppm( 55)/      10.0000745d0/
      data      gp2m( 55)/       6.1761092d0/
      data      hspm( 55)/       0.4647853d0/
!
!       Data for Element  56:       Barium
!
      data      ussm( 56)/     -10.1125345d0/
      data      uppm( 56)/      -8.2347224d0/
      data    betasm( 56)/      -9.9994459d0/
      data    betapm( 56)/      -9.6197255d0/
      data       zsm( 56)/       1.9765973d0/
      data       zpm( 56)/       1.3157348d0/
      data      alpm( 56)/       0.8594840d0/
      data      gssm( 56)/       4.8486178d0/
      data      gspm( 56)/       4.5659982d0/
      data      gppm( 56)/       5.0937708d0/
      data      gp2m( 56)/       5.2125824d0/
      data      hspm( 56)/       0.5237082d0/
!
!       Data for Element  78:     Platinum
!
      data      ussm( 78)/     -80.0844280d0/
      data      uppm( 78)/     -52.9422696d0/
      data      uddm( 78)/    -111.3718762d0/
      data    betasm( 78)/      -4.7115964d0/
      data    betapm( 78)/      -3.0319706d0/
      data    betadm( 78)/     -10.3250349d0/
      data       zsm( 78)/       1.8655763d0/
      data       zpm( 78)/       1.9475781d0/
      data       zdm( 78)/       2.8552253d0/
      data      zsnm( 78)/       1.8965538d0/
      data      zpnm( 78)/       1.6684682d0/
      data      zdnm( 78)/       2.7361239d0/
      data      alpm( 78)/       0.4999907d0/
      data      pocm( 78)/       1.3743794d0/
!
!       Data for Element  80:      Mercury
!
      data      ussm( 80)/     -19.8095740d0/
      data      uppm( 80)/     -13.1025300d0/
      data    betasm( 80)/      -0.4045250d0/
      data    betapm( 80)/      -6.2066830d0/
      data       zsm( 80)/       2.2181840d0/
      data       zpm( 80)/       2.0650380d0/
      data      alpm( 80)/       1.3356410d0/
      data      gssm( 80)/      10.8000000d0/
      data      gspm( 80)/       9.3000000d0/
      data      gppm( 80)/      14.3000000d0/
      data      gp2m( 80)/      13.5000000d0/
      data      hspm( 80)/       1.3000000d0/
!
!       Data for Element  81:     Thallium
!
      data      ussm( 81)/     -29.7009655d0/
      data      uppm( 81)/     -29.5937539d0/
      data    betasm( 81)/      -4.9667442d0/
      data    betapm( 81)/      -7.7616060d0/
      data       zsm( 81)/       4.0000447d0/
      data       zpm( 81)/       1.8076332d0/
      data      alpm( 81)/       1.3116968d0/
      data      gssm( 81)/       8.8675337d0/
      data      gspm( 81)/      12.1148290d0/
      data      gppm( 81)/      10.6532769d0/
      data      gp2m( 81)/      13.5333191d0/
      data      hspm( 81)/       0.5997565d0/
      data guesm1( 81,1)/      -0.7940727d0/
      data guesm2( 81,1)/       0.9999962d0/
      data guesm3( 81,1)/       0.4999732d0/
!
!       Data for Element  82:         Lead
!
      data      ussm( 82)/     -47.3196920d0/
      data      uppm( 82)/     -28.8475600d0/
      data    betasm( 82)/      -8.0423870d0/
      data    betapm( 82)/      -3.0000000d0/
      data       zsm( 82)/       2.4982860d0/
      data       zpm( 82)/       2.0820710d0/
      data      alpm( 82)/       1.7283330d0/
      data      gssm( 82)/       9.8000000d0/
      data      gspm( 82)/       8.3000000d0/
      data      gppm( 82)/       7.3000000d0/
      data      gp2m( 82)/       6.5000000d0/
      data      hspm( 82)/       1.3000000d0/
!
!       Data for Element  83:      Bismuth
!
      data      ussm( 83)/     -53.5827147d0/
      data      uppm( 83)/     -39.4572213d0/
      data    betasm( 83)/      -9.0000249d0/
      data    betapm( 83)/      -1.9830269d0/
      data       zsm( 83)/       2.6772255d0/
      data       zpm( 83)/       0.6936864d0/
      data      alpm( 83)/       5.7660628d0/
      data      gssm( 83)/       8.3702778d0/
      data      gspm( 83)/       7.7974668d0/
      data      gppm( 83)/       9.8303165d0/
      data      gp2m( 83)/       8.9291355d0/
      data      hspm( 83)/       0.5999908d0/
      data guesm1( 83,1)/      -0.1281535d0/
      data guesm2( 83,1)/       3.0003211d0/
      data guesm3( 83,1)/       1.7993215d0/
!
!       Data for Element  85:     Astatine
!
      data      alpm( 85)/       3.0000000d0/
      data      gssm( 85)/      10.0000000d0/
!
!       Data for Element  87:     Francium
!
      data      alpm( 87)/       3.0000000d0/
      data      gssm( 87)/      10.0000000d0/
!
!       Data for Element  90:      Thorium
!
      data      ussm( 90)/     -40.5682920d0/
      data      uppm( 90)/     -28.0891870d0/
      data    betasm( 90)/      -4.2562180d0/
      data    betapm( 90)/      -4.2562180d0/
      data       zsm( 90)/       1.4353060d0/
      data       zpm( 90)/       1.4353060d0/
      data      alpm( 90)/       2.1961078d0/
      data      gssm( 90)/       9.8200000d0/
      data      gspm( 90)/       8.3600000d0/
      data      gppm( 90)/       7.3100000d0/
      data      gp2m( 90)/       6.5400000d0/
      data      hspm( 90)/       1.3200000d0/
!
!       Data for Element 100:   3+ Sparkle
!
      data      alpm(100)/       1.5000000d0/
!
!       Data for Element 101:   3- Sparkle
!
      data      alpm(101)/       1.5000000d0/
!
!       Data for Element 102:  Capped bond
!
      data      ussm(102)/     -11.9062760d0/
      data    betasm(102)/-9999999.0000000d0/
      data       zsm(102)/       4.0000000d0/
      data       zpm(102)/       0.3000000d0/
      data       zdm(102)/       0.3000000d0/
      data      alpm(102)/       2.5441341d0/
      data      gssm(102)/      12.8480000d0/
      data      hspm(102)/       0.1000000d0/
!
!       Data for Element 103:   ++ Sparkle
!
      data      alpm(103)/       1.5000000d0/
!
!       Data for Element 104:    + Sparkle
!
      data      alpm(104)/       1.5000000d0/
!
!       Data for Element 105:   -- Sparkle
!
      data      alpm(105)/       1.5000000d0/
!
!       Data for Element 106:    - Sparkle
!
      data      alpm(106)/       1.5000000d0/
      end module Parameters_for_MNDO_C
