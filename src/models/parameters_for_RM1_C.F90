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

  module Parameters_for_RM1_C
    double precision, dimension(107) :: ussRM1, uppRM1, uddRM1, zsRM1, zpRM1, zdRM1, betasRM1, &
    betapRM1, betadRM1, gssRM1, gspRM1, gppRM1, gp2RM1, hspRM1, alpRM1, zsnRM1, zpnRM1, zdnRM1, &
    poc_RM1, f0sdRM1, g2sdRM1
    character (len=80), dimension(107) :: refRM1
    double precision, dimension(107,4) :: guess1RM1, guess2RM1, guess3RM1
!
!                    Data for Element  1         Hydrogen
!
      data     ussRM1(  1)/     -11.9606770D0/
      data   betasRM1(  1)/      -5.7654447D0/
      data      zsRM1(  1)/       1.0826737D0/
      data     alpRM1(  1)/       3.0683595D0/
      data     gssRM1(  1)/      13.9832130D0/
      data guess1RM1(  1,1)/       0.1028888D0/
      data guess2RM1(  1,1)/       5.9017227D0/
      data guess3RM1(  1,1)/       1.1750118D0/
      data guess1RM1(  1,2)/       0.0645745D0/
      data guess2RM1(  1,2)/       6.4178567D0/
      data guess3RM1(  1,2)/       1.9384448D0/
      data guess1RM1(  1,3)/      -0.0356739D0/
      data guess2RM1(  1,3)/       2.8047313D0/
      data guess3RM1(  1,3)/       1.6365524D0/
!
!                    Data for Element  6           Carbon
!
      data     ussRM1(  6)/     -51.7255603D0/
      data     uppRM1(  6)/     -39.4072894D0/
      data   betasRM1(  6)/     -15.4593243D0/
      data   betapRM1(  6)/      -8.2360864D0/
      data      zsRM1(  6)/       1.8501880D0/
      data      zpRM1(  6)/       1.7683009D0/
      data     alpRM1(  6)/       2.7928208D0/
      data     gssRM1(  6)/      13.0531244D0/
      data     gspRM1(  6)/      11.3347939D0/
      data     gppRM1(  6)/      10.9511374D0/
      data     gp2RM1(  6)/       9.7239510D0/
      data     hspRM1(  6)/       1.5521513D0/
      data guess1RM1(  6,1)/       0.0746227D0/
      data guess2RM1(  6,1)/       5.7392160D0/
      data guess3RM1(  6,1)/       1.0439698D0/
      data guess1RM1(  6,2)/       0.0117705D0/
      data guess2RM1(  6,2)/       6.9240173D0/
      data guess3RM1(  6,2)/       1.6615957D0/
      data guess1RM1(  6,3)/       0.0372066D0/
      data guess2RM1(  6,3)/       6.2615894D0/
      data guess3RM1(  6,3)/       1.6315872D0/
      data guess1RM1(  6,4)/      -0.0027066D0/
      data guess2RM1(  6,4)/       9.0000373D0/
      data guess3RM1(  6,4)/       2.7955790D0/
!
!                    Data for Element  7         Nitrogen
!
      data     ussRM1(  7)/     -70.8512372D0/
      data     uppRM1(  7)/     -57.9773092D0/
      data   betasRM1(  7)/     -20.8712455D0/
      data   betapRM1(  7)/     -16.6717185D0/
      data      zsRM1(  7)/       2.3744716D0/
      data      zpRM1(  7)/       1.9781257D0/
      data     alpRM1(  7)/       2.9642254D0/
      data     gssRM1(  7)/      13.0873623D0/
      data     gspRM1(  7)/      13.2122683D0/
      data     gppRM1(  7)/      13.6992432D0/
      data     gp2RM1(  7)/      11.9410395D0/
      data     hspRM1(  7)/       5.0000085D0/
      data guess1RM1(  7,1)/       0.0607338D0/
      data guess2RM1(  7,1)/       4.5889295D0/
      data guess3RM1(  7,1)/       1.3787388D0/
      data guess1RM1(  7,2)/       0.0243856D0/
      data guess2RM1(  7,2)/       4.6273052D0/
      data guess3RM1(  7,2)/       2.0837070D0/
      data guess1RM1(  7,3)/      -0.0228343D0/
      data guess2RM1(  7,3)/       2.0527466D0/
      data guess3RM1(  7,3)/       1.8676382D0/
!
!                    Data for Element  8           Oxygen
!
      data     ussRM1(  8)/     -96.9494807D0/
      data     uppRM1(  8)/     -77.8909298D0/
      data   betasRM1(  8)/     -29.8510121D0/
      data   betapRM1(  8)/     -29.1510131D0/
      data      zsRM1(  8)/       3.1793691D0/
      data      zpRM1(  8)/       2.5536191D0/
      data     alpRM1(  8)/       4.1719672D0/
      data     gssRM1(  8)/      14.0024279D0/
      data     gspRM1(  8)/      14.9562504D0/
      data     gppRM1(  8)/      14.1451514D0/
      data     gp2RM1(  8)/      12.7032550D0/
      data     hspRM1(  8)/       3.9321716D0/
      data guess1RM1(  8,1)/       0.2309355D0/
      data guess2RM1(  8,1)/       5.2182874D0/
      data guess3RM1(  8,1)/       0.9036356D0/
      data guess1RM1(  8,2)/       0.0585987D0/
      data guess2RM1(  8,2)/       7.4293293D0/
      data guess3RM1(  8,2)/       1.5175461D0/
!
!                    Data for Element  9         Fluorine
!
      data     ussRM1(  9)/    -134.1836959D0/
      data     uppRM1(  9)/    -107.8466092D0/
      data   betasRM1(  9)/     -70.0000051D0/
      data   betapRM1(  9)/     -32.6798271D0/
      data      zsRM1(  9)/       4.4033791D0/
      data      zpRM1(  9)/       2.6484156D0/
      data     alpRM1(  9)/       6.0000006D0/
      data     gssRM1(  9)/      16.7209132D0/
      data     gspRM1(  9)/      16.7614263D0/
      data     gppRM1(  9)/      15.2258103D0/
      data     gp2RM1(  9)/      14.8657868D0/
      data     hspRM1(  9)/       1.9976617D0/
      data guess1RM1(  9,1)/       0.4030203D0/
      data guess2RM1(  9,1)/       7.2044196D0/
      data guess3RM1(  9,1)/       0.8165301D0/
      data guess1RM1(  9,2)/       0.0708583D0/
      data guess2RM1(  9,2)/       9.0000156D0/
      data guess3RM1(  9,2)/       1.4380238D0/
!
!                    Data for Element 15       Phosphorus
!
      data     ussRM1( 15)/     -41.8153318D0/
      data     uppRM1( 15)/     -34.3834253D0/
      data   betasRM1( 15)/      -6.1351497D0/
      data   betapRM1( 15)/      -5.9444213D0/
      data      zsRM1( 15)/       2.1224012D0/
      data      zpRM1( 15)/       1.7432795D0/
      data     alpRM1( 15)/       1.9099329D0/
      data     gssRM1( 15)/      11.0805926D0/
      data     gspRM1( 15)/       5.6833920D0/
      data     gppRM1( 15)/       7.6041756D0/
      data     gp2RM1( 15)/       7.4026518D0/
      data     hspRM1( 15)/       1.1618179D0/
      data guess1RM1( 15,1)/      -0.4106347D0/
      data guess2RM1( 15,1)/       6.0875283D0/
      data guess3RM1( 15,1)/       1.3165026D0/
      data guess1RM1( 15,2)/      -0.1629929D0/
      data guess2RM1( 15,2)/       7.0947260D0/
      data guess3RM1( 15,2)/       1.9072132D0/
      data guess1RM1( 15,3)/      -0.0488713D0/
      data guess2RM1( 15,3)/       8.9997931D0/
      data guess3RM1( 15,3)/       2.6585778D0/
!
!                    Data for Element 16           Sulfur
!
      data     ussRM1( 16)/     -55.1677512D0/
      data     uppRM1( 16)/     -46.5293042D0/
      data   betasRM1( 16)/      -1.9591072D0/
      data   betapRM1( 16)/      -8.7743065D0/
      data      zsRM1( 16)/       2.1334431D0/
      data      zpRM1( 16)/       1.8746065D0/
      data     alpRM1( 16)/       2.4401564D0/
      data     gssRM1( 16)/      12.4882841D0/
      data     gspRM1( 16)/       8.5691057D0/
      data     gppRM1( 16)/       8.5230117D0/
      data     gp2RM1( 16)/       7.6686330D0/
      data     hspRM1( 16)/       3.8897893D0/
      data guess1RM1( 16,1)/      -0.7460106D0/
      data guess2RM1( 16,1)/       4.8103800D0/
      data guess3RM1( 16,1)/       0.5938013D0/
      data guess1RM1( 16,2)/      -0.0651929D0/
      data guess2RM1( 16,2)/       7.2076086D0/
      data guess3RM1( 16,2)/       1.2949201D0/
      data guess1RM1( 16,3)/      -0.0065598D0/
      data guess2RM1( 16,3)/       9.0000018D0/
      data guess3RM1( 16,3)/       1.8006015D0/
!
!                    Data for Element 17         Chlorine
!
      data     ussRM1( 17)/    -118.4730692D0/
      data     uppRM1( 17)/     -76.3533034D0/
      data   betasRM1( 17)/     -19.9243043D0/
      data   betapRM1( 17)/     -11.5293520D0/
      data      zsRM1( 17)/       3.8649107D0/
      data      zpRM1( 17)/       1.8959314D0/
      data     alpRM1( 17)/       3.6935883D0/
      data     gssRM1( 17)/      15.3602310D0/
      data     gspRM1( 17)/      13.3067117D0/
      data     gppRM1( 17)/      12.5650264D0/
      data     gp2RM1( 17)/       9.6639708D0/
      data     hspRM1( 17)/       1.7648990D0/
      data guess1RM1( 17,1)/       0.1294711D0/
      data guess2RM1( 17,1)/       2.9772442D0/
      data guess3RM1( 17,1)/       1.4674978D0/
      data guess1RM1( 17,2)/       0.0028890D0/
      data guess2RM1( 17,2)/       7.0982759D0/
      data guess3RM1( 17,2)/       2.5000272D0/

!
!                    Data for Element 35          Bromine
!
      data     ussRM1( 35)/    -113.4839818D0/
      data     uppRM1( 35)/     -76.1872002D0/
      data   betasRM1( 35)/      -1.3413984D0/
      data   betapRM1( 35)/      -8.2022599D0/
      data      zsRM1( 35)/       5.7315721D0/
      data      zpRM1( 35)/       2.0314758D0/
      data     alpRM1( 35)/       2.8671053D0/
      data     gssRM1( 35)/      17.1156307D0/
      data     gspRM1( 35)/      15.6241925D0/
      data     gppRM1( 35)/      10.7354629D0/
      data     gp2RM1( 35)/       8.8605620D0/
      data     hspRM1( 35)/       2.2351276D0/
      data guess1RM1( 35,1)/       0.9868994D0/
      data guess2RM1( 35,1)/       4.2848419D0/
      data guess3RM1( 35,1)/       2.0001970D0/
      data guess1RM1( 35,2)/      -0.9273125D0/
      data guess2RM1( 35,2)/       4.5400591D0/
      data guess3RM1( 35,2)/       2.0161770D0/
!
!                    Data for Element 53           Iodine
!
      data     ussRM1( 53)/     -74.8999784D0/
      data     uppRM1( 53)/     -51.4102380D0/
      data   betasRM1( 53)/      -4.1931615D0/
      data   betapRM1( 53)/      -4.4003841D0/
      data      zsRM1( 53)/       2.5300375D0/
      data      zpRM1( 53)/       2.3173868D0/
      data     alpRM1( 53)/       2.1415709D0/
      data     gssRM1( 53)/      19.9997413D0/
      data     gspRM1( 53)/       7.6895767D0/
      data     gppRM1( 53)/       7.3048834D0/
      data     gp2RM1( 53)/       6.8542461D0/
      data     hspRM1( 53)/       1.4160294D0/
      data guess1RM1( 53,1)/      -0.0814772D0/
      data guess2RM1( 53,1)/       1.5606507D0/
      data guess3RM1( 53,1)/       2.0000206D0/
      data guess1RM1( 53,2)/       0.0591499D0/
      data guess2RM1( 53,2)/       5.7611127D0/
      data guess3RM1( 53,2)/       2.2048880D0/
!
!                    Data for Element  57        Lanthanum
!
      data     ussRM1( 57)/       -14.680434D0/
      data     uppRM1( 57)/        -6.734739D0/
      data     uddRM1( 57)/       -20.489967D0/
      data   betasRM1( 57)/        -7.669555D0/
      data   betapRM1( 57)/         0.477696D0/
      data   betadRM1( 57)/        -3.711477D0/
      data      zsRM1( 57)/         1.272677D0/
      data      zpRM1( 57)/         1.423276D0/
      data      zdRM1( 57)/         1.410369D0/
      data     zsnRM1( 57)/         0.784530D0/
      data     zpnRM1( 57)/         1.506612D0/
      data     zdnRM1( 57)/         1.172062D0/
      data     alpRM1( 57)/         1.284047D0/
      data     gssRM1( 57)/         2.984540D0/
      data     gspRM1( 57)/         3.510670D0/
      data     gppRM1( 57)/         6.287652D0/
      data     gp2RM1( 57)/         5.453445D0/
      data     hspRM1( 57)/         0.284330D0/
      data    poc_RM1( 57)/         1.875176D0/
      data    f0sdRM1( 57)/         7.720813D0/
      data    g2sdRM1( 57)/         3.916745D0/
      data guess1RM1( 57,1)/         0.628559D0/
      data guess2RM1( 57,1)/         7.860849D0/
      data guess3RM1( 57,1)/         1.304475D0/
      data guess1RM1( 57,2)/         0.081642D0/
      data guess2RM1( 57,2)/        10.346858D0/
      data guess3RM1( 57,2)/         3.247040D0/
!
!                    Data for Element  58           Cerium
!
      data     ussRM1( 58)/       -14.719389D0/
      data     uppRM1( 58)/        -7.689429D0/
      data     uddRM1( 58)/       -20.451577D0/
      data   betasRM1( 58)/        -7.668787D0/
      data   betapRM1( 58)/         0.444442D0/
      data   betadRM1( 58)/        -3.744938D0/
      data      zsRM1( 58)/         1.281028D0/
      data      zpRM1( 58)/         1.425366D0/
      data      zdRM1( 58)/         1.412866D0/
      data     zsnRM1( 58)/         0.819413D0/
      data     zpnRM1( 58)/         1.430349D0/
      data     zdnRM1( 58)/         1.193220D0/
      data     alpRM1( 58)/         1.286233D0/
      data     gssRM1( 58)/         3.117245D0/
      data     gspRM1( 58)/         3.637115D0/
      data     gppRM1( 58)/         5.969378D0/
      data     gp2RM1( 58)/         5.177398D0/
      data     hspRM1( 58)/         0.401838D0/
      data    poc_RM1( 58)/         1.875084D0/
      data    f0sdRM1( 58)/         7.715129D0/
      data    g2sdRM1( 58)/         3.918293D0/
      data guess1RM1( 58,1)/         0.661961D0/
      data guess2RM1( 58,1)/         7.890258D0/
      data guess3RM1( 58,1)/         1.262827D0/
      data guess1RM1( 58,2)/         0.075992D0/
      data guess2RM1( 58,2)/        10.316109D0/
      data guess3RM1( 58,2)/         3.247817D0/
!
!                    Data for Element  59     Praseodymium
!
      data     ussRM1( 59)/       -14.524081D0/
      data     uppRM1( 59)/        -7.056827D0/
      data     uddRM1( 59)/       -20.689328D0/
      data   betasRM1( 59)/        -7.947993D0/
      data   betapRM1( 59)/         0.853816D0/
      data   betadRM1( 59)/        -3.830293D0/
      data      zsRM1( 59)/         1.538039D0/
      data      zpRM1( 59)/         1.581647D0/
      data      zdRM1( 59)/         1.374904D0/
      data     zsnRM1( 59)/         0.781944D0/
      data     zpnRM1( 59)/         1.292982D0/
      data     zdnRM1( 59)/         0.989605D0/
      data     alpRM1( 59)/         1.280603D0/
      data     gssRM1( 59)/         2.974704D0/
      data     gspRM1( 59)/         3.449344D0/
      data     gppRM1( 59)/         5.396095D0/
      data     gp2RM1( 59)/         4.680174D0/
      data     hspRM1( 59)/         0.444736D0/
      data    poc_RM1( 59)/         1.847416D0/
      data    f0sdRM1( 59)/         7.618301D0/
      data    g2sdRM1( 59)/         3.965863D0/
      data guess1RM1( 59,1)/         0.453379D0/
      data guess2RM1( 59,1)/         7.823199D0/
      data guess3RM1( 59,1)/         1.565167D0/
      data guess1RM1( 59,2)/         0.010394D0/
      data guess2RM1( 59,2)/        10.288365D0/
      data guess3RM1( 59,2)/         3.268703D0/
!
!                    Data for Element  60        Neodymium
!
      data     ussRM1( 60)/       -14.633586D0/
      data     uppRM1( 60)/        -7.047863D0/
      data     uddRM1( 60)/       -19.641391D0/
      data   betasRM1( 60)/        -7.929373D0/
      data   betapRM1( 60)/         0.964369D0/
      data   betadRM1( 60)/        -3.812351D0/
      data      zsRM1( 60)/         1.458290D0/
      data      zpRM1( 60)/         1.570516D0/
      data      zdRM1( 60)/         1.513561D0/
      data     zsnRM1( 60)/         0.809460D0/
      data     zpnRM1( 60)/         1.441356D0/
      data     zdnRM1( 60)/         0.997255D0/
      data     alpRM1( 60)/         1.289770D0/
      data     gssRM1( 60)/         3.079380D0/
      data     gspRM1( 60)/         3.600006D0/
      data     gppRM1( 60)/         6.015314D0/
      data     gp2RM1( 60)/         5.217239D0/
      data     hspRM1( 60)/         0.374325D0/
      data    poc_RM1( 60)/         1.713827D0/
      data    f0sdRM1( 60)/         7.673391D0/
      data    g2sdRM1( 60)/         4.134872D0/
      data guess1RM1( 60,1)/         0.396892D0/
      data guess2RM1( 60,1)/         7.706685D0/
      data guess3RM1( 60,1)/         1.589180D0/
      data guess1RM1( 60,2)/         0.027000D0/
      data guess2RM1( 60,2)/        10.322645D0/
      data guess3RM1( 60,2)/         3.236907D0/
!
!                    Data for Element  61       Promethium
!
      data     ussRM1( 61)/       -15.493641D0/
      data     uppRM1( 61)/        -7.101212D0/
      data     uddRM1( 61)/       -19.423422D0/
      data   betasRM1( 61)/        -7.982257D0/
      data   betapRM1( 61)/         0.977908D0/
      data   betadRM1( 61)/        -3.925565D0/
      data      zsRM1( 61)/         1.065536D0/
      data      zpRM1( 61)/         1.846925D0/
      data      zdRM1( 61)/         1.424049D0/
      data     zsnRM1( 61)/         0.746406D0/
      data     zpnRM1( 61)/         1.467790D0/
      data     zdnRM1( 61)/         0.769147D0/
      data     alpRM1( 61)/         1.393272D0/
      data     gssRM1( 61)/         2.839509D0/
      data     gspRM1( 61)/         3.345263D0/
      data     gppRM1( 61)/         6.125633D0/
      data     gp2RM1( 61)/         5.312922D0/
      data     hspRM1( 61)/         0.248937D0/
      data    poc_RM1( 61)/         1.756461D0/
      data    f0sdRM1( 61)/         7.985070D0/
      data    g2sdRM1( 61)/         4.129014D0/
      data guess1RM1( 61,1)/         0.347564D0/
      data guess2RM1( 61,1)/         7.683601D0/
      data guess3RM1( 61,1)/         1.649495D0/
      data guess1RM1( 61,2)/         0.065302D0/
      data guess2RM1( 61,2)/         9.917302D0/
      data guess3RM1( 61,2)/         3.117827D0/
!
!                    Data for Element  62         Samarium
!
      data     ussRM1( 62)/       -16.509324D0/
      data     uppRM1( 62)/        -7.018030D0/
      data     uddRM1( 62)/       -18.873252D0/
      data   betasRM1( 62)/        -8.077902D0/
      data   betapRM1( 62)/         0.939330D0/
      data   betadRM1( 62)/        -4.019389D0/
      data      zsRM1( 62)/         1.293914D0/
      data      zpRM1( 62)/         1.738656D0/
      data      zdRM1( 62)/         1.521378D0/
      data     zsnRM1( 62)/         0.663480D0/
      data     zpnRM1( 62)/         1.492303D0/
      data     zdnRM1( 62)/         0.853154D0/
      data     alpRM1( 62)/         1.307705D0/
      data     gssRM1( 62)/         2.524039D0/
      data     gspRM1( 62)/         2.992065D0/
      data     gppRM1( 62)/         6.227938D0/
      data     gp2RM1( 62)/         5.401653D0/
      data     hspRM1( 62)/         0.130304D0/
      data    poc_RM1( 62)/         1.739423D0/
      data    f0sdRM1( 62)/         7.999828D0/
      data    g2sdRM1( 62)/         4.163729D0/
      data guess1RM1( 62,1)/         0.253884D0/
      data guess2RM1( 62,1)/         7.713826D0/
      data guess3RM1( 62,1)/         1.731093D0/
      data guess1RM1( 62,2)/         0.037754D0/
      data guess2RM1( 62,2)/        10.013365D0/
      data guess3RM1( 62,2)/         3.210086D0/
!
!                    Data for Element  63         Europium
!
      data     ussRM1( 63)/       -18.544097D0/
      data     uppRM1( 63)/        -6.936941D0/
      data     uddRM1( 63)/       -19.870881D0/
      data   betasRM1( 63)/        -8.168076D0/
      data   betapRM1( 63)/         0.880709D0/
      data   betadRM1( 63)/        -3.987470D0/
      data      zsRM1( 63)/         1.350342D0/
      data      zpRM1( 63)/         1.733714D0/
      data      zdRM1( 63)/         1.494122D0/
      data     zsnRM1( 63)/         0.679880D0/
      data     zpnRM1( 63)/         1.578577D0/
      data     zdnRM1( 63)/         0.798561D0/
      data     alpRM1( 63)/         1.334371D0/
      data     gssRM1( 63)/         2.586427D0/
      data     gspRM1( 63)/         3.068954D0/
      data     gppRM1( 63)/         6.587990D0/
      data     gp2RM1( 63)/         5.713936D0/
      data     hspRM1( 63)/         0.116134D0/
      data    poc_RM1( 63)/         1.758060D0/
      data    f0sdRM1( 63)/         8.017234D0/
      data    g2sdRM1( 63)/         3.095150D0/
      data guess1RM1( 63,1)/         0.233304D0/
      data guess2RM1( 63,1)/         7.759708D0/
      data guess3RM1( 63,1)/         1.747035D0/
      data guess1RM1( 63,2)/         0.019498D0/
      data guess2RM1( 63,2)/         9.971389D0/
      data guess3RM1( 63,2)/         3.229312D0/
!
!                    Data for Element  64       Gadolinium
!
      data     ussRM1( 64)/       -19.908768D0/
      data     uppRM1( 64)/        -7.758304D0/
      data     uddRM1( 64)/       -18.965077D0/
      data   betasRM1( 64)/        -7.588699D0/
      data   betapRM1( 64)/        -2.041090D0/
      data   betadRM1( 64)/        -4.281267D0/
      data      zsRM1( 64)/         1.272776D0/
      data      zpRM1( 64)/         1.908122D0/
      data      zdRM1( 64)/         1.515905D0/
      data     zsnRM1( 64)/         0.986492D0/
      data     zpnRM1( 64)/         2.043026D0/
      data     zdnRM1( 64)/         0.976186D0/
      data     alpRM1( 64)/         1.296754D0/
      data     gssRM1( 64)/         3.752853D0/
      data     gspRM1( 64)/         4.433961D0/
      data     gppRM1( 64)/         8.526309D0/
      data     gp2RM1( 64)/         7.395090D0/
      data     hspRM1( 64)/         0.271406D0/
      data    poc_RM1( 64)/         1.521880D0/
      data    f0sdRM1( 64)/         8.204767D0/
      data    g2sdRM1( 64)/         1.639570D0/
      data guess1RM1( 64,1)/         1.030249D0/
      data guess2RM1( 64,1)/         7.375972D0/
      data guess3RM1( 64,1)/         1.549473D0/
      data guess1RM1( 64,2)/         0.035904D0/
      data guess2RM1( 64,2)/         7.596826D0/
      data guess3RM1( 64,2)/         3.036203D0/
!
!                    Data for Element  65          Terbium
!
      data     ussRM1( 65)/       -20.926991D0/
      data     uppRM1( 65)/        -7.752425D0/
      data     uddRM1( 65)/       -19.971117D0/
      data   betasRM1( 65)/        -7.575097D0/
      data   betapRM1( 65)/        -2.056425D0/
      data   betadRM1( 65)/        -4.373081D0/
      data      zsRM1( 65)/         1.210052D0/
      data      zpRM1( 65)/         1.921514D0/
      data      zdRM1( 65)/         1.528123D0/
      data     zsnRM1( 65)/         1.206260D0/
      data     zpnRM1( 65)/         2.068824D0/
      data     zdnRM1( 65)/         0.858750D0/
      data     alpRM1( 65)/         1.298316D0/
      data     gssRM1( 65)/         4.588904D0/
      data     gspRM1( 65)/         5.344161D0/
      data     gppRM1( 65)/         8.633975D0/
      data     gp2RM1( 65)/         7.488472D0/
      data     hspRM1( 65)/         0.621955D0/
      data    poc_RM1( 65)/         1.525065D0/
      data    f0sdRM1( 65)/         8.212703D0/
      data    g2sdRM1( 65)/         1.336917D0/
      data guess1RM1( 65,1)/         1.037032D0/
      data guess2RM1( 65,1)/         7.599123D0/
      data guess3RM1( 65,1)/         1.570196D0/
      data guess1RM1( 65,2)/         0.046496D0/
      data guess2RM1( 65,2)/         7.618297D0/
      data guess3RM1( 65,2)/         3.128177D0/
!
!                    Data for Element  66       Dysprosium
!
      data     ussRM1( 66)/       -20.926240D0/
      data     uppRM1( 66)/        -7.667306D0/
      data     uddRM1( 66)/       -17.940815D0/
      data   betasRM1( 66)/        -7.606705D0/
      data   betapRM1( 66)/         1.961734D0/
      data   betadRM1( 66)/        -4.368527D0/
      data      zsRM1( 66)/         1.295275D0/
      data      zpRM1( 66)/         1.912107D0/
      data      zdRM1( 66)/         1.413397D0/
      data     zsnRM1( 66)/         1.372366D0/
      data     zpnRM1( 66)/         1.074073D0/
      data     zdnRM1( 66)/         0.819144D0/
      data     alpRM1( 66)/         1.348259D0/
      data     gssRM1( 66)/         5.220812D0/
      data     gspRM1( 66)/         4.498122D0/
      data     gppRM1( 66)/         4.482505D0/
      data     gp2RM1( 66)/         3.887794D0/
      data     hspRM1( 66)/         0.973576D0/
      data    poc_RM1( 66)/         1.625055D0/
      data    f0sdRM1( 66)/         8.305431D0/
      data    g2sdRM1( 66)/         1.310365D0/
      data guess1RM1( 66,1)/         1.130715D0/
      data guess2RM1( 66,1)/         7.711956D0/
      data guess3RM1( 66,1)/         1.536658D0/
      data guess1RM1( 66,2)/         0.068456D0/
      data guess2RM1( 66,2)/         7.506540D0/
      data guess3RM1( 66,2)/         3.234171D0/
!
!                    Data for Element  67          Holmium
!
      data     ussRM1( 67)/       -22.057459D0/
      data     uppRM1( 67)/        -7.595638D0/
      data     uddRM1( 67)/       -18.000406D0/
      data   betasRM1( 67)/        -5.645226D0/
      data   betapRM1( 67)/         0.006537D0/
      data   betadRM1( 67)/        -4.312899D0/
      data      zsRM1( 67)/         1.330550D0/
      data      zpRM1( 67)/         1.779559D0/
      data      zdRM1( 67)/         1.536524D0/
      data     zsnRM1( 67)/         1.498038D0/
      data     zpnRM1( 67)/         1.967497D0/
      data     zdnRM1( 67)/         0.663021D0/
      data     alpRM1( 67)/         1.330075D0/
      data     gssRM1( 67)/         5.698899D0/
      data     gspRM1( 67)/         6.321991D0/
      data     gppRM1( 67)/         8.211100D0/
      data     gp2RM1( 67)/         7.121701D0/
      data     hspRM1( 67)/         1.317056D0/
      data    poc_RM1( 67)/         1.719560D0/
      data    f0sdRM1( 67)/         8.240569D0/
      data    g2sdRM1( 67)/         1.245432D0/
      data guess1RM1( 67,1)/         1.090708D0/
      data guess2RM1( 67,1)/         7.571516D0/
      data guess3RM1( 67,1)/         1.490954D0/
      data guess1RM1( 67,2)/         0.001419D0/
      data guess2RM1( 67,2)/         7.799696D0/
      data guess3RM1( 67,2)/         3.254251D0/
!
!                    Data for Element  68           Erbium
!
      data     ussRM1( 68)/       -21.978399D0/
      data     uppRM1( 68)/        -7.607850D0/
      data     uddRM1( 68)/       -17.976841D0/
      data   betasRM1( 68)/        -5.634710D0/
      data   betapRM1( 68)/        -0.018972D0/
      data   betadRM1( 68)/        -4.250679D0/
      data      zsRM1( 68)/         1.347757D0/
      data      zpRM1( 68)/         1.806481D0/
      data      zdRM1( 68)/         1.466189D0/
      data     zsnRM1( 68)/         1.446757D0/
      data     zpnRM1( 68)/         1.973883D0/
      data     zdnRM1( 68)/         0.650461D0/
      data     alpRM1( 68)/         1.320103D0/
      data     gssRM1( 68)/         5.503813D0/
      data     gspRM1( 68)/         6.164998D0/
      data     gppRM1( 68)/         8.237750D0/
      data     gp2RM1( 68)/         7.144816D0/
      data     hspRM1( 68)/         1.210092D0/
      data    poc_RM1( 68)/         2.717136D0/
      data    f0sdRM1( 68)/         8.257327D0/
      data    g2sdRM1( 68)/         1.248745D0/
      data guess1RM1( 68,1)/         1.174177D0/
      data guess2RM1( 68,1)/         7.583252D0/
      data guess3RM1( 68,1)/         1.503549D0/
      data guess1RM1( 68,2)/         0.008646D0/
      data guess2RM1( 68,2)/         7.813788D0/
      data guess3RM1( 68,2)/         3.233597D0/
!
!                    Data for Element  69          Thulium
!
      data     ussRM1( 69)/       -21.890870D0/
      data     uppRM1( 69)/        -7.252807D0/
      data     uddRM1( 69)/       -18.183885D0/
      data   betasRM1( 69)/        -5.480593D0/
      data   betapRM1( 69)/         0.078860D0/
      data   betadRM1( 69)/        -4.321791D0/
      data      zsRM1( 69)/         1.369147D0/
      data      zpRM1( 69)/         1.674365D0/
      data      zdRM1( 69)/         1.714394D0/
      data     zsnRM1( 69)/         1.238387D0/
      data     zpnRM1( 69)/         1.821235D0/
      data     zdnRM1( 69)/         0.956205D0/
      data     alpRM1( 69)/         1.266431D0/
      data     gssRM1( 69)/         4.711122D0/
      data     gspRM1( 69)/         5.362973D0/
      data     gppRM1( 69)/         7.600691D0/
      data     gp2RM1( 69)/         6.592277D0/
      data     hspRM1( 69)/         0.914909D0/
      data    poc_RM1( 69)/         2.763482D0/
      data    f0sdRM1( 69)/         8.326130D0/
      data    g2sdRM1( 69)/         1.474432D0/
      data guess1RM1( 69,1)/         1.345257D0/
      data guess2RM1( 69,1)/         7.850614D0/
      data guess3RM1( 69,1)/         1.256951D0/
      data guess1RM1( 69,2)/         0.012733D0/
      data guess2RM1( 69,2)/         7.564753D0/
      data guess3RM1( 69,2)/         2.885614D0/
!
!                    Data for Element  70        Ytterbium
!
      data     ussRM1( 70)/       -21.983455D0/
      data     uppRM1( 70)/        -7.652815D0/
      data     uddRM1( 70)/       -18.071899D0/
      data   betasRM1( 70)/        -5.532949D0/
      data   betapRM1( 70)/        -0.086913D0/
      data   betadRM1( 70)/        -4.143079D0/
      data      zsRM1( 70)/         1.239808D0/
      data      zpRM1( 70)/         1.849144D0/
      data      zdRM1( 70)/         1.485378D0/
      data     zsnRM1( 70)/         1.568094D0/
      data     zpnRM1( 70)/         1.854817D0/
      data     zdnRM1( 70)/         0.749402D0/
      data     alpRM1( 70)/         1.306335D0/
      data     gssRM1( 70)/         5.965409D0/
      data     gspRM1( 70)/         6.406787D0/
      data     gppRM1( 70)/         7.740843D0/
      data     gp2RM1( 70)/         6.713836D0/
      data     hspRM1( 70)/         1.510228D0/
      data    poc_RM1( 70)/         2.608576D0/
      data    f0sdRM1( 70)/         8.369310D0/
      data    g2sdRM1( 70)/         1.261632D0/
      data guess1RM1( 70,1)/         1.319492D0/
      data guess2RM1( 70,1)/         7.587562D0/
      data guess3RM1( 70,1)/         1.517605D0/
      data guess1RM1( 70,2)/         0.025848D0/
      data guess2RM1( 70,2)/         7.882701D0/
      data guess3RM1( 70,2)/         3.232953D0/
!
!                    Data for Element  71         Lutetium
!
      data     ussRM1( 71)/       -22.032739D0/
      data     uppRM1( 71)/        -7.542271D0/
      data     uddRM1( 71)/       -18.222119D0/
      data   betasRM1( 71)/        -5.527432D0/
      data   betapRM1( 71)/        -0.244868D0/
      data   betadRM1( 71)/        -4.214400D0/
      data      zsRM1( 71)/         1.425302D0/
      data      zpRM1( 71)/         1.790353D0/
      data      zdRM1( 71)/         1.642603D0/
      data     zsnRM1( 71)/         1.475988D0/
      data     zpnRM1( 71)/         2.136483D0/
      data     zdnRM1( 71)/         0.660000D0/
      data     alpRM1( 71)/         1.434498D0/
      data     gssRM1( 71)/         5.615013D0/
      data     gspRM1( 71)/         6.372097D0/
      data     gppRM1( 71)/         8.916340D0/
      data     gp2RM1( 71)/         7.733375D0/
      data     hspRM1( 71)/         1.122642D0/
      data    poc_RM1( 71)/         2.209660D0/
      data    f0sdRM1( 71)/         8.171496D0/
      data    g2sdRM1( 71)/         1.086967D0/
      data guess1RM1( 71,1)/         0.771212D0/
      data guess2RM1( 71,1)/         7.664851D0/
      data guess3RM1( 71,1)/         1.726925D0/
      data guess1RM1( 71,2)/         0.011350D0/
      data guess2RM1( 71,2)/         7.870821D0/
      data guess3RM1( 71,2)/         3.465402D0/
!
!                    Data for Element102      Capped bond
!
      data     ussRM1(102)/     -11.9062760D0/
      data   betasRM1(102)/-9999999.0000000D0/
      data      zsRM1(102)/       4.0000000D0/
      data      zpRM1(102)/       0.3000000D0/
      data     alpRM1(102)/       2.5441341D0/
      data     gssRM1(102)/      12.8480000D0/
      data     hspRM1(102)/       0.1000000D0/
!
!       Data for Element 100:   3+ Sparkle
!
      data      alpRM1(100)/       1.5000000d0/
!
!       Data for Element 101:   3- Sparkle
!
      data      alpRM1(101)/       1.5000000d0/
!
!       Data for Element 103:   ++ Sparkle
!
      data      alpRM1(103)/       1.5000000d0/
!
!       Data for Element 104:    + Sparkle
!
      data      alpRM1(104)/       1.5000000d0/
!
!       Data for Element 105:   -- Sparkle
!
      data      alpRM1(105)/       1.5000000d0/
!
!       Data for Element 106:    - Sparkle
!
      data      alpRM1(106)/       1.5000000d0/
 end module Parameters_for_RM1_C
