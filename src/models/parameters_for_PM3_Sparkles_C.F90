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

      module Parameters_for_PM3_Sparkles_C
      double precision, dimension(107) :: alpPM3sp, gssPM3sp
      double precision, dimension(107,4) :: guesPM3sp1, guesPM3sp2, guesPM3sp3
!
!       Data for Element  57:    Lanthanum
!
      data    alpPM3sp( 57)/          2.0790327441D0/
      data    gssPM3sp(57) /          55.7063938597D0/
      data    guesPM3sp1( 57,1)/       0.9478615182D0/
      data    guesPM3sp2( 57,1)/       7.2366011717D0/
      data    guesPM3sp3( 57,1)/       1.8543786298D0/
      data    guesPM3sp1( 57,2)/       0.3177671904D0/
      data    guesPM3sp2( 57,2)/       8.5224006292D0/
      data    guesPM3sp3( 57,2)/       3.0752474724D0/

!
!       Data for Element  58:       Cerium
!
      data    alpPM3sp( 58)/           2.5665085968D0/
      data    gssPM3sp(58)  /         58.5701153062D0/
      data    guesPM3sp1( 58,1)/       1.8026688761D0/
      data    guesPM3sp2( 58,1)/       7.5971870028D0/
      data    guesPM3sp3( 58,1)/       1.8009003439D0/
      data    guesPM3sp1( 58,2)/       0.1319892158D0/
      data    guesPM3sp2( 58,2)/       9.6116040841D0/
      data    guesPM3sp3( 58,2)/       3.0613741124D0/

!
!       Data for Element  59: Praseodymium
!
      data    alpPM3sp( 59)/           2.6002012752D0/
      data    gssPM3sp(59)  /         58.9488361614D0/
      data    guesPM3sp1( 59,1)/       1.7432011456D0/
      data    guesPM3sp2( 59,1)/       8.1376077067D0/
      data    guesPM3sp3( 59,1)/       1.8084052715D0/
      data    guesPM3sp1( 59,2)/       0.0112202332D0/
      data    guesPM3sp2( 59,2)/       8.7388958526D0/
      data    guesPM3sp3( 59,2)/       2.9806781047D0/

!
!       Data for Element  60:    Neodymium
!
      data    alpPM3sp( 60)/           4.7057677595D0/
      data    gssPM3sp(60)  /         57.4944898977D0/
      data    guesPM3sp1( 60,1)/       1.0715972265D0/
      data    guesPM3sp2( 60,1)/       6.9565346287D0/
      data    guesPM3sp3( 60,1)/       1.7812099249D0/
      data    guesPM3sp1( 60,2)/       0.0886417116D0/
      data    guesPM3sp2( 60,2)/      10.8664473398D0/
      data    guesPM3sp3( 60,2)/       3.0992613820D0/

!
!       Data for Element  61:   Promethium
!
      data    alpPM3sp( 61)/           3.1490918074D0/
      data    gssPM3sp(61)  /         59.2924444913D0/
      data    guesPM3sp1( 61,1)/       1.6572814674D0/
      data    guesPM3sp2( 61,1)/       9.2529413759D0/
      data    guesPM3sp3( 61,1)/       1.7412637448D0/
      data    guesPM3sp1( 61,2)/       0.1851223683D0/
      data    guesPM3sp2( 61,2)/       7.4186533283D0/
      data    guesPM3sp3( 61,2)/       3.0623727738D0/

!
!       Data for Element  62:     Samarium
!
      data    alpPM3sp( 62)/           3.6813938335D0/
      data    gssPM3sp(62)  /         54.8086404668D0/
      data    guesPM3sp1( 62,1)/       0.7706615984D0/
      data    guesPM3sp2( 62,1)/       6.6020324700D0/
      data    guesPM3sp3( 62,1)/       1.7636673188D0/
      data    guesPM3sp1( 62,2)/       0.0936188340D0/
      data    guesPM3sp2( 62,2)/       9.3136737687D0/
      data    guesPM3sp3( 62,2)/       2.9879390071D0/

!
!       Data for Element  63:     Europium
!
      data    alpPM3sp( 63)/           2.1398139884D0/
      data    gssPM3sp(63) /          55.5863246694D0/
      data    guesPM3sp1( 63,1)/       0.6101627168D0/
      data    guesPM3sp2( 63,1)/       7.1373146362D0/
      data    guesPM3sp3( 63,1)/       1.7807085112D0/
      data    guesPM3sp1( 63,2)/       0.3415714636D0/
      data    guesPM3sp2( 63,2)/       9.1732778046D0/
      data    guesPM3sp3( 63,2)/       3.0121099267D0/

!
!       Data for Element  64:   Gadolinium
!
      data    alpPM3sp( 64)/           3.6813938335D0/
      data    gssPM3sp(64) /          54.8086404668D0/
      data    guesPM3sp1( 64,1)/       0.7706615984D0/
      data    guesPM3sp2( 64,1)/       7.5453068267D0/
      data    guesPM3sp3( 64,1)/       1.7636673188D0/
      data    guesPM3sp1( 64,2)/       0.0936188340D0/
      data    guesPM3sp2( 64,2)/       8.2224517067D0/
      data    guesPM3sp3( 64,2)/       2.9879390071D0/
!
!       Data for Element  65:      Terbium
!
      data    alpPM3sp( 65)/           2.8245126194D0/
      data    gssPM3sp(65) /          56.2564137683D0/
      data    guesPM3sp1( 65,1)/       1.3428294115D0/
      data    guesPM3sp2( 65,1)/       7.5782265384D0/
      data    guesPM3sp3( 65,1)/       1.7181508908D0/
      data    guesPM3sp1( 65,2)/       0.2651000290D0/
      data    guesPM3sp2( 65,2)/       6.4476118233D0/
      data    guesPM3sp3( 65,2)/       2.9952711306D0/

!
!       Data for Element  66:   Dysprosium
!
      data    alpPM3sp( 66)/           2.4630183002D0/
      data    gssPM3sp(66)  /         55.7563629021D0/
      data    guesPM3sp1( 66,1)/       1.3287435702D0/
      data    guesPM3sp2( 66,1)/       7.9816784235D0/
      data    guesPM3sp3( 66,1)/       1.7080927380D0/
      data    guesPM3sp1( 66,2)/       0.3332412720D0/
      data    guesPM3sp2( 66,2)/       9.7816147381D0/
      data    guesPM3sp3( 66,2)/       2.9323938165D0/

!
!       Data for Element  67:      Holmium
!
      data    alpPM3sp( 67)/           3.7240820504D0/
      data    gssPM3sp(67)  /         58.0161995449D0/
      data    guesPM3sp1( 67,1)/       1.0370243923D0/
      data    guesPM3sp2( 67,1)/       8.7236664457D0/
      data    guesPM3sp3( 67,1)/       1.7412327388D0/
      data    guesPM3sp1( 67,2)/       0.5175396118D0/
      data    guesPM3sp2( 67,2)/      10.6247030904D0/
      data    guesPM3sp3( 67,2)/       3.0090983082D0/

!
!       Data for Element  68:       Erbium
!
      data    alpPM3sp( 68)/           3.6824732031D0/
      data    gssPM3sp(68)  /         58.0537274879D0/
      data    guesPM3sp1( 68,1)/       0.5361755749D0/
      data    guesPM3sp2( 68,1)/       8.7390296713D0/
      data    guesPM3sp3( 68,1)/       1.7856436299D0/
      data    guesPM3sp1( 68,2)/       0.0776867913D0/
      data    guesPM3sp2( 68,2)/       8.6267120701D0/
      data    guesPM3sp3( 68,2)/       2.9875549001D0/

!
!       Data for Element  69:      Thulium
!
      data    alpPM3sp( 69)/           2.7914853437D0/
      data    gssPM3sp(69)  /         55.8262787893D0/
      data    guesPM3sp1( 69,1)/       1.0428337088D0/
      data    guesPM3sp2( 69,1)/       7.6146043282D0/
      data    guesPM3sp3( 69,1)/       1.7354274982D0/
      data    guesPM3sp1( 69,2)/       0.4554143900D0/
      data    guesPM3sp2( 69,2)/       8.5730264221D0/
      data    guesPM3sp3( 69,2)/       2.9246089938D0/

!
!       Data for Element  70:    Ytterbium
!
      data    alpPM3sp( 70)/           3.9448800504D0/
      data    gssPM3sp(70)  /         55.8087638754D0/
      data    guesPM3sp1( 70,1)/       1.4143797707D0/
      data    guesPM3sp2( 70,1)/       8.3678869738D0/
      data    guesPM3sp3( 70,1)/       1.7015172798D0/
      data    guesPM3sp1( 70,2)/       0.1969092096D0/
      data    guesPM3sp2( 70,2)/       7.5001690472D0/
      data    guesPM3sp3( 70,2)/       2.8840018666D0/

!
!       Data for Element  71:     Lutetium
!
      data    alpPM3sp( 71)/           4.0585822429D0/
      data    gssPM3sp(71)  /         56.6177268844D0/
      data    guesPM3sp1( 71,1)/       1.1216536112D0/
      data    guesPM3sp2( 71,1)/       7.9408016011D0/
      data    guesPM3sp3( 71,1)/       1.7060269731D0/
      data    guesPM3sp1( 71,2)/       0.3305655101D0/
      data    guesPM3sp2( 71,2)/       7.3013810219D0/
      data    guesPM3sp3( 71,2)/       2.9090465178D0/

      end module Parameters_for_PM3_Sparkles_C
