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

      module Parameters_for_AM1_Sparkles_C
      double precision, dimension(107) :: alpam1sp, gssam1sp
      double precision, dimension(107,4) :: guesam1sp1, guesam1sp2, guesam1sp3
!
!       Data for Element  57:    Lanthanum
!
      data      alpam1sp( 57)/       2.1879021d0/
      data    guesam1sp1( 57,1)/       1.3207809d0/
      data    guesam1sp2( 57,1)/       7.1394307d0/
      data    guesam1sp3( 57,1)/       1.8503282d0/
      data    guesam1sp1( 57,2)/       0.3425778d0/
      data    guesam1sp2( 57,2)/       8.7780632d0/
      data    guesam1sp3( 57,2)/       3.1678964d0/
      data    gssam1sp(57) /          55.7344864002d0/
!
!       Data for Element  58:       Cerium
!
      data      alpam1sp( 58)/       2.6637770d0/
      data    guesam1sp1( 58,1)/       1.7507655d0/
      data    guesam1sp2( 58,1)/       7.6163181d0/
      data    guesam1sp3( 58,1)/       1.8064853d0/
      data    guesam1sp1( 58,2)/       0.0093401d0/
      data    guesam1sp2( 58,2)/       8.7664931d0/
      data    guesam1sp3( 58,2)/       3.2008171d0/
      data    gssam1sp(58)  /         58.7223887052d0/
!
!       Data for Element  59: Praseodymium
!
      data      alpam1sp( 59)/       2.6104230d0/
      data    guesam1sp1( 59,1)/       1.7515391d0/
      data    guesam1sp2( 59,1)/       7.6039743d0/
      data    guesam1sp3( 59,1)/       1.8084677d0/
      data    guesam1sp1( 59,2)/       0.0097057d0/
      data    guesam1sp2( 59,2)/       8.7264195d0/
      data    guesam1sp3( 59,2)/       2.9111890d0/
      data    gssam1sp(59)  /         58.9017644267d0/
!
!       Data for Element  60:    Neodymium
!
      data      alpam1sp( 60)/       4.5002951d0/
      data    guesam1sp1( 60,1)/       1.1206946d0/
      data    guesam1sp2( 60,1)/       6.8295606d0/
      data    guesam1sp3( 60,1)/       1.7859050d0/
      data    guesam1sp1( 60,2)/       0.1070369d0/
      data    guesam1sp2( 60,2)/      10.7894805d0/
      data    guesam1sp3( 60,2)/       3.1628661d0/
      data    gssam1sp(60)  /         57.6242766015d0/
!
!       Data for Element  61:   Promethium
!
      data      alpam1sp( 61)/       3.1059834d0/
      data    guesam1sp1( 61,1)/       1.7347671d0/
      data    guesam1sp2( 61,1)/       9.2464226d0/
      data    guesam1sp3( 61,1)/       1.7533419d0/
      data    guesam1sp1( 61,2)/       0.2571017d0/
      data    guesam1sp2( 61,2)/       7.8793445d0/
      data    guesam1sp3( 61,2)/       3.0498163d0/
      data    gssam1sp(61)  /         59.4249705519d0/
!
!       Data for Element  62:     Samarium
!
      data      alpam1sp( 62)/       4.1758509d0/
      data    guesam1sp1( 62,1)/       0.9592885d0/
      data    guesam1sp2( 62,1)/       6.4799924d0/
      data    guesam1sp3( 62,1)/       1.7381402d0/
      data    guesam1sp1( 62,2)/       0.0261004d0/
      data    guesam1sp2( 62,2)/       9.7391952d0/
      data    guesam1sp3( 62,2)/       2.8881177d0/
      data    gssam1sp(62)  /         56.9935144820d0/
!
!       Data for Element  63:     Europium
!
      data      alpam1sp( 63)/       2.1247189d0/
      data    guesam1sp1( 63,1)/       0.5695122d0/
      data    guesam1sp2( 63,1)/       7.4680208d0/
      data    guesam1sp3( 63,1)/       1.7319730d0/
      data    guesam1sp1( 63,2)/       0.3286619d0/
      data    guesam1sp2( 63,2)/       7.8009780d0/
      data    guesam1sp3( 63,2)/       2.9641285d0/
      data    gssam1sp(63) /          55.6059122033D0/
!
!       Data for Element  64:   Gadolinium
!
      data      alpam1sp( 64)/       3.6525485d0/
      data    guesam1sp1( 64,1)/       0.7013512d0/
      data    guesam1sp2( 64,1)/       7.5454483d0/
      data    guesam1sp3( 64,1)/       1.7761953d0/
      data    guesam1sp1( 64,2)/       0.1293094d0/
      data    guesam1sp2( 64,2)/       8.3437991d0/
      data    guesam1sp3( 64,2)/       3.0110320d0/
      data    gssam1sp(64) /          55.7083247618D0/
!
!       Data for Element  65:      Terbium
!
      data      alpam1sp( 65)/       2.3418889d0/
      data    guesam1sp1( 65,1)/       0.7734458d0/
      data    guesam1sp2( 65,1)/       7.6510526d0/
      data    guesam1sp3( 65,1)/       1.7033464d0/
      data    guesam1sp1( 65,2)/       0.3936233d0/
      data    guesam1sp2( 65,2)/       7.9261457d0/
      data    guesam1sp3( 65,2)/       3.0132951d0/
      data    gssam1sp(65) /          55.7245956904D0/
!
!       Data for Element  66:   Dysprosium
!
      data      alpam1sp( 66)/       2.4164925d0/
      data    guesam1sp1( 66,1)/       1.0385214d0/
      data    guesam1sp2( 66,1)/       8.0016177d0/
      data    guesam1sp3( 66,1)/       1.7161372d0/
      data    guesam1sp1( 66,2)/       0.3018082d0/
      data    guesam1sp2( 66,2)/       8.7917319d0/
      data    guesam1sp3( 66,2)/       2.9953821d0/
      data    gssam1sp(66)  /         55.7676495431d0/
!
!       Data for Element  67:      Holmium
!
      data      alpam1sp( 67)/       3.5558256d0/
      data    guesam1sp1( 67,1)/       0.9819488d0/
      data    guesam1sp2( 67,1)/       8.9158637d0/
      data    guesam1sp3( 67,1)/       1.7513929d0/
      data    guesam1sp1( 67,2)/       0.5454297d0/
      data    guesam1sp2( 67,2)/       8.6358144d0/
      data    guesam1sp3( 67,2)/       3.0176829d0/
      data    gssam1sp(67)  /         58.0375569983d0/
!
!       Data for Element  68:       Erbium
!
      data      alpam1sp( 68)/       3.6568233d0/
      data    guesam1sp1( 68,1)/       0.7029402d0/
      data    guesam1sp2( 68,1)/       8.7235010d0/
      data    guesam1sp3( 68,1)/       1.7746085d0/
      data    guesam1sp1( 68,2)/       0.1321262d0/
      data    guesam1sp2( 68,2)/       8.3498076d0/
      data    guesam1sp3( 68,2)/       3.0114807d0/
      data    gssam1sp(68)  /         58.0489423317d0/
!
!       Data for Element  69:      Thulium
!
      data      alpam1sp( 69)/       2.6969697d0/
      data    guesam1sp1( 69,1)/       0.9252929d0/
      data    guesam1sp2( 69,1)/       7.8636142d0/
      data    guesam1sp3( 69,1)/       1.7113099d0/
      data    guesam1sp1( 69,2)/       0.4125806d0/
      data    guesam1sp2( 69,2)/       8.5211941d0/
      data    guesam1sp3( 69,2)/       2.9817875d0/
      data    gssam1sp(69)  /         55.8847434446d0/
!
!       Data for Element  70:    Ytterbium
!
      data      alpam1sp( 70)/       4.0022391d0/
      data    guesam1sp1( 70,1)/       1.0231211d0/
      data    guesam1sp2( 70,1)/       8.3969198d0/
      data    guesam1sp3( 70,1)/       1.7046767d0/
      data    guesam1sp1( 70,2)/       0.3351965d0/
      data    guesam1sp2( 70,2)/       7.2817207d0/
      data    guesam1sp3( 70,2)/       2.9140122d0/
      data    gssam1sp(70)  /         56.1788359906d0/
!
!       Data for Element  71:     Lutetium
!
      data      alpam1sp( 71)/       4.0203424d0/
      data    guesam1sp1( 71,1)/       1.0381639d0/
      data    guesam1sp2( 71,1)/       8.4911797d0/
      data    guesam1sp3( 71,1)/       1.7034421d0/
      data    guesam1sp1( 71,2)/       0.3342233d0/
      data    guesam1sp2( 71,2)/       7.2729947d0/
      data    guesam1sp3( 71,2)/       2.9153096d0/
      data    gssam1sp(71)  /         56.1751741742d0/

   end module Parameters_for_AM1_Sparkles_C
