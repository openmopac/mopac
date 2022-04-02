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

  module Parameters_for_PM6_Sparkles_C
    double precision, dimension(107) ::  gss6sp, alp6sp
    double precision, dimension(107,2) :: gues6sp1, gues6sp2, gues6sp3
!
!       Data for Element  57:                       Lanthanum
!
      data     alp6sp( 57)/         2.0955474333D0/
      data     gss6sp( 57)/        55.7614959637D0/
      data gues6sp1( 57,1)/         0.9198962192D0/
      data gues6sp2( 57,1)/         7.1956586116D0/
      data gues6sp3( 57,1)/         1.8688421745D0/
      data gues6sp1( 57,2)/         0.3395617280D0/
      data gues6sp2( 57,2)/         8.5194840290D0/
      data gues6sp3( 57,2)/         3.1236739454D0/

!                    Data for Element  58           Cerium
!
      data     alp6sp( 58)/         2.1249588196D0/
      data     gss6sp( 58)/        58.8260906171D0/
      data gues6sp1( 58,1)/         1.7329167054D0/
      data gues6sp2( 58,1)/         7.4140930052D0/
      data gues6sp3( 58,1)/         1.7149281546D0/
      data gues6sp1( 58,2)/         0.0764294472D0/
      data gues6sp2( 58,2)/         8.4974750829D0/
      data gues6sp3( 58,2)/         3.0778367381D0/

!
!       Data for Element  59:                       Praseodymium
!
      data     alp6sp( 59)/         2.4693712260D0/
      data     gss6sp( 59)/        58.3343604075D0/
      data gues6sp1( 59,1)/         2.8321015232D0/
      data gues6sp2( 59,1)/         7.1195524904D0/
      data gues6sp3( 59,1)/         1.6208236553D0/
      data gues6sp1( 59,2)/         0.0541169724D0/
      data gues6sp2( 59,2)/         7.8230014741D0/
      data gues6sp3( 59,2)/         3.1133411960D0/

!
!       Data for Element  60:                       Neodymium
!
      data     alp6sp( 60)/          4.1738480733D0/
      data     gss6sp( 60)/         57.6974644976D0/
      data gues6sp1( 60,1)/          1.1507966088D0/
      data gues6sp2( 60,1)/          6.4949658378D0/
      data gues6sp3( 60,1)/          1.5653255583D0/
      data gues6sp1( 60,2)/          0.1889516026D0/
      data gues6sp2( 60,2)/         10.9231117908D0/
      data gues6sp3( 60,2)/          3.0169407149D0/

!
!       Data for Element  61:                       Promethium
!
      data     alp6sp( 61)/          3.0374070006D0/
      data     gss6sp( 61)/         59.5665725491D0/
      data gues6sp1( 61,1)/          1.8134017776D0/
      data gues6sp2( 61,1)/          9.0994056545D0/
      data gues6sp3( 61,1)/          1.6148716177D0/
      data gues6sp1( 61,2)/          0.2759193756D0/
      data gues6sp2( 61,2)/          7.2120871121D0/
      data gues6sp3( 61,2)/          3.0287226366D0/

!
!       Data for Element  62:                       Samarium
!
      data     alp6sp( 62)/          4.0858458124D0/
      data     gss6sp( 62)/         56.8573294165D0/
      data gues6sp1( 62,1)/          1.5645679440D0/
      data gues6sp2( 62,1)/          6.4255324886D0/
      data gues6sp3( 62,1)/          1.4885991013D0/
      data gues6sp1( 62,2)/          0.1021969444D0/
      data gues6sp2( 62,2)/          9.4102061689D0/
      data gues6sp3( 62,2)/          3.1094973204D0/

!
!       Data for Element  63:                       Europium
!
      data     alp6sp( 63)/          2.0467722838D0/
      data     gss6sp( 63)/         55.6632255486D0/
      data gues6sp1( 63,1)/          0.2712333739D0/
      data gues6sp2( 63,1)/          7.3743656586D0/
      data gues6sp3( 63,1)/          1.7955662564D0/
      data gues6sp1( 63,2)/          0.3493713916D0/
      data gues6sp2( 63,2)/          7.7881047906D0/
      data gues6sp3( 63,2)/          2.9632616015D0/

!
!       Data for Element  64:                       Gadolinium
!
      data     alp6sp( 64)/          2.1346333468D0/
      data     gss6sp( 64)/         56.8944696903D0/
      data gues6sp1( 64,1)/          0.2517865588D0/
      data gues6sp2( 64,1)/          8.7505991931D0/
      data gues6sp3( 64,1)/          1.7313405711D0/
      data gues6sp1( 64,2)/          0.1221903028D0/
      data gues6sp2( 64,2)/          7.4979582981D0/
      data gues6sp3( 64,2)/          2.9344373061D0/

!
!       Data for Element  65:                       Terbium
!
      data     alp6sp( 65)/          2.5139941133D0/
      data     gss6sp( 65)/         55.2205687662D0/
      data gues6sp1( 65,1)/          0.5222813525D0/
      data gues6sp2( 65,1)/          7.9527873623D0/
      data gues6sp3( 65,1)/          1.7550018623D0/
      data gues6sp1( 65,2)/          0.3099626210D0/
      data gues6sp2( 65,2)/          6.6812787003D0/
      data gues6sp3( 65,2)/          2.9759649920D0/

!
!       Data for Element  66:                       Dysprosium
!
      data     alp6sp( 66)/          2.5510632015D0/
      data     gss6sp( 66)/         55.8786332882D0/
      data gues6sp1( 66,1)/          1.1809130487D0/
      data gues6sp2( 66,1)/          8.9849822704D0/
      data gues6sp3( 66,1)/          1.6756952638D0/
      data gues6sp1( 66,2)/          0.4066395540D0/
      data gues6sp2( 66,2)/          8.9799453811D0/
      data gues6sp3( 66,2)/          2.9787279400D0/

!
!       Data for Element  67:                       Holmium
!
      data     alp6sp( 67)/          3.4814819284D0/
      data     gss6sp( 67)/         56.0010800433D0/
      data gues6sp1( 67,1)/          0.3389541104D0/
      data gues6sp2( 67,1)/          8.1824420999D0/
      data gues6sp3( 67,1)/          1.6446707189D0/
      data gues6sp1( 67,2)/          0.1333201849D0/
      data gues6sp2( 67,2)/          8.7112042124D0/
      data gues6sp3( 67,2)/          2.9809112221D0/

!
!       Data for Element  68:                       Erbium
!
      data     alp6sp( 68)/          3.6603230421D0/
      data     gss6sp( 68)/         58.4870426290D0/
      data gues6sp1( 68,1)/          0.4687052850D0/
      data gues6sp2( 68,1)/          9.3819581436D0/
      data gues6sp3( 68,1)/          1.7306657473D0/
      data gues6sp1( 68,2)/          0.2107436963D0/
      data gues6sp2( 68,2)/          8.4256041357D0/
      data gues6sp3( 68,2)/          2.7714227710D0/

!
!       Data for Element  69:                       Thulium
!
      data     alp6sp( 69)/          2.3042905227D0/
      data     gss6sp( 69)/         56.3699484190D0/
      data gues6sp1( 69,1)/          0.7757838661D0/
      data gues6sp2( 69,1)/          8.3570694122D0/
      data gues6sp3( 69,1)/          1.6489766048D0/
      data gues6sp1( 69,2)/          0.2905574744D0/
      data gues6sp2( 69,2)/          7.6933919381D0/
      data gues6sp3( 69,2)/          2.9316355211D0/

!
!       Data for Element  70:                       Ytterbium
!
      data     alp6sp( 70)/          4.2104920412D0/
      data     gss6sp( 70)/         56.3592753390D0/
      data gues6sp1( 70,1)/          1.0542080628D0/
      data gues6sp2( 70,1)/          8.5454710978D0/
      data gues6sp3( 70,1)/          1.4993488570D0/
      data gues6sp1( 70,2)/          0.1983232376D0/
      data gues6sp2( 70,2)/          8.4702758246D0/
      data gues6sp3( 70,2)/          2.8575372636D0/

!
!       Data for Element  71:                       Lutetium
!
      data     alp6sp( 71)/          3.2076166779D0/
      data     gss6sp( 71)/         56.2871833842D0/
      data gues6sp1( 71,1)/          0.6496316332D0/
      data gues6sp2( 71,1)/          9.2468459960D0/
      data gues6sp3( 71,1)/          1.5344631779D0/
      data gues6sp1( 71,2)/          0.2355401411D0/
      data gues6sp2( 71,2)/          7.3208383861D0/
      data gues6sp3( 71,2)/          2.9270906286D0/
    end module Parameters_for_PM6_Sparkles_C
