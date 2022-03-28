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

      module Parameters_for_MNDOD_C
      double precision, dimension(107) :: ussd, uppd, zsd, zpd, betasd, betapd, &
        alpd, gssd, gppd, gspd, gp2d, hspd, uddd, zdd, betadd, zsnd, zpnd, zdnd&
        , poc_d
!
!                    Data for Element  1         Hydrogen
!
      data     ussd(  1)/     -11.9062760D0/
      data   betasd(  1)/      -6.9890640D0/
      data      zsd(  1)/       1.3319670D0/
      data     zsnd(  1)/       0.0000000D0/
      data     alpd(  1)/       2.5441341D0/
      data     gssd(  1)/      12.8480000D0/
!
!                    Data for Element  2           Helium
!
      data     ussd(  2)/     -35.4931170D0/
      data     uppd(  2)/       9.9998433D0/
      data   betasd(  2)/     -17.5545020D0/
      data   betapd(  2)/     -26.0526064D0/
      data      zsd(  2)/       1.7710761D0/
      data      zpd(  2)/       6.9018258D0/
      data     zsnd(  2)/       0.0000000D0/
      data     alpd(  2)/       5.7813931D0/
      data     gssd(  2)/       8.0590135D0/
      data     gspd(  2)/      11.0152229D0/
      data     gppd(  2)/      19.3843339D0/
      data     gp2d(  2)/      11.3517224D0/
      data     hspd(  2)/       0.5321205D0/
!
!                    Data for Element  3          Lithium
!
      data     ussd(  3)/      -4.8578570D0/
      data     uppd(  3)/      -2.0266084D0/
      data   betasd(  3)/      -0.1904126D0/
      data   betapd(  3)/      -1.6081355D0/
      data      zsd(  3)/       0.4296141D0/
      data      zpd(  3)/       0.7554884D0/
      data     zsnd(  3)/       0.0000000D0/
      data     alpd(  3)/       1.2083987D0/
      data     gssd(  3)/       7.5947069D0/
      data     gspd(  3)/       6.7259856D0/
      data     gppd(  3)/       8.6596829D0/
      data     gp2d(  3)/       3.8714751D0/
      data     hspd(  3)/       5.0003381D0/
!
!                    Data for Element  4        Beryllium
!
      data     ussd(  4)/     -16.6023780D0/
      data     uppd(  4)/     -10.7037710D0/
      data   betasd(  4)/      -4.0170960D0/
      data   betapd(  4)/      -4.0170960D0/
      data      zsd(  4)/       1.0042100D0/
      data      zpd(  4)/       1.0042100D0/
      data     zsnd(  4)/       0.0000000D0/
      data     alpd(  4)/       1.6694340D0/
      data     gssd(  4)/       9.0000000D0/
      data     gspd(  4)/       7.4300000D0/
      data     gppd(  4)/       6.9700000D0/
      data     gp2d(  4)/       6.2200000D0/
      data     hspd(  4)/       1.2800000D0/
!
!                    Data for Element  5            Boron
!
      data     ussd(  5)/     -34.5471300D0/
      data     uppd(  5)/     -23.1216900D0/
      data   betasd(  5)/      -8.2520540D0/
      data   betapd(  5)/      -8.2520540D0/
      data      zsd(  5)/       1.5068010D0/
      data      zpd(  5)/       1.5068010D0/
      data     zsnd(  5)/       0.0000000D0/
      data     alpd(  5)/       2.1349930D0/
      data     gssd(  5)/      10.5900000D0/
      data     gspd(  5)/       9.5600000D0/
      data     gppd(  5)/       8.8600000D0/
      data     gp2d(  5)/       7.8600000D0/
      data     hspd(  5)/       1.8100000D0/
!
!                    Data for Element  6           Carbon
!
      data     ussd(  6)/     -52.2797450D0/
      data     uppd(  6)/     -39.2055580D0/
      data   betasd(  6)/     -18.9850440D0/
      data   betapd(  6)/      -7.9341220D0/
      data      zsd(  6)/       1.7875370D0/
      data      zpd(  6)/       1.7875370D0/
      data     zsnd(  6)/       0.0000000D0/
      data     alpd(  6)/       2.5463800D0/
      data     gssd(  6)/      12.2300000D0/
      data     gspd(  6)/      11.4700000D0/
      data     gppd(  6)/      11.0800000D0/
      data     gp2d(  6)/       9.8400000D0/
      data     hspd(  6)/       2.4300000D0/
!
!                    Data for Element  7         Nitrogen
!
      data     ussd(  7)/     -71.9321220D0/
      data     uppd(  7)/     -57.1723190D0/
      data   betasd(  7)/     -20.4957580D0/
      data   betapd(  7)/     -20.4957580D0/
      data      zsd(  7)/       2.2556140D0/
      data      zpd(  7)/       2.2556140D0/
      data     zsnd(  7)/       0.0000000D0/
      data     alpd(  7)/       2.8613420D0/
      data     gssd(  7)/      13.5900000D0/
      data     gspd(  7)/      12.6600000D0/
      data     gppd(  7)/      12.9800000D0/
      data     gp2d(  7)/      11.5900000D0/
      data     hspd(  7)/       3.1400000D0/
!
!                    Data for Element  8           Oxygen
!
      data     ussd(  8)/     -99.6443090D0/
      data     uppd(  8)/     -77.7974720D0/
      data   betasd(  8)/     -32.6880820D0/
      data   betapd(  8)/     -32.6880820D0/
      data      zsd(  8)/       2.6999050D0/
      data      zpd(  8)/       2.6999050D0/
      data     zsnd(  8)/       0.0000000D0/
      data     alpd(  8)/       3.1606040D0/
      data     gssd(  8)/      15.4200000D0/
      data     gspd(  8)/      14.4800000D0/
      data     gppd(  8)/      14.5200000D0/
      data     gp2d(  8)/      12.9800000D0/
      data     hspd(  8)/       3.9400000D0/
!
!                    Data for Element  9         Fluorine
!
      data     ussd(  9)/    -131.0715480D0/
      data     uppd(  9)/    -105.7821370D0/
      data   betasd(  9)/     -48.2904660D0/
      data   betapd(  9)/     -36.5085400D0/
      data      zsd(  9)/       2.8484870D0/
      data      zpd(  9)/       2.8484870D0/
      data     zsnd(  9)/       0.0000000D0/
      data     alpd(  9)/       3.4196606D0/
      data     gssd(  9)/      16.9200000D0/
      data     gspd(  9)/      17.2500000D0/
      data     gppd(  9)/      16.7100000D0/
      data     gp2d(  9)/      14.9100000D0/
      data     hspd(  9)/       4.8300000D0/
!
!                    Data for Element 10             Neon
!
      data     ussd( 10)/       9.6554736D0/
      data     uppd( 10)/     -71.1466474D0/
      data   betasd( 10)/      -0.1512657D0/
      data   betapd( 10)/     -23.8470187D0/
      data      zsd( 10)/       5.9998745D0/
      data      zpd( 10)/       4.1752600D0/
      data     zsnd( 10)/       0.0000000D0/
      data     alpd( 10)/       2.5377255D0/
      data     gssd( 10)/       0.4989374D0/
      data     gspd( 10)/      10.1592044D0/
      data     gppd( 10)/      18.9449976D0/
      data     gp2d( 10)/       8.5648589D0/
      data     hspd( 10)/       0.2999991D0/
!
!                    Data for Element 11           Sodium
!
      data     ussd( 11)/      -5.2010000D0/
      data     uppd( 11)/      -2.7125732D0/
      data   betasd( 11)/      -1.0873817D0/
      data   betapd( 11)/      -0.4862394D0/
      data      zsd( 11)/       0.9875083D0/
      data      zpd( 11)/       0.8933498D0/
      data     zsnd( 11)/       0.6541126D0/
      data     zpnd( 11)/       0.5644087D0/
      data     alpd( 11)/       1.1701020D0/
      data     gssd( 11)/       4.5944448D0/
      data     gspd( 11)/       4.1475740D0/
      data     gppd( 11)/       4.2991976D0/
      data     gp2d( 11)/       3.7969573D0/
      data     hspd( 11)/       0.5344087D0/
      data    poc_d( 11)/       1.5305533D0/
!
!                    Data for Element 12        Magnesium
!
      data     ussd( 12)/     -15.0970000D0/
      data     uppd( 12)/     -10.6500000D0/
      data   betasd( 12)/      -1.8958835D0/
      data   betapd( 12)/      -2.1410894D0/
      data      zsd( 12)/       1.4489045D0/
      data      zpd( 12)/       0.9529300D0/
      data     zsnd( 12)/       1.0500000D0/
      data     zpnd( 12)/       0.9252719D0/
      data     alpd( 12)/       1.6214698D0/
      data     gssd( 12)/       7.3751326D0/
      data     gspd( 12)/       6.8889074D0/
      data     gppd( 12)/       7.0479538D0/
      data     gp2d( 12)/       6.2245987D0/
      data     hspd( 12)/       0.7267339D0/
      data    poc_d( 12)/       1.3507762D0/
!
!                    Data for Element 13        Aluminium
!
      data     ussd( 13)/     -24.0179291D0/
      data     uppd( 13)/     -20.7959797D0/
      data     uddd( 13)/      -5.2208274D0/
      data   betasd( 13)/      -7.1018585D0/
      data   betapd( 13)/      -2.3180962D0/
      data   betadd( 13)/      -3.3563854D0/
      data      zsd( 13)/       1.7940227D0/
      data      zpd( 13)/       1.3713092D0/
      data      zdd( 13)/       0.8059113D0/
      data     zsnd( 13)/       1.2224927D0/
      data     zpnd( 13)/       1.0929199D0/
      data     zdnd( 13)/       0.8003828D0/
      data     alpd( 13)/       1.4430168D0/
      data     gssd( 13)/       8.5867102D0/
      data     gspd( 13)/       7.6646931D0/
      data     gppd( 13)/       8.3249572D0/
      data     gp2d( 13)/       7.3524202D0/
      data     hspd( 13)/       0.5440129D0/
      data    poc_d( 13)/       1.5844252D0/
!
!                    Data for Element 14          Silicon
!
      data     ussd( 14)/     -36.0515300D0/
      data     uppd( 14)/     -27.5356910D0/
      data     uddd( 14)/     -14.6774390D0/
      data   betasd( 14)/      -8.2107342D0/
      data   betapd( 14)/      -4.8846203D0/
      data   betadd( 14)/      -2.6080115D0/
      data      zsd( 14)/       1.9156546D0/
      data      zpd( 14)/       1.6816113D0/
      data     zdd( 14)/       0.9667717D0/
      data     zsnd( 14)/       1.5292918D0/
      data     zpnd( 14)/       0.9762808D0/
      data     zdnd( 14)/       0.9381644D0/
      data     alpd( 14)/       1.6600693D0/
      data     gssd( 14)/      10.7416470D0/
      data     gspd( 14)/       7.5606664D0/
      data     gppd( 14)/       7.4364969D0/
      data     gp2d( 14)/       6.5677515D0/
      data     hspd( 14)/       0.8775388D0/
      data    poc_d( 14)/       1.2665655D0/
!
!                    Data for Element 15       Phosphorus
!
      data     ussd( 15)/     -47.0555290D0/
      data     uppd( 15)/     -38.0670590D0/
      data     uddd( 15)/     -23.6915970D0/
      data   betasd( 15)/      -8.9021043D0/
      data   betapd( 15)/      -9.3861108D0/
      data   betadd( 15)/      -2.0917008D0/
      data      zsd( 15)/       2.2664629D0/
      data      zpd( 15)/       1.9400149D0/
      data     zdd( 15)/       1.1001090D0/
      data     zsnd( 15)/       1.6343761D0/
      data     zpnd( 15)/       1.0829117D0/
      data     zdnd( 15)/       1.0065147D0/
      data     alpd( 15)/       1.8525512D0/
      data     gssd( 15)/      11.4797530D0/
      data     gspd( 15)/       8.5575691D0/
      data     gppd( 15)/       8.2487228D0/
      data     gp2d( 15)/       7.2850917D0/
      data     hspd( 15)/       2.1078044D0/
      data    poc_d( 15)/       1.1851300D0/
!
!                    Data for Element 16           Sulfur
!
      data     ussd( 16)/     -56.8891280D0/
      data     uppd( 16)/     -47.2747450D0/
      data     uddd( 16)/     -25.0951180D0/
      data   betasd( 16)/     -10.9995450D0/
      data   betapd( 16)/     -12.2154370D0/
      data   betadd( 16)/      -1.8806695D0/
      data      zsd( 16)/       2.2258505D0/
      data      zpd( 16)/       2.0997056D0/
      data     zdd( 16)/       1.2314725D0/
      data     zsnd( 16)/       1.7363914D0/
      data     zpnd( 16)/       1.1211817D0/
      data     zdnd( 16)/       1.0508467D0/
      data     alpd( 16)/       2.0230595D0/
      data     gssd( 16)/      12.1963020D0/
      data     gspd( 16)/       8.8539009D0/
      data     gppd( 16)/       8.5402324D0/
      data     gp2d( 16)/       7.5425465D0/
      data     hspd( 16)/       2.6463523D0/
      data    poc_d( 16)/       1.1155021D0/
!
!                    Data for Element 17         Chlorine
!
      data     ussd( 17)/     -69.6229728D0/
      data     uppd( 17)/     -59.1007290D0/
      data     uddd( 17)/     -36.6745732D0/
      data   betasd( 17)/      -6.0372917D0/
      data   betapd( 17)/     -19.1833850D0/
      data   betadd( 17)/      -1.8777820D0/
      data      zsd( 17)/       2.5616106D0/
      data      zpd( 17)/       2.3893380D0/
      data     zdd( 17)/       1.2513978D0/
      data     zsnd( 17)/       1.8808755D0/
      data     zpnd( 17)/       1.1810423D0/
      data     zdnd( 17)/       1.1406155D0/
      data     alpd( 17)/       2.1803002D0/
      data     gssd( 17)/      13.2111485D0/
      data     gspd( 17)/       9.4194951D0/
      data     gppd( 17)/       8.9962003D0/
      data     gp2d( 17)/       7.9452475D0/
      data     hspd( 17)/       3.0814986D0/
      data    poc_d( 17)/       1.0298121D0/
!
!                    Data for Element 30             Zinc
!
      data     ussd( 30)/     -18.0230000D0/
      data     uppd( 30)/     -12.2421658D0/
      data   betasd( 30)/      -5.0172608D0/
      data   betapd( 30)/      -0.7120597D0/
      data      zsd( 30)/       1.7315035D0/
      data      zpd( 30)/       1.3935830D0/
      data     zsnd( 30)/       1.5660000D0/
      data     zpnd( 30)/       0.8628398D0/
      data     alpd( 30)/       1.5176370D0/
      data     gssd( 30)/       8.5607284D0/
      data     gspd( 30)/       7.4900360D0/
      data     gppd( 30)/       5.1396483D0/
      data     gp2d( 30)/       4.5054031D0/
      data     hspd( 30)/       0.5329461D0/
      data    poc_d( 30)/       1.5892339D0/
!
!                    Data for Element 35          Bromine
!
      data     ussd( 35)/     -65.4027779D0/
      data     uppd( 35)/     -54.5537535D0/
      data     uddd( 35)/     -13.7280993D0/
      data   betasd( 35)/      -8.3149761D0/
      data   betapd( 35)/     -10.5070414D0/
      data   betadd( 35)/      -0.9625993D0/
      data      zsd( 35)/       2.5905410D0/
      data      zpd( 35)/       2.3308565D0/
      data     zdd( 35)/       1.3573612D0/
      data     zsnd( 35)/       2.2358154D0/
      data     zpnd( 35)/       1.4329265D0/
      data     zdnd( 35)/       1.2425783D0/
      data     alpd( 35)/       2.0910500D0/
      data     gssd( 35)/      12.2223555D0/
      data     gspd( 35)/       8.2637201D0/
      data     gppd( 35)/       8.5354644D0/
      data     gp2d( 35)/       7.4821671D0/
      data     hspd( 35)/       2.7495223D0/
      data    poc_d( 35)/       1.1131242D0/
!
!                    Data for Element 48          Cadmium
!
      data     ussd( 48)/     -16.9697000D0/
      data     uppd( 48)/     -12.4009648D0/
      data   betasd( 48)/      -2.7715438D0/
      data   betapd( 48)/      -1.8056502D0/
      data      zsd( 48)/       1.7488056D0/
      data      zpd( 48)/       1.5632147D0/
      data     zsnd( 48)/       1.7631484D0/
      data     zpnd( 48)/       1.5255190D0/
      data     alpd( 48)/       1.4246133D0/
      data     gssd( 48)/       7.9044344D0/
      data     gspd( 48)/       7.5157069D0/
      data     gppd( 48)/       7.4799999D0/
      data     gp2d( 48)/       6.5186642D0/
      data     hspd( 48)/       0.6367444D0/
      data    poc_d( 48)/       1.7211858D0/
!
!                    Data for Element 53           Iodine
!
      data     ussd( 53)/     -62.7653526D0/
      data     uppd( 53)/     -50.2921157D0/
      data     uddd( 53)/     -12.2483050D0/
      data   betasd( 53)/     -10.6994867D0/
      data   betapd( 53)/      -4.9411778D0/
      data   betadd( 53)/      -2.3504610D0/
      data      zsd( 53)/       2.7565432D0/
      data      zpd( 53)/       2.2530795D0/
      data     zdd( 53)/       1.5023351D0/
      data     zsnd( 53)/       2.6724110D0/
      data     zpnd( 53)/       1.5722987D0/
      data     zdnd( 53)/       1.2588480D0/
      data     alpd( 53)/       1.9061744D0/
      data     gssd( 53)/      11.9807820D0/
      data     gspd( 53)/       7.8559019D0/
      data     gppd( 53)/       7.7093723D0/
      data     gp2d( 53)/       6.7185573D0/
      data     hspd( 53)/       2.0714746D0/
      data    poc_d( 53)/       1.1355686D0/
!
!                    Data for Element 80          Mercury
!
      data     ussd( 80)/     -18.8156490D0/
      data     uppd( 80)/     -13.3971135D0/
      data   betasd( 80)/      -2.2187224D0/
      data   betapd( 80)/      -2.9097857D0/
      data      zsd( 80)/       2.3331076D0/
      data      zpd( 80)/       1.7083107D0/
      data     zsnd( 80)/       2.1860001D0/
      data     zpnd( 80)/       1.7050046D0/
      data     alpd( 80)/       1.3822417D0/
      data     gssd( 80)/       8.3156495D0/
      data     gspd( 80)/       8.2121730D0/
      data     gppd( 80)/       7.1152588D0/
      data     gp2d( 80)/       6.1712498D0/
      data     hspd( 80)/       0.8359410D0/
      data    poc_d( 80)/       1.6360718D0/
      end module Parameters_for_MNDOD_C
