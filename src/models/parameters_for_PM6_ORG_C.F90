! Molecular Orbital PACkage (MOPAC)
! Copyright 2021 Virginia Polytechnic Institute and State University
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!    http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

  module Parameters_for_PM6_ORG_C
    double precision, dimension(107) :: uss_org, upp_org, udd_org, zs_org, zp_org, zd_org, betas_org, &
    betap_org, betad_org, gss_org, gsp_org, gpp_org, gp2_org, hsp_org, polvo_org, poc__org, &
    zsn_org, zpn_org, zdn_org, f0sd_org, g2sd_org, alp_org, &
    CPE_Zet_org, CPE_Z0_org, CPE_B_org, CPE_Xlo_org, CPE_Xhi_org
    double precision :: v_par_org(60)
    double precision, dimension(107,4) :: gues_org1, gues_org2, gues_org3
!
!                    Data for Element   1         Hydrogen
!
      data     uss_org(  1)/     -12.55594289D0/
      data   betas_org(  1)/      -8.78462079D0/
      data      zs_org(  1)/       1.22600545D0/
      data     gss_org(  1)/      16.51220819D0/
      data   polvo_org(  1)/       0.26211400D0/
      data gues_org1(  1,1)/      -0.00418845D0/
      data gues_org2(  1,1)/       4.53239266D0/
      data gues_org3(  1,1)/       2.02900589D0/
!
!                    Data for Element   6           Carbon
!
      data     uss_org(  6)/     -49.97894304D0/
      data     upp_org(  6)/     -40.45461533D0/
      data   betas_org(  6)/     -11.51865824D0/
      data   betap_org(  6)/      -8.49404612D0/
      data      zs_org(  6)/       2.05935067D0/
      data      zp_org(  6)/       1.70277718D0/
      data     gss_org(  6)/      13.65002631D0/
      data     gsp_org(  6)/      11.89475833D0/
      data     gpp_org(  6)/      11.68748408D0/
      data     gp2_org(  6)/       9.87880661D0/
      data     hsp_org(  6)/       0.86511600D0/
      data   polvo_org(  6)/       0.48507100D0/
      data gues_org1(  6,1)/       0.02037941D0/
      data gues_org2(  6,1)/       5.86041774D0/
      data gues_org3(  6,1)/       1.61433014D0/
!
!                    Data for Element   7         Nitrogen
!
      data     uss_org(  7)/     -58.46660524D0/
      data     upp_org(  7)/     -44.93993591D0/
      data   betas_org(  7)/     -22.38999733D0/
      data   betap_org(  7)/     -16.84861487D0/
      data      zs_org(  7)/       2.26444238D0/
      data      zp_org(  7)/       2.15161849D0/
      data     gss_org(  7)/       8.49212330D0/
      data     gsp_org(  7)/       9.35824093D0/
      data     gpp_org(  7)/      11.15955253D0/
      data     gp2_org(  7)/       9.83500604D0/
      data     hsp_org(  7)/       4.37603155D0/
      data   polvo_org(  7)/       0.20474300D0/
      data gues_org1(  7,1)/       0.02761821D0/
      data gues_org2(  7,1)/       2.15056844D0/
      data gues_org3(  7,1)/       1.57143677D0/
!
!                    Data for Element   8           Oxygen
!
      data     uss_org(  8)/     -92.29368783D0/
      data     upp_org(  8)/     -68.63750130D0/
      data   betas_org(  8)/     -62.84670430D0/
      data   betap_org(  8)/     -17.60167625D0/
      data      zs_org(  8)/       4.84901421D0/
      data      zp_org(  8)/       2.22296191D0/
      data     gss_org(  8)/      15.79259941D0/
      data     gsp_org(  8)/      14.90978302D0/
      data     gpp_org(  8)/      13.21455827D0/
      data     gp2_org(  8)/      11.40301293D0/
      data     hsp_org(  8)/       8.87691264D0/
      data   polvo_org(  8)/       0.15430100D0/
      data gues_org1(  8,1)/      -0.00991525D0/
      data gues_org2(  8,1)/       2.33011528D0/
      data gues_org3(  8,1)/       2.21804989D0/
!
!                    Data for Element   9         Fluorine
!
      data     uss_org(  9)/    -137.49084253D0/
      data     upp_org(  9)/     -98.60644823D0/
      data   betas_org(  9)/     -76.00423918D0/
      data   betap_org(  9)/     -27.48652426D0/
      data      zs_org(  9)/       6.88765955D0/
      data      zp_org(  9)/       2.72983077D0/
      data     gss_org(  9)/       8.83106915D0/
      data     gsp_org(  9)/      20.76644115D0/
      data     gpp_org(  9)/      11.64269222D0/
      data     gp2_org(  9)/      12.12267483D0/
      data     hsp_org(  9)/       6.70622637D0/
      data   polvo_org(  9)/       0.19961100D0/
      data gues_org1(  9,1)/      -0.00784070D0/
      data gues_org2(  9,1)/       4.54137015D0/
      data gues_org3(  9,1)/       1.56706840D0/
!
!                    Data for Element  11           Sodium
!
      data     uss_org( 11)/      -4.95173405D0/
      data     upp_org( 11)/      -2.85030259D0/
      data   betas_org( 11)/      -0.10258029D0/
      data   betap_org( 11)/      -5.30239870D0/
      data      zs_org( 11)/       0.89075745D0/
      data      zp_org( 11)/       0.84312659D0/
      data     gss_org( 11)/       5.50902119D0/
      data     gsp_org( 11)/       6.41933250D0/
      data     gpp_org( 11)/      12.04894719D0/
      data     gp2_org( 11)/      14.01273067D0/
      data     hsp_org( 11)/       1.90534477D0/
      data gues_org1( 11,1)/      -0.20214451D0/
      data gues_org2( 11,1)/       3.59261867D0/
      data gues_org3( 11,1)/       2.03002270D0/
!
!                    Data for Element  12        Magnesium
!
      data     uss_org( 12)/     -11.96334304D0/
      data     upp_org( 12)/     -11.69107988D0/
      data   betas_org( 12)/      -5.01124807D0/
      data   betap_org( 12)/      -1.14370715D0/
      data      zs_org( 12)/       1.17743010D0/
      data      zp_org( 12)/       1.14720276D0/
      data     gss_org( 12)/       4.80798289D0/
      data     gsp_org( 12)/       7.18851244D0/
      data     gpp_org( 12)/      11.14377866D0/
      data     gp2_org( 12)/       6.66068989D0/
      data     hsp_org( 12)/       0.11765466D0/
!
!                    Data for Element  15       Phosphorus
!
      data     uss_org( 15)/     -48.75635103D0/
      data     upp_org( 15)/     -39.47147070D0/
      data     udd_org( 15)/      -7.09971722D0/
      data   betas_org( 15)/      -9.08367246D0/
      data   betap_org( 15)/     -12.08969359D0/
      data   betad_org( 15)/     -22.22987996D0/
      data      zs_org( 15)/       1.97679584D0/
      data      zp_org( 15)/       1.81471099D0/
      data      zd_org( 15)/       1.38655362D0/
      data     zsn_org( 15)/       4.21416730D0/
      data     zpn_org( 15)/       1.16554745D0/
      data     zdn_org( 15)/       7.95024299D0/
      data     gss_org( 15)/       8.39197936D0/
      data     gsp_org( 15)/       8.59076849D0/
      data     gpp_org( 15)/       8.89666101D0/
      data     gp2_org( 15)/       7.84071729D0/
      data     hsp_org( 15)/       1.76922744D0/
      data   polvo_org( 15)/       2.31461000D0/
      data gues_org1( 15,1)/      -0.00588080D0/
      data gues_org2( 15,1)/       5.47368452D0/
      data gues_org3( 15,1)/       2.51534669D0/
!
!                    Data for Element  16           Sulfur
!
      data     uss_org( 16)/     -52.84384720D0/
      data     upp_org( 16)/     -36.96356507D0/
      data     udd_org( 16)/     -46.28679594D0/
      data   betas_org( 16)/     -13.46369626D0/
      data   betap_org( 16)/     -10.31738995D0/
      data   betad_org( 16)/     -13.24600677D0/
      data      zs_org( 16)/       1.83414635D0/
      data      zp_org( 16)/       2.01054933D0/
      data      zd_org( 16)/       3.28706036D0/
      data     zsn_org( 16)/       2.10174879D0/
      data     zpn_org( 16)/       0.64664107D0/
      data     zdn_org( 16)/       1.75166002D0/
      data     gss_org( 16)/       8.84550473D0/
      data     gsp_org( 16)/       6.39148567D0/
      data     gpp_org( 16)/       7.26625310D0/
      data     gp2_org( 16)/       6.61407354D0/
      data     hsp_org( 16)/       4.97868817D0/
      data   polvo_org( 16)/       1.45331000D0/
      data gues_org1( 16,1)/      -0.34813009D0/
      data gues_org2( 16,1)/       0.56686034D0/
      data gues_org3( 16,1)/       0.70475091D0/
!
!                    Data for Element  17         Chlorine
!
      data     uss_org( 17)/     -66.26070535D0/
      data     upp_org( 17)/     -62.26023080D0/
      data     udd_org( 17)/     -43.70187019D0/
      data   betas_org( 17)/       1.37962037D0/
      data   betap_org( 17)/     -12.15978875D0/
      data   betad_org( 17)/      -5.34992632D0/
      data      zs_org( 17)/       2.25268102D0/
      data      zp_org( 17)/       2.08518895D0/
      data      zd_org( 17)/       1.21815745D0/
      data     zsn_org( 17)/       1.16718094D0/
      data     zpn_org( 17)/       2.37662942D0/
      data     zdn_org( 17)/       8.52113127D0/
      data     gss_org( 17)/      11.42538436D0/
      data     gsp_org( 17)/       9.47705295D0/
      data     gpp_org( 17)/      10.42961146D0/
      data     gp2_org( 17)/       8.59200761D0/
      data     hsp_org( 17)/       3.13144591D0/
      data   polvo_org( 17)/       1.23621000D0/
      data gues_org1( 17,1)/      -0.06645192D0/
      data gues_org2( 17,1)/       0.99319979D0/
      data gues_org3( 17,1)/       1.80535831D0/
!
!                    Data for Element  19        Potassium
!
      data     uss_org( 19)/      -4.48852441D0/
      data     upp_org( 19)/      -2.55854938D0/
      data   betas_org( 19)/       2.83518288D0/
      data   betap_org( 19)/      -5.45250990D0/
      data      zs_org( 19)/      20.79243083D0/
      data      zp_org( 19)/       1.18451908D0/
      data     gss_org( 19)/       8.89475202D0/
      data     gsp_org( 19)/       7.91937666D0/
      data     gpp_org( 19)/       5.82876802D0/
      data     gp2_org( 19)/      14.55099039D0/
      data     hsp_org( 19)/       0.20023681D0/
      data gues_org1( 19,1)/      -2.46232435D0/
      data gues_org2( 19,1)/       7.26344700D0/
      data gues_org3( 19,1)/       0.39859888D0/
!
!                    Data for Element  20          Calcium
!
      data     uss_org( 20)/     -11.30040321D0/
      data     upp_org( 20)/     -11.24709193D0/
      data   betas_org( 20)/     -13.95200982D0/
      data   betap_org( 20)/      -1.03687937D0/
      data      zs_org( 20)/       1.13633888D0/
      data      zp_org( 20)/       1.36391661D0/
      data     gss_org( 20)/       5.27215378D0/
      data     gsp_org( 20)/       8.01812014D0/
      data     gpp_org( 20)/      10.40294168D0/
      data     gp2_org( 20)/       9.61343941D0/
      data     hsp_org( 20)/       0.57750655D0/
      data gues_org1( 20,1)/      -0.05481494D0/
      data gues_org2( 20,1)/       3.06482279D0/
      data gues_org3( 20,1)/       2.55875837D0/
!
!                    Data for Element  26             Iron
!
      data     uss_org( 26)/     -63.11378305D0/
      data     upp_org( 26)/     -55.17003417D0/
      data     udd_org( 26)/     -92.79635767D0/
      data   betas_org( 26)/     -11.89674547D0/
      data   betap_org( 26)/      -6.15606271D0/
      data   betad_org( 26)/      -4.53813034D0/
      data      zs_org( 26)/       0.88051290D0/
      data      zp_org( 26)/       3.86416580D0/
      data      zd_org( 26)/       1.47143813D0/
      data     zsn_org( 26)/       0.87572725D0/
      data     zpn_org( 26)/       1.18049638D0/
      data     zdn_org( 26)/       1.97150826D0/
      data     gss_org( 26)/       4.78751275D0/
      data     gsp_org( 26)/       5.38521983D0/
      data     gpp_org( 26)/       7.03218047D0/
      data     gp2_org( 26)/       6.16439214D0/
      data     hsp_org( 26)/       1.11548030D0/
      data    poc__org( 26)/       1.47728348D0/
      data    f0sd_org( 26)/       8.08211779D0/
      data    g2sd_org( 26)/       2.28344696D0/
!
!                    Data for Element  27           Cobalt
!
      data     uss_org( 27)/     -22.47045814D0/
      data     upp_org( 27)/       9.61831967D0/
      data     udd_org( 27)/     -28.60773321D0/
      data   betas_org( 27)/      -5.66619572D0/
      data   betap_org( 27)/      -4.26404529D0/
      data   betad_org( 27)/      -2.93416835D0/
      data      zs_org( 27)/       1.56886152D0/
      data      zp_org( 27)/       1.54218792D0/
      data      zd_org( 27)/       1.49767547D0/
      data     zsn_org( 27)/       0.54524491D0/
      data     zpn_org( 27)/       1.75774350D0/
      data     zdn_org( 27)/       0.37542882D0/
      data     gss_org( 27)/       2.98079905D0/
      data     gsp_org( 27)/       3.69684549D0/
      data     gpp_org( 27)/      10.47082373D0/
      data     gp2_org( 27)/       9.17869838D0/
      data     hsp_org( 27)/       0.08217134D0/
      data    f0sd_org( 27)/       1.30529985D0/
      data    g2sd_org( 27)/       1.76269441D0/
!
!                    Data for Element  30             Zinc
!
      data     uss_org( 30)/     -18.23240645D0/
      data     upp_org( 30)/     -11.46772746D0/
      data   betas_org( 30)/     -14.24021703D0/
      data   betap_org( 30)/       1.36108212D0/
      data      zs_org( 30)/       1.42059168D0/
      data      zp_org( 30)/       1.74700301D0/
      data     gss_org( 30)/       9.12311506D0/
      data     gsp_org( 30)/       6.12338449D0/
      data     gpp_org( 30)/      15.82159114D0/
      data     gp2_org( 30)/       7.49082352D0/
      data     hsp_org( 30)/       0.00000010D0/
!
!                    Data for Element  34         Selenium
!
      data     uss_org( 34)/     -45.41143600D0/
      data     upp_org( 34)/     -31.48645829D0/
      data   betas_org( 34)/     -18.63044575D0/
      data   betap_org( 34)/      -6.48230673D0/
      data      zs_org( 34)/       2.00088084D0/
      data      zp_org( 34)/       1.83033932D0/
      data     gss_org( 34)/       5.35353800D0/
      data     gsp_org( 34)/       4.63343348D0/
      data     gpp_org( 34)/       6.92575107D0/
      data     gp2_org( 34)/       6.25703133D0/
      data     hsp_org( 34)/       5.23183164D0/
!
!                    Data for Element  35          Bromine
!
      data     uss_org( 35)/     -58.73303903D0/
      data     upp_org( 35)/     -47.96135197D0/
      data     udd_org( 35)/      11.95670453D0/
      data   betas_org( 35)/     -28.11981180D0/
      data   betap_org( 35)/     -12.08697290D0/
      data   betad_org( 35)/     -12.28177003D0/
      data      zs_org( 35)/       3.31096356D0/
      data      zp_org( 35)/       2.31454673D0/
      data      zd_org( 35)/       2.05905955D0/
      data     zsn_org( 35)/       1.27562087D0/
      data     zpn_org( 35)/       6.89411237D0/
      data     zdn_org( 35)/       1.64387242D0/
      data     gss_org( 35)/      12.24747271D0/
      data     gsp_org( 35)/       7.50538389D0/
      data     gpp_org( 35)/       9.63134303D0/
      data     gp2_org( 35)/       7.99752064D0/
      data     hsp_org( 35)/      10.02797159D0/
      data   polvo_org( 35)/       2.14242000D0/
      data gues_org1( 35,1)/       0.00058420D0/
      data gues_org2( 35,1)/       0.87506795D0/
      data gues_org3( 35,1)/       3.35644944D0/
!
!                    Data for Element  53           Iodine
!
      data     uss_org( 53)/     -60.49284738D0/
      data     upp_org( 53)/     -57.78472794D0/
      data     udd_org( 53)/     -49.99411990D0/
      data   betas_org( 53)/     -35.12059789D0/
      data   betap_org( 53)/      -8.66134604D0/
      data   betad_org( 53)/      -4.80443540D0/
      data      zs_org( 53)/       4.20089050D0/
      data      zp_org( 53)/       2.09377419D0/
      data      zd_org( 53)/       1.57418387D0/
      data     zsn_org( 53)/       7.31124291D0/
      data     zpn_org( 53)/       7.18205916D0/
      data     zdn_org( 53)/       1.85901047D0/
      data     gss_org( 53)/       6.82894450D0/
      data     gsp_org( 53)/       8.63341975D0/
      data     gpp_org( 53)/      10.12095945D0/
      data     gp2_org( 53)/       7.86186368D0/
      data     hsp_org( 53)/       1.38482851D0/
      data   polvo_org( 53)/       3.82316000D0/
      data gues_org1( 53,1)/      -0.87065695D0/
      data gues_org2( 53,1)/       2.81130957D0/
      data gues_org3( 53,1)/       0.50815367D0/
!
!                    Data for Element  85         Astatine
!
      data     alp_org( 85)/       3.00000000D0/
      data     gss_org( 85)/      10.00000000D0/
!
!                    Data for Element  87         Francium
!
      data     alp_org( 87)/       3.00000000D0/
      data     gss_org( 87)/      10.00000000D0/
!
!                    Data for Element 100       3+ Sparkle
!
      data     alp_org(100)/       1.50000000D0/
!
!                    Data for Element 101       3- Sparkle
!
      data     alp_org(101)/       1.50000000D0/
!
!                    Data for Element 102      Capped bond
!
      data   betas_org(102)/-9999999.00000000D0/
      data      zs_org(102)/       4.00000000D0/
      data     gss_org(102)/      12.84800000D0/
!
!                    Data for Element 103       ++ Sparkle
!
      data     alp_org(103)/       1.50000000D0/
!
!                    Data for Element 104        + Sparkle
!
      data     alp_org(104)/       1.50000000D0/
!
!                    Data for Element 105       -- Sparkle
!
      data     alp_org(105)/       1.50000000D0/
!
!                    Data for Element 106        - Sparkle
!
      data     alp_org(106)/       1.50000000D0/
!
!
!                     Global parameters
!
!
      data  v_par_org(1)/  9.27846500d0/  ! Used in ccrep C-C triple bonds.
      data  v_par_org(2)/  5.98375200d0/  ! Used in ccrep for exponent correction of C-C triple bonds.
      data  v_par_org(3)/  1.88549346d0/  ! O-H scalar                 correction.
      data  v_par_org(4)/  1.29663856d0/  ! O-H       exponent         correction.
      data  v_par_org(5)/  1.12465557d0/  ! O-H               offset   correction.
      data  v_par_org(7)/  0.69935551d0/  ! Used in dftd3 to set "s6"  in D3H4
      data  v_par_org(8)/ 22.90643249d0/  ! Used in dftd3 to set "alp" in D3H4
      data  v_par_org(9)/  0.63309647d0/  ! Used in dftd3 to set "rs6" in D3H4
      data  v_par_org(11)/ 3.16312327d0/  ! C-H       exponent         correction.
      data  v_par_org(12)/ 1.85191308d0/  ! C-H               offset   correction.
      data  v_par_org(13)/ 0.03076443d0/  ! C-C scalar                 correction.
      data  v_par_org(14)/ 5.89911268d0/  ! C-C       exponent         correction.
      data  v_par_org(15)/ 2.54117311d0/  ! C-C               offset   correction.
      data  v_par_org(16)/ 1.95072008d0/  ! H-H scalar                 correction.
      data  v_par_org(17)/ 9.93668160d0/  ! H-H       exponent         correction.
      data  v_par_org(18)/ 1.73448743d0/  ! H-H               offset   correction.
      data  v_par_org(19)/ 0.97353953d0/  ! C-H scalar                 correction.
      data  v_par_org(20)/ 0.11796327d0/  ! O-C scalar                 correction.
      data  v_par_org(21)/ 4.71122155d0/  ! O-C       exponent         correction.
      data  v_par_org(22)/ 2.33136344d0/  ! O-C               offset   correction.
      data  v_par_org(23)/-0.01163876d0/  ! S-O scalar                 correction.
      data  v_par_org(24)/ 7.15357690d0/  ! S-O       exponent         correction.
      data  v_par_org(25)/ 2.75878894d0/  ! S-O               offset   correction.
      data  v_par_org(26)/ 0.06759169d0/  ! O-N scalar                 correction.
      data  v_par_org(27)/ 7.83635303d0/  ! O-N       exponent         correction.
      data  v_par_org(28)/ 2.39364886d0/  ! O-N               offset   correction.
      data  v_par_org(29)/ 0.06917316d0/  ! F-H scalar                 correction.
      data  v_par_org(30)/ 5.86395123d0/  ! F-H       exponent         correction.
      data  v_par_org(31)/ 2.08820430d0/  ! F-H               offset   correction.
      data  v_par_org(32)/ 0.22070756d0/  ! O-O scalar                 correction.
      data  v_par_org(33)/ 6.65530337d0/  ! O-O       exponent         correction.
      data  v_par_org(34)/ 2.28968120d0/  ! O-O               offset   correction.
      data  v_par_org(35)/ 0.07450271d0/  ! N-C scalar                 correction.
      data  v_par_org(36)/ 7.59405828d0/  ! N-C       exponent         correction.
      data  v_par_org(37)/ 2.83228505d0/  ! N-C               offset   correction.
      data  v_par_org(38)/ 0.99964027d0/  ! N-H scalar                 correction.
      data  v_par_org(39)/ 4.85532734d0/  ! N-H       exponent         correction.
      data  v_par_org(40)/ 1.78032801d0/  ! N-H               offset   correction.
      data  v_par_org(41)/ 0.15560543d0/  ! S-N scalar                 correction.
      data  v_par_org(42)/ 4.94373969d0/  ! S-N       exponent         correction.
      data  v_par_org(43)/ 3.15741465d0/  ! S-N               offset   correction.
      data  v_par_org(44)/ 0.22855618d0/  ! S-H scalar                 correction.
      data  v_par_org(45)/ 6.08398751d0/  ! S-H       exponent         correction.
      data  v_par_org(46)/ 2.44074236d0/  ! S-H               offset   correction.
  contains
  subroutine alpb_and_xfac_pm6_ORG
    use parameters_C, only : xfac, alpb
 !
      alpb( 1, 1) =     5.50167002d0 !    Hydrogen -     Hydrogen
      xfac( 1, 1) =     4.96873740d0 !    Hydrogen -     Hydrogen
 !
      alpb( 6, 1) =     1.09926251d0 !      Carbon -     Hydrogen
      xfac( 6, 1) =     0.20244255d0 !      Carbon -     Hydrogen
      alpb( 6, 6) =     2.49174033d0 !      Carbon -       Carbon
      xfac( 6, 6) =     0.75827957d0 !      Carbon -       Carbon
 !
      alpb( 7, 1) =     1.18019929d0 !    Nitrogen -     Hydrogen
      xfac( 7, 1) =     0.19459121d0 !    Nitrogen -     Hydrogen
      alpb( 7, 6) =     2.91019769d0 !    Nitrogen -       Carbon
      xfac( 7, 6) =     1.19919594d0 !    Nitrogen -       Carbon
      alpb( 7, 7) =     3.26594309d0 !    Nitrogen -     Nitrogen
      xfac( 7, 7) =     1.56001466d0 !    Nitrogen -     Nitrogen
 !
      alpb( 8, 1) =     1.60341671d0 !      Oxygen -     Hydrogen
      xfac( 8, 1) =     0.14777972d0 !      Oxygen -     Hydrogen
      alpb( 8, 6) =     2.83076261d0 !      Oxygen -       Carbon
      xfac( 8, 6) =     0.76483641d0 !      Oxygen -       Carbon
      alpb( 8, 7) =     3.22126658d0 !      Oxygen -     Nitrogen
      xfac( 8, 7) =     1.02400611d0 !      Oxygen -     Nitrogen
      alpb( 8, 8) =     3.22832635d0 !      Oxygen -       Oxygen
      xfac( 8, 8) =     0.75857042d0 !      Oxygen -       Oxygen
 !
      alpb( 9, 1) =     3.00908696d0 !    Fluorine -     Hydrogen
      xfac( 9, 1) =     0.57468838d0 !    Fluorine -     Hydrogen
      alpb( 9, 6) =     2.99751365d0 !    Fluorine -       Carbon
      xfac( 9, 6) =     0.76790658d0 !    Fluorine -       Carbon
      alpb( 9, 7) =     3.70135239d0 !    Fluorine -     Nitrogen
      xfac( 9, 7) =     1.83973176d0 !    Fluorine -     Nitrogen
      alpb( 9, 8) =     3.04821245d0 !    Fluorine -       Oxygen
      xfac( 9, 8) =     0.62903055d0 !    Fluorine -       Oxygen
      alpb( 9, 9) =     4.01574430d0 !    Fluorine -     Fluorine
      xfac( 9, 9) =     2.09551890d0 !    Fluorine -     Fluorine
 !
      alpb(11, 1) =     1.27239598d0 !      Sodium -     Hydrogen
      xfac(11, 1) =     1.79556240d0 !      Sodium -     Hydrogen
      alpb(11, 6) =     0.72742814d0 !      Sodium -       Carbon
      xfac(11, 6) =     0.33826248d0 !      Sodium -       Carbon
      alpb(11, 7) =     2.03351539d0 !      Sodium -     Nitrogen
      xfac(11, 7) =     7.11947177d0 !      Sodium -     Nitrogen
      alpb(11, 8) =     1.07417163d0 !      Sodium -       Oxygen
      xfac(11, 8) =     0.46593554d0 !      Sodium -       Oxygen
      alpb(11, 9) =     1.50126936d0 !      Sodium -     Fluorine
      xfac(11, 9) =     0.48884039d0 !      Sodium -     Fluorine
      alpb(11,11) =     0.49603378d0 !      Sodium -       Sodium
      xfac(11,11) =     0.52981333d0 !      Sodium -       Sodium
 !
      alpb(12, 1) =     1.74026385d0 !   Magnesium -     Hydrogen
      xfac(12, 1) =     1.64641220d0 !   Magnesium -     Hydrogen
      alpb(12, 6) =     1.91332801d0 !   Magnesium -       Carbon
      xfac(12, 6) =     1.56070554d0 !   Magnesium -       Carbon
      alpb(12, 7) =     1.24194361d0 !   Magnesium -     Nitrogen
      xfac(12, 7) =     0.42400581d0 !   Magnesium -     Nitrogen
      alpb(12, 8) =     1.64306809d0 !   Magnesium -       Oxygen
      xfac(12, 8) =     0.64270914d0 !   Magnesium -       Oxygen
      alpb(12, 9) =     4.10107545d0 !   Magnesium -     Fluorine
      xfac(12, 9) =    16.89931519d0 !   Magnesium -     Fluorine
      alpb(12,11) =     1.50677300d0 !   Magnesium -       Sodium
      xfac(12,11) =     8.67561900d0 !   Magnesium -       Sodium
      alpb(12,12) =     0.25961712d0 !   Magnesium -    Magnesium
      xfac(12,12) =     0.17210861d0 !   Magnesium -    Magnesium
 !
      alpb(15, 1) =     1.77701283d0 !  Phosphorus -     Hydrogen
      xfac(15, 1) =     0.92231547d0 !  Phosphorus -     Hydrogen
      alpb(15, 6) =     1.94273921d0 !  Phosphorus -       Carbon
      xfac(15, 6) =     0.86185262d0 !  Phosphorus -       Carbon
      alpb(15, 7) =     2.08046079d0 !  Phosphorus -     Nitrogen
      xfac(15, 7) =     0.85814293d0 !  Phosphorus -     Nitrogen
      alpb(15, 8) =     2.22358855d0 !  Phosphorus -       Oxygen
      xfac(15, 8) =     0.77229960d0 !  Phosphorus -       Oxygen
      alpb(15, 9) =     2.40049800d0 !  Phosphorus -     Fluorine
      xfac(15, 9) =     0.68791787d0 !  Phosphorus -     Fluorine
      alpb(15,11) =     1.50032000d0 !  Phosphorus -       Sodium
      xfac(15,11) =     2.83709500d0 !  Phosphorus -       Sodium
      alpb(15,12) =     1.33601450d0 !  Phosphorus -    Magnesium
      xfac(15,12) =     1.19491696d0 !  Phosphorus -    Magnesium
      alpb(15,13) =     1.98072700d0 !  Phosphorus -     Aluminum
      xfac(15,13) =     5.05081600d0 !  Phosphorus -     Aluminum
      alpb(15,14) =     3.31346600d0 !  Phosphorus -      Silicon
      xfac(15,14) =    13.23912100d0 !  Phosphorus -      Silicon
      alpb(15,15) =     1.31423390d0 !  Phosphorus -   Phosphorus
      xfac(15,15) =     0.53298471d0 !  Phosphorus -   Phosphorus
 !
      alpb(16, 1) =     1.99051064d0 !      Sulfur -     Hydrogen
      xfac(16, 1) =     0.87147940d0 !      Sulfur -     Hydrogen
      alpb(16, 6) =     2.04551959d0 !      Sulfur -       Carbon
      xfac(16, 6) =     0.81385984d0 !      Sulfur -       Carbon
      alpb(16, 7) =     2.38647451d0 !      Sulfur -     Nitrogen
      xfac(16, 7) =     1.24545127d0 !      Sulfur -     Nitrogen
      alpb(16, 8) =     2.00418133d0 !      Sulfur -       Oxygen
      xfac(16, 8) =     0.57116134d0 !      Sulfur -       Oxygen
      alpb(16, 9) =     2.02072554d0 !      Sulfur -     Fluorine
      xfac(16, 9) =     0.53705603d0 !      Sulfur -     Fluorine
      alpb(16,11) =     1.24987406d0 !      Sulfur -       Sodium
      xfac(16,11) =     0.96273055d0 !      Sulfur -       Sodium
      alpb(16,12) =     1.65688125d0 !      Sulfur -    Magnesium
      xfac(16,12) =     1.64556653d0 !      Sulfur -    Magnesium
      alpb(16,13) =     1.97670500d0 !      Sulfur -     Aluminum
      xfac(16,13) =     2.34738400d0 !      Sulfur -     Aluminum
      alpb(16,14) =     1.88591600d0 !      Sulfur -      Silicon
      xfac(16,14) =     0.87665800d0 !      Sulfur -      Silicon
      alpb(16,15) =     1.72980793d0 !      Sulfur -   Phosphorus
      xfac(16,15) =     0.92899705d0 !      Sulfur -   Phosphorus
      alpb(16,16) =     1.82481904d0 !      Sulfur -       Sulfur
      xfac(16,16) =     0.95994800d0 !      Sulfur -       Sulfur
 !
      alpb(17, 1) =     1.77202317d0 !    Chlorine -     Hydrogen
      xfac(17, 1) =     0.35128616d0 !    Chlorine -     Hydrogen
      alpb(17, 6) =     1.77000810d0 !    Chlorine -       Carbon
      xfac(17, 6) =     0.31362264d0 !    Chlorine -       Carbon
      alpb(17, 7) =     1.88009782d0 !    Chlorine -     Nitrogen
      xfac(17, 7) =     0.37268446d0 !    Chlorine -     Nitrogen
      alpb(17, 8) =     1.72179335d0 !    Chlorine -       Oxygen
      xfac(17, 8) =     0.23922164d0 !    Chlorine -       Oxygen
      alpb(17, 9) =     1.75531626d0 !    Chlorine -     Fluorine
      xfac(17, 9) =     0.20571654d0 !    Chlorine -     Fluorine
      alpb(17,11) =     1.81642900d0 !    Chlorine -       Sodium
      xfac(17,11) =     1.35789400d0 !    Chlorine -       Sodium
      alpb(17,12) =     1.92937720d0 !    Chlorine -    Magnesium
      xfac(17,12) =     1.41507295d0 !    Chlorine -    Magnesium
      alpb(17,13) =     2.12593900d0 !    Chlorine -     Aluminum
      xfac(17,13) =     2.15345100d0 !    Chlorine -     Aluminum
      alpb(17,14) =     1.68497800d0 !    Chlorine -      Silicon
      xfac(17,14) =     0.51300000d0 !    Chlorine -      Silicon
      alpb(17,15) =     1.52465837d0 !    Chlorine -   Phosphorus
      xfac(17,15) =     0.40850293d0 !    Chlorine -   Phosphorus
      alpb(17,16) =     1.63061304d0 !    Chlorine -       Sulfur
      xfac(17,16) =     0.48323910d0 !    Chlorine -       Sulfur
      alpb(17,17) =     1.68012133d0 !    Chlorine -     Chlorine
      xfac(17,17) =     0.38845132d0 !    Chlorine -     Chlorine
 !
      alpb(19, 1) =     0.92695277d0 !   Potassium -     Hydrogen
      xfac(19, 1) =     0.99027410d0 !   Potassium -     Hydrogen
      alpb(19, 6) =     1.51170545d0 !   Potassium -       Carbon
      xfac(19, 6) =     1.50180357d0 !   Potassium -       Carbon
      alpb(19, 7) =     1.52806612d0 !   Potassium -     Nitrogen
      xfac(19, 7) =     1.85579457d0 !   Potassium -     Nitrogen
      alpb(19, 8) =     1.52565876d0 !   Potassium -       Oxygen
      xfac(19, 8) =     1.09688509d0 !   Potassium -       Oxygen
      alpb(19, 9) =     2.89600597d0 !   Potassium -     Fluorine
      xfac(19, 9) =     6.96328886d0 !   Potassium -     Fluorine
      alpb(19,11) =     0.88430700d0 !   Potassium -       Sodium
      xfac(19,11) =     5.56302700d0 !   Potassium -       Sodium
      alpb(19,12) =     0.88481000d0 !   Potassium -    Magnesium
      xfac(19,12) =     3.29050200d0 !   Potassium -    Magnesium
      alpb(19,13) =     1.97607600d0 !   Potassium -     Aluminum
      xfac(19,13) =    29.94470800d0 !   Potassium -     Aluminum
      alpb(19,14) =     1.67593000d0 !   Potassium -      Silicon
      xfac(19,14) =     8.27920000d0 !   Potassium -      Silicon
      alpb(19,15) =     1.44373800d0 !   Potassium -   Phosphorus
      xfac(19,15) =     4.47538400d0 !   Potassium -   Phosphorus
      alpb(19,16) =     2.34162192d0 !   Potassium -       Sulfur
      xfac(19,16) =    28.74098661d0 !   Potassium -       Sulfur
      alpb(19,17) =     1.56138261d0 !   Potassium -     Chlorine
      xfac(19,17) =     1.76932814d0 !   Potassium -     Chlorine
      alpb(19,18) =     2.30280300d0 !   Potassium -        Argon
      xfac(19,18) =     9.71050800d0 !   Potassium -        Argon
      alpb(19,19) =     1.14821680d0 !   Potassium -    Potassium
      xfac(19,19) =     4.60699839d0 !   Potassium -    Potassium
 !
      alpb(20, 1) =     1.69147617d0 !     Calcium -     Hydrogen
      xfac(20, 1) =     4.09412903d0 !     Calcium -     Hydrogen
      alpb(20, 5) =     1.70001000d0 !     Calcium -        Boron
      xfac(20, 5) =     1.70001000d0 !     Calcium -        Boron
      alpb(20, 6) =     1.05066825d0 !     Calcium -       Carbon
      xfac(20, 6) =     0.32236564d0 !     Calcium -       Carbon
      alpb(20, 7) =     1.11272184d0 !     Calcium -     Nitrogen
      xfac(20, 7) =     0.36106406d0 !     Calcium -     Nitrogen
      alpb(20, 8) =     1.00032461d0 !     Calcium -       Oxygen
      xfac(20, 8) =     0.16930809d0 !     Calcium -       Oxygen
      alpb(20, 9) =     3.85224570d0 !     Calcium -     Fluorine
      xfac(20, 9) =    10.98269271d0 !     Calcium -     Fluorine
      alpb(20,11) =     3.10710400d0 !     Calcium -       Sodium
      xfac(20,11) =     9.65750900d0 !     Calcium -       Sodium
      alpb(20,12) =     2.29980000d0 !     Calcium -    Magnesium
      xfac(20,12) =     8.59980000d0 !     Calcium -    Magnesium
      alpb(20,14) =     1.21878800d0 !     Calcium -      Silicon
      xfac(20,14) =     0.33623300d0 !     Calcium -      Silicon
      alpb(20,15) =     1.02414200d0 !     Calcium -   Phosphorus
      xfac(20,15) =     0.41084000d0 !     Calcium -   Phosphorus
      alpb(20,16) =     0.54371148d0 !     Calcium -       Sulfur
      xfac(20,16) =     0.20200309d0 !     Calcium -       Sulfur
      alpb(20,17) =     0.81998389d0 !     Calcium -     Chlorine
      xfac(20,17) =     0.16380876d0 !     Calcium -     Chlorine
      alpb(20,18) =     1.03488100d0 !     Calcium -        Argon
      xfac(20,18) =     0.29107200d0 !     Calcium -        Argon
      alpb(20,19) =     1.11920000d0 !     Calcium -    Potassium
      xfac(20,19) =     1.24032000d0 !     Calcium -    Potassium
      alpb(20,20) =     1.14828793d0 !     Calcium -      Calcium
      xfac(20,20) =    29.47195357d0 !     Calcium -      Calcium
 !
      alpb(26, 1) =     0.79193667d0 !        Iron -     Hydrogen
      xfac(26, 1) =     0.05838526d0 !        Iron -     Hydrogen
      alpb(26, 6) =     3.99017141d0 !        Iron -       Carbon
      xfac(26, 6) =     5.33514976d0 !        Iron -       Carbon
      alpb(26, 7) =     4.02796928d0 !        Iron -     Nitrogen
      xfac(26, 7) =     5.63401184d0 !        Iron -     Nitrogen
      alpb(26, 8) =     3.23501120d0 !        Iron -       Oxygen
      xfac(26, 8) =     1.61877513d0 !        Iron -       Oxygen
      alpb(26, 9) =     3.27029865d0 !        Iron -     Fluorine
      xfac(26, 9) =     1.57958147d0 !        Iron -     Fluorine
      alpb(26,15) =     1.42597278d0 !        Iron -   Phosphorus
      xfac(26,15) =     0.59759172d0 !        Iron -   Phosphorus
      alpb(26,16) =     2.39883900d0 !        Iron -       Sulfur
      xfac(26,16) =     1.26086296d0 !        Iron -       Sulfur
      alpb(26,17) =     1.96664930d0 !        Iron -     Chlorine
      xfac(26,17) =     0.31495551d0 !        Iron -     Chlorine
      alpb(26,19) =     2.00000000d0 !        Iron -    Potassium
      xfac(26,19) =     6.00000000d0 !        Iron -    Potassium
      alpb(26,26) =     3.56505839d0 !        Iron -         Iron
      xfac(26,26) =    24.92536659d0 !        Iron -         Iron
 !
      alpb(27, 1) =     2.34165801d0 !      Cobalt -     Hydrogen
      xfac(27, 1) =     1.94845346d0 !      Cobalt -     Hydrogen
      alpb(27, 5) =     3.20000000d0 !      Cobalt -        Boron
      xfac(27, 5) =     1.00000000d0 !      Cobalt -        Boron
      alpb(27, 6) =     3.55962148d0 !      Cobalt -       Carbon
      xfac(27, 6) =     5.07784400d0 !      Cobalt -       Carbon
      alpb(27, 7) =     3.10684864d0 !      Cobalt -     Nitrogen
      xfac(27, 7) =     2.57511613d0 !      Cobalt -     Nitrogen
      alpb(27, 8) =     3.12463554d0 !      Cobalt -       Oxygen
      xfac(27, 8) =     2.40385387d0 !      Cobalt -       Oxygen
      alpb(27, 9) =     3.44800954d0 !      Cobalt -     Fluorine
      xfac(27, 9) =     3.53698290d0 !      Cobalt -     Fluorine
      alpb(27,15) =     1.15939687d0 !      Cobalt -   Phosphorus
      xfac(27,15) =     0.05506725d0 !      Cobalt -   Phosphorus
      alpb(27,16) =     2.11118270d0 !      Cobalt -       Sulfur
      xfac(27,16) =     1.15338065d0 !      Cobalt -       Sulfur
      alpb(27,17) =     1.95752618d0 !      Cobalt -     Chlorine
      xfac(27,17) =     0.34957665d0 !      Cobalt -     Chlorine
      alpb(27,27) =     3.43844100d0 !      Cobalt -       Cobalt
      xfac(27,27) =     3.99349531d0 !      Cobalt -       Cobalt
 !
      alpb(30, 1) =     2.09471205d0 !        Zinc -     Hydrogen
      xfac(30, 1) =     3.90267492d0 !        Zinc -     Hydrogen
      alpb(30, 6) =     1.61726178d0 !        Zinc -       Carbon
      xfac(30, 6) =     0.69422546d0 !        Zinc -       Carbon
      alpb(30, 7) =     1.53345218d0 !        Zinc -     Nitrogen
      xfac(30, 7) =     0.58920321d0 !        Zinc -     Nitrogen
      alpb(30, 8) =     1.93125665d0 !        Zinc -       Oxygen
      xfac(30, 8) =     0.87718124d0 !        Zinc -       Oxygen
      alpb(30, 9) =     2.41002100d0 !        Zinc -     Fluorine
      xfac(30, 9) =     1.22554500d0 !        Zinc -     Fluorine
      alpb(30,15) =     1.22048000d0 !        Zinc -   Phosphorus
      xfac(30,15) =     0.58153000d0 !        Zinc -   Phosphorus
      alpb(30,16) =     1.37971732d0 !        Zinc -       Sulfur
      xfac(30,16) =     0.77139403d0 !        Zinc -       Sulfur
      alpb(30,17) =     1.31543174d0 !        Zinc -     Chlorine
      xfac(30,17) =     0.35258456d0 !        Zinc -     Chlorine
      alpb(30,20) =     1.11918000d0 !        Zinc -      Calcium
      xfac(30,20) =     1.24029000d0 !        Zinc -      Calcium
      alpb(30,30) =     0.63690130d0 !        Zinc -         Zinc
      xfac(30,30) =     0.20744744d0 !        Zinc -         Zinc
 !
      alpb(34, 1) =     1.48506586d0 !    Selenium -     Hydrogen
      xfac(34, 1) =     0.47959118d0 !    Selenium -     Hydrogen
      alpb(34, 6) =     2.03352801d0 !    Selenium -       Carbon
      xfac(34, 6) =     0.53659777d0 !    Selenium -       Carbon
      alpb(34, 7) =     1.73268900d0 !    Selenium -     Nitrogen
      xfac(34, 7) =     0.22936373d0 !    Selenium -     Nitrogen
      alpb(34, 8) =     1.97750067d0 !    Selenium -       Oxygen
      xfac(34, 8) =     0.30794027d0 !    Selenium -       Oxygen
      alpb(34, 9) =     1.91312772d0 !    Selenium -     Fluorine
      xfac(34, 9) =     0.20007031d0 !    Selenium -     Fluorine
      alpb(34,15) =     1.14004774d0 !    Selenium -   Phosphorus
      xfac(34,15) =     0.26921054d0 !    Selenium -   Phosphorus
      alpb(34,16) =     1.39091004d0 !    Selenium -       Sulfur
      xfac(34,16) =     0.46848253d0 !    Selenium -       Sulfur
      alpb(34,17) =     1.36119988d0 !    Selenium -     Chlorine
      xfac(34,17) =     0.23519901d0 !    Selenium -     Chlorine
      alpb(34,27) =     2.52345000d0 !    Selenium -       Cobalt
      xfac(34,27) =     2.20241000d0 !    Selenium -       Cobalt
      alpb(34,30) =     1.08896205d0 !    Selenium -         Zinc
      xfac(34,30) =     0.59959280d0 !    Selenium -         Zinc
      alpb(34,34) =     1.06068552d0 !    Selenium -     Selenium
      xfac(34,34) =     0.15281377d0 !    Selenium -     Selenium
 !
      alpb(35, 1) =     2.19042544d0 !     Bromine -     Hydrogen
      xfac(35, 1) =     0.95776729d0 !     Bromine -     Hydrogen
      alpb(35, 6) =     2.28719180d0 !     Bromine -       Carbon
      xfac(35, 6) =     1.00009158d0 !     Bromine -       Carbon
      alpb(35, 7) =     3.97658112d0 !     Bromine -     Nitrogen
      xfac(35, 7) =    27.73297436d0 !     Bromine -     Nitrogen
      alpb(35, 8) =     2.87622656d0 !     Bromine -       Oxygen
      xfac(35, 8) =     1.81002542d0 !     Bromine -       Oxygen
      alpb(35, 9) =     2.75170253d0 !     Bromine -     Fluorine
      xfac(35, 9) =     0.99530599d0 !     Bromine -     Fluorine
      alpb(35,11) =     1.62221800d0 !     Bromine -       Sodium
      xfac(35,11) =     1.75293700d0 !     Bromine -       Sodium
      alpb(35,12) =     1.76836414d0 !     Bromine -    Magnesium
      xfac(35,12) =     2.24092055d0 !     Bromine -    Magnesium
      alpb(35,15) =     1.73435576d0 !     Bromine -   Phosphorus
      xfac(35,15) =     0.86422011d0 !     Bromine -   Phosphorus
      alpb(35,16) =     2.24818777d0 !     Bromine -       Sulfur
      xfac(35,16) =     1.83250480d0 !     Bromine -       Sulfur
      alpb(35,17) =     1.47106975d0 !     Bromine -     Chlorine
      xfac(35,17) =     0.20035006d0 !     Bromine -     Chlorine
      alpb(35,19) =     1.58047134d0 !     Bromine -    Potassium
      xfac(35,19) =     4.76810610d0 !     Bromine -    Potassium
      alpb(35,20) =     1.04553736d0 !     Bromine -      Calcium
      xfac(35,20) =     0.53031339d0 !     Bromine -      Calcium
      alpb(35,26) =     3.16463530d0 !     Bromine -         Iron
      xfac(35,26) =     4.41099184d0 !     Bromine -         Iron
      alpb(35,27) =     0.42328296d0 !     Bromine -       Cobalt
      xfac(35,27) =     0.03509636d0 !     Bromine -       Cobalt
      alpb(35,30) =     1.47179767d0 !     Bromine -         Zinc
      xfac(35,30) =     0.85333446d0 !     Bromine -         Zinc
      alpb(35,34) =     2.60066785d0 !     Bromine -     Selenium
      xfac(35,34) =     5.81873823d0 !     Bromine -     Selenium
      alpb(35,35) =     2.24012626d0 !     Bromine -      Bromine
      xfac(35,35) =     1.64245682d0 !     Bromine -      Bromine
 !
      alpb(53, 1) =     2.03084264d0 !      Iodine -     Hydrogen
      xfac(53, 1) =     1.56972456d0 !      Iodine -     Hydrogen
      alpb(53, 5) =     2.66760500d0 !      Iodine -        Boron
      xfac(53, 5) =     3.16138500d0 !      Iodine -        Boron
      alpb(53, 6) =     1.83816299d0 !      Iodine -       Carbon
      xfac(53, 6) =     0.85052273d0 !      Iodine -       Carbon
      alpb(53, 7) =     1.88024651d0 !      Iodine -     Nitrogen
      xfac(53, 7) =     0.71355484d0 !      Iodine -     Nitrogen
      alpb(53, 8) =     1.89231884d0 !      Iodine -       Oxygen
      xfac(53, 8) =     0.49996437d0 !      Iodine -       Oxygen
      alpb(53, 9) =     2.08961050d0 !      Iodine -     Fluorine
      xfac(53, 9) =     0.49882184d0 !      Iodine -     Fluorine
      alpb(53,11) =     1.40309000d0 !      Iodine -       Sodium
      xfac(53,11) =     1.98611200d0 !      Iodine -       Sodium
      alpb(53,12) =     2.04513700d0 !      Iodine -    Magnesium
      xfac(53,12) =     3.27691400d0 !      Iodine -    Magnesium
      alpb(53,15) =     1.66521677d0 !      Iodine -   Phosphorus
      xfac(53,15) =     2.28682057d0 !      Iodine -   Phosphorus
      alpb(53,16) =     1.73982165d0 !      Iodine -       Sulfur
      xfac(53,16) =     0.82461290d0 !      Iodine -       Sulfur
      alpb(53,17) =     1.76148881d0 !      Iodine -     Chlorine
      xfac(53,17) =     0.49999855d0 !      Iodine -     Chlorine
      alpb(53,19) =     1.09923911d0 !      Iodine -    Potassium
      xfac(53,19) =     2.79051582d0 !      Iodine -    Potassium
      alpb(53,20) =     0.82340461d0 !      Iodine -      Calcium
      xfac(53,20) =     0.77638091d0 !      Iodine -      Calcium
      alpb(53,26) =     1.48387437d0 !      Iodine -         Iron
      xfac(53,26) =     0.31060459d0 !      Iodine -         Iron
      alpb(53,27) =     1.06507565d0 !      Iodine -       Cobalt
      xfac(53,27) =     0.15484143d0 !      Iodine -       Cobalt
      alpb(53,30) =     1.21426495d0 !      Iodine -         Zinc
      xfac(53,30) =     1.28262860d0 !      Iodine -         Zinc
      alpb(53,35) =     1.43409914d0 !      Iodine -      Bromine
      xfac(53,35) =     0.49975667d0 !      Iodine -      Bromine
      alpb(53,53) =     1.24208422d0 !      Iodine -       Iodine
      xfac(53,53) =     0.49512749d0 !      Iodine -       Iodine
    end subroutine alpb_and_xfac_pm6_ORG
  end module Parameters_for_PM6_ORG_C
