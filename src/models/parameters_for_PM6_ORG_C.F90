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
      data     uss_org(  1)/     -11.51804056D0/
      data   betas_org(  1)/      -9.02609740D0/
      data      zs_org(  1)/       1.25406043D0/
      data     gss_org(  1)/      15.86493302D0/
      data   polvo_org(  1)/       0.26211400D0/
      data gues_org1(  1,1)/      -0.00065824D0/
      data gues_org2(  1,1)/       3.59078585D0/
      data gues_org3(  1,1)/       2.35426666D0/
!
!                    Data for Element   2           Helium
!
      data     uss_org(  2)/     -31.77096900D0/
      data     upp_org(  2)/      -5.85638200D0/
      data   betas_org(  2)/     -58.90377400D0/
      data   betap_org(  2)/     -37.03997400D0/
      data      zs_org(  2)/       3.31320400D0/
      data      zp_org(  2)/       3.65713300D0/
      data     gss_org(  2)/       9.44529900D0/
      data     gsp_org(  2)/      11.20141900D0/
      data     gpp_org(  2)/       9.21454800D0/
      data     gp2_org(  2)/      13.04611500D0/
      data     hsp_org(  2)/       0.29995400D0/
!
!                    Data for Element   3          Lithium
!
      data     uss_org(  3)/      -4.70991200D0/
      data     upp_org(  3)/      -2.72258100D0/
      data   betas_org(  3)/      -2.28394600D0/
      data   betap_org(  3)/      -7.53557300D0/
      data      zs_org(  3)/       0.98104100D0/
      data      zp_org(  3)/       2.95344500D0/
      data     gss_org(  3)/      11.03590700D0/
      data     gsp_org(  3)/      19.99864700D0/
      data     gpp_org(  3)/      11.54365000D0/
      data     gp2_org(  3)/       9.05903600D0/
      data     hsp_org(  3)/       1.64188600D0/
!
!                    Data for Element   4        Beryllium
!
      data     uss_org(  4)/     -16.36031500D0/
      data     upp_org(  4)/     -16.33921600D0/
      data   betas_org(  4)/      -3.19954900D0/
      data   betap_org(  4)/      -4.45192000D0/
      data      zs_org(  4)/       1.21253900D0/
      data      zp_org(  4)/       1.27648700D0/
      data     gss_org(  4)/       7.55280400D0/
      data     gsp_org(  4)/      10.20314600D0/
      data     gpp_org(  4)/      12.86215300D0/
      data     gp2_org(  4)/      13.60285800D0/
      data     hsp_org(  4)/       1.50145200D0/
      data gues_org1(  4,1)/       0.16418000D0/
      data gues_org2(  4,1)/       1.70482800D0/
      data gues_org3(  4,1)/       1.78559100D0/
!
!                    Data for Element   5            Boron
!
      data     uss_org(  5)/     -25.96767900D0/
      data     upp_org(  5)/     -19.11586400D0/
      data   betas_org(  5)/      -4.95970600D0/
      data   betap_org(  5)/      -4.65675300D0/
      data      zs_org(  5)/       1.63417400D0/
      data      zp_org(  5)/       1.47919500D0/
      data     gss_org(  5)/       8.17934100D0/
      data     gsp_org(  5)/       7.29402100D0/
      data     gpp_org(  5)/       7.82939500D0/
      data     gp2_org(  5)/       6.40107200D0/
      data     hsp_org(  5)/       1.25284500D0/
!
!                    Data for Element   6           Carbon
!
      data     uss_org(  6)/     -50.21748281D0/
      data     upp_org(  6)/     -40.24936947D0/
      data   betas_org(  6)/     -13.19290435D0/
      data   betap_org(  6)/      -8.12535240D0/
      data      zs_org(  6)/       2.06257083D0/
      data      zp_org(  6)/       1.68590132D0/
      data     gss_org(  6)/      12.54662491D0/
      data     gsp_org(  6)/      12.05342723D0/
      data     gpp_org(  6)/      11.51643967D0/
      data     gp2_org(  6)/       9.58468368D0/
      data     hsp_org(  6)/       0.77526172D0/
      data   polvo_org(  6)/       0.48507100D0/
      data gues_org1(  6,1)/       0.02622155D0/
      data gues_org2(  6,1)/       5.74556494D0/
      data gues_org3(  6,1)/       1.55433422D0/
!
!                    Data for Element   7         Nitrogen
!
      data     uss_org(  7)/     -59.39873932D0/
      data     upp_org(  7)/     -46.29419110D0/
      data   betas_org(  7)/     -21.25148749D0/
      data   betap_org(  7)/     -16.80968497D0/
      data      zs_org(  7)/       2.24595554D0/
      data      zp_org(  7)/       2.10239125D0/
      data     gss_org(  7)/      10.36135750D0/
      data     gsp_org(  7)/      10.45048800D0/
      data     gpp_org(  7)/      12.06395930D0/
      data     gp2_org(  7)/      10.44042562D0/
      data     hsp_org(  7)/       7.54400703D0/
      data   polvo_org(  7)/       0.20474300D0/
      data gues_org1(  7,1)/       0.02396931D0/
      data gues_org2(  7,1)/       2.64054835D0/
      data gues_org3(  7,1)/       1.63598591D0/
!
!                    Data for Element   8           Oxygen
!
      data     uss_org(  8)/     -92.44667207D0/
      data     upp_org(  8)/     -68.88994307D0/
      data   betas_org(  8)/     -63.94814293D0/
      data   betap_org(  8)/     -19.65680718D0/
      data      zs_org(  8)/       5.46670031D0/
      data      zp_org(  8)/       2.36483491D0/
      data     gss_org(  8)/      15.44708890D0/
      data     gsp_org(  8)/      15.62015275D0/
      data     gpp_org(  8)/      12.61015274D0/
      data     gp2_org(  8)/      10.67382310D0/
      data     hsp_org(  8)/       7.13812989D0/
      data   polvo_org(  8)/       0.15430100D0/
      data gues_org1(  8,1)/      -0.00407258D0/
      data gues_org2(  8,1)/       0.92192926D0/
      data gues_org3(  8,1)/       1.87329535D0/
!
!                    Data for Element   9         Fluorine
!
      data     uss_org(  9)/    -137.46462929D0/
      data     upp_org(  9)/     -98.88187919D0/
      data   betas_org(  9)/     -76.27545211D0/
      data   betap_org(  9)/     -27.95305669D0/
      data      zs_org(  9)/       7.11585740D0/
      data      zp_org(  9)/       2.77565793D0/
      data     gss_org(  9)/       9.11675979D0/
      data     gsp_org(  9)/      20.42916477D0/
      data     gpp_org(  9)/      10.61179864D0/
      data     gp2_org(  9)/      12.12433246D0/
      data     hsp_org(  9)/       6.24805955D0/
      data   polvo_org(  9)/       0.19961100D0/
      data gues_org1(  9,1)/      -0.00904697D0/
      data gues_org2(  9,1)/       3.93358548D0/
      data gues_org3(  9,1)/       1.87267339D0/
!
!                    Data for Element  10             Neon
!
      data     uss_org( 10)/      -2.97872900D0/
      data     upp_org( 10)/     -85.44111800D0/
      data   betas_org( 10)/     -69.79347500D0/
      data   betap_org( 10)/     -33.26196200D0/
      data      zs_org( 10)/       6.00014800D0/
      data      zp_org( 10)/       3.83452800D0/
      data     gss_org( 10)/      19.99957400D0/
      data     gsp_org( 10)/      16.89695100D0/
      data     gpp_org( 10)/       8.96356000D0/
      data     gp2_org( 10)/      16.02779900D0/
      data     hsp_org( 10)/       1.77928000D0/
!
!                    Data for Element  11           Sodium
!
      data     uss_org( 11)/      -4.83260447D0/
      data     upp_org( 11)/      -2.73144016D0/
      data   betas_org( 11)/      -0.19109424D0/
      data   betap_org( 11)/      -5.89251999D0/
      data      zs_org( 11)/       0.76390472D0/
      data      zp_org( 11)/       0.67900479D0/
      data     gss_org( 11)/       5.05187143D0/
      data     gsp_org( 11)/       6.58150063D0/
      data     gpp_org( 11)/      12.12326108D0/
      data     gp2_org( 11)/      13.88371923D0/
      data     hsp_org( 11)/       1.67398150D0/
      data gues_org1( 11,1)/      -0.17685775D0/
      data gues_org2( 11,1)/       3.48357540D0/
      data gues_org3( 11,1)/       1.88145523D0/
!
!                    Data for Element  12        Magnesium
!
      data     uss_org( 12)/     -12.58764527D0/
      data     upp_org( 12)/     -13.64783717D0/
      data   betas_org( 12)/      -3.06349561D0/
      data   betap_org( 12)/      -1.61715918D0/
      data      zs_org( 12)/       1.36341647D0/
      data      zp_org( 12)/       1.49515190D0/
      data     gss_org( 12)/       4.79939024D0/
      data     gsp_org( 12)/       6.80385889D0/
      data     gpp_org( 12)/       9.71740128D0/
      data     gp2_org( 12)/       8.37871633D0/
      data     hsp_org( 12)/       0.63119158D0/
!
!                    Data for Element  13         Aluminum
!
      data     uss_org( 13)/     -24.54677800D0/
      data     upp_org( 13)/     -20.10443400D0/
      data     udd_org( 13)/       8.00439400D0/
      data   betas_org( 13)/     -18.37522900D0/
      data   betap_org( 13)/      -9.38270000D0/
      data   betad_org( 13)/     -20.84047400D0/
      data      zs_org( 13)/       2.36426400D0/
      data      zp_org( 13)/       1.74910200D0/
      data      zd_org( 13)/       1.26938400D0/
      data     zsn_org( 13)/       4.74234100D0/
      data     zpn_org( 13)/       4.66962600D0/
      data     zdn_org( 13)/       7.13113800D0/
      data     alp_org( 13)/       0.96879800D0/
      data     gss_org( 13)/       6.65215500D0/
      data     gsp_org( 13)/       7.45943500D0/
      data     gpp_org( 13)/       7.66885700D0/
      data     gp2_org( 13)/       6.67329900D0/
      data     hsp_org( 13)/       0.43506000D0/
      data gues_org1( 13,1)/       1.00222200D0/
      data gues_org2( 13,1)/       1.51740000D0/
      data gues_org3( 13,1)/       0.65910100D0/
!
!                    Data for Element  14          Silicon
!
      data     uss_org( 14)/     -27.35805800D0/
      data     upp_org( 14)/     -20.49057800D0/
      data     udd_org( 14)/     -22.75190000D0/
      data   betas_org( 14)/      -8.68690900D0/
      data   betap_org( 14)/      -1.85648200D0/
      data   betad_org( 14)/      -6.36062700D0/
      data      zs_org( 14)/       1.75274100D0/
      data      zp_org( 14)/       1.19841300D0/
      data      zd_org( 14)/       2.12859300D0/
      data     zsn_org( 14)/       8.38811100D0/
      data     zpn_org( 14)/       1.84304800D0/
      data     zdn_org( 14)/       0.70860000D0/
      data     gss_org( 14)/       5.19480500D0/
      data     gsp_org( 14)/       5.09053400D0/
      data     gpp_org( 14)/       5.18515000D0/
      data     gp2_org( 14)/       4.76977500D0/
      data     hsp_org( 14)/       1.42501200D0/
      data   polvo_org( 14)/       1.88611000D0/
      data gues_org1( 14,1)/       0.20857100D0/
      data gues_org2( 14,1)/       6.00048300D0/
      data gues_org3( 14,1)/       1.18524500D0/
!
!                    Data for Element  15       Phosphorus
!
      data     uss_org( 15)/     -47.86742323D0/
      data     upp_org( 15)/     -39.20419043D0/
      data     udd_org( 15)/      -7.37596714D0/
      data   betas_org( 15)/      -9.80027923D0/
      data   betap_org( 15)/     -12.28361682D0/
      data   betad_org( 15)/     -22.86588191D0/
      data      zs_org( 15)/       2.04739172D0/
      data      zp_org( 15)/       1.82112954D0/
      data      zd_org( 15)/       1.27956078D0/
      data     zsn_org( 15)/       3.26632817D0/
      data     zpn_org( 15)/       1.32382922D0/
      data     zdn_org( 15)/       7.00631004D0/
      data     gss_org( 15)/       8.62371066D0/
      data     gsp_org( 15)/       8.30807318D0/
      data     gpp_org( 15)/       8.70036212D0/
      data     gp2_org( 15)/       7.86508651D0/
      data     hsp_org( 15)/       1.84450770D0/
      data   polvo_org( 15)/       2.31461000D0/
      data gues_org1( 15,1)/      -0.01765766D0/
      data gues_org2( 15,1)/       4.32287106D0/
      data gues_org3( 15,1)/       2.47010191D0/
!
!                    Data for Element  16           Sulfur
!
      data     uss_org( 16)/     -52.81945777D0/
      data     upp_org( 16)/     -36.87335599D0/
      data     udd_org( 16)/     -46.14171619D0/
      data   betas_org( 16)/     -13.40463493D0/
      data   betap_org( 16)/      -9.26738586D0/
      data   betad_org( 16)/     -12.98139394D0/
      data      zs_org( 16)/       1.84968731D0/
      data      zp_org( 16)/       1.94369437D0/
      data      zd_org( 16)/       2.88026287D0/
      data     zsn_org( 16)/       1.57989965D0/
      data     zpn_org( 16)/       0.69868617D0/
      data     zdn_org( 16)/       2.30608926D0/
      data     gss_org( 16)/       8.86729470D0/
      data     gsp_org( 16)/       6.19182172D0/
      data     gpp_org( 16)/       7.55957479D0/
      data     gp2_org( 16)/       6.72917030D0/
      data     hsp_org( 16)/       5.33308599D0/
      data   polvo_org( 16)/       1.45331000D0/
      data gues_org1( 16,1)/      -0.11203409D0/
      data gues_org2( 16,1)/       1.02612340D0/
      data gues_org3( 16,1)/       1.61268145D0/
!
!                    Data for Element  17         Chlorine
!
      data     uss_org( 17)/     -66.46128618D0/
      data     upp_org( 17)/     -62.24794228D0/
      data     udd_org( 17)/     -43.46840998D0/
      data   betas_org( 17)/       1.24470145D0/
      data   betap_org( 17)/     -12.27639617D0/
      data   betad_org( 17)/      -4.93751724D0/
      data      zs_org( 17)/       2.37132118D0/
      data      zp_org( 17)/       2.09554394D0/
      data      zd_org( 17)/       1.21389223D0/
      data     zsn_org( 17)/       0.76291736D0/
      data     zpn_org( 17)/       2.19897547D0/
      data     zdn_org( 17)/       8.62093236D0/
      data     gss_org( 17)/      10.99267343D0/
      data     gsp_org( 17)/       9.50390305D0/
      data     gpp_org( 17)/       9.62518424D0/
      data     gp2_org( 17)/       8.38049254D0/
      data     hsp_org( 17)/       2.32744386D0/
      data   polvo_org( 17)/       1.23621000D0/
      data gues_org1( 17,1)/      -0.07781945D0/
      data gues_org2( 17,1)/       0.98368004D0/
      data gues_org3( 17,1)/       1.82890575D0/
!
!                    Data for Element  18            Argon
!
      data     uss_org( 18)/      -7.79793100D0/
      data     upp_org( 18)/     -83.21148700D0/
      data   betas_org( 18)/      -8.83984200D0/
      data   betap_org( 18)/     -28.42730300D0/
      data      zs_org( 18)/       6.00027200D0/
      data      zp_org( 18)/       5.94917000D0/
      data     gss_org( 18)/      17.85877600D0/
      data     gsp_org( 18)/       4.16845100D0/
      data     gpp_org( 18)/      11.85250000D0/
      data     gp2_org( 18)/      15.66954300D0/
      data     hsp_org( 18)/       4.57454900D0/
!
!                    Data for Element  19        Potassium
!
      data     uss_org( 19)/      -3.52807031D0/
      data     upp_org( 19)/       0.24886472D0/
      data   betas_org( 19)/      -5.70393556D0/
      data   betap_org( 19)/      -8.97131501D0/
      data      zs_org( 19)/       7.76683045D0/
      data      zp_org( 19)/       1.16903706D0/
      data     gss_org( 19)/       2.25960503D0/
      data     gsp_org( 19)/       6.68884024D0/
      data     gpp_org( 19)/       5.49931070D0/
      data     gp2_org( 19)/      23.41230341D0/
      data     hsp_org( 19)/       0.00000010D0/
      data gues_org1( 19,1)/       0.23946479D0/
      data gues_org2( 19,1)/       7.66204824D0/
      data gues_org3( 19,1)/       2.07367453D0/
!
!                    Data for Element  20          Calcium
!
      data     uss_org( 20)/     -11.47553748D0/
      data     upp_org( 20)/     -11.34698800D0/
      data   betas_org( 20)/     -13.58930464D0/
      data   betap_org( 20)/      -0.77954077D0/
      data      zs_org( 20)/       1.20624088D0/
      data      zp_org( 20)/       1.77240549D0/
      data     gss_org( 20)/       5.36170799D0/
      data     gsp_org( 20)/       7.01307557D0/
      data     gpp_org( 20)/      10.22493218D0/
      data     gp2_org( 20)/       8.89140215D0/
      data     hsp_org( 20)/       0.69547543D0/
      data gues_org1( 20,1)/      -0.08613230D0/
      data gues_org2( 20,1)/       2.25318729D0/
      data gues_org3( 20,1)/       2.41309577D0/
!
!                    Data for Element  21         Scandium
!
      data     uss_org( 21)/     -15.54446100D0/
      data     upp_org( 21)/     -18.64629500D0/
      data     udd_org( 21)/     -16.06944400D0/
      data   betas_org( 21)/      -8.62094400D0/
      data   betap_org( 21)/       3.07594800D0/
      data   betad_org( 21)/      -9.76866100D0/
      data      zs_org( 21)/       1.40246900D0/
      data      zp_org( 21)/       1.34519600D0/
      data      zd_org( 21)/       1.85901200D0/
      data     zsn_org( 21)/       0.84841800D0/
      data     zpn_org( 21)/       2.45172900D0/
      data     zdn_org( 21)/       0.78937200D0/
      data     alp_org( 21)/       0.81655600D0/
      data     gss_org( 21)/       4.63821534D0/
      data     gsp_org( 21)/       5.73916362D0/
      data     gpp_org( 21)/      14.60487239D0/
      data     gp2_org( 21)/      12.80259529D0/
      data     hsp_org( 21)/       0.19383456D0/
      data    poc__org( 21)/       3.17373400D0/
      data    f0sd_org( 21)/       4.79831300D0/
      data    g2sd_org( 21)/       5.38013600D0/
!
!                    Data for Element  22         Titanium
!
      data     uss_org( 22)/     -25.50797300D0/
      data     upp_org( 22)/     -17.26090900D0/
      data     udd_org( 22)/     -23.80948600D0/
      data   betas_org( 22)/       3.38914200D0/
      data   betap_org( 22)/      -3.35535000D0/
      data   betad_org( 22)/      -1.84282900D0/
      data      zs_org( 22)/       5.32477700D0/
      data      zp_org( 22)/       1.16406800D0/
      data      zd_org( 22)/       1.41828000D0/
      data     zsn_org( 22)/       1.04590400D0/
      data     zpn_org( 22)/       1.07684400D0/
      data     zdn_org( 22)/       0.71794500D0/
      data     gss_org( 22)/       5.71785132D0/
      data     gsp_org( 22)/       5.80001498D0/
      data     gpp_org( 22)/       6.41472577D0/
      data     gp2_org( 22)/       5.62313287D0/
      data     hsp_org( 22)/       1.40373164D0/
      data    f0sd_org( 22)/       6.56056200D0/
      data    g2sd_org( 22)/       3.39623500D0/
!
!                    Data for Element  23         Vanadium
!
      data     uss_org( 23)/     -32.16227600D0/
      data     upp_org( 23)/     -21.57250100D0/
      data     udd_org( 23)/     -34.50624500D0/
      data   betas_org( 23)/      -1.21133000D0/
      data   betap_org( 23)/       0.74074600D0/
      data   betad_org( 23)/       3.15366900D0/
      data      zs_org( 23)/       1.97433000D0/
      data      zp_org( 23)/       1.06310600D0/
      data      zd_org( 23)/       1.39480600D0/
      data     zsn_org( 23)/       1.09442600D0/
      data     zpn_org( 23)/       0.75537800D0/
      data     zdn_org( 23)/       1.09936700D0/
      data     gss_org( 23)/       5.98311618D0/
      data     gsp_org( 23)/       4.73676949D0/
      data     gpp_org( 23)/       4.49976294D0/
      data     gp2_org( 23)/       3.94448115D0/
      data     hsp_org( 23)/       0.90110519D0/
      data    f0sd_org( 23)/       6.81002100D0/
      data    g2sd_org( 23)/       1.83140700D0/
!
!                    Data for Element  24         Chromium
!
      data     uss_org( 24)/     -34.86433900D0/
      data     upp_org( 24)/     -26.97861500D0/
      data     udd_org( 24)/     -54.43103600D0/
      data   betas_org( 24)/      -5.12261500D0/
      data   betap_org( 24)/       3.92671100D0/
      data   betad_org( 24)/      -4.23055000D0/
      data      zs_org( 24)/       3.28346000D0/
      data      zp_org( 24)/       1.02939400D0/
      data      zd_org( 24)/       1.62311900D0/
      data     zsn_org( 24)/       1.61985300D0/
      data     zpn_org( 24)/       0.84826600D0/
      data     zdn_org( 24)/       1.40501500D0/
      data     gss_org( 24)/       8.85557242D0/
      data     gsp_org( 24)/       5.58863066D0/
      data     gpp_org( 24)/       5.05309383D0/
      data     gp2_org( 24)/       4.42952965D0/
      data     hsp_org( 24)/       0.64803936D0/
      data    f0sd_org( 24)/       6.15013600D0/
      data    g2sd_org( 24)/       2.00030000D0/
!
!                    Data for Element  25        Manganese
!
      data     uss_org( 25)/     -51.46000000D0/
      data     upp_org( 25)/     -37.54399000D0/
      data     udd_org( 25)/     -47.65537000D0/
      data   betas_org( 25)/      -4.18529000D0/
      data   betap_org( 25)/      -3.47963000D0/
      data   betad_org( 25)/     -13.47319000D0/
      data      zs_org( 25)/       2.13168000D0/
      data      zp_org( 25)/       1.52588000D0/
      data      zd_org( 25)/       2.60780000D0/
      data     zsn_org( 25)/       1.13245000D0/
      data     zpn_org( 25)/       1.39074000D0/
      data     zdn_org( 25)/       0.96255000D0/
      data     gss_org( 25)/       6.19098954D0/
      data     gsp_org( 25)/       6.75742701D0/
      data     gpp_org( 25)/       8.28459435D0/
      data     gp2_org( 25)/       7.26225507D0/
      data     hsp_org( 25)/       1.52051823D0/
      data    f0sd_org( 25)/       7.69092000D0/
      data    g2sd_org( 25)/       1.10533000D0/
!
!                    Data for Element  26             Iron
!
      data     uss_org( 26)/     -63.58803126D0/
      data     upp_org( 26)/     -56.75246939D0/
      data     udd_org( 26)/     -86.87882209D0/
      data   betas_org( 26)/      -7.00398257D0/
      data   betap_org( 26)/      -6.50486868D0/
      data   betad_org( 26)/      -5.50434298D0/
      data      zs_org( 26)/       1.28927823D0/
      data      zp_org( 26)/       4.34287794D0/
      data      zd_org( 26)/       1.38098203D0/
      data     zsn_org( 26)/       1.19929712D0/
      data     zpn_org( 26)/       1.19788000D0/
      data     zdn_org( 26)/       1.67720640D0/
      data     gss_org( 26)/       6.55643598D0/
      data     gsp_org( 26)/       6.55255788D0/
      data     gpp_org( 26)/       7.13573341D0/
      data     gp2_org( 26)/       6.25516639D0/
      data     hsp_org( 26)/       1.58823536D0/
      data    poc__org( 26)/       1.45706296D0/
      data    f0sd_org( 26)/       7.54776190D0/
      data    g2sd_org( 26)/       1.39266710D0/
!
!                    Data for Element  27           Cobalt
!
      data     uss_org( 27)/     -21.23671841D0/
      data     upp_org( 27)/       9.84676450D0/
      data     udd_org( 27)/     -28.26385032D0/
      data   betas_org( 27)/      -5.78517148D0/
      data   betap_org( 27)/      -3.49202346D0/
      data   betad_org( 27)/      -2.68166071D0/
      data      zs_org( 27)/       1.33517899D0/
      data      zp_org( 27)/       1.25921162D0/
      data      zd_org( 27)/       1.70611919D0/
      data     zsn_org( 27)/       0.53706986D0/
      data     zpn_org( 27)/       1.67470244D0/
      data     zdn_org( 27)/       0.36428989D0/
      data     gss_org( 27)/       2.93610657D0/
      data     gsp_org( 27)/       3.63927114D0/
      data     gpp_org( 27)/       9.97614966D0/
      data     gp2_org( 27)/       8.74506830D0/
      data     hsp_org( 27)/       0.09223269D0/
      data    f0sd_org( 27)/       1.28975559D0/
      data    g2sd_org( 27)/       1.76910311D0/
!
!                    Data for Element  28           Nickel
!
      data     uss_org( 28)/     -47.62024700D0/
      data     upp_org( 28)/     -32.87840800D0/
      data     udd_org( 28)/     -93.02639500D0/
      data   betas_org( 28)/      -9.15152100D0/
      data   betap_org( 28)/      -8.08669600D0/
      data   betad_org( 28)/      -8.65591000D0/
      data      zs_org( 28)/       1.59182800D0/
      data      zp_org( 28)/       2.30473900D0/
      data      zd_org( 28)/       2.51476100D0/
      data     zsn_org( 28)/       0.74647000D0/
      data     zpn_org( 28)/       0.75332700D0/
      data     zdn_org( 28)/       1.46134500D0/
      data     alp_org( 28)/       2.89496000D0/
      data     gss_org( 28)/       4.08087594D0/
      data     gsp_org( 28)/       4.09945168D0/
      data     gpp_org( 28)/       4.48754520D0/
      data     gp2_org( 28)/       3.93377111D0/
      data     hsp_org( 28)/       0.99349774D0/
      data    poc__org( 28)/       1.58697900D0/
      data    f0sd_org( 28)/       4.65166400D0/
      data    g2sd_org( 28)/       1.88050200D0/
!
!                    Data for Element  29           Copper
!
      data     uss_org( 29)/     -97.00220500D0/
      data     upp_org( 29)/      -1.00000000D0/
      data     udd_org( 29)/    -110.44259200D0/
      data   betas_org( 29)/      -9.36950800D0/
      data   betap_org( 29)/      -0.10000000D0/
      data   betad_org( 29)/     -16.98209200D0/
      data      zs_org( 29)/       1.66909600D0/
      data      zp_org( 29)/       3.00000000D0/
      data      zd_org( 29)/       2.73499000D0/
      data     zsn_org( 29)/       1.89959800D0/
      data     zpn_org( 29)/       3.00000000D0/
      data     zdn_org( 29)/       1.48431700D0/
      data     gss_org( 29)/      10.38491002D0/
      data     gsp_org( 29)/      12.14536058D0/
      data     gpp_org( 29)/      17.87090545D0/
      data     gp2_org( 29)/      15.66559186D0/
      data     hsp_org( 29)/       2.03739400D0/
      data    f0sd_org( 29)/       9.84880700D0/
      data    g2sd_org( 29)/       9.84757700D0/
!
!                    Data for Element  30             Zinc
!
      data     uss_org( 30)/     -18.08806266D0/
      data     upp_org( 30)/      -9.67029382D0/
      data   betas_org( 30)/     -13.10860471D0/
      data   betap_org( 30)/       1.85795151D0/
      data      zs_org( 30)/       1.59949465D0/
      data      zp_org( 30)/       1.60377982D0/
      data     gss_org( 30)/       8.79420586D0/
      data     gsp_org( 30)/       4.78518617D0/
      data     gpp_org( 30)/      14.37082405D0/
      data     gp2_org( 30)/       6.06638427D0/
      data     hsp_org( 30)/       0.16164719D0/
!
!                    Data for Element  31          Gallium
!
      data     uss_org( 31)/     -30.60022600D0/
      data     upp_org( 31)/     -21.03242500D0/
      data   betas_org( 31)/     -10.80832000D0/
      data   betap_org( 31)/      -4.18550000D0/
      data      zs_org( 31)/       2.33906700D0/
      data      zp_org( 31)/       1.72959200D0/
      data     gss_org( 31)/      10.35488500D0/
      data     gsp_org( 31)/       7.99367400D0/
      data     gpp_org( 31)/       6.09018400D0/
      data     gp2_org( 31)/       6.29922600D0/
      data     hsp_org( 31)/       1.29597400D0/
!
!                    Data for Element  32        Germanium
!
      data     uss_org( 32)/     -32.74733800D0/
      data     upp_org( 32)/     -24.70901600D0/
      data   betas_org( 32)/     -14.85429700D0/
      data   betap_org( 32)/      -2.59126000D0/
      data      zs_org( 32)/       2.54607300D0/
      data      zp_org( 32)/       1.70913000D0/
      data     gss_org( 32)/       7.51830100D0/
      data     gsp_org( 32)/       6.59444300D0/
      data     gpp_org( 32)/       6.06680100D0/
      data     gp2_org( 32)/       5.30594700D0/
      data     hsp_org( 32)/       0.29074200D0/
!
!                    Data for Element  33          Arsenic
!
      data     uss_org( 33)/     -37.95696500D0/
      data     upp_org( 33)/     -38.45370100D0/
      data     udd_org( 33)/     -30.28265800D0/
      data   betas_org( 33)/     -11.96372500D0/
      data   betap_org( 33)/      -7.34007300D0/
      data   betad_org( 33)/       3.75300500D0/
      data      zs_org( 33)/       2.92617100D0/
      data      zp_org( 33)/       1.76519100D0/
      data      zd_org( 33)/       1.39214200D0/
      data     zsn_org( 33)/       2.00654300D0/
      data     zpn_org( 33)/       3.31683200D0/
      data     zdn_org( 33)/       4.65344000D0/
      data     gss_org( 33)/       6.66503000D0/
      data     gsp_org( 33)/       6.21386700D0/
      data     gpp_org( 33)/       9.31083600D0/
      data     gp2_org( 33)/       8.71254200D0/
      data     hsp_org( 33)/       0.28066200D0/
!
!                    Data for Element  34         Selenium
!
      data     uss_org( 34)/     -37.28591885D0/
      data     upp_org( 34)/     -33.59240043D0/
      data   betas_org( 34)/      -6.24107962D0/
      data   betap_org( 34)/      -5.03066509D0/
      data      zs_org( 34)/       2.14361983D0/
      data      zp_org( 34)/       1.83381324D0/
      data     gss_org( 34)/       5.24788145D0/
      data     gsp_org( 34)/       3.97949564D0/
      data     gpp_org( 34)/       5.90315008D0/
      data     gp2_org( 34)/       5.90414302D0/
      data     hsp_org( 34)/       0.53718015D0/
!
!                    Data for Element  35          Bromine
!
      data     uss_org( 35)/     -58.74432465D0/
      data     upp_org( 35)/     -48.18471406D0/
      data     udd_org( 35)/      11.72958259D0/
      data   betas_org( 35)/     -28.01073019D0/
      data   betap_org( 35)/     -12.42223445D0/
      data   betad_org( 35)/     -12.37839387D0/
      data      zs_org( 35)/       3.43573316D0/
      data      zp_org( 35)/       2.31387649D0/
      data      zd_org( 35)/       2.25611310D0/
      data     zsn_org( 35)/       0.98875815D0/
      data     zpn_org( 35)/       6.76856856D0/
      data     zdn_org( 35)/       1.74468307D0/
      data     gss_org( 35)/      12.80570147D0/
      data     gsp_org( 35)/       6.74334750D0/
      data     gpp_org( 35)/       9.47123235D0/
      data     gp2_org( 35)/       8.23187235D0/
      data     hsp_org( 35)/       9.71494721D0/
      data   polvo_org( 35)/       2.14242000D0/
      data gues_org1( 35,1)/      -0.00113945D0/
      data gues_org2( 35,1)/       0.87257533D0/
      data gues_org3( 35,1)/       3.42097381D0/
!
!                    Data for Element  36          Krypton
!
      data     uss_org( 36)/       8.53538400D0/
      data     upp_org( 36)/     -80.48432100D0/
      data   betas_org( 36)/      -2.72708800D0/
      data   betap_org( 36)/     -16.14295100D0/
      data      zs_org( 36)/       1.31224800D0/
      data      zp_org( 36)/       4.49137100D0/
      data     gss_org( 36)/      19.99985700D0/
      data     gsp_org( 36)/       1.17530400D0/
      data     gpp_org( 36)/       9.17478400D0/
      data     gp2_org( 36)/      14.92694800D0/
      data     hsp_org( 36)/       0.29986700D0/
!
!                    Data for Element  37         Rubidium
!
      data     uss_org( 37)/      -3.63650500D0/
      data     upp_org( 37)/      -2.50067100D0/
      data   betas_org( 37)/       9.99874400D0/
      data   betap_org( 37)/       1.34300400D0/
      data      zs_org( 37)/       5.51014500D0/
      data      zp_org( 37)/       1.33517000D0/
      data     gss_org( 37)/       6.68082400D0/
      data     gsp_org( 37)/      20.00109800D0/
      data     gpp_org( 37)/       5.06887400D0/
      data     gp2_org( 37)/       2.74786000D0/
      data     hsp_org( 37)/       3.60283400D0/
!
!                    Data for Element  38        Strontium
!
      data     uss_org( 38)/     -10.42767100D0/
      data     upp_org( 38)/      -9.94375100D0/
      data   betas_org( 38)/      -6.25310800D0/
      data   betap_org( 38)/      -9.84449800D0/
      data      zs_org( 38)/       2.19730300D0/
      data      zp_org( 38)/       1.73013700D0/
      data     gss_org( 38)/       4.60366400D0/
      data     gsp_org( 38)/       5.71606900D0/
      data     gpp_org( 38)/       7.33462000D0/
      data     gp2_org( 38)/       7.44308800D0/
      data     hsp_org( 38)/       0.83152700D0/
      data gues_org1( 38,1)/      -0.01294800D0/
      data gues_org2( 38,1)/       6.00012600D0/
      data gues_org3( 38,1)/       3.01196400D0/
!
!                    Data for Element  39          Yttrium
!
      data     uss_org( 39)/     -14.24780900D0/
      data     upp_org( 39)/     -14.81714000D0/
      data     udd_org( 39)/     -16.39430200D0/
      data   betas_org( 39)/       0.34333600D0/
      data   betap_org( 39)/      -3.18080700D0/
      data   betad_org( 39)/      -4.50895700D0/
      data      zs_org( 39)/       0.59336800D0/
      data      zp_org( 39)/       1.49042200D0/
      data      zd_org( 39)/       1.65089300D0/
      data     zsn_org( 39)/       0.90261100D0/
      data     zpn_org( 39)/       1.48440000D0/
      data     zdn_org( 39)/       1.38423800D0/
      data     alp_org( 39)/       0.50072700D0/
      data     gss_org( 39)/       4.04673328D0/
      data     gsp_org( 39)/       4.72627720D0/
      data     gpp_org( 39)/       7.27875241D0/
      data     gp2_org( 39)/       6.34328112D0/
      data     hsp_org( 39)/       0.67922805D0/
      data    poc__org( 39)/       2.77370300D0/
      data    f0sd_org( 39)/       4.97271600D0/
      data    g2sd_org( 39)/       5.01636400D0/
!
!                    Data for Element  40        Zirconium
!
      data     uss_org( 40)/     -20.00888400D0/
      data     upp_org( 40)/     -14.55969200D0/
      data     udd_org( 40)/     -21.30265700D0/
      data   betas_org( 40)/       9.55195200D0/
      data   betap_org( 40)/      -4.55191500D0/
      data   betad_org( 40)/      -3.21327400D0/
      data      zs_org( 40)/       1.69259000D0/
      data      zp_org( 40)/       1.69491600D0/
      data      zd_org( 40)/       1.56739200D0/
      data     zsn_org( 40)/       1.18910900D0/
      data     zpn_org( 40)/       0.80909200D0/
      data     zdn_org( 40)/       1.19024900D0/
      data     gss_org( 40)/       5.33120797D0/
      data     gsp_org( 40)/       4.15057919D0/
      data     gpp_org( 40)/       3.96738099D0/
      data     gp2_org( 40)/       3.45748990D0/
      data     hsp_org( 40)/       0.74367621D0/
      data    f0sd_org( 40)/       5.01070400D0/
      data    g2sd_org( 40)/       2.94365200D0/
!
!                    Data for Element  41          Niobium
!
      data     uss_org( 41)/     -31.26929800D0/
      data     upp_org( 41)/     -20.15127700D0/
      data     udd_org( 41)/     -35.89311600D0/
      data   betas_org( 41)/     -12.04524400D0/
      data   betap_org( 41)/       1.46576200D0/
      data   betad_org( 41)/      -5.92016000D0/
      data      zs_org( 41)/       2.35556200D0/
      data      zp_org( 41)/       1.38690700D0/
      data      zd_org( 41)/       1.97732400D0/
      data     zsn_org( 41)/       1.49075400D0/
      data     zpn_org( 41)/       0.89276000D0/
      data     zdn_org( 41)/       1.44383700D0/
      data     alp_org( 41)/       0.84397400D0/
      data     gss_org( 41)/       6.68359218D0/
      data     gsp_org( 41)/       4.68533917D0/
      data     gpp_org( 41)/       4.37764686D0/
      data     gp2_org( 41)/       3.81502806D0/
      data     hsp_org( 41)/       0.65067949D0/
      data    f0sd_org( 41)/       6.55067400D0/
      data    g2sd_org( 41)/       1.06557700D0/
!
!                    Data for Element  42       Molybdenum
!
      data     uss_org( 42)/     -53.46772800D0/
      data     upp_org( 42)/     -35.29195100D0/
      data     udd_org( 42)/     -55.83697700D0/
      data   betas_org( 42)/      -0.18934400D0/
      data   betap_org( 42)/       7.01776200D0/
      data   betad_org( 42)/     -10.94112600D0/
      data      zs_org( 42)/       1.06042900D0/
      data      zp_org( 42)/       1.35041200D0/
      data      zd_org( 42)/       1.82715200D0/
      data     zsn_org( 42)/       1.91299500D0/
      data     zpn_org( 42)/       1.35505500D0/
      data     zdn_org( 42)/       1.87623100D0/
      data     gss_org( 42)/       8.57665210D0/
      data     gsp_org( 42)/       6.88829273D0/
      data     gpp_org( 42)/       6.64450947D0/
      data     gp2_org( 42)/       5.79055161D0/
      data     hsp_org( 42)/       1.31736764D0/
      data    f0sd_org( 42)/      10.00060800D0/
      data    g2sd_org( 42)/       1.21675200D0/
!
!                    Data for Element  43       Technetium
!
      data     uss_org( 43)/     -41.85029200D0/
      data     upp_org( 43)/     -34.91029300D0/
      data     udd_org( 43)/     -45.53041200D0/
      data   betas_org( 43)/      -2.79102400D0/
      data   betap_org( 43)/      -8.08669700D0/
      data   betad_org( 43)/      -5.72433500D0/
      data      zs_org( 43)/       1.95624500D0/
      data      zp_org( 43)/       6.00629900D0/
      data      zd_org( 43)/       1.76736000D0/
      data     zsn_org( 43)/       1.41103300D0/
      data     zpn_org( 43)/       1.14131300D0/
      data     zdn_org( 43)/       1.15931200D0/
      data     gss_org( 43)/       6.32617395D0/
      data     gsp_org( 43)/       5.58713806D0/
      data     gpp_org( 43)/       5.59642600D0/
      data     gp2_org( 43)/       4.87716869D0/
      data     hsp_org( 43)/       1.25898912D0/
      data    f0sd_org( 43)/       5.43488600D0/
      data    g2sd_org( 43)/       1.10687500D0/
!
!                    Data for Element  44        Ruthenium
!
      data     uss_org( 44)/     -44.90152100D0/
      data     upp_org( 44)/     -41.42440900D0/
      data     udd_org( 44)/     -37.93451400D0/
      data   betas_org( 44)/     -12.85950800D0/
      data   betap_org( 44)/      -8.47551800D0/
      data   betad_org( 44)/      -3.83079700D0/
      data      zs_org( 44)/       1.45919500D0/
      data      zp_org( 44)/       5.53720100D0/
      data      zd_org( 44)/       2.09316400D0/
      data     zsn_org( 44)/       0.98444900D0/
      data     zpn_org( 44)/       4.58661300D0/
      data     zdn_org( 44)/       0.76533200D0/
      data     gss_org( 44)/       4.41364279D0/
      data     gsp_org( 44)/       5.35699582D0/
      data     gpp_org( 44)/      22.49044761D0/
      data     gp2_org( 44)/      19.59995666D0/
      data     hsp_org( 44)/       0.00805809D0/
      data    f0sd_org( 44)/       5.91740400D0/
      data    g2sd_org( 44)/       5.85973800D0/
!
!                    Data for Element  45          Rhodium
!
      data     uss_org( 45)/     -20.51375600D0/
      data     upp_org( 45)/     -40.04543100D0/
      data     udd_org( 45)/     -35.81849200D0/
      data   betas_org( 45)/      -8.22214100D0/
      data   betap_org( 45)/     -15.55669100D0/
      data   betad_org( 45)/     -13.39618200D0/
      data      zs_org( 45)/       1.32491900D0/
      data      zp_org( 45)/       4.30611100D0/
      data      zd_org( 45)/       2.90140600D0/
      data     zsn_org( 45)/       0.80992300D0/
      data     zpn_org( 45)/       6.89825900D0/
      data     zdn_org( 45)/       0.64313400D0/
      data     gss_org( 45)/       3.63117927D0/
      data     gsp_org( 45)/       4.40781970D0/
      data     gpp_org( 45)/      33.82559912D0/
      data     gp2_org( 45)/      29.47830511D0/
      data     hsp_org( 45)/       0.00009171D0/
      data    f0sd_org( 45)/       1.77549700D0/
      data    g2sd_org( 45)/       1.85157100D0/
!
!                    Data for Element  46        Palladium
!
      data     uss_org( 46)/     -76.14019600D0/
      data     upp_org( 46)/     -21.07336200D0/
      data     udd_org( 46)/     -85.32530100D0/
      data   betas_org( 46)/      -8.03824500D0/
      data   betap_org( 46)/       0.74003700D0/
      data   betad_org( 46)/      -2.39449800D0/
      data      zs_org( 46)/       1.65850300D0/
      data      zp_org( 46)/       1.15671800D0/
      data      zd_org( 46)/       2.21986100D0/
      data     zsn_org( 46)/       1.79408500D0/
      data     zpn_org( 46)/       6.15877800D0/
      data     zdn_org( 46)/       1.63091300D0/
      data     gss_org( 46)/       8.04353534D0/
      data     gsp_org( 46)/       9.75504166D0/
      data     gpp_org( 46)/      30.19955553D0/
      data     gp2_org( 46)/      26.31828364D0/
      data     hsp_org( 46)/       0.08612136D0/
      data    f0sd_org( 46)/       8.00444700D0/
      data    g2sd_org( 46)/       2.61314800D0/
!
!                    Data for Element  47           Silver
!
      data     uss_org( 47)/     -25.48413700D0/
      data     upp_org( 47)/     -36.11602300D0/
      data     udd_org( 47)/     -35.66827200D0/
      data   betas_org( 47)/      -6.12962300D0/
      data   betap_org( 47)/       1.00411500D0/
      data   betad_org( 47)/     -69.23834700D0/
      data      zs_org( 47)/       1.99400400D0/
      data      zp_org( 47)/       0.68181700D0/
      data      zd_org( 47)/       6.00732800D0/
      data     zsn_org( 47)/       0.69551400D0/
      data     zpn_org( 47)/       4.72994900D0/
      data     zdn_org( 47)/       0.50652200D0/
      data     gss_org( 47)/       3.11824213D0/
      data     gsp_org( 47)/       3.78515183D0/
      data     gpp_org( 47)/      23.19329540D0/
      data     gp2_org( 47)/      20.21247387D0/
      data     hsp_org( 47)/       0.00043152D0/
      data    f0sd_org( 47)/       1.93832700D0/
      data    g2sd_org( 47)/       1.07190100D0/
!
!                    Data for Element  48          Cadmium
!
      data     uss_org( 48)/     -14.64579200D0/
      data     upp_org( 48)/      -9.31866400D0/
      data   betas_org( 48)/     -11.61318300D0/
      data   betap_org( 48)/       1.66317800D0/
      data      zs_org( 48)/       1.38410800D0/
      data      zp_org( 48)/       1.95741300D0/
      data     gss_org( 48)/       6.67728400D0/
      data     gsp_org( 48)/       5.95337300D0/
      data     gpp_org( 48)/      18.72984300D0/
      data     gp2_org( 48)/       9.91745200D0/
      data     hsp_org( 48)/       0.82519200D0/
!
!                    Data for Element  49           Indium
!
      data     uss_org( 49)/     -28.33924600D0/
      data     upp_org( 49)/     -23.37387500D0/
      data   betas_org( 49)/      -1.98237600D0/
      data   betap_org( 49)/      -3.33029400D0/
      data      zs_org( 49)/       2.02308700D0/
      data      zp_org( 49)/       2.10661800D0/
      data     gss_org( 49)/       9.90609100D0/
      data     gsp_org( 49)/      10.52006000D0/
      data     gpp_org( 49)/       4.82600600D0/
      data     gp2_org( 49)/       7.90656300D0/
      data     hsp_org( 49)/       3.50029900D0/
!
!                    Data for Element  50              Tin
!
      data     uss_org( 50)/     -29.88821700D0/
      data     upp_org( 50)/     -22.15695400D0/
      data   betas_org( 50)/      -8.62108700D0/
      data   betap_org( 50)/      -4.98975200D0/
      data      zs_org( 50)/       2.38394100D0/
      data      zp_org( 50)/       2.05790800D0/
      data     gss_org( 50)/       8.26965500D0/
      data     gsp_org( 50)/       5.01334900D0/
      data     gpp_org( 50)/       6.58487400D0/
      data     gp2_org( 50)/       5.85515900D0/
      data     hsp_org( 50)/       0.53121200D0/
      data gues_org1( 50,1)/      -1.00458700D0/
      data gues_org2( 50,1)/       4.70625200D0/
      data gues_org3( 50,1)/       1.18021800D0/
!
!                    Data for Element  51         Antimony
!
      data     uss_org( 51)/     -41.68887900D0/
      data     upp_org( 51)/     -39.54118000D0/
      data     udd_org( 51)/      -6.58166300D0/
      data   betas_org( 51)/      -7.47232200D0/
      data   betap_org( 51)/      -5.94075000D0/
      data   betad_org( 51)/      -3.97910800D0/
      data      zs_org( 51)/       2.39117800D0/
      data      zp_org( 51)/       1.77300600D0/
      data      zd_org( 51)/       2.46559000D0/
      data     zsn_org( 51)/       5.99359100D0/
      data     zpn_org( 51)/       6.14508600D0/
      data     zdn_org( 51)/       5.70403100D0/
      data     gss_org( 51)/      10.58883200D0/
      data     gsp_org( 51)/       7.31002300D0/
      data     gpp_org( 51)/       9.28160900D0/
      data     gp2_org( 51)/       8.95408100D0/
      data     hsp_org( 51)/       0.77911200D0/
!
!                    Data for Element  52        Tellurium
!
      data     uss_org( 52)/    -114.73331600D0/
      data     upp_org( 52)/     -50.09638900D0/
      data   betas_org( 52)/     -70.00106200D0/
      data   betap_org( 52)/      -6.15164200D0/
      data      zs_org( 52)/       2.76986200D0/
      data      zp_org( 52)/       1.73131900D0/
      data     gss_org( 52)/       7.03062600D0/
      data     gsp_org( 52)/      12.60138900D0/
      data     gpp_org( 52)/       7.88347900D0/
      data     gp2_org( 52)/       6.97316300D0/
      data     hsp_org( 52)/       5.00082600D0/
!
!                    Data for Element  53           Iodine
!
      data     uss_org( 53)/     -61.73780525D0/
      data     upp_org( 53)/     -54.91097768D0/
      data     udd_org( 53)/     -65.06725189D0/
      data   betas_org( 53)/     -31.10448815D0/
      data   betap_org( 53)/     -10.41671111D0/
      data   betad_org( 53)/      -4.83328024D0/
      data      zs_org( 53)/       3.80693109D0/
      data      zp_org( 53)/       2.10056853D0/
      data      zd_org( 53)/       1.66121671D0/
      data     zsn_org( 53)/       7.75109882D0/
      data     zpn_org( 53)/       8.70927258D0/
      data     zdn_org( 53)/       2.35542871D0/
      data     gss_org( 53)/       3.45548116D0/
      data     gsp_org( 53)/       9.51671899D0/
      data     gpp_org( 53)/       8.92941421D0/
      data     gp2_org( 53)/       7.08111703D0/
      data     hsp_org( 53)/       2.82920092D0/
      data   polvo_org( 53)/       3.82316000D0/
      data gues_org1( 53,1)/       0.24432378D0/
      data gues_org2( 53,1)/       1.89315782D0/
      data gues_org3( 53,1)/       1.25922210D0/
!
!                    Data for Element  54            Xenon
!
      data     uss_org( 54)/     -18.27022700D0/
      data     upp_org( 54)/    -167.16306300D0/
      data   betas_org( 54)/      -3.98062200D0/
      data   betap_org( 54)/     -38.82279200D0/
      data      zs_org( 54)/       2.75978700D0/
      data      zp_org( 54)/       1.97744600D0/
      data     gss_org( 54)/      20.00025200D0/
      data     gsp_org( 54)/       4.17590200D0/
      data     gpp_org( 54)/       2.30578700D0/
      data     gp2_org( 54)/       4.06322000D0/
      data     hsp_org( 54)/       4.41884300D0/
!
!                    Data for Element  55           Cesium
!
      data     uss_org( 55)/      -3.74860900D0/
      data     upp_org( 55)/      -2.34810900D0/
      data   betas_org( 55)/       2.28783800D0/
      data   betap_org( 55)/      -5.90807100D0/
      data      zs_org( 55)/       5.95600800D0/
      data      zp_org( 55)/       1.61948500D0/
      data     gss_org( 55)/       6.46475100D0/
      data     gsp_org( 55)/       4.00450100D0/
      data     gpp_org( 55)/      13.77539000D0/
      data     gp2_org( 55)/      12.91253700D0/
      data     hsp_org( 55)/       1.02692800D0/
!
!                    Data for Element  56           Barium
!
      data     uss_org( 56)/      -9.30698500D0/
      data     upp_org( 56)/      -8.82671300D0/
      data   betas_org( 56)/      10.00312500D0/
      data   betap_org( 56)/      -6.33516000D0/
      data      zs_org( 56)/       1.39537900D0/
      data      zp_org( 56)/       1.43013900D0/
      data     gss_org( 56)/       3.60082300D0/
      data     gsp_org( 56)/       4.74057900D0/
      data     gpp_org( 56)/       3.34516600D0/
      data     gp2_org( 56)/       3.14278300D0/
      data     hsp_org( 56)/       0.92942900D0/
!
!                    Data for Element  57        Lanthanum
!
      data     uss_org( 57)/     -19.64195300D0/
      data     upp_org( 57)/     -22.05943100D0/
      data     udd_org( 57)/     -22.63898600D0/
      data   betas_org( 57)/       0.79672700D0/
      data   betap_org( 57)/     -10.85605600D0/
      data   betad_org( 57)/      -0.48492200D0/
      data      zs_org( 57)/       2.67378000D0/
      data      zp_org( 57)/       1.24819200D0/
      data      zd_org( 57)/       1.68856200D0/
      data     zsn_org( 57)/       1.61778400D0/
      data     zpn_org( 57)/       4.33162000D0/
      data     zdn_org( 57)/       2.28573800D0/
      data     alp_org( 57)/       5.94044300D0/
      data     gss_org( 57)/       6.15444012D0/
      data     gsp_org( 57)/       7.32270386D0/
      data     gpp_org( 57)/      18.07746454D0/
      data     gp2_org( 57)/      15.67905690D0/
      data     hsp_org( 57)/       0.13860106D0/
      data    poc__org( 57)/       2.51170100D0/
      data    f0sd_org( 57)/       8.85685800D0/
      data    g2sd_org( 57)/       7.92558500D0/
!
!                    Data for Element  71         Lutetium
!
      data     uss_org( 71)/     -15.95499400D0/
      data     upp_org( 71)/     -11.60621300D0/
      data     udd_org( 71)/     -13.05005600D0/
      data   betas_org( 71)/      -5.59077800D0/
      data   betap_org( 71)/      -0.93767900D0/
      data   betad_org( 71)/      -7.73775200D0/
      data      zs_org( 71)/       5.47174100D0/
      data      zp_org( 71)/       1.71229600D0/
      data      zd_org( 71)/       2.22589200D0/
      data     zsn_org( 71)/       1.63233500D0/
      data     zpn_org( 71)/       4.03312800D0/
      data     zdn_org( 71)/       0.92199900D0/
      data     gss_org( 71)/       6.20979563D0/
      data     gsp_org( 71)/       7.37910233D0/
      data     gpp_org( 71)/      16.83174619D0/
      data     gp2_org( 71)/      14.59861286D0/
      data     hsp_org( 71)/       0.20900816D0/
      data    poc__org( 71)/       2.74326200D0/
      data    f0sd_org( 71)/       3.92492700D0/
      data    g2sd_org( 71)/       1.00094600D0/
!
!                    Data for Element  72          Hafnium
!
      data     uss_org( 72)/     -22.37514000D0/
      data     upp_org( 72)/     -13.08167000D0/
      data     udd_org( 72)/     -20.63774100D0/
      data   betas_org( 72)/      -5.36635100D0/
      data   betap_org( 72)/     -21.55011900D0/
      data   betad_org( 72)/      -3.88444300D0/
      data      zs_org( 72)/       3.08534400D0/
      data      zp_org( 72)/       1.57581900D0/
      data      zd_org( 72)/       1.84084000D0/
      data     zsn_org( 72)/       0.94692700D0/
      data     zpn_org( 72)/       3.53891100D0/
      data     zdn_org( 72)/       0.94028300D0/
      data     gss_org( 72)/       3.60233846D0/
      data     gsp_org( 72)/       4.29372878D0/
      data     gpp_org( 72)/      14.76919446D0/
      data     gp2_org( 72)/      12.80970790D0/
      data     hsp_org( 72)/       0.01102838D0/
      data    f0sd_org( 72)/       4.84290000D0/
      data    g2sd_org( 72)/       4.38610100D0/
!
!                    Data for Element  73         Tantalum
!
      data     uss_org( 73)/     -39.00998400D0/
      data     upp_org( 73)/       1.16397500D0/
      data     udd_org( 73)/     -43.26631500D0/
      data   betas_org( 73)/     -17.19960500D0/
      data   betap_org( 73)/      -5.81883900D0/
      data   betad_org( 73)/      -9.81679400D0/
      data      zs_org( 73)/       4.57808700D0/
      data      zp_org( 73)/       4.84124400D0/
      data      zd_org( 73)/       1.83824900D0/
      data     zsn_org( 73)/       1.74136700D0/
      data     zpn_org( 73)/       3.43015700D0/
      data     zdn_org( 73)/       2.31119800D0/
      data     gss_org( 73)/       6.62457962D0/
      data     gsp_org( 73)/       7.80532092D0/
      data     gpp_org( 73)/      14.31532349D0/
      data     gp2_org( 73)/      12.41605376D0/
      data     hsp_org( 73)/       0.57726330D0/
      data    f0sd_org( 73)/       8.54442700D0/
      data    g2sd_org( 73)/       2.07425400D0/
!
!                    Data for Element  74         Tungsten
!
      data     uss_org( 74)/     -44.52495000D0/
      data     upp_org( 74)/     -40.01150000D0/
      data     udd_org( 74)/     -46.49041000D0/
      data   betas_org( 74)/     -16.94646000D0/
      data   betap_org( 74)/       5.62317000D0/
      data   betad_org( 74)/      -2.94734000D0/
      data      zs_org( 74)/       2.66456000D0/
      data      zp_org( 74)/       1.62401000D0/
      data      zd_org( 74)/       1.79440000D0/
      data     zsn_org( 74)/       1.49886000D0/
      data     zpn_org( 74)/       1.96590000D0/
      data     zdn_org( 74)/       1.87645000D0/
      data     gss_org( 74)/       5.70202457D0/
      data     gsp_org( 74)/       6.32314460D0/
      data     gpp_org( 74)/       8.20443334D0/
      data     gp2_org( 74)/       7.11591921D0/
      data     hsp_org( 74)/       1.31991216D0/
      data    f0sd_org( 74)/       7.78818000D0/
      data    g2sd_org( 74)/       1.68494000D0/
!
!                    Data for Element  75          Rhenium
!
      data     uss_org( 75)/     -41.29134200D0/
      data     upp_org( 75)/     -35.08959200D0/
      data     udd_org( 75)/     -44.17898500D0/
      data   betas_org( 75)/       3.83007500D0/
      data   betap_org( 75)/      -1.63853000D0/
      data   betad_org( 75)/      -1.41441100D0/
      data      zs_org( 75)/       2.41183900D0/
      data      zp_org( 75)/       1.81535100D0/
      data      zd_org( 75)/       2.52276600D0/
      data     zsn_org( 75)/       1.68082300D0/
      data     zpn_org( 75)/       1.33121800D0/
      data     zdn_org( 75)/       1.49062300D0/
      data     gss_org( 75)/       6.39425566D0/
      data     gsp_org( 75)/       5.55557053D0/
      data     gpp_org( 75)/       5.55566882D0/
      data     gp2_org( 75)/       4.81857660D0/
      data     hsp_org( 75)/       1.22091294D0/
      data    f0sd_org( 75)/       5.44281800D0/
      data    g2sd_org( 75)/       2.37627900D0/
!
!                    Data for Element  76           Osmium
!
      data     uss_org( 76)/     -26.43408000D0/
      data     upp_org( 76)/     -48.73950000D0/
      data     udd_org( 76)/     -55.83788000D0/
      data   betas_org( 76)/     -12.50873000D0/
      data   betap_org( 76)/       0.84688000D0/
      data   betad_org( 76)/       5.16436000D0/
      data      zs_org( 76)/       3.03100000D0/
      data      zp_org( 76)/       1.59396000D0/
      data      zd_org( 76)/       1.77557000D0/
      data     zsn_org( 76)/       1.84470000D0/
      data     zpn_org( 76)/       1.56422000D0/
      data     zdn_org( 76)/       1.77001000D0/
      data     gss_org( 76)/       7.01768325D0/
      data     gsp_org( 76)/       6.38419982D0/
      data     gpp_org( 76)/       6.52807300D0/
      data     gp2_org( 76)/       5.66196813D0/
      data     hsp_org( 76)/       1.50892580D0/
      data    f0sd_org( 76)/       2.02117000D0/
      data    g2sd_org( 76)/       1.39213000D0/
!
!                    Data for Element  77          Iridium
!
      data     uss_org( 77)/     -29.70397400D0/
      data     upp_org( 77)/     -38.21092400D0/
      data     udd_org( 77)/     -32.53820200D0/
      data   betas_org( 77)/     -10.94342700D0/
      data   betap_org( 77)/       2.90888000D0/
      data   betad_org( 77)/      -3.79173100D0/
      data      zs_org( 77)/       1.50090700D0/
      data      zp_org( 77)/       4.10637300D0/
      data      zd_org( 77)/       2.67604700D0/
      data     zsn_org( 77)/       0.92724600D0/
      data     zpn_org( 77)/       3.19189200D0/
      data     zdn_org( 77)/       0.66200700D0/
      data     gss_org( 77)/       3.52746719D0/
      data     gsp_org( 77)/       4.20381991D0/
      data     gpp_org( 77)/      13.32095485D0/
      data     gp2_org( 77)/      11.55361188D0/
      data     hsp_org( 77)/       0.01850052D0/
      data    f0sd_org( 77)/       2.62717000D0/
      data    g2sd_org( 77)/       2.99602900D0/
!
!                    Data for Element  78         Platinum
!
      data     uss_org( 78)/     -73.51617300D0/
      data     upp_org( 78)/     -68.32005600D0/
      data     udd_org( 78)/     -76.59887300D0/
      data   betas_org( 78)/       1.15141800D0/
      data   betap_org( 78)/       3.29869400D0/
      data   betad_org( 78)/     -18.04473700D0/
      data      zs_org( 78)/       2.30126400D0/
      data      zp_org( 78)/       1.66240400D0/
      data      zd_org( 78)/       3.16885200D0/
      data     zsn_org( 78)/       2.27069900D0/
      data     zpn_org( 78)/       1.94989600D0/
      data     zdn_org( 78)/       1.71385600D0/
      data     gss_org( 78)/       8.63828609D0/
      data     gsp_org( 78)/       7.92225366D0/
      data     gpp_org( 78)/       8.13764268D0/
      data     gp2_org( 78)/       7.05798993D0/
      data     hsp_org( 78)/       1.89261714D0/
      data    f0sd_org( 78)/       7.09859100D0/
      data    g2sd_org( 78)/       4.48418300D0/
!
!                    Data for Element  79             Gold
!
      data     uss_org( 79)/     -95.04184600D0/
      data     upp_org( 79)/     -63.89015800D0/
      data     udd_org( 79)/     -88.06608700D0/
      data   betas_org( 79)/      -7.47962500D0/
      data   betap_org( 79)/       3.66435600D0/
      data   betad_org( 79)/     -61.71546800D0/
      data      zs_org( 79)/       1.81416900D0/
      data      zp_org( 79)/       1.61865700D0/
      data      zd_org( 79)/       5.05316700D0/
      data     zsn_org( 79)/       2.44468000D0/
      data     zpn_org( 79)/       7.01499000D0/
      data     zdn_org( 79)/       1.77708900D0/
      data     gss_org( 79)/       9.30015173D0/
      data     gsp_org( 79)/      11.07344344D0/
      data     gpp_org( 79)/      29.27616758D0/
      data     gp2_org( 79)/      25.39198438D0/
      data     hsp_org( 79)/       0.14438427D0/
      data    f0sd_org( 79)/       8.82725700D0/
      data    g2sd_org( 79)/       4.91562500D0/
!
!                    Data for Element  80          Mercury
!
      data     uss_org( 80)/     -17.60873200D0/
      data     upp_org( 80)/     -18.36941700D0/
      data   betas_org( 80)/      -3.04523900D0/
      data   betap_org( 80)/      -5.69355600D0/
      data      zs_org( 80)/       2.10489600D0/
      data      zp_org( 80)/       1.51629300D0/
      data     gss_org( 80)/       6.37282200D0/
      data     gsp_org( 80)/      10.14317600D0/
      data     gpp_org( 80)/      10.39739300D0/
      data     gp2_org( 80)/      14.79405600D0/
      data     hsp_org( 80)/       0.92612800D0/
!
!                    Data for Element  81         Thallium
!
      data     uss_org( 81)/     -29.51862100D0/
      data     upp_org( 81)/     -29.82690700D0/
      data   betas_org( 81)/      -7.23017000D0/
      data   betap_org( 81)/      -7.57554400D0/
      data      zs_org( 81)/       3.33588300D0/
      data      zp_org( 81)/       1.76614100D0/
      data     gss_org( 81)/       5.01511800D0/
      data     gsp_org( 81)/      13.93204900D0/
      data     gpp_org( 81)/      10.49555100D0/
      data     gp2_org( 81)/      10.52619800D0/
      data     hsp_org( 81)/       0.29376000D0/
!
!                    Data for Element  82             Lead
!
      data     uss_org( 82)/     -35.03814500D0/
      data     upp_org( 82)/     -25.41340100D0/
      data   betas_org( 82)/      -8.32379200D0/
      data   betap_org( 82)/      -2.23789100D0/
      data      zs_org( 82)/       2.36890100D0/
      data      zp_org( 82)/       1.68524600D0/
      data     gss_org( 82)/       5.25412800D0/
      data     gsp_org( 82)/       7.06101600D0/
      data     gpp_org( 82)/       6.81855100D0/
      data     gp2_org( 82)/       5.60301900D0/
      data     hsp_org( 82)/       1.01881900D0/
      data gues_org1( 82,1)/      -0.23946300D0/
      data gues_org2( 82,1)/       5.44433800D0/
      data gues_org3( 82,1)/       1.61368200D0/
!
!                    Data for Element  83          Bismuth
!
      data     uss_org( 83)/     -42.40917700D0/
      data     upp_org( 83)/     -36.39374600D0/
      data   betas_org( 83)/     -34.95157800D0/
      data   betap_org( 83)/      -7.35906000D0/
      data      zs_org( 83)/       3.70237700D0/
      data      zp_org( 83)/       1.87232700D0/
      data     gss_org( 83)/       5.85180300D0/
      data     gsp_org( 83)/       6.79058300D0/
      data     gpp_org( 83)/       8.38944200D0/
      data     gp2_org( 83)/       7.72421900D0/
      data     hsp_org( 83)/       0.29560600D0/
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
!                    Data for Element  90          Thorium
!
      data     uss_org( 90)/     -40.56829200D0/
      data     upp_org( 90)/     -28.08918700D0/
      data   betas_org( 90)/      -4.25621800D0/
      data   betap_org( 90)/      -4.25621800D0/
      data      zs_org( 90)/       1.43530600D0/
      data      zp_org( 90)/       1.43530600D0/
      data     gss_org( 90)/       9.82000000D0/
      data     gsp_org( 90)/       8.36000000D0/
      data     gpp_org( 90)/       7.31000000D0/
      data     gp2_org( 90)/       6.54000000D0/
      data     hsp_org( 90)/       1.32000000D0/
!
!                    Data for Element  97        Berkelium
!
      data gues_org1( 97,1)/       1.48000000D0/
      data gues_org2( 97,1)/       0.96000000D0/
      data gues_org1( 97,2)/       1.56000000D0/
      data gues_org2( 97,2)/       0.76000000D0/
      data gues_org1( 97,3)/       1.55000000D0/
      data gues_org2( 97,3)/       0.85000000D0/
!
!                    Data for Element  98          Mithril
!
      data     uss_org( 98)/      -3.00000000D0/
      data   betas_org( 98)/     -99.00000000D0/
      data      zs_org( 98)/       2.00000000D0/
      data     gss_org( 98)/      12.00000000D0/
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
      data  v_par_org(1)/  9.27846500d0/  ! Used in ccrep for scalar correction of C-C triple bonds.
      data  v_par_org(2)/  5.98375200d0/  ! Used in ccrep for exponent correction of C-C triple bonds.
      data  v_par_org(3)/  0.80182891d0/  ! Used in ccrep for scalar correction of O-H term.
      data  v_par_org(4)/  5.06345816d0/  ! Used in ccrep for exponent correction of O-H term.
      data  v_par_org(5)/  1.71474407d0/  ! Used in ccrep for offset correction of O-H term.
      data  v_par_org(7)/  0.52610950d0/  ! Used in dftd3 to set "s6"  in D3H4
      data  v_par_org(8)/ 21.75352387d0/  ! Used in dftd3 to set "alp" in D3H4
      data  v_par_org(9)/  0.68676676d0/  ! Used in dftd3 to set "rs6" in D3H4
      data  v_par_org(11)/ 1.93360969d0/  ! Used in ccrep for exponent correction of C-H term.
      data  v_par_org(12)/ 1.70212501d0/  ! Used in ccrep for offset correction of C-H term.
      data  v_par_org(13)/ 0.00491400d0/  ! Used in ccrep for scalar correction of C-C term.
      data  v_par_org(14)/ 4.79512600d0/  ! Used in ccrep for exponent correction of C-C term.
      data  v_par_org(15)/ 3.13542300d0/  ! Used in ccrep for offset correction of C-C term.
      data  v_par_org(16)/ 2.59803786d0/  ! Used in ccrep for scalar correction of H-H term.
      data  v_par_org(17)/ 6.27278958d0/  ! Used in ccrep for exponent correction of H-H term.
      data  v_par_org(18)/ 1.44426992d0/  ! Used in ccrep for offset correction of H-H term.
      data  v_par_org(19)/ 0.76974534d0/  ! Used in ccrep for scalar correction of C-H term.
      data  v_par_org(20)/ 0.07499093d0/  ! Used in ccrep for scalar correction of C-O term.
      data  v_par_org(21)/15.81072083d0/  ! Used in ccrep for exponent correction of C-O term.
      data  v_par_org(22)/ 2.52666184d0/  ! Used in ccrep for offset correction of C-O term.
      data  v_par_org(23)/ 0.18520409d0/  ! Used in ccrep for scalar correction of S-O term.
      data  v_par_org(24)/ 5.00121981d0/  ! Used in ccrep for exponent correction of S-O term.
      data  v_par_org(25)/ 2.42896530d0/  ! Used in ccrep for offset correction of S-O term.
      data  v_par_org(26)/ 0.04567638d0/  ! Used in ccrep for scalar correction of O-N term.
      data  v_par_org(27)/ 6.14689078d0/  ! Used in ccrep for exponent correction of O-N term.
      data  v_par_org(28)/ 2.10448237d0/  ! Used in ccrep for offset correction of O-N term.
      data  v_par_org(29)/ 0.06917316d0/  ! Used in ccrep for scalar correction of F-H term.
      data  v_par_org(30)/ 5.86395123d0/  ! Used in ccrep for exponent correction of F-H term.
      data  v_par_org(31)/ 2.08820430d0/  ! Used in ccrep for offset correction of F-H term.
  contains
  subroutine alpb_and_xfac_pm6_ORG
    use parameters_C, only : xfac, alpb
 !
      alpb( 1, 1) =     3.86761310d0 !    Hydrogen -     Hydrogen
      xfac( 1, 1) =     2.57164124d0 !    Hydrogen -     Hydrogen
 !
      alpb( 2, 1) =     2.98988100d0 !      Helium -     Hydrogen
      xfac( 2, 1) =     2.37119900d0 !      Helium -     Hydrogen
      alpb( 2, 2) =     3.78355900d0 !      Helium -       Helium
      xfac( 2, 2) =     3.45090000d0 !      Helium -       Helium
 !
      alpb( 3, 1) =     2.13626500d0 !     Lithium -     Hydrogen
      xfac( 3, 1) =     2.19198500d0 !     Lithium -     Hydrogen
      alpb( 3, 2) =     3.11240300d0 !     Lithium -       Helium
      xfac( 3, 2) =     9.27367600d0 !     Lithium -       Helium
      alpb( 3, 3) =     4.71467400d0 !     Lithium -      Lithium
      xfac( 3, 3) =    16.11638400d0 !     Lithium -      Lithium
 !
      alpb( 4, 1) =     2.47541800d0 !   Beryllium -     Hydrogen
      xfac( 4, 1) =     2.56283100d0 !   Beryllium -     Hydrogen
      alpb( 4, 2) =     3.30670200d0 !   Beryllium -       Helium
      xfac( 4, 2) =    12.54487800d0 !   Beryllium -       Helium
      alpb( 4, 3) =     2.23672800d0 !   Beryllium -      Lithium
      xfac( 4, 3) =     3.28716500d0 !   Beryllium -      Lithium
      alpb( 4, 4) =     1.49990700d0 !   Beryllium -    Beryllium
      xfac( 4, 4) =     0.23863300d0 !   Beryllium -    Beryllium
 !
      alpb( 5, 1) =     2.61523100d0 !       Boron -     Hydrogen
      xfac( 5, 1) =     1.32139400d0 !       Boron -     Hydrogen
      alpb( 5, 2) =     3.16314000d0 !       Boron -       Helium
      xfac( 5, 2) =     1.97417000d0 !       Boron -       Helium
      alpb( 5, 3) =     3.75939700d0 !       Boron -      Lithium
      xfac( 5, 3) =     7.88601800d0 !       Boron -      Lithium
      alpb( 5, 4) =     1.88899800d0 !       Boron -    Beryllium
      xfac( 5, 4) =     1.15179200d0 !       Boron -    Beryllium
      alpb( 5, 5) =     3.31862400d0 !       Boron -        Boron
      xfac( 5, 5) =     3.59361900d0 !       Boron -        Boron
 !
      alpb( 6, 1) =     0.99864555d0 !      Carbon -     Hydrogen
      xfac( 6, 1) =     0.20175360d0 !      Carbon -     Hydrogen
      alpb( 6, 2) =     3.04270500d0 !      Carbon -       Helium
      xfac( 6, 2) =     3.21397100d0 !      Carbon -       Helium
      alpb( 6, 3) =     3.24187400d0 !      Carbon -      Lithium
      xfac( 6, 3) =    16.18000200d0 !      Carbon -      Lithium
      alpb( 6, 4) =     4.21288200d0 !      Carbon -    Beryllium
      xfac( 6, 4) =    25.03587900d0 !      Carbon -    Beryllium
      alpb( 6, 5) =     2.91900700d0 !      Carbon -        Boron
      xfac( 6, 5) =     1.87485900d0 !      Carbon -        Boron
      alpb( 6, 6) =     2.57231839d0 !      Carbon -       Carbon
      xfac( 6, 6) =     0.88750253d0 !      Carbon -       Carbon
 !
      alpb( 7, 1) =     0.91440427d0 !    Nitrogen -     Hydrogen
      xfac( 7, 1) =     0.16692722d0 !    Nitrogen -     Hydrogen
      alpb( 7, 2) =     2.81433900d0 !    Nitrogen -       Helium
      xfac( 7, 2) =     1.07786100d0 !    Nitrogen -       Helium
      alpb( 7, 3) =     2.64062300d0 !    Nitrogen -      Lithium
      xfac( 7, 3) =     2.82340300d0 !    Nitrogen -      Lithium
      alpb( 7, 4) =     2.58089500d0 !    Nitrogen -    Beryllium
      xfac( 7, 4) =     1.74060500d0 !    Nitrogen -    Beryllium
      alpb( 7, 5) =     2.47700400d0 !    Nitrogen -        Boron
      xfac( 7, 5) =     0.95288200d0 !    Nitrogen -        Boron
      alpb( 7, 6) =     2.75419777d0 !    Nitrogen -       Carbon
      xfac( 7, 6) =     0.99512536d0 !    Nitrogen -       Carbon
      alpb( 7, 7) =     3.02637026d0 !    Nitrogen -     Nitrogen
      xfac( 7, 7) =     1.11212861d0 !    Nitrogen -     Nitrogen
 !
      alpb( 8, 1) =     1.35972470d0 !      Oxygen -     Hydrogen
      xfac( 8, 1) =     0.15119107d0 !      Oxygen -     Hydrogen
      alpb( 8, 2) =     3.65377500d0 !      Oxygen -       Helium
      xfac( 8, 2) =     6.68452500d0 !      Oxygen -       Helium
      alpb( 8, 3) =     2.58444200d0 !      Oxygen -      Lithium
      xfac( 8, 3) =     1.96859800d0 !      Oxygen -      Lithium
      alpb( 8, 4) =     3.05186700d0 !      Oxygen -    Beryllium
      xfac( 8, 4) =     3.21815500d0 !      Oxygen -    Beryllium
      alpb( 8, 5) =     2.69535100d0 !      Oxygen -        Boron
      xfac( 8, 5) =     1.26980100d0 !      Oxygen -        Boron
      alpb( 8, 6) =     2.96304092d0 !      Oxygen -       Carbon
      xfac( 8, 6) =     0.92566298d0 !      Oxygen -       Carbon
      alpb( 8, 7) =     3.16879398d0 !      Oxygen -     Nitrogen
      xfac( 8, 7) =     0.95083801d0 !      Oxygen -     Nitrogen
      alpb( 8, 8) =     2.70452075d0 !      Oxygen -       Oxygen
      xfac( 8, 8) =     0.43220577d0 !      Oxygen -       Oxygen
 !
      alpb( 9, 1) =     3.24193702d0 !    Fluorine -     Hydrogen
      xfac( 9, 1) =     0.72809348d0 !    Fluorine -     Hydrogen
      alpb( 9, 2) =     2.85654300d0 !    Fluorine -       Helium
      xfac( 9, 2) =     0.74510700d0 !    Fluorine -       Helium
      alpb( 9, 3) =     3.04390100d0 !    Fluorine -      Lithium
      xfac( 9, 3) =     1.97598500d0 !    Fluorine -      Lithium
      alpb( 9, 4) =     3.72692300d0 !    Fluorine -    Beryllium
      xfac( 9, 4) =     3.88299300d0 !    Fluorine -    Beryllium
      alpb( 9, 5) =     2.82383700d0 !    Fluorine -        Boron
      xfac( 9, 5) =     0.86276100d0 !    Fluorine -        Boron
      alpb( 9, 6) =     3.03392715d0 !    Fluorine -       Carbon
      xfac( 9, 6) =     0.78711266d0 !    Fluorine -       Carbon
      alpb( 9, 7) =     3.57433335d0 !    Fluorine -     Nitrogen
      xfac( 9, 7) =     1.48863535d0 !    Fluorine -     Nitrogen
      alpb( 9, 8) =     3.42309777d0 !    Fluorine -       Oxygen
      xfac( 9, 8) =     1.02434853d0 !    Fluorine -       Oxygen
      alpb( 9, 9) =     3.72590569d0 !    Fluorine -     Fluorine
      xfac( 9, 9) =     1.38536051d0 !    Fluorine -     Fluorine
 !
      alpb(10, 1) =     5.99968000d0 !        Neon -     Hydrogen
      xfac(10, 1) =     5.53502100d0 !        Neon -     Hydrogen
      alpb(10, 2) =     3.67775800d0 !        Neon -       Helium
      xfac(10, 2) =     1.96092400d0 !        Neon -       Helium
      alpb(10, 3) =     2.19366600d0 !        Neon -      Lithium
      xfac(10, 3) =     0.70495800d0 !        Neon -      Lithium
      alpb(10, 4) =     1.31658800d0 !        Neon -    Beryllium
      xfac(10, 4) =     0.39262800d0 !        Neon -    Beryllium
      alpb(10, 5) =     2.75619000d0 !        Neon -        Boron
      xfac(10, 5) =     2.76414000d0 !        Neon -        Boron
      alpb(10, 6) =     3.44118800d0 !        Neon -       Carbon
      xfac(10, 6) =     5.46878000d0 !        Neon -       Carbon
      alpb(10, 7) =     4.42637000d0 !        Neon -     Nitrogen
      xfac(10, 7) =    29.99960900d0 !        Neon -     Nitrogen
      alpb(10, 8) =     2.88958700d0 !        Neon -       Oxygen
      xfac(10, 8) =     0.76389900d0 !        Neon -       Oxygen
      alpb(10, 9) =     3.67561100d0 !        Neon -     Fluorine
      xfac(10, 9) =     2.70675400d0 !        Neon -     Fluorine
      alpb(10,10) =     3.97456700d0 !        Neon -         Neon
      xfac(10,10) =     2.79483000d0 !        Neon -         Neon
 !
      alpb(11, 1) =     0.51527161d0 !      Sodium -     Hydrogen
      xfac(11, 1) =     0.37076659d0 !      Sodium -     Hydrogen
      alpb(11, 2) =     1.70302900d0 !      Sodium -       Helium
      xfac(11, 2) =     4.28251700d0 !      Sodium -       Helium
      alpb(11, 3) =     1.13682830d0 !      Sodium -      Lithium
      xfac(11, 3) =     3.87369093d0 !      Sodium -      Lithium
      alpb(11, 4) =     1.25548000d0 !      Sodium -    Beryllium
      xfac(11, 4) =     3.12162000d0 !      Sodium -    Beryllium
      alpb(11, 5) =     3.32368357d0 !      Sodium -        Boron
      xfac(11, 5) =     3.02021725d0 !      Sodium -        Boron
      alpb(11, 6) =     0.47556435d0 !      Sodium -       Carbon
      xfac(11, 6) =     0.18435042d0 !      Sodium -       Carbon
      alpb(11, 7) =     2.03351686d0 !      Sodium -     Nitrogen
      xfac(11, 7) =     7.12402246d0 !      Sodium -     Nitrogen
      alpb(11, 8) =     1.11150192d0 !      Sodium -       Oxygen
      xfac(11, 8) =     0.33426088d0 !      Sodium -       Oxygen
      alpb(11, 9) =     1.67842921d0 !      Sodium -     Fluorine
      xfac(11, 9) =     0.42914701d0 !      Sodium -     Fluorine
      alpb(11,10) =     1.60620354d0 !      Sodium -         Neon
      xfac(11,10) =     1.41555315d0 !      Sodium -         Neon
      alpb(11,11) =     0.19848982d0 !      Sodium -       Sodium
      xfac(11,11) =     0.36490548d0 !      Sodium -       Sodium
 !
      alpb(12, 1) =     1.95625469d0 !   Magnesium -     Hydrogen
      xfac(12, 1) =     2.62124083d0 !   Magnesium -     Hydrogen
      alpb(12, 2) =     2.21060300d0 !   Magnesium -       Helium
      xfac(12, 2) =     3.72585000d0 !   Magnesium -       Helium
      alpb(12, 3) =     1.18438000d0 !   Magnesium -      Lithium
      xfac(12, 3) =     2.49025000d0 !   Magnesium -      Lithium
      alpb(12, 4) =     1.55759100d0 !   Magnesium -    Beryllium
      xfac(12, 4) =     2.06639200d0 !   Magnesium -    Beryllium
      alpb(12, 5) =     2.52744100d0 !   Magnesium -        Boron
      xfac(12, 5) =     6.14670100d0 !   Magnesium -        Boron
      alpb(12, 6) =     2.78987067d0 !   Magnesium -       Carbon
      xfac(12, 6) =     7.96134202d0 !   Magnesium -       Carbon
      alpb(12, 7) =     2.04885974d0 !   Magnesium -     Nitrogen
      xfac(12, 7) =     1.79455363d0 !   Magnesium -     Nitrogen
      alpb(12, 8) =     2.74355805d0 !   Magnesium -       Oxygen
      xfac(12, 8) =     3.88906745d0 !   Magnesium -       Oxygen
      alpb(12, 9) =     3.97448921d0 !   Magnesium -     Fluorine
      xfac(12, 9) =    18.30661744d0 !   Magnesium -     Fluorine
      alpb(12,10) =     2.03167600d0 !   Magnesium -         Neon
      xfac(12,10) =     1.21485900d0 !   Magnesium -         Neon
      alpb(12,11) =     1.50677300d0 !   Magnesium -       Sodium
      xfac(12,11) =     8.67561900d0 !   Magnesium -       Sodium
      alpb(12,12) =     1.35402211d0 !   Magnesium -    Magnesium
      xfac(12,12) =     1.74715966d0 !   Magnesium -    Magnesium
 !
      alpb(13, 1) =     2.02599600d0 !    Aluminum -     Hydrogen
      xfac(13, 1) =     2.95837900d0 !    Aluminum -     Hydrogen
      alpb(13, 2) =     2.25583000d0 !    Aluminum -       Helium
      xfac(13, 2) =     2.70140000d0 !    Aluminum -       Helium
      alpb(13, 3) =     1.58159300d0 !    Aluminum -      Lithium
      xfac(13, 3) =     1.10681900d0 !    Aluminum -      Lithium
      alpb(13, 4) =     1.93823700d0 !    Aluminum -    Beryllium
      xfac(13, 4) =     5.03721400d0 !    Aluminum -    Beryllium
      alpb(13, 5) =     2.05956900d0 !    Aluminum -        Boron
      xfac(13, 5) =     2.74147900d0 !    Aluminum -        Boron
      alpb(13, 6) =     2.26744000d0 !    Aluminum -       Carbon
      xfac(13, 6) =     2.92805600d0 !    Aluminum -       Carbon
      alpb(13, 7) =     2.00975400d0 !    Aluminum -     Nitrogen
      xfac(13, 7) =     1.34520200d0 !    Aluminum -     Nitrogen
      alpb(13, 8) =     2.49866000d0 !    Aluminum -       Oxygen
      xfac(13, 8) =     2.13139600d0 !    Aluminum -       Oxygen
      alpb(13, 9) =     3.08425800d0 !    Aluminum -     Fluorine
      xfac(13, 9) =     1.97563500d0 !    Aluminum -     Fluorine
      alpb(13,10) =     2.44786900d0 !    Aluminum -         Neon
      xfac(13,10) =     1.70920000d0 !    Aluminum -         Neon
      alpb(13,11) =     1.20287100d0 !    Aluminum -       Sodium
      xfac(13,11) =     2.07184700d0 !    Aluminum -       Sodium
      alpb(13,12) =     1.97253000d0 !    Aluminum -    Magnesium
      xfac(13,12) =    13.47244300d0 !    Aluminum -    Magnesium
      alpb(13,13) =     1.38771400d0 !    Aluminum -     Aluminum
      xfac(13,13) =     2.13920000d0 !    Aluminum -     Aluminum
 !
      alpb(14, 1) =     1.89695000d0 !     Silicon -     Hydrogen
      xfac(14, 1) =     0.92419600d0 !     Silicon -     Hydrogen
      alpb(14, 2) =     2.04049800d0 !     Silicon -       Helium
      xfac(14, 2) =     1.85358300d0 !     Silicon -       Helium
      alpb(14, 3) =     1.78960900d0 !     Silicon -      Lithium
      xfac(14, 3) =     3.09079100d0 !     Silicon -      Lithium
      alpb(14, 4) =     1.26313200d0 !     Silicon -    Beryllium
      xfac(14, 4) =     0.62343300d0 !     Silicon -    Beryllium
      alpb(14, 5) =     1.98265300d0 !     Silicon -        Boron
      xfac(14, 5) =     1.02828700d0 !     Silicon -        Boron
      alpb(14, 6) =     1.98449800d0 !     Silicon -       Carbon
      xfac(14, 6) =     0.78574500d0 !     Silicon -       Carbon
      alpb(14, 7) =     1.81898800d0 !     Silicon -     Nitrogen
      xfac(14, 7) =     0.59297200d0 !     Silicon -     Nitrogen
      alpb(14, 8) =     1.92360000d0 !     Silicon -       Oxygen
      xfac(14, 8) =     0.75109500d0 !     Silicon -       Oxygen
      alpb(14, 9) =     2.13102800d0 !     Silicon -     Fluorine
      xfac(14, 9) =     0.54351600d0 !     Silicon -     Fluorine
      alpb(14,10) =     2.86778400d0 !     Silicon -         Neon
      xfac(14,10) =    14.37867600d0 !     Silicon -         Neon
      alpb(14,11) =     2.00761500d0 !     Silicon -       Sodium
      xfac(14,11) =     9.23764400d0 !     Silicon -       Sodium
      alpb(14,12) =     3.13974900d0 !     Silicon -    Magnesium
      xfac(14,12) =    29.99452000d0 !     Silicon -    Magnesium
      alpb(14,13) =     1.90000000d0 !     Silicon -     Aluminum
      xfac(14,13) =     2.00000000d0 !     Silicon -     Aluminum
      alpb(14,14) =     1.32900000d0 !     Silicon -      Silicon
      xfac(14,14) =     0.27347700d0 !     Silicon -      Silicon
 !
      alpb(15, 1) =     1.63770846d0 !  Phosphorus -     Hydrogen
      xfac(15, 1) =     0.81812089d0 !  Phosphorus -     Hydrogen
      alpb(15, 2) =     2.09315800d0 !  Phosphorus -       Helium
      xfac(15, 2) =     1.49021800d0 !  Phosphorus -       Helium
      alpb(15, 3) =     1.39454400d0 !  Phosphorus -      Lithium
      xfac(15, 3) =     1.12295000d0 !  Phosphorus -      Lithium
      alpb(15, 4) =     1.80007000d0 !  Phosphorus -    Beryllium
      xfac(15, 4) =     1.68483100d0 !  Phosphorus -    Beryllium
      alpb(15, 5) =     1.92316800d0 !  Phosphorus -        Boron
      xfac(15, 5) =     1.45088600d0 !  Phosphorus -        Boron
      alpb(15, 6) =     1.79751022d0 !  Phosphorus -       Carbon
      xfac(15, 6) =     0.76432343d0 !  Phosphorus -       Carbon
      alpb(15, 7) =     1.92145817d0 !  Phosphorus -     Nitrogen
      xfac(15, 7) =     0.70889722d0 !  Phosphorus -     Nitrogen
      alpb(15, 8) =     2.24757842d0 !  Phosphorus -       Oxygen
      xfac(15, 8) =     0.81247514d0 !  Phosphorus -       Oxygen
      alpb(15, 9) =     2.24534369d0 !  Phosphorus -     Fluorine
      xfac(15, 9) =     0.54703350d0 !  Phosphorus -     Fluorine
      alpb(15,10) =     2.21903600d0 !  Phosphorus -         Neon
      xfac(15,10) =     0.77495400d0 !  Phosphorus -         Neon
      alpb(15,11) =     1.50032000d0 !  Phosphorus -       Sodium
      xfac(15,11) =     2.83709500d0 !  Phosphorus -       Sodium
      alpb(15,12) =     1.22015396d0 !  Phosphorus -    Magnesium
      xfac(15,12) =     1.23302173d0 !  Phosphorus -    Magnesium
      alpb(15,13) =     1.98072700d0 !  Phosphorus -     Aluminum
      xfac(15,13) =     5.05081600d0 !  Phosphorus -     Aluminum
      alpb(15,14) =     3.31346600d0 !  Phosphorus -      Silicon
      xfac(15,14) =    13.23912100d0 !  Phosphorus -      Silicon
      alpb(15,15) =     1.22886287d0 !  Phosphorus -   Phosphorus
      xfac(15,15) =     0.54734266d0 !  Phosphorus -   Phosphorus
 !
      alpb(16, 1) =     1.79726324d0 !      Sulfur -     Hydrogen
      xfac(16, 1) =     0.61549662d0 !      Sulfur -     Hydrogen
      alpb(16, 2) =     1.95914900d0 !      Sulfur -       Helium
      xfac(16, 2) =     0.43761800d0 !      Sulfur -       Helium
      alpb(16, 3) =     2.29427500d0 !      Sulfur -      Lithium
      xfac(16, 3) =     2.64250200d0 !      Sulfur -      Lithium
      alpb(16, 4) =     2.78173600d0 !      Sulfur -    Beryllium
      xfac(16, 4) =     3.79156500d0 !      Sulfur -    Beryllium
      alpb(16, 5) =     2.40369600d0 !      Sulfur -        Boron
      xfac(16, 5) =     1.12539400d0 !      Sulfur -        Boron
      alpb(16, 6) =     1.95864281d0 !      Sulfur -       Carbon
      xfac(16, 6) =     0.60266499d0 !      Sulfur -       Carbon
      alpb(16, 7) =     1.97016067d0 !      Sulfur -     Nitrogen
      xfac(16, 7) =     0.55284773d0 !      Sulfur -     Nitrogen
      alpb(16, 8) =     2.24141697d0 !      Sulfur -       Oxygen
      xfac(16, 8) =     0.63787106d0 !      Sulfur -       Oxygen
      alpb(16, 9) =     1.90325343d0 !      Sulfur -     Fluorine
      xfac(16, 9) =     0.33697482d0 !      Sulfur -     Fluorine
      alpb(16,10) =     2.78705800d0 !      Sulfur -         Neon
      xfac(16,10) =     3.29616000d0 !      Sulfur -         Neon
      alpb(16,11) =     0.57939050d0 !      Sulfur -       Sodium
      xfac(16,11) =     0.25066746d0 !      Sulfur -       Sodium
      alpb(16,12) =     2.17697792d0 !      Sulfur -    Magnesium
      xfac(16,12) =     3.78289700d0 !      Sulfur -    Magnesium
      alpb(16,13) =     1.97670500d0 !      Sulfur -     Aluminum
      xfac(16,13) =     2.34738400d0 !      Sulfur -     Aluminum
      alpb(16,14) =     1.88591600d0 !      Sulfur -      Silicon
      xfac(16,14) =     0.87665800d0 !      Sulfur -      Silicon
      alpb(16,15) =     1.56670499d0 !      Sulfur -   Phosphorus
      xfac(16,15) =     0.66196577d0 !      Sulfur -   Phosphorus
      alpb(16,16) =     1.67096453d0 !      Sulfur -       Sulfur
      xfac(16,16) =     0.59310014d0 !      Sulfur -       Sulfur
 !
      alpb(17, 1) =     1.70789841d0 !    Chlorine -     Hydrogen
      xfac(17, 1) =     0.36046457d0 !    Chlorine -     Hydrogen
      alpb(17, 2) =     1.67167700d0 !    Chlorine -       Helium
      xfac(17, 2) =     0.27296400d0 !    Chlorine -       Helium
      alpb(17, 3) =     2.78300100d0 !    Chlorine -      Lithium
      xfac(17, 3) =     4.22779400d0 !    Chlorine -      Lithium
      alpb(17, 4) =     2.82267600d0 !    Chlorine -    Beryllium
      xfac(17, 4) =     2.50727500d0 !    Chlorine -    Beryllium
      alpb(17, 5) =     2.25932300d0 !    Chlorine -        Boron
      xfac(17, 5) =     0.82212900d0 !    Chlorine -        Boron
      alpb(17, 6) =     1.69392115d0 !    Chlorine -       Carbon
      xfac(17, 6) =     0.29258663d0 !    Chlorine -       Carbon
      alpb(17, 7) =     1.73250222d0 !    Chlorine -     Nitrogen
      xfac(17, 7) =     0.30394869d0 !    Chlorine -     Nitrogen
      alpb(17, 8) =     1.74133731d0 !    Chlorine -       Oxygen
      xfac(17, 8) =     0.25638696d0 !    Chlorine -       Oxygen
      alpb(17, 9) =     1.66575401d0 !    Chlorine -     Fluorine
      xfac(17, 9) =     0.18656760d0 !    Chlorine -     Fluorine
      alpb(17,10) =     1.70315100d0 !    Chlorine -         Neon
      xfac(17,10) =     0.12513300d0 !    Chlorine -         Neon
      alpb(17,11) =     1.81642900d0 !    Chlorine -       Sodium
      xfac(17,11) =     1.35789400d0 !    Chlorine -       Sodium
      alpb(17,12) =     2.32537421d0 !    Chlorine -    Magnesium
      xfac(17,12) =     3.15279932d0 !    Chlorine -    Magnesium
      alpb(17,13) =     2.12593900d0 !    Chlorine -     Aluminum
      xfac(17,13) =     2.15345100d0 !    Chlorine -     Aluminum
      alpb(17,14) =     1.68497800d0 !    Chlorine -      Silicon
      xfac(17,14) =     0.51300000d0 !    Chlorine -      Silicon
      alpb(17,15) =     1.41706691d0 !    Chlorine -   Phosphorus
      xfac(17,15) =     0.37570275d0 !    Chlorine -   Phosphorus
      alpb(17,16) =     1.55645719d0 !    Chlorine -       Sulfur
      xfac(17,16) =     0.39562708d0 !    Chlorine -       Sulfur
      alpb(17,17) =     1.60269643d0 !    Chlorine -     Chlorine
      xfac(17,17) =     0.37280371d0 !    Chlorine -     Chlorine
 !
      alpb(18, 1) =     4.05616700d0 !       Argon -     Hydrogen
      xfac(18, 1) =     3.93344500d0 !       Argon -     Hydrogen
      alpb(18, 2) =     2.71656200d0 !       Argon -       Helium
      xfac(18, 2) =     1.17721100d0 !       Argon -       Helium
      alpb(18, 3) =     3.12289500d0 !       Argon -      Lithium
      xfac(18, 3) =     3.36291000d0 !       Argon -      Lithium
      alpb(18, 4) =     3.04400700d0 !       Argon -    Beryllium
      xfac(18, 4) =     2.75549200d0 !       Argon -    Beryllium
      alpb(18, 5) =     2.41547100d0 !       Argon -        Boron
      xfac(18, 5) =     1.93158600d0 !       Argon -        Boron
      alpb(18, 6) =     1.47130900d0 !       Argon -       Carbon
      xfac(18, 6) =     0.12230900d0 !       Argon -       Carbon
      alpb(18, 7) =     2.32680500d0 !       Argon -     Nitrogen
      xfac(18, 7) =     0.56258100d0 !       Argon -     Nitrogen
      alpb(18, 8) =     2.24067300d0 !       Argon -       Oxygen
      xfac(18, 8) =     0.35579500d0 !       Argon -       Oxygen
      alpb(18, 9) =     3.92065800d0 !       Argon -     Fluorine
      xfac(18, 9) =     9.26971500d0 !       Argon -     Fluorine
      alpb(18,10) =     2.96374700d0 !       Argon -         Neon
      xfac(18,10) =     1.30469700d0 !       Argon -         Neon
      alpb(18,11) =     2.16767700d0 !       Argon -       Sodium
      xfac(18,11) =     3.39813800d0 !       Argon -       Sodium
      alpb(18,12) =     2.09266400d0 !       Argon -    Magnesium
      xfac(18,12) =     1.97063800d0 !       Argon -    Magnesium
      alpb(18,13) =     2.64516500d0 !       Argon -     Aluminum
      xfac(18,13) =     1.85200900d0 !       Argon -     Aluminum
      alpb(18,14) =     1.78035000d0 !       Argon -      Silicon
      xfac(18,14) =     1.06789000d0 !       Argon -      Silicon
      alpb(18,15) =     4.37251600d0 !       Argon -   Phosphorus
      xfac(18,15) =     0.17101400d0 !       Argon -   Phosphorus
      alpb(18,16) =     2.04939800d0 !       Argon -       Sulfur
      xfac(18,16) =     0.65376900d0 !       Argon -       Sulfur
      alpb(18,17) =     2.55444900d0 !       Argon -     Chlorine
      xfac(18,17) =     2.25609400d0 !       Argon -     Chlorine
      alpb(18,18) =     2.30643200d0 !       Argon -        Argon
      xfac(18,18) =     0.97269900d0 !       Argon -        Argon
 !
      alpb(19, 1) =     1.26382808d0 !   Potassium -     Hydrogen
      xfac(19, 1) =     3.92136315d0 !   Potassium -     Hydrogen
      alpb(19, 2) =     1.41850100d0 !   Potassium -       Helium
      xfac(19, 2) =     2.89504500d0 !   Potassium -       Helium
      alpb(19, 3) =     1.03648700d0 !   Potassium -      Lithium
      xfac(19, 3) =     4.37456700d0 !   Potassium -      Lithium
      alpb(19, 4) =     1.93188800d0 !   Potassium -    Beryllium
      xfac(19, 4) =     6.73222100d0 !   Potassium -    Beryllium
      alpb(19, 5) =     2.03176800d0 !   Potassium -        Boron
      xfac(19, 5) =     8.90054100d0 !   Potassium -        Boron
      alpb(19, 6) =     1.50263562d0 !   Potassium -       Carbon
      xfac(19, 6) =     2.61589246d0 !   Potassium -       Carbon
      alpb(19, 7) =     1.78547490d0 !   Potassium -     Nitrogen
      xfac(19, 7) =     2.96745329d0 !   Potassium -     Nitrogen
      alpb(19, 8) =     1.02003042d0 !   Potassium -       Oxygen
      xfac(19, 8) =     0.28228579d0 !   Potassium -       Oxygen
      alpb(19, 9) =     3.51882956d0 !   Potassium -     Fluorine
      xfac(19, 9) =     6.83779848d0 !   Potassium -     Fluorine
      alpb(19,10) =     1.13802100d0 !   Potassium -         Neon
      xfac(19,10) =     0.23399500d0 !   Potassium -         Neon
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
      alpb(19,16) =     2.20715908d0 !   Potassium -       Sulfur
      xfac(19,16) =    29.40272366d0 !   Potassium -       Sulfur
      alpb(19,17) =     1.07718412d0 !   Potassium -     Chlorine
      xfac(19,17) =     0.66754058d0 !   Potassium -     Chlorine
      alpb(19,18) =     2.30280300d0 !   Potassium -        Argon
      xfac(19,18) =     9.71050800d0 !   Potassium -        Argon
      alpb(19,19) =     1.62207850d0 !   Potassium -    Potassium
      xfac(19,19) =     5.77267253d0 !   Potassium -    Potassium
 !
      alpb(20, 1) =     1.75280570d0 !     Calcium -     Hydrogen
      xfac(20, 1) =     4.44935744d0 !     Calcium -     Hydrogen
      alpb(20, 2) =     1.71984700d0 !     Calcium -       Helium
      xfac(20, 2) =     2.91385200d0 !     Calcium -       Helium
      alpb(20, 5) =     1.70001000d0 !     Calcium -        Boron
      xfac(20, 5) =     1.70001000d0 !     Calcium -        Boron
      alpb(20, 6) =     1.35726415d0 !     Calcium -       Carbon
      xfac(20, 6) =     0.62878644d0 !     Calcium -       Carbon
      alpb(20, 7) =     1.66602298d0 !     Calcium -     Nitrogen
      xfac(20, 7) =     1.17109711d0 !     Calcium -     Nitrogen
      alpb(20, 8) =     1.71094766d0 !     Calcium -       Oxygen
      xfac(20, 8) =     0.67980812d0 !     Calcium -       Oxygen
      alpb(20, 9) =     3.52091284d0 !     Calcium -     Fluorine
      xfac(20, 9) =    10.98966965d0 !     Calcium -     Fluorine
      alpb(20,10) =     0.95453000d0 !     Calcium -         Neon
      xfac(20,10) =     0.33258600d0 !     Calcium -         Neon
      alpb(20,11) =     3.10710400d0 !     Calcium -       Sodium
      xfac(20,11) =     9.65750900d0 !     Calcium -       Sodium
      alpb(20,12) =     2.29980000d0 !     Calcium -    Magnesium
      xfac(20,12) =     8.59980000d0 !     Calcium -    Magnesium
      alpb(20,13) =     1.61256500d0 !     Calcium -     Aluminum
      xfac(20,13) =     4.18855500d0 !     Calcium -     Aluminum
      alpb(20,14) =     1.21878800d0 !     Calcium -      Silicon
      xfac(20,14) =     0.33623300d0 !     Calcium -      Silicon
      alpb(20,15) =     1.02414200d0 !     Calcium -   Phosphorus
      xfac(20,15) =     0.41084000d0 !     Calcium -   Phosphorus
      alpb(20,16) =     0.75537722d0 !     Calcium -       Sulfur
      xfac(20,16) =     0.27652269d0 !     Calcium -       Sulfur
      alpb(20,17) =     1.03798455d0 !     Calcium -     Chlorine
      xfac(20,17) =     0.31117199d0 !     Calcium -     Chlorine
      alpb(20,18) =     1.03488100d0 !     Calcium -        Argon
      xfac(20,18) =     0.29107200d0 !     Calcium -        Argon
      alpb(20,19) =     1.11920000d0 !     Calcium -    Potassium
      xfac(20,19) =     1.24032000d0 !     Calcium -    Potassium
      alpb(20,20) =     1.16372040d0 !     Calcium -      Calcium
      xfac(20,20) =    29.49299918d0 !     Calcium -      Calcium
 !
      alpb(21, 1) =     1.17948500d0 !    Scandium -     Hydrogen
      xfac(21, 1) =     0.35119900d0 !    Scandium -     Hydrogen
      alpb(21, 6) =     2.63049000d0 !    Scandium -       Carbon
      xfac(21, 6) =     8.60805200d0 !    Scandium -       Carbon
      alpb(21, 7) =     2.27000400d0 !    Scandium -     Nitrogen
      xfac(21, 7) =     3.23188100d0 !    Scandium -     Nitrogen
      alpb(21, 8) =     2.25651600d0 !    Scandium -       Oxygen
      xfac(21, 8) =     3.05867200d0 !    Scandium -       Oxygen
      alpb(21, 9) =     3.10798500d0 !    Scandium -     Fluorine
      xfac(21, 9) =     7.25234700d0 !    Scandium -     Fluorine
      alpb(21,13) =     1.00355000d0 !    Scandium -     Aluminum
      xfac(21,13) =     0.50062000d0 !    Scandium -     Aluminum
      alpb(21,14) =     2.01687000d0 !    Scandium -      Silicon
      xfac(21,14) =     3.21907000d0 !    Scandium -      Silicon
      alpb(21,15) =     0.86816500d0 !    Scandium -   Phosphorus
      xfac(21,15) =     0.62674900d0 !    Scandium -   Phosphorus
      alpb(21,16) =     0.42293900d0 !    Scandium -       Sulfur
      xfac(21,16) =     0.21185000d0 !    Scandium -       Sulfur
      alpb(21,17) =     2.14147400d0 !    Scandium -     Chlorine
      xfac(21,17) =     2.99612900d0 !    Scandium -     Chlorine
      alpb(21,21) =     1.13283800d0 !    Scandium -     Scandium
      xfac(21,21) =     2.59816600d0 !    Scandium -     Scandium
 !
      alpb(22, 1) =     0.83266900d0 !    Titanium -     Hydrogen
      xfac(22, 1) =     0.14372200d0 !    Titanium -     Hydrogen
      alpb(22, 5) =     1.62871000d0 !    Titanium -        Boron
      xfac(22, 5) =     0.64936000d0 !    Titanium -        Boron
      alpb(22, 6) =     1.59797300d0 !    Titanium -       Carbon
      xfac(22, 6) =     0.41670600d0 !    Titanium -       Carbon
      alpb(22, 7) =     1.67868600d0 !    Titanium -     Nitrogen
      xfac(22, 7) =     0.54546100d0 !    Titanium -     Nitrogen
      alpb(22, 8) =     1.78911800d0 !    Titanium -       Oxygen
      xfac(22, 8) =     0.79948600d0 !    Titanium -       Oxygen
      alpb(22, 9) =     2.30708700d0 !    Titanium -     Fluorine
      xfac(22, 9) =     1.08574200d0 !    Titanium -     Fluorine
      alpb(22,12) =     1.91134000d0 !    Titanium -    Magnesium
      xfac(22,12) =     4.33024000d0 !    Titanium -    Magnesium
      alpb(22,13) =     1.36948600d0 !    Titanium -     Aluminum
      xfac(22,13) =     2.09184100d0 !    Titanium -     Aluminum
      alpb(22,14) =     2.85603800d0 !    Titanium -      Silicon
      xfac(22,14) =     6.77381500d0 !    Titanium -      Silicon
      alpb(22,15) =     2.15192900d0 !    Titanium -   Phosphorus
      xfac(22,15) =     4.15050000d0 !    Titanium -   Phosphorus
      alpb(22,16) =     1.84643900d0 !    Titanium -       Sulfur
      xfac(22,16) =     0.94378400d0 !    Titanium -       Sulfur
      alpb(22,17) =     1.46103400d0 !    Titanium -     Chlorine
      xfac(22,17) =     0.33329700d0 !    Titanium -     Chlorine
      alpb(22,20) =     2.00000000d0 !    Titanium -      Calcium
      xfac(22,20) =     4.10914100d0 !    Titanium -      Calcium
      alpb(22,22) =     2.64859700d0 !    Titanium -     Titanium
      xfac(22,22) =     2.00000000d0 !    Titanium -     Titanium
 !
      alpb(23, 1) =     1.28013300d0 !    Vanadium -     Hydrogen
      xfac(23, 1) =     0.10520400d0 !    Vanadium -     Hydrogen
      alpb(23, 6) =     2.78985500d0 !    Vanadium -       Carbon
      xfac(23, 6) =     1.93876000d0 !    Vanadium -       Carbon
      alpb(23, 7) =     1.60754000d0 !    Vanadium -     Nitrogen
      xfac(23, 7) =     0.27672500d0 !    Vanadium -     Nitrogen
      alpb(23, 8) =     1.62397300d0 !    Vanadium -       Oxygen
      xfac(23, 8) =     0.41531200d0 !    Vanadium -       Oxygen
      alpb(23, 9) =     1.82516000d0 !    Vanadium -     Fluorine
      xfac(23, 9) =     0.34281500d0 !    Vanadium -     Fluorine
      alpb(23,11) =     2.55101000d0 !    Vanadium -       Sodium
      xfac(23,11) =     8.27602000d0 !    Vanadium -       Sodium
      alpb(23,15) =     2.54915400d0 !    Vanadium -   Phosphorus
      xfac(23,15) =     6.25062400d0 !    Vanadium -   Phosphorus
      alpb(23,16) =     2.70412400d0 !    Vanadium -       Sulfur
      xfac(23,16) =     2.03503900d0 !    Vanadium -       Sulfur
      alpb(23,17) =     1.68852900d0 !    Vanadium -     Chlorine
      xfac(23,17) =     0.24365700d0 !    Vanadium -     Chlorine
      alpb(23,19) =     4.52136000d0 !    Vanadium -    Potassium
      xfac(23,19) =     2.02659000d0 !    Vanadium -    Potassium
      alpb(23,23) =     4.83239100d0 !    Vanadium -     Vanadium
      xfac(23,23) =    10.77989200d0 !    Vanadium -     Vanadium
 !
      alpb(24, 1) =     0.88266100d0 !    Chromium -     Hydrogen
      xfac(24, 1) =     0.04446900d0 !    Chromium -     Hydrogen
      alpb(24, 6) =     3.65675400d0 !    Chromium -       Carbon
      xfac(24, 6) =     6.11018700d0 !    Chromium -       Carbon
      alpb(24, 7) =     3.02918600d0 !    Chromium -     Nitrogen
      xfac(24, 7) =     1.92032400d0 !    Chromium -     Nitrogen
      alpb(24, 8) =     2.50000000d0 !    Chromium -       Oxygen
      xfac(24, 8) =     1.05551100d0 !    Chromium -       Oxygen
      alpb(24, 9) =     2.71652100d0 !    Chromium -     Fluorine
      xfac(24, 9) =     0.73760700d0 !    Chromium -     Fluorine
      alpb(24,11) =     2.29505600d0 !    Chromium -       Sodium
      xfac(24,11) =     8.36427400d0 !    Chromium -       Sodium
      alpb(24,14) =     1.86076000d0 !    Chromium -      Silicon
      xfac(24,14) =     1.02911000d0 !    Chromium -      Silicon
      alpb(24,15) =     1.69538300d0 !    Chromium -   Phosphorus
      xfac(24,15) =     0.60017700d0 !    Chromium -   Phosphorus
      alpb(24,16) =     2.26097800d0 !    Chromium -       Sulfur
      xfac(24,16) =     0.55033400d0 !    Chromium -       Sulfur
      alpb(24,17) =     2.15261800d0 !    Chromium -     Chlorine
      xfac(24,17) =     0.36907300d0 !    Chromium -     Chlorine
      alpb(24,19) =     2.00000000d0 !    Chromium -    Potassium
      xfac(24,19) =     2.00000000d0 !    Chromium -    Potassium
      alpb(24,24) =     4.65541900d0 !    Chromium -     Chromium
      xfac(24,24) =    10.31860700d0 !    Chromium -     Chromium
 !
      alpb(25, 1) =     2.30994000d0 !   Manganese -     Hydrogen
      xfac(25, 1) =     1.26921000d0 !   Manganese -     Hydrogen
      alpb(25, 6) =     3.00075000d0 !   Manganese -       Carbon
      xfac(25, 6) =     2.58311000d0 !   Manganese -       Carbon
      alpb(25, 7) =     2.92147000d0 !   Manganese -     Nitrogen
      xfac(25, 7) =     1.95675000d0 !   Manganese -     Nitrogen
      alpb(25, 8) =     2.57754000d0 !   Manganese -       Oxygen
      xfac(25, 8) =     1.28562000d0 !   Manganese -       Oxygen
      alpb(25, 9) =     2.79195000d0 !   Manganese -     Fluorine
      xfac(25, 9) =     1.11307000d0 !   Manganese -     Fluorine
      alpb(25,13) =     1.76836000d0 !   Manganese -     Aluminum
      xfac(25,13) =     1.04079000d0 !   Manganese -     Aluminum
      alpb(25,14) =     1.93795900d0 !   Manganese -      Silicon
      xfac(25,14) =     0.95058000d0 !   Manganese -      Silicon
      alpb(25,15) =     1.94702000d0 !   Manganese -   Phosphorus
      xfac(25,15) =     1.13032000d0 !   Manganese -   Phosphorus
      alpb(25,16) =     2.48251000d0 !   Manganese -       Sulfur
      xfac(25,16) =     1.61265000d0 !   Manganese -       Sulfur
      alpb(25,17) =     1.65701000d0 !   Manganese -     Chlorine
      xfac(25,17) =     0.20185000d0 !   Manganese -     Chlorine
      alpb(25,20) =     1.49144000d0 !   Manganese -      Calcium
      xfac(25,20) =     0.62018000d0 !   Manganese -      Calcium
      alpb(25,25) =     2.66542000d0 !   Manganese -    Manganese
      xfac(25,25) =     2.46004000d0 !   Manganese -    Manganese
 !
      alpb(26, 1) =     1.40315740d0 !        Iron -     Hydrogen
      xfac(26, 1) =     0.29604424d0 !        Iron -     Hydrogen
      alpb(26, 6) =     3.14575186d0 !        Iron -       Carbon
      xfac(26, 6) =     2.62934167d0 !        Iron -       Carbon
      alpb(26, 7) =     3.37932899d0 !        Iron -     Nitrogen
      xfac(26, 7) =     2.09260038d0 !        Iron -     Nitrogen
      alpb(26, 8) =     3.98069771d0 !        Iron -       Oxygen
      xfac(26, 8) =     5.69691446d0 !        Iron -       Oxygen
      alpb(26, 9) =     3.35802534d0 !        Iron -     Fluorine
      xfac(26, 9) =     1.53997079d0 !        Iron -     Fluorine
      alpb(26,15) =     1.42597541d0 !        Iron -   Phosphorus
      xfac(26,15) =     0.59756052d0 !        Iron -   Phosphorus
      alpb(26,16) =     1.06819377d0 !        Iron -       Sulfur
      xfac(26,16) =     0.11563478d0 !        Iron -       Sulfur
      alpb(26,17) =     1.84765876d0 !        Iron -     Chlorine
      xfac(26,17) =     0.21314306d0 !        Iron -     Chlorine
      alpb(26,19) =     2.00000000d0 !        Iron -    Potassium
      xfac(26,19) =     6.00000000d0 !        Iron -    Potassium
      alpb(26,26) =     4.84293092d0 !        Iron -         Iron
      xfac(26,26) =    24.87954318d0 !        Iron -         Iron
 !
      alpb(27, 1) =     2.17108831d0 !      Cobalt -     Hydrogen
      xfac(27, 1) =     0.76194642d0 !      Cobalt -     Hydrogen
      alpb(27, 5) =     3.20000000d0 !      Cobalt -        Boron
      xfac(27, 5) =     1.00000000d0 !      Cobalt -        Boron
      alpb(27, 6) =     3.58818631d0 !      Cobalt -       Carbon
      xfac(27, 6) =     4.27873559d0 !      Cobalt -       Carbon
      alpb(27, 7) =     3.09432046d0 !      Cobalt -     Nitrogen
      xfac(27, 7) =     1.96214158d0 !      Cobalt -     Nitrogen
      alpb(27, 8) =     3.21946846d0 !      Cobalt -       Oxygen
      xfac(27, 8) =     2.17460434d0 !      Cobalt -       Oxygen
      alpb(27, 9) =     3.55245306d0 !      Cobalt -     Fluorine
      xfac(27, 9) =     3.58467454d0 !      Cobalt -     Fluorine
      alpb(27,14) =     2.46980500d0 !      Cobalt -      Silicon
      xfac(27,14) =     1.09024000d0 !      Cobalt -      Silicon
      alpb(27,15) =     1.15791919d0 !      Cobalt -   Phosphorus
      xfac(27,15) =     0.05332018d0 !      Cobalt -   Phosphorus
      alpb(27,16) =     2.33607740d0 !      Cobalt -       Sulfur
      xfac(27,16) =     1.34289480d0 !      Cobalt -       Sulfur
      alpb(27,17) =     1.98024559d0 !      Cobalt -     Chlorine
      xfac(27,17) =     0.18647888d0 !      Cobalt -     Chlorine
      alpb(27,27) =     3.44775733d0 !      Cobalt -       Cobalt
      xfac(27,27) =     3.99903614d0 !      Cobalt -       Cobalt
 !
      alpb(28, 1) =     2.63528000d0 !      Nickel -     Hydrogen
      xfac(28, 1) =     1.76312400d0 !      Nickel -     Hydrogen
      alpb(28, 6) =     4.28551300d0 !      Nickel -       Carbon
      xfac(28, 6) =     7.13332400d0 !      Nickel -       Carbon
      alpb(28, 7) =     3.84521500d0 !      Nickel -     Nitrogen
      xfac(28, 7) =     4.28680000d0 !      Nickel -     Nitrogen
      alpb(28, 8) =     2.93723200d0 !      Nickel -       Oxygen
      xfac(28, 8) =     0.88594200d0 !      Nickel -       Oxygen
      alpb(28, 9) =     3.44024100d0 !      Nickel -     Fluorine
      xfac(28, 9) =     1.08820800d0 !      Nickel -     Fluorine
      alpb(28,14) =     2.06888100d0 !      Nickel -      Silicon
      xfac(28,14) =     0.93864600d0 !      Nickel -      Silicon
      alpb(28,15) =     3.26028300d0 !      Nickel -   Phosphorus
      xfac(28,15) =     5.05972700d0 !      Nickel -   Phosphorus
      alpb(28,16) =     2.00275200d0 !      Nickel -       Sulfur
      xfac(28,16) =     0.27485200d0 !      Nickel -       Sulfur
      alpb(28,17) =     2.20051200d0 !      Nickel -     Chlorine
      xfac(28,17) =     0.20231300d0 !      Nickel -     Chlorine
      alpb(28,28) =     1.09796000d0 !      Nickel -       Nickel
      xfac(28,28) =     0.03547400d0 !      Nickel -       Nickel
 !
      alpb(29, 1) =     2.33535900d0 !      Copper -     Hydrogen
      xfac(29, 1) =     0.60359100d0 !      Copper -     Hydrogen
      alpb(29, 6) =     4.63877300d0 !      Copper -       Carbon
      xfac(29, 6) =     7.06779400d0 !      Copper -       Carbon
      alpb(29, 7) =     4.21433700d0 !      Copper -     Nitrogen
      xfac(29, 7) =     3.22866700d0 !      Copper -     Nitrogen
      alpb(29, 8) =     3.95995100d0 !      Copper -       Oxygen
      xfac(29, 8) =     2.00000000d0 !      Copper -       Oxygen
      alpb(29, 9) =     4.47883200d0 !      Copper -     Fluorine
      xfac(29, 9) =     1.28210800d0 !      Copper -     Fluorine
      alpb(29,15) =     0.21064000d0 !      Copper -   Phosphorus
      xfac(29,15) =     0.02012600d0 !      Copper -   Phosphorus
      alpb(29,16) =     0.27311200d0 !      Copper -       Sulfur
      xfac(29,16) =     0.00524800d0 !      Copper -       Sulfur
      alpb(29,17) =     2.77653100d0 !      Copper -     Chlorine
      xfac(29,17) =     0.13906500d0 !      Copper -     Chlorine
      alpb(29,29) =     3.61684600d0 !      Copper -       Copper
      xfac(29,29) =     5.18437600d0 !      Copper -       Copper
 !
      alpb(30, 1) =     1.91652557d0 !        Zinc -     Hydrogen
      xfac(30, 1) =     2.67852675d0 !        Zinc -     Hydrogen
      alpb(30, 6) =     1.91622954d0 !        Zinc -       Carbon
      xfac(30, 6) =     1.45805800d0 !        Zinc -       Carbon
      alpb(30, 7) =     1.67976867d0 !        Zinc -     Nitrogen
      xfac(30, 7) =     0.78960720d0 !        Zinc -     Nitrogen
      alpb(30, 8) =     2.19645641d0 !        Zinc -       Oxygen
      xfac(30, 8) =     1.49579306d0 !        Zinc -       Oxygen
      alpb(30, 9) =     2.41002100d0 !        Zinc -     Fluorine
      xfac(30, 9) =     1.22554500d0 !        Zinc -     Fluorine
      alpb(30,14) =     1.41074308d0 !        Zinc -      Silicon
      xfac(30,14) =     2.81658103d0 !        Zinc -      Silicon
      alpb(30,15) =     1.22048000d0 !        Zinc -   Phosphorus
      xfac(30,15) =     0.58153000d0 !        Zinc -   Phosphorus
      alpb(30,16) =     1.46484141d0 !        Zinc -       Sulfur
      xfac(30,16) =     0.82618049d0 !        Zinc -       Sulfur
      alpb(30,17) =     1.57853485d0 !        Zinc -     Chlorine
      xfac(30,17) =     0.71210113d0 !        Zinc -     Chlorine
      alpb(30,20) =     1.11918000d0 !        Zinc -      Calcium
      xfac(30,20) =     1.24029000d0 !        Zinc -      Calcium
      alpb(30,30) =     0.73736738d0 !        Zinc -         Zinc
      xfac(30,30) =     0.19303833d0 !        Zinc -         Zinc
 !
      alpb(31, 1) =     1.84735000d0 !     Gallium -     Hydrogen
      xfac(31, 1) =     1.38665200d0 !     Gallium -     Hydrogen
      alpb(31, 6) =     2.32541000d0 !     Gallium -       Carbon
      xfac(31, 6) =     1.96299000d0 !     Gallium -       Carbon
      alpb(31, 7) =     2.12182000d0 !     Gallium -     Nitrogen
      xfac(31, 7) =     1.18833800d0 !     Gallium -     Nitrogen
      alpb(31, 8) =     2.34834700d0 !     Gallium -       Oxygen
      xfac(31, 8) =     1.52364400d0 !     Gallium -       Oxygen
      alpb(31, 9) =     2.67986900d0 !     Gallium -     Fluorine
      xfac(31, 9) =     1.41694200d0 !     Gallium -     Fluorine
      alpb(31,14) =     1.91378000d0 !     Gallium -      Silicon
      xfac(31,14) =     1.00229000d0 !     Gallium -      Silicon
      alpb(31,15) =     2.97965000d0 !     Gallium -   Phosphorus
      xfac(31,15) =     0.50000000d0 !     Gallium -   Phosphorus
      alpb(31,16) =     2.23210800d0 !     Gallium -       Sulfur
      xfac(31,16) =     2.45628400d0 !     Gallium -       Sulfur
      alpb(31,17) =     2.02471000d0 !     Gallium -     Chlorine
      xfac(31,17) =     1.18666100d0 !     Gallium -     Chlorine
      alpb(31,31) =     1.33464300d0 !     Gallium -      Gallium
      xfac(31,31) =     1.19839400d0 !     Gallium -      Gallium
 !
      alpb(32, 1) =     2.20679300d0 !   Germanium -     Hydrogen
      xfac(32, 1) =     1.73322600d0 !   Germanium -     Hydrogen
      alpb(32, 6) =     2.25746900d0 !   Germanium -       Carbon
      xfac(32, 6) =     1.29751000d0 !   Germanium -       Carbon
      alpb(32, 7) =     1.98822600d0 !   Germanium -     Nitrogen
      xfac(32, 7) =     0.63750600d0 !   Germanium -     Nitrogen
      alpb(32, 8) =     2.13941300d0 !   Germanium -       Oxygen
      xfac(32, 8) =     0.82696400d0 !   Germanium -       Oxygen
      alpb(32, 9) =     2.38477700d0 !   Germanium -     Fluorine
      xfac(32, 9) =     0.65197700d0 !   Germanium -     Fluorine
      alpb(32,14) =     0.29972100d0 !   Germanium -      Silicon
      xfac(32,14) =     0.17868000d0 !   Germanium -      Silicon
      alpb(32,15) =     2.46929100d0 !   Germanium -   Phosphorus
      xfac(32,15) =     5.61634900d0 !   Germanium -   Phosphorus
      alpb(32,16) =     2.02458800d0 !   Germanium -       Sulfur
      xfac(32,16) =     1.16095700d0 !   Germanium -       Sulfur
      alpb(32,17) =     1.77122800d0 !   Germanium -     Chlorine
      xfac(32,17) =     0.54523900d0 !   Germanium -     Chlorine
      alpb(32,25) =     2.38283400d0 !   Germanium -    Manganese
      xfac(32,25) =     2.25515100d0 !   Germanium -    Manganese
      alpb(32,27) =     2.85261000d0 !   Germanium -       Cobalt
      xfac(32,27) =     2.15185000d0 !   Germanium -       Cobalt
      alpb(32,32) =     2.01900000d0 !   Germanium -    Germanium
      xfac(32,32) =     3.02300000d0 !   Germanium -    Germanium
 !
      alpb(33, 1) =     1.99352700d0 !     Arsenic -     Hydrogen
      xfac(33, 1) =     1.09058900d0 !     Arsenic -     Hydrogen
      alpb(33, 6) =     1.85506900d0 !     Arsenic -       Carbon
      xfac(33, 6) =     0.57909800d0 !     Arsenic -       Carbon
      alpb(33, 7) =     1.49654300d0 !     Arsenic -     Nitrogen
      xfac(33, 7) =     0.27333700d0 !     Arsenic -     Nitrogen
      alpb(33, 8) =     2.00395000d0 !     Arsenic -       Oxygen
      xfac(33, 8) =     0.70161400d0 !     Arsenic -       Oxygen
      alpb(33, 9) =     2.01258300d0 !     Arsenic -     Fluorine
      xfac(33, 9) =     0.40262800d0 !     Arsenic -     Fluorine
      alpb(33,13) =     1.15278600d0 !     Arsenic -     Aluminum
      xfac(33,13) =     1.00358000d0 !     Arsenic -     Aluminum
      alpb(33,14) =     1.91560000d0 !     Arsenic -      Silicon
      xfac(33,14) =     1.43070600d0 !     Arsenic -      Silicon
      alpb(33,16) =     1.95436800d0 !     Arsenic -       Sulfur
      xfac(33,16) =     1.03378400d0 !     Arsenic -       Sulfur
      alpb(33,17) =     1.69107000d0 !     Arsenic -     Chlorine
      xfac(33,17) =     0.45443300d0 !     Arsenic -     Chlorine
      alpb(33,22) =     1.93291100d0 !     Arsenic -     Titanium
      xfac(33,22) =     1.58131700d0 !     Arsenic -     Titanium
      alpb(33,27) =     3.36814000d0 !     Arsenic -       Cobalt
      xfac(33,27) =     1.67524000d0 !     Arsenic -       Cobalt
      alpb(33,30) =     1.45913000d0 !     Arsenic -         Zinc
      xfac(33,30) =     3.15657100d0 !     Arsenic -         Zinc
      alpb(33,31) =     1.73097700d0 !     Arsenic -      Gallium
      xfac(33,31) =     1.68629800d0 !     Arsenic -      Gallium
      alpb(33,33) =     1.58826400d0 !     Arsenic -      Arsenic
      xfac(33,33) =     0.73730700d0 !     Arsenic -      Arsenic
 !
      alpb(34, 1) =     1.91168436d0 !    Selenium -     Hydrogen
      xfac(34, 1) =     0.61962704d0 !    Selenium -     Hydrogen
      alpb(34, 6) =     1.99722384d0 !    Selenium -       Carbon
      xfac(34, 6) =     0.40774439d0 !    Selenium -       Carbon
      alpb(34, 7) =     1.42047963d0 !    Selenium -     Nitrogen
      xfac(34, 7) =     0.14955363d0 !    Selenium -     Nitrogen
      alpb(34, 8) =     2.55201064d0 !    Selenium -       Oxygen
      xfac(34, 8) =     0.70452922d0 !    Selenium -       Oxygen
      alpb(34, 9) =     1.81182997d0 !    Selenium -     Fluorine
      xfac(34, 9) =     0.18328066d0 !    Selenium -     Fluorine
      alpb(34,14) =     1.52981700d0 !    Selenium -      Silicon
      xfac(34,14) =     0.51822700d0 !    Selenium -      Silicon
      alpb(34,15) =     1.09337576d0 !    Selenium -   Phosphorus
      xfac(34,15) =     0.24826135d0 !    Selenium -   Phosphorus
      alpb(34,16) =     1.40169306d0 !    Selenium -       Sulfur
      xfac(34,16) =     0.45309356d0 !    Selenium -       Sulfur
      alpb(34,17) =     1.39667869d0 !    Selenium -     Chlorine
      xfac(34,17) =     0.19163082d0 !    Selenium -     Chlorine
      alpb(34,25) =     2.64803800d0 !    Selenium -    Manganese
      xfac(34,25) =     2.18072000d0 !    Selenium -    Manganese
      alpb(34,27) =     2.52345000d0 !    Selenium -       Cobalt
      xfac(34,27) =     2.20241000d0 !    Selenium -       Cobalt
      alpb(34,30) =     1.18241520d0 !    Selenium -         Zinc
      xfac(34,30) =     0.53118428d0 !    Selenium -         Zinc
      alpb(34,32) =     2.66905700d0 !    Selenium -    Germanium
      xfac(34,32) =     5.87205100d0 !    Selenium -    Germanium
      alpb(34,33) =     1.66528000d0 !    Selenium -      Arsenic
      xfac(34,33) =     0.71126100d0 !    Selenium -      Arsenic
      alpb(34,34) =     1.04604879d0 !    Selenium -     Selenium
      xfac(34,34) =     0.08287186d0 !    Selenium -     Selenium
 !
      alpb(35, 1) =     2.32330237d0 !     Bromine -     Hydrogen
      xfac(35, 1) =     1.16256718d0 !     Bromine -     Hydrogen
      alpb(35, 2) =     2.12827500d0 !     Bromine -       Helium
      xfac(35, 2) =     1.06204300d0 !     Bromine -       Helium
      alpb(35, 3) =     2.07444100d0 !     Bromine -      Lithium
      xfac(35, 3) =     1.85886600d0 !     Bromine -      Lithium
      alpb(35, 4) =     2.36714600d0 !     Bromine -    Beryllium
      xfac(35, 4) =     1.94093300d0 !     Bromine -    Beryllium
      alpb(35, 5) =     2.30789000d0 !     Bromine -        Boron
      xfac(35, 5) =     1.22642000d0 !     Bromine -        Boron
      alpb(35, 6) =     2.28910440d0 !     Bromine -       Carbon
      xfac(35, 6) =     1.00298419d0 !     Bromine -       Carbon
      alpb(35, 7) =     4.02853661d0 !     Bromine -     Nitrogen
      xfac(35, 7) =    27.79959282d0 !     Bromine -     Nitrogen
      alpb(35, 8) =     3.17944225d0 !     Bromine -       Oxygen
      xfac(35, 8) =     2.68460222d0 !     Bromine -       Oxygen
      alpb(35, 9) =     2.88334772d0 !     Bromine -     Fluorine
      xfac(35, 9) =     1.12242172d0 !     Bromine -     Fluorine
      alpb(35,10) =     2.46417200d0 !     Bromine -         Neon
      xfac(35,10) =     1.00615900d0 !     Bromine -         Neon
      alpb(35,11) =     1.62221800d0 !     Bromine -       Sodium
      xfac(35,11) =     1.75293700d0 !     Bromine -       Sodium
      alpb(35,12) =     2.12290506d0 !     Bromine -    Magnesium
      xfac(35,12) =     4.26099410d0 !     Bromine -    Magnesium
      alpb(35,13) =     1.89414100d0 !     Bromine -     Aluminum
      xfac(35,13) =     2.35713000d0 !     Bromine -     Aluminum
      alpb(35,14) =     1.57082500d0 !     Bromine -      Silicon
      xfac(35,14) =     0.58951100d0 !     Bromine -      Silicon
      alpb(35,15) =     1.61508016d0 !     Bromine -   Phosphorus
      xfac(35,15) =     0.73579343d0 !     Bromine -   Phosphorus
      alpb(35,16) =     2.27133539d0 !     Bromine -       Sulfur
      xfac(35,16) =     1.76382187d0 !     Bromine -       Sulfur
      alpb(35,17) =     1.56815388d0 !     Bromine -     Chlorine
      xfac(35,17) =     0.30249893d0 !     Bromine -     Chlorine
      alpb(35,18) =     2.45080100d0 !     Bromine -        Argon
      xfac(35,18) =     3.26266800d0 !     Bromine -        Argon
      alpb(35,19) =     1.42409156d0 !     Bromine -    Potassium
      xfac(35,19) =     4.68569482d0 !     Bromine -    Potassium
      alpb(35,20) =     1.31539226d0 !     Bromine -      Calcium
      xfac(35,20) =     1.03360637d0 !     Bromine -      Calcium
      alpb(35,21) =     1.79348600d0 !     Bromine -     Scandium
      xfac(35,21) =     2.09825100d0 !     Bromine -     Scandium
      alpb(35,22) =     1.67484700d0 !     Bromine -     Titanium
      xfac(35,22) =     0.88343400d0 !     Bromine -     Titanium
      alpb(35,23) =     1.90290400d0 !     Bromine -     Vanadium
      xfac(35,23) =     0.61269800d0 !     Bromine -     Vanadium
      alpb(35,24) =     1.56602800d0 !     Bromine -     Chromium
      xfac(35,24) =     0.21785300d0 !     Bromine -     Chromium
      alpb(35,25) =     2.28382000d0 !     Bromine -    Manganese
      xfac(35,25) =     1.18358000d0 !     Bromine -    Manganese
      alpb(35,26) =     2.58821794d0 !     Bromine -         Iron
      xfac(35,26) =     1.95368018d0 !     Bromine -         Iron
      alpb(35,27) =     0.51635973d0 !     Bromine -       Cobalt
      xfac(35,27) =     0.03830792d0 !     Bromine -       Cobalt
      alpb(35,28) =     2.77213600d0 !     Bromine -       Nickel
      xfac(35,28) =     0.63214500d0 !     Bromine -       Nickel
      alpb(35,29) =     5.82640700d0 !     Bromine -       Copper
      xfac(35,29) =     0.76851700d0 !     Bromine -       Copper
      alpb(35,30) =     1.42314736d0 !     Bromine -         Zinc
      xfac(35,30) =     0.93939095d0 !     Bromine -         Zinc
      alpb(35,31) =     1.81910500d0 !     Bromine -      Gallium
      xfac(35,31) =     1.26103600d0 !     Bromine -      Gallium
      alpb(35,32) =     1.60236600d0 !     Bromine -    Germanium
      xfac(35,32) =     0.62773700d0 !     Bromine -    Germanium
      alpb(35,33) =     1.52017000d0 !     Bromine -      Arsenic
      xfac(35,33) =     0.51415300d0 !     Bromine -      Arsenic
      alpb(35,34) =     2.16491691d0 !     Bromine -     Selenium
      xfac(35,34) =     2.15588763d0 !     Bromine -     Selenium
      alpb(35,35) =     2.17381413d0 !     Bromine -      Bromine
      xfac(35,35) =     1.21049397d0 !     Bromine -      Bromine
 !
      alpb(36, 1) =     3.77045300d0 !     Krypton -     Hydrogen
      xfac(36, 1) =     5.12589700d0 !     Krypton -     Hydrogen
      alpb(36, 2) =     1.99694300d0 !     Krypton -       Helium
      xfac(36, 2) =     0.62770100d0 !     Krypton -       Helium
      alpb(36, 3) =     3.31456200d0 !     Krypton -      Lithium
      xfac(36, 3) =     8.75869700d0 !     Krypton -      Lithium
      alpb(36, 4) =     3.25304800d0 !     Krypton -    Beryllium
      xfac(36, 4) =    10.23779600d0 !     Krypton -    Beryllium
      alpb(36, 5) =     2.36316900d0 !     Krypton -        Boron
      xfac(36, 5) =     2.94678100d0 !     Krypton -        Boron
      alpb(36, 6) =     2.07673800d0 !     Krypton -       Carbon
      xfac(36, 6) =     0.65262300d0 !     Krypton -       Carbon
      alpb(36, 7) =     1.64405200d0 !     Krypton -     Nitrogen
      xfac(36, 7) =     0.19960600d0 !     Krypton -     Nitrogen
      alpb(36, 8) =     0.29230000d0 !     Krypton -       Oxygen
      xfac(36, 8) =     0.00673300d0 !     Krypton -       Oxygen
      alpb(36, 9) =     3.45232100d0 !     Krypton -     Fluorine
      xfac(36, 9) =     4.13440700d0 !     Krypton -     Fluorine
      alpb(36,10) =     2.81367900d0 !     Krypton -         Neon
      xfac(36,10) =     1.43372200d0 !     Krypton -         Neon
      alpb(36,11) =     2.29428745d0 !     Krypton -       Sodium
      xfac(36,11) =     8.35794790d0 !     Krypton -       Sodium
      alpb(36,12) =     1.39148700d0 !     Krypton -    Magnesium
      xfac(36,12) =     0.88843600d0 !     Krypton -    Magnesium
      alpb(36,13) =     2.46713100d0 !     Krypton -     Aluminum
      xfac(36,13) =     5.09171600d0 !     Krypton -     Aluminum
      alpb(36,14) =     1.76410000d0 !     Krypton -      Silicon
      xfac(36,14) =     0.55425000d0 !     Krypton -      Silicon
      alpb(36,17) =     1.88497400d0 !     Krypton -     Chlorine
      xfac(36,17) =     0.52021700d0 !     Krypton -     Chlorine
      alpb(36,18) =     1.99512500d0 !     Krypton -        Argon
      xfac(36,18) =     0.55487400d0 !     Krypton -        Argon
      alpb(36,19) =     2.18248700d0 !     Krypton -    Potassium
      xfac(36,19) =     8.60978200d0 !     Krypton -    Potassium
      alpb(36,20) =     1.30519700d0 !     Krypton -      Calcium
      xfac(36,20) =     0.87889100d0 !     Krypton -      Calcium
      alpb(36,35) =     1.52900600d0 !     Krypton -      Bromine
      xfac(36,35) =     0.30809800d0 !     Krypton -      Bromine
      alpb(36,36) =     1.13531900d0 !     Krypton -      Krypton
      xfac(36,36) =     0.05209900d0 !     Krypton -      Krypton
 !
      alpb(37, 1) =     2.44355600d0 !    Rubidium -     Hydrogen
      xfac(37, 1) =    29.86163200d0 !    Rubidium -     Hydrogen
      alpb(37, 2) =     1.27074100d0 !    Rubidium -       Helium
      xfac(37, 2) =     1.86258500d0 !    Rubidium -       Helium
      alpb(37, 5) =     5.53223900d0 !    Rubidium -        Boron
      xfac(37, 5) =     9.04049300d0 !    Rubidium -        Boron
      alpb(37, 6) =     2.76583000d0 !    Rubidium -       Carbon
      xfac(37, 6) =    29.97403100d0 !    Rubidium -       Carbon
      alpb(37, 7) =     0.76104700d0 !    Rubidium -     Nitrogen
      xfac(37, 7) =     0.02463600d0 !    Rubidium -     Nitrogen
      alpb(37, 8) =     1.33490800d0 !    Rubidium -       Oxygen
      xfac(37, 8) =     1.12535000d0 !    Rubidium -       Oxygen
      alpb(37, 9) =     3.63812200d0 !    Rubidium -     Fluorine
      xfac(37, 9) =    28.81527800d0 !    Rubidium -     Fluorine
      alpb(37,10) =     2.26759100d0 !    Rubidium -         Neon
      xfac(37,10) =     7.73656300d0 !    Rubidium -         Neon
      alpb(37,13) =     0.79877400d0 !    Rubidium -     Aluminum
      xfac(37,13) =     2.99245700d0 !    Rubidium -     Aluminum
      alpb(37,16) =     1.30318400d0 !    Rubidium -       Sulfur
      xfac(37,16) =     0.96441100d0 !    Rubidium -       Sulfur
      alpb(37,17) =     2.27441100d0 !    Rubidium -     Chlorine
      xfac(37,17) =    10.38448600d0 !    Rubidium -     Chlorine
      alpb(37,18) =     2.51097700d0 !    Rubidium -        Argon
      xfac(37,18) =    18.43332900d0 !    Rubidium -        Argon
      alpb(37,35) =     1.79776600d0 !    Rubidium -      Bromine
      xfac(37,35) =     5.17621400d0 !    Rubidium -      Bromine
      alpb(37,36) =     2.26875300d0 !    Rubidium -      Krypton
      xfac(37,36) =    15.30750300d0 !    Rubidium -      Krypton
      alpb(37,37) =     1.18081800d0 !    Rubidium -     Rubidium
      xfac(37,37) =    20.14761000d0 !    Rubidium -     Rubidium
 !
      alpb(38, 1) =     2.10591400d0 !   Strontium -     Hydrogen
      xfac(38, 1) =    12.97331600d0 !   Strontium -     Hydrogen
      alpb(38, 6) =     1.98668800d0 !   Strontium -       Carbon
      xfac(38, 6) =     6.65465700d0 !   Strontium -       Carbon
      alpb(38, 7) =     2.18362900d0 !   Strontium -     Nitrogen
      xfac(38, 7) =     6.85386600d0 !   Strontium -     Nitrogen
      alpb(38, 8) =     2.13839900d0 !   Strontium -       Oxygen
      xfac(38, 8) =     3.56139600d0 !   Strontium -       Oxygen
      alpb(38, 9) =     3.05066600d0 !   Strontium -     Fluorine
      xfac(38, 9) =    10.97170500d0 !   Strontium -     Fluorine
      alpb(38,14) =     2.96978000d0 !   Strontium -      Silicon
      xfac(38,14) =     2.76475000d0 !   Strontium -      Silicon
      alpb(38,15) =     2.78915000d0 !   Strontium -   Phosphorus
      xfac(38,15) =     2.55210000d0 !   Strontium -   Phosphorus
      alpb(38,16) =     1.59810600d0 !   Strontium -       Sulfur
      xfac(38,16) =     3.12960300d0 !   Strontium -       Sulfur
      alpb(38,17) =     1.85419000d0 !   Strontium -     Chlorine
      xfac(38,17) =     3.78395500d0 !   Strontium -     Chlorine
      alpb(38,22) =     2.88003000d0 !   Strontium -     Titanium
      xfac(38,22) =     2.81725000d0 !   Strontium -     Titanium
      alpb(38,35) =     1.52431600d0 !   Strontium -      Bromine
      xfac(38,35) =     2.76656700d0 !   Strontium -      Bromine
      alpb(38,38) =     1.00004000d0 !   Strontium -    Strontium
      xfac(38,38) =     5.37212000d0 !   Strontium -    Strontium
 !
      alpb(39, 1) =     1.18905300d0 !     Yttrium -     Hydrogen
      xfac(39, 1) =     0.61239900d0 !     Yttrium -     Hydrogen
      alpb(39, 6) =     1.33609400d0 !     Yttrium -       Carbon
      xfac(39, 6) =     0.50430600d0 !     Yttrium -       Carbon
      alpb(39, 7) =     1.77879600d0 !     Yttrium -     Nitrogen
      xfac(39, 7) =     1.62790300d0 !     Yttrium -     Nitrogen
      alpb(39, 8) =     1.85103000d0 !     Yttrium -       Oxygen
      xfac(39, 8) =     1.74292200d0 !     Yttrium -       Oxygen
      alpb(39, 9) =     2.64804600d0 !     Yttrium -     Fluorine
      xfac(39, 9) =     4.43380900d0 !     Yttrium -     Fluorine
      alpb(39,13) =     1.00350000d0 !     Yttrium -     Aluminum
      xfac(39,13) =     0.50067000d0 !     Yttrium -     Aluminum
      alpb(39,14) =     2.01682000d0 !     Yttrium -      Silicon
      xfac(39,14) =     3.21903000d0 !     Yttrium -      Silicon
      alpb(39,15) =     0.95445000d0 !     Yttrium -   Phosphorus
      xfac(39,15) =     0.54166000d0 !     Yttrium -   Phosphorus
      alpb(39,16) =     0.97168800d0 !     Yttrium -       Sulfur
      xfac(39,16) =     0.31822200d0 !     Yttrium -       Sulfur
      alpb(39,17) =     1.63015200d0 !     Yttrium -     Chlorine
      xfac(39,17) =     1.15495900d0 !     Yttrium -     Chlorine
      alpb(39,35) =     1.40120800d0 !     Yttrium -      Bromine
      xfac(39,35) =     1.05431600d0 !     Yttrium -      Bromine
      alpb(39,39) =     1.01268100d0 !     Yttrium -      Yttrium
      xfac(39,39) =     1.69172500d0 !     Yttrium -      Yttrium
 !
      alpb(40, 1) =     1.37970300d0 !   Zirconium -     Hydrogen
      xfac(40, 1) =     0.59373200d0 !   Zirconium -     Hydrogen
      alpb(40, 6) =     2.02942700d0 !   Zirconium -       Carbon
      xfac(40, 6) =     1.99918200d0 !   Zirconium -       Carbon
      alpb(40, 7) =     1.70708300d0 !   Zirconium -     Nitrogen
      xfac(40, 7) =     0.99504500d0 !   Zirconium -     Nitrogen
      alpb(40, 8) =     1.70957000d0 !   Zirconium -       Oxygen
      xfac(40, 8) =     1.05752500d0 !   Zirconium -       Oxygen
      alpb(40, 9) =     1.90092500d0 !   Zirconium -     Fluorine
      xfac(40, 9) =     0.86114200d0 !   Zirconium -     Fluorine
      alpb(40,13) =     1.27062000d0 !   Zirconium -     Aluminum
      xfac(40,13) =     0.87406000d0 !   Zirconium -     Aluminum
      alpb(40,14) =     1.75083300d0 !   Zirconium -      Silicon
      xfac(40,14) =     1.72334300d0 !   Zirconium -      Silicon
      alpb(40,15) =     1.09185800d0 !   Zirconium -   Phosphorus
      xfac(40,15) =     0.74837600d0 !   Zirconium -   Phosphorus
      alpb(40,16) =     2.12976100d0 !   Zirconium -       Sulfur
      xfac(40,16) =     2.42932400d0 !   Zirconium -       Sulfur
      alpb(40,17) =     1.32883500d0 !   Zirconium -     Chlorine
      xfac(40,17) =     0.44309900d0 !   Zirconium -     Chlorine
      alpb(40,35) =     1.44686800d0 !   Zirconium -      Bromine
      xfac(40,35) =     0.85890900d0 !   Zirconium -      Bromine
      alpb(40,40) =     3.86596800d0 !   Zirconium -    Zirconium
      xfac(40,40) =     3.07777300d0 !   Zirconium -    Zirconium
 !
      alpb(41, 1) =     2.50591200d0 !     Niobium -     Hydrogen
      xfac(41, 1) =     3.60377900d0 !     Niobium -     Hydrogen
      alpb(41, 6) =     2.62101200d0 !     Niobium -       Carbon
      xfac(41, 6) =     4.57548100d0 !     Niobium -       Carbon
      alpb(41, 7) =     2.02386300d0 !     Niobium -     Nitrogen
      xfac(41, 7) =     1.21358700d0 !     Niobium -     Nitrogen
      alpb(41, 8) =     2.04948900d0 !     Niobium -       Oxygen
      xfac(41, 8) =     1.18471900d0 !     Niobium -       Oxygen
      alpb(41, 9) =     3.00315700d0 !     Niobium -     Fluorine
      xfac(41, 9) =     3.66368200d0 !     Niobium -     Fluorine
      alpb(41,11) =     2.55101000d0 !     Niobium -       Sodium
      xfac(41,11) =     8.27602000d0 !     Niobium -       Sodium
      alpb(41,15) =     2.22160800d0 !     Niobium -   Phosphorus
      xfac(41,15) =     6.20150700d0 !     Niobium -   Phosphorus
      alpb(41,16) =     2.24948200d0 !     Niobium -       Sulfur
      xfac(41,16) =     2.46002000d0 !     Niobium -       Sulfur
      alpb(41,17) =     2.21527500d0 !     Niobium -     Chlorine
      xfac(41,17) =     1.89155700d0 !     Niobium -     Chlorine
      alpb(41,19) =     4.52136000d0 !     Niobium -    Potassium
      xfac(41,19) =     2.02659000d0 !     Niobium -    Potassium
      alpb(41,35) =     2.00667800d0 !     Niobium -      Bromine
      xfac(41,35) =     1.92126900d0 !     Niobium -      Bromine
      alpb(41,41) =     1.72794100d0 !     Niobium -      Niobium
      xfac(41,41) =     2.12238800d0 !     Niobium -      Niobium
 !
      alpb(42, 1) =     2.03574800d0 !  Molybdenum -     Hydrogen
      xfac(42, 1) =     0.93468600d0 !  Molybdenum -     Hydrogen
      alpb(42, 6) =     2.19867200d0 !  Molybdenum -       Carbon
      xfac(42, 6) =     1.19074200d0 !  Molybdenum -       Carbon
      alpb(42, 7) =     1.86947500d0 !  Molybdenum -     Nitrogen
      xfac(42, 7) =     0.60826800d0 !  Molybdenum -     Nitrogen
      alpb(42, 8) =     1.75542400d0 !  Molybdenum -       Oxygen
      xfac(42, 8) =     0.51126700d0 !  Molybdenum -       Oxygen
      alpb(42, 9) =     2.20259300d0 !  Molybdenum -     Fluorine
      xfac(42, 9) =     0.61042900d0 !  Molybdenum -     Fluorine
      alpb(42,11) =     2.44077000d0 !  Molybdenum -       Sodium
      xfac(42,11) =     8.28655000d0 !  Molybdenum -       Sodium
      alpb(42,15) =     1.85044100d0 !  Molybdenum -   Phosphorus
      xfac(42,15) =     1.52284600d0 !  Molybdenum -   Phosphorus
      alpb(42,16) =     1.93965800d0 !  Molybdenum -       Sulfur
      xfac(42,16) =     0.83042800d0 !  Molybdenum -       Sulfur
      alpb(42,17) =     1.78336200d0 !  Molybdenum -     Chlorine
      xfac(42,17) =     0.47432500d0 !  Molybdenum -     Chlorine
      alpb(42,19) =     3.93942000d0 !  Molybdenum -    Potassium
      xfac(42,19) =     2.14239000d0 !  Molybdenum -    Potassium
      alpb(42,24) =     2.67461600d0 !  Molybdenum -     Chromium
      xfac(42,24) =     1.74194300d0 !  Molybdenum -     Chromium
      alpb(42,35) =     1.28333400d0 !  Molybdenum -      Bromine
      xfac(42,35) =     0.22591800d0 !  Molybdenum -      Bromine
      alpb(42,42) =     2.03425400d0 !  Molybdenum -   Molybdenum
      xfac(42,42) =     0.62646200d0 !  Molybdenum -   Molybdenum
 !
      alpb(43, 1) =     2.83034500d0 !  Technetium -     Hydrogen
      xfac(43, 1) =     6.31033400d0 !  Technetium -     Hydrogen
      alpb(43, 6) =     3.19832600d0 !  Technetium -       Carbon
      xfac(43, 6) =     3.97243900d0 !  Technetium -       Carbon
      alpb(43, 7) =     2.31541700d0 !  Technetium -     Nitrogen
      xfac(43, 7) =     0.72713000d0 !  Technetium -     Nitrogen
      alpb(43, 8) =     2.40519000d0 !  Technetium -       Oxygen
      xfac(43, 8) =     1.02461600d0 !  Technetium -       Oxygen
      alpb(43, 9) =     3.60481500d0 !  Technetium -     Fluorine
      xfac(43, 9) =     5.81178400d0 !  Technetium -     Fluorine
      alpb(43,16) =     2.46340100d0 !  Technetium -       Sulfur
      xfac(43,16) =     1.49650200d0 !  Technetium -       Sulfur
      alpb(43,17) =     2.57204300d0 !  Technetium -     Chlorine
      xfac(43,17) =     1.65158300d0 !  Technetium -     Chlorine
      alpb(43,32) =     2.85282000d0 !  Technetium -    Germanium
      xfac(43,32) =     2.15206000d0 !  Technetium -    Germanium
      alpb(43,34) =     2.52366000d0 !  Technetium -     Selenium
      xfac(43,34) =     2.20262000d0 !  Technetium -     Selenium
      alpb(43,35) =     2.82826400d0 !  Technetium -      Bromine
      xfac(43,35) =     3.82013000d0 !  Technetium -      Bromine
 !
      alpb(44, 1) =     2.89289900d0 !   Ruthenium -     Hydrogen
      xfac(44, 1) =     7.13797600d0 !   Ruthenium -     Hydrogen
      alpb(44, 6) =     2.78483300d0 !   Ruthenium -       Carbon
      xfac(44, 6) =     1.13493600d0 !   Ruthenium -       Carbon
      alpb(44, 7) =     3.05550400d0 !   Ruthenium -     Nitrogen
      xfac(44, 7) =     2.33409400d0 !   Ruthenium -     Nitrogen
      alpb(44, 8) =     3.13494000d0 !   Ruthenium -       Oxygen
      xfac(44, 8) =     2.97627900d0 !   Ruthenium -       Oxygen
      alpb(44, 9) =     3.87871100d0 !   Ruthenium -     Fluorine
      xfac(44, 9) =     6.94712800d0 !   Ruthenium -     Fluorine
      alpb(44,14) =     2.77591000d0 !   Ruthenium -      Silicon
      xfac(44,14) =     0.84943000d0 !   Ruthenium -      Silicon
      alpb(44,15) =     0.29891600d0 !   Ruthenium -   Phosphorus
      xfac(44,15) =     0.05697400d0 !   Ruthenium -   Phosphorus
      alpb(44,16) =     2.50807600d0 !   Ruthenium -       Sulfur
      xfac(44,16) =     1.00668300d0 !   Ruthenium -       Sulfur
      alpb(44,17) =     1.75988300d0 !   Ruthenium -     Chlorine
      xfac(44,17) =     0.12658600d0 !   Ruthenium -     Chlorine
      alpb(44,32) =     2.85232000d0 !   Ruthenium -    Germanium
      xfac(44,32) =     2.15156000d0 !   Ruthenium -    Germanium
      alpb(44,34) =     2.52316000d0 !   Ruthenium -     Selenium
      xfac(44,34) =     2.20212000d0 !   Ruthenium -     Selenium
      alpb(44,35) =     2.58473500d0 !   Ruthenium -      Bromine
      xfac(44,35) =     0.65988100d0 !   Ruthenium -      Bromine
      alpb(44,44) =     0.57205600d0 !   Ruthenium -    Ruthenium
      xfac(44,44) =     0.09780500d0 !   Ruthenium -    Ruthenium
 !
      alpb(45, 1) =     3.10416500d0 !     Rhodium -     Hydrogen
      xfac(45, 1) =     2.30610700d0 !     Rhodium -     Hydrogen
      alpb(45, 6) =     3.41599100d0 !     Rhodium -       Carbon
      xfac(45, 6) =     3.48807900d0 !     Rhodium -       Carbon
      alpb(45, 7) =     3.58546200d0 !     Rhodium -     Nitrogen
      xfac(45, 7) =     4.00094700d0 !     Rhodium -     Nitrogen
      alpb(45, 8) =     3.92783000d0 !     Rhodium -       Oxygen
      xfac(45, 8) =    10.29867600d0 !     Rhodium -       Oxygen
      alpb(45, 9) =     4.05165400d0 !     Rhodium -     Fluorine
      xfac(45, 9) =     9.06538400d0 !     Rhodium -     Fluorine
      alpb(45,14) =     2.77649000d0 !     Rhodium -      Silicon
      xfac(45,14) =     0.85001000d0 !     Rhodium -      Silicon
      alpb(45,15) =     2.33460700d0 !     Rhodium -   Phosphorus
      xfac(45,15) =     1.03814100d0 !     Rhodium -   Phosphorus
      alpb(45,16) =     3.15400600d0 !     Rhodium -       Sulfur
      xfac(45,16) =     4.81641000d0 !     Rhodium -       Sulfur
      alpb(45,17) =     3.30013000d0 !     Rhodium -     Chlorine
      xfac(45,17) =     3.58686500d0 !     Rhodium -     Chlorine
      alpb(45,32) =     2.85290000d0 !     Rhodium -    Germanium
      xfac(45,32) =     2.15214000d0 !     Rhodium -    Germanium
      alpb(45,34) =     2.52374000d0 !     Rhodium -     Selenium
      xfac(45,34) =     2.20270000d0 !     Rhodium -     Selenium
      alpb(45,35) =     2.92808200d0 !     Rhodium -      Bromine
      xfac(45,35) =     1.51014900d0 !     Rhodium -      Bromine
      alpb(45,45) =     2.49732800d0 !     Rhodium -      Rhodium
      xfac(45,45) =     2.07011400d0 !     Rhodium -      Rhodium
 !
      alpb(46, 1) =     2.18376100d0 !   Palladium -     Hydrogen
      xfac(46, 1) =     0.44326900d0 !   Palladium -     Hydrogen
      alpb(46, 6) =     4.77719200d0 !   Palladium -       Carbon
      xfac(46, 6) =     9.85371500d0 !   Palladium -       Carbon
      alpb(46, 7) =     2.32804600d0 !   Palladium -     Nitrogen
      xfac(46, 7) =     0.24970300d0 !   Palladium -     Nitrogen
      alpb(46, 8) =     2.15486700d0 !   Palladium -       Oxygen
      xfac(46, 8) =     0.21640300d0 !   Palladium -       Oxygen
      alpb(46, 9) =     4.23731200d0 !   Palladium -     Fluorine
      xfac(46, 9) =     6.94531200d0 !   Palladium -     Fluorine
      alpb(46,13) =     1.57272000d0 !   Palladium -     Aluminum
      xfac(46,13) =     1.05729000d0 !   Palladium -     Aluminum
      alpb(46,14) =     2.94820000d0 !   Palladium -      Silicon
      xfac(46,14) =     2.22510400d0 !   Palladium -      Silicon
      alpb(46,15) =     0.80363000d0 !   Palladium -   Phosphorus
      xfac(46,15) =     0.04501700d0 !   Palladium -   Phosphorus
      alpb(46,16) =     2.17780100d0 !   Palladium -       Sulfur
      xfac(46,16) =     0.25522900d0 !   Palladium -       Sulfur
      alpb(46,17) =     3.87124300d0 !   Palladium -     Chlorine
      xfac(46,17) =     2.96989100d0 !   Palladium -     Chlorine
      alpb(46,35) =     5.99487900d0 !   Palladium -      Bromine
      xfac(46,35) =     4.63805100d0 !   Palladium -      Bromine
      alpb(46,46) =     1.06437500d0 !   Palladium -    Palladium
      xfac(46,46) =     0.05195600d0 !   Palladium -    Palladium
 !
      alpb(47, 1) =     2.89593600d0 !      Silver -     Hydrogen
      xfac(47, 1) =     1.99516800d0 !      Silver -     Hydrogen
      alpb(47, 6) =     4.40433600d0 !      Silver -       Carbon
      xfac(47, 6) =    11.33545600d0 !      Silver -       Carbon
      alpb(47, 7) =     4.65987100d0 !      Silver -     Nitrogen
      xfac(47, 7) =    19.80371000d0 !      Silver -     Nitrogen
      alpb(47, 8) =     1.89387400d0 !      Silver -       Oxygen
      xfac(47, 8) =     0.16566100d0 !      Silver -       Oxygen
      alpb(47, 9) =     4.62842300d0 !      Silver -     Fluorine
      xfac(47, 9) =    12.69588400d0 !      Silver -     Fluorine
      alpb(47,13) =     1.92880000d0 !      Silver -     Aluminum
      xfac(47,13) =     0.89651400d0 !      Silver -     Aluminum
      alpb(47,15) =     6.00000600d0 !      Silver -   Phosphorus
      xfac(47,15) =     0.04993200d0 !      Silver -   Phosphorus
      alpb(47,16) =     3.65312100d0 !      Silver -       Sulfur
      xfac(47,16) =    11.18802200d0 !      Silver -       Sulfur
      alpb(47,17) =     4.44117600d0 !      Silver -     Chlorine
      xfac(47,17) =    23.76545900d0 !      Silver -     Chlorine
      alpb(47,35) =     3.67749100d0 !      Silver -      Bromine
      xfac(47,35) =     1.71436900d0 !      Silver -      Bromine
      alpb(47,47) =     2.12764500d0 !      Silver -       Silver
      xfac(47,47) =     0.55774200d0 !      Silver -       Silver
 !
      alpb(48, 1) =     2.62874800d0 !     Cadmium -     Hydrogen
      xfac(48, 1) =    11.91420100d0 !     Cadmium -     Hydrogen
      alpb(48, 6) =     1.42567800d0 !     Cadmium -       Carbon
      xfac(48, 6) =     0.60344100d0 !     Cadmium -       Carbon
      alpb(48, 7) =     0.97042300d0 !     Cadmium -     Nitrogen
      xfac(48, 7) =     0.18066300d0 !     Cadmium -     Nitrogen
      alpb(48, 8) =     1.69667300d0 !     Cadmium -       Oxygen
      xfac(48, 8) =     0.92614600d0 !     Cadmium -       Oxygen
      alpb(48, 9) =     2.31213500d0 !     Cadmium -     Fluorine
      xfac(48, 9) =     1.35366500d0 !     Cadmium -     Fluorine
      alpb(48,14) =     1.37122500d0 !     Cadmium -      Silicon
      xfac(48,14) =     2.25334600d0 !     Cadmium -      Silicon
      alpb(48,16) =     1.18220200d0 !     Cadmium -       Sulfur
      xfac(48,16) =     0.36138900d0 !     Cadmium -       Sulfur
      alpb(48,17) =     0.94354700d0 !     Cadmium -     Chlorine
      xfac(48,17) =     0.14042400d0 !     Cadmium -     Chlorine
      alpb(48,35) =     1.00145100d0 !     Cadmium -      Bromine
      xfac(48,35) =     0.27226700d0 !     Cadmium -      Bromine
      alpb(48,48) =     1.56404400d0 !     Cadmium -      Cadmium
      xfac(48,48) =    18.61799900d0 !     Cadmium -      Cadmium
 !
      alpb(49, 1) =     3.06414400d0 !      Indium -     Hydrogen
      xfac(49, 1) =    14.97529300d0 !      Indium -     Hydrogen
      alpb(49, 6) =     2.18927200d0 !      Indium -       Carbon
      xfac(49, 6) =     2.18738500d0 !      Indium -       Carbon
      alpb(49, 7) =     2.46986800d0 !      Indium -     Nitrogen
      xfac(49, 7) =     3.36999300d0 !      Indium -     Nitrogen
      alpb(49, 8) =     2.66209500d0 !      Indium -       Oxygen
      xfac(49, 8) =     4.12858300d0 !      Indium -       Oxygen
      alpb(49, 9) =     2.94879700d0 !      Indium -     Fluorine
      xfac(49, 9) =     3.70101600d0 !      Indium -     Fluorine
      alpb(49,16) =     2.54213100d0 !      Indium -       Sulfur
      xfac(49,16) =     6.34110500d0 !      Indium -       Sulfur
      alpb(49,17) =     2.23340500d0 !      Indium -     Chlorine
      xfac(49,17) =     2.38855200d0 !      Indium -     Chlorine
      alpb(49,31) =     1.62887000d0 !      Indium -      Gallium
      xfac(49,31) =     2.42198700d0 !      Indium -      Gallium
      alpb(49,33) =     2.29955200d0 !      Indium -      Arsenic
      xfac(49,33) =     6.20835000d0 !      Indium -      Arsenic
      alpb(49,34) =     1.90657200d0 !      Indium -     Selenium
      xfac(49,34) =     2.31932300d0 !      Indium -     Selenium
      alpb(49,35) =     2.25795700d0 !      Indium -      Bromine
      xfac(49,35) =     3.72859800d0 !      Indium -      Bromine
      alpb(49,49) =     2.07324100d0 !      Indium -       Indium
      xfac(49,49) =     8.06349100d0 !      Indium -       Indium
 !
      alpb(50, 1) =     2.64891000d0 !         Tin -     Hydrogen
      xfac(50, 1) =     6.53516200d0 !         Tin -     Hydrogen
      alpb(50, 6) =     2.44053800d0 !         Tin -       Carbon
      xfac(50, 6) =     3.37435500d0 !         Tin -       Carbon
      alpb(50, 7) =     2.08558900d0 !         Tin -     Nitrogen
      xfac(50, 7) =     1.39190000d0 !         Tin -     Nitrogen
      alpb(50, 8) =     2.72726000d0 !         Tin -       Oxygen
      xfac(50, 8) =     4.37401700d0 !         Tin -       Oxygen
      alpb(50, 9) =     3.72428600d0 !         Tin -     Fluorine
      xfac(50, 9) =    18.59866400d0 !         Tin -     Fluorine
      alpb(50,16) =     2.13154200d0 !         Tin -       Sulfur
      xfac(50,16) =     2.31487000d0 !         Tin -       Sulfur
      alpb(50,17) =     1.77152200d0 !         Tin -     Chlorine
      xfac(50,17) =     0.80778200d0 !         Tin -     Chlorine
      alpb(50,32) =     2.52463300d0 !         Tin -    Germanium
      xfac(50,32) =    12.34341100d0 !         Tin -    Germanium
      alpb(50,34) =     2.12737700d0 !         Tin -     Selenium
      xfac(50,34) =     3.06188500d0 !         Tin -     Selenium
      alpb(50,35) =     1.53508900d0 !         Tin -      Bromine
      xfac(50,35) =     0.66879800d0 !         Tin -      Bromine
      alpb(50,50) =     0.92100000d0 !         Tin -          Tin
      xfac(50,50) =     0.28700000d0 !         Tin -          Tin
 !
      alpb(51, 1) =     1.57127200d0 !    Antimony -     Hydrogen
      xfac(51, 1) =     0.79534300d0 !    Antimony -     Hydrogen
      alpb(51, 6) =     1.69620600d0 !    Antimony -       Carbon
      xfac(51, 6) =     0.57921200d0 !    Antimony -       Carbon
      alpb(51, 7) =     0.67611500d0 !    Antimony -     Nitrogen
      xfac(51, 7) =     0.08206500d0 !    Antimony -     Nitrogen
      alpb(51, 8) =     1.84638400d0 !    Antimony -       Oxygen
      xfac(51, 8) =     0.63423400d0 !    Antimony -       Oxygen
      alpb(51, 9) =     2.18292200d0 !    Antimony -     Fluorine
      xfac(51, 9) =     0.65027700d0 !    Antimony -     Fluorine
      alpb(51,13) =     1.42264100d0 !    Antimony -     Aluminum
      xfac(51,13) =     1.61669000d0 !    Antimony -     Aluminum
      alpb(51,14) =     2.68659000d0 !    Antimony -      Silicon
      xfac(51,14) =     8.71374900d0 !    Antimony -      Silicon
      alpb(51,16) =     1.41883700d0 !    Antimony -       Sulfur
      xfac(51,16) =     0.39696900d0 !    Antimony -       Sulfur
      alpb(51,17) =     1.11728700d0 !    Antimony -     Chlorine
      xfac(51,17) =     0.15647500d0 !    Antimony -     Chlorine
      alpb(51,25) =     2.40032000d0 !    Antimony -    Manganese
      xfac(51,25) =     2.23671000d0 !    Antimony -    Manganese
      alpb(51,27) =     2.20463000d0 !    Antimony -       Cobalt
      xfac(51,27) =     2.27605000d0 !    Antimony -       Cobalt
      alpb(51,35) =     1.06391600d0 !    Antimony -      Bromine
      xfac(51,35) =     0.19804400d0 !    Antimony -      Bromine
      alpb(51,43) =     2.20485000d0 !    Antimony -   Technetium
      xfac(51,43) =     2.27626000d0 !    Antimony -   Technetium
      alpb(51,44) =     2.20435000d0 !    Antimony -    Ruthenium
      xfac(51,44) =     2.27576000d0 !    Antimony -    Ruthenium
      alpb(51,45) =     2.20493000d0 !    Antimony -      Rhodium
      xfac(51,45) =     2.27634000d0 !    Antimony -      Rhodium
      alpb(51,49) =     2.14193300d0 !    Antimony -       Indium
      xfac(51,49) =     6.66080100d0 !    Antimony -       Indium
      alpb(51,51) =     1.34853500d0 !    Antimony -     Antimony
      xfac(51,51) =     0.72488500d0 !    Antimony -     Antimony
 !
      alpb(52, 1) =     2.03913000d0 !   Tellurium -     Hydrogen
      xfac(52, 1) =     1.80767900d0 !   Tellurium -     Hydrogen
      alpb(52, 6) =     1.99281600d0 !   Tellurium -       Carbon
      xfac(52, 6) =     0.97049400d0 !   Tellurium -       Carbon
      alpb(52, 7) =     1.72226900d0 !   Tellurium -     Nitrogen
      xfac(52, 7) =     0.35859300d0 !   Tellurium -     Nitrogen
      alpb(52, 8) =     1.85306400d0 !   Tellurium -       Oxygen
      xfac(52, 8) =     0.38292600d0 !   Tellurium -       Oxygen
      alpb(52, 9) =     1.99857600d0 !   Tellurium -     Fluorine
      xfac(52, 9) =     0.20082200d0 !   Tellurium -     Fluorine
      alpb(52,13) =     1.38754100d0 !   Tellurium -     Aluminum
      xfac(52,13) =     2.10681200d0 !   Tellurium -     Aluminum
      alpb(52,15) =     1.45371800d0 !   Tellurium -   Phosphorus
      xfac(52,15) =     1.10928900d0 !   Tellurium -   Phosphorus
      alpb(52,16) =     1.83017000d0 !   Tellurium -       Sulfur
      xfac(52,16) =     0.94392500d0 !   Tellurium -       Sulfur
      alpb(52,17) =     1.30026000d0 !   Tellurium -     Chlorine
      xfac(52,17) =     0.28547800d0 !   Tellurium -     Chlorine
      alpb(52,30) =     1.17728821d0 !   Tellurium -         Zinc
      xfac(52,30) =     1.76233190d0 !   Tellurium -         Zinc
      alpb(52,32) =     2.34237200d0 !   Tellurium -    Germanium
      xfac(52,32) =     7.01904900d0 !   Tellurium -    Germanium
      alpb(52,33) =     1.18925300d0 !   Tellurium -      Arsenic
      xfac(52,33) =     0.68577400d0 !   Tellurium -      Arsenic
      alpb(52,34) =     1.56600800d0 !   Tellurium -     Selenium
      xfac(52,34) =     1.18782600d0 !   Tellurium -     Selenium
      alpb(52,35) =     1.25094000d0 !   Tellurium -      Bromine
      xfac(52,35) =     0.39420200d0 !   Tellurium -      Bromine
      alpb(52,48) =     1.30726200d0 !   Tellurium -      Cadmium
      xfac(52,48) =     1.08591900d0 !   Tellurium -      Cadmium
      alpb(52,49) =     1.54098800d0 !   Tellurium -       Indium
      xfac(52,49) =     2.03958200d0 !   Tellurium -       Indium
      alpb(52,50) =     1.76394100d0 !   Tellurium -          Tin
      xfac(52,50) =     2.95197600d0 !   Tellurium -          Tin
      alpb(52,52) =     1.16497800d0 !   Tellurium -    Tellurium
      xfac(52,52) =     0.64248600d0 !   Tellurium -    Tellurium
 !
      alpb(53, 1) =     1.93208961d0 !      Iodine -     Hydrogen
      xfac(53, 1) =     1.45911551d0 !      Iodine -     Hydrogen
      alpb(53, 2) =     2.17298400d0 !      Iodine -       Helium
      xfac(53, 2) =     1.63072100d0 !      Iodine -       Helium
      alpb(53, 3) =     2.12125100d0 !      Iodine -      Lithium
      xfac(53, 3) =     4.16859900d0 !      Iodine -      Lithium
      alpb(53, 4) =     2.28802300d0 !      Iodine -    Beryllium
      xfac(53, 4) =     2.35189800d0 !      Iodine -    Beryllium
      alpb(53, 5) =     2.66760500d0 !      Iodine -        Boron
      xfac(53, 5) =     3.16138500d0 !      Iodine -        Boron
      alpb(53, 6) =     1.74356688d0 !      Iodine -       Carbon
      xfac(53, 6) =     0.70734286d0 !      Iodine -       Carbon
      alpb(53, 7) =     1.70203519d0 !      Iodine -     Nitrogen
      xfac(53, 7) =     0.43292581d0 !      Iodine -     Nitrogen
      alpb(53, 8) =     2.26759981d0 !      Iodine -       Oxygen
      xfac(53, 8) =     1.24722053d0 !      Iodine -       Oxygen
      alpb(53, 9) =     1.38492211d0 !      Iodine -     Fluorine
      xfac(53, 9) =     0.04921958d0 !      Iodine -     Fluorine
      alpb(53,10) =     2.41441500d0 !      Iodine -         Neon
      xfac(53,10) =     1.50356800d0 !      Iodine -         Neon
      alpb(53,11) =     1.40309000d0 !      Iodine -       Sodium
      xfac(53,11) =     1.98611200d0 !      Iodine -       Sodium
      alpb(53,12) =     2.04513700d0 !      Iodine -    Magnesium
      xfac(53,12) =     3.27691400d0 !      Iodine -    Magnesium
      alpb(53,13) =     1.81606800d0 !      Iodine -     Aluminum
      xfac(53,13) =     2.92908000d0 !      Iodine -     Aluminum
      alpb(53,14) =     1.55957900d0 !      Iodine -      Silicon
      xfac(53,14) =     0.70029900d0 !      Iodine -      Silicon
      alpb(53,15) =     1.61394058d0 !      Iodine -   Phosphorus
      xfac(53,15) =     2.29228806d0 !      Iodine -   Phosphorus
      alpb(53,16) =     1.74862139d0 !      Iodine -       Sulfur
      xfac(53,16) =     1.22089065d0 !      Iodine -       Sulfur
      alpb(53,17) =     0.92363033d0 !      Iodine -     Chlorine
      xfac(53,17) =     0.06592143d0 !      Iodine -     Chlorine
      alpb(53,18) =     1.57658700d0 !      Iodine -        Argon
      xfac(53,18) =     0.30536700d0 !      Iodine -        Argon
      alpb(53,19) =     0.95633287d0 !      Iodine -    Potassium
      xfac(53,19) =     2.85564328d0 !      Iodine -    Potassium
      alpb(53,20) =     1.13608897d0 !      Iodine -      Calcium
      xfac(53,20) =     1.55760575d0 !      Iodine -      Calcium
      alpb(53,21) =     1.81488400d0 !      Iodine -     Scandium
      xfac(53,21) =     3.11428200d0 !      Iodine -     Scandium
      alpb(53,22) =     1.93346900d0 !      Iodine -     Titanium
      xfac(53,22) =     2.42674700d0 !      Iodine -     Titanium
      alpb(53,23) =     2.68352000d0 !      Iodine -     Vanadium
      xfac(53,23) =     6.19811200d0 !      Iodine -     Vanadium
      alpb(53,24) =     2.63422400d0 !      Iodine -     Chromium
      xfac(53,24) =     2.59859000d0 !      Iodine -     Chromium
      alpb(53,25) =     2.26660000d0 !      Iodine -    Manganese
      xfac(53,25) =     1.19341000d0 !      Iodine -    Manganese
      alpb(53,26) =     0.82428922d0 !      Iodine -         Iron
      xfac(53,26) =     0.09398336d0 !      Iodine -         Iron
      alpb(53,27) =     0.78587225d0 !      Iodine -       Cobalt
      xfac(53,27) =     0.07897620d0 !      Iodine -       Cobalt
      alpb(53,28) =     1.08534300d0 !      Iodine -       Nickel
      xfac(53,28) =     0.01745900d0 !      Iodine -       Nickel
      alpb(53,29) =     0.83430500d0 !      Iodine -       Copper
      xfac(53,29) =     0.00678100d0 !      Iodine -       Copper
      alpb(53,30) =     1.34495474d0 !      Iodine -         Zinc
      xfac(53,30) =     2.03209735d0 !      Iodine -         Zinc
      alpb(53,31) =     1.67172900d0 !      Iodine -      Gallium
      xfac(53,31) =     1.25216800d0 !      Iodine -      Gallium
      alpb(53,32) =     1.81742500d0 !      Iodine -    Germanium
      xfac(53,32) =     1.32326700d0 !      Iodine -    Germanium
      alpb(53,33) =     1.24526200d0 !      Iodine -      Arsenic
      xfac(53,33) =     0.31082400d0 !      Iodine -      Arsenic
      alpb(53,35) =     1.47849596d0 !      Iodine -      Bromine
      xfac(53,35) =     0.52925630d0 !      Iodine -      Bromine
      alpb(53,36) =     1.23857400d0 !      Iodine -      Krypton
      xfac(53,36) =     0.20113600d0 !      Iodine -      Krypton
      alpb(53,37) =     1.43267500d0 !      Iodine -     Rubidium
      xfac(53,37) =     4.09244600d0 !      Iodine -     Rubidium
      alpb(53,38) =     1.26204200d0 !      Iodine -    Strontium
      xfac(53,38) =     2.10394100d0 !      Iodine -    Strontium
      alpb(53,39) =     1.27911000d0 !      Iodine -      Yttrium
      xfac(53,39) =     1.02140200d0 !      Iodine -      Yttrium
      alpb(53,40) =     1.99518200d0 !      Iodine -    Zirconium
      xfac(53,40) =     4.51394300d0 !      Iodine -    Zirconium
      alpb(53,41) =     1.96725100d0 !      Iodine -      Niobium
      xfac(53,41) =     2.39929800d0 !      Iodine -      Niobium
      alpb(53,42) =     0.94846100d0 !      Iodine -   Molybdenum
      xfac(53,42) =     0.12469500d0 !      Iodine -   Molybdenum
      alpb(53,43) =     1.29231200d0 !      Iodine -   Technetium
      xfac(53,43) =     0.11059400d0 !      Iodine -   Technetium
      alpb(53,44) =     3.95320300d0 !      Iodine -    Ruthenium
      xfac(53,44) =     7.83771000d0 !      Iodine -    Ruthenium
      alpb(53,45) =     3.70817000d0 !      Iodine -      Rhodium
      xfac(53,45) =     2.35794400d0 !      Iodine -      Rhodium
      alpb(53,46) =     5.14454400d0 !      Iodine -    Palladium
      xfac(53,46) =     3.52201700d0 !      Iodine -    Palladium
      alpb(53,47) =     2.59316100d0 !      Iodine -       Silver
      xfac(53,47) =     0.04890400d0 !      Iodine -       Silver
      alpb(53,48) =     0.99623800d0 !      Iodine -      Cadmium
      xfac(53,48) =     0.39678400d0 !      Iodine -      Cadmium
      alpb(53,49) =     2.35175800d0 !      Iodine -       Indium
      xfac(53,49) =     5.94782100d0 !      Iodine -       Indium
      alpb(53,50) =     1.85563300d0 !      Iodine -          Tin
      xfac(53,50) =     1.78316300d0 !      Iodine -          Tin
      alpb(53,51) =     1.15531500d0 !      Iodine -     Antimony
      xfac(53,51) =     0.31819000d0 !      Iodine -     Antimony
      alpb(53,52) =     1.49395100d0 !      Iodine -    Tellurium
      xfac(53,52) =     1.10111600d0 !      Iodine -    Tellurium
      alpb(53,53) =     1.13078270d0 !      Iodine -       Iodine
      xfac(53,53) =     0.39946436d0 !      Iodine -       Iodine
 !
      alpb(54, 1) =     1.35686100d0 !       Xenon -     Hydrogen
      xfac(54, 1) =     0.70101600d0 !       Xenon -     Hydrogen
      alpb(54, 2) =     2.49783200d0 !       Xenon -       Helium
      xfac(54, 2) =     2.59947100d0 !       Xenon -       Helium
      alpb(54, 3) =     2.46689500d0 !       Xenon -      Lithium
      xfac(54, 3) =     4.58208100d0 !       Xenon -      Lithium
      alpb(54, 4) =     6.00000300d0 !       Xenon -    Beryllium
      xfac(54, 4) =     0.66052500d0 !       Xenon -    Beryllium
      alpb(54, 5) =     5.05195700d0 !       Xenon -        Boron
      xfac(54, 5) =     1.10061200d0 !       Xenon -        Boron
      alpb(54, 6) =     1.70444000d0 !       Xenon -       Carbon
      xfac(54, 6) =     0.82672700d0 !       Xenon -       Carbon
      alpb(54, 7) =     1.93295200d0 !       Xenon -     Nitrogen
      xfac(54, 7) =     0.92562400d0 !       Xenon -     Nitrogen
      alpb(54, 8) =     0.83923300d0 !       Xenon -       Oxygen
      xfac(54, 8) =     0.03535600d0 !       Xenon -       Oxygen
      alpb(54, 9) =     1.12881200d0 !       Xenon -     Fluorine
      xfac(54, 9) =     0.06501100d0 !       Xenon -     Fluorine
      alpb(54,10) =     1.33020200d0 !       Xenon -         Neon
      xfac(54,10) =     0.29386200d0 !       Xenon -         Neon
      alpb(54,11) =     2.12828932d0 !       Xenon -       Sodium
      xfac(54,11) =     8.36612980d0 !       Xenon -       Sodium
      alpb(54,12) =     2.69841400d0 !       Xenon -    Magnesium
      xfac(54,12) =     9.72357200d0 !       Xenon -    Magnesium
      alpb(54,13) =     2.41203900d0 !       Xenon -     Aluminum
      xfac(54,13) =     7.40446500d0 !       Xenon -     Aluminum
      alpb(54,14) =     3.08706000d0 !       Xenon -      Silicon
      xfac(54,14) =    16.09200000d0 !       Xenon -      Silicon
      alpb(54,17) =     1.54639600d0 !       Xenon -     Chlorine
      xfac(54,17) =     0.46375800d0 !       Xenon -     Chlorine
      alpb(54,18) =     0.59152000d0 !       Xenon -        Argon
      xfac(54,18) =     0.04926600d0 !       Xenon -        Argon
      alpb(54,19) =     1.17125000d0 !       Xenon -    Potassium
      xfac(54,19) =     1.22488900d0 !       Xenon -    Potassium
      alpb(54,20) =     1.51065300d0 !       Xenon -      Calcium
      xfac(54,20) =     1.71712100d0 !       Xenon -      Calcium
      alpb(54,35) =     1.43961800d0 !       Xenon -      Bromine
      xfac(54,35) =     0.47511600d0 !       Xenon -      Bromine
      alpb(54,36) =     0.55156100d0 !       Xenon -      Krypton
      xfac(54,36) =     0.04979300d0 !       Xenon -      Krypton
      alpb(54,37) =     1.08782300d0 !       Xenon -     Rubidium
      xfac(54,37) =     0.97496500d0 !       Xenon -     Rubidium
      alpb(54,53) =     0.79915500d0 !       Xenon -       Iodine
      xfac(54,53) =     0.11209000d0 !       Xenon -       Iodine
      alpb(54,54) =     1.24476200d0 !       Xenon -        Xenon
      xfac(54,54) =     0.34447400d0 !       Xenon -        Xenon
 !
      alpb(55, 1) =     0.26488200d0 !      Cesium -     Hydrogen
      xfac(55, 1) =     0.09690100d0 !      Cesium -     Hydrogen
      alpb(55, 5) =     1.48711000d0 !      Cesium -        Boron
      xfac(55, 5) =    10.39261000d0 !      Cesium -        Boron
      alpb(55, 6) =     2.14710400d0 !      Cesium -       Carbon
      xfac(55, 6) =    24.51462300d0 !      Cesium -       Carbon
      alpb(55, 7) =     2.44653200d0 !      Cesium -     Nitrogen
      xfac(55, 7) =    29.71107700d0 !      Cesium -     Nitrogen
      alpb(55, 8) =     2.08513900d0 !      Cesium -       Oxygen
      xfac(55, 8) =     8.17684300d0 !      Cesium -       Oxygen
      alpb(55, 9) =     2.83410000d0 !      Cesium -     Fluorine
      xfac(55, 9) =    22.23341600d0 !      Cesium -     Fluorine
      alpb(55,15) =     2.92495300d0 !      Cesium -   Phosphorus
      xfac(55,15) =     0.50651200d0 !      Cesium -   Phosphorus
      alpb(55,16) =     0.28941200d0 !      Cesium -       Sulfur
      xfac(55,16) =     0.09174300d0 !      Cesium -       Sulfur
      alpb(55,17) =     1.67366300d0 !      Cesium -     Chlorine
      xfac(55,17) =     4.53196500d0 !      Cesium -     Chlorine
      alpb(55,35) =     1.16718900d0 !      Cesium -      Bromine
      xfac(55,35) =     1.65842700d0 !      Cesium -      Bromine
      alpb(55,53) =     0.91956200d0 !      Cesium -       Iodine
      xfac(55,53) =     1.07217800d0 !      Cesium -       Iodine
      alpb(55,55) =     1.17084300d0 !      Cesium -       Cesium
      xfac(55,55) =    25.32005500d0 !      Cesium -       Cesium
 !
      alpb(56, 1) =     6.00013500d0 !      Barium -     Hydrogen
      xfac(56, 1) =     2.04000400d0 !      Barium -     Hydrogen
      alpb(56, 6) =     0.77062600d0 !      Barium -       Carbon
      xfac(56, 6) =     0.11979300d0 !      Barium -       Carbon
      alpb(56, 7) =     1.14823300d0 !      Barium -     Nitrogen
      xfac(56, 7) =     0.20793400d0 !      Barium -     Nitrogen
      alpb(56, 8) =     1.28301800d0 !      Barium -       Oxygen
      xfac(56, 8) =     0.34894500d0 !      Barium -       Oxygen
      alpb(56, 9) =     3.00061800d0 !      Barium -     Fluorine
      xfac(56, 9) =     5.57525500d0 !      Barium -     Fluorine
      alpb(56,13) =     2.10592400d0 !      Barium -     Aluminum
      xfac(56,13) =     9.53909900d0 !      Barium -     Aluminum
      alpb(56,14) =     1.24042000d0 !      Barium -      Silicon
      xfac(56,14) =     1.21266000d0 !      Barium -      Silicon
      alpb(56,16) =     0.70518800d0 !      Barium -       Sulfur
      xfac(56,16) =     0.21538600d0 !      Barium -       Sulfur
      alpb(56,17) =     1.07104400d0 !      Barium -     Chlorine
      xfac(56,17) =     0.16017700d0 !      Barium -     Chlorine
      alpb(56,22) =     2.17604000d0 !      Barium -     Titanium
      xfac(56,22) =     9.49353000d0 !      Barium -     Titanium
      alpb(56,35) =     1.19034600d0 !      Barium -      Bromine
      xfac(56,35) =     0.82879400d0 !      Barium -      Bromine
      alpb(56,53) =     0.98252800d0 !      Barium -       Iodine
      xfac(56,53) =     0.83559700d0 !      Barium -       Iodine
      alpb(56,56) =     0.33926900d0 !      Barium -       Barium
      xfac(56,56) =     0.35618600d0 !      Barium -       Barium
 !
      alpb(57, 1) =     0.83366700d0 !   Lanthanum -     Hydrogen
      xfac(57, 1) =     0.62350100d0 !   Lanthanum -     Hydrogen
      alpb(57, 6) =     0.60486900d0 !   Lanthanum -       Carbon
      xfac(57, 6) =     0.10864900d0 !   Lanthanum -       Carbon
      alpb(57, 7) =     0.75888100d0 !   Lanthanum -     Nitrogen
      xfac(57, 7) =     0.10477800d0 !   Lanthanum -     Nitrogen
      alpb(57, 8) =     1.31833300d0 !   Lanthanum -       Oxygen
      xfac(57, 8) =     0.55795700d0 !   Lanthanum -       Oxygen
      alpb(57, 9) =     2.37933500d0 !   Lanthanum -     Fluorine
      xfac(57, 9) =     2.40190300d0 !   Lanthanum -     Fluorine
      alpb(57,13) =     1.00351000d0 !   Lanthanum -     Aluminum
      xfac(57,13) =     0.50054000d0 !   Lanthanum -     Aluminum
      alpb(57,14) =     2.01682000d0 !   Lanthanum -      Silicon
      xfac(57,14) =     3.21903000d0 !   Lanthanum -      Silicon
      alpb(57,15) =     0.95445000d0 !   Lanthanum -   Phosphorus
      xfac(57,15) =     0.54166000d0 !   Lanthanum -   Phosphorus
      alpb(57,16) =     1.83412900d0 !   Lanthanum -       Sulfur
      xfac(57,16) =     2.68241200d0 !   Lanthanum -       Sulfur
      alpb(57,17) =     0.99375300d0 !   Lanthanum -     Chlorine
      xfac(57,17) =     0.23020300d0 !   Lanthanum -     Chlorine
      alpb(57,35) =     0.75818400d0 !   Lanthanum -      Bromine
      xfac(57,35) =     0.23858200d0 !   Lanthanum -      Bromine
      alpb(57,53) =     0.59266600d0 !   Lanthanum -       Iodine
      xfac(57,53) =     0.22688300d0 !   Lanthanum -       Iodine
      alpb(57,57) =     4.24806700d0 !   Lanthanum -    Lanthanum
      xfac(57,57) =     5.17516200d0 !   Lanthanum -    Lanthanum
 !
      alpb(64, 1) =     0.39087000d0 !  Gadolinium -     Hydrogen
      xfac(64, 1) =     0.13581000d0 !  Gadolinium -     Hydrogen
      alpb(64, 6) =     0.44687000d0 !  Gadolinium -       Carbon
      xfac(64, 6) =     0.05304000d0 !  Gadolinium -       Carbon
      alpb(64, 7) =     1.15941000d0 !  Gadolinium -     Nitrogen
      xfac(64, 7) =     0.20505000d0 !  Gadolinium -     Nitrogen
      alpb(64, 8) =     0.86204000d0 !  Gadolinium -       Oxygen
      xfac(64, 8) =     0.17580000d0 !  Gadolinium -       Oxygen
      alpb(64, 9) =     1.49798000d0 !  Gadolinium -     Fluorine
      xfac(64, 9) =     0.33463000d0 !  Gadolinium -     Fluorine
      alpb(64,13) =     1.00351000d0 !  Gadolinium -     Aluminum
      xfac(64,13) =     0.50054000d0 !  Gadolinium -     Aluminum
      alpb(64,14) =     2.01682000d0 !  Gadolinium -      Silicon
      xfac(64,14) =     3.21903000d0 !  Gadolinium -      Silicon
      alpb(64,15) =     0.95445000d0 !  Gadolinium -   Phosphorus
      xfac(64,15) =     0.54166000d0 !  Gadolinium -   Phosphorus
      alpb(64,16) =     2.00393000d0 !  Gadolinium -       Sulfur
      xfac(64,16) =     2.65540000d0 !  Gadolinium -       Sulfur
      alpb(64,17) =     0.80681000d0 !  Gadolinium -     Chlorine
      xfac(64,17) =     0.08997000d0 !  Gadolinium -     Chlorine
      alpb(64,35) =     0.71581000d0 !  Gadolinium -      Bromine
      xfac(64,35) =     0.24074000d0 !  Gadolinium -      Bromine
      alpb(64,53) =     0.58536000d0 !  Gadolinium -       Iodine
      xfac(64,53) =     0.27824000d0 !  Gadolinium -       Iodine
      alpb(64,64) =     3.34818000d0 !  Gadolinium -   Gadolinium
      xfac(64,64) =     2.67040000d0 !  Gadolinium -   Gadolinium
 !
      alpb(71, 1) =     1.41579000d0 !    Lutetium -     Hydrogen
      xfac(71, 1) =     0.78792000d0 !    Lutetium -     Hydrogen
      alpb(71, 6) =     2.31281300d0 !    Lutetium -       Carbon
      xfac(71, 6) =     4.45382500d0 !    Lutetium -       Carbon
      alpb(71, 7) =     2.14130200d0 !    Lutetium -     Nitrogen
      xfac(71, 7) =     2.86082800d0 !    Lutetium -     Nitrogen
      alpb(71, 8) =     2.19248600d0 !    Lutetium -       Oxygen
      xfac(71, 8) =     2.91707600d0 !    Lutetium -       Oxygen
      alpb(71,15) =     5.61882000d0 !    Lutetium -   Phosphorus
      xfac(71,15) =     0.50000000d0 !    Lutetium -   Phosphorus
      alpb(71,17) =     2.75363600d0 !    Lutetium -     Chlorine
      xfac(71,17) =    12.75709900d0 !    Lutetium -     Chlorine
      alpb(71,35) =     2.32261800d0 !    Lutetium -      Bromine
      xfac(71,35) =     8.64827400d0 !    Lutetium -      Bromine
      alpb(71,53) =     2.24834800d0 !    Lutetium -       Iodine
      xfac(71,53) =    10.08231500d0 !    Lutetium -       Iodine
 !
      alpb(72, 1) =     1.42378800d0 !     Hafnium -     Hydrogen
      xfac(72, 1) =     3.42731200d0 !     Hafnium -     Hydrogen
      alpb(72, 5) =     1.63350000d0 !     Hafnium -        Boron
      xfac(72, 5) =     0.65927000d0 !     Hafnium -        Boron
      alpb(72, 6) =     1.00219400d0 !     Hafnium -       Carbon
      xfac(72, 6) =     0.37857900d0 !     Hafnium -       Carbon
      alpb(72, 7) =     1.33241000d0 !     Hafnium -     Nitrogen
      xfac(72, 7) =     0.65579500d0 !     Hafnium -     Nitrogen
      alpb(72, 8) =     1.63328900d0 !     Hafnium -       Oxygen
      xfac(72, 8) =     1.03471800d0 !     Hafnium -       Oxygen
      alpb(72, 9) =     2.29080300d0 !     Hafnium -     Fluorine
      xfac(72, 9) =     1.67933500d0 !     Hafnium -     Fluorine
      alpb(72,12) =     1.91135000d0 !     Hafnium -    Magnesium
      xfac(72,12) =     4.33025000d0 !     Hafnium -    Magnesium
      alpb(72,13) =     0.94915000d0 !     Hafnium -     Aluminum
      xfac(72,13) =     0.62252000d0 !     Hafnium -     Aluminum
      alpb(72,14) =     2.18930000d0 !     Hafnium -      Silicon
      xfac(72,14) =     3.38230000d0 !     Hafnium -      Silicon
      alpb(72,15) =     1.23122000d0 !     Hafnium -   Phosphorus
      xfac(72,15) =     0.50553000d0 !     Hafnium -   Phosphorus
      alpb(72,16) =     2.32711000d0 !     Hafnium -       Sulfur
      xfac(72,16) =     1.66676000d0 !     Hafnium -       Sulfur
      alpb(72,17) =     1.29711700d0 !     Hafnium -     Chlorine
      xfac(72,17) =     0.70642100d0 !     Hafnium -     Chlorine
      alpb(72,20) =     2.05450000d0 !     Hafnium -      Calcium
      xfac(72,20) =     4.31951000d0 !     Hafnium -      Calcium
      alpb(72,33) =     1.79950000d0 !     Hafnium -      Arsenic
      xfac(72,33) =     1.28082000d0 !     Hafnium -      Arsenic
      alpb(72,35) =     1.09075900d0 !     Hafnium -      Bromine
      xfac(72,35) =     0.69245600d0 !     Hafnium -      Bromine
      alpb(72,53) =     1.01409600d0 !     Hafnium -       Iodine
      xfac(72,53) =     0.82094800d0 !     Hafnium -       Iodine
      alpb(72,56) =     2.26483000d0 !     Hafnium -       Barium
      xfac(72,56) =     9.02252000d0 !     Hafnium -       Barium
      alpb(72,72) =     0.54414400d0 !     Hafnium -      Hafnium
      xfac(72,72) =     1.05891100d0 !     Hafnium -      Hafnium
 !
      alpb(73, 1) =     2.28801400d0 !    Tantalum -     Hydrogen
      xfac(73, 1) =     2.82766900d0 !    Tantalum -     Hydrogen
      alpb(73, 6) =     1.83894900d0 !    Tantalum -       Carbon
      xfac(73, 6) =     0.84743900d0 !    Tantalum -       Carbon
      alpb(73, 7) =     2.05367900d0 !    Tantalum -     Nitrogen
      xfac(73, 7) =     1.01546100d0 !    Tantalum -     Nitrogen
      alpb(73, 8) =     2.41262900d0 !    Tantalum -       Oxygen
      xfac(73, 8) =     1.75108300d0 !    Tantalum -       Oxygen
      alpb(73, 9) =     3.10739000d0 !    Tantalum -     Fluorine
      xfac(73, 9) =     3.14652000d0 !    Tantalum -     Fluorine
      alpb(73,11) =     2.55112000d0 !    Tantalum -       Sodium
      xfac(73,11) =     8.27613000d0 !    Tantalum -       Sodium
      alpb(73,15) =     2.51380000d0 !    Tantalum -   Phosphorus
      xfac(73,15) =     6.26188000d0 !    Tantalum -   Phosphorus
      alpb(73,16) =     2.24672300d0 !    Tantalum -       Sulfur
      xfac(73,16) =     2.97598000d0 !    Tantalum -       Sulfur
      alpb(73,17) =     1.60880500d0 !    Tantalum -     Chlorine
      xfac(73,17) =     0.51641300d0 !    Tantalum -     Chlorine
      alpb(73,19) =     4.52147000d0 !    Tantalum -    Potassium
      xfac(73,19) =     2.02670000d0 !    Tantalum -    Potassium
      alpb(73,35) =     1.64037600d0 !    Tantalum -      Bromine
      xfac(73,35) =     0.79144500d0 !    Tantalum -      Bromine
      alpb(73,53) =     2.40105300d0 !    Tantalum -       Iodine
      xfac(73,53) =     6.55155100d0 !    Tantalum -       Iodine
      alpb(73,73) =     2.08286300d0 !    Tantalum -     Tantalum
      xfac(73,73) =    10.98705300d0 !    Tantalum -     Tantalum
 !
      alpb(74, 1) =     2.13088000d0 !    Tungsten -     Hydrogen
      xfac(74, 1) =     1.83227000d0 !    Tungsten -     Hydrogen
      alpb(74, 6) =     2.09748000d0 !    Tungsten -       Carbon
      xfac(74, 6) =     1.16077000d0 !    Tungsten -       Carbon
      alpb(74, 7) =     1.59604000d0 !    Tungsten -     Nitrogen
      xfac(74, 7) =     0.47835000d0 !    Tungsten -     Nitrogen
      alpb(74, 8) =     1.35902000d0 !    Tungsten -       Oxygen
      xfac(74, 8) =     0.34901000d0 !    Tungsten -       Oxygen
      alpb(74, 9) =     1.44605000d0 !    Tungsten -     Fluorine
      xfac(74, 9) =     0.21389000d0 !    Tungsten -     Fluorine
      alpb(74,11) =     2.55103000d0 !    Tungsten -       Sodium
      xfac(74,11) =     8.27604000d0 !    Tungsten -       Sodium
      alpb(74,15) =     2.33806000d0 !    Tungsten -   Phosphorus
      xfac(74,15) =     5.95386000d0 !    Tungsten -   Phosphorus
      alpb(74,16) =     1.54257000d0 !    Tungsten -       Sulfur
      xfac(74,16) =     0.48863000d0 !    Tungsten -       Sulfur
      alpb(74,17) =     1.31069000d0 !    Tungsten -     Chlorine
      xfac(74,17) =     0.27800000d0 !    Tungsten -     Chlorine
      alpb(74,19) =     4.52138000d0 !    Tungsten -    Potassium
      xfac(74,19) =     2.02661000d0 !    Tungsten -    Potassium
      alpb(74,35) =     1.29326000d0 !    Tungsten -      Bromine
      xfac(74,35) =     0.37239000d0 !    Tungsten -      Bromine
      alpb(74,53) =     1.57357000d0 !    Tungsten -       Iodine
      xfac(74,53) =     1.07737000d0 !    Tungsten -       Iodine
      alpb(74,74) =     2.94087000d0 !    Tungsten -     Tungsten
      xfac(74,74) =     7.47139000d0 !    Tungsten -     Tungsten
 !
      alpb(75, 1) =     1.63450000d0 !     Rhenium -     Hydrogen
      xfac(75, 1) =     0.34589400d0 !     Rhenium -     Hydrogen
      alpb(75, 6) =     2.30628500d0 !     Rhenium -       Carbon
      xfac(75, 6) =     0.69068700d0 !     Rhenium -       Carbon
      alpb(75, 7) =     1.91833200d0 !     Rhenium -     Nitrogen
      xfac(75, 7) =     0.44521300d0 !     Rhenium -     Nitrogen
      alpb(75, 8) =     1.96774700d0 !     Rhenium -       Oxygen
      xfac(75, 8) =     0.63596000d0 !     Rhenium -       Oxygen
      alpb(75, 9) =     2.15421900d0 !     Rhenium -     Fluorine
      xfac(75, 9) =     0.53596600d0 !     Rhenium -     Fluorine
      alpb(75,14) =     2.77593000d0 !     Rhenium -      Silicon
      xfac(75,14) =     0.84945000d0 !     Rhenium -      Silicon
      alpb(75,15) =     1.80416800d0 !     Rhenium -   Phosphorus
      xfac(75,15) =     0.96694200d0 !     Rhenium -   Phosphorus
      alpb(75,16) =     1.08391900d0 !     Rhenium -       Sulfur
      xfac(75,16) =     0.06887400d0 !     Rhenium -       Sulfur
      alpb(75,17) =     1.43387500d0 !     Rhenium -     Chlorine
      xfac(75,17) =     0.14631900d0 !     Rhenium -     Chlorine
      alpb(75,32) =     2.85234000d0 !     Rhenium -    Germanium
      xfac(75,32) =     2.15158000d0 !     Rhenium -    Germanium
      alpb(75,34) =     2.52317000d0 !     Rhenium -     Selenium
      xfac(75,34) =     2.20214000d0 !     Rhenium -     Selenium
      alpb(75,35) =     1.60306000d0 !     Rhenium -      Bromine
      xfac(75,35) =     0.28752800d0 !     Rhenium -      Bromine
      alpb(75,51) =     2.20436000d0 !     Rhenium -     Antimony
      xfac(75,51) =     2.27578000d0 !     Rhenium -     Antimony
      alpb(75,53) =     2.61011900d0 !     Rhenium -       Iodine
      xfac(75,53) =     3.55928600d0 !     Rhenium -       Iodine
      alpb(75,75) =     6.00025800d0 !     Rhenium -      Rhenium
      xfac(75,75) =     4.48885200d0 !     Rhenium -      Rhenium
 !
      alpb(76, 1) =     3.40418000d0 !      Osmium -     Hydrogen
      xfac(76, 1) =     4.39387000d0 !      Osmium -     Hydrogen
      alpb(76, 6) =     2.33650000d0 !      Osmium -       Carbon
      xfac(76, 6) =     0.49841000d0 !      Osmium -       Carbon
      alpb(76, 7) =     1.14309000d0 !      Osmium -     Nitrogen
      xfac(76, 7) =     0.08087000d0 !      Osmium -     Nitrogen
      alpb(76, 8) =     1.35036000d0 !      Osmium -       Oxygen
      xfac(76, 8) =     0.18430000d0 !      Osmium -       Oxygen
      alpb(76, 9) =     1.50762000d0 !      Osmium -     Fluorine
      xfac(76, 9) =     0.14005000d0 !      Osmium -     Fluorine
      alpb(76,11) =     2.55074000d0 !      Osmium -       Sodium
      xfac(76,11) =     8.27575000d0 !      Osmium -       Sodium
      alpb(76,15) =     2.83609000d0 !      Osmium -   Phosphorus
      xfac(76,15) =     6.05830000d0 !      Osmium -   Phosphorus
      alpb(76,16) =     2.80950000d0 !      Osmium -       Sulfur
      xfac(76,16) =     4.18605000d0 !      Osmium -       Sulfur
      alpb(76,17) =     1.83307000d0 !      Osmium -     Chlorine
      xfac(76,17) =     0.32792000d0 !      Osmium -     Chlorine
      alpb(76,19) =     4.52109000d0 !      Osmium -    Potassium
      xfac(76,19) =     2.02632000d0 !      Osmium -    Potassium
      alpb(76,35) =     1.76688000d0 !      Osmium -      Bromine
      xfac(76,35) =     0.38243000d0 !      Osmium -      Bromine
      alpb(76,53) =     2.20376000d0 !      Osmium -       Iodine
      xfac(76,53) =     2.19919000d0 !      Osmium -       Iodine
      alpb(76,76) =     2.02163000d0 !      Osmium -       Osmium
      xfac(76,76) =     0.83044000d0 !      Osmium -       Osmium
 !
      alpb(77, 1) =     1.03390000d0 !     Iridium -     Hydrogen
      xfac(77, 1) =     0.05804700d0 !     Iridium -     Hydrogen
      alpb(77, 6) =     1.69029500d0 !     Iridium -       Carbon
      xfac(77, 6) =     0.11504700d0 !     Iridium -       Carbon
      alpb(77, 7) =     3.93450800d0 !     Iridium -     Nitrogen
      xfac(77, 7) =     8.51864000d0 !     Iridium -     Nitrogen
      alpb(77, 8) =     3.74827200d0 !     Iridium -       Oxygen
      xfac(77, 8) =     9.62540200d0 !     Iridium -       Oxygen
      alpb(77, 9) =     2.98279900d0 !     Iridium -     Fluorine
      xfac(77, 9) =     1.49963900d0 !     Iridium -     Fluorine
      alpb(77,11) =     2.55082000d0 !     Iridium -       Sodium
      xfac(77,11) =     8.27583000d0 !     Iridium -       Sodium
      alpb(77,15) =     2.71406000d0 !     Iridium -   Phosphorus
      xfac(77,15) =     6.28467000d0 !     Iridium -   Phosphorus
      alpb(77,16) =     3.20483400d0 !     Iridium -       Sulfur
      xfac(77,16) =     4.13573200d0 !     Iridium -       Sulfur
      alpb(77,17) =     2.00977000d0 !     Iridium -     Chlorine
      xfac(77,17) =     0.25891600d0 !     Iridium -     Chlorine
      alpb(77,19) =     4.52117000d0 !     Iridium -    Potassium
      xfac(77,19) =     2.02640000d0 !     Iridium -    Potassium
      alpb(77,35) =     2.03814200d0 !     Iridium -      Bromine
      xfac(77,35) =     0.17187900d0 !     Iridium -      Bromine
      alpb(77,53) =     3.41091400d0 !     Iridium -       Iodine
      xfac(77,53) =     1.49714800d0 !     Iridium -       Iodine
      alpb(77,77) =     5.77166300d0 !     Iridium -      Iridium
      xfac(77,77) =    11.17519300d0 !     Iridium -      Iridium
 !
      alpb(78, 1) =     4.00119800d0 !    Platinum -     Hydrogen
      xfac(78, 1) =     8.92401500d0 !    Platinum -     Hydrogen
      alpb(78, 6) =     3.30672200d0 !    Platinum -       Carbon
      xfac(78, 6) =     3.49340300d0 !    Platinum -       Carbon
      alpb(78, 7) =     2.30792300d0 !    Platinum -     Nitrogen
      xfac(78, 7) =     0.54073000d0 !    Platinum -     Nitrogen
      alpb(78, 8) =     2.11056300d0 !    Platinum -       Oxygen
      xfac(78, 8) =     0.48775600d0 !    Platinum -       Oxygen
      alpb(78, 9) =     3.71444100d0 !    Platinum -     Fluorine
      xfac(78, 9) =     5.61701400d0 !    Platinum -     Fluorine
      alpb(78,13) =     1.57236000d0 !    Platinum -     Aluminum
      xfac(78,13) =     1.05693000d0 !    Platinum -     Aluminum
      alpb(78,14) =     0.99999000d0 !    Platinum -      Silicon
      xfac(78,14) =     0.09999000d0 !    Platinum -      Silicon
      alpb(78,15) =     1.40323900d0 !    Platinum -   Phosphorus
      xfac(78,15) =     0.23371200d0 !    Platinum -   Phosphorus
      alpb(78,16) =     2.79150000d0 !    Platinum -       Sulfur
      xfac(78,16) =     2.22426300d0 !    Platinum -       Sulfur
      alpb(78,17) =     2.10852600d0 !    Platinum -     Chlorine
      xfac(78,17) =     0.34100100d0 !    Platinum -     Chlorine
      alpb(78,35) =     2.18530700d0 !    Platinum -      Bromine
      xfac(78,35) =     0.52036100d0 !    Platinum -      Bromine
      alpb(78,53) =     3.07733800d0 !    Platinum -       Iodine
      xfac(78,53) =     4.60124800d0 !    Platinum -       Iodine
      alpb(78,78) =     3.40427600d0 !    Platinum -     Platinum
      xfac(78,78) =     9.01025200d0 !    Platinum -     Platinum
 !
      alpb(79, 1) =     3.36904100d0 !        Gold -     Hydrogen
      xfac(79, 1) =     2.60528300d0 !        Gold -     Hydrogen
      alpb(79, 6) =     4.58001600d0 !        Gold -       Carbon
      xfac(79, 6) =    21.48563400d0 !        Gold -       Carbon
      alpb(79, 7) =     2.13809500d0 !        Gold -     Nitrogen
      xfac(79, 7) =     0.22205900d0 !        Gold -     Nitrogen
      alpb(79, 8) =     1.54876300d0 !        Gold -       Oxygen
      xfac(79, 8) =     0.07719200d0 !        Gold -       Oxygen
      alpb(79, 9) =     4.45314500d0 !        Gold -     Fluorine
      xfac(79, 9) =     9.59438400d0 !        Gold -     Fluorine
      alpb(79,13) =     1.57257000d0 !        Gold -     Aluminum
      xfac(79,13) =     1.05714000d0 !        Gold -     Aluminum
      alpb(79,15) =     1.61871300d0 !        Gold -   Phosphorus
      xfac(79,15) =     0.06700100d0 !        Gold -   Phosphorus
      alpb(79,16) =     4.30623800d0 !        Gold -       Sulfur
      xfac(79,16) =    21.61914500d0 !        Gold -       Sulfur
      alpb(79,17) =     3.53941400d0 !        Gold -     Chlorine
      xfac(79,17) =     2.25770200d0 !        Gold -     Chlorine
      alpb(79,35) =     0.58191100d0 !        Gold -      Bromine
      xfac(79,35) =     0.00423700d0 !        Gold -      Bromine
      alpb(79,53) =     0.57791600d0 !        Gold -       Iodine
      xfac(79,53) =     0.00881600d0 !        Gold -       Iodine
      alpb(79,79) =     0.90316200d0 !        Gold -         Gold
      xfac(79,79) =     0.01309100d0 !        Gold -         Gold
 !
      alpb(80, 1) =     1.13658700d0 !     Mercury -     Hydrogen
      xfac(80, 1) =     0.79939900d0 !     Mercury -     Hydrogen
      alpb(80, 6) =     0.79581600d0 !     Mercury -       Carbon
      xfac(80, 6) =     0.14712800d0 !     Mercury -       Carbon
      alpb(80, 7) =     0.33215200d0 !     Mercury -     Nitrogen
      xfac(80, 7) =     0.05024000d0 !     Mercury -     Nitrogen
      alpb(80, 8) =     1.05214500d0 !     Mercury -       Oxygen
      xfac(80, 8) =     0.24072000d0 !     Mercury -       Oxygen
      alpb(80, 9) =     1.24057200d0 !     Mercury -     Fluorine
      xfac(80, 9) =     0.11382700d0 !     Mercury -     Fluorine
      alpb(80,14) =     2.77086000d0 !     Mercury -      Silicon
      xfac(80,14) =     3.68074000d0 !     Mercury -      Silicon
      alpb(80,15) =     0.60860400d0 !     Mercury -   Phosphorus
      xfac(80,15) =     0.21495100d0 !     Mercury -   Phosphorus
      alpb(80,16) =     1.04168200d0 !     Mercury -       Sulfur
      xfac(80,16) =     0.34738300d0 !     Mercury -       Sulfur
      alpb(80,17) =     0.43073100d0 !     Mercury -     Chlorine
      xfac(80,17) =     0.05366000d0 !     Mercury -     Chlorine
      alpb(80,22) =     3.41463000d0 !     Mercury -     Titanium
      xfac(80,22) =     2.95720000d0 !     Mercury -     Titanium
      alpb(80,35) =     0.63871700d0 !     Mercury -      Bromine
      xfac(80,35) =     0.17236300d0 !     Mercury -      Bromine
      alpb(80,52) =     0.29150000d0 !     Mercury -    Tellurium
      xfac(80,52) =     0.21273200d0 !     Mercury -    Tellurium
      alpb(80,53) =     0.75816200d0 !     Mercury -       Iodine
      xfac(80,53) =     0.34205800d0 !     Mercury -       Iodine
      alpb(80,80) =     0.47441300d0 !     Mercury -      Mercury
      xfac(80,80) =     0.42327600d0 !     Mercury -      Mercury
 !
      alpb(81, 1) =     0.67365800d0 !    Thallium -     Hydrogen
      xfac(81, 1) =     0.13820500d0 !    Thallium -     Hydrogen
      alpb(81, 5) =     1.52834700d0 !    Thallium -        Boron
      xfac(81, 5) =    10.50433800d0 !    Thallium -        Boron
      alpb(81, 6) =     1.39034500d0 !    Thallium -       Carbon
      xfac(81, 6) =     0.58289500d0 !    Thallium -       Carbon
      alpb(81, 7) =     0.98233500d0 !    Thallium -     Nitrogen
      xfac(81, 7) =     0.15881200d0 !    Thallium -     Nitrogen
      alpb(81, 8) =     1.55006800d0 !    Thallium -       Oxygen
      xfac(81, 8) =     0.63690600d0 !    Thallium -       Oxygen
      alpb(81, 9) =     1.46951600d0 !    Thallium -     Fluorine
      xfac(81, 9) =     0.22616600d0 !    Thallium -     Fluorine
      alpb(81,16) =     0.99485100d0 !    Thallium -       Sulfur
      xfac(81,16) =     0.30342600d0 !    Thallium -       Sulfur
      alpb(81,17) =     0.84619300d0 !    Thallium -     Chlorine
      xfac(81,17) =     0.16203700d0 !    Thallium -     Chlorine
      alpb(81,35) =     0.87441900d0 !    Thallium -      Bromine
      xfac(81,35) =     0.29683600d0 !    Thallium -      Bromine
      alpb(81,53) =     0.90201200d0 !    Thallium -       Iodine
      xfac(81,53) =     0.43003300d0 !    Thallium -       Iodine
      alpb(81,81) =     1.19168400d0 !    Thallium -     Thallium
      xfac(81,81) =     9.53512700d0 !    Thallium -     Thallium
 !
      alpb(82, 1) =     1.52267600d0 !        Lead -     Hydrogen
      xfac(82, 1) =     0.84009600d0 !        Lead -     Hydrogen
      alpb(82, 3) =     1.00181000d0 !        Lead -      Lithium
      xfac(82, 3) =     1.28506400d0 !        Lead -      Lithium
      alpb(82, 5) =     0.91119700d0 !        Lead -        Boron
      xfac(82, 5) =     1.13815700d0 !        Lead -        Boron
      alpb(82, 6) =     1.52559300d0 !        Lead -       Carbon
      xfac(82, 6) =     0.40465600d0 !        Lead -       Carbon
      alpb(82, 7) =     1.31739400d0 !        Lead -     Nitrogen
      xfac(82, 7) =     0.33578700d0 !        Lead -     Nitrogen
      alpb(82, 8) =     1.76321000d0 !        Lead -       Oxygen
      xfac(82, 8) =     0.78250600d0 !        Lead -       Oxygen
      alpb(82, 9) =     3.28890200d0 !        Lead -     Fluorine
      xfac(82, 9) =     8.36856200d0 !        Lead -     Fluorine
      alpb(82,15) =     4.51680000d0 !        Lead -   Phosphorus
      xfac(82,15) =     5.03320000d0 !        Lead -   Phosphorus
      alpb(82,16) =     1.02751900d0 !        Lead -       Sulfur
      xfac(82,16) =     0.17515000d0 !        Lead -       Sulfur
      alpb(82,17) =     1.09412300d0 !        Lead -     Chlorine
      xfac(82,17) =     0.16481400d0 !        Lead -     Chlorine
      alpb(82,23) =     1.50000000d0 !        Lead -     Vanadium
      xfac(82,23) =     1.00000000d0 !        Lead -     Vanadium
      alpb(82,24) =     1.86076000d0 !        Lead -     Chromium
      xfac(82,24) =     1.02911000d0 !        Lead -     Chromium
      alpb(82,30) =     1.50000000d0 !        Lead -         Zinc
      xfac(82,30) =     1.00000000d0 !        Lead -         Zinc
      alpb(82,34) =     2.00000000d0 !        Lead -     Selenium
      xfac(82,34) =     0.11119500d0 !        Lead -     Selenium
      alpb(82,35) =     0.86555000d0 !        Lead -      Bromine
      xfac(82,35) =     0.14822900d0 !        Lead -      Bromine
      alpb(82,41) =     1.50000000d0 !        Lead -      Niobium
      xfac(82,41) =     1.00000000d0 !        Lead -      Niobium
      alpb(82,42) =     2.00000000d0 !        Lead -   Molybdenum
      xfac(82,42) =     5.00000000d0 !        Lead -   Molybdenum
      alpb(82,52) =     1.00255900d0 !        Lead -    Tellurium
      xfac(82,52) =     0.80904200d0 !        Lead -    Tellurium
      alpb(82,53) =     0.98347400d0 !        Lead -       Iodine
      xfac(82,53) =     0.26742600d0 !        Lead -       Iodine
      alpb(82,82) =     1.88176400d0 !        Lead -         Lead
      xfac(82,82) =     2.36234300d0 !        Lead -         Lead
 !
      alpb(83, 1) =     1.67990500d0 !     Bismuth -     Hydrogen
      xfac(83, 1) =     1.39746200d0 !     Bismuth -     Hydrogen
      alpb(83, 3) =     0.34014000d0 !     Bismuth -      Lithium
      xfac(83, 3) =     0.69532000d0 !     Bismuth -      Lithium
      alpb(83, 6) =     1.53402500d0 !     Bismuth -       Carbon
      xfac(83, 6) =     0.57617900d0 !     Bismuth -       Carbon
      alpb(83, 7) =     1.14387600d0 !     Bismuth -     Nitrogen
      xfac(83, 7) =     0.15273800d0 !     Bismuth -     Nitrogen
      alpb(83, 8) =     1.55329700d0 !     Bismuth -       Oxygen
      xfac(83, 8) =     0.33304200d0 !     Bismuth -       Oxygen
      alpb(83, 9) =     2.35540000d0 !     Bismuth -     Fluorine
      xfac(83, 9) =     1.03532400d0 !     Bismuth -     Fluorine
      alpb(83,16) =     1.46687900d0 !     Bismuth -       Sulfur
      xfac(83,16) =     0.62099700d0 !     Bismuth -       Sulfur
      alpb(83,17) =     1.27297500d0 !     Bismuth -     Chlorine
      xfac(83,17) =     0.32687100d0 !     Bismuth -     Chlorine
      alpb(83,34) =     1.34474600d0 !     Bismuth -     Selenium
      xfac(83,34) =     0.65120800d0 !     Bismuth -     Selenium
      alpb(83,35) =     1.14623300d0 !     Bismuth -      Bromine
      xfac(83,35) =     0.38117000d0 !     Bismuth -      Bromine
      alpb(83,53) =     1.30217100d0 !     Bismuth -       Iodine
      xfac(83,53) =     0.86237700d0 !     Bismuth -       Iodine
      alpb(83,83) =     1.07406400d0 !     Bismuth -      Bismuth
      xfac(83,83) =     1.16821400d0 !     Bismuth -      Bismuth
 !
      alpb(87, 7) =     2.21881000d0 !    Francium -     Nitrogen
      xfac(87, 7) =     1.01263000d0 !    Francium -     Nitrogen
      alpb(87, 9) =     2.21881000d0 !    Francium -     Fluorine
      xfac(87, 9) =     1.01263000d0 !    Francium -     Fluorine
      alpb(87,17) =     1.57966000d0 !    Francium -     Chlorine
      xfac(87,17) =     0.76156000d0 !    Francium -     Chlorine
      alpb(87,87) =     1.57966000d0 !    Francium -     Francium
      xfac(87,87) =     0.76156000d0 !    Francium -     Francium
    end subroutine alpb_and_xfac_pm6_ORG
  end module Parameters_for_PM6_ORG_C
