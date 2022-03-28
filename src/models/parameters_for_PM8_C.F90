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

  module Parameters_for_PM8_C
    double precision, dimension(107) :: uss8, upp8, udd8, zs8, zp8, zd8, betas8, &
    betap8, betad8, gss8, gsp8, gpp8, gp28, hsp8, polvo8, poc_8, &
    zsn8, zpn8, zdn8, f0sd8, g2sd8, alp8, &
    CPE_Zet8, CPE_Z08, CPE_B8, CPE_Xlo8, CPE_Xhi8
    double precision :: v_par8(60)
    double precision, dimension(107,4) :: gues81, gues82, gues83
!
!                    Data for Element   1         Hydrogen
!
      data     uss8(  1)/       -11.246958D0/
      data   betas8(  1)/        -8.352984D0/
      data      zs8(  1)/         1.268641D0/
      data     gss8(  1)/        14.448686D0/
      data   polvo8(  1)/         0.262114D0/
      data CPE_Zet8(  1)/         1.342400D0/
      data  CPE_Z08(  1)/         2.255100D0/
      data   CPE_B8(  1)/         0.856600D0/
      data CPE_Xlo8(  1)/         0.379600D0/
      data CPE_Xhi8(  1)/         1.379600D0/
      data gues81(  1,1)/         0.024184D0/
      data gues82(  1,1)/         3.055953D0/
      data gues83(  1,1)/         1.786011D0/
!
!                    Data for Element   2           Helium
!
      data     uss8(  2)/       -31.770969D0/
      data     upp8(  2)/        -5.856382D0/
      data   betas8(  2)/       -58.903774D0/
      data   betap8(  2)/       -37.039974D0/
      data      zs8(  2)/         3.313204D0/
      data      zp8(  2)/         3.657133D0/
      data     gss8(  2)/         9.445299D0/
      data     gsp8(  2)/        11.201419D0/
      data     gpp8(  2)/         9.214548D0/
      data     gp28(  2)/        13.046115D0/
      data     hsp8(  2)/         0.299954D0/
!
!                    Data for Element   3          Lithium
!
      data     uss8(  3)/        -4.709912D0/
      data     upp8(  3)/        -2.722581D0/
      data   betas8(  3)/        -2.283946D0/
      data   betap8(  3)/        -7.535573D0/
      data      zs8(  3)/         0.981041D0/
      data      zp8(  3)/         2.953445D0/
      data     gss8(  3)/        11.035907D0/
      data     gsp8(  3)/        19.998647D0/
      data     gpp8(  3)/        11.543650D0/
      data     gp28(  3)/         9.059036D0/
      data     hsp8(  3)/         1.641886D0/
!
!                    Data for Element   4        Beryllium
!
      data     uss8(  4)/       -16.360315D0/
      data     upp8(  4)/       -16.339216D0/
      data   betas8(  4)/        -3.199549D0/
      data   betap8(  4)/        -4.451920D0/
      data      zs8(  4)/         1.212539D0/
      data      zp8(  4)/         1.276487D0/
      data     gss8(  4)/         7.552804D0/
      data     gsp8(  4)/        10.203146D0/
      data     gpp8(  4)/        12.862153D0/
      data     gp28(  4)/        13.602858D0/
      data     hsp8(  4)/         1.501452D0/
      data gues81(  4,1)/         0.164180D0/
      data gues82(  4,1)/         1.704828D0/
      data gues83(  4,1)/         1.785591D0/
!
!                    Data for Element   5            Boron
!
      data     uss8(  5)/       -25.967679D0/
      data     upp8(  5)/       -19.115864D0/
      data   betas8(  5)/        -4.959706D0/
      data   betap8(  5)/        -4.656753D0/
      data      zs8(  5)/         1.634174D0/
      data      zp8(  5)/         1.479195D0/
      data     gss8(  5)/         8.179341D0/
      data     gsp8(  5)/         7.294021D0/
      data     gpp8(  5)/         7.829395D0/
      data     gp28(  5)/         6.401072D0/
      data     hsp8(  5)/         1.252845D0/
!
!                    Data for Element   6           Carbon
!
      data     uss8(  6)/       -51.089653D0/
      data     upp8(  6)/       -39.937920D0/
      data   betas8(  6)/       -15.385236D0/
      data   betap8(  6)/        -7.471929D0/
      data      zs8(  6)/         2.047558D0/
      data      zp8(  6)/         1.702841D0/
      data     gss8(  6)/        13.335519D0/
      data     gsp8(  6)/        11.528134D0/
      data     gpp8(  6)/        10.778326D0/
      data     gp28(  6)/         9.486212D0/
      data     hsp8(  6)/         0.717322D0/
      data   polvo8(  6)/         0.485071D0/
      data CPE_Zet8(  6)/         1.167040D0/
      data  CPE_Z08(  6)/         1.178300D0/
      data   CPE_B8(  6)/         0.004800D0/
      data CPE_Xlo8(  6)/         1.086200D0/
      data CPE_Xhi8(  6)/         2.353000D0/
      data gues81(  6,1)/         0.046302D0/
      data gues82(  6,1)/         2.100206D0/
      data gues83(  6,1)/         1.333959D0/
!
!                    Data for Element   7         Nitrogen
!
      data     uss8(  7)/       -57.784823D0/
      data     upp8(  7)/       -49.893036D0/
      data   betas8(  7)/       -17.979377D0/
      data   betap8(  7)/       -15.055017D0/
      data      zs8(  7)/         2.380406D0/
      data      zp8(  7)/         1.999246D0/
      data     gss8(  7)/        12.357026D0/
      data     gsp8(  7)/         9.636190D0/
      data     gpp8(  7)/        12.570756D0/
      data     gp28(  7)/        10.576425D0/
      data     hsp8(  7)/         2.871545D0/
      data   polvo8(  7)/         0.204743D0/
      data CPE_Zet8(  7)/         1.378880D0/
      data  CPE_Z08(  7)/         2.029200D0/
      data   CPE_B8(  7)/         0.323800D0/
      data CPE_Xlo8(  7)/         1.651100D0/
      data CPE_Xhi8(  7)/         2.292100D0/
      data gues81(  7,1)/        -0.001436D0/
      data gues82(  7,1)/         0.495196D0/
      data gues83(  7,1)/         1.704857D0/
!
!                    Data for Element   8           Oxygen
!
      data     uss8(  8)/       -91.678761D0/
      data     upp8(  8)/       -70.460949D0/
      data   betas8(  8)/       -65.635137D0/
      data   betap8(  8)/       -21.622604D0/
      data      zs8(  8)/         5.421751D0/
      data      zp8(  8)/         2.270960D0/
      data     gss8(  8)/        11.304042D0/
      data     gsp8(  8)/        15.807424D0/
      data     gpp8(  8)/        13.618205D0/
      data     gp28(  8)/        10.332765D0/
      data     hsp8(  8)/         5.010801D0/
      data   polvo8(  8)/         0.154301D0/
      data CPE_Zet8(  8)/         1.585280D0/
      data  CPE_Z08(  8)/         4.322700D0/
      data   CPE_B8(  8)/         0.045100D0/
      data CPE_Xlo8(  8)/         3.483200D0/
      data CPE_Xhi8(  8)/         3.605000D0/
      data gues81(  8,1)/        -0.017771D0/
      data gues82(  8,1)/         3.058310D0/
      data gues83(  8,1)/         1.896435D0/
!
!                    Data for Element   9         Fluorine
!
      data     uss8(  9)/      -140.225626D0/
      data     upp8(  9)/       -98.778044D0/
      data   betas8(  9)/       -69.922593D0/
      data   betap8(  9)/       -30.448165D0/
      data      zs8(  9)/         6.043849D0/
      data      zp8(  9)/         2.906722D0/
      data     gss8(  9)/        12.446818D0/
      data     gsp8(  9)/        18.496082D0/
      data     gpp8(  9)/         8.417366D0/
      data     gp28(  9)/        12.179816D0/
      data     hsp8(  9)/         2.604382D0/
      data   polvo8(  9)/         0.199611D0/
      data gues81(  9,1)/        -0.010792D0/
      data gues82(  9,1)/         6.004648D0/
      data gues83(  9,1)/         1.847724D0/
!
!                    Data for Element  10             Neon
!
      data     uss8( 10)/        -2.978729D0/
      data     upp8( 10)/       -85.441118D0/
      data   betas8( 10)/       -69.793475D0/
      data   betap8( 10)/       -33.261962D0/
      data      zs8( 10)/         6.000148D0/
      data      zp8( 10)/         3.834528D0/
      data     gss8( 10)/        19.999574D0/
      data     gsp8( 10)/        16.896951D0/
      data     gpp8( 10)/         8.963560D0/
      data     gp28( 10)/        16.027799D0/
      data     hsp8( 10)/         1.779280D0/
!
!                    Data for Element  11           Sodium
!
      data     uss8( 11)/        -4.537153D0/
      data     upp8( 11)/        -2.433015D0/
      data   betas8( 11)/         0.244853D0/
      data   betap8( 11)/         0.491998D0/
      data      zs8( 11)/         0.686327D0/
      data      zp8( 11)/         0.950068D0/
      data     gss8( 11)/         4.059972D0/
      data     gsp8( 11)/         7.061183D0/
      data     gpp8( 11)/         9.283540D0/
      data     gp28( 11)/        17.034978D0/
      data     hsp8( 11)/         0.640715D0/
      data gues81( 11,1)/        -1.026036D0/
      data gues82( 11,1)/         2.014506D0/
      data gues83( 11,1)/         1.271202D0/
!
!                    Data for Element  12        Magnesium
!
      data     uss8( 12)/       -14.574226D0/
      data     upp8( 12)/        -7.583850D0/
      data   betas8( 12)/        -9.604932D0/
      data   betap8( 12)/         3.416908D0/
      data      zs8( 12)/         1.310830D0/
      data      zp8( 12)/         1.388897D0/
      data     gss8( 12)/         7.115328D0/
      data     gsp8( 12)/         3.253024D0/
      data     gpp8( 12)/         4.737311D0/
      data     gp28( 12)/         8.428485D0/
      data     hsp8( 12)/         0.877379D0/
!
!                    Data for Element  13         Aluminum
!
      data     uss8( 13)/       -24.546778D0/
      data     upp8( 13)/       -20.104434D0/
      data     udd8( 13)/         8.004394D0/
      data   betas8( 13)/       -18.375229D0/
      data   betap8( 13)/        -9.382700D0/
      data   betad8( 13)/       -20.840474D0/
      data      zs8( 13)/         2.364264D0/
      data      zp8( 13)/         1.749102D0/
      data      zd8( 13)/         1.269384D0/
      data     zsn8( 13)/         4.742341D0/
      data     zpn8( 13)/         4.669626D0/
      data     zdn8( 13)/         7.131138D0/
      data     alp8( 13)/         0.968798D0/
      data     gss8( 13)/         6.652155D0/
      data     gsp8( 13)/         7.459435D0/
      data     gpp8( 13)/         7.668857D0/
      data     gp28( 13)/         6.673299D0/
      data     hsp8( 13)/         0.435060D0/
      data gues81( 13,1)/         1.002222D0/
      data gues82( 13,1)/         1.517400D0/
      data gues83( 13,1)/         0.659101D0/
!
!                    Data for Element  14          Silicon
!
      data     uss8( 14)/       -27.358058D0/
      data     upp8( 14)/       -20.490578D0/
      data     udd8( 14)/       -22.751900D0/
      data   betas8( 14)/        -8.686909D0/
      data   betap8( 14)/        -1.856482D0/
      data   betad8( 14)/        -6.360627D0/
      data      zs8( 14)/         1.752741D0/
      data      zp8( 14)/         1.198413D0/
      data      zd8( 14)/         2.128593D0/
      data     zsn8( 14)/         8.388111D0/
      data     zpn8( 14)/         1.843048D0/
      data     zdn8( 14)/         0.708600D0/
      data     gss8( 14)/         5.194805D0/
      data     gsp8( 14)/         5.090534D0/
      data     gpp8( 14)/         5.185150D0/
      data     gp28( 14)/         4.769775D0/
      data     hsp8( 14)/         1.425012D0/
      data   polvo8( 14)/         1.886110D0/
      data gues81( 14,1)/         0.208571D0/
      data gues82( 14,1)/         6.000483D0/
      data gues83( 14,1)/         1.185245D0/
!
!                    Data for Element  15       Phosphorus
!
      data     uss8( 15)/       -48.729905D0/
      data     upp8( 15)/       -40.354689D0/
      data     udd8( 15)/        -7.349246D0/
      data   betas8( 15)/       -14.583780D0/
      data   betap8( 15)/       -11.744725D0/
      data   betad8( 15)/       -20.099893D0/
      data      zs8( 15)/         2.158033D0/
      data      zp8( 15)/         1.805343D0/
      data      zd8( 15)/         1.230358D0/
      data     zsn8( 15)/         6.042706D0/
      data     zpn8( 15)/         2.376473D0/
      data     zdn8( 15)/         7.147750D0/
      data     gss8( 15)/         8.758856D0/
      data     gsp8( 15)/         8.483679D0/
      data     gpp8( 15)/         8.662754D0/
      data     gp28( 15)/         7.734264D0/
      data     hsp8( 15)/         0.871681D0/
      data   polvo8( 15)/         2.314610D0/
      data gues81( 15,1)/        -0.034320D0/
      data gues82( 15,1)/         6.001394D0/
      data gues83( 15,1)/         2.296737D0/
!
!                    Data for Element  16           Sulfur
!
      data     uss8( 16)/       -47.530706D0/
      data     upp8( 16)/       -39.191045D0/
      data     udd8( 16)/       -46.306944D0/
      data   betas8( 16)/       -13.827440D0/
      data   betap8( 16)/        -7.664613D0/
      data   betad8( 16)/        -9.986172D0/
      data      zs8( 16)/         2.192844D0/
      data      zp8( 16)/         1.841078D0/
      data      zd8( 16)/         3.109401D0/
      data     zsn8( 16)/         0.479722D0/
      data     zpn8( 16)/         1.015507D0/
      data     zdn8( 16)/         4.317470D0/
      data     gss8( 16)/         9.170350D0/
      data     gsp8( 16)/         5.944296D0/
      data     gpp8( 16)/         8.165473D0/
      data     gp28( 16)/         7.301878D0/
      data     hsp8( 16)/         5.005404D0/
      data   polvo8( 16)/         1.453310D0/
      data gues81( 16,1)/        -0.036928D0/
      data gues82( 16,1)/         1.795067D0/
      data gues83( 16,1)/         2.082618D0/
!
!                    Data for Element  17         Chlorine
!
      data     uss8( 17)/       -61.389930D0/
      data     upp8( 17)/       -54.482801D0/
      data     udd8( 17)/       -38.258155D0/
      data   betas8( 17)/        -2.367988D0/
      data   betap8( 17)/       -13.802139D0/
      data   betad8( 17)/        -4.037751D0/
      data      zs8( 17)/         2.637050D0/
      data      zp8( 17)/         2.118146D0/
      data      zd8( 17)/         1.324033D0/
      data     zsn8( 17)/         0.956297D0/
      data     zpn8( 17)/         2.464067D0/
      data     zdn8( 17)/         6.410325D0/
      data     gss8( 17)/        11.142654D0/
      data     gsp8( 17)/         7.487881D0/
      data     gpp8( 17)/         9.551886D0/
      data     gp28( 17)/         8.128436D0/
      data     hsp8( 17)/         5.004267D0/
      data   polvo8( 17)/         1.236210D0/
      data gues81( 17,1)/        -0.013213D0/
      data gues82( 17,1)/         3.687022D0/
      data gues83( 17,1)/         2.544635D0/
!
!                    Data for Element  18            Argon
!
      data     uss8( 18)/        -7.797931D0/
      data     upp8( 18)/       -83.211487D0/
      data   betas8( 18)/        -8.839842D0/
      data   betap8( 18)/       -28.427303D0/
      data      zs8( 18)/         6.000272D0/
      data      zp8( 18)/         5.949170D0/
      data     gss8( 18)/        17.858776D0/
      data     gsp8( 18)/         4.168451D0/
      data     gpp8( 18)/        11.852500D0/
      data     gp28( 18)/        15.669543D0/
      data     hsp8( 18)/         4.574549D0/
!
!                    Data for Element  19        Potassium
!
      data     uss8( 19)/        -3.801108D0/
      data     upp8( 19)/        -3.339656D0/
      data   betas8( 19)/        -8.755195D0/
      data   betap8( 19)/        -1.788061D0/
      data      zs8( 19)/         6.000478D0/
      data      zp8( 19)/         1.127503D0/
      data     gss8( 19)/         3.369251D0/
      data     gsp8( 19)/         6.129351D0/
      data     gpp8( 19)/         0.999505D0/
      data     gp28( 19)/        18.999148D0/
      data     hsp8( 19)/         0.300325D0/
      data gues81( 19,1)/         0.157519D0/
      data gues82( 19,1)/         6.000566D0/
      data gues83( 19,1)/         2.047539D0/
!
!                    Data for Element  20          Calcium
!
      data     uss8( 20)/       -10.770058D0/
      data     upp8( 20)/        -9.754177D0/
      data   betas8( 20)/        -4.343881D0/
      data   betap8( 20)/        -1.296612D0/
      data      zs8( 20)/         1.528258D0/
      data      zp8( 20)/         2.060094D0/
      data     gss8( 20)/         5.725773D0/
      data     gsp8( 20)/         4.781065D0/
      data     gpp8( 20)/         7.172103D0/
      data     gp28( 20)/         7.431876D0/
      data     hsp8( 20)/         1.240572D0/
      data gues81( 20,1)/        -0.025275D0/
      data gues82( 20,1)/         0.500017D0/
      data gues83( 20,1)/         2.329051D0/
!
!                    Data for Element  21         Scandium
!
      data     uss8( 21)/       -15.544461D0/
      data     upp8( 21)/       -18.646295D0/
      data     udd8( 21)/       -16.069444D0/
      data   betas8( 21)/        -8.620944D0/
      data   betap8( 21)/         3.075948D0/
      data   betad8( 21)/        -9.768661D0/
      data      zs8( 21)/         1.402469D0/
      data      zp8( 21)/         1.345196D0/
      data      zd8( 21)/         1.859012D0/
      data     zsn8( 21)/         0.848418D0/
      data     zpn8( 21)/         2.451729D0/
      data     zdn8( 21)/         0.789372D0/
      data     alp8( 21)/         0.816556D0/
      data     gss8( 21)/         4.638215D0/
      data     gsp8( 21)/         5.739164D0/
      data     gpp8( 21)/        14.604872D0/
      data     gp28( 21)/        12.802595D0/
      data     hsp8( 21)/         0.193835D0/
      data    poc_8( 21)/         3.173734D0/
      data    f0sd8( 21)/         4.798313D0/
      data    g2sd8( 21)/         5.380136D0/
!
!                    Data for Element  22         Titanium
!
      data     uss8( 22)/       -25.507973D0/
      data     upp8( 22)/       -17.260909D0/
      data     udd8( 22)/       -23.809486D0/
      data   betas8( 22)/         3.389142D0/
      data   betap8( 22)/        -3.355350D0/
      data   betad8( 22)/        -1.842829D0/
      data      zs8( 22)/         5.324777D0/
      data      zp8( 22)/         1.164068D0/
      data      zd8( 22)/         1.418280D0/
      data     zsn8( 22)/         1.045904D0/
      data     zpn8( 22)/         1.076844D0/
      data     zdn8( 22)/         0.717945D0/
      data     gss8( 22)/         5.717851D0/
      data     gsp8( 22)/         5.800015D0/
      data     gpp8( 22)/         6.414726D0/
      data     gp28( 22)/         5.623133D0/
      data     hsp8( 22)/         1.403732D0/
      data    f0sd8( 22)/         6.560562D0/
      data    g2sd8( 22)/         3.396235D0/
!
!                    Data for Element  23         Vanadium
!
      data     uss8( 23)/       -32.162276D0/
      data     upp8( 23)/       -21.572501D0/
      data     udd8( 23)/       -34.506245D0/
      data   betas8( 23)/        -1.211330D0/
      data   betap8( 23)/         0.740746D0/
      data   betad8( 23)/         3.153669D0/
      data      zs8( 23)/         1.974330D0/
      data      zp8( 23)/         1.063106D0/
      data      zd8( 23)/         1.394806D0/
      data     zsn8( 23)/         1.094426D0/
      data     zpn8( 23)/         0.755378D0/
      data     zdn8( 23)/         1.099367D0/
      data     gss8( 23)/         5.983116D0/
      data     gsp8( 23)/         4.736769D0/
      data     gpp8( 23)/         4.499763D0/
      data     gp28( 23)/         3.944481D0/
      data     hsp8( 23)/         0.901105D0/
      data    f0sd8( 23)/         6.810021D0/
      data    g2sd8( 23)/         1.831407D0/
!
!                    Data for Element  24         Chromium
!
      data     uss8( 24)/       -34.864339D0/
      data     upp8( 24)/       -26.978615D0/
      data     udd8( 24)/       -54.431036D0/
      data   betas8( 24)/        -5.122615D0/
      data   betap8( 24)/         3.926711D0/
      data   betad8( 24)/        -4.230550D0/
      data      zs8( 24)/         3.283460D0/
      data      zp8( 24)/         1.029394D0/
      data      zd8( 24)/         1.623119D0/
      data     zsn8( 24)/         1.619853D0/
      data     zpn8( 24)/         0.848266D0/
      data     zdn8( 24)/         1.405015D0/
      data     gss8( 24)/         8.855572D0/
      data     gsp8( 24)/         5.588631D0/
      data     gpp8( 24)/         5.053094D0/
      data     gp28( 24)/         4.429530D0/
      data     hsp8( 24)/         0.648039D0/
      data    f0sd8( 24)/         6.150136D0/
      data    g2sd8( 24)/         2.000300D0/
!
!                    Data for Element  25        Manganese
!
      data     uss8( 25)/       -51.460000D0/
      data     upp8( 25)/       -37.543990D0/
      data     udd8( 25)/       -47.655370D0/
      data   betas8( 25)/        -4.185290D0/
      data   betap8( 25)/        -3.479630D0/
      data   betad8( 25)/       -13.473190D0/
      data      zs8( 25)/         2.131680D0/
      data      zp8( 25)/         1.525880D0/
      data      zd8( 25)/         2.607800D0/
      data     zsn8( 25)/         1.132450D0/
      data     zpn8( 25)/         1.390740D0/
      data     zdn8( 25)/         0.962550D0/
      data     gss8( 25)/         6.190990D0/
      data     gsp8( 25)/         6.757427D0/
      data     gpp8( 25)/         8.284594D0/
      data     gp28( 25)/         7.262255D0/
      data     hsp8( 25)/         1.520518D0/
      data    f0sd8( 25)/         7.690920D0/
      data    g2sd8( 25)/         1.105330D0/
!
!                    Data for Element  26             Iron
!
      data     uss8( 26)/       -70.515047D0/
      data     upp8( 26)/       -62.963069D0/
      data     udd8( 26)/      -103.631790D0/
      data   betas8( 26)/         8.027621D0/
      data   betap8( 26)/        -1.125760D0/
      data   betad8( 26)/        -3.507531D0/
      data      zs8( 26)/         1.479150D0/
      data      zp8( 26)/         6.002246D0/
      data      zd8( 26)/         1.080747D0/
      data     zsn8( 26)/         1.459152D0/
      data     zpn8( 26)/         1.392614D0/
      data     zdn8( 26)/         2.161909D0/
      data     gss8( 26)/         7.977036D0/
      data     gsp8( 26)/         7.786867D0/
      data     gpp8( 26)/         8.295758D0/
      data     gp28( 26)/         7.272041D0/
      data     hsp8( 26)/         1.880189D0/
      data    poc_8( 26)/         1.272092D0/
      data    f0sd8( 26)/         9.300165D0/
      data    g2sd8( 26)/         1.601345D0/
!
!                    Data for Element  27           Cobalt
!
      data     uss8( 27)/       -21.039413D0/
      data     upp8( 27)/        10.000000D0/
      data     udd8( 27)/       -28.068971D0/
      data   betas8( 27)/        -8.992062D0/
      data   betap8( 27)/        -0.100000D0/
      data   betad8( 27)/        -2.481509D0/
      data      zs8( 27)/         1.166613D0/
      data      zp8( 27)/         3.000000D0/
      data      zd8( 27)/         1.860218D0/
      data     zsn8( 27)/         0.519518D0/
      data     zpn8( 27)/         1.000000D0/
      data     zdn8( 27)/         0.352115D0/
      data     gss8( 27)/         2.840152D0/
      data     gsp8( 27)/         3.425933D0/
      data     gpp8( 27)/         5.956968D0/
      data     gp28( 27)/         5.221864D0/
      data     hsp8( 27)/         0.390087D0/
      data    f0sd8( 27)/         1.446283D0/
      data    g2sd8( 27)/         1.680225D0/
!
!                    Data for Element  28           Nickel
!
      data     uss8( 28)/       -47.620247D0/
      data     upp8( 28)/       -32.878408D0/
      data     udd8( 28)/       -93.026395D0/
      data   betas8( 28)/        -9.151521D0/
      data   betap8( 28)/        -8.086696D0/
      data   betad8( 28)/        -8.655910D0/
      data      zs8( 28)/         1.591828D0/
      data      zp8( 28)/         2.304739D0/
      data      zd8( 28)/         2.514761D0/
      data     zsn8( 28)/         0.746470D0/
      data     zpn8( 28)/         0.753327D0/
      data     zdn8( 28)/         1.461345D0/
      data     alp8( 28)/         2.894960D0/
      data     gss8( 28)/         4.080876D0/
      data     gsp8( 28)/         4.099452D0/
      data     gpp8( 28)/         4.487545D0/
      data     gp28( 28)/         3.933771D0/
      data     hsp8( 28)/         0.993498D0/
      data    poc_8( 28)/         1.586979D0/
      data    f0sd8( 28)/         4.651664D0/
      data    g2sd8( 28)/         1.880502D0/
!
!                    Data for Element  29           Copper
!
      data     uss8( 29)/       -97.002205D0/
      data     upp8( 29)/        -1.000000D0/
      data     udd8( 29)/      -110.442592D0/
      data   betas8( 29)/        -9.369508D0/
      data   betap8( 29)/        -0.100000D0/
      data   betad8( 29)/       -16.982092D0/
      data      zs8( 29)/         1.669096D0/
      data      zp8( 29)/         3.000000D0/
      data      zd8( 29)/         2.734990D0/
      data     zsn8( 29)/         1.899598D0/
      data     zpn8( 29)/         3.000000D0/
      data     zdn8( 29)/         1.484317D0/
      data     gss8( 29)/        10.384910D0/
      data     gsp8( 29)/        12.145361D0/
      data     gpp8( 29)/        17.870905D0/
      data     gp28( 29)/        15.665592D0/
      data     hsp8( 29)/         2.037394D0/
      data    f0sd8( 29)/         9.848807D0/
      data    g2sd8( 29)/         9.847577D0/
!
!                    Data for Element  30             Zinc
!
      data     uss8( 30)/       -18.040862D0/
      data     upp8( 30)/        -7.834895D0/
      data   betas8( 30)/       -13.276583D0/
      data   betap8( 30)/         1.479642D0/
      data      zs8( 30)/         1.512875D0/
      data      zp8( 30)/         1.789482D0/
      data     gss8( 30)/         8.707424D0/
      data     gsp8( 30)/         3.436116D0/
      data     gpp8( 30)/        20.000041D0/
      data     gp28( 30)/         6.782785D0/
      data     hsp8( 30)/         0.662036D0/
!
!                    Data for Element  31          Gallium
!
      data     uss8( 31)/       -30.600226D0/
      data     upp8( 31)/       -21.032425D0/
      data   betas8( 31)/       -10.808320D0/
      data   betap8( 31)/        -4.185500D0/
      data      zs8( 31)/         2.339067D0/
      data      zp8( 31)/         1.729592D0/
      data     gss8( 31)/        10.354885D0/
      data     gsp8( 31)/         7.993674D0/
      data     gpp8( 31)/         6.090184D0/
      data     gp28( 31)/         6.299226D0/
      data     hsp8( 31)/         1.295974D0/
!
!                    Data for Element  32        Germanium
!
      data     uss8( 32)/       -32.747338D0/
      data     upp8( 32)/       -24.709016D0/
      data   betas8( 32)/       -14.854297D0/
      data   betap8( 32)/        -2.591260D0/
      data      zs8( 32)/         2.546073D0/
      data      zp8( 32)/         1.709130D0/
      data     gss8( 32)/         7.518301D0/
      data     gsp8( 32)/         6.594443D0/
      data     gpp8( 32)/         6.066801D0/
      data     gp28( 32)/         5.305947D0/
      data     hsp8( 32)/         0.290742D0/
!
!                    Data for Element  33          Arsenic
!
      data     uss8( 33)/       -37.956965D0/
      data     upp8( 33)/       -38.453701D0/
      data     udd8( 33)/       -30.282658D0/
      data   betas8( 33)/       -11.963725D0/
      data   betap8( 33)/        -7.340073D0/
      data   betad8( 33)/         3.753005D0/
      data      zs8( 33)/         2.926171D0/
      data      zp8( 33)/         1.765191D0/
      data      zd8( 33)/         1.392142D0/
      data     zsn8( 33)/         2.006543D0/
      data     zpn8( 33)/         3.316832D0/
      data     zdn8( 33)/         4.653440D0/
      data     gss8( 33)/         6.665030D0/
      data     gsp8( 33)/         6.213867D0/
      data     gpp8( 33)/         9.310836D0/
      data     gp28( 33)/         8.712542D0/
      data     hsp8( 33)/         0.280662D0/
!
!                    Data for Element  34         Selenium
!
      data     uss8( 34)/       -32.671088D0/
      data     upp8( 34)/       -32.522220D0/
      data   betas8( 34)/         2.636001D0/
      data   betap8( 34)/        -9.557700D0/
      data      zs8( 34)/         2.512366D0/
      data      zp8( 34)/         2.007576D0/
      data     gss8( 34)/         5.522356D0/
      data     gsp8( 34)/         2.907562D0/
      data     gpp8( 34)/         8.042391D0/
      data     gp28( 34)/         6.735106D0/
      data     hsp8( 34)/         3.095789D0/
!
!                    Data for Element  35          Bromine
!
      data     uss8( 35)/       -45.834364D0/
      data     upp8( 35)/       -50.293675D0/
      data     udd8( 35)/         7.086738D0/
      data   betas8( 35)/       -32.131665D0/
      data   betap8( 35)/        -9.514484D0/
      data   betad8( 35)/        -9.839124D0/
      data      zs8( 35)/         4.670684D0/
      data      zp8( 35)/         2.035626D0/
      data      zd8( 35)/         1.521031D0/
      data     zsn8( 35)/         3.094777D0/
      data     zpn8( 35)/         3.065764D0/
      data     zdn8( 35)/         2.820003D0/
      data     gss8( 35)/         7.616791D0/
      data     gsp8( 35)/         5.010425D0/
      data     gpp8( 35)/         9.649216D0/
      data     gp28( 35)/         8.343792D0/
      data     hsp8( 35)/         4.996553D0/
      data   polvo8( 35)/         2.142420D0/
      data gues81( 35,1)/        -0.004996D0/
      data gues82( 35,1)/         6.001292D0/
      data gues83( 35,1)/         2.895153D0/
!
!                    Data for Element  36          Krypton
!
      data     uss8( 36)/         8.535384D0/
      data     upp8( 36)/       -80.484321D0/
      data   betas8( 36)/        -2.727088D0/
      data   betap8( 36)/       -16.142951D0/
      data      zs8( 36)/         1.312248D0/
      data      zp8( 36)/         4.491371D0/
      data     gss8( 36)/        19.999857D0/
      data     gsp8( 36)/         1.175304D0/
      data     gpp8( 36)/         9.174784D0/
      data     gp28( 36)/        14.926948D0/
      data     hsp8( 36)/         0.299867D0/
!
!                    Data for Element  37         Rubidium
!
      data     uss8( 37)/        -3.636505D0/
      data     upp8( 37)/        -2.500671D0/
      data   betas8( 37)/         9.998744D0/
      data   betap8( 37)/         1.343004D0/
      data      zs8( 37)/         5.510145D0/
      data      zp8( 37)/         1.335170D0/
      data     gss8( 37)/         6.680824D0/
      data     gsp8( 37)/        20.001098D0/
      data     gpp8( 37)/         5.068874D0/
      data     gp28( 37)/         2.747860D0/
      data     hsp8( 37)/         3.602834D0/
!
!                    Data for Element  38        Strontium
!
      data     uss8( 38)/       -10.427671D0/
      data     upp8( 38)/        -9.943751D0/
      data   betas8( 38)/        -6.253108D0/
      data   betap8( 38)/        -9.844498D0/
      data      zs8( 38)/         2.197303D0/
      data      zp8( 38)/         1.730137D0/
      data     gss8( 38)/         4.603664D0/
      data     gsp8( 38)/         5.716069D0/
      data     gpp8( 38)/         7.334620D0/
      data     gp28( 38)/         7.443088D0/
      data     hsp8( 38)/         0.831527D0/
      data gues81( 38,1)/        -0.012948D0/
      data gues82( 38,1)/         6.000126D0/
      data gues83( 38,1)/         3.011964D0/
!
!                    Data for Element  39          Yttrium
!
      data     uss8( 39)/       -14.247809D0/
      data     upp8( 39)/       -14.817140D0/
      data     udd8( 39)/       -16.394302D0/
      data   betas8( 39)/         0.343336D0/
      data   betap8( 39)/        -3.180807D0/
      data   betad8( 39)/        -4.508957D0/
      data      zs8( 39)/         0.593368D0/
      data      zp8( 39)/         1.490422D0/
      data      zd8( 39)/         1.650893D0/
      data     zsn8( 39)/         0.902611D0/
      data     zpn8( 39)/         1.484400D0/
      data     zdn8( 39)/         1.384238D0/
      data     alp8( 39)/         0.500727D0/
      data     gss8( 39)/         4.046733D0/
      data     gsp8( 39)/         4.726277D0/
      data     gpp8( 39)/         7.278752D0/
      data     gp28( 39)/         6.343281D0/
      data     hsp8( 39)/         0.679228D0/
      data    poc_8( 39)/         2.773703D0/
      data    f0sd8( 39)/         4.972716D0/
      data    g2sd8( 39)/         5.016364D0/
!
!                    Data for Element  40        Zirconium
!
      data     uss8( 40)/       -20.008884D0/
      data     upp8( 40)/       -14.559692D0/
      data     udd8( 40)/       -21.302657D0/
      data   betas8( 40)/         9.551952D0/
      data   betap8( 40)/        -4.551915D0/
      data   betad8( 40)/        -3.213274D0/
      data      zs8( 40)/         1.692590D0/
      data      zp8( 40)/         1.694916D0/
      data      zd8( 40)/         1.567392D0/
      data     zsn8( 40)/         1.189109D0/
      data     zpn8( 40)/         0.809092D0/
      data     zdn8( 40)/         1.190249D0/
      data     gss8( 40)/         5.331208D0/
      data     gsp8( 40)/         4.150579D0/
      data     gpp8( 40)/         3.967381D0/
      data     gp28( 40)/         3.457490D0/
      data     hsp8( 40)/         0.743676D0/
      data    f0sd8( 40)/         5.010704D0/
      data    g2sd8( 40)/         2.943652D0/
!
!                    Data for Element  41          Niobium
!
      data     uss8( 41)/       -31.269298D0/
      data     upp8( 41)/       -20.151277D0/
      data     udd8( 41)/       -35.893116D0/
      data   betas8( 41)/       -12.045244D0/
      data   betap8( 41)/         1.465762D0/
      data   betad8( 41)/        -5.920160D0/
      data      zs8( 41)/         2.355562D0/
      data      zp8( 41)/         1.386907D0/
      data      zd8( 41)/         1.977324D0/
      data     zsn8( 41)/         1.490754D0/
      data     zpn8( 41)/         0.892760D0/
      data     zdn8( 41)/         1.443837D0/
      data     alp8( 41)/         0.843974D0/
      data     gss8( 41)/         6.683592D0/
      data     gsp8( 41)/         4.685339D0/
      data     gpp8( 41)/         4.377647D0/
      data     gp28( 41)/         3.815028D0/
      data     hsp8( 41)/         0.650679D0/
      data    f0sd8( 41)/         6.550674D0/
      data    g2sd8( 41)/         1.065577D0/
!
!                    Data for Element  42       Molybdenum
!
      data     uss8( 42)/       -53.467728D0/
      data     upp8( 42)/       -35.291951D0/
      data     udd8( 42)/       -55.836977D0/
      data   betas8( 42)/        -0.189344D0/
      data   betap8( 42)/         7.017762D0/
      data   betad8( 42)/       -10.941126D0/
      data      zs8( 42)/         1.060429D0/
      data      zp8( 42)/         1.350412D0/
      data      zd8( 42)/         1.827152D0/
      data     zsn8( 42)/         1.912995D0/
      data     zpn8( 42)/         1.355055D0/
      data     zdn8( 42)/         1.876231D0/
      data     gss8( 42)/         8.576652D0/
      data     gsp8( 42)/         6.888293D0/
      data     gpp8( 42)/         6.644509D0/
      data     gp28( 42)/         5.790552D0/
      data     hsp8( 42)/         1.317368D0/
      data    f0sd8( 42)/        10.000608D0/
      data    g2sd8( 42)/         1.216752D0/
!
!                    Data for Element  43       Technetium
!
      data     uss8( 43)/       -41.850292D0/
      data     upp8( 43)/       -34.910293D0/
      data     udd8( 43)/       -45.530412D0/
      data   betas8( 43)/        -2.791024D0/
      data   betap8( 43)/        -8.086697D0/
      data   betad8( 43)/        -5.724335D0/
      data      zs8( 43)/         1.956245D0/
      data      zp8( 43)/         6.006299D0/
      data      zd8( 43)/         1.767360D0/
      data     zsn8( 43)/         1.411033D0/
      data     zpn8( 43)/         1.141313D0/
      data     zdn8( 43)/         1.159312D0/
      data     gss8( 43)/         6.326174D0/
      data     gsp8( 43)/         5.587138D0/
      data     gpp8( 43)/         5.596426D0/
      data     gp28( 43)/         4.877169D0/
      data     hsp8( 43)/         1.258989D0/
      data    f0sd8( 43)/         5.434886D0/
      data    g2sd8( 43)/         1.106875D0/
!
!                    Data for Element  44        Ruthenium
!
      data     uss8( 44)/       -44.901521D0/
      data     upp8( 44)/       -41.424409D0/
      data     udd8( 44)/       -37.934514D0/
      data   betas8( 44)/       -12.859508D0/
      data   betap8( 44)/        -8.475518D0/
      data   betad8( 44)/        -3.830797D0/
      data      zs8( 44)/         1.459195D0/
      data      zp8( 44)/         5.537201D0/
      data      zd8( 44)/         2.093164D0/
      data     zsn8( 44)/         0.984449D0/
      data     zpn8( 44)/         4.586613D0/
      data     zdn8( 44)/         0.765332D0/
      data     gss8( 44)/         4.413643D0/
      data     gsp8( 44)/         5.356996D0/
      data     gpp8( 44)/        22.490448D0/
      data     gp28( 44)/        19.599957D0/
      data     hsp8( 44)/         0.008058D0/
      data    f0sd8( 44)/         5.917404D0/
      data    g2sd8( 44)/         5.859738D0/
!
!                    Data for Element  45          Rhodium
!
      data     uss8( 45)/       -20.513756D0/
      data     upp8( 45)/       -40.045431D0/
      data     udd8( 45)/       -35.818492D0/
      data   betas8( 45)/        -8.222141D0/
      data   betap8( 45)/       -15.556691D0/
      data   betad8( 45)/       -13.396182D0/
      data      zs8( 45)/         1.324919D0/
      data      zp8( 45)/         4.306111D0/
      data      zd8( 45)/         2.901406D0/
      data     zsn8( 45)/         0.809923D0/
      data     zpn8( 45)/         6.898259D0/
      data     zdn8( 45)/         0.643134D0/
      data     gss8( 45)/         3.631179D0/
      data     gsp8( 45)/         4.407820D0/
      data     gpp8( 45)/        33.825599D0/
      data     gp28( 45)/        29.478305D0/
      data     hsp8( 45)/         0.000092D0/
      data    f0sd8( 45)/         1.775497D0/
      data    g2sd8( 45)/         1.851571D0/
!
!                    Data for Element  46        Palladium
!
      data     uss8( 46)/       -76.140196D0/
      data     upp8( 46)/       -21.073362D0/
      data     udd8( 46)/       -85.325301D0/
      data   betas8( 46)/        -8.038245D0/
      data   betap8( 46)/         0.740037D0/
      data   betad8( 46)/        -2.394498D0/
      data      zs8( 46)/         1.658503D0/
      data      zp8( 46)/         1.156718D0/
      data      zd8( 46)/         2.219861D0/
      data     zsn8( 46)/         1.794085D0/
      data     zpn8( 46)/         6.158778D0/
      data     zdn8( 46)/         1.630913D0/
      data     gss8( 46)/         8.043535D0/
      data     gsp8( 46)/         9.755042D0/
      data     gpp8( 46)/        30.199556D0/
      data     gp28( 46)/        26.318284D0/
      data     hsp8( 46)/         0.086121D0/
      data    f0sd8( 46)/         8.004447D0/
      data    g2sd8( 46)/         2.613148D0/
!
!                    Data for Element  47           Silver
!
      data     uss8( 47)/       -25.484137D0/
      data     upp8( 47)/       -36.116023D0/
      data     udd8( 47)/       -35.668272D0/
      data   betas8( 47)/        -6.129623D0/
      data   betap8( 47)/         1.004115D0/
      data   betad8( 47)/       -69.238347D0/
      data      zs8( 47)/         1.994004D0/
      data      zp8( 47)/         0.681817D0/
      data      zd8( 47)/         6.007328D0/
      data     zsn8( 47)/         0.695514D0/
      data     zpn8( 47)/         4.729949D0/
      data     zdn8( 47)/         0.506522D0/
      data     gss8( 47)/         3.118242D0/
      data     gsp8( 47)/         3.785152D0/
      data     gpp8( 47)/        23.193295D0/
      data     gp28( 47)/        20.212474D0/
      data     hsp8( 47)/         0.000432D0/
      data    f0sd8( 47)/         1.938327D0/
      data    g2sd8( 47)/         1.071901D0/
!
!                    Data for Element  48          Cadmium
!
      data     uss8( 48)/       -14.645792D0/
      data     upp8( 48)/        -9.318664D0/
      data   betas8( 48)/       -11.613183D0/
      data   betap8( 48)/         1.663178D0/
      data      zs8( 48)/         1.384108D0/
      data      zp8( 48)/         1.957413D0/
      data     gss8( 48)/         6.677284D0/
      data     gsp8( 48)/         5.953373D0/
      data     gpp8( 48)/        18.729843D0/
      data     gp28( 48)/         9.917452D0/
      data     hsp8( 48)/         0.825192D0/
!
!                    Data for Element  49           Indium
!
      data     uss8( 49)/       -28.339246D0/
      data     upp8( 49)/       -23.373875D0/
      data   betas8( 49)/        -1.982376D0/
      data   betap8( 49)/        -3.330294D0/
      data      zs8( 49)/         2.023087D0/
      data      zp8( 49)/         2.106618D0/
      data     gss8( 49)/         9.906091D0/
      data     gsp8( 49)/        10.520060D0/
      data     gpp8( 49)/         4.826006D0/
      data     gp28( 49)/         7.906563D0/
      data     hsp8( 49)/         3.500299D0/
!
!                    Data for Element  50              Tin
!
      data     uss8( 50)/       -29.888217D0/
      data     upp8( 50)/       -22.156954D0/
      data   betas8( 50)/        -8.621087D0/
      data   betap8( 50)/        -4.989752D0/
      data      zs8( 50)/         2.383941D0/
      data      zp8( 50)/         2.057908D0/
      data     gss8( 50)/         8.269655D0/
      data     gsp8( 50)/         5.013349D0/
      data     gpp8( 50)/         6.584874D0/
      data     gp28( 50)/         5.855159D0/
      data     hsp8( 50)/         0.531212D0/
      data gues81( 50,1)/        -1.004587D0/
      data gues82( 50,1)/         4.706252D0/
      data gues83( 50,1)/         1.180218D0/
!
!                    Data for Element  51         Antimony
!
      data     uss8( 51)/       -41.688879D0/
      data     upp8( 51)/       -39.541180D0/
      data     udd8( 51)/        -6.581663D0/
      data   betas8( 51)/        -7.472322D0/
      data   betap8( 51)/        -5.940750D0/
      data   betad8( 51)/        -3.979108D0/
      data      zs8( 51)/         2.391178D0/
      data      zp8( 51)/         1.773006D0/
      data      zd8( 51)/         2.465590D0/
      data     zsn8( 51)/         5.993591D0/
      data     zpn8( 51)/         6.145086D0/
      data     zdn8( 51)/         5.704031D0/
      data     gss8( 51)/        10.588832D0/
      data     gsp8( 51)/         7.310023D0/
      data     gpp8( 51)/         9.281609D0/
      data     gp28( 51)/         8.954081D0/
      data     hsp8( 51)/         0.779112D0/
!
!                    Data for Element  52        Tellurium
!
      data     uss8( 52)/      -114.733316D0/
      data     upp8( 52)/       -50.096389D0/
      data   betas8( 52)/       -70.001062D0/
      data   betap8( 52)/        -6.151642D0/
      data      zs8( 52)/         2.769862D0/
      data      zp8( 52)/         1.731319D0/
      data     gss8( 52)/         7.030626D0/
      data     gsp8( 52)/        12.601389D0/
      data     gpp8( 52)/         7.883479D0/
      data     gp28( 52)/         6.973163D0/
      data     hsp8( 52)/         5.000826D0/
!
!                    Data for Element  53           Iodine
!
      data     uss8( 53)/       -59.973232D0/
      data     upp8( 53)/       -56.459835D0/
      data     udd8( 53)/       -28.822603D0/
      data   betas8( 53)/       -30.522481D0/
      data   betap8( 53)/        -5.942120D0/
      data   betad8( 53)/        -7.676107D0/
      data      zs8( 53)/         4.498653D0/
      data      zp8( 53)/         1.917072D0/
      data      zd8( 53)/         1.875175D0/
      data     zsn8( 53)/         9.135244D0/
      data     zpn8( 53)/         6.888191D0/
      data     zdn8( 53)/         3.791523D0/
      data     gss8( 53)/         7.234759D0/
      data     gsp8( 53)/         9.154406D0/
      data     gpp8( 53)/         9.877466D0/
      data     gp28( 53)/         8.035916D0/
      data     hsp8( 53)/         5.004215D0/
      data   polvo8( 53)/         3.823160D0/
      data gues81( 53,1)/        -0.035519D0/
      data gues82( 53,1)/         1.744389D0/
      data gues83( 53,1)/         1.223844D0/
!
!                    Data for Element  54            Xenon
!
      data     uss8( 54)/       -18.270227D0/
      data     upp8( 54)/      -167.163063D0/
      data   betas8( 54)/        -3.980622D0/
      data   betap8( 54)/       -38.822792D0/
      data      zs8( 54)/         2.759787D0/
      data      zp8( 54)/         1.977446D0/
      data     gss8( 54)/        20.000252D0/
      data     gsp8( 54)/         4.175902D0/
      data     gpp8( 54)/         2.305787D0/
      data     gp28( 54)/         4.063220D0/
      data     hsp8( 54)/         4.418843D0/
!
!                    Data for Element  55           Cesium
!
      data     uss8( 55)/        -3.748609D0/
      data     upp8( 55)/        -2.348109D0/
      data   betas8( 55)/         2.287838D0/
      data   betap8( 55)/        -5.908071D0/
      data      zs8( 55)/         5.956008D0/
      data      zp8( 55)/         1.619485D0/
      data     gss8( 55)/         6.464751D0/
      data     gsp8( 55)/         4.004501D0/
      data     gpp8( 55)/        13.775390D0/
      data     gp28( 55)/        12.912537D0/
      data     hsp8( 55)/         1.026928D0/
!
!                    Data for Element  56           Barium
!
      data     uss8( 56)/        -9.306985D0/
      data     upp8( 56)/        -8.826713D0/
      data   betas8( 56)/        10.003125D0/
      data   betap8( 56)/        -6.335160D0/
      data      zs8( 56)/         1.395379D0/
      data      zp8( 56)/         1.430139D0/
      data     gss8( 56)/         3.600823D0/
      data     gsp8( 56)/         4.740579D0/
      data     gpp8( 56)/         3.345166D0/
      data     gp28( 56)/         3.142783D0/
      data     hsp8( 56)/         0.929429D0/
!
!                    Data for Element  57        Lanthanum
!
      data     uss8( 57)/       -19.641953D0/
      data     upp8( 57)/       -22.059431D0/
      data     udd8( 57)/       -22.638986D0/
      data   betas8( 57)/         0.796727D0/
      data   betap8( 57)/       -10.856056D0/
      data   betad8( 57)/        -0.484922D0/
      data      zs8( 57)/         2.673780D0/
      data      zp8( 57)/         1.248192D0/
      data      zd8( 57)/         1.688562D0/
      data     zsn8( 57)/         1.617784D0/
      data     zpn8( 57)/         4.331620D0/
      data     zdn8( 57)/         2.285738D0/
      data     alp8( 57)/         5.940443D0/
      data     gss8( 57)/         6.154440D0/
      data     gsp8( 57)/         7.322704D0/
      data     gpp8( 57)/        18.077465D0/
      data     gp28( 57)/        15.679057D0/
      data     hsp8( 57)/         0.138601D0/
      data    poc_8( 57)/         2.511701D0/
      data    f0sd8( 57)/         8.856858D0/
      data    g2sd8( 57)/         7.925585D0/
!
!                    Data for Element  71         Lutetium
!
      data     uss8( 71)/       -15.954994D0/
      data     upp8( 71)/       -11.606213D0/
      data     udd8( 71)/       -13.050056D0/
      data   betas8( 71)/        -5.590778D0/
      data   betap8( 71)/        -0.937679D0/
      data   betad8( 71)/        -7.737752D0/
      data      zs8( 71)/         5.471741D0/
      data      zp8( 71)/         1.712296D0/
      data      zd8( 71)/         2.225892D0/
      data     zsn8( 71)/         1.632335D0/
      data     zpn8( 71)/         4.033128D0/
      data     zdn8( 71)/         0.921999D0/
      data     gss8( 71)/         6.209796D0/
      data     gsp8( 71)/         7.379102D0/
      data     gpp8( 71)/        16.831746D0/
      data     gp28( 71)/        14.598613D0/
      data     hsp8( 71)/         0.209008D0/
      data    poc_8( 71)/         2.743262D0/
      data    f0sd8( 71)/         3.924927D0/
      data    g2sd8( 71)/         1.000946D0/
!
!                    Data for Element  72          Hafnium
!
      data     uss8( 72)/       -22.375140D0/
      data     upp8( 72)/       -13.081670D0/
      data     udd8( 72)/       -20.637741D0/
      data   betas8( 72)/        -5.366351D0/
      data   betap8( 72)/       -21.550119D0/
      data   betad8( 72)/        -3.884443D0/
      data      zs8( 72)/         3.085344D0/
      data      zp8( 72)/         1.575819D0/
      data      zd8( 72)/         1.840840D0/
      data     zsn8( 72)/         0.946927D0/
      data     zpn8( 72)/         3.538911D0/
      data     zdn8( 72)/         0.940283D0/
      data     gss8( 72)/         3.602338D0/
      data     gsp8( 72)/         4.293729D0/
      data     gpp8( 72)/        14.769194D0/
      data     gp28( 72)/        12.809708D0/
      data     hsp8( 72)/         0.011028D0/
      data    f0sd8( 72)/         4.842900D0/
      data    g2sd8( 72)/         4.386101D0/
!
!                    Data for Element  73         Tantalum
!
      data     uss8( 73)/       -39.009984D0/
      data     upp8( 73)/         1.163975D0/
      data     udd8( 73)/       -43.266315D0/
      data   betas8( 73)/       -17.199605D0/
      data   betap8( 73)/        -5.818839D0/
      data   betad8( 73)/        -9.816794D0/
      data      zs8( 73)/         4.578087D0/
      data      zp8( 73)/         4.841244D0/
      data      zd8( 73)/         1.838249D0/
      data     zsn8( 73)/         1.741367D0/
      data     zpn8( 73)/         3.430157D0/
      data     zdn8( 73)/         2.311198D0/
      data     gss8( 73)/         6.624580D0/
      data     gsp8( 73)/         7.805321D0/
      data     gpp8( 73)/        14.315323D0/
      data     gp28( 73)/        12.416054D0/
      data     hsp8( 73)/         0.577263D0/
      data    f0sd8( 73)/         8.544427D0/
      data    g2sd8( 73)/         2.074254D0/
!
!                    Data for Element  74         Tungsten
!
      data     uss8( 74)/       -44.524950D0/
      data     upp8( 74)/       -40.011500D0/
      data     udd8( 74)/       -46.490410D0/
      data   betas8( 74)/       -16.946460D0/
      data   betap8( 74)/         5.623170D0/
      data   betad8( 74)/        -2.947340D0/
      data      zs8( 74)/         2.664560D0/
      data      zp8( 74)/         1.624010D0/
      data      zd8( 74)/         1.794400D0/
      data     zsn8( 74)/         1.498860D0/
      data     zpn8( 74)/         1.965900D0/
      data     zdn8( 74)/         1.876450D0/
      data     gss8( 74)/         5.702025D0/
      data     gsp8( 74)/         6.323145D0/
      data     gpp8( 74)/         8.204433D0/
      data     gp28( 74)/         7.115919D0/
      data     hsp8( 74)/         1.319912D0/
      data    f0sd8( 74)/         7.788180D0/
      data    g2sd8( 74)/         1.684940D0/
!
!                    Data for Element  75          Rhenium
!
      data     uss8( 75)/       -41.291342D0/
      data     upp8( 75)/       -35.089592D0/
      data     udd8( 75)/       -44.178985D0/
      data   betas8( 75)/         3.830075D0/
      data   betap8( 75)/        -1.638530D0/
      data   betad8( 75)/        -1.414411D0/
      data      zs8( 75)/         2.411839D0/
      data      zp8( 75)/         1.815351D0/
      data      zd8( 75)/         2.522766D0/
      data     zsn8( 75)/         1.680823D0/
      data     zpn8( 75)/         1.331218D0/
      data     zdn8( 75)/         1.490623D0/
      data     gss8( 75)/         6.394256D0/
      data     gsp8( 75)/         5.555571D0/
      data     gpp8( 75)/         5.555669D0/
      data     gp28( 75)/         4.818577D0/
      data     hsp8( 75)/         1.220913D0/
      data    f0sd8( 75)/         5.442818D0/
      data    g2sd8( 75)/         2.376279D0/
!
!                    Data for Element  76           Osmium
!
      data     uss8( 76)/       -26.434080D0/
      data     upp8( 76)/       -48.739500D0/
      data     udd8( 76)/       -55.837880D0/
      data   betas8( 76)/       -12.508730D0/
      data   betap8( 76)/         0.846880D0/
      data   betad8( 76)/         5.164360D0/
      data      zs8( 76)/         3.031000D0/
      data      zp8( 76)/         1.593960D0/
      data      zd8( 76)/         1.775570D0/
      data     zsn8( 76)/         1.844700D0/
      data     zpn8( 76)/         1.564220D0/
      data     zdn8( 76)/         1.770010D0/
      data     gss8( 76)/         7.017683D0/
      data     gsp8( 76)/         6.384200D0/
      data     gpp8( 76)/         6.528073D0/
      data     gp28( 76)/         5.661968D0/
      data     hsp8( 76)/         1.508926D0/
      data    f0sd8( 76)/         2.021170D0/
      data    g2sd8( 76)/         1.392130D0/
!
!                    Data for Element  77          Iridium
!
      data     uss8( 77)/       -29.703974D0/
      data     upp8( 77)/       -38.210924D0/
      data     udd8( 77)/       -32.538202D0/
      data   betas8( 77)/       -10.943427D0/
      data   betap8( 77)/         2.908880D0/
      data   betad8( 77)/        -3.791731D0/
      data      zs8( 77)/         1.500907D0/
      data      zp8( 77)/         4.106373D0/
      data      zd8( 77)/         2.676047D0/
      data     zsn8( 77)/         0.927246D0/
      data     zpn8( 77)/         3.191892D0/
      data     zdn8( 77)/         0.662007D0/
      data     gss8( 77)/         3.527467D0/
      data     gsp8( 77)/         4.203820D0/
      data     gpp8( 77)/        13.320955D0/
      data     gp28( 77)/        11.553612D0/
      data     hsp8( 77)/         0.018501D0/
      data    f0sd8( 77)/         2.627170D0/
      data    g2sd8( 77)/         2.996029D0/
!
!                    Data for Element  78         Platinum
!
      data     uss8( 78)/       -73.516173D0/
      data     upp8( 78)/       -68.320056D0/
      data     udd8( 78)/       -76.598873D0/
      data   betas8( 78)/         1.151418D0/
      data   betap8( 78)/         3.298694D0/
      data   betad8( 78)/       -18.044737D0/
      data      zs8( 78)/         2.301264D0/
      data      zp8( 78)/         1.662404D0/
      data      zd8( 78)/         3.168852D0/
      data     zsn8( 78)/         2.270699D0/
      data     zpn8( 78)/         1.949896D0/
      data     zdn8( 78)/         1.713856D0/
      data     gss8( 78)/         8.638286D0/
      data     gsp8( 78)/         7.922254D0/
      data     gpp8( 78)/         8.137643D0/
      data     gp28( 78)/         7.057990D0/
      data     hsp8( 78)/         1.892617D0/
      data    f0sd8( 78)/         7.098591D0/
      data    g2sd8( 78)/         4.484183D0/
!
!                    Data for Element  79             Gold
!
      data     uss8( 79)/       -95.041846D0/
      data     upp8( 79)/       -63.890158D0/
      data     udd8( 79)/       -88.066087D0/
      data   betas8( 79)/        -7.479625D0/
      data   betap8( 79)/         3.664356D0/
      data   betad8( 79)/       -61.715468D0/
      data      zs8( 79)/         1.814169D0/
      data      zp8( 79)/         1.618657D0/
      data      zd8( 79)/         5.053167D0/
      data     zsn8( 79)/         2.444680D0/
      data     zpn8( 79)/         7.014990D0/
      data     zdn8( 79)/         1.777089D0/
      data     gss8( 79)/         9.300152D0/
      data     gsp8( 79)/        11.073443D0/
      data     gpp8( 79)/        29.276168D0/
      data     gp28( 79)/        25.391984D0/
      data     hsp8( 79)/         0.144384D0/
      data    f0sd8( 79)/         8.827257D0/
      data    g2sd8( 79)/         4.915625D0/
!
!                    Data for Element  80          Mercury
!
      data     uss8( 80)/       -17.608732D0/
      data     upp8( 80)/       -18.369417D0/
      data   betas8( 80)/        -3.045239D0/
      data   betap8( 80)/        -5.693556D0/
      data      zs8( 80)/         2.104896D0/
      data      zp8( 80)/         1.516293D0/
      data     gss8( 80)/         6.372822D0/
      data     gsp8( 80)/        10.143176D0/
      data     gpp8( 80)/        10.397393D0/
      data     gp28( 80)/        14.794056D0/
      data     hsp8( 80)/         0.926128D0/
!
!                    Data for Element  81         Thallium
!
      data     uss8( 81)/       -29.518621D0/
      data     upp8( 81)/       -29.826907D0/
      data   betas8( 81)/        -7.230170D0/
      data   betap8( 81)/        -7.575544D0/
      data      zs8( 81)/         3.335883D0/
      data      zp8( 81)/         1.766141D0/
      data     gss8( 81)/         5.015118D0/
      data     gsp8( 81)/        13.932049D0/
      data     gpp8( 81)/        10.495551D0/
      data     gp28( 81)/        10.526198D0/
      data     hsp8( 81)/         0.293760D0/
!
!                    Data for Element  82             Lead
!
      data     uss8( 82)/       -35.038145D0/
      data     upp8( 82)/       -25.413401D0/
      data   betas8( 82)/        -8.323792D0/
      data   betap8( 82)/        -2.237891D0/
      data      zs8( 82)/         2.368901D0/
      data      zp8( 82)/         1.685246D0/
      data     gss8( 82)/         5.254128D0/
      data     gsp8( 82)/         7.061016D0/
      data     gpp8( 82)/         6.818551D0/
      data     gp28( 82)/         5.603019D0/
      data     hsp8( 82)/         1.018819D0/
      data gues81( 82,1)/        -0.239463D0/
      data gues82( 82,1)/         5.444338D0/
      data gues83( 82,1)/         1.613682D0/
!
!                    Data for Element  83          Bismuth
!
      data     uss8( 83)/       -42.409177D0/
      data     upp8( 83)/       -36.393746D0/
      data   betas8( 83)/       -34.951578D0/
      data   betap8( 83)/        -7.359060D0/
      data      zs8( 83)/         3.702377D0/
      data      zp8( 83)/         1.872327D0/
      data     gss8( 83)/         5.851803D0/
      data     gsp8( 83)/         6.790583D0/
      data     gpp8( 83)/         8.389442D0/
      data     gp28( 83)/         7.724219D0/
      data     hsp8( 83)/         0.295606D0/
!
!                    Data for Element  85         Astatine
!
      data     alp8( 85)/         3.000000D0/
      data     gss8( 85)/        10.000000D0/
!
!                    Data for Element  87         Francium
!
      data     alp8( 87)/         3.000000D0/
      data     gss8( 87)/        10.000000D0/
!
!                    Data for Element  90          Thorium
!
      data     uss8( 90)/       -40.568292D0/
      data     upp8( 90)/       -28.089187D0/
      data   betas8( 90)/        -4.256218D0/
      data   betap8( 90)/        -4.256218D0/
      data      zs8( 90)/         1.435306D0/
      data      zp8( 90)/         1.435306D0/
      data     gss8( 90)/         9.820000D0/
      data     gsp8( 90)/         8.360000D0/
      data     gpp8( 90)/         7.310000D0/
      data     gp28( 90)/         6.540000D0/
      data     hsp8( 90)/         1.320000D0/
!
!                    Data for Element  97        Berkelium
!
      data gues81( 97,1)/         1.480000D0/
      data gues82( 97,1)/         0.960000D0/
      data gues81( 97,2)/         1.560000D0/
      data gues82( 97,2)/         0.760000D0/
      data gues81( 97,3)/         1.550000D0/
      data gues82( 97,3)/         0.850000D0/
!
!                    Data for Element  98          Mithril
!
      data     uss8( 98)/        -3.000000D0/
      data   betas8( 98)/       -99.000000D0/
      data      zs8( 98)/         2.000000D0/
      data     gss8( 98)/        12.000000D0/
!
!                    Data for Element 100       3+ Sparkle
!
      data     alp8(100)/         1.500000D0/
!
!                    Data for Element 101       3- Sparkle
!
      data     alp8(101)/         1.500000D0/
!
!                    Data for Element 102      Capped bond
!
      data   betas8(102)/  -9999999.000000D0/
      data      zs8(102)/         4.000000D0/
      data     gss8(102)/        12.848000D0/
!
!                    Data for Element 103       ++ Sparkle
!
      data     alp8(103)/         1.500000D0/
!
!                    Data for Element 104        + Sparkle
!
      data     alp8(104)/         1.500000D0/
!
!                    Data for Element 105       -- Sparkle
!
      data     alp8(105)/         1.500000D0/
!
!                    Data for Element 106        - Sparkle
!
      data     alp8(106)/         1.500000D0/
!
!
!                     Global parameters
!
!
      data   v_par8(1)/    9.278465d0/  ! Used in ccrep for scalar correction of C-C triple bonds.
      data   v_par8(2)/    5.983752d0/  ! Used in ccrep for exponent correction of C-C triple bonds.
      data   v_par8(7)/    1.000000d0/  ! Used in dftd3 to set "s6"  in D3H4
      data   v_par8(8)/   14.000000d0/  ! Used in dftd3 to set "alp" in D3H4
      data   v_par8(9)/    1.560000d0/  ! Used in dftd3 to set "rs6" in D3H4
  contains
  subroutine alpb_and_xfac_pm8
    use parameters_C, only : xfac, alpb
 !
      alpb( 1, 1) =     3.540942d0 !    Hydrogen -     Hydrogen
      xfac( 1, 1) =     2.243587d0 !    Hydrogen -     Hydrogen
 !
      alpb( 2, 1) =     2.989881d0 !      Helium -     Hydrogen
      xfac( 2, 1) =     2.371199d0 !      Helium -     Hydrogen
      alpb( 2, 2) =     3.783559d0 !      Helium -       Helium
      xfac( 2, 2) =     3.450900d0 !      Helium -       Helium
 !
      alpb( 3, 1) =     2.136265d0 !     Lithium -     Hydrogen
      xfac( 3, 1) =     2.191985d0 !     Lithium -     Hydrogen
      alpb( 3, 2) =     3.112403d0 !     Lithium -       Helium
      xfac( 3, 2) =     9.273676d0 !     Lithium -       Helium
      alpb( 3, 3) =     4.714674d0 !     Lithium -      Lithium
      xfac( 3, 3) =    16.116384d0 !     Lithium -      Lithium
 !
      alpb( 4, 1) =     2.475418d0 !   Beryllium -     Hydrogen
      xfac( 4, 1) =     2.562831d0 !   Beryllium -     Hydrogen
      alpb( 4, 2) =     3.306702d0 !   Beryllium -       Helium
      xfac( 4, 2) =    12.544878d0 !   Beryllium -       Helium
      alpb( 4, 3) =     2.236728d0 !   Beryllium -      Lithium
      xfac( 4, 3) =     3.287165d0 !   Beryllium -      Lithium
      alpb( 4, 4) =     1.499907d0 !   Beryllium -    Beryllium
      xfac( 4, 4) =     0.238633d0 !   Beryllium -    Beryllium
 !
      alpb( 5, 1) =     2.615231d0 !       Boron -     Hydrogen
      xfac( 5, 1) =     1.321394d0 !       Boron -     Hydrogen
      alpb( 5, 2) =     3.163140d0 !       Boron -       Helium
      xfac( 5, 2) =     1.974170d0 !       Boron -       Helium
      alpb( 5, 3) =     3.759397d0 !       Boron -      Lithium
      xfac( 5, 3) =     7.886018d0 !       Boron -      Lithium
      alpb( 5, 4) =     1.888998d0 !       Boron -    Beryllium
      xfac( 5, 4) =     1.151792d0 !       Boron -    Beryllium
      alpb( 5, 5) =     3.318624d0 !       Boron -        Boron
      xfac( 5, 5) =     3.593619d0 !       Boron -        Boron
 !
      alpb( 6, 1) =     1.027806d0 !      Carbon -     Hydrogen
      xfac( 6, 1) =     0.216506d0 !      Carbon -     Hydrogen
      alpb( 6, 2) =     3.042705d0 !      Carbon -       Helium
      xfac( 6, 2) =     3.213971d0 !      Carbon -       Helium
      alpb( 6, 3) =     3.241874d0 !      Carbon -      Lithium
      xfac( 6, 3) =    16.180002d0 !      Carbon -      Lithium
      alpb( 6, 4) =     4.212882d0 !      Carbon -    Beryllium
      xfac( 6, 4) =    25.035879d0 !      Carbon -    Beryllium
      alpb( 6, 5) =     2.919007d0 !      Carbon -        Boron
      xfac( 6, 5) =     1.874859d0 !      Carbon -        Boron
      alpb( 6, 6) =     2.613713d0 !      Carbon -       Carbon
      xfac( 6, 6) =     0.813510d0 !      Carbon -       Carbon
 !
      alpb( 7, 1) =     0.969406d0 !    Nitrogen -     Hydrogen
      xfac( 7, 1) =     0.175506d0 !    Nitrogen -     Hydrogen
      alpb( 7, 2) =     2.814339d0 !    Nitrogen -       Helium
      xfac( 7, 2) =     1.077861d0 !    Nitrogen -       Helium
      alpb( 7, 3) =     2.640623d0 !    Nitrogen -      Lithium
      xfac( 7, 3) =     2.823403d0 !    Nitrogen -      Lithium
      alpb( 7, 4) =     2.580895d0 !    Nitrogen -    Beryllium
      xfac( 7, 4) =     1.740605d0 !    Nitrogen -    Beryllium
      alpb( 7, 5) =     2.477004d0 !    Nitrogen -        Boron
      xfac( 7, 5) =     0.952882d0 !    Nitrogen -        Boron
      alpb( 7, 6) =     2.686108d0 !    Nitrogen -       Carbon
      xfac( 7, 6) =     0.859949d0 !    Nitrogen -       Carbon
      alpb( 7, 7) =     2.574502d0 !    Nitrogen -     Nitrogen
      xfac( 7, 7) =     0.675313d0 !    Nitrogen -     Nitrogen
 !
      alpb( 8, 1) =     1.260942d0 !      Oxygen -     Hydrogen
      xfac( 8, 1) =     0.192295d0 !      Oxygen -     Hydrogen
      alpb( 8, 2) =     3.653775d0 !      Oxygen -       Helium
      xfac( 8, 2) =     6.684525d0 !      Oxygen -       Helium
      alpb( 8, 3) =     2.584442d0 !      Oxygen -      Lithium
      xfac( 8, 3) =     1.968598d0 !      Oxygen -      Lithium
      alpb( 8, 4) =     3.051867d0 !      Oxygen -    Beryllium
      xfac( 8, 4) =     3.218155d0 !      Oxygen -    Beryllium
      alpb( 8, 5) =     2.695351d0 !      Oxygen -        Boron
      xfac( 8, 5) =     1.269801d0 !      Oxygen -        Boron
      alpb( 8, 6) =     2.889607d0 !      Oxygen -       Carbon
      xfac( 8, 6) =     0.990211d0 !      Oxygen -       Carbon
      alpb( 8, 7) =     2.784292d0 !      Oxygen -     Nitrogen
      xfac( 8, 7) =     0.764756d0 !      Oxygen -     Nitrogen
      alpb( 8, 8) =     2.623998d0 !      Oxygen -       Oxygen
      xfac( 8, 8) =     0.535112d0 !      Oxygen -       Oxygen
 !
      alpb( 9, 1) =     3.136740d0 !    Fluorine -     Hydrogen
      xfac( 9, 1) =     0.815802d0 !    Fluorine -     Hydrogen
      alpb( 9, 2) =     2.856543d0 !    Fluorine -       Helium
      xfac( 9, 2) =     0.745107d0 !    Fluorine -       Helium
      alpb( 9, 3) =     3.043901d0 !    Fluorine -      Lithium
      xfac( 9, 3) =     1.975985d0 !    Fluorine -      Lithium
      alpb( 9, 4) =     3.726923d0 !    Fluorine -    Beryllium
      xfac( 9, 4) =     3.882993d0 !    Fluorine -    Beryllium
      alpb( 9, 5) =     2.823837d0 !    Fluorine -        Boron
      xfac( 9, 5) =     0.862761d0 !    Fluorine -        Boron
      alpb( 9, 6) =     3.027600d0 !    Fluorine -       Carbon
      xfac( 9, 6) =     0.732968d0 !    Fluorine -       Carbon
      alpb( 9, 7) =     2.856646d0 !    Fluorine -     Nitrogen
      xfac( 9, 7) =     0.635854d0 !    Fluorine -     Nitrogen
      alpb( 9, 8) =     3.015444d0 !    Fluorine -       Oxygen
      xfac( 9, 8) =     0.674251d0 !    Fluorine -       Oxygen
      alpb( 9, 9) =     3.175759d0 !    Fluorine -     Fluorine
      xfac( 9, 9) =     0.681343d0 !    Fluorine -     Fluorine
 !
      alpb(10, 1) =     5.999680d0 !        Neon -     Hydrogen
      xfac(10, 1) =     5.535021d0 !        Neon -     Hydrogen
      alpb(10, 2) =     3.677758d0 !        Neon -       Helium
      xfac(10, 2) =     1.960924d0 !        Neon -       Helium
      alpb(10, 3) =     2.193666d0 !        Neon -      Lithium
      xfac(10, 3) =     0.704958d0 !        Neon -      Lithium
      alpb(10, 4) =     1.316588d0 !        Neon -    Beryllium
      xfac(10, 4) =     0.392628d0 !        Neon -    Beryllium
      alpb(10, 5) =     2.756190d0 !        Neon -        Boron
      xfac(10, 5) =     2.764140d0 !        Neon -        Boron
      alpb(10, 6) =     3.441188d0 !        Neon -       Carbon
      xfac(10, 6) =     5.468780d0 !        Neon -       Carbon
      alpb(10, 7) =     4.426370d0 !        Neon -     Nitrogen
      xfac(10, 7) =    29.999609d0 !        Neon -     Nitrogen
      alpb(10, 8) =     2.889587d0 !        Neon -       Oxygen
      xfac(10, 8) =     0.763899d0 !        Neon -       Oxygen
      alpb(10, 9) =     3.675611d0 !        Neon -     Fluorine
      xfac(10, 9) =     2.706754d0 !        Neon -     Fluorine
      alpb(10,10) =     3.974567d0 !        Neon -         Neon
      xfac(10,10) =     2.794830d0 !        Neon -         Neon
 !
      alpb(11, 1) =     0.500326d0 !      Sodium -     Hydrogen
      xfac(11, 1) =     0.207831d0 !      Sodium -     Hydrogen
      alpb(11, 2) =     1.703029d0 !      Sodium -       Helium
      xfac(11, 2) =     4.282517d0 !      Sodium -       Helium
      alpb(11, 3) =     1.267299d0 !      Sodium -      Lithium
      xfac(11, 3) =     0.881482d0 !      Sodium -      Lithium
      alpb(11, 4) =     1.255480d0 !      Sodium -    Beryllium
      xfac(11, 4) =     3.121620d0 !      Sodium -    Beryllium
      alpb(11, 5) =     1.569961d0 !      Sodium -        Boron
      xfac(11, 5) =     3.188608d0 !      Sodium -        Boron
      alpb(11, 6) =     2.196050d0 !      Sodium -       Carbon
      xfac(11, 6) =     4.520429d0 !      Sodium -       Carbon
      alpb(11, 7) =     2.494384d0 !      Sodium -     Nitrogen
      xfac(11, 7) =     8.586387d0 !      Sodium -     Nitrogen
      alpb(11, 8) =     1.981449d0 !      Sodium -       Oxygen
      xfac(11, 8) =     3.270079d0 !      Sodium -       Oxygen
      alpb(11, 9) =     2.619551d0 !      Sodium -     Fluorine
      xfac(11, 9) =     7.047351d0 !      Sodium -     Fluorine
      alpb(11,10) =     1.774236d0 !      Sodium -         Neon
      xfac(11,10) =     1.343037d0 !      Sodium -         Neon
      alpb(11,11) =     0.446435d0 !      Sodium -       Sodium
      xfac(11,11) =     0.287137d0 !      Sodium -       Sodium
 !
      alpb(12, 1) =     2.651594d0 !   Magnesium -     Hydrogen
      xfac(12, 1) =     7.758237d0 !   Magnesium -     Hydrogen
      alpb(12, 2) =     2.210603d0 !   Magnesium -       Helium
      xfac(12, 2) =     3.725850d0 !   Magnesium -       Helium
      alpb(12, 3) =     1.184380d0 !   Magnesium -      Lithium
      xfac(12, 3) =     2.490250d0 !   Magnesium -      Lithium
      alpb(12, 4) =     1.557591d0 !   Magnesium -    Beryllium
      xfac(12, 4) =     2.066392d0 !   Magnesium -    Beryllium
      alpb(12, 5) =     2.527441d0 !   Magnesium -        Boron
      xfac(12, 5) =     6.146701d0 !   Magnesium -        Boron
      alpb(12, 6) =     3.040946d0 !   Magnesium -       Carbon
      xfac(12, 6) =    10.517690d0 !   Magnesium -       Carbon
      alpb(12, 7) =     2.079125d0 !   Magnesium -     Nitrogen
      xfac(12, 7) =     1.208075d0 !   Magnesium -     Nitrogen
      alpb(12, 8) =     2.251520d0 !   Magnesium -       Oxygen
      xfac(12, 8) =     1.535734d0 !   Magnesium -       Oxygen
      alpb(12, 9) =     3.362208d0 !   Magnesium -     Fluorine
      xfac(12, 9) =     5.859023d0 !   Magnesium -     Fluorine
      alpb(12,10) =     2.031676d0 !   Magnesium -         Neon
      xfac(12,10) =     1.214859d0 !   Magnesium -         Neon
      alpb(12,11) =     1.506773d0 !   Magnesium -       Sodium
      xfac(12,11) =     8.675619d0 !   Magnesium -       Sodium
      alpb(12,12) =     1.093573d0 !   Magnesium -    Magnesium
      xfac(12,12) =     0.465645d0 !   Magnesium -    Magnesium
 !
      alpb(13, 1) =     2.025996d0 !    Aluminum -     Hydrogen
      xfac(13, 1) =     2.958379d0 !    Aluminum -     Hydrogen
      alpb(13, 2) =     2.255830d0 !    Aluminum -       Helium
      xfac(13, 2) =     2.701400d0 !    Aluminum -       Helium
      alpb(13, 3) =     1.581593d0 !    Aluminum -      Lithium
      xfac(13, 3) =     1.106819d0 !    Aluminum -      Lithium
      alpb(13, 4) =     1.938237d0 !    Aluminum -    Beryllium
      xfac(13, 4) =     5.037214d0 !    Aluminum -    Beryllium
      alpb(13, 5) =     2.059569d0 !    Aluminum -        Boron
      xfac(13, 5) =     2.741479d0 !    Aluminum -        Boron
      alpb(13, 6) =     2.267440d0 !    Aluminum -       Carbon
      xfac(13, 6) =     2.928056d0 !    Aluminum -       Carbon
      alpb(13, 7) =     2.009754d0 !    Aluminum -     Nitrogen
      xfac(13, 7) =     1.345202d0 !    Aluminum -     Nitrogen
      alpb(13, 8) =     2.498660d0 !    Aluminum -       Oxygen
      xfac(13, 8) =     2.131396d0 !    Aluminum -       Oxygen
      alpb(13, 9) =     3.084258d0 !    Aluminum -     Fluorine
      xfac(13, 9) =     1.975635d0 !    Aluminum -     Fluorine
      alpb(13,10) =     2.447869d0 !    Aluminum -         Neon
      xfac(13,10) =     1.709200d0 !    Aluminum -         Neon
      alpb(13,11) =     1.202871d0 !    Aluminum -       Sodium
      xfac(13,11) =     2.071847d0 !    Aluminum -       Sodium
      alpb(13,12) =     1.972530d0 !    Aluminum -    Magnesium
      xfac(13,12) =    13.472443d0 !    Aluminum -    Magnesium
      alpb(13,13) =     1.387714d0 !    Aluminum -     Aluminum
      xfac(13,13) =     2.139200d0 !    Aluminum -     Aluminum
 !
      alpb(14, 1) =     1.896950d0 !     Silicon -     Hydrogen
      xfac(14, 1) =     0.924196d0 !     Silicon -     Hydrogen
      alpb(14, 2) =     2.040498d0 !     Silicon -       Helium
      xfac(14, 2) =     1.853583d0 !     Silicon -       Helium
      alpb(14, 3) =     1.789609d0 !     Silicon -      Lithium
      xfac(14, 3) =     3.090791d0 !     Silicon -      Lithium
      alpb(14, 4) =     1.263132d0 !     Silicon -    Beryllium
      xfac(14, 4) =     0.623433d0 !     Silicon -    Beryllium
      alpb(14, 5) =     1.982653d0 !     Silicon -        Boron
      xfac(14, 5) =     1.028287d0 !     Silicon -        Boron
      alpb(14, 6) =     1.984498d0 !     Silicon -       Carbon
      xfac(14, 6) =     0.785745d0 !     Silicon -       Carbon
      alpb(14, 7) =     1.818988d0 !     Silicon -     Nitrogen
      xfac(14, 7) =     0.592972d0 !     Silicon -     Nitrogen
      alpb(14, 8) =     1.923600d0 !     Silicon -       Oxygen
      xfac(14, 8) =     0.751095d0 !     Silicon -       Oxygen
      alpb(14, 9) =     2.131028d0 !     Silicon -     Fluorine
      xfac(14, 9) =     0.543516d0 !     Silicon -     Fluorine
      alpb(14,10) =     2.867784d0 !     Silicon -         Neon
      xfac(14,10) =    14.378676d0 !     Silicon -         Neon
      alpb(14,11) =     2.007615d0 !     Silicon -       Sodium
      xfac(14,11) =     9.237644d0 !     Silicon -       Sodium
      alpb(14,12) =     3.139749d0 !     Silicon -    Magnesium
      xfac(14,12) =    29.994520d0 !     Silicon -    Magnesium
      alpb(14,13) =     1.900000d0 !     Silicon -     Aluminum
      xfac(14,13) =     2.000000d0 !     Silicon -     Aluminum
      alpb(14,14) =     1.329000d0 !     Silicon -      Silicon
      xfac(14,14) =     0.273477d0 !     Silicon -      Silicon
 !
      alpb(15, 1) =     1.926537d0 !  Phosphorus -     Hydrogen
      xfac(15, 1) =     1.234986d0 !  Phosphorus -     Hydrogen
      alpb(15, 2) =     2.093158d0 !  Phosphorus -       Helium
      xfac(15, 2) =     1.490218d0 !  Phosphorus -       Helium
      alpb(15, 3) =     1.394544d0 !  Phosphorus -      Lithium
      xfac(15, 3) =     1.122950d0 !  Phosphorus -      Lithium
      alpb(15, 4) =     1.800070d0 !  Phosphorus -    Beryllium
      xfac(15, 4) =     1.684831d0 !  Phosphorus -    Beryllium
      alpb(15, 5) =     1.923168d0 !  Phosphorus -        Boron
      xfac(15, 5) =     1.450886d0 !  Phosphorus -        Boron
      alpb(15, 6) =     1.994653d0 !  Phosphorus -       Carbon
      xfac(15, 6) =     0.979512d0 !  Phosphorus -       Carbon
      alpb(15, 7) =     2.147042d0 !  Phosphorus -     Nitrogen
      xfac(15, 7) =     0.972154d0 !  Phosphorus -     Nitrogen
      alpb(15, 8) =     2.220768d0 !  Phosphorus -       Oxygen
      xfac(15, 8) =     0.878705d0 !  Phosphorus -       Oxygen
      alpb(15, 9) =     2.234356d0 !  Phosphorus -     Fluorine
      xfac(15, 9) =     0.514575d0 !  Phosphorus -     Fluorine
      alpb(15,10) =     2.219036d0 !  Phosphorus -         Neon
      xfac(15,10) =     0.774954d0 !  Phosphorus -         Neon
      alpb(15,11) =     1.500320d0 !  Phosphorus -       Sodium
      xfac(15,11) =     2.837095d0 !  Phosphorus -       Sodium
      alpb(15,12) =     1.383773d0 !  Phosphorus -    Magnesium
      xfac(15,12) =     1.177881d0 !  Phosphorus -    Magnesium
      alpb(15,13) =     1.980727d0 !  Phosphorus -     Aluminum
      xfac(15,13) =     5.050816d0 !  Phosphorus -     Aluminum
      alpb(15,14) =     3.313466d0 !  Phosphorus -      Silicon
      xfac(15,14) =    13.239121d0 !  Phosphorus -      Silicon
      alpb(15,15) =     1.505792d0 !  Phosphorus -   Phosphorus
      xfac(15,15) =     0.902501d0 !  Phosphorus -   Phosphorus
 !
      alpb(16, 1) =     2.215975d0 !      Sulfur -     Hydrogen
      xfac(16, 1) =     0.849712d0 !      Sulfur -     Hydrogen
      alpb(16, 2) =     1.959149d0 !      Sulfur -       Helium
      xfac(16, 2) =     0.437618d0 !      Sulfur -       Helium
      alpb(16, 3) =     2.294275d0 !      Sulfur -      Lithium
      xfac(16, 3) =     2.642502d0 !      Sulfur -      Lithium
      alpb(16, 4) =     2.781736d0 !      Sulfur -    Beryllium
      xfac(16, 4) =     3.791565d0 !      Sulfur -    Beryllium
      alpb(16, 5) =     2.403696d0 !      Sulfur -        Boron
      xfac(16, 5) =     1.125394d0 !      Sulfur -        Boron
      alpb(16, 6) =     2.210305d0 !      Sulfur -       Carbon
      xfac(16, 6) =     0.666849d0 !      Sulfur -       Carbon
      alpb(16, 7) =     2.289990d0 !      Sulfur -     Nitrogen
      xfac(16, 7) =     0.738710d0 !      Sulfur -     Nitrogen
      alpb(16, 8) =     2.383289d0 !      Sulfur -       Oxygen
      xfac(16, 8) =     0.747215d0 !      Sulfur -       Oxygen
      alpb(16, 9) =     2.187186d0 !      Sulfur -     Fluorine
      xfac(16, 9) =     0.375251d0 !      Sulfur -     Fluorine
      alpb(16,10) =     2.787058d0 !      Sulfur -         Neon
      xfac(16,10) =     3.296160d0 !      Sulfur -         Neon
      alpb(16,11) =     1.400850d0 !      Sulfur -       Sodium
      xfac(16,11) =     0.852434d0 !      Sulfur -       Sodium
      alpb(16,12) =     1.500163d0 !      Sulfur -    Magnesium
      xfac(16,12) =     0.500748d0 !      Sulfur -    Magnesium
      alpb(16,13) =     1.976705d0 !      Sulfur -     Aluminum
      xfac(16,13) =     2.347384d0 !      Sulfur -     Aluminum
      alpb(16,14) =     1.885916d0 !      Sulfur -      Silicon
      xfac(16,14) =     0.876658d0 !      Sulfur -      Silicon
      alpb(16,15) =     1.595325d0 !      Sulfur -   Phosphorus
      xfac(16,15) =     0.562266d0 !      Sulfur -   Phosphorus
      alpb(16,16) =     1.794556d0 !      Sulfur -       Sulfur
      xfac(16,16) =     0.473856d0 !      Sulfur -       Sulfur
 !
      alpb(17, 1) =     2.402886d0 !    Chlorine -     Hydrogen
      xfac(17, 1) =     0.754831d0 !    Chlorine -     Hydrogen
      alpb(17, 2) =     1.671677d0 !    Chlorine -       Helium
      xfac(17, 2) =     0.272964d0 !    Chlorine -       Helium
      alpb(17, 3) =     2.783001d0 !    Chlorine -      Lithium
      xfac(17, 3) =     4.227794d0 !    Chlorine -      Lithium
      alpb(17, 4) =     2.822676d0 !    Chlorine -    Beryllium
      xfac(17, 4) =     2.507275d0 !    Chlorine -    Beryllium
      alpb(17, 5) =     2.259323d0 !    Chlorine -        Boron
      xfac(17, 5) =     0.822129d0 !    Chlorine -        Boron
      alpb(17, 6) =     2.162197d0 !    Chlorine -       Carbon
      xfac(17, 6) =     0.515787d0 !    Chlorine -       Carbon
      alpb(17, 7) =     2.172134d0 !    Chlorine -     Nitrogen
      xfac(17, 7) =     0.520745d0 !    Chlorine -     Nitrogen
      alpb(17, 8) =     2.323236d0 !    Chlorine -       Oxygen
      xfac(17, 8) =     0.585510d0 !    Chlorine -       Oxygen
      alpb(17, 9) =     2.313270d0 !    Chlorine -     Fluorine
      xfac(17, 9) =     0.411124d0 !    Chlorine -     Fluorine
      alpb(17,10) =     1.703151d0 !    Chlorine -         Neon
      xfac(17,10) =     0.125133d0 !    Chlorine -         Neon
      alpb(17,11) =     1.816429d0 !    Chlorine -       Sodium
      xfac(17,11) =     1.357894d0 !    Chlorine -       Sodium
      alpb(17,12) =     2.391806d0 !    Chlorine -    Magnesium
      xfac(17,12) =     2.430856d0 !    Chlorine -    Magnesium
      alpb(17,13) =     2.125939d0 !    Chlorine -     Aluminum
      xfac(17,13) =     2.153451d0 !    Chlorine -     Aluminum
      alpb(17,14) =     1.684978d0 !    Chlorine -      Silicon
      xfac(17,14) =     0.513000d0 !    Chlorine -      Silicon
      alpb(17,15) =     1.468306d0 !    Chlorine -   Phosphorus
      xfac(17,15) =     0.352361d0 !    Chlorine -   Phosphorus
      alpb(17,16) =     1.715435d0 !    Chlorine -       Sulfur
      xfac(17,16) =     0.356971d0 !    Chlorine -       Sulfur
      alpb(17,17) =     1.823239d0 !    Chlorine -     Chlorine
      xfac(17,17) =     0.332919d0 !    Chlorine -     Chlorine
 !
      alpb(18, 1) =     4.056167d0 !       Argon -     Hydrogen
      xfac(18, 1) =     3.933445d0 !       Argon -     Hydrogen
      alpb(18, 2) =     2.716562d0 !       Argon -       Helium
      xfac(18, 2) =     1.177211d0 !       Argon -       Helium
      alpb(18, 3) =     3.122895d0 !       Argon -      Lithium
      xfac(18, 3) =     3.362910d0 !       Argon -      Lithium
      alpb(18, 4) =     3.044007d0 !       Argon -    Beryllium
      xfac(18, 4) =     2.755492d0 !       Argon -    Beryllium
      alpb(18, 5) =     2.415471d0 !       Argon -        Boron
      xfac(18, 5) =     1.931586d0 !       Argon -        Boron
      alpb(18, 6) =     1.471309d0 !       Argon -       Carbon
      xfac(18, 6) =     0.122309d0 !       Argon -       Carbon
      alpb(18, 7) =     2.326805d0 !       Argon -     Nitrogen
      xfac(18, 7) =     0.562581d0 !       Argon -     Nitrogen
      alpb(18, 8) =     2.240673d0 !       Argon -       Oxygen
      xfac(18, 8) =     0.355795d0 !       Argon -       Oxygen
      alpb(18, 9) =     3.920658d0 !       Argon -     Fluorine
      xfac(18, 9) =     9.269715d0 !       Argon -     Fluorine
      alpb(18,10) =     2.963747d0 !       Argon -         Neon
      xfac(18,10) =     1.304697d0 !       Argon -         Neon
      alpb(18,11) =     2.167677d0 !       Argon -       Sodium
      xfac(18,11) =     3.398138d0 !       Argon -       Sodium
      alpb(18,12) =     2.092664d0 !       Argon -    Magnesium
      xfac(18,12) =     1.970638d0 !       Argon -    Magnesium
      alpb(18,13) =     2.645165d0 !       Argon -     Aluminum
      xfac(18,13) =     1.852009d0 !       Argon -     Aluminum
      alpb(18,14) =     1.780350d0 !       Argon -      Silicon
      xfac(18,14) =     1.067890d0 !       Argon -      Silicon
      alpb(18,15) =     4.372516d0 !       Argon -   Phosphorus
      xfac(18,15) =     0.171014d0 !       Argon -   Phosphorus
      alpb(18,16) =     2.049398d0 !       Argon -       Sulfur
      xfac(18,16) =     0.653769d0 !       Argon -       Sulfur
      alpb(18,17) =     2.554449d0 !       Argon -     Chlorine
      xfac(18,17) =     2.256094d0 !       Argon -     Chlorine
      alpb(18,18) =     2.306432d0 !       Argon -        Argon
      xfac(18,18) =     0.972699d0 !       Argon -        Argon
 !
      alpb(19, 1) =     0.648173d0 !   Potassium -     Hydrogen
      xfac(19, 1) =     0.369340d0 !   Potassium -     Hydrogen
      alpb(19, 2) =     1.418501d0 !   Potassium -       Helium
      xfac(19, 2) =     2.895045d0 !   Potassium -       Helium
      alpb(19, 3) =     1.036487d0 !   Potassium -      Lithium
      xfac(19, 3) =     4.374567d0 !   Potassium -      Lithium
      alpb(19, 4) =     1.931888d0 !   Potassium -    Beryllium
      xfac(19, 4) =     6.732221d0 !   Potassium -    Beryllium
      alpb(19, 5) =     2.031768d0 !   Potassium -        Boron
      xfac(19, 5) =     8.900541d0 !   Potassium -        Boron
      alpb(19, 6) =     2.241757d0 !   Potassium -       Carbon
      xfac(19, 6) =    10.317987d0 !   Potassium -       Carbon
      alpb(19, 7) =     2.325859d0 !   Potassium -     Nitrogen
      xfac(19, 7) =     7.977707d0 !   Potassium -     Nitrogen
      alpb(19, 8) =     1.508571d0 !   Potassium -       Oxygen
      xfac(19, 8) =     1.012275d0 !   Potassium -       Oxygen
      alpb(19, 9) =     3.182817d0 !   Potassium -     Fluorine
      xfac(19, 9) =     6.592971d0 !   Potassium -     Fluorine
      alpb(19,10) =     1.138021d0 !   Potassium -         Neon
      xfac(19,10) =     0.233995d0 !   Potassium -         Neon
      alpb(19,11) =     0.884307d0 !   Potassium -       Sodium
      xfac(19,11) =     5.563027d0 !   Potassium -       Sodium
      alpb(19,12) =     0.884810d0 !   Potassium -    Magnesium
      xfac(19,12) =     3.290502d0 !   Potassium -    Magnesium
      alpb(19,13) =     1.976076d0 !   Potassium -     Aluminum
      xfac(19,13) =    29.944708d0 !   Potassium -     Aluminum
      alpb(19,14) =     1.675930d0 !   Potassium -      Silicon
      xfac(19,14) =     8.279200d0 !   Potassium -      Silicon
      alpb(19,15) =     1.443738d0 !   Potassium -   Phosphorus
      xfac(19,15) =     4.475384d0 !   Potassium -   Phosphorus
      alpb(19,16) =     2.512156d0 !   Potassium -       Sulfur
      xfac(19,16) =    29.528951d0 !   Potassium -       Sulfur
      alpb(19,17) =     1.622163d0 !   Potassium -     Chlorine
      xfac(19,17) =     1.231481d0 !   Potassium -     Chlorine
      alpb(19,18) =     2.302803d0 !   Potassium -        Argon
      xfac(19,18) =     9.710508d0 !   Potassium -        Argon
      alpb(19,19) =     1.435514d0 !   Potassium -    Potassium
      xfac(19,19) =     5.934329d0 !   Potassium -    Potassium
 !
      alpb(20, 1) =     2.141859d0 !     Calcium -     Hydrogen
      xfac(20, 1) =     7.728606d0 !     Calcium -     Hydrogen
      alpb(20, 2) =     1.719847d0 !     Calcium -       Helium
      xfac(20, 2) =     2.913852d0 !     Calcium -       Helium
      alpb(20, 5) =     1.700010d0 !     Calcium -        Boron
      xfac(20, 5) =     1.700010d0 !     Calcium -        Boron
      alpb(20, 6) =     1.035305d0 !     Calcium -       Carbon
      xfac(20, 6) =     0.148450d0 !     Calcium -       Carbon
      alpb(20, 7) =     2.386600d0 !     Calcium -     Nitrogen
      xfac(20, 7) =     2.988074d0 !     Calcium -     Nitrogen
      alpb(20, 8) =     3.263897d0 !     Calcium -       Oxygen
      xfac(20, 8) =    17.028946d0 !     Calcium -       Oxygen
      alpb(20, 9) =     2.645053d0 !     Calcium -     Fluorine
      xfac(20, 9) =     3.482821d0 !     Calcium -     Fluorine
      alpb(20,10) =     0.954530d0 !     Calcium -         Neon
      xfac(20,10) =     0.332586d0 !     Calcium -         Neon
      alpb(20,11) =     3.107104d0 !     Calcium -       Sodium
      xfac(20,11) =     9.657509d0 !     Calcium -       Sodium
      alpb(20,12) =     2.299800d0 !     Calcium -    Magnesium
      xfac(20,12) =     8.599800d0 !     Calcium -    Magnesium
      alpb(20,13) =     1.612565d0 !     Calcium -     Aluminum
      xfac(20,13) =     4.188555d0 !     Calcium -     Aluminum
      alpb(20,14) =     1.218788d0 !     Calcium -      Silicon
      xfac(20,14) =     0.336233d0 !     Calcium -      Silicon
      alpb(20,15) =     1.024142d0 !     Calcium -   Phosphorus
      xfac(20,15) =     0.410840d0 !     Calcium -   Phosphorus
      alpb(20,16) =     0.958171d0 !     Calcium -       Sulfur
      xfac(20,16) =     0.325739d0 !     Calcium -       Sulfur
      alpb(20,17) =     2.383391d0 !     Calcium -     Chlorine
      xfac(20,17) =     5.956144d0 !     Calcium -     Chlorine
      alpb(20,18) =     1.034881d0 !     Calcium -        Argon
      xfac(20,18) =     0.291072d0 !     Calcium -        Argon
      alpb(20,19) =     1.119200d0 !     Calcium -    Potassium
      xfac(20,19) =     1.240320d0 !     Calcium -    Potassium
      alpb(20,20) =     1.889674d0 !     Calcium -      Calcium
      xfac(20,20) =    30.003591d0 !     Calcium -      Calcium
 !
      alpb(21, 1) =     1.179485d0 !    Scandium -     Hydrogen
      xfac(21, 1) =     0.351199d0 !    Scandium -     Hydrogen
      alpb(21, 6) =     2.630490d0 !    Scandium -       Carbon
      xfac(21, 6) =     8.608052d0 !    Scandium -       Carbon
      alpb(21, 7) =     2.270004d0 !    Scandium -     Nitrogen
      xfac(21, 7) =     3.231881d0 !    Scandium -     Nitrogen
      alpb(21, 8) =     2.256516d0 !    Scandium -       Oxygen
      xfac(21, 8) =     3.058672d0 !    Scandium -       Oxygen
      alpb(21, 9) =     3.107985d0 !    Scandium -     Fluorine
      xfac(21, 9) =     7.252347d0 !    Scandium -     Fluorine
      alpb(21,13) =     1.003550d0 !    Scandium -     Aluminum
      xfac(21,13) =     0.500620d0 !    Scandium -     Aluminum
      alpb(21,14) =     2.016870d0 !    Scandium -      Silicon
      xfac(21,14) =     3.219070d0 !    Scandium -      Silicon
      alpb(21,15) =     0.868165d0 !    Scandium -   Phosphorus
      xfac(21,15) =     0.626749d0 !    Scandium -   Phosphorus
      alpb(21,16) =     0.422939d0 !    Scandium -       Sulfur
      xfac(21,16) =     0.211850d0 !    Scandium -       Sulfur
      alpb(21,17) =     2.141474d0 !    Scandium -     Chlorine
      xfac(21,17) =     2.996129d0 !    Scandium -     Chlorine
      alpb(21,21) =     1.132838d0 !    Scandium -     Scandium
      xfac(21,21) =     2.598166d0 !    Scandium -     Scandium
 !
      alpb(22, 1) =     0.832669d0 !    Titanium -     Hydrogen
      xfac(22, 1) =     0.143722d0 !    Titanium -     Hydrogen
      alpb(22, 5) =     1.628710d0 !    Titanium -        Boron
      xfac(22, 5) =     0.649360d0 !    Titanium -        Boron
      alpb(22, 6) =     1.597973d0 !    Titanium -       Carbon
      xfac(22, 6) =     0.416706d0 !    Titanium -       Carbon
      alpb(22, 7) =     1.678686d0 !    Titanium -     Nitrogen
      xfac(22, 7) =     0.545461d0 !    Titanium -     Nitrogen
      alpb(22, 8) =     1.789118d0 !    Titanium -       Oxygen
      xfac(22, 8) =     0.799486d0 !    Titanium -       Oxygen
      alpb(22, 9) =     2.307087d0 !    Titanium -     Fluorine
      xfac(22, 9) =     1.085742d0 !    Titanium -     Fluorine
      alpb(22,12) =     1.911340d0 !    Titanium -    Magnesium
      xfac(22,12) =     4.330240d0 !    Titanium -    Magnesium
      alpb(22,13) =     1.369486d0 !    Titanium -     Aluminum
      xfac(22,13) =     2.091841d0 !    Titanium -     Aluminum
      alpb(22,14) =     2.856038d0 !    Titanium -      Silicon
      xfac(22,14) =     6.773815d0 !    Titanium -      Silicon
      alpb(22,15) =     2.151929d0 !    Titanium -   Phosphorus
      xfac(22,15) =     4.150500d0 !    Titanium -   Phosphorus
      alpb(22,16) =     1.846439d0 !    Titanium -       Sulfur
      xfac(22,16) =     0.943784d0 !    Titanium -       Sulfur
      alpb(22,17) =     1.461034d0 !    Titanium -     Chlorine
      xfac(22,17) =     0.333297d0 !    Titanium -     Chlorine
      alpb(22,20) =     2.000000d0 !    Titanium -      Calcium
      xfac(22,20) =     4.109141d0 !    Titanium -      Calcium
      alpb(22,22) =     2.648597d0 !    Titanium -     Titanium
      xfac(22,22) =     2.000000d0 !    Titanium -     Titanium
 !
      alpb(23, 1) =     1.280133d0 !    Vanadium -     Hydrogen
      xfac(23, 1) =     0.105204d0 !    Vanadium -     Hydrogen
      alpb(23, 6) =     2.789855d0 !    Vanadium -       Carbon
      xfac(23, 6) =     1.938760d0 !    Vanadium -       Carbon
      alpb(23, 7) =     1.607540d0 !    Vanadium -     Nitrogen
      xfac(23, 7) =     0.276725d0 !    Vanadium -     Nitrogen
      alpb(23, 8) =     1.623973d0 !    Vanadium -       Oxygen
      xfac(23, 8) =     0.415312d0 !    Vanadium -       Oxygen
      alpb(23, 9) =     1.825160d0 !    Vanadium -     Fluorine
      xfac(23, 9) =     0.342815d0 !    Vanadium -     Fluorine
      alpb(23,11) =     2.551010d0 !    Vanadium -       Sodium
      xfac(23,11) =     8.276020d0 !    Vanadium -       Sodium
      alpb(23,15) =     2.549154d0 !    Vanadium -   Phosphorus
      xfac(23,15) =     6.250624d0 !    Vanadium -   Phosphorus
      alpb(23,16) =     2.704124d0 !    Vanadium -       Sulfur
      xfac(23,16) =     2.035039d0 !    Vanadium -       Sulfur
      alpb(23,17) =     1.688529d0 !    Vanadium -     Chlorine
      xfac(23,17) =     0.243657d0 !    Vanadium -     Chlorine
      alpb(23,19) =     4.521360d0 !    Vanadium -    Potassium
      xfac(23,19) =     2.026590d0 !    Vanadium -    Potassium
      alpb(23,23) =     4.832391d0 !    Vanadium -     Vanadium
      xfac(23,23) =    10.779892d0 !    Vanadium -     Vanadium
 !
      alpb(24, 1) =     0.882661d0 !    Chromium -     Hydrogen
      xfac(24, 1) =     0.044469d0 !    Chromium -     Hydrogen
      alpb(24, 6) =     3.656754d0 !    Chromium -       Carbon
      xfac(24, 6) =     6.110187d0 !    Chromium -       Carbon
      alpb(24, 7) =     3.029186d0 !    Chromium -     Nitrogen
      xfac(24, 7) =     1.920324d0 !    Chromium -     Nitrogen
      alpb(24, 8) =     2.500000d0 !    Chromium -       Oxygen
      xfac(24, 8) =     1.055511d0 !    Chromium -       Oxygen
      alpb(24, 9) =     2.716521d0 !    Chromium -     Fluorine
      xfac(24, 9) =     0.737607d0 !    Chromium -     Fluorine
      alpb(24,11) =     2.295056d0 !    Chromium -       Sodium
      xfac(24,11) =     8.364274d0 !    Chromium -       Sodium
      alpb(24,14) =     1.860760d0 !    Chromium -      Silicon
      xfac(24,14) =     1.029110d0 !    Chromium -      Silicon
      alpb(24,15) =     1.695383d0 !    Chromium -   Phosphorus
      xfac(24,15) =     0.600177d0 !    Chromium -   Phosphorus
      alpb(24,16) =     2.260978d0 !    Chromium -       Sulfur
      xfac(24,16) =     0.550334d0 !    Chromium -       Sulfur
      alpb(24,17) =     2.152618d0 !    Chromium -     Chlorine
      xfac(24,17) =     0.369073d0 !    Chromium -     Chlorine
      alpb(24,19) =     2.000000d0 !    Chromium -    Potassium
      xfac(24,19) =     2.000000d0 !    Chromium -    Potassium
      alpb(24,24) =     4.655419d0 !    Chromium -     Chromium
      xfac(24,24) =    10.318607d0 !    Chromium -     Chromium
 !
      alpb(25, 1) =     2.309940d0 !   Manganese -     Hydrogen
      xfac(25, 1) =     1.269210d0 !   Manganese -     Hydrogen
      alpb(25, 6) =     3.000750d0 !   Manganese -       Carbon
      xfac(25, 6) =     2.583110d0 !   Manganese -       Carbon
      alpb(25, 7) =     2.921470d0 !   Manganese -     Nitrogen
      xfac(25, 7) =     1.956750d0 !   Manganese -     Nitrogen
      alpb(25, 8) =     2.577540d0 !   Manganese -       Oxygen
      xfac(25, 8) =     1.285620d0 !   Manganese -       Oxygen
      alpb(25, 9) =     2.791950d0 !   Manganese -     Fluorine
      xfac(25, 9) =     1.113070d0 !   Manganese -     Fluorine
      alpb(25,13) =     1.768360d0 !   Manganese -     Aluminum
      xfac(25,13) =     1.040790d0 !   Manganese -     Aluminum
      alpb(25,14) =     1.937959d0 !   Manganese -      Silicon
      xfac(25,14) =     0.950580d0 !   Manganese -      Silicon
      alpb(25,15) =     1.947020d0 !   Manganese -   Phosphorus
      xfac(25,15) =     1.130320d0 !   Manganese -   Phosphorus
      alpb(25,16) =     2.482510d0 !   Manganese -       Sulfur
      xfac(25,16) =     1.612650d0 !   Manganese -       Sulfur
      alpb(25,17) =     1.657010d0 !   Manganese -     Chlorine
      xfac(25,17) =     0.201850d0 !   Manganese -     Chlorine
      alpb(25,20) =     1.491440d0 !   Manganese -      Calcium
      xfac(25,20) =     0.620180d0 !   Manganese -      Calcium
      alpb(25,25) =     2.665420d0 !   Manganese -    Manganese
      xfac(25,25) =     2.460040d0 !   Manganese -    Manganese
 !
      alpb(26, 1) =     0.854488d0 !        Iron -     Hydrogen
      xfac(26, 1) =     0.025195d0 !        Iron -     Hydrogen
      alpb(26, 6) =     3.991343d0 !        Iron -       Carbon
      xfac(26, 6) =     0.366835d0 !        Iron -       Carbon
      alpb(26, 7) =     2.500486d0 !        Iron -     Nitrogen
      xfac(26, 7) =     0.155342d0 !        Iron -     Nitrogen
      alpb(26, 8) =     1.726313d0 !        Iron -       Oxygen
      xfac(26, 8) =     0.136422d0 !        Iron -       Oxygen
      alpb(26, 9) =     4.294707d0 !        Iron -     Fluorine
      xfac(26, 9) =     3.657350d0 !        Iron -     Fluorine
      alpb(26,15) =     2.567534d0 !        Iron -   Phosphorus
      xfac(26,15) =     0.431291d0 !        Iron -   Phosphorus
      alpb(26,16) =     0.988991d0 !        Iron -       Sulfur
      xfac(26,16) =     0.033478d0 !        Iron -       Sulfur
      alpb(26,17) =     1.229793d0 !        Iron -     Chlorine
      xfac(26,17) =     0.019473d0 !        Iron -     Chlorine
      alpb(26,19) =     2.000000d0 !        Iron -    Potassium
      xfac(26,19) =     6.000000d0 !        Iron -    Potassium
      alpb(26,26) =     2.720785d0 !        Iron -         Iron
      xfac(26,26) =     1.846890d0 !        Iron -         Iron
 !
      alpb(27, 1) =     2.966518d0 !      Cobalt -     Hydrogen
      xfac(27, 1) =     2.472465d0 !      Cobalt -     Hydrogen
      alpb(27, 5) =     3.200000d0 !      Cobalt -        Boron
      xfac(27, 5) =     1.000000d0 !      Cobalt -        Boron
      alpb(27, 6) =     3.716233d0 !      Cobalt -       Carbon
      xfac(27, 6) =     2.123930d0 !      Cobalt -       Carbon
      alpb(27, 7) =     3.618638d0 !      Cobalt -     Nitrogen
      xfac(27, 7) =     2.653836d0 !      Cobalt -     Nitrogen
      alpb(27, 8) =     3.726911d0 !      Cobalt -       Oxygen
      xfac(27, 8) =     5.252022d0 !      Cobalt -       Oxygen
      alpb(27, 9) =     3.956347d0 !      Cobalt -     Fluorine
      xfac(27, 9) =     4.585030d0 !      Cobalt -     Fluorine
      alpb(27,14) =     2.469805d0 !      Cobalt -      Silicon
      xfac(27,14) =     1.090240d0 !      Cobalt -      Silicon
      alpb(27,15) =     1.152505d0 !      Cobalt -   Phosphorus
      xfac(27,15) =     0.105936d0 !      Cobalt -   Phosphorus
      alpb(27,16) =     2.429255d0 !      Cobalt -       Sulfur
      xfac(27,16) =     0.436707d0 !      Cobalt -       Sulfur
      alpb(27,17) =     3.217497d0 !      Cobalt -     Chlorine
      xfac(27,17) =     1.033414d0 !      Cobalt -     Chlorine
      alpb(27,27) =     3.288166d0 !      Cobalt -       Cobalt
      xfac(27,27) =     3.919618d0 !      Cobalt -       Cobalt
 !
      alpb(28, 1) =     2.635280d0 !      Nickel -     Hydrogen
      xfac(28, 1) =     1.763124d0 !      Nickel -     Hydrogen
      alpb(28, 6) =     4.285513d0 !      Nickel -       Carbon
      xfac(28, 6) =     7.133324d0 !      Nickel -       Carbon
      alpb(28, 7) =     3.845215d0 !      Nickel -     Nitrogen
      xfac(28, 7) =     4.286800d0 !      Nickel -     Nitrogen
      alpb(28, 8) =     2.937232d0 !      Nickel -       Oxygen
      xfac(28, 8) =     0.885942d0 !      Nickel -       Oxygen
      alpb(28, 9) =     3.440241d0 !      Nickel -     Fluorine
      xfac(28, 9) =     1.088208d0 !      Nickel -     Fluorine
      alpb(28,14) =     2.068881d0 !      Nickel -      Silicon
      xfac(28,14) =     0.938646d0 !      Nickel -      Silicon
      alpb(28,15) =     3.260283d0 !      Nickel -   Phosphorus
      xfac(28,15) =     5.059727d0 !      Nickel -   Phosphorus
      alpb(28,16) =     2.002752d0 !      Nickel -       Sulfur
      xfac(28,16) =     0.274852d0 !      Nickel -       Sulfur
      alpb(28,17) =     2.200512d0 !      Nickel -     Chlorine
      xfac(28,17) =     0.202313d0 !      Nickel -     Chlorine
      alpb(28,28) =     1.097960d0 !      Nickel -       Nickel
      xfac(28,28) =     0.035474d0 !      Nickel -       Nickel
 !
      alpb(29, 1) =     2.335359d0 !      Copper -     Hydrogen
      xfac(29, 1) =     0.603591d0 !      Copper -     Hydrogen
      alpb(29, 6) =     4.638773d0 !      Copper -       Carbon
      xfac(29, 6) =     7.067794d0 !      Copper -       Carbon
      alpb(29, 7) =     4.214337d0 !      Copper -     Nitrogen
      xfac(29, 7) =     3.228667d0 !      Copper -     Nitrogen
      alpb(29, 8) =     3.959951d0 !      Copper -       Oxygen
      xfac(29, 8) =     2.000000d0 !      Copper -       Oxygen
      alpb(29, 9) =     4.478832d0 !      Copper -     Fluorine
      xfac(29, 9) =     1.282108d0 !      Copper -     Fluorine
      alpb(29,15) =     0.210640d0 !      Copper -   Phosphorus
      xfac(29,15) =     0.020126d0 !      Copper -   Phosphorus
      alpb(29,16) =     0.273112d0 !      Copper -       Sulfur
      xfac(29,16) =     0.005248d0 !      Copper -       Sulfur
      alpb(29,17) =     2.776531d0 !      Copper -     Chlorine
      xfac(29,17) =     0.139065d0 !      Copper -     Chlorine
      alpb(29,29) =     3.616846d0 !      Copper -       Copper
      xfac(29,29) =     5.184376d0 !      Copper -       Copper
 !
      alpb(30, 1) =     1.987891d0 !        Zinc -     Hydrogen
      xfac(30, 1) =     3.109193d0 !        Zinc -     Hydrogen
      alpb(30, 6) =     1.802327d0 !        Zinc -       Carbon
      xfac(30, 6) =     0.991465d0 !        Zinc -       Carbon
      alpb(30, 7) =     1.844579d0 !        Zinc -     Nitrogen
      xfac(30, 7) =     0.952476d0 !        Zinc -     Nitrogen
      alpb(30, 8) =     2.335054d0 !        Zinc -       Oxygen
      xfac(30, 8) =     2.265313d0 !        Zinc -       Oxygen
      alpb(30, 9) =     2.410021d0 !        Zinc -     Fluorine
      xfac(30, 9) =     1.225545d0 !        Zinc -     Fluorine
      alpb(30,14) =     1.832058d0 !        Zinc -      Silicon
      xfac(30,14) =     3.783905d0 !        Zinc -      Silicon
      alpb(30,15) =     1.220480d0 !        Zinc -   Phosphorus
      xfac(30,15) =     0.581530d0 !        Zinc -   Phosphorus
      alpb(30,16) =     1.455000d0 !        Zinc -       Sulfur
      xfac(30,16) =     0.648000d0 !        Zinc -       Sulfur
      alpb(30,17) =     1.625176d0 !        Zinc -     Chlorine
      xfac(30,17) =     0.721351d0 !        Zinc -     Chlorine
      alpb(30,20) =     1.119180d0 !        Zinc -      Calcium
      xfac(30,20) =     1.240290d0 !        Zinc -      Calcium
      alpb(30,30) =     0.929000d0 !        Zinc -         Zinc
      xfac(30,30) =     0.465000d0 !        Zinc -         Zinc
 !
      alpb(31, 1) =     1.847350d0 !     Gallium -     Hydrogen
      xfac(31, 1) =     1.386652d0 !     Gallium -     Hydrogen
      alpb(31, 6) =     2.325410d0 !     Gallium -       Carbon
      xfac(31, 6) =     1.962990d0 !     Gallium -       Carbon
      alpb(31, 7) =     2.121820d0 !     Gallium -     Nitrogen
      xfac(31, 7) =     1.188338d0 !     Gallium -     Nitrogen
      alpb(31, 8) =     2.348347d0 !     Gallium -       Oxygen
      xfac(31, 8) =     1.523644d0 !     Gallium -       Oxygen
      alpb(31, 9) =     2.679869d0 !     Gallium -     Fluorine
      xfac(31, 9) =     1.416942d0 !     Gallium -     Fluorine
      alpb(31,14) =     1.913780d0 !     Gallium -      Silicon
      xfac(31,14) =     1.002290d0 !     Gallium -      Silicon
      alpb(31,15) =     2.979650d0 !     Gallium -   Phosphorus
      xfac(31,15) =     0.500000d0 !     Gallium -   Phosphorus
      alpb(31,16) =     2.232108d0 !     Gallium -       Sulfur
      xfac(31,16) =     2.456284d0 !     Gallium -       Sulfur
      alpb(31,17) =     2.024710d0 !     Gallium -     Chlorine
      xfac(31,17) =     1.186661d0 !     Gallium -     Chlorine
      alpb(31,31) =     1.334643d0 !     Gallium -      Gallium
      xfac(31,31) =     1.198394d0 !     Gallium -      Gallium
 !
      alpb(32, 1) =     2.206793d0 !   Germanium -     Hydrogen
      xfac(32, 1) =     1.733226d0 !   Germanium -     Hydrogen
      alpb(32, 6) =     2.257469d0 !   Germanium -       Carbon
      xfac(32, 6) =     1.297510d0 !   Germanium -       Carbon
      alpb(32, 7) =     1.988226d0 !   Germanium -     Nitrogen
      xfac(32, 7) =     0.637506d0 !   Germanium -     Nitrogen
      alpb(32, 8) =     2.139413d0 !   Germanium -       Oxygen
      xfac(32, 8) =     0.826964d0 !   Germanium -       Oxygen
      alpb(32, 9) =     2.384777d0 !   Germanium -     Fluorine
      xfac(32, 9) =     0.651977d0 !   Germanium -     Fluorine
      alpb(32,14) =     0.299721d0 !   Germanium -      Silicon
      xfac(32,14) =     0.178680d0 !   Germanium -      Silicon
      alpb(32,15) =     2.469291d0 !   Germanium -   Phosphorus
      xfac(32,15) =     5.616349d0 !   Germanium -   Phosphorus
      alpb(32,16) =     2.024588d0 !   Germanium -       Sulfur
      xfac(32,16) =     1.160957d0 !   Germanium -       Sulfur
      alpb(32,17) =     1.771228d0 !   Germanium -     Chlorine
      xfac(32,17) =     0.545239d0 !   Germanium -     Chlorine
      alpb(32,25) =     2.382834d0 !   Germanium -    Manganese
      xfac(32,25) =     2.255151d0 !   Germanium -    Manganese
      alpb(32,27) =     2.852610d0 !   Germanium -       Cobalt
      xfac(32,27) =     2.151850d0 !   Germanium -       Cobalt
      alpb(32,32) =     2.019000d0 !   Germanium -    Germanium
      xfac(32,32) =     3.023000d0 !   Germanium -    Germanium
 !
      alpb(33, 1) =     1.993527d0 !     Arsenic -     Hydrogen
      xfac(33, 1) =     1.090589d0 !     Arsenic -     Hydrogen
      alpb(33, 6) =     1.855069d0 !     Arsenic -       Carbon
      xfac(33, 6) =     0.579098d0 !     Arsenic -       Carbon
      alpb(33, 7) =     1.496543d0 !     Arsenic -     Nitrogen
      xfac(33, 7) =     0.273337d0 !     Arsenic -     Nitrogen
      alpb(33, 8) =     2.003950d0 !     Arsenic -       Oxygen
      xfac(33, 8) =     0.701614d0 !     Arsenic -       Oxygen
      alpb(33, 9) =     2.012583d0 !     Arsenic -     Fluorine
      xfac(33, 9) =     0.402628d0 !     Arsenic -     Fluorine
      alpb(33,13) =     1.152786d0 !     Arsenic -     Aluminum
      xfac(33,13) =     1.003580d0 !     Arsenic -     Aluminum
      alpb(33,14) =     1.915600d0 !     Arsenic -      Silicon
      xfac(33,14) =     1.430706d0 !     Arsenic -      Silicon
      alpb(33,16) =     1.954368d0 !     Arsenic -       Sulfur
      xfac(33,16) =     1.033784d0 !     Arsenic -       Sulfur
      alpb(33,17) =     1.691070d0 !     Arsenic -     Chlorine
      xfac(33,17) =     0.454433d0 !     Arsenic -     Chlorine
      alpb(33,22) =     1.932911d0 !     Arsenic -     Titanium
      xfac(33,22) =     1.581317d0 !     Arsenic -     Titanium
      alpb(33,27) =     3.368140d0 !     Arsenic -       Cobalt
      xfac(33,27) =     1.675240d0 !     Arsenic -       Cobalt
      alpb(33,30) =     1.459130d0 !     Arsenic -         Zinc
      xfac(33,30) =     3.156571d0 !     Arsenic -         Zinc
      alpb(33,31) =     1.730977d0 !     Arsenic -      Gallium
      xfac(33,31) =     1.686298d0 !     Arsenic -      Gallium
      alpb(33,33) =     1.588264d0 !     Arsenic -      Arsenic
      xfac(33,33) =     0.737307d0 !     Arsenic -      Arsenic
 !
      alpb(34, 1) =     2.035068d0 !    Selenium -     Hydrogen
      xfac(34, 1) =     0.847998d0 !    Selenium -     Hydrogen
      alpb(34, 6) =     2.387118d0 !    Selenium -       Carbon
      xfac(34, 6) =     1.114787d0 !    Selenium -       Carbon
      alpb(34, 7) =     1.937764d0 !    Selenium -     Nitrogen
      xfac(34, 7) =     0.482840d0 !    Selenium -     Nitrogen
      alpb(34, 8) =     2.484263d0 !    Selenium -       Oxygen
      xfac(34, 8) =     0.955161d0 !    Selenium -       Oxygen
      alpb(34, 9) =     2.302180d0 !    Selenium -     Fluorine
      xfac(34, 9) =     0.444806d0 !    Selenium -     Fluorine
      alpb(34,14) =     1.529817d0 !    Selenium -      Silicon
      xfac(34,14) =     0.518227d0 !    Selenium -      Silicon
      alpb(34,15) =     1.048183d0 !    Selenium -   Phosphorus
      xfac(34,15) =     0.292052d0 !    Selenium -   Phosphorus
      alpb(34,16) =     1.479606d0 !    Selenium -       Sulfur
      xfac(34,16) =     0.391721d0 !    Selenium -       Sulfur
      alpb(34,17) =     2.128861d0 !    Selenium -     Chlorine
      xfac(34,17) =     0.981067d0 !    Selenium -     Chlorine
      alpb(34,25) =     2.648038d0 !    Selenium -    Manganese
      xfac(34,25) =     2.180720d0 !    Selenium -    Manganese
      alpb(34,27) =     2.523450d0 !    Selenium -       Cobalt
      xfac(34,27) =     2.202410d0 !    Selenium -       Cobalt
      alpb(34,30) =     1.186242d0 !    Selenium -         Zinc
      xfac(34,30) =     0.511594d0 !    Selenium -         Zinc
      alpb(34,32) =     2.669057d0 !    Selenium -    Germanium
      xfac(34,32) =     5.872051d0 !    Selenium -    Germanium
      alpb(34,33) =     1.665280d0 !    Selenium -      Arsenic
      xfac(34,33) =     0.711261d0 !    Selenium -      Arsenic
      alpb(34,34) =     1.795894d0 !    Selenium -     Selenium
      xfac(34,34) =     0.821823d0 !    Selenium -     Selenium
 !
      alpb(35, 1) =     2.192803d0 !     Bromine -     Hydrogen
      xfac(35, 1) =     0.850378d0 !     Bromine -     Hydrogen
      alpb(35, 2) =     2.128275d0 !     Bromine -       Helium
      xfac(35, 2) =     1.062043d0 !     Bromine -       Helium
      alpb(35, 3) =     2.074441d0 !     Bromine -      Lithium
      xfac(35, 3) =     1.858866d0 !     Bromine -      Lithium
      alpb(35, 4) =     2.367146d0 !     Bromine -    Beryllium
      xfac(35, 4) =     1.940933d0 !     Bromine -    Beryllium
      alpb(35, 5) =     2.307890d0 !     Bromine -        Boron
      xfac(35, 5) =     1.226420d0 !     Bromine -        Boron
      alpb(35, 6) =     2.015086d0 !     Bromine -       Carbon
      xfac(35, 6) =     0.570686d0 !     Bromine -       Carbon
      alpb(35, 7) =     4.224901d0 !     Bromine -     Nitrogen
      xfac(35, 7) =    30.000133d0 !     Bromine -     Nitrogen
      alpb(35, 8) =     2.283046d0 !     Bromine -       Oxygen
      xfac(35, 8) =     0.706584d0 !     Bromine -       Oxygen
      alpb(35, 9) =     2.031765d0 !     Bromine -     Fluorine
      xfac(35, 9) =     0.293500d0 !     Bromine -     Fluorine
      alpb(35,10) =     2.464172d0 !     Bromine -         Neon
      xfac(35,10) =     1.006159d0 !     Bromine -         Neon
      alpb(35,11) =     1.622218d0 !     Bromine -       Sodium
      xfac(35,11) =     1.752937d0 !     Bromine -       Sodium
      alpb(35,12) =     2.195697d0 !     Bromine -    Magnesium
      xfac(35,12) =     2.916280d0 !     Bromine -    Magnesium
      alpb(35,13) =     1.894141d0 !     Bromine -     Aluminum
      xfac(35,13) =     2.357130d0 !     Bromine -     Aluminum
      alpb(35,14) =     1.570825d0 !     Bromine -      Silicon
      xfac(35,14) =     0.589511d0 !     Bromine -      Silicon
      alpb(35,15) =     1.402139d0 !     Bromine -   Phosphorus
      xfac(35,15) =     0.456521d0 !     Bromine -   Phosphorus
      alpb(35,16) =     1.509874d0 !     Bromine -       Sulfur
      xfac(35,16) =     0.286688d0 !     Bromine -       Sulfur
      alpb(35,17) =     1.710331d0 !     Bromine -     Chlorine
      xfac(35,17) =     0.389238d0 !     Bromine -     Chlorine
      alpb(35,18) =     2.450801d0 !     Bromine -        Argon
      xfac(35,18) =     3.262668d0 !     Bromine -        Argon
      alpb(35,19) =     1.616093d0 !     Bromine -    Potassium
      xfac(35,19) =     3.322795d0 !     Bromine -    Potassium
      alpb(35,20) =     2.078405d0 !     Bromine -      Calcium
      xfac(35,20) =     4.052910d0 !     Bromine -      Calcium
      alpb(35,21) =     1.793486d0 !     Bromine -     Scandium
      xfac(35,21) =     2.098251d0 !     Bromine -     Scandium
      alpb(35,22) =     1.674847d0 !     Bromine -     Titanium
      xfac(35,22) =     0.883434d0 !     Bromine -     Titanium
      alpb(35,23) =     1.902904d0 !     Bromine -     Vanadium
      xfac(35,23) =     0.612698d0 !     Bromine -     Vanadium
      alpb(35,24) =     1.566028d0 !     Bromine -     Chromium
      xfac(35,24) =     0.217853d0 !     Bromine -     Chromium
      alpb(35,25) =     2.283820d0 !     Bromine -    Manganese
      xfac(35,25) =     1.183580d0 !     Bromine -    Manganese
      alpb(35,26) =     3.641782d0 !     Bromine -         Iron
      xfac(35,26) =     6.061921d0 !     Bromine -         Iron
      alpb(35,27) =     2.632688d0 !     Bromine -       Cobalt
      xfac(35,27) =     0.425148d0 !     Bromine -       Cobalt
      alpb(35,28) =     2.772136d0 !     Bromine -       Nickel
      xfac(35,28) =     0.632145d0 !     Bromine -       Nickel
      alpb(35,29) =     5.826407d0 !     Bromine -       Copper
      xfac(35,29) =     0.768517d0 !     Bromine -       Copper
      alpb(35,30) =     1.416120d0 !     Bromine -         Zinc
      xfac(35,30) =     0.747027d0 !     Bromine -         Zinc
      alpb(35,31) =     1.819105d0 !     Bromine -      Gallium
      xfac(35,31) =     1.261036d0 !     Bromine -      Gallium
      alpb(35,32) =     1.602366d0 !     Bromine -    Germanium
      xfac(35,32) =     0.627737d0 !     Bromine -    Germanium
      alpb(35,33) =     1.520170d0 !     Bromine -      Arsenic
      xfac(35,33) =     0.514153d0 !     Bromine -      Arsenic
      alpb(35,34) =     1.483713d0 !     Bromine -     Selenium
      xfac(35,34) =     0.319342d0 !     Bromine -     Selenium
      alpb(35,35) =     1.758146d0 !     Bromine -      Bromine
      xfac(35,35) =     0.615308d0 !     Bromine -      Bromine
 !
      alpb(36, 1) =     3.770453d0 !     Krypton -     Hydrogen
      xfac(36, 1) =     5.125897d0 !     Krypton -     Hydrogen
      alpb(36, 2) =     1.996943d0 !     Krypton -       Helium
      xfac(36, 2) =     0.627701d0 !     Krypton -       Helium
      alpb(36, 3) =     3.314562d0 !     Krypton -      Lithium
      xfac(36, 3) =     8.758697d0 !     Krypton -      Lithium
      alpb(36, 4) =     3.253048d0 !     Krypton -    Beryllium
      xfac(36, 4) =    10.237796d0 !     Krypton -    Beryllium
      alpb(36, 5) =     2.363169d0 !     Krypton -        Boron
      xfac(36, 5) =     2.946781d0 !     Krypton -        Boron
      alpb(36, 6) =     2.076738d0 !     Krypton -       Carbon
      xfac(36, 6) =     0.652623d0 !     Krypton -       Carbon
      alpb(36, 7) =     1.644052d0 !     Krypton -     Nitrogen
      xfac(36, 7) =     0.199606d0 !     Krypton -     Nitrogen
      alpb(36, 8) =     0.292300d0 !     Krypton -       Oxygen
      xfac(36, 8) =     0.006733d0 !     Krypton -       Oxygen
      alpb(36, 9) =     3.452321d0 !     Krypton -     Fluorine
      xfac(36, 9) =     4.134407d0 !     Krypton -     Fluorine
      alpb(36,10) =     2.813679d0 !     Krypton -         Neon
      xfac(36,10) =     1.433722d0 !     Krypton -         Neon
      alpb(36,11) =     2.480598d0 !     Krypton -       Sodium
      xfac(36,11) =     8.354448d0 !     Krypton -       Sodium
      alpb(36,12) =     1.391487d0 !     Krypton -    Magnesium
      xfac(36,12) =     0.888436d0 !     Krypton -    Magnesium
      alpb(36,13) =     2.467131d0 !     Krypton -     Aluminum
      xfac(36,13) =     5.091716d0 !     Krypton -     Aluminum
      alpb(36,14) =     1.764100d0 !     Krypton -      Silicon
      xfac(36,14) =     0.554250d0 !     Krypton -      Silicon
      alpb(36,17) =     1.884974d0 !     Krypton -     Chlorine
      xfac(36,17) =     0.520217d0 !     Krypton -     Chlorine
      alpb(36,18) =     1.995125d0 !     Krypton -        Argon
      xfac(36,18) =     0.554874d0 !     Krypton -        Argon
      alpb(36,19) =     2.182487d0 !     Krypton -    Potassium
      xfac(36,19) =     8.609782d0 !     Krypton -    Potassium
      alpb(36,20) =     1.305197d0 !     Krypton -      Calcium
      xfac(36,20) =     0.878891d0 !     Krypton -      Calcium
      alpb(36,35) =     1.529006d0 !     Krypton -      Bromine
      xfac(36,35) =     0.308098d0 !     Krypton -      Bromine
      alpb(36,36) =     1.135319d0 !     Krypton -      Krypton
      xfac(36,36) =     0.052099d0 !     Krypton -      Krypton
 !
      alpb(37, 1) =     2.443556d0 !    Rubidium -     Hydrogen
      xfac(37, 1) =    29.861632d0 !    Rubidium -     Hydrogen
      alpb(37, 2) =     1.270741d0 !    Rubidium -       Helium
      xfac(37, 2) =     1.862585d0 !    Rubidium -       Helium
      alpb(37, 5) =     5.532239d0 !    Rubidium -        Boron
      xfac(37, 5) =     9.040493d0 !    Rubidium -        Boron
      alpb(37, 6) =     2.765830d0 !    Rubidium -       Carbon
      xfac(37, 6) =    29.974031d0 !    Rubidium -       Carbon
      alpb(37, 7) =     0.761047d0 !    Rubidium -     Nitrogen
      xfac(37, 7) =     0.024636d0 !    Rubidium -     Nitrogen
      alpb(37, 8) =     1.334908d0 !    Rubidium -       Oxygen
      xfac(37, 8) =     1.125350d0 !    Rubidium -       Oxygen
      alpb(37, 9) =     3.638122d0 !    Rubidium -     Fluorine
      xfac(37, 9) =    28.815278d0 !    Rubidium -     Fluorine
      alpb(37,10) =     2.267591d0 !    Rubidium -         Neon
      xfac(37,10) =     7.736563d0 !    Rubidium -         Neon
      alpb(37,13) =     0.798774d0 !    Rubidium -     Aluminum
      xfac(37,13) =     2.992457d0 !    Rubidium -     Aluminum
      alpb(37,16) =     1.303184d0 !    Rubidium -       Sulfur
      xfac(37,16) =     0.964411d0 !    Rubidium -       Sulfur
      alpb(37,17) =     2.274411d0 !    Rubidium -     Chlorine
      xfac(37,17) =    10.384486d0 !    Rubidium -     Chlorine
      alpb(37,18) =     2.510977d0 !    Rubidium -        Argon
      xfac(37,18) =    18.433329d0 !    Rubidium -        Argon
      alpb(37,35) =     1.797766d0 !    Rubidium -      Bromine
      xfac(37,35) =     5.176214d0 !    Rubidium -      Bromine
      alpb(37,36) =     2.268753d0 !    Rubidium -      Krypton
      xfac(37,36) =    15.307503d0 !    Rubidium -      Krypton
      alpb(37,37) =     1.180818d0 !    Rubidium -     Rubidium
      xfac(37,37) =    20.147610d0 !    Rubidium -     Rubidium
 !
      alpb(38, 1) =     2.105914d0 !   Strontium -     Hydrogen
      xfac(38, 1) =    12.973316d0 !   Strontium -     Hydrogen
      alpb(38, 6) =     1.986688d0 !   Strontium -       Carbon
      xfac(38, 6) =     6.654657d0 !   Strontium -       Carbon
      alpb(38, 7) =     2.183629d0 !   Strontium -     Nitrogen
      xfac(38, 7) =     6.853866d0 !   Strontium -     Nitrogen
      alpb(38, 8) =     2.138399d0 !   Strontium -       Oxygen
      xfac(38, 8) =     3.561396d0 !   Strontium -       Oxygen
      alpb(38, 9) =     3.050666d0 !   Strontium -     Fluorine
      xfac(38, 9) =    10.971705d0 !   Strontium -     Fluorine
      alpb(38,14) =     2.969780d0 !   Strontium -      Silicon
      xfac(38,14) =     2.764750d0 !   Strontium -      Silicon
      alpb(38,15) =     2.789150d0 !   Strontium -   Phosphorus
      xfac(38,15) =     2.552100d0 !   Strontium -   Phosphorus
      alpb(38,16) =     1.598106d0 !   Strontium -       Sulfur
      xfac(38,16) =     3.129603d0 !   Strontium -       Sulfur
      alpb(38,17) =     1.854190d0 !   Strontium -     Chlorine
      xfac(38,17) =     3.783955d0 !   Strontium -     Chlorine
      alpb(38,22) =     2.880030d0 !   Strontium -     Titanium
      xfac(38,22) =     2.817250d0 !   Strontium -     Titanium
      alpb(38,35) =     1.524316d0 !   Strontium -      Bromine
      xfac(38,35) =     2.766567d0 !   Strontium -      Bromine
      alpb(38,38) =     1.000040d0 !   Strontium -    Strontium
      xfac(38,38) =     5.372120d0 !   Strontium -    Strontium
 !
      alpb(39, 1) =     1.189053d0 !     Yttrium -     Hydrogen
      xfac(39, 1) =     0.612399d0 !     Yttrium -     Hydrogen
      alpb(39, 6) =     1.336094d0 !     Yttrium -       Carbon
      xfac(39, 6) =     0.504306d0 !     Yttrium -       Carbon
      alpb(39, 7) =     1.778796d0 !     Yttrium -     Nitrogen
      xfac(39, 7) =     1.627903d0 !     Yttrium -     Nitrogen
      alpb(39, 8) =     1.851030d0 !     Yttrium -       Oxygen
      xfac(39, 8) =     1.742922d0 !     Yttrium -       Oxygen
      alpb(39, 9) =     2.648046d0 !     Yttrium -     Fluorine
      xfac(39, 9) =     4.433809d0 !     Yttrium -     Fluorine
      alpb(39,13) =     1.003500d0 !     Yttrium -     Aluminum
      xfac(39,13) =     0.500670d0 !     Yttrium -     Aluminum
      alpb(39,14) =     2.016820d0 !     Yttrium -      Silicon
      xfac(39,14) =     3.219030d0 !     Yttrium -      Silicon
      alpb(39,15) =     0.954450d0 !     Yttrium -   Phosphorus
      xfac(39,15) =     0.541660d0 !     Yttrium -   Phosphorus
      alpb(39,16) =     0.971688d0 !     Yttrium -       Sulfur
      xfac(39,16) =     0.318222d0 !     Yttrium -       Sulfur
      alpb(39,17) =     1.630152d0 !     Yttrium -     Chlorine
      xfac(39,17) =     1.154959d0 !     Yttrium -     Chlorine
      alpb(39,35) =     1.401208d0 !     Yttrium -      Bromine
      xfac(39,35) =     1.054316d0 !     Yttrium -      Bromine
      alpb(39,39) =     1.012681d0 !     Yttrium -      Yttrium
      xfac(39,39) =     1.691725d0 !     Yttrium -      Yttrium
 !
      alpb(40, 1) =     1.379703d0 !   Zirconium -     Hydrogen
      xfac(40, 1) =     0.593732d0 !   Zirconium -     Hydrogen
      alpb(40, 6) =     2.029427d0 !   Zirconium -       Carbon
      xfac(40, 6) =     1.999182d0 !   Zirconium -       Carbon
      alpb(40, 7) =     1.707083d0 !   Zirconium -     Nitrogen
      xfac(40, 7) =     0.995045d0 !   Zirconium -     Nitrogen
      alpb(40, 8) =     1.709570d0 !   Zirconium -       Oxygen
      xfac(40, 8) =     1.057525d0 !   Zirconium -       Oxygen
      alpb(40, 9) =     1.900925d0 !   Zirconium -     Fluorine
      xfac(40, 9) =     0.861142d0 !   Zirconium -     Fluorine
      alpb(40,13) =     1.270620d0 !   Zirconium -     Aluminum
      xfac(40,13) =     0.874060d0 !   Zirconium -     Aluminum
      alpb(40,14) =     1.750833d0 !   Zirconium -      Silicon
      xfac(40,14) =     1.723343d0 !   Zirconium -      Silicon
      alpb(40,15) =     1.091858d0 !   Zirconium -   Phosphorus
      xfac(40,15) =     0.748376d0 !   Zirconium -   Phosphorus
      alpb(40,16) =     2.129761d0 !   Zirconium -       Sulfur
      xfac(40,16) =     2.429324d0 !   Zirconium -       Sulfur
      alpb(40,17) =     1.328835d0 !   Zirconium -     Chlorine
      xfac(40,17) =     0.443099d0 !   Zirconium -     Chlorine
      alpb(40,35) =     1.446868d0 !   Zirconium -      Bromine
      xfac(40,35) =     0.858909d0 !   Zirconium -      Bromine
      alpb(40,40) =     3.865968d0 !   Zirconium -    Zirconium
      xfac(40,40) =     3.077773d0 !   Zirconium -    Zirconium
 !
      alpb(41, 1) =     2.505912d0 !     Niobium -     Hydrogen
      xfac(41, 1) =     3.603779d0 !     Niobium -     Hydrogen
      alpb(41, 6) =     2.621012d0 !     Niobium -       Carbon
      xfac(41, 6) =     4.575481d0 !     Niobium -       Carbon
      alpb(41, 7) =     2.023863d0 !     Niobium -     Nitrogen
      xfac(41, 7) =     1.213587d0 !     Niobium -     Nitrogen
      alpb(41, 8) =     2.049489d0 !     Niobium -       Oxygen
      xfac(41, 8) =     1.184719d0 !     Niobium -       Oxygen
      alpb(41, 9) =     3.003157d0 !     Niobium -     Fluorine
      xfac(41, 9) =     3.663682d0 !     Niobium -     Fluorine
      alpb(41,11) =     2.551010d0 !     Niobium -       Sodium
      xfac(41,11) =     8.276020d0 !     Niobium -       Sodium
      alpb(41,15) =     2.221608d0 !     Niobium -   Phosphorus
      xfac(41,15) =     6.201507d0 !     Niobium -   Phosphorus
      alpb(41,16) =     2.249482d0 !     Niobium -       Sulfur
      xfac(41,16) =     2.460020d0 !     Niobium -       Sulfur
      alpb(41,17) =     2.215275d0 !     Niobium -     Chlorine
      xfac(41,17) =     1.891557d0 !     Niobium -     Chlorine
      alpb(41,19) =     4.521360d0 !     Niobium -    Potassium
      xfac(41,19) =     2.026590d0 !     Niobium -    Potassium
      alpb(41,35) =     2.006678d0 !     Niobium -      Bromine
      xfac(41,35) =     1.921269d0 !     Niobium -      Bromine
      alpb(41,41) =     1.727941d0 !     Niobium -      Niobium
      xfac(41,41) =     2.122388d0 !     Niobium -      Niobium
 !
      alpb(42, 1) =     2.035748d0 !  Molybdenum -     Hydrogen
      xfac(42, 1) =     0.934686d0 !  Molybdenum -     Hydrogen
      alpb(42, 6) =     2.198672d0 !  Molybdenum -       Carbon
      xfac(42, 6) =     1.190742d0 !  Molybdenum -       Carbon
      alpb(42, 7) =     1.869475d0 !  Molybdenum -     Nitrogen
      xfac(42, 7) =     0.608268d0 !  Molybdenum -     Nitrogen
      alpb(42, 8) =     1.755424d0 !  Molybdenum -       Oxygen
      xfac(42, 8) =     0.511267d0 !  Molybdenum -       Oxygen
      alpb(42, 9) =     2.202593d0 !  Molybdenum -     Fluorine
      xfac(42, 9) =     0.610429d0 !  Molybdenum -     Fluorine
      alpb(42,11) =     2.440770d0 !  Molybdenum -       Sodium
      xfac(42,11) =     8.286550d0 !  Molybdenum -       Sodium
      alpb(42,15) =     1.850441d0 !  Molybdenum -   Phosphorus
      xfac(42,15) =     1.522846d0 !  Molybdenum -   Phosphorus
      alpb(42,16) =     1.939658d0 !  Molybdenum -       Sulfur
      xfac(42,16) =     0.830428d0 !  Molybdenum -       Sulfur
      alpb(42,17) =     1.783362d0 !  Molybdenum -     Chlorine
      xfac(42,17) =     0.474325d0 !  Molybdenum -     Chlorine
      alpb(42,19) =     3.939420d0 !  Molybdenum -    Potassium
      xfac(42,19) =     2.142390d0 !  Molybdenum -    Potassium
      alpb(42,24) =     2.674616d0 !  Molybdenum -     Chromium
      xfac(42,24) =     1.741943d0 !  Molybdenum -     Chromium
      alpb(42,35) =     1.283334d0 !  Molybdenum -      Bromine
      xfac(42,35) =     0.225918d0 !  Molybdenum -      Bromine
      alpb(42,42) =     2.034254d0 !  Molybdenum -   Molybdenum
      xfac(42,42) =     0.626462d0 !  Molybdenum -   Molybdenum
 !
      alpb(43, 1) =     2.830345d0 !  Technetium -     Hydrogen
      xfac(43, 1) =     6.310334d0 !  Technetium -     Hydrogen
      alpb(43, 6) =     3.198326d0 !  Technetium -       Carbon
      xfac(43, 6) =     3.972439d0 !  Technetium -       Carbon
      alpb(43, 7) =     2.315417d0 !  Technetium -     Nitrogen
      xfac(43, 7) =     0.727130d0 !  Technetium -     Nitrogen
      alpb(43, 8) =     2.405190d0 !  Technetium -       Oxygen
      xfac(43, 8) =     1.024616d0 !  Technetium -       Oxygen
      alpb(43, 9) =     3.604815d0 !  Technetium -     Fluorine
      xfac(43, 9) =     5.811784d0 !  Technetium -     Fluorine
      alpb(43,16) =     2.463401d0 !  Technetium -       Sulfur
      xfac(43,16) =     1.496502d0 !  Technetium -       Sulfur
      alpb(43,17) =     2.572043d0 !  Technetium -     Chlorine
      xfac(43,17) =     1.651583d0 !  Technetium -     Chlorine
      alpb(43,32) =     2.852820d0 !  Technetium -    Germanium
      xfac(43,32) =     2.152060d0 !  Technetium -    Germanium
      alpb(43,34) =     2.523660d0 !  Technetium -     Selenium
      xfac(43,34) =     2.202620d0 !  Technetium -     Selenium
      alpb(43,35) =     2.828264d0 !  Technetium -      Bromine
      xfac(43,35) =     3.820130d0 !  Technetium -      Bromine
 !
      alpb(44, 1) =     2.892899d0 !   Ruthenium -     Hydrogen
      xfac(44, 1) =     7.137976d0 !   Ruthenium -     Hydrogen
      alpb(44, 6) =     2.784833d0 !   Ruthenium -       Carbon
      xfac(44, 6) =     1.134936d0 !   Ruthenium -       Carbon
      alpb(44, 7) =     3.055504d0 !   Ruthenium -     Nitrogen
      xfac(44, 7) =     2.334094d0 !   Ruthenium -     Nitrogen
      alpb(44, 8) =     3.134940d0 !   Ruthenium -       Oxygen
      xfac(44, 8) =     2.976279d0 !   Ruthenium -       Oxygen
      alpb(44, 9) =     3.878711d0 !   Ruthenium -     Fluorine
      xfac(44, 9) =     6.947128d0 !   Ruthenium -     Fluorine
      alpb(44,14) =     2.775910d0 !   Ruthenium -      Silicon
      xfac(44,14) =     0.849430d0 !   Ruthenium -      Silicon
      alpb(44,15) =     0.298916d0 !   Ruthenium -   Phosphorus
      xfac(44,15) =     0.056974d0 !   Ruthenium -   Phosphorus
      alpb(44,16) =     2.508076d0 !   Ruthenium -       Sulfur
      xfac(44,16) =     1.006683d0 !   Ruthenium -       Sulfur
      alpb(44,17) =     1.759883d0 !   Ruthenium -     Chlorine
      xfac(44,17) =     0.126586d0 !   Ruthenium -     Chlorine
      alpb(44,32) =     2.852320d0 !   Ruthenium -    Germanium
      xfac(44,32) =     2.151560d0 !   Ruthenium -    Germanium
      alpb(44,34) =     2.523160d0 !   Ruthenium -     Selenium
      xfac(44,34) =     2.202120d0 !   Ruthenium -     Selenium
      alpb(44,35) =     2.584735d0 !   Ruthenium -      Bromine
      xfac(44,35) =     0.659881d0 !   Ruthenium -      Bromine
      alpb(44,44) =     0.572056d0 !   Ruthenium -    Ruthenium
      xfac(44,44) =     0.097805d0 !   Ruthenium -    Ruthenium
 !
      alpb(45, 1) =     3.104165d0 !     Rhodium -     Hydrogen
      xfac(45, 1) =     2.306107d0 !     Rhodium -     Hydrogen
      alpb(45, 6) =     3.415991d0 !     Rhodium -       Carbon
      xfac(45, 6) =     3.488079d0 !     Rhodium -       Carbon
      alpb(45, 7) =     3.585462d0 !     Rhodium -     Nitrogen
      xfac(45, 7) =     4.000947d0 !     Rhodium -     Nitrogen
      alpb(45, 8) =     3.927830d0 !     Rhodium -       Oxygen
      xfac(45, 8) =    10.298676d0 !     Rhodium -       Oxygen
      alpb(45, 9) =     4.051654d0 !     Rhodium -     Fluorine
      xfac(45, 9) =     9.065384d0 !     Rhodium -     Fluorine
      alpb(45,14) =     2.776490d0 !     Rhodium -      Silicon
      xfac(45,14) =     0.850010d0 !     Rhodium -      Silicon
      alpb(45,15) =     2.334607d0 !     Rhodium -   Phosphorus
      xfac(45,15) =     1.038141d0 !     Rhodium -   Phosphorus
      alpb(45,16) =     3.154006d0 !     Rhodium -       Sulfur
      xfac(45,16) =     4.816410d0 !     Rhodium -       Sulfur
      alpb(45,17) =     3.300130d0 !     Rhodium -     Chlorine
      xfac(45,17) =     3.586865d0 !     Rhodium -     Chlorine
      alpb(45,32) =     2.852900d0 !     Rhodium -    Germanium
      xfac(45,32) =     2.152140d0 !     Rhodium -    Germanium
      alpb(45,34) =     2.523740d0 !     Rhodium -     Selenium
      xfac(45,34) =     2.202700d0 !     Rhodium -     Selenium
      alpb(45,35) =     2.928082d0 !     Rhodium -      Bromine
      xfac(45,35) =     1.510149d0 !     Rhodium -      Bromine
      alpb(45,45) =     2.497328d0 !     Rhodium -      Rhodium
      xfac(45,45) =     2.070114d0 !     Rhodium -      Rhodium
 !
      alpb(46, 1) =     2.183761d0 !   Palladium -     Hydrogen
      xfac(46, 1) =     0.443269d0 !   Palladium -     Hydrogen
      alpb(46, 6) =     4.777192d0 !   Palladium -       Carbon
      xfac(46, 6) =     9.853715d0 !   Palladium -       Carbon
      alpb(46, 7) =     2.328046d0 !   Palladium -     Nitrogen
      xfac(46, 7) =     0.249703d0 !   Palladium -     Nitrogen
      alpb(46, 8) =     2.154867d0 !   Palladium -       Oxygen
      xfac(46, 8) =     0.216403d0 !   Palladium -       Oxygen
      alpb(46, 9) =     4.237312d0 !   Palladium -     Fluorine
      xfac(46, 9) =     6.945312d0 !   Palladium -     Fluorine
      alpb(46,13) =     1.572720d0 !   Palladium -     Aluminum
      xfac(46,13) =     1.057290d0 !   Palladium -     Aluminum
      alpb(46,14) =     2.948200d0 !   Palladium -      Silicon
      xfac(46,14) =     2.225104d0 !   Palladium -      Silicon
      alpb(46,15) =     0.803630d0 !   Palladium -   Phosphorus
      xfac(46,15) =     0.045017d0 !   Palladium -   Phosphorus
      alpb(46,16) =     2.177801d0 !   Palladium -       Sulfur
      xfac(46,16) =     0.255229d0 !   Palladium -       Sulfur
      alpb(46,17) =     3.871243d0 !   Palladium -     Chlorine
      xfac(46,17) =     2.969891d0 !   Palladium -     Chlorine
      alpb(46,35) =     5.994879d0 !   Palladium -      Bromine
      xfac(46,35) =     4.638051d0 !   Palladium -      Bromine
      alpb(46,46) =     1.064375d0 !   Palladium -    Palladium
      xfac(46,46) =     0.051956d0 !   Palladium -    Palladium
 !
      alpb(47, 1) =     2.895936d0 !      Silver -     Hydrogen
      xfac(47, 1) =     1.995168d0 !      Silver -     Hydrogen
      alpb(47, 6) =     4.404336d0 !      Silver -       Carbon
      xfac(47, 6) =    11.335456d0 !      Silver -       Carbon
      alpb(47, 7) =     4.659871d0 !      Silver -     Nitrogen
      xfac(47, 7) =    19.803710d0 !      Silver -     Nitrogen
      alpb(47, 8) =     1.893874d0 !      Silver -       Oxygen
      xfac(47, 8) =     0.165661d0 !      Silver -       Oxygen
      alpb(47, 9) =     4.628423d0 !      Silver -     Fluorine
      xfac(47, 9) =    12.695884d0 !      Silver -     Fluorine
      alpb(47,13) =     1.928800d0 !      Silver -     Aluminum
      xfac(47,13) =     0.896514d0 !      Silver -     Aluminum
      alpb(47,15) =     6.000006d0 !      Silver -   Phosphorus
      xfac(47,15) =     0.049932d0 !      Silver -   Phosphorus
      alpb(47,16) =     3.653121d0 !      Silver -       Sulfur
      xfac(47,16) =    11.188022d0 !      Silver -       Sulfur
      alpb(47,17) =     4.441176d0 !      Silver -     Chlorine
      xfac(47,17) =    23.765459d0 !      Silver -     Chlorine
      alpb(47,35) =     3.677491d0 !      Silver -      Bromine
      xfac(47,35) =     1.714369d0 !      Silver -      Bromine
      alpb(47,47) =     2.127645d0 !      Silver -       Silver
      xfac(47,47) =     0.557742d0 !      Silver -       Silver
 !
      alpb(48, 1) =     2.628748d0 !     Cadmium -     Hydrogen
      xfac(48, 1) =    11.914201d0 !     Cadmium -     Hydrogen
      alpb(48, 6) =     1.425678d0 !     Cadmium -       Carbon
      xfac(48, 6) =     0.603441d0 !     Cadmium -       Carbon
      alpb(48, 7) =     0.970423d0 !     Cadmium -     Nitrogen
      xfac(48, 7) =     0.180663d0 !     Cadmium -     Nitrogen
      alpb(48, 8) =     1.696673d0 !     Cadmium -       Oxygen
      xfac(48, 8) =     0.926146d0 !     Cadmium -       Oxygen
      alpb(48, 9) =     2.312135d0 !     Cadmium -     Fluorine
      xfac(48, 9) =     1.353665d0 !     Cadmium -     Fluorine
      alpb(48,14) =     1.371225d0 !     Cadmium -      Silicon
      xfac(48,14) =     2.253346d0 !     Cadmium -      Silicon
      alpb(48,16) =     1.182202d0 !     Cadmium -       Sulfur
      xfac(48,16) =     0.361389d0 !     Cadmium -       Sulfur
      alpb(48,17) =     0.943547d0 !     Cadmium -     Chlorine
      xfac(48,17) =     0.140424d0 !     Cadmium -     Chlorine
      alpb(48,35) =     1.001451d0 !     Cadmium -      Bromine
      xfac(48,35) =     0.272267d0 !     Cadmium -      Bromine
      alpb(48,48) =     1.564044d0 !     Cadmium -      Cadmium
      xfac(48,48) =    18.617999d0 !     Cadmium -      Cadmium
 !
      alpb(49, 1) =     3.064144d0 !      Indium -     Hydrogen
      xfac(49, 1) =    14.975293d0 !      Indium -     Hydrogen
      alpb(49, 6) =     2.189272d0 !      Indium -       Carbon
      xfac(49, 6) =     2.187385d0 !      Indium -       Carbon
      alpb(49, 7) =     2.469868d0 !      Indium -     Nitrogen
      xfac(49, 7) =     3.369993d0 !      Indium -     Nitrogen
      alpb(49, 8) =     2.662095d0 !      Indium -       Oxygen
      xfac(49, 8) =     4.128583d0 !      Indium -       Oxygen
      alpb(49, 9) =     2.948797d0 !      Indium -     Fluorine
      xfac(49, 9) =     3.701016d0 !      Indium -     Fluorine
      alpb(49,16) =     2.542131d0 !      Indium -       Sulfur
      xfac(49,16) =     6.341105d0 !      Indium -       Sulfur
      alpb(49,17) =     2.233405d0 !      Indium -     Chlorine
      xfac(49,17) =     2.388552d0 !      Indium -     Chlorine
      alpb(49,31) =     1.628870d0 !      Indium -      Gallium
      xfac(49,31) =     2.421987d0 !      Indium -      Gallium
      alpb(49,33) =     2.299552d0 !      Indium -      Arsenic
      xfac(49,33) =     6.208350d0 !      Indium -      Arsenic
      alpb(49,34) =     1.906572d0 !      Indium -     Selenium
      xfac(49,34) =     2.319323d0 !      Indium -     Selenium
      alpb(49,35) =     2.257957d0 !      Indium -      Bromine
      xfac(49,35) =     3.728598d0 !      Indium -      Bromine
      alpb(49,49) =     2.073241d0 !      Indium -       Indium
      xfac(49,49) =     8.063491d0 !      Indium -       Indium
 !
      alpb(50, 1) =     2.648910d0 !         Tin -     Hydrogen
      xfac(50, 1) =     6.535162d0 !         Tin -     Hydrogen
      alpb(50, 6) =     2.440538d0 !         Tin -       Carbon
      xfac(50, 6) =     3.374355d0 !         Tin -       Carbon
      alpb(50, 7) =     2.085589d0 !         Tin -     Nitrogen
      xfac(50, 7) =     1.391900d0 !         Tin -     Nitrogen
      alpb(50, 8) =     2.727260d0 !         Tin -       Oxygen
      xfac(50, 8) =     4.374017d0 !         Tin -       Oxygen
      alpb(50, 9) =     3.724286d0 !         Tin -     Fluorine
      xfac(50, 9) =    18.598664d0 !         Tin -     Fluorine
      alpb(50,16) =     2.131542d0 !         Tin -       Sulfur
      xfac(50,16) =     2.314870d0 !         Tin -       Sulfur
      alpb(50,17) =     1.771522d0 !         Tin -     Chlorine
      xfac(50,17) =     0.807782d0 !         Tin -     Chlorine
      alpb(50,32) =     2.524633d0 !         Tin -    Germanium
      xfac(50,32) =    12.343411d0 !         Tin -    Germanium
      alpb(50,34) =     2.127377d0 !         Tin -     Selenium
      xfac(50,34) =     3.061885d0 !         Tin -     Selenium
      alpb(50,35) =     1.535089d0 !         Tin -      Bromine
      xfac(50,35) =     0.668798d0 !         Tin -      Bromine
      alpb(50,50) =     0.921000d0 !         Tin -          Tin
      xfac(50,50) =     0.287000d0 !         Tin -          Tin
 !
      alpb(51, 1) =     1.571272d0 !    Antimony -     Hydrogen
      xfac(51, 1) =     0.795343d0 !    Antimony -     Hydrogen
      alpb(51, 6) =     1.696206d0 !    Antimony -       Carbon
      xfac(51, 6) =     0.579212d0 !    Antimony -       Carbon
      alpb(51, 7) =     0.676115d0 !    Antimony -     Nitrogen
      xfac(51, 7) =     0.082065d0 !    Antimony -     Nitrogen
      alpb(51, 8) =     1.846384d0 !    Antimony -       Oxygen
      xfac(51, 8) =     0.634234d0 !    Antimony -       Oxygen
      alpb(51, 9) =     2.182922d0 !    Antimony -     Fluorine
      xfac(51, 9) =     0.650277d0 !    Antimony -     Fluorine
      alpb(51,13) =     1.422641d0 !    Antimony -     Aluminum
      xfac(51,13) =     1.616690d0 !    Antimony -     Aluminum
      alpb(51,14) =     2.686590d0 !    Antimony -      Silicon
      xfac(51,14) =     8.713749d0 !    Antimony -      Silicon
      alpb(51,16) =     1.418837d0 !    Antimony -       Sulfur
      xfac(51,16) =     0.396969d0 !    Antimony -       Sulfur
      alpb(51,17) =     1.117287d0 !    Antimony -     Chlorine
      xfac(51,17) =     0.156475d0 !    Antimony -     Chlorine
      alpb(51,25) =     2.400320d0 !    Antimony -    Manganese
      xfac(51,25) =     2.236710d0 !    Antimony -    Manganese
      alpb(51,27) =     2.204630d0 !    Antimony -       Cobalt
      xfac(51,27) =     2.276050d0 !    Antimony -       Cobalt
      alpb(51,35) =     1.063916d0 !    Antimony -      Bromine
      xfac(51,35) =     0.198044d0 !    Antimony -      Bromine
      alpb(51,43) =     2.204850d0 !    Antimony -   Technetium
      xfac(51,43) =     2.276260d0 !    Antimony -   Technetium
      alpb(51,44) =     2.204350d0 !    Antimony -    Ruthenium
      xfac(51,44) =     2.275760d0 !    Antimony -    Ruthenium
      alpb(51,45) =     2.204930d0 !    Antimony -      Rhodium
      xfac(51,45) =     2.276340d0 !    Antimony -      Rhodium
      alpb(51,49) =     2.141933d0 !    Antimony -       Indium
      xfac(51,49) =     6.660801d0 !    Antimony -       Indium
      alpb(51,51) =     1.348535d0 !    Antimony -     Antimony
      xfac(51,51) =     0.724885d0 !    Antimony -     Antimony
 !
      alpb(52, 1) =     2.039130d0 !   Tellurium -     Hydrogen
      xfac(52, 1) =     1.807679d0 !   Tellurium -     Hydrogen
      alpb(52, 6) =     1.992816d0 !   Tellurium -       Carbon
      xfac(52, 6) =     0.970494d0 !   Tellurium -       Carbon
      alpb(52, 7) =     1.722269d0 !   Tellurium -     Nitrogen
      xfac(52, 7) =     0.358593d0 !   Tellurium -     Nitrogen
      alpb(52, 8) =     1.853064d0 !   Tellurium -       Oxygen
      xfac(52, 8) =     0.382926d0 !   Tellurium -       Oxygen
      alpb(52, 9) =     1.998576d0 !   Tellurium -     Fluorine
      xfac(52, 9) =     0.200822d0 !   Tellurium -     Fluorine
      alpb(52,13) =     1.387541d0 !   Tellurium -     Aluminum
      xfac(52,13) =     2.106812d0 !   Tellurium -     Aluminum
      alpb(52,15) =     1.453718d0 !   Tellurium -   Phosphorus
      xfac(52,15) =     1.109289d0 !   Tellurium -   Phosphorus
      alpb(52,16) =     1.830170d0 !   Tellurium -       Sulfur
      xfac(52,16) =     0.943925d0 !   Tellurium -       Sulfur
      alpb(52,17) =     1.300260d0 !   Tellurium -     Chlorine
      xfac(52,17) =     0.285478d0 !   Tellurium -     Chlorine
      alpb(52,30) =     1.218929d0 !   Tellurium -         Zinc
      xfac(52,30) =     1.756070d0 !   Tellurium -         Zinc
      alpb(52,32) =     2.342372d0 !   Tellurium -    Germanium
      xfac(52,32) =     7.019049d0 !   Tellurium -    Germanium
      alpb(52,33) =     1.189253d0 !   Tellurium -      Arsenic
      xfac(52,33) =     0.685774d0 !   Tellurium -      Arsenic
      alpb(52,34) =     1.566008d0 !   Tellurium -     Selenium
      xfac(52,34) =     1.187826d0 !   Tellurium -     Selenium
      alpb(52,35) =     1.250940d0 !   Tellurium -      Bromine
      xfac(52,35) =     0.394202d0 !   Tellurium -      Bromine
      alpb(52,48) =     1.307262d0 !   Tellurium -      Cadmium
      xfac(52,48) =     1.085919d0 !   Tellurium -      Cadmium
      alpb(52,49) =     1.540988d0 !   Tellurium -       Indium
      xfac(52,49) =     2.039582d0 !   Tellurium -       Indium
      alpb(52,50) =     1.763941d0 !   Tellurium -          Tin
      xfac(52,50) =     2.951976d0 !   Tellurium -          Tin
      alpb(52,52) =     1.164978d0 !   Tellurium -    Tellurium
      xfac(52,52) =     0.642486d0 !   Tellurium -    Tellurium
 !
      alpb(53, 1) =     2.139913d0 !      Iodine -     Hydrogen
      xfac(53, 1) =     0.981898d0 !      Iodine -     Hydrogen
      alpb(53, 2) =     2.172984d0 !      Iodine -       Helium
      xfac(53, 2) =     1.630721d0 !      Iodine -       Helium
      alpb(53, 3) =     2.121251d0 !      Iodine -      Lithium
      xfac(53, 3) =     4.168599d0 !      Iodine -      Lithium
      alpb(53, 4) =     2.288023d0 !      Iodine -    Beryllium
      xfac(53, 4) =     2.351898d0 !      Iodine -    Beryllium
      alpb(53, 5) =     2.667605d0 !      Iodine -        Boron
      xfac(53, 5) =     3.161385d0 !      Iodine -        Boron
      alpb(53, 6) =     2.068710d0 !      Iodine -       Carbon
      xfac(53, 6) =     0.810156d0 !      Iodine -       Carbon
      alpb(53, 7) =     1.677518d0 !      Iodine -     Nitrogen
      xfac(53, 7) =     0.264903d0 !      Iodine -     Nitrogen
      alpb(53, 8) =     2.288919d0 !      Iodine -       Oxygen
      xfac(53, 8) =     0.866204d0 !      Iodine -       Oxygen
      alpb(53, 9) =     2.203580d0 !      Iodine -     Fluorine
      xfac(53, 9) =     0.392425d0 !      Iodine -     Fluorine
      alpb(53,10) =     2.414415d0 !      Iodine -         Neon
      xfac(53,10) =     1.503568d0 !      Iodine -         Neon
      alpb(53,11) =     1.403090d0 !      Iodine -       Sodium
      xfac(53,11) =     1.986112d0 !      Iodine -       Sodium
      alpb(53,12) =     2.045137d0 !      Iodine -    Magnesium
      xfac(53,12) =     3.276914d0 !      Iodine -    Magnesium
      alpb(53,13) =     1.816068d0 !      Iodine -     Aluminum
      xfac(53,13) =     2.929080d0 !      Iodine -     Aluminum
      alpb(53,14) =     1.559579d0 !      Iodine -      Silicon
      xfac(53,14) =     0.700299d0 !      Iodine -      Silicon
      alpb(53,15) =     2.131593d0 !      Iodine -   Phosphorus
      xfac(53,15) =     3.047207d0 !      Iodine -   Phosphorus
      alpb(53,16) =     1.855110d0 !      Iodine -       Sulfur
      xfac(53,16) =     0.709929d0 !      Iodine -       Sulfur
      alpb(53,17) =     1.574161d0 !      Iodine -     Chlorine
      xfac(53,17) =     0.310474d0 !      Iodine -     Chlorine
      alpb(53,18) =     1.576587d0 !      Iodine -        Argon
      xfac(53,18) =     0.305367d0 !      Iodine -        Argon
      alpb(53,19) =     1.539714d0 !      Iodine -    Potassium
      xfac(53,19) =     4.824353d0 !      Iodine -    Potassium
      alpb(53,20) =     2.196490d0 !      Iodine -      Calcium
      xfac(53,20) =     7.689921d0 !      Iodine -      Calcium
      alpb(53,21) =     1.814884d0 !      Iodine -     Scandium
      xfac(53,21) =     3.114282d0 !      Iodine -     Scandium
      alpb(53,22) =     1.933469d0 !      Iodine -     Titanium
      xfac(53,22) =     2.426747d0 !      Iodine -     Titanium
      alpb(53,23) =     2.683520d0 !      Iodine -     Vanadium
      xfac(53,23) =     6.198112d0 !      Iodine -     Vanadium
      alpb(53,24) =     2.634224d0 !      Iodine -     Chromium
      xfac(53,24) =     2.598590d0 !      Iodine -     Chromium
      alpb(53,25) =     2.266600d0 !      Iodine -    Manganese
      xfac(53,25) =     1.193410d0 !      Iodine -    Manganese
      alpb(53,26) =     1.912829d0 !      Iodine -         Iron
      xfac(53,26) =     0.532622d0 !      Iodine -         Iron
      alpb(53,27) =     3.235204d0 !      Iodine -       Cobalt
      xfac(53,27) =     1.105239d0 !      Iodine -       Cobalt
      alpb(53,28) =     1.085343d0 !      Iodine -       Nickel
      xfac(53,28) =     0.017459d0 !      Iodine -       Nickel
      alpb(53,29) =     0.834305d0 !      Iodine -       Copper
      xfac(53,29) =     0.006781d0 !      Iodine -       Copper
      alpb(53,30) =     1.394762d0 !      Iodine -         Zinc
      xfac(53,30) =     0.976607d0 !      Iodine -         Zinc
      alpb(53,31) =     1.671729d0 !      Iodine -      Gallium
      xfac(53,31) =     1.252168d0 !      Iodine -      Gallium
      alpb(53,32) =     1.817425d0 !      Iodine -    Germanium
      xfac(53,32) =     1.323267d0 !      Iodine -    Germanium
      alpb(53,33) =     1.245262d0 !      Iodine -      Arsenic
      xfac(53,33) =     0.310824d0 !      Iodine -      Arsenic
      alpb(53,35) =     1.579376d0 !      Iodine -      Bromine
      xfac(53,35) =     0.483054d0 !      Iodine -      Bromine
      alpb(53,36) =     1.238574d0 !      Iodine -      Krypton
      xfac(53,36) =     0.201136d0 !      Iodine -      Krypton
      alpb(53,37) =     1.432675d0 !      Iodine -     Rubidium
      xfac(53,37) =     4.092446d0 !      Iodine -     Rubidium
      alpb(53,38) =     1.262042d0 !      Iodine -    Strontium
      xfac(53,38) =     2.103941d0 !      Iodine -    Strontium
      alpb(53,39) =     1.279110d0 !      Iodine -      Yttrium
      xfac(53,39) =     1.021402d0 !      Iodine -      Yttrium
      alpb(53,40) =     1.995182d0 !      Iodine -    Zirconium
      xfac(53,40) =     4.513943d0 !      Iodine -    Zirconium
      alpb(53,41) =     1.967251d0 !      Iodine -      Niobium
      xfac(53,41) =     2.399298d0 !      Iodine -      Niobium
      alpb(53,42) =     0.948461d0 !      Iodine -   Molybdenum
      xfac(53,42) =     0.124695d0 !      Iodine -   Molybdenum
      alpb(53,43) =     1.292312d0 !      Iodine -   Technetium
      xfac(53,43) =     0.110594d0 !      Iodine -   Technetium
      alpb(53,44) =     3.953203d0 !      Iodine -    Ruthenium
      xfac(53,44) =     7.837710d0 !      Iodine -    Ruthenium
      alpb(53,45) =     3.708170d0 !      Iodine -      Rhodium
      xfac(53,45) =     2.357944d0 !      Iodine -      Rhodium
      alpb(53,46) =     5.144544d0 !      Iodine -    Palladium
      xfac(53,46) =     3.522017d0 !      Iodine -    Palladium
      alpb(53,47) =     2.593161d0 !      Iodine -       Silver
      xfac(53,47) =     0.048904d0 !      Iodine -       Silver
      alpb(53,48) =     0.996238d0 !      Iodine -      Cadmium
      xfac(53,48) =     0.396784d0 !      Iodine -      Cadmium
      alpb(53,49) =     2.351758d0 !      Iodine -       Indium
      xfac(53,49) =     5.947821d0 !      Iodine -       Indium
      alpb(53,50) =     1.855633d0 !      Iodine -          Tin
      xfac(53,50) =     1.783163d0 !      Iodine -          Tin
      alpb(53,51) =     1.155315d0 !      Iodine -     Antimony
      xfac(53,51) =     0.318190d0 !      Iodine -     Antimony
      alpb(53,52) =     1.493951d0 !      Iodine -    Tellurium
      xfac(53,52) =     1.101116d0 !      Iodine -    Tellurium
      alpb(53,53) =     1.519925d0 !      Iodine -       Iodine
      xfac(53,53) =     0.510542d0 !      Iodine -       Iodine
 !
      alpb(54, 1) =     1.356861d0 !       Xenon -     Hydrogen
      xfac(54, 1) =     0.701016d0 !       Xenon -     Hydrogen
      alpb(54, 2) =     2.497832d0 !       Xenon -       Helium
      xfac(54, 2) =     2.599471d0 !       Xenon -       Helium
      alpb(54, 3) =     2.466895d0 !       Xenon -      Lithium
      xfac(54, 3) =     4.582081d0 !       Xenon -      Lithium
      alpb(54, 4) =     6.000003d0 !       Xenon -    Beryllium
      xfac(54, 4) =     0.660525d0 !       Xenon -    Beryllium
      alpb(54, 5) =     5.051957d0 !       Xenon -        Boron
      xfac(54, 5) =     1.100612d0 !       Xenon -        Boron
      alpb(54, 6) =     1.704440d0 !       Xenon -       Carbon
      xfac(54, 6) =     0.826727d0 !       Xenon -       Carbon
      alpb(54, 7) =     1.932952d0 !       Xenon -     Nitrogen
      xfac(54, 7) =     0.925624d0 !       Xenon -     Nitrogen
      alpb(54, 8) =     0.839233d0 !       Xenon -       Oxygen
      xfac(54, 8) =     0.035356d0 !       Xenon -       Oxygen
      alpb(54, 9) =     1.128812d0 !       Xenon -     Fluorine
      xfac(54, 9) =     0.065011d0 !       Xenon -     Fluorine
      alpb(54,10) =     1.330202d0 !       Xenon -         Neon
      xfac(54,10) =     0.293862d0 !       Xenon -         Neon
      alpb(54,11) =     2.103003d0 !       Xenon -       Sodium
      xfac(54,11) =     8.368204d0 !       Xenon -       Sodium
      alpb(54,12) =     2.698414d0 !       Xenon -    Magnesium
      xfac(54,12) =     9.723572d0 !       Xenon -    Magnesium
      alpb(54,13) =     2.412039d0 !       Xenon -     Aluminum
      xfac(54,13) =     7.404465d0 !       Xenon -     Aluminum
      alpb(54,14) =     3.087060d0 !       Xenon -      Silicon
      xfac(54,14) =    16.092000d0 !       Xenon -      Silicon
      alpb(54,17) =     1.546396d0 !       Xenon -     Chlorine
      xfac(54,17) =     0.463758d0 !       Xenon -     Chlorine
      alpb(54,18) =     0.591520d0 !       Xenon -        Argon
      xfac(54,18) =     0.049266d0 !       Xenon -        Argon
      alpb(54,19) =     1.171250d0 !       Xenon -    Potassium
      xfac(54,19) =     1.224889d0 !       Xenon -    Potassium
      alpb(54,20) =     1.510653d0 !       Xenon -      Calcium
      xfac(54,20) =     1.717121d0 !       Xenon -      Calcium
      alpb(54,35) =     1.439618d0 !       Xenon -      Bromine
      xfac(54,35) =     0.475116d0 !       Xenon -      Bromine
      alpb(54,36) =     0.551561d0 !       Xenon -      Krypton
      xfac(54,36) =     0.049793d0 !       Xenon -      Krypton
      alpb(54,37) =     1.087823d0 !       Xenon -     Rubidium
      xfac(54,37) =     0.974965d0 !       Xenon -     Rubidium
      alpb(54,53) =     0.799155d0 !       Xenon -       Iodine
      xfac(54,53) =     0.112090d0 !       Xenon -       Iodine
      alpb(54,54) =     1.244762d0 !       Xenon -        Xenon
      xfac(54,54) =     0.344474d0 !       Xenon -        Xenon
 !
      alpb(55, 1) =     0.264882d0 !      Cesium -     Hydrogen
      xfac(55, 1) =     0.096901d0 !      Cesium -     Hydrogen
      alpb(55, 5) =     1.487110d0 !      Cesium -        Boron
      xfac(55, 5) =    10.392610d0 !      Cesium -        Boron
      alpb(55, 6) =     2.147104d0 !      Cesium -       Carbon
      xfac(55, 6) =    24.514623d0 !      Cesium -       Carbon
      alpb(55, 7) =     2.446532d0 !      Cesium -     Nitrogen
      xfac(55, 7) =    29.711077d0 !      Cesium -     Nitrogen
      alpb(55, 8) =     2.085139d0 !      Cesium -       Oxygen
      xfac(55, 8) =     8.176843d0 !      Cesium -       Oxygen
      alpb(55, 9) =     2.834100d0 !      Cesium -     Fluorine
      xfac(55, 9) =    22.233416d0 !      Cesium -     Fluorine
      alpb(55,15) =     2.924953d0 !      Cesium -   Phosphorus
      xfac(55,15) =     0.506512d0 !      Cesium -   Phosphorus
      alpb(55,16) =     0.289412d0 !      Cesium -       Sulfur
      xfac(55,16) =     0.091743d0 !      Cesium -       Sulfur
      alpb(55,17) =     1.673663d0 !      Cesium -     Chlorine
      xfac(55,17) =     4.531965d0 !      Cesium -     Chlorine
      alpb(55,35) =     1.167189d0 !      Cesium -      Bromine
      xfac(55,35) =     1.658427d0 !      Cesium -      Bromine
      alpb(55,53) =     0.919562d0 !      Cesium -       Iodine
      xfac(55,53) =     1.072178d0 !      Cesium -       Iodine
      alpb(55,55) =     1.170843d0 !      Cesium -       Cesium
      xfac(55,55) =    25.320055d0 !      Cesium -       Cesium
 !
      alpb(56, 1) =     6.000135d0 !      Barium -     Hydrogen
      xfac(56, 1) =     2.040004d0 !      Barium -     Hydrogen
      alpb(56, 6) =     0.770626d0 !      Barium -       Carbon
      xfac(56, 6) =     0.119793d0 !      Barium -       Carbon
      alpb(56, 7) =     1.148233d0 !      Barium -     Nitrogen
      xfac(56, 7) =     0.207934d0 !      Barium -     Nitrogen
      alpb(56, 8) =     1.283018d0 !      Barium -       Oxygen
      xfac(56, 8) =     0.348945d0 !      Barium -       Oxygen
      alpb(56, 9) =     3.000618d0 !      Barium -     Fluorine
      xfac(56, 9) =     5.575255d0 !      Barium -     Fluorine
      alpb(56,13) =     2.105924d0 !      Barium -     Aluminum
      xfac(56,13) =     9.539099d0 !      Barium -     Aluminum
      alpb(56,14) =     1.240420d0 !      Barium -      Silicon
      xfac(56,14) =     1.212660d0 !      Barium -      Silicon
      alpb(56,16) =     0.705188d0 !      Barium -       Sulfur
      xfac(56,16) =     0.215386d0 !      Barium -       Sulfur
      alpb(56,17) =     1.071044d0 !      Barium -     Chlorine
      xfac(56,17) =     0.160177d0 !      Barium -     Chlorine
      alpb(56,22) =     2.176040d0 !      Barium -     Titanium
      xfac(56,22) =     9.493530d0 !      Barium -     Titanium
      alpb(56,35) =     1.190346d0 !      Barium -      Bromine
      xfac(56,35) =     0.828794d0 !      Barium -      Bromine
      alpb(56,53) =     0.982528d0 !      Barium -       Iodine
      xfac(56,53) =     0.835597d0 !      Barium -       Iodine
      alpb(56,56) =     0.339269d0 !      Barium -       Barium
      xfac(56,56) =     0.356186d0 !      Barium -       Barium
 !
      alpb(57, 1) =     0.833667d0 !   Lanthanum -     Hydrogen
      xfac(57, 1) =     0.623501d0 !   Lanthanum -     Hydrogen
      alpb(57, 6) =     0.604869d0 !   Lanthanum -       Carbon
      xfac(57, 6) =     0.108649d0 !   Lanthanum -       Carbon
      alpb(57, 7) =     0.758881d0 !   Lanthanum -     Nitrogen
      xfac(57, 7) =     0.104778d0 !   Lanthanum -     Nitrogen
      alpb(57, 8) =     1.318333d0 !   Lanthanum -       Oxygen
      xfac(57, 8) =     0.557957d0 !   Lanthanum -       Oxygen
      alpb(57, 9) =     2.379335d0 !   Lanthanum -     Fluorine
      xfac(57, 9) =     2.401903d0 !   Lanthanum -     Fluorine
      alpb(57,13) =     1.003510d0 !   Lanthanum -     Aluminum
      xfac(57,13) =     0.500540d0 !   Lanthanum -     Aluminum
      alpb(57,14) =     2.016820d0 !   Lanthanum -      Silicon
      xfac(57,14) =     3.219030d0 !   Lanthanum -      Silicon
      alpb(57,15) =     0.954450d0 !   Lanthanum -   Phosphorus
      xfac(57,15) =     0.541660d0 !   Lanthanum -   Phosphorus
      alpb(57,16) =     1.834129d0 !   Lanthanum -       Sulfur
      xfac(57,16) =     2.682412d0 !   Lanthanum -       Sulfur
      alpb(57,17) =     0.993753d0 !   Lanthanum -     Chlorine
      xfac(57,17) =     0.230203d0 !   Lanthanum -     Chlorine
      alpb(57,35) =     0.758184d0 !   Lanthanum -      Bromine
      xfac(57,35) =     0.238582d0 !   Lanthanum -      Bromine
      alpb(57,53) =     0.592666d0 !   Lanthanum -       Iodine
      xfac(57,53) =     0.226883d0 !   Lanthanum -       Iodine
      alpb(57,57) =     4.248067d0 !   Lanthanum -    Lanthanum
      xfac(57,57) =     5.175162d0 !   Lanthanum -    Lanthanum
 !
      alpb(64, 1) =     0.390870d0 !  Gadolinium -     Hydrogen
      xfac(64, 1) =     0.135810d0 !  Gadolinium -     Hydrogen
      alpb(64, 6) =     0.446870d0 !  Gadolinium -       Carbon
      xfac(64, 6) =     0.053040d0 !  Gadolinium -       Carbon
      alpb(64, 7) =     1.159410d0 !  Gadolinium -     Nitrogen
      xfac(64, 7) =     0.205050d0 !  Gadolinium -     Nitrogen
      alpb(64, 8) =     0.862040d0 !  Gadolinium -       Oxygen
      xfac(64, 8) =     0.175800d0 !  Gadolinium -       Oxygen
      alpb(64, 9) =     1.497980d0 !  Gadolinium -     Fluorine
      xfac(64, 9) =     0.334630d0 !  Gadolinium -     Fluorine
      alpb(64,13) =     1.003510d0 !  Gadolinium -     Aluminum
      xfac(64,13) =     0.500540d0 !  Gadolinium -     Aluminum
      alpb(64,14) =     2.016820d0 !  Gadolinium -      Silicon
      xfac(64,14) =     3.219030d0 !  Gadolinium -      Silicon
      alpb(64,15) =     0.954450d0 !  Gadolinium -   Phosphorus
      xfac(64,15) =     0.541660d0 !  Gadolinium -   Phosphorus
      alpb(64,16) =     2.003930d0 !  Gadolinium -       Sulfur
      xfac(64,16) =     2.655400d0 !  Gadolinium -       Sulfur
      alpb(64,17) =     0.806810d0 !  Gadolinium -     Chlorine
      xfac(64,17) =     0.089970d0 !  Gadolinium -     Chlorine
      alpb(64,35) =     0.715810d0 !  Gadolinium -      Bromine
      xfac(64,35) =     0.240740d0 !  Gadolinium -      Bromine
      alpb(64,53) =     0.585360d0 !  Gadolinium -       Iodine
      xfac(64,53) =     0.278240d0 !  Gadolinium -       Iodine
      alpb(64,64) =     3.348180d0 !  Gadolinium -   Gadolinium
      xfac(64,64) =     2.670400d0 !  Gadolinium -   Gadolinium
 !
      alpb(71, 1) =     1.415790d0 !    Lutetium -     Hydrogen
      xfac(71, 1) =     0.787920d0 !    Lutetium -     Hydrogen
      alpb(71, 6) =     2.312813d0 !    Lutetium -       Carbon
      xfac(71, 6) =     4.453825d0 !    Lutetium -       Carbon
      alpb(71, 7) =     2.141302d0 !    Lutetium -     Nitrogen
      xfac(71, 7) =     2.860828d0 !    Lutetium -     Nitrogen
      alpb(71, 8) =     2.192486d0 !    Lutetium -       Oxygen
      xfac(71, 8) =     2.917076d0 !    Lutetium -       Oxygen
      alpb(71,15) =     5.618820d0 !    Lutetium -   Phosphorus
      xfac(71,15) =     0.500000d0 !    Lutetium -   Phosphorus
      alpb(71,17) =     2.753636d0 !    Lutetium -     Chlorine
      xfac(71,17) =    12.757099d0 !    Lutetium -     Chlorine
      alpb(71,35) =     2.322618d0 !    Lutetium -      Bromine
      xfac(71,35) =     8.648274d0 !    Lutetium -      Bromine
      alpb(71,53) =     2.248348d0 !    Lutetium -       Iodine
      xfac(71,53) =    10.082315d0 !    Lutetium -       Iodine
 !
      alpb(72, 1) =     1.423788d0 !     Hafnium -     Hydrogen
      xfac(72, 1) =     3.427312d0 !     Hafnium -     Hydrogen
      alpb(72, 5) =     1.633500d0 !     Hafnium -        Boron
      xfac(72, 5) =     0.659270d0 !     Hafnium -        Boron
      alpb(72, 6) =     1.002194d0 !     Hafnium -       Carbon
      xfac(72, 6) =     0.378579d0 !     Hafnium -       Carbon
      alpb(72, 7) =     1.332410d0 !     Hafnium -     Nitrogen
      xfac(72, 7) =     0.655795d0 !     Hafnium -     Nitrogen
      alpb(72, 8) =     1.633289d0 !     Hafnium -       Oxygen
      xfac(72, 8) =     1.034718d0 !     Hafnium -       Oxygen
      alpb(72, 9) =     2.290803d0 !     Hafnium -     Fluorine
      xfac(72, 9) =     1.679335d0 !     Hafnium -     Fluorine
      alpb(72,12) =     1.911350d0 !     Hafnium -    Magnesium
      xfac(72,12) =     4.330250d0 !     Hafnium -    Magnesium
      alpb(72,13) =     0.949150d0 !     Hafnium -     Aluminum
      xfac(72,13) =     0.622520d0 !     Hafnium -     Aluminum
      alpb(72,14) =     2.189300d0 !     Hafnium -      Silicon
      xfac(72,14) =     3.382300d0 !     Hafnium -      Silicon
      alpb(72,15) =     1.231220d0 !     Hafnium -   Phosphorus
      xfac(72,15) =     0.505530d0 !     Hafnium -   Phosphorus
      alpb(72,16) =     2.327110d0 !     Hafnium -       Sulfur
      xfac(72,16) =     1.666760d0 !     Hafnium -       Sulfur
      alpb(72,17) =     1.297117d0 !     Hafnium -     Chlorine
      xfac(72,17) =     0.706421d0 !     Hafnium -     Chlorine
      alpb(72,20) =     2.054500d0 !     Hafnium -      Calcium
      xfac(72,20) =     4.319510d0 !     Hafnium -      Calcium
      alpb(72,33) =     1.799500d0 !     Hafnium -      Arsenic
      xfac(72,33) =     1.280820d0 !     Hafnium -      Arsenic
      alpb(72,35) =     1.090759d0 !     Hafnium -      Bromine
      xfac(72,35) =     0.692456d0 !     Hafnium -      Bromine
      alpb(72,53) =     1.014096d0 !     Hafnium -       Iodine
      xfac(72,53) =     0.820948d0 !     Hafnium -       Iodine
      alpb(72,56) =     2.264830d0 !     Hafnium -       Barium
      xfac(72,56) =     9.022520d0 !     Hafnium -       Barium
      alpb(72,72) =     0.544144d0 !     Hafnium -      Hafnium
      xfac(72,72) =     1.058911d0 !     Hafnium -      Hafnium
 !
      alpb(73, 1) =     2.288014d0 !    Tantalum -     Hydrogen
      xfac(73, 1) =     2.827669d0 !    Tantalum -     Hydrogen
      alpb(73, 6) =     1.838949d0 !    Tantalum -       Carbon
      xfac(73, 6) =     0.847439d0 !    Tantalum -       Carbon
      alpb(73, 7) =     2.053679d0 !    Tantalum -     Nitrogen
      xfac(73, 7) =     1.015461d0 !    Tantalum -     Nitrogen
      alpb(73, 8) =     2.412629d0 !    Tantalum -       Oxygen
      xfac(73, 8) =     1.751083d0 !    Tantalum -       Oxygen
      alpb(73, 9) =     3.107390d0 !    Tantalum -     Fluorine
      xfac(73, 9) =     3.146520d0 !    Tantalum -     Fluorine
      alpb(73,11) =     2.551120d0 !    Tantalum -       Sodium
      xfac(73,11) =     8.276130d0 !    Tantalum -       Sodium
      alpb(73,15) =     2.513800d0 !    Tantalum -   Phosphorus
      xfac(73,15) =     6.261880d0 !    Tantalum -   Phosphorus
      alpb(73,16) =     2.246723d0 !    Tantalum -       Sulfur
      xfac(73,16) =     2.975980d0 !    Tantalum -       Sulfur
      alpb(73,17) =     1.608805d0 !    Tantalum -     Chlorine
      xfac(73,17) =     0.516413d0 !    Tantalum -     Chlorine
      alpb(73,19) =     4.521470d0 !    Tantalum -    Potassium
      xfac(73,19) =     2.026700d0 !    Tantalum -    Potassium
      alpb(73,35) =     1.640376d0 !    Tantalum -      Bromine
      xfac(73,35) =     0.791445d0 !    Tantalum -      Bromine
      alpb(73,53) =     2.401053d0 !    Tantalum -       Iodine
      xfac(73,53) =     6.551551d0 !    Tantalum -       Iodine
      alpb(73,73) =     2.082863d0 !    Tantalum -     Tantalum
      xfac(73,73) =    10.987053d0 !    Tantalum -     Tantalum
 !
      alpb(74, 1) =     2.130880d0 !    Tungsten -     Hydrogen
      xfac(74, 1) =     1.832270d0 !    Tungsten -     Hydrogen
      alpb(74, 6) =     2.097480d0 !    Tungsten -       Carbon
      xfac(74, 6) =     1.160770d0 !    Tungsten -       Carbon
      alpb(74, 7) =     1.596040d0 !    Tungsten -     Nitrogen
      xfac(74, 7) =     0.478350d0 !    Tungsten -     Nitrogen
      alpb(74, 8) =     1.359020d0 !    Tungsten -       Oxygen
      xfac(74, 8) =     0.349010d0 !    Tungsten -       Oxygen
      alpb(74, 9) =     1.446050d0 !    Tungsten -     Fluorine
      xfac(74, 9) =     0.213890d0 !    Tungsten -     Fluorine
      alpb(74,11) =     2.551030d0 !    Tungsten -       Sodium
      xfac(74,11) =     8.276040d0 !    Tungsten -       Sodium
      alpb(74,15) =     2.338060d0 !    Tungsten -   Phosphorus
      xfac(74,15) =     5.953860d0 !    Tungsten -   Phosphorus
      alpb(74,16) =     1.542570d0 !    Tungsten -       Sulfur
      xfac(74,16) =     0.488630d0 !    Tungsten -       Sulfur
      alpb(74,17) =     1.310690d0 !    Tungsten -     Chlorine
      xfac(74,17) =     0.278000d0 !    Tungsten -     Chlorine
      alpb(74,19) =     4.521380d0 !    Tungsten -    Potassium
      xfac(74,19) =     2.026610d0 !    Tungsten -    Potassium
      alpb(74,35) =     1.293260d0 !    Tungsten -      Bromine
      xfac(74,35) =     0.372390d0 !    Tungsten -      Bromine
      alpb(74,53) =     1.573570d0 !    Tungsten -       Iodine
      xfac(74,53) =     1.077370d0 !    Tungsten -       Iodine
      alpb(74,74) =     2.940870d0 !    Tungsten -     Tungsten
      xfac(74,74) =     7.471390d0 !    Tungsten -     Tungsten
 !
      alpb(75, 1) =     1.634500d0 !     Rhenium -     Hydrogen
      xfac(75, 1) =     0.345894d0 !     Rhenium -     Hydrogen
      alpb(75, 6) =     2.306285d0 !     Rhenium -       Carbon
      xfac(75, 6) =     0.690687d0 !     Rhenium -       Carbon
      alpb(75, 7) =     1.918332d0 !     Rhenium -     Nitrogen
      xfac(75, 7) =     0.445213d0 !     Rhenium -     Nitrogen
      alpb(75, 8) =     1.967747d0 !     Rhenium -       Oxygen
      xfac(75, 8) =     0.635960d0 !     Rhenium -       Oxygen
      alpb(75, 9) =     2.154219d0 !     Rhenium -     Fluorine
      xfac(75, 9) =     0.535966d0 !     Rhenium -     Fluorine
      alpb(75,14) =     2.775930d0 !     Rhenium -      Silicon
      xfac(75,14) =     0.849450d0 !     Rhenium -      Silicon
      alpb(75,15) =     1.804168d0 !     Rhenium -   Phosphorus
      xfac(75,15) =     0.966942d0 !     Rhenium -   Phosphorus
      alpb(75,16) =     1.083919d0 !     Rhenium -       Sulfur
      xfac(75,16) =     0.068874d0 !     Rhenium -       Sulfur
      alpb(75,17) =     1.433875d0 !     Rhenium -     Chlorine
      xfac(75,17) =     0.146319d0 !     Rhenium -     Chlorine
      alpb(75,32) =     2.852340d0 !     Rhenium -    Germanium
      xfac(75,32) =     2.151580d0 !     Rhenium -    Germanium
      alpb(75,34) =     2.523170d0 !     Rhenium -     Selenium
      xfac(75,34) =     2.202140d0 !     Rhenium -     Selenium
      alpb(75,35) =     1.603060d0 !     Rhenium -      Bromine
      xfac(75,35) =     0.287528d0 !     Rhenium -      Bromine
      alpb(75,51) =     2.204360d0 !     Rhenium -     Antimony
      xfac(75,51) =     2.275780d0 !     Rhenium -     Antimony
      alpb(75,53) =     2.610119d0 !     Rhenium -       Iodine
      xfac(75,53) =     3.559286d0 !     Rhenium -       Iodine
      alpb(75,75) =     6.000258d0 !     Rhenium -      Rhenium
      xfac(75,75) =     4.488852d0 !     Rhenium -      Rhenium
 !
      alpb(76, 1) =     3.404180d0 !      Osmium -     Hydrogen
      xfac(76, 1) =     4.393870d0 !      Osmium -     Hydrogen
      alpb(76, 6) =     2.336500d0 !      Osmium -       Carbon
      xfac(76, 6) =     0.498410d0 !      Osmium -       Carbon
      alpb(76, 7) =     1.143090d0 !      Osmium -     Nitrogen
      xfac(76, 7) =     0.080870d0 !      Osmium -     Nitrogen
      alpb(76, 8) =     1.350360d0 !      Osmium -       Oxygen
      xfac(76, 8) =     0.184300d0 !      Osmium -       Oxygen
      alpb(76, 9) =     1.507620d0 !      Osmium -     Fluorine
      xfac(76, 9) =     0.140050d0 !      Osmium -     Fluorine
      alpb(76,11) =     2.550740d0 !      Osmium -       Sodium
      xfac(76,11) =     8.275750d0 !      Osmium -       Sodium
      alpb(76,15) =     2.836090d0 !      Osmium -   Phosphorus
      xfac(76,15) =     6.058300d0 !      Osmium -   Phosphorus
      alpb(76,16) =     2.809500d0 !      Osmium -       Sulfur
      xfac(76,16) =     4.186050d0 !      Osmium -       Sulfur
      alpb(76,17) =     1.833070d0 !      Osmium -     Chlorine
      xfac(76,17) =     0.327920d0 !      Osmium -     Chlorine
      alpb(76,19) =     4.521090d0 !      Osmium -    Potassium
      xfac(76,19) =     2.026320d0 !      Osmium -    Potassium
      alpb(76,35) =     1.766880d0 !      Osmium -      Bromine
      xfac(76,35) =     0.382430d0 !      Osmium -      Bromine
      alpb(76,53) =     2.203760d0 !      Osmium -       Iodine
      xfac(76,53) =     2.199190d0 !      Osmium -       Iodine
      alpb(76,76) =     2.021630d0 !      Osmium -       Osmium
      xfac(76,76) =     0.830440d0 !      Osmium -       Osmium
 !
      alpb(77, 1) =     1.033900d0 !     Iridium -     Hydrogen
      xfac(77, 1) =     0.058047d0 !     Iridium -     Hydrogen
      alpb(77, 6) =     1.690295d0 !     Iridium -       Carbon
      xfac(77, 6) =     0.115047d0 !     Iridium -       Carbon
      alpb(77, 7) =     3.934508d0 !     Iridium -     Nitrogen
      xfac(77, 7) =     8.518640d0 !     Iridium -     Nitrogen
      alpb(77, 8) =     3.748272d0 !     Iridium -       Oxygen
      xfac(77, 8) =     9.625402d0 !     Iridium -       Oxygen
      alpb(77, 9) =     2.982799d0 !     Iridium -     Fluorine
      xfac(77, 9) =     1.499639d0 !     Iridium -     Fluorine
      alpb(77,11) =     2.550820d0 !     Iridium -       Sodium
      xfac(77,11) =     8.275830d0 !     Iridium -       Sodium
      alpb(77,15) =     2.714060d0 !     Iridium -   Phosphorus
      xfac(77,15) =     6.284670d0 !     Iridium -   Phosphorus
      alpb(77,16) =     3.204834d0 !     Iridium -       Sulfur
      xfac(77,16) =     4.135732d0 !     Iridium -       Sulfur
      alpb(77,17) =     2.009770d0 !     Iridium -     Chlorine
      xfac(77,17) =     0.258916d0 !     Iridium -     Chlorine
      alpb(77,19) =     4.521170d0 !     Iridium -    Potassium
      xfac(77,19) =     2.026400d0 !     Iridium -    Potassium
      alpb(77,35) =     2.038142d0 !     Iridium -      Bromine
      xfac(77,35) =     0.171879d0 !     Iridium -      Bromine
      alpb(77,53) =     3.410914d0 !     Iridium -       Iodine
      xfac(77,53) =     1.497148d0 !     Iridium -       Iodine
      alpb(77,77) =     5.771663d0 !     Iridium -      Iridium
      xfac(77,77) =    11.175193d0 !     Iridium -      Iridium
 !
      alpb(78, 1) =     4.001198d0 !    Platinum -     Hydrogen
      xfac(78, 1) =     8.924015d0 !    Platinum -     Hydrogen
      alpb(78, 6) =     3.306722d0 !    Platinum -       Carbon
      xfac(78, 6) =     3.493403d0 !    Platinum -       Carbon
      alpb(78, 7) =     2.307923d0 !    Platinum -     Nitrogen
      xfac(78, 7) =     0.540730d0 !    Platinum -     Nitrogen
      alpb(78, 8) =     2.110563d0 !    Platinum -       Oxygen
      xfac(78, 8) =     0.487756d0 !    Platinum -       Oxygen
      alpb(78, 9) =     3.714441d0 !    Platinum -     Fluorine
      xfac(78, 9) =     5.617014d0 !    Platinum -     Fluorine
      alpb(78,13) =     1.572360d0 !    Platinum -     Aluminum
      xfac(78,13) =     1.056930d0 !    Platinum -     Aluminum
      alpb(78,14) =     0.999990d0 !    Platinum -      Silicon
      xfac(78,14) =     0.099990d0 !    Platinum -      Silicon
      alpb(78,15) =     1.403239d0 !    Platinum -   Phosphorus
      xfac(78,15) =     0.233712d0 !    Platinum -   Phosphorus
      alpb(78,16) =     2.791500d0 !    Platinum -       Sulfur
      xfac(78,16) =     2.224263d0 !    Platinum -       Sulfur
      alpb(78,17) =     2.108526d0 !    Platinum -     Chlorine
      xfac(78,17) =     0.341001d0 !    Platinum -     Chlorine
      alpb(78,35) =     2.185307d0 !    Platinum -      Bromine
      xfac(78,35) =     0.520361d0 !    Platinum -      Bromine
      alpb(78,53) =     3.077338d0 !    Platinum -       Iodine
      xfac(78,53) =     4.601248d0 !    Platinum -       Iodine
      alpb(78,78) =     3.404276d0 !    Platinum -     Platinum
      xfac(78,78) =     9.010252d0 !    Platinum -     Platinum
 !
      alpb(79, 1) =     3.369041d0 !        Gold -     Hydrogen
      xfac(79, 1) =     2.605283d0 !        Gold -     Hydrogen
      alpb(79, 6) =     4.580016d0 !        Gold -       Carbon
      xfac(79, 6) =    21.485634d0 !        Gold -       Carbon
      alpb(79, 7) =     2.138095d0 !        Gold -     Nitrogen
      xfac(79, 7) =     0.222059d0 !        Gold -     Nitrogen
      alpb(79, 8) =     1.548763d0 !        Gold -       Oxygen
      xfac(79, 8) =     0.077192d0 !        Gold -       Oxygen
      alpb(79, 9) =     4.453145d0 !        Gold -     Fluorine
      xfac(79, 9) =     9.594384d0 !        Gold -     Fluorine
      alpb(79,13) =     1.572570d0 !        Gold -     Aluminum
      xfac(79,13) =     1.057140d0 !        Gold -     Aluminum
      alpb(79,15) =     1.618713d0 !        Gold -   Phosphorus
      xfac(79,15) =     0.067001d0 !        Gold -   Phosphorus
      alpb(79,16) =     4.306238d0 !        Gold -       Sulfur
      xfac(79,16) =    21.619145d0 !        Gold -       Sulfur
      alpb(79,17) =     3.539414d0 !        Gold -     Chlorine
      xfac(79,17) =     2.257702d0 !        Gold -     Chlorine
      alpb(79,35) =     0.581911d0 !        Gold -      Bromine
      xfac(79,35) =     0.004237d0 !        Gold -      Bromine
      alpb(79,53) =     0.577916d0 !        Gold -       Iodine
      xfac(79,53) =     0.008816d0 !        Gold -       Iodine
      alpb(79,79) =     0.903162d0 !        Gold -         Gold
      xfac(79,79) =     0.013091d0 !        Gold -         Gold
 !
      alpb(80, 1) =     1.136587d0 !     Mercury -     Hydrogen
      xfac(80, 1) =     0.799399d0 !     Mercury -     Hydrogen
      alpb(80, 6) =     0.795816d0 !     Mercury -       Carbon
      xfac(80, 6) =     0.147128d0 !     Mercury -       Carbon
      alpb(80, 7) =     0.332152d0 !     Mercury -     Nitrogen
      xfac(80, 7) =     0.050240d0 !     Mercury -     Nitrogen
      alpb(80, 8) =     1.052145d0 !     Mercury -       Oxygen
      xfac(80, 8) =     0.240720d0 !     Mercury -       Oxygen
      alpb(80, 9) =     1.240572d0 !     Mercury -     Fluorine
      xfac(80, 9) =     0.113827d0 !     Mercury -     Fluorine
      alpb(80,14) =     2.770860d0 !     Mercury -      Silicon
      xfac(80,14) =     3.680740d0 !     Mercury -      Silicon
      alpb(80,15) =     0.608604d0 !     Mercury -   Phosphorus
      xfac(80,15) =     0.214951d0 !     Mercury -   Phosphorus
      alpb(80,16) =     1.041682d0 !     Mercury -       Sulfur
      xfac(80,16) =     0.347383d0 !     Mercury -       Sulfur
      alpb(80,17) =     0.430731d0 !     Mercury -     Chlorine
      xfac(80,17) =     0.053660d0 !     Mercury -     Chlorine
      alpb(80,22) =     3.414630d0 !     Mercury -     Titanium
      xfac(80,22) =     2.957200d0 !     Mercury -     Titanium
      alpb(80,35) =     0.638717d0 !     Mercury -      Bromine
      xfac(80,35) =     0.172363d0 !     Mercury -      Bromine
      alpb(80,52) =     0.291500d0 !     Mercury -    Tellurium
      xfac(80,52) =     0.212732d0 !     Mercury -    Tellurium
      alpb(80,53) =     0.758162d0 !     Mercury -       Iodine
      xfac(80,53) =     0.342058d0 !     Mercury -       Iodine
      alpb(80,80) =     0.474413d0 !     Mercury -      Mercury
      xfac(80,80) =     0.423276d0 !     Mercury -      Mercury
 !
      alpb(81, 1) =     0.673658d0 !    Thallium -     Hydrogen
      xfac(81, 1) =     0.138205d0 !    Thallium -     Hydrogen
      alpb(81, 5) =     1.528347d0 !    Thallium -        Boron
      xfac(81, 5) =    10.504338d0 !    Thallium -        Boron
      alpb(81, 6) =     1.390345d0 !    Thallium -       Carbon
      xfac(81, 6) =     0.582895d0 !    Thallium -       Carbon
      alpb(81, 7) =     0.982335d0 !    Thallium -     Nitrogen
      xfac(81, 7) =     0.158812d0 !    Thallium -     Nitrogen
      alpb(81, 8) =     1.550068d0 !    Thallium -       Oxygen
      xfac(81, 8) =     0.636906d0 !    Thallium -       Oxygen
      alpb(81, 9) =     1.469516d0 !    Thallium -     Fluorine
      xfac(81, 9) =     0.226166d0 !    Thallium -     Fluorine
      alpb(81,16) =     0.994851d0 !    Thallium -       Sulfur
      xfac(81,16) =     0.303426d0 !    Thallium -       Sulfur
      alpb(81,17) =     0.846193d0 !    Thallium -     Chlorine
      xfac(81,17) =     0.162037d0 !    Thallium -     Chlorine
      alpb(81,35) =     0.874419d0 !    Thallium -      Bromine
      xfac(81,35) =     0.296836d0 !    Thallium -      Bromine
      alpb(81,53) =     0.902012d0 !    Thallium -       Iodine
      xfac(81,53) =     0.430033d0 !    Thallium -       Iodine
      alpb(81,81) =     1.191684d0 !    Thallium -     Thallium
      xfac(81,81) =     9.535127d0 !    Thallium -     Thallium
 !
      alpb(82, 1) =     1.522676d0 !        Lead -     Hydrogen
      xfac(82, 1) =     0.840096d0 !        Lead -     Hydrogen
      alpb(82, 3) =     1.001810d0 !        Lead -      Lithium
      xfac(82, 3) =     1.285064d0 !        Lead -      Lithium
      alpb(82, 5) =     0.911197d0 !        Lead -        Boron
      xfac(82, 5) =     1.138157d0 !        Lead -        Boron
      alpb(82, 6) =     1.525593d0 !        Lead -       Carbon
      xfac(82, 6) =     0.404656d0 !        Lead -       Carbon
      alpb(82, 7) =     1.317394d0 !        Lead -     Nitrogen
      xfac(82, 7) =     0.335787d0 !        Lead -     Nitrogen
      alpb(82, 8) =     1.763210d0 !        Lead -       Oxygen
      xfac(82, 8) =     0.782506d0 !        Lead -       Oxygen
      alpb(82, 9) =     3.288902d0 !        Lead -     Fluorine
      xfac(82, 9) =     8.368562d0 !        Lead -     Fluorine
      alpb(82,15) =     4.516800d0 !        Lead -   Phosphorus
      xfac(82,15) =     5.033200d0 !        Lead -   Phosphorus
      alpb(82,16) =     1.027519d0 !        Lead -       Sulfur
      xfac(82,16) =     0.175150d0 !        Lead -       Sulfur
      alpb(82,17) =     1.094123d0 !        Lead -     Chlorine
      xfac(82,17) =     0.164814d0 !        Lead -     Chlorine
      alpb(82,23) =     1.500000d0 !        Lead -     Vanadium
      xfac(82,23) =     1.000000d0 !        Lead -     Vanadium
      alpb(82,24) =     1.860760d0 !        Lead -     Chromium
      xfac(82,24) =     1.029110d0 !        Lead -     Chromium
      alpb(82,30) =     1.500000d0 !        Lead -         Zinc
      xfac(82,30) =     1.000000d0 !        Lead -         Zinc
      alpb(82,34) =     2.000000d0 !        Lead -     Selenium
      xfac(82,34) =     0.111195d0 !        Lead -     Selenium
      alpb(82,35) =     0.865550d0 !        Lead -      Bromine
      xfac(82,35) =     0.148229d0 !        Lead -      Bromine
      alpb(82,41) =     1.500000d0 !        Lead -      Niobium
      xfac(82,41) =     1.000000d0 !        Lead -      Niobium
      alpb(82,42) =     2.000000d0 !        Lead -   Molybdenum
      xfac(82,42) =     5.000000d0 !        Lead -   Molybdenum
      alpb(82,52) =     1.002559d0 !        Lead -    Tellurium
      xfac(82,52) =     0.809042d0 !        Lead -    Tellurium
      alpb(82,53) =     0.983474d0 !        Lead -       Iodine
      xfac(82,53) =     0.267426d0 !        Lead -       Iodine
      alpb(82,82) =     1.881764d0 !        Lead -         Lead
      xfac(82,82) =     2.362343d0 !        Lead -         Lead
 !
      alpb(83, 1) =     1.679905d0 !     Bismuth -     Hydrogen
      xfac(83, 1) =     1.397462d0 !     Bismuth -     Hydrogen
      alpb(83, 3) =     0.340140d0 !     Bismuth -      Lithium
      xfac(83, 3) =     0.695320d0 !     Bismuth -      Lithium
      alpb(83, 6) =     1.534025d0 !     Bismuth -       Carbon
      xfac(83, 6) =     0.576179d0 !     Bismuth -       Carbon
      alpb(83, 7) =     1.143876d0 !     Bismuth -     Nitrogen
      xfac(83, 7) =     0.152738d0 !     Bismuth -     Nitrogen
      alpb(83, 8) =     1.553297d0 !     Bismuth -       Oxygen
      xfac(83, 8) =     0.333042d0 !     Bismuth -       Oxygen
      alpb(83, 9) =     2.355400d0 !     Bismuth -     Fluorine
      xfac(83, 9) =     1.035324d0 !     Bismuth -     Fluorine
      alpb(83,16) =     1.466879d0 !     Bismuth -       Sulfur
      xfac(83,16) =     0.620997d0 !     Bismuth -       Sulfur
      alpb(83,17) =     1.272975d0 !     Bismuth -     Chlorine
      xfac(83,17) =     0.326871d0 !     Bismuth -     Chlorine
      alpb(83,34) =     1.344746d0 !     Bismuth -     Selenium
      xfac(83,34) =     0.651208d0 !     Bismuth -     Selenium
      alpb(83,35) =     1.146233d0 !     Bismuth -      Bromine
      xfac(83,35) =     0.381170d0 !     Bismuth -      Bromine
      alpb(83,53) =     1.302171d0 !     Bismuth -       Iodine
      xfac(83,53) =     0.862377d0 !     Bismuth -       Iodine
      alpb(83,83) =     1.074064d0 !     Bismuth -      Bismuth
      xfac(83,83) =     1.168214d0 !     Bismuth -      Bismuth
 !
      alpb(87, 7) =     2.218810d0 !    Francium -     Nitrogen
      xfac(87, 7) =     1.012630d0 !    Francium -     Nitrogen
      alpb(87, 9) =     2.218810d0 !    Francium -     Fluorine
      xfac(87, 9) =     1.012630d0 !    Francium -     Fluorine
      alpb(87,17) =     1.579660d0 !    Francium -     Chlorine
      xfac(87,17) =     0.761560d0 !    Francium -     Chlorine
      alpb(87,87) =     1.579660d0 !    Francium -     Francium
      xfac(87,87) =     0.761560d0 !    Francium -     Francium
    end subroutine alpb_and_xfac_pm8
  end module Parameters_for_PM8_C
