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

  module Parameters_for_PM7_C
    double precision, dimension(107) :: uss7, upp7, udd7, zs7, zp7, zd7, betas7, &
    betap7, betad7, gss7, gsp7, gpp7, gp27, hsp7, polvo7, poc_7, &
    zsn7, zpn7, zdn7, f0sd7, g2sd7, alp7, &
    CPE_Zet7, CPE_Z07, CPE_B7, CPE_Xlo7, CPE_Xhi7
    double precision :: v_par7(60)
    double precision, dimension(107,4) :: gues71, gues72, gues73
!
!                    Data for Element   1         Hydrogen
!
      data     uss7(  1)/       -11.070112D0/
      data   betas7(  1)/        -8.389745D0/
      data      zs7(  1)/         1.260237D0/
      data     gss7(  1)/        14.149656D0/
      data   polvo7(  1)/         0.229769D0/
      data gues71(  1,1)/         0.177854D0/
      data gues72(  1,1)/         1.428710D0/
      data gues73(  1,1)/         0.991324D0/
!
!                    Data for Element   2           Helium
!
      data     uss7(  2)/       -31.770969D0/
      data     upp7(  2)/        -5.856382D0/
      data   betas7(  2)/       -58.903774D0/
      data   betap7(  2)/       -37.039974D0/
      data      zs7(  2)/         3.313204D0/
      data      zp7(  2)/         3.657133D0/
      data     gss7(  2)/         9.445299D0/
      data     gsp7(  2)/        11.201419D0/
      data     gpp7(  2)/         9.214548D0/
      data     gp27(  2)/        13.046115D0/
      data     hsp7(  2)/         0.299954D0/
!
!                    Data for Element   3          Lithium
!
      data     uss7(  3)/        -4.804124D0/
      data     upp7(  3)/        -2.450842D0/
      data   betas7(  3)/        -2.082310D0/
      data   betap7(  3)/       -27.085547D0/
      data      zs7(  3)/         0.804974D0/
      data      zp7(  3)/         6.027530D0/
      data     gss7(  3)/         9.175811D0/
      data     gsp7(  3)/        16.614419D0/
      data     gpp7(  3)/        14.193195D0/
      data     gp27(  3)/        11.289123D0/
      data     hsp7(  3)/         3.533317D0/
!
!                    Data for Element   4        Beryllium
!
      data     uss7(  4)/       -17.427477D0/
      data     upp7(  4)/       -14.843910D0/
      data   betas7(  4)/        -3.965129D0/
      data   betap7(  4)/        -8.623194D0/
      data      zs7(  4)/         1.036199D0/
      data      zp7(  4)/         1.764629D0/
      data     gss7(  4)/         9.590009D0/
      data     gsp7(  4)/         9.878338D0/
      data     gpp7(  4)/         8.195145D0/
      data     gp27(  4)/        10.136549D0/
      data     hsp7(  4)/         0.966381D0/
!
!                    Data for Element   5            Boron
!
      data     uss7(  5)/       -26.613990D0/
      data     upp7(  5)/       -23.468278D0/
      data   betas7(  5)/        -7.509528D0/
      data   betap7(  5)/        -3.775165D0/
      data      zs7(  5)/         1.560481D0/
      data      zp7(  5)/         1.449712D0/
      data     gss7(  5)/         6.418667D0/
      data     gsp7(  5)/        10.191243D0/
      data     gpp7(  5)/         5.675076D0/
      data     gp27(  5)/         6.689156D0/
      data     hsp7(  5)/         0.942923D0/
!
!                    Data for Element   6           Carbon
!
      data     uss7(  6)/       -51.372620D0/
      data     upp7(  6)/       -40.135421D0/
      data   betas7(  6)/       -14.414930D0/
      data   betap7(  6)/        -7.893717D0/
      data      zs7(  6)/         1.942244D0/
      data      zp7(  6)/         1.708723D0/
      data     gss7(  6)/        12.347323D0/
      data     gsp7(  6)/        11.932801D0/
      data     gpp7(  6)/        10.452226D0/
      data     gp27(  6)/         9.385498D0/
      data     hsp7(  6)/         0.802634D0/
      data   polvo7(  6)/         0.935524D0/
      data gues71(  6,1)/         0.045888D0/
      data gues72(  6,1)/         5.037055D0/
      data gues73(  6,1)/         1.588715D0/
!
!                    Data for Element   7         Nitrogen
!
      data     uss7(  7)/       -61.623260D0/
      data     upp7(  7)/       -48.949816D0/
      data   betas7(  7)/       -22.110958D0/
      data   betap7(  7)/       -15.464433D0/
      data      zs7(  7)/         2.354344D0/
      data      zp7(  7)/         2.028288D0/
      data     gss7(  7)/        11.881002D0/
      data     gsp7(  7)/         9.572117D0/
      data     gpp7(  7)/        12.210140D0/
      data     gp27(  7)/        10.299947D0/
      data     hsp7(  7)/         2.977433D0/
      data   polvo7(  7)/         0.551046D0/
      data gues71(  7,1)/         0.015143D0/
      data gues72(  7,1)/         4.734148D0/
      data gues73(  7,1)/         1.516714D0/
!
!                    Data for Element   8           Oxygen
!
      data     uss7(  8)/       -96.087736D0/
      data     upp7(  8)/       -71.087738D0/
      data   betas7(  8)/       -67.776339D0/
      data   betap7(  8)/       -20.267981D0/
      data      zs7(  8)/         5.972309D0/
      data      zp7(  8)/         2.349017D0/
      data     gss7(  8)/        14.955121D0/
      data     gsp7(  8)/        16.088521D0/
      data     gpp7(  8)/        12.403337D0/
      data     gp27(  8)/        10.499706D0/
      data     hsp7(  8)/         5.028656D0/
      data   polvo7(  8)/         0.312052D0/
      data gues71(  8,1)/        -0.016243D0/
      data gues72(  8,1)/         1.871970D0/
      data gues73(  8,1)/         1.844360D0/
!
!                    Data for Element   9         Fluorine
!
      data     uss7(  9)/      -137.395658D0/
      data     upp7(  9)/       -98.051239D0/
      data   betas7(  9)/       -69.725688D0/
      data   betap7(  9)/       -29.746124D0/
      data      zs7(  9)/         6.070030D0/
      data      zp7(  9)/         2.930631D0/
      data     gss7(  9)/        13.745134D0/
      data     gsp7(  9)/        17.991219D0/
      data     gpp7(  9)/         9.188912D0/
      data     gp27(  9)/        12.322503D0/
      data     hsp7(  9)/         2.906187D0/
      data   polvo7(  9)/         0.121273D0/
!
!                    Data for Element  10             Neon
!
      data     uss7( 10)/        -2.978729D0/
      data     upp7( 10)/       -85.441118D0/
      data   betas7( 10)/       -69.793475D0/
      data   betap7( 10)/       -33.261962D0/
      data      zs7( 10)/         6.000148D0/
      data      zp7( 10)/         3.834528D0/
      data     gss7( 10)/        19.999574D0/
      data     gsp7( 10)/        16.896951D0/
      data     gpp7( 10)/         8.963560D0/
      data     gp27( 10)/        16.027799D0/
      data     hsp7( 10)/         1.779280D0/
!
!                    Data for Element  11           Sodium
!
      data     uss7( 11)/        -5.815476D0/
      data     upp7( 11)/        -3.731003D0/
      data   betas7( 11)/         8.483380D0/
      data   betap7( 11)/        -5.735680D0/
      data      zs7( 11)/         1.666701D0/
      data      zp7( 11)/         1.397571D0/
      data     gss7( 11)/        20.011368D0/
      data     gsp7( 11)/        20.020053D0/
      data     gpp7( 11)/        12.820792D0/
      data     gp27( 11)/        19.015416D0/
      data     hsp7( 11)/         5.020547D0/
!
!                    Data for Element  12        Magnesium
!
      data     uss7( 12)/       -14.858681D0/
      data     upp7( 12)/       -12.451227D0/
      data   betas7( 12)/       -12.576970D0/
      data   betap7( 12)/        -0.702739D0/
      data      zs7( 12)/         1.170297D0/
      data      zp7( 12)/         1.840439D0/
      data     gss7( 12)/         7.480635D0/
      data     gsp7( 12)/         9.602125D0/
      data     gpp7( 12)/         8.869755D0/
      data     gp27( 12)/         6.241718D0/
      data     hsp7( 12)/         0.992746D0/
!
!                    Data for Element  13         Aluminum
!
      data     uss7( 13)/       -32.518856D0/
      data     upp7( 13)/       -24.873064D0/
      data     udd7( 13)/       -31.418925D0/
      data   betas7( 13)/         6.109627D0/
      data   betap7( 13)/        -2.986557D0/
      data   betad7( 13)/       -28.937998D0/
      data      zs7( 13)/         1.232599D0/
      data      zp7( 13)/         1.219336D0/
      data      zd7( 13)/         1.617502D0/
      data     zsn7( 13)/         2.346908D0/
      data     zpn7( 13)/         1.529050D0/
      data     zdn7( 13)/         3.682742D0/
      data     alp7( 13)/         5.341685D0/
      data     gss7( 13)/        10.347944D0/
      data     gsp7( 13)/         9.180517D0/
      data     gpp7( 13)/         7.181623D0/
      data     gp27( 13)/         5.626787D0/
      data     hsp7( 13)/         0.920418D0/
!
!                    Data for Element  14          Silicon
!
      data     uss7( 14)/       -41.586357D0/
      data     upp7( 14)/       -36.694055D0/
      data     udd7( 14)/       -16.775635D0/
      data   betas7( 14)/       -10.755885D0/
      data   betap7( 14)/        -3.922152D0/
      data   betad7( 14)/        -4.736877D0/
      data      zs7( 14)/         1.433994D0/
      data      zp7( 14)/         1.671776D0/
      data      zd7( 14)/         1.221915D0/
      data     zsn7( 14)/         2.002570D0/
      data     zpn7( 14)/         0.818377D0/
      data     zdn7( 14)/         2.591238D0/
      data     gss7( 14)/         8.159128D0/
      data     gsp7( 14)/        11.213512D0/
      data     gpp7( 14)/         8.521933D0/
      data     gp27( 14)/         8.493112D0/
      data     hsp7( 14)/         0.959479D0/
      data   polvo7( 14)/         3.453040D0/
!
!                    Data for Element  15       Phosphorus
!
      data     uss7( 15)/       -66.928393D0/
      data     upp7( 15)/       -30.435889D0/
      data     udd7( 15)/        -9.320605D0/
      data   betas7( 15)/       -45.578128D0/
      data   betap7( 15)/       -12.487910D0/
      data   betad7( 15)/       -34.182622D0/
      data      zs7( 15)/         2.257933D0/
      data      zp7( 15)/         1.555172D0/
      data      zd7( 15)/         1.235995D0/
      data     zsn7( 15)/         4.925330D0/
      data     zpn7( 15)/         1.040649D0/
      data     zdn7( 15)/        12.110811D0/
      data     gss7( 15)/         3.279518D0/
      data     gsp7( 15)/         7.561540D0/
      data     gpp7( 15)/         5.922736D0/
      data     gp27( 15)/         3.990430D0/
      data     hsp7( 15)/         1.699233D0/
      data   polvo7( 15)/         2.858410D0/
!
!                    Data for Element  16           Sulfur
!
      data     uss7( 16)/       -51.157757D0/
      data     upp7( 16)/       -40.352643D0/
      data     udd7( 16)/       -48.529935D0/
      data   betas7( 16)/       -11.422550D0/
      data   betap7( 16)/        -7.191896D0/
      data   betad7( 16)/       -10.695329D0/
      data      zs7( 16)/         2.046153D0/
      data      zp7( 16)/         1.807678D0/
      data      zd7( 16)/         3.510309D0/
      data     zsn7( 16)/         1.131343D0/
      data     zpn7( 16)/         0.823803D0/
      data     zdn7( 16)/         2.296065D0/
      data     gss7( 16)/         8.728478D0/
      data     gsp7( 16)/         6.483871D0/
      data     gpp7( 16)/         7.357401D0/
      data     gp27( 16)/         6.875448D0/
      data     hsp7( 16)/         3.012199D0/
      data   polvo7( 16)/         1.726420D0/
!
!                    Data for Element  17         Chlorine
!
      data     uss7( 17)/       -68.847511D0/
      data     upp7( 17)/       -56.857991D0/
      data     udd7( 17)/       -49.258143D0/
      data   betas7( 17)/        -2.893931D0/
      data   betap7( 17)/       -13.528255D0/
      data   betad7( 17)/         1.888153D0/
      data      zs7( 17)/         2.223076D0/
      data      zp7( 17)/         2.264466D0/
      data      zd7( 17)/         0.949994D0/
      data     zsn7( 17)/         1.992900D0/
      data     zpn7( 17)/         1.874460D0/
      data     zdn7( 17)/         5.469221D0/
      data     gss7( 17)/        11.540894D0/
      data     gsp7( 17)/         9.032229D0/
      data     gpp7( 17)/         8.457262D0/
      data     gp27( 17)/         7.988115D0/
      data     hsp7( 17)/         5.000141D0/
      data   polvo7( 17)/         1.617180D0/
!
!                    Data for Element  18            Argon
!
      data     uss7( 18)/        -7.797931D0/
      data     upp7( 18)/       -83.211487D0/
      data   betas7( 18)/        -8.839842D0/
      data   betap7( 18)/       -28.427303D0/
      data      zs7( 18)/         6.000272D0/
      data      zp7( 18)/         5.949170D0/
      data     gss7( 18)/        17.858776D0/
      data     gsp7( 18)/         4.168451D0/
      data     gpp7( 18)/        11.852500D0/
      data     gp27( 18)/        15.669543D0/
      data     hsp7( 18)/         4.574549D0/
!
!                    Data for Element  19        Potassium
!
      data     uss7( 19)/        -4.888065D0/
      data     upp7( 19)/        -3.763457D0/
      data   betas7( 19)/        10.013029D0/
      data   betap7( 19)/        -2.882668D0/
      data      zs7( 19)/         5.422018D0/
      data      zp7( 19)/         1.471023D0/
      data     gss7( 19)/        19.497974D0/
      data     gsp7( 19)/         4.674636D0/
      data     gpp7( 19)/         4.339481D0/
      data     gp27( 19)/         5.981455D0/
      data     hsp7( 19)/         1.092988D0/
!
!                    Data for Element  20          Calcium
!
      data     uss7( 20)/       -13.503503D0/
      data     upp7( 20)/       -10.559344D0/
      data   betas7( 20)/       -11.696053D0/
      data   betap7( 20)/         4.968210D0/
      data      zs7( 20)/         1.477988D0/
      data      zp7( 20)/         2.220194D0/
      data     gss7( 20)/         7.914200D0/
      data     gsp7( 20)/         6.712903D0/
      data     gpp7( 20)/         4.997910D0/
      data     gp27( 20)/         4.995881D0/
      data     hsp7( 20)/         1.170905D0/
!
!                    Data for Element  21         Scandium
!
      data     uss7( 21)/       -19.383239D0/
      data     upp7( 21)/       -15.936628D0/
      data     udd7( 21)/       -20.365590D0/
      data   betas7( 21)/       -16.127750D0/
      data   betap7( 21)/        -4.714646D0/
      data   betad7( 21)/        -8.631714D0/
      data      zs7( 21)/         1.794897D0/
      data      zp7( 21)/         2.174934D0/
      data      zd7( 21)/         5.992860D0/
      data     zsn7( 21)/         1.314009D0/
      data     zpn7( 21)/         1.020629D0/
      data     zdn7( 21)/         1.437857D0/
      data     alp7( 21)/         0.991198D0/
      data     gss7( 21)/         7.183554D0/
      data     gsp7( 21)/         6.188166D0/
      data     gpp7( 21)/         6.079855D0/
      data     gp27( 21)/         5.329586D0/
      data     hsp7( 21)/         1.340355D0/
      data    poc_7( 21)/         1.070880D0/
      data    f0sd7( 21)/         8.096837D0/
      data    g2sd7( 21)/         3.531412D0/
!
!                    Data for Element  22         Titanium
!
      data     uss7( 22)/       -26.608414D0/
      data     upp7( 22)/       -23.616842D0/
      data     udd7( 22)/       -28.876758D0/
      data   betas7( 22)/        -5.411644D0/
      data   betap7( 22)/        -4.838856D0/
      data   betad7( 22)/         0.774574D0/
      data      zs7( 22)/         1.448579D0/
      data      zp7( 22)/         1.940695D0/
      data      zd7( 22)/         1.093648D0/
      data     zsn7( 22)/         1.078295D0/
      data     zpn7( 22)/         4.663707D0/
      data     zdn7( 22)/         0.954258D0/
      data     gss7( 22)/         5.894930D0/
      data     gsp7( 22)/         7.330203D0/
      data     gpp7( 22)/        27.781556D0/
      data     gp27( 22)/        24.353243D0/
      data     hsp7( 22)/         0.044555D0/
      data    f0sd7( 22)/         6.384127D0/
      data    g2sd7( 22)/         3.488564D0/
!
!                    Data for Element  23         Vanadium
!
      data     uss7( 23)/       -32.598954D0/
      data     upp7( 23)/       -20.496422D0/
      data     udd7( 23)/       -43.169867D0/
      data   betas7( 23)/        -4.628385D0/
      data   betap7( 23)/        -3.039568D0/
      data   betad7( 23)/        -3.704203D0/
      data      zs7( 23)/         6.051795D0/
      data      zp7( 23)/         2.249871D0/
      data      zd7( 23)/         1.087345D0/
      data     zsn7( 23)/         1.215500D0/
      data     zpn7( 23)/         0.877260D0/
      data     zdn7( 23)/         1.512555D0/
      data     gss7( 23)/         6.645015D0/
      data     gsp7( 23)/         5.436952D0/
      data     gpp7( 23)/         5.225810D0/
      data     gp27( 23)/         4.580932D0/
      data     hsp7( 23)/         1.092636D0/
      data    f0sd7( 23)/         6.560730D0/
      data    g2sd7( 23)/         1.196816D0/
!
!                    Data for Element  24         Chromium
!
      data     uss7( 24)/       -41.077064D0/
      data     upp7( 24)/       -19.350873D0/
      data     udd7( 24)/       -80.190851D0/
      data   betas7( 24)/       -13.781066D0/
      data   betap7( 24)/         0.735757D0/
      data   betad7( 24)/        -6.372908D0/
      data      zs7( 24)/         2.838413D0/
      data      zp7( 24)/         1.379560D0/
      data      zd7( 24)/         1.188729D0/
      data     zsn7( 24)/         2.174521D0/
      data     zpn7( 24)/         4.770642D0/
      data     zdn7( 24)/         2.141579D0/
      data     gss7( 24)/        11.887886D0/
      data     gsp7( 24)/        14.518298D0/
      data     gpp7( 24)/        28.418564D0/
      data     gp27( 24)/        24.911643D0/
      data     hsp7( 24)/         1.187458D0/
      data    f0sd7( 24)/         7.511007D0/
      data    g2sd7( 24)/         2.622589D0/
!
!                    Data for Element  25        Manganese
!
      data     uss7( 25)/       -42.374682D0/
      data     upp7( 25)/       -18.304981D0/
      data     udd7( 25)/       -54.430991D0/
      data   betas7( 25)/       -19.986721D0/
      data   betap7( 25)/       -51.153604D0/
      data   betad7( 25)/       -28.049908D0/
      data      zs7( 25)/         1.666440D0/
      data      zp7( 25)/         2.078735D0/
      data      zd7( 25)/         2.897070D0/
      data     zsn7( 25)/         1.299761D0/
      data     zpn7( 25)/         4.059245D0/
      data     zdn7( 25)/         1.146085D0/
      data     gss7( 25)/         7.105662D0/
      data     gsp7( 25)/         8.807648D0/
      data     gpp7( 25)/        24.180795D0/
      data     gp27( 25)/        21.196825D0/
      data     hsp7( 25)/         0.221872D0/
      data    f0sd7( 25)/         4.784190D0/
      data    g2sd7( 25)/         2.008311D0/
!
!                    Data for Element  26             Iron
!
      data     uss7( 26)/       -74.715611D0/
      data     upp7( 26)/       -56.758188D0/
      data     udd7( 26)/       -90.918476D0/
      data   betas7( 26)/        -4.365430D0/
      data   betap7( 26)/        -4.256080D0/
      data   betad7( 26)/       -12.531631D0/
      data      zs7( 26)/         1.157576D0/
      data      zp7( 26)/         2.737621D0/
      data      zd7( 26)/         1.860792D0/
      data     zsn7( 26)/         2.223065D0/
      data     zpn7( 26)/         1.314405D0/
      data     zdn7( 26)/         1.769722D0/
      data     gss7( 26)/        12.153271D0/
      data     gsp7( 26)/         8.511068D0/
      data     gpp7( 26)/         7.829869D0/
      data     gp27( 26)/         6.863644D0/
      data     hsp7( 26)/         1.267977D0/
      data    poc_7( 26)/         0.993526D0/
      data    f0sd7( 26)/         9.314037D0/
      data    g2sd7( 26)/         1.970401D0/
!
!                    Data for Element  27           Cobalt
!
      data     uss7( 27)/       -37.720682D0/
      data     upp7( 27)/        -0.230340D0/
      data     udd7( 27)/       -85.185900D0/
      data   betas7( 27)/       -11.175136D0/
      data   betap7( 27)/       -18.331339D0/
      data   betad7( 27)/        -5.935777D0/
      data      zs7( 27)/         1.789441D0/
      data      zp7( 27)/         1.531664D0/
      data      zd7( 27)/         1.951497D0/
      data     zsn7( 27)/         1.710796D0/
      data     zpn7( 27)/         0.928007D0/
      data     zdn7( 27)/         1.563753D0/
      data     gss7( 27)/         9.352749D0/
      data     gsp7( 27)/         6.087093D0/
      data     gpp7( 27)/         5.528108D0/
      data     gp27( 27)/         4.845926D0/
      data     hsp7( 27)/         0.763102D0/
      data    poc_7( 27)/         1.433458D0/
      data    f0sd7( 27)/         3.045500D0/
      data    g2sd7( 27)/         1.015102D0/
!
!                    Data for Element  28           Nickel
!
      data     uss7( 28)/       -55.503570D0/
      data     upp7( 28)/       -30.601744D0/
      data     udd7( 28)/       -68.610896D0/
      data   betas7( 28)/       -15.417178D0/
      data   betap7( 28)/       -21.305796D0/
      data   betad7( 28)/        -4.094535D0/
      data      zs7( 28)/         1.708340D0/
      data      zp7( 28)/         2.000099D0/
      data      zd7( 28)/         5.698724D0/
      data     zsn7( 28)/         1.177087D0/
      data     zpn7( 28)/         1.013217D0/
      data     zdn7( 28)/         1.017987D0/
      data     gss7( 28)/         6.435016D0/
      data     gsp7( 28)/         5.921995D0/
      data     gpp7( 28)/         6.035702D0/
      data     gp27( 28)/         5.290881D0/
      data     hsp7( 28)/         1.379687D0/
      data    poc_7( 28)/         2.208500D0/
      data    f0sd7( 28)/         5.492550D0/
      data    g2sd7( 28)/         2.469437D0/
!
!                    Data for Element  29           Copper
!
      data     uss7( 29)/       -55.174441D0/
      data     upp7( 29)/         3.200458D0/
      data     udd7( 29)/      -118.258961D0/
      data   betas7( 29)/       -11.801588D0/
      data   betap7( 29)/       -37.165178D0/
      data   betad7( 29)/       -14.652492D0/
      data      zs7( 29)/         1.735325D0/
      data      zp7( 29)/         3.219976D0/
      data      zd7( 29)/         6.013523D0/
      data     zsn7( 29)/         2.419271D0/
      data     zpn7( 29)/         0.302125D0/
      data     zdn7( 29)/         1.678203D0/
      data     gss7( 29)/        13.225910D0/
      data     gsp7( 29)/         2.055274D0/
      data     gpp7( 29)/         1.799749D0/
      data     gp27( 29)/         1.577656D0/
      data     hsp7( 29)/         0.000420D0/
      data    f0sd7( 29)/         5.160900D0/
      data    g2sd7( 29)/         2.792359D0/
!
!                    Data for Element  30             Zinc
!
      data     uss7( 30)/       -16.700035D0/
      data     upp7( 30)/       -14.844247D0/
      data   betas7( 30)/       -16.770975D0/
      data   betap7( 30)/         2.907797D0/
      data      zs7( 30)/         1.560140D0/
      data      zp7( 30)/         1.915631D0/
      data     gss7( 30)/         6.421475D0/
      data     gsp7( 30)/        10.243652D0/
      data     gpp7( 30)/        20.001326D0/
      data     gp27( 30)/        16.126802D0/
      data     hsp7( 30)/         0.983644D0/
!
!                    Data for Element  31          Gallium
!
      data     uss7( 31)/       -30.812730D0/
      data     upp7( 31)/       -22.498885D0/
      data   betas7( 31)/       -15.082480D0/
      data   betap7( 31)/        -0.938845D0/
      data      zs7( 31)/         1.913326D0/
      data      zp7( 31)/         1.811217D0/
      data     gss7( 31)/         9.436450D0/
      data     gsp7( 31)/         9.189262D0/
      data     gpp7( 31)/         5.480436D0/
      data     gp27( 31)/         6.991064D0/
      data     hsp7( 31)/         0.970992D0/
!
!                    Data for Element  32        Germanium
!
      data     uss7( 32)/       -35.694620D0/
      data     upp7( 32)/       -29.273804D0/
      data   betas7( 32)/       -18.071730D0/
      data   betap7( 32)/        -1.563157D0/
      data      zs7( 32)/         2.762845D0/
      data      zp7( 32)/         1.531131D0/
      data     gss7( 32)/         4.991616D0/
      data     gsp7( 32)/         9.108444D0/
      data     gpp7( 32)/         6.693916D0/
      data     gp27( 32)/         5.914950D0/
      data     hsp7( 32)/         0.801760D0/
!
!                    Data for Element  33          Arsenic
!
      data     uss7( 33)/       -41.523302D0/
      data     upp7( 33)/       -36.959219D0/
      data     udd7( 33)/       -35.859071D0/
      data   betas7( 33)/       -17.440295D0/
      data   betap7( 33)/        -6.603566D0/
      data   betad7( 33)/        -2.445554D0/
      data      zs7( 33)/         3.213850D0/
      data      zp7( 33)/         1.628384D0/
      data      zd7( 33)/         3.314358D0/
      data     zsn7( 33)/         0.916221D0/
      data     zpn7( 33)/         1.115722D0/
      data     zdn7( 33)/         2.137809D0/
      data     gss7( 33)/         8.088350D0/
      data     gsp7( 33)/         7.457692D0/
      data     gpp7( 33)/         8.918517D0/
      data     gp27( 33)/         7.102449D0/
      data     hsp7( 33)/         0.969630D0/
!
!                    Data for Element  34         Selenium
!
      data     uss7( 34)/       -47.218303D0/
      data     upp7( 34)/       -35.588294D0/
      data   betas7( 34)/        -9.614158D0/
      data   betap7( 34)/        -6.121302D0/
      data      zs7( 34)/         2.751130D0/
      data      zp7( 34)/         1.901764D0/
      data     gss7( 34)/         4.895424D0/
      data     gsp7( 34)/         6.792977D0/
      data     gpp7( 34)/         5.775063D0/
      data     gp27( 34)/         5.578480D0/
      data     hsp7( 34)/         3.152775D0/
!
!                    Data for Element  35          Bromine
!
      data     uss7( 35)/       -49.141354D0/
      data     upp7( 35)/       -48.274409D0/
      data     udd7( 35)/         2.677328D0/
      data   betas7( 35)/       -32.458894D0/
      data   betap7( 35)/       -10.270309D0/
      data   betad7( 35)/       -19.977175D0/
      data      zs7( 35)/         3.725480D0/
      data      zp7( 35)/         2.242318D0/
      data      zd7( 35)/         1.591034D0/
      data     zsn7( 35)/        10.522069D0/
      data     zpn7( 35)/         9.531017D0/
      data     zdn7( 35)/         5.776829D0/
      data     gss7( 35)/         8.131665D0/
      data     gsp7( 35)/         4.285572D0/
      data     gpp7( 35)/         8.056519D0/
      data     gp27( 35)/         7.520115D0/
      data     hsp7( 35)/         1.567275D0/
      data   polvo7( 35)/         2.645090D0/
!
!                    Data for Element  36          Krypton
!
      data     uss7( 36)/         8.535384D0/
      data     upp7( 36)/       -80.484321D0/
      data   betas7( 36)/        -2.727088D0/
      data   betap7( 36)/       -16.142951D0/
      data      zs7( 36)/         1.312248D0/
      data      zp7( 36)/         4.491371D0/
      data     gss7( 36)/        19.999857D0/
      data     gsp7( 36)/         1.175304D0/
      data     gpp7( 36)/         9.174784D0/
      data     gp27( 36)/        14.926948D0/
      data     hsp7( 36)/         0.299867D0/
!
!                    Data for Element  37         Rubidium
!
      data     uss7( 37)/        -4.120962D0/
      data     upp7( 37)/        -0.633300D0/
      data   betas7( 37)/        -8.442947D0/
      data   betap7( 37)/         4.853952D0/
      data      zs7( 37)/         1.314831D0/
      data      zp7( 37)/         6.015581D0/
      data     gss7( 37)/        11.892047D0/
      data     gsp7( 37)/         3.477383D0/
      data     gpp7( 37)/         6.000901D0/
      data     gp27( 37)/         6.008182D0/
      data     hsp7( 37)/         0.998126D0/
!
!                    Data for Element  38        Strontium
!
      data     uss7( 38)/       -10.693066D0/
      data     upp7( 38)/        -8.539218D0/
      data   betas7( 38)/        -4.904378D0/
      data   betap7( 38)/         8.809297D0/
      data      zs7( 38)/         2.092264D0/
      data      zp7( 38)/         3.314082D0/
      data     gss7( 38)/         6.494973D0/
      data     gsp7( 38)/         4.045506D0/
      data     gpp7( 38)/         2.547611D0/
      data     gp27( 38)/         4.121201D0/
      data     hsp7( 38)/         1.102790D0/
!
!                    Data for Element  39          Yttrium
!
      data     uss7( 39)/       -17.035117D0/
      data     upp7( 39)/       -16.168689D0/
      data     udd7( 39)/       -16.354811D0/
      data   betas7( 39)/       -10.513848D0/
      data   betap7( 39)/       -11.341408D0/
      data   betad7( 39)/       -10.701025D0/
      data      zs7( 39)/         1.605083D0/
      data      zp7( 39)/         2.131069D0/
      data      zd7( 39)/         6.021645D0/
      data     zsn7( 39)/         1.186263D0/
      data     zpn7( 39)/         2.244351D0/
      data     zdn7( 39)/         0.911477D0/
      data     gss7( 39)/         5.318448D0/
      data     gsp7( 39)/         6.318665D0/
      data     gpp7( 39)/        11.005171D0/
      data     gp27( 39)/         9.590777D0/
      data     hsp7( 39)/         0.637634D0/
      data    poc_7( 39)/         2.019557D0/
      data    f0sd7( 39)/         6.855595D0/
      data    g2sd7( 39)/         5.889822D0/
!
!                    Data for Element  40        Zirconium
!
      data     uss7( 40)/       -18.679203D0/
      data     upp7( 40)/        -0.049727D0/
      data     udd7( 40)/       -23.299157D0/
      data   betas7( 40)/         1.688872D0/
      data   betap7( 40)/         2.330045D0/
      data   betad7( 40)/        -4.552268D0/
      data      zs7( 40)/         1.373517D0/
      data      zp7( 40)/         1.141705D0/
      data      zd7( 40)/         1.618769D0/
      data     zsn7( 40)/         1.082243D0/
      data     zpn7( 40)/         2.978817D0/
      data     zdn7( 40)/         1.417227D0/
      data     gss7( 40)/         4.852089D0/
      data     gsp7( 40)/         5.870299D0/
      data     gpp7( 40)/        14.606623D0/
      data     gp27( 40)/        12.729368D0/
      data     hsp7( 40)/         0.151363D0/
      data    f0sd7( 40)/         4.737734D0/
      data    g2sd7( 40)/         2.141620D0/
!
!                    Data for Element  41          Niobium
!
      data     uss7( 41)/       -33.497110D0/
      data     upp7( 41)/       -34.762698D0/
      data     udd7( 41)/       -44.819149D0/
      data   betas7( 41)/       -23.566737D0/
      data   betap7( 41)/        -1.623945D0/
      data   betad7( 41)/        -7.421668D0/
      data      zs7( 41)/         2.761686D0/
      data      zp7( 41)/         5.999062D0/
      data      zd7( 41)/         1.611677D0/
      data     zsn7( 41)/         1.429235D0/
      data     zpn7( 41)/         2.911794D0/
      data     zdn7( 41)/         1.950434D0/
      data     gss7( 41)/         6.407780D0/
      data     gsp7( 41)/         7.659444D0/
      data     gpp7( 41)/        14.277976D0/
      data     gp27( 41)/        12.442959D0/
      data     hsp7( 41)/         0.619230D0/
      data    f0sd7( 41)/         6.393769D0/
      data    g2sd7( 41)/         1.759636D0/
!
!                    Data for Element  42       Molybdenum
!
      data     uss7( 42)/       -51.662768D0/
      data     upp7( 42)/        46.059429D0/
      data     udd7( 42)/       -57.269405D0/
      data   betas7( 42)/         6.685073D0/
      data   betap7( 42)/         5.485123D0/
      data   betad7( 42)/       -13.146960D0/
      data      zs7( 42)/         1.595399D0/
      data      zp7( 42)/         1.426575D0/
      data      zd7( 42)/         1.787748D0/
      data     zsn7( 42)/         1.903541D0/
      data     zpn7( 42)/         1.592195D0/
      data     zdn7( 42)/         1.889678D0/
      data     gss7( 42)/         8.534266D0/
      data     gsp7( 42)/         7.704937D0/
      data     gpp7( 42)/         7.807325D0/
      data     gp27( 42)/         6.803921D0/
      data     hsp7( 42)/         1.787407D0/
      data    f0sd7( 42)/         9.654475D0/
      data    g2sd7( 42)/         2.314954D0/
!
!                    Data for Element  43       Technetium
!
      data     uss7( 43)/       -48.916740D0/
      data     upp7( 43)/       -21.908166D0/
      data     udd7( 43)/       -53.807590D0/
      data   betas7( 43)/       -17.096185D0/
      data   betap7( 43)/       -17.740652D0/
      data   betad7( 43)/        -7.241592D0/
      data      zs7( 43)/         2.104672D0/
      data      zp7( 43)/         2.669984D0/
      data      zd7( 43)/         3.030496D0/
      data     zsn7( 43)/         2.061082D0/
      data     zpn7( 43)/         0.888524D0/
      data     zdn7( 43)/         1.575315D0/
      data     gss7( 43)/         9.240580D0/
      data     gsp7( 43)/         4.795964D0/
      data     gpp7( 43)/         4.356876D0/
      data     gp27( 43)/         3.796926D0/
      data     hsp7( 43)/         0.248174D0/
      data    f0sd7( 43)/         7.521148D0/
      data    g2sd7( 43)/         3.106149D0/
!
!                    Data for Element  44        Ruthenium
!
      data     uss7( 44)/       -41.151429D0/
      data     upp7( 44)/       -42.965344D0/
      data     udd7( 44)/       -45.714719D0/
      data   betas7( 44)/        -4.989393D0/
      data   betap7( 44)/       -10.778690D0/
      data   betad7( 44)/         1.216566D0/
      data      zs7( 44)/         1.605646D0/
      data      zp7( 44)/         4.580820D0/
      data      zd7( 44)/         1.244578D0/
      data     zsn7( 44)/         1.172546D0/
      data     zpn7( 44)/         1.373361D0/
      data     zdn7( 44)/         1.018114D0/
      data     gss7( 44)/         5.256950D0/
      data     gsp7( 44)/         5.631870D0/
      data     gpp7( 44)/         6.734273D0/
      data     gp27( 44)/         5.868779D0/
      data     hsp7( 44)/         1.326657D0/
      data    f0sd7( 44)/         4.898881D0/
      data    g2sd7( 44)/         2.648488D0/
!
!                    Data for Element  45          Rhodium
!
      data     uss7( 45)/       -24.613157D0/
      data     upp7( 45)/         6.621039D0/
      data     udd7( 45)/       -81.764165D0/
      data   betas7( 45)/        -9.488908D0/
      data   betap7( 45)/        -6.699556D0/
      data   betad7( 45)/        -7.997845D0/
      data      zs7( 45)/         1.591465D0/
      data      zp7( 45)/         4.546046D0/
      data      zd7( 45)/         2.685918D0/
      data     zsn7( 45)/         2.079986D0/
      data     zpn7( 45)/         9.641003D0/
      data     zdn7( 45)/         1.787794D0/
      data     gss7( 45)/         9.325333D0/
      data     gsp7( 45)/        11.318440D0/
      data     gpp7( 45)/        47.274639D0/
      data     gp27( 45)/        41.198863D0/
      data     hsp7( 45)/         0.017585D0/
      data    f0sd7( 45)/         2.230584D0/
      data    g2sd7( 45)/         1.492841D0/
!
!                    Data for Element  46        Palladium
!
      data     uss7( 46)/       -90.670356D0/
      data     upp7( 46)/        45.018147D0/
      data     udd7( 46)/       -94.618031D0/
      data   betas7( 46)/       -18.862423D0/
      data   betap7( 46)/       -18.107010D0/
      data   betad7( 46)/        -3.592862D0/
      data      zs7( 46)/         5.790768D0/
      data      zp7( 46)/         2.169788D0/
      data      zd7( 46)/         1.327661D0/
      data     zsn7( 46)/         1.985663D0/
      data     zpn7( 46)/         0.621281D0/
      data     zdn7( 46)/         1.768258D0/
      data     gss7( 46)/         8.902449D0/
      data     gsp7( 46)/         3.376439D0/
      data     gpp7( 46)/         3.046450D0/
      data     gp27( 46)/         2.654918D0/
      data     hsp7( 46)/         0.043028D0/
      data    f0sd7( 46)/         9.251409D0/
      data    g2sd7( 46)/         1.948722D0/
!
!                    Data for Element  47           Silver
!
      data     uss7( 47)/       -92.280499D0/
      data     upp7( 47)/        29.229985D0/
      data     udd7( 47)/       -82.344865D0/
      data   betas7( 47)/        -9.850776D0/
      data   betap7( 47)/       -29.894728D0/
      data   betad7( 47)/       -63.636331D0/
      data      zs7( 47)/         1.793032D0/
      data      zp7( 47)/         2.528721D0/
      data      zd7( 47)/         3.524808D0/
      data     zsn7( 47)/         1.619764D0/
      data     zpn7( 47)/         0.439729D0/
      data     zdn7( 47)/         1.210202D0/
      data     gss7( 47)/         7.261991D0/
      data     gsp7( 47)/         2.391732D0/
      data     gpp7( 47)/         2.156210D0/
      data     gp27( 47)/         1.879092D0/
      data     hsp7( 47)/         0.014435D0/
      data    f0sd7( 47)/         8.987758D0/
      data    g2sd7( 47)/         4.716654D0/
!
!                    Data for Element  48          Cadmium
!
      data     uss7( 48)/       -18.127987D0/
      data     upp7( 48)/       -13.777839D0/
      data   betas7( 48)/       -23.781665D0/
      data   betap7( 48)/       -11.892060D0/
      data      zs7( 48)/         3.670047D0/
      data      zp7( 48)/         1.857036D0/
      data     gss7( 48)/         8.904816D0/
      data     gsp7( 48)/         9.232666D0/
      data     gpp7( 48)/        11.103045D0/
      data     gp27( 48)/        10.905897D0/
      data     hsp7( 48)/         0.981926D0/
!
!                    Data for Element  49           Indium
!
      data     uss7( 49)/       -26.891944D0/
      data     upp7( 49)/       -28.519053D0/
      data   betas7( 49)/        -0.447307D0/
      data   betap7( 49)/        -4.269337D0/
      data      zs7( 49)/         1.902085D0/
      data      zp7( 49)/         1.940127D0/
      data     gss7( 49)/         6.493621D0/
      data     gsp7( 49)/        12.576468D0/
      data     gpp7( 49)/        10.282533D0/
      data     gp27( 49)/        10.903195D0/
      data     hsp7( 49)/         2.133796D0/
!
!                    Data for Element  50              Tin
!
      data     uss7( 50)/       -33.880164D0/
      data     upp7( 50)/       -39.128186D0/
      data   betas7( 50)/         0.443105D0/
      data   betap7( 50)/        -8.486074D0/
      data      zs7( 50)/         1.959238D0/
      data      zp7( 50)/         1.976146D0/
      data     gss7( 50)/         6.196917D0/
      data     gsp7( 50)/        10.595744D0/
      data     gpp7( 50)/        14.691065D0/
      data     gp27( 50)/        13.501111D0/
      data     hsp7( 50)/         1.234523D0/
!
!                    Data for Element  51         Antimony
!
      data     uss7( 51)/       -42.835901D0/
      data     upp7( 51)/       -19.996258D0/
      data     udd7( 51)/       -20.317174D0/
      data   betas7( 51)/       -13.037071D0/
      data   betap7( 51)/        -6.166480D0/
      data   betad7( 51)/        -9.740725D0/
      data      zs7( 51)/         1.998600D0/
      data      zp7( 51)/         1.887062D0/
      data      zd7( 51)/         1.475516D0/
      data     zsn7( 51)/         2.179206D0/
      data     zpn7( 51)/         0.862318D0/
      data     zdn7( 51)/         4.147596D0/
      data     gss7( 51)/         9.994149D0/
      data     gsp7( 51)/         1.434008D0/
      data     gpp7( 51)/         7.208157D0/
      data     gp27( 51)/         6.212730D0/
      data     hsp7( 51)/         3.566032D0/
!
!                    Data for Element  52        Tellurium
!
      data     uss7( 52)/       -97.416118D0/
      data     upp7( 52)/       -50.000552D0/
      data   betas7( 52)/       -70.028904D0/
      data   betap7( 52)/       -11.183348D0/
      data      zs7( 52)/         3.024819D0/
      data      zp7( 52)/         2.598283D0/
      data     gss7( 52)/        18.350494D0/
      data     gsp7( 52)/        11.255114D0/
      data     gpp7( 52)/         8.695261D0/
      data     gp27( 52)/         7.622556D0/
      data     hsp7( 52)/         3.626912D0/
!
!                    Data for Element  53           Iodine
!
      data     uss7( 53)/       -63.618928D0/
      data     upp7( 53)/       -45.760969D0/
      data     udd7( 53)/         6.109119D0/
      data   betas7( 53)/       -37.373318D0/
      data   betap7( 53)/       -10.174482D0/
      data   betad7( 53)/       -11.807267D0/
      data      zs7( 53)/         3.316202D0/
      data      zp7( 53)/         2.449124D0/
      data      zd7( 53)/         1.716121D0/
      data     zsn7( 53)/         4.000764D0/
      data     zpn7( 53)/         3.993847D0/
      data     zdn7( 53)/         3.946706D0/
      data     gss7( 53)/         7.658717D0/
      data     gsp7( 53)/         8.237228D0/
      data     gpp7( 53)/         5.667030D0/
      data     gp27( 53)/         5.661068D0/
      data     hsp7( 53)/         2.688576D0/
      data   polvo7( 53)/         4.400520D0/
!
!                    Data for Element  54            Xenon
!
      data     uss7( 54)/       -18.964330D0/
      data     upp7( 54)/      -108.181436D0/
      data   betas7( 54)/        -2.718707D0/
      data   betap7( 54)/       -44.936370D0/
      data      zs7( 54)/         3.208788D0/
      data      zp7( 54)/         2.727979D0/
      data     gss7( 54)/        17.906443D0/
      data     gsp7( 54)/         4.106228D0/
      data     gpp7( 54)/         1.716979D0/
      data     gp27( 54)/        18.971469D0/
      data     hsp7( 54)/         4.990194D0/
!
!                    Data for Element  55           Cesium
!
      data     uss7( 55)/        -3.996308D0/
      data     upp7( 55)/        -2.569885D0/
      data   betas7( 55)/       -11.167340D0/
      data   betap7( 55)/         9.691485D0/
      data      zs7( 55)/         1.776064D0/
      data      zp7( 55)/         6.025310D0/
      data     gss7( 55)/        18.164131D0/
      data     gsp7( 55)/         6.920824D0/
      data     gpp7( 55)/        16.792426D0/
      data     gp27( 55)/         8.175881D0/
      data     hsp7( 55)/         4.590034D0/
!
!                    Data for Element  56           Barium
!
      data     uss7( 56)/       -11.571532D0/
      data     upp7( 56)/        -9.917993D0/
      data   betas7( 56)/       -10.914737D0/
      data   betap7( 56)/         9.727920D0/
      data      zs7( 56)/         1.750490D0/
      data      zp7( 56)/         1.968788D0/
      data     gss7( 56)/         7.843618D0/
      data     gsp7( 56)/        19.900648D0/
      data     gpp7( 56)/        20.004643D0/
      data     gp27( 56)/        19.020523D0/
      data     hsp7( 56)/         0.979914D0/
!
!                    Data for Element  57        Lanthanum
!
      data     uss7( 57)/       -15.586927D0/
      data     upp7( 57)/        58.477136D0/
      data     udd7( 57)/       -19.818759D0/
      data   betas7( 57)/       -18.460416D0/
      data   betap7( 57)/       -19.708547D0/
      data   betad7( 57)/         0.849478D0/
      data      zs7( 57)/         3.398968D0/
      data      zp7( 57)/         1.811983D0/
      data      zd7( 57)/         1.894574D0/
      data     zsn7( 57)/         1.187188D0/
      data     zpn7( 57)/         2.542482D0/
      data     zdn7( 57)/         2.306744D0/
      data     gss7( 57)/         4.516349D0/
      data     gsp7( 57)/         5.344112D0/
      data     gpp7( 57)/        10.610725D0/
      data     gp27( 57)/         9.202959D0/
      data     hsp7( 57)/         0.286106D0/
      data    poc_7( 57)/         1.846287D0/
      data    f0sd7( 57)/         7.808849D0/
      data    g2sd7( 57)/         5.958952D0/
!
!                    Data for Element  71         Lutetium
!
      data     uss7( 71)/       -21.914035D0/
      data     upp7( 71)/        54.132176D0/
      data     udd7( 71)/       -24.661582D0/
      data   betas7( 71)/       -26.143720D0/
      data   betap7( 71)/        -9.506888D0/
      data   betad7( 71)/         3.472080D0/
      data      zs7( 71)/         2.327039D0/
      data      zp7( 71)/         6.000335D0/
      data      zd7( 71)/         1.208414D0/
      data     zsn7( 71)/         0.449170D0/
      data     zpn7( 71)/         2.469444D0/
      data     zdn7( 71)/         2.216418D0/
      data     gss7( 71)/         1.708751D0/
      data     gsp7( 71)/         2.037073D0/
      data     gpp7( 71)/        10.305910D0/
      data     gp27( 71)/         8.938585D0/
      data     hsp7( 71)/         0.000293D0/
      data    poc_7( 71)/         5.824175D0/
      data    f0sd7( 71)/         9.700135D0/
      data    g2sd7( 71)/         6.013887D0/
!
!                    Data for Element  72          Hafnium
!
      data     uss7( 72)/       -25.690382D0/
      data     upp7( 72)/        -9.479410D0/
      data     udd7( 72)/       -39.077741D0/
      data   betas7( 72)/        -4.866355D0/
      data   betap7( 72)/       -21.264221D0/
      data   betad7( 72)/       -12.878794D0/
      data      zs7( 72)/         2.854938D0/
      data      zp7( 72)/         3.079458D0/
      data      zd7( 72)/         2.067146D0/
      data     zsn7( 72)/         3.099683D0/
      data     zpn7( 72)/         3.333027D0/
      data     zdn7( 72)/         3.025020D0/
      data     gss7( 72)/        11.791941D0/
      data     gsp7( 72)/        12.198754D0/
      data     gpp7( 72)/        13.909964D0/
      data     gp27( 72)/        12.064475D0/
      data     hsp7( 72)/         3.057466D0/
      data    f0sd7( 72)/         4.020384D0/
      data    g2sd7( 72)/         4.323408D0/
!
!                    Data for Element  73         Tantalum
!
      data     uss7( 73)/       -34.075891D0/
      data     upp7( 73)/        -5.504664D0/
      data     udd7( 73)/       -35.650460D0/
      data   betas7( 73)/       -15.943219D0/
      data   betap7( 73)/         8.985389D0/
      data   betad7( 73)/       -11.508162D0/
      data      zs7( 73)/         4.116264D0/
      data      zp7( 73)/         3.380936D0/
      data      zd7( 73)/         1.755408D0/
      data     zsn7( 73)/         1.011432D0/
      data     zpn7( 73)/         2.139168D0/
      data     zdn7( 73)/         1.685479D0/
      data     gss7( 73)/         3.847731D0/
      data     gsp7( 73)/         4.550506D0/
      data     gpp7( 73)/         8.927545D0/
      data     gp27( 73)/         7.743093D0/
      data     hsp7( 73)/         0.256277D0/
      data    f0sd7( 73)/         7.257766D0/
      data    g2sd7( 73)/         1.619809D0/
!
!                    Data for Element  74         Tungsten
!
      data     uss7( 74)/       -52.048404D0/
      data     upp7( 74)/       -39.590059D0/
      data     udd7( 74)/       -53.556920D0/
      data   betas7( 74)/       -63.148771D0/
      data   betap7( 74)/        -2.737119D0/
      data   betad7( 74)/         1.132748D0/
      data      zs7( 74)/         3.881177D0/
      data      zp7( 74)/         2.044717D0/
      data      zd7( 74)/         1.928901D0/
      data     zsn7( 74)/         3.461491D0/
      data     zpn7( 74)/         1.904387D0/
      data     zdn7( 74)/         2.180340D0/
      data     gss7( 74)/        13.168346D0/
      data     gsp7( 74)/         8.485482D0/
      data     gpp7( 74)/         7.947717D0/
      data     gp27( 74)/         6.893262D0/
      data     hsp7( 74)/         0.826801D0/
      data    f0sd7( 74)/        10.002668D0/
      data    g2sd7( 74)/         3.417555D0/
!
!                    Data for Element  75          Rhenium
!
      data     uss7( 75)/       -41.679545D0/
      data     upp7( 75)/        43.429894D0/
      data     udd7( 75)/       -54.761512D0/
      data   betas7( 75)/         8.467231D0/
      data   betap7( 75)/        -6.468335D0/
      data   betad7( 75)/       -11.136390D0/
      data      zs7( 75)/         2.452162D0/
      data      zp7( 75)/         1.583194D0/
      data      zd7( 75)/         2.414839D0/
      data     zsn7( 75)/         2.433415D0/
      data     zpn7( 75)/         0.838026D0/
      data     zdn7( 75)/         1.921708D0/
      data     gss7( 75)/         9.257297D0/
      data     gsp7( 75)/         3.796290D0/
      data     gpp7( 75)/         3.497395D0/
      data     gp27( 75)/         3.033382D0/
      data     hsp7( 75)/         0.046330D0/
      data    f0sd7( 75)/         5.229236D0/
      data    g2sd7( 75)/         1.821985D0/
!
!                    Data for Element  76           Osmium
!
      data     uss7( 76)/       -65.963764D0/
      data     upp7( 76)/        37.736568D0/
      data     udd7( 76)/       -89.718816D0/
      data   betas7( 76)/       -43.486712D0/
      data   betap7( 76)/       -25.607006D0/
      data   betad7( 76)/        -1.430819D0/
      data      zs7( 76)/         3.094808D0/
      data      zp7( 76)/         2.845232D0/
      data      zd7( 76)/         1.986395D0/
      data     zsn7( 76)/         2.613281D0/
      data     zpn7( 76)/         2.062936D0/
      data     zdn7( 76)/         2.944917D0/
      data     gss7( 76)/         9.941551D0/
      data     gsp7( 76)/         8.617668D0/
      data     gpp7( 76)/         8.609401D0/
      data     gp27( 76)/         7.467158D0/
      data     hsp7( 76)/         1.886033D0/
      data    f0sd7( 76)/         8.758980D0/
      data    g2sd7( 76)/         4.717871D0/
!
!                    Data for Element  77          Iridium
!
      data     uss7( 77)/       -40.856798D0/
      data     upp7( 77)/        -2.270208D0/
      data     udd7( 77)/       -68.020682D0/
      data   betas7( 77)/       -11.770307D0/
      data   betap7( 77)/       -13.487742D0/
      data   betad7( 77)/        -5.642629D0/
      data      zs7( 77)/         1.924564D0/
      data      zp7( 77)/         3.510744D0/
      data      zd7( 77)/         2.437796D0/
      data     zsn7( 77)/         2.108777D0/
      data     zpn7( 77)/         0.618406D0/
      data     zdn7( 77)/         1.826929D0/
      data     gss7( 77)/         8.022296D0/
      data     gsp7( 77)/         2.803574D0/
      data     gpp7( 77)/         2.580839D0/
      data     gp27( 77)/         2.238429D0/
      data     hsp7( 77)/         0.013100D0/
      data    f0sd7( 77)/         3.726074D0/
      data    g2sd7( 77)/         2.747207D0/
!
!                    Data for Element  78         Platinum
!
      data     uss7( 78)/       -55.878758D0/
      data     upp7( 78)/        52.660706D0/
      data     udd7( 78)/       -92.789895D0/
      data   betas7( 78)/       -10.270452D0/
      data   betap7( 78)/        10.016048D0/
      data   betad7( 78)/        -7.705919D0/
      data      zs7( 78)/         2.922551D0/
      data      zp7( 78)/         0.725689D0/
      data      zd7( 78)/         2.158085D0/
      data     zsn7( 78)/         3.083320D0/
      data     zpn7( 78)/        19.427280D0/
      data     zdn7( 78)/         2.233704D0/
      data     gss7( 78)/        11.729692D0/
      data     gsp7( 78)/        13.983535D0/
      data     gpp7( 78)/        81.077279D0/
      data     gp27( 78)/        70.320441D0/
      data     hsp7( 78)/         0.000643D0/
      data    f0sd7( 78)/         4.725137D0/
      data    g2sd7( 78)/         3.459127D0/
!
!                    Data for Element  79             Gold
!
      data     uss7( 79)/       -94.841695D0/
      data     upp7( 79)/       -61.195249D0/
      data     udd7( 79)/      -114.242383D0/
      data   betas7( 79)/       -13.460355D0/
      data   betap7( 79)/       -24.921790D0/
      data   betad7( 79)/       -63.835796D0/
      data      zs7( 79)/         1.904923D0/
      data      zp7( 79)/         2.408005D0/
      data      zd7( 79)/         4.377691D0/
      data     zsn7( 79)/         2.228930D0/
      data     zpn7( 79)/         4.555019D0/
      data     zdn7( 79)/         2.406645D0/
      data     gss7( 79)/         8.479387D0/
      data     gsp7( 79)/        10.011584D0/
      data     gpp7( 79)/        19.009792D0/
      data     gp27( 79)/        16.487689D0/
      data     hsp7( 79)/         0.645273D0/
      data    f0sd7( 79)/         9.054233D0/
      data    g2sd7( 79)/         5.690708D0/
!
!                    Data for Element  80          Mercury
!
      data     uss7( 80)/       -18.205464D0/
      data     upp7( 80)/       -15.280873D0/
      data   betas7( 80)/        -9.919880D0/
      data   betap7( 80)/         2.134696D0/
      data      zs7( 80)/         2.575831D0/
      data      zp7( 80)/         1.955505D0/
      data     gss7( 80)/         7.687301D0/
      data     gsp7( 80)/         8.355361D0/
      data     gpp7( 80)/         5.135379D0/
      data     gp27( 80)/         9.529761D0/
      data     hsp7( 80)/         0.954727D0/
!
!                    Data for Element  81         Thallium
!
      data     uss7( 81)/       -31.112183D0/
      data     upp7( 81)/       -18.547083D0/
      data     udd7( 81)/        10.294626D0/
      data   betas7( 81)/        -2.456566D0/
      data   betap7( 81)/        -4.949902D0/
      data   betad7( 81)/         0.079835D0/
      data      zs7( 81)/         1.903342D0/
      data      zp7( 81)/         2.838647D0/
      data      zd7( 81)/         5.015677D0/
      data     gss7( 81)/        11.438997D0/
      data     gsp7( 81)/         6.598450D0/
      data     gpp7( 81)/         6.054580D0/
      data     gp27( 81)/         5.507624D0/
      data     hsp7( 81)/         0.894757D0/
!
!                    Data for Element  82             Lead
!
      data     uss7( 82)/       -39.446347D0/
      data     upp7( 82)/       -29.348301D0/
      data     udd7( 82)/       -72.584748D0/
      data   betas7( 82)/       -64.174888D0/
      data   betap7( 82)/        -4.631384D0/
      data   betad7( 82)/        -5.319005D0/
      data      zs7( 82)/         4.706006D0/
      data      zp7( 82)/         2.591455D0/
      data     gss7( 82)/         8.368048D0/
      data     gsp7( 82)/         8.606930D0/
      data     gpp7( 82)/         6.431147D0/
      data     gp27( 82)/         6.550076D0/
      data     hsp7( 82)/         0.984819D0/
!
!                    Data for Element  83          Bismuth
!
      data     uss7( 83)/       -36.561343D0/
      data     upp7( 83)/       -30.823167D0/
      data     udd7( 83)/       -19.667431D0/
      data   betas7( 83)/       -63.673960D0/
      data   betap7( 83)/        -6.931981D0/
      data   betad7( 83)/        -8.868066D0/
      data      zs7( 83)/         5.465413D0/
      data      zp7( 83)/         2.037481D0/
      data      zd7( 83)/         2.855400D0/
      data     zsn7( 83)/         4.275828D0/
      data     zpn7( 83)/         3.018252D0/
      data     zdn7( 83)/         4.889868D0/
      data     gss7( 83)/         3.438678D0/
      data     gsp7( 83)/         3.987429D0/
      data     gpp7( 83)/         8.221978D0/
      data     gp27( 83)/         8.183927D0/
      data     hsp7( 83)/         1.610989D0/
!
!                    Data for Element  85         Astatine
!
      data     alp7( 85)/         3.000000D0/
      data     gss7( 85)/        10.000000D0/
!
!                    Data for Element  87         Francium
!
      data     alp7( 87)/         3.000000D0/
      data     gss7( 87)/        10.000000D0/
!
!                    Data for Element  90          Thorium
!
      data     uss7( 90)/       -40.568292D0/
      data     upp7( 90)/       -28.089187D0/
      data   betas7( 90)/        -4.256218D0/
      data   betap7( 90)/        -4.256218D0/
      data      zs7( 90)/         1.435306D0/
      data      zp7( 90)/         1.435306D0/
      data     gss7( 90)/         9.820000D0/
      data     gsp7( 90)/         8.360000D0/
      data     gpp7( 90)/         7.310000D0/
      data     gp27( 90)/         6.540000D0/
      data     hsp7( 90)/         1.320000D0/
!
!                    Data for Element  97        Berkelium
!
      data gues71( 97,1)/         1.185538D0/
      data gues72( 97,1)/         1.167439D0/
      data gues71( 97,2)/         1.373356D0/
      data gues72( 97,2)/         0.786508D0/
      data gues71( 97,3)/         1.450618D0/
      data gues72( 97,3)/         2.001315D0/
!
!                    Data for Element  98          Mithril
!
      data     uss7( 98)/        -3.000000D0/
      data   betas7( 98)/       -99.000000D0/
      data      zs7( 98)/         2.000000D0/
      data     gss7( 98)/        12.000000D0/
!
!                    Data for Element 100       3+ Sparkle
!
      data     alp7(100)/         1.500000D0/
!
!                    Data for Element 101       3- Sparkle
!
      data     alp7(101)/         1.500000D0/
!
!                    Data for Element 102      Capped bond
!
      data   betas7(102)/  -9999999.000000D0/
      data      zs7(102)/         4.000000D0/
      data     gss7(102)/        12.848000D0/
!
!                    Data for Element 103       ++ Sparkle
!
      data     alp7(103)/         1.500000D0/
!
!                    Data for Element 104        + Sparkle
!
      data     alp7(104)/         1.500000D0/
!
!                    Data for Element 105       -- Sparkle
!
      data     alp7(105)/         1.500000D0/
!
!                    Data for Element 106        - Sparkle
!
      data     alp7(106)/         1.500000D0/
!
!
!                     Global parameters
!
!
      data   v_par7(1)/    8.947612d0/  ! Used in ccrep for scalar correction of C-C triple bonds.
      data   v_par7(2)/    6.024265d0/  ! Used in ccrep for exponent correction of C-C triple bonds.
      data   v_par7(3)/   -0.012037d0/  ! Used in ccrep for scalar correction of O-H term.
      data   v_par7(4)/    0.701333d0/  ! Used in ccrep for exponent correction of C-C triple bonds.
      data   v_par7(7)/    0.900000d0/  ! Used in dftd3 to set "s6"  in D3H4
      data   v_par7(8)/   14.000000d0/  ! Used in dftd3 to set "alp" in D3H4
      data   v_par7(9)/    1.561000d0/  ! Used in dftd3 to set "rs6" in D3H4
  contains
  subroutine alpb_and_xfac_pm7
    use parameters_C, only : xfac, alpb
      xfac = 0.d0
      alpb = 0.d0
 !
      alpb( 1, 1) =     4.051163d0 !    Hydrogen -     Hydrogen
      xfac( 1, 1) =     2.845627d0 !    Hydrogen -     Hydrogen
 !
      alpb( 2, 1) =     2.989881d0 !      Helium -     Hydrogen
      xfac( 2, 1) =     2.371199d0 !      Helium -     Hydrogen
      alpb( 2, 2) =     3.783559d0 !      Helium -       Helium
      xfac( 2, 2) =     3.450900d0 !      Helium -       Helium
 !
      alpb( 3, 1) =     1.265105d0 !     Lithium -     Hydrogen
      xfac( 3, 1) =     0.488118d0 !     Lithium -     Hydrogen
      alpb( 3, 2) =     2.982569d0 !     Lithium -       Helium
      xfac( 3, 2) =     8.316732d0 !     Lithium -       Helium
      alpb( 3, 3) =     3.213216d0 !     Lithium -      Lithium
      xfac( 3, 3) =    16.832394d0 !     Lithium -      Lithium
 !
      alpb( 4, 1) =     2.854611d0 !   Beryllium -     Hydrogen
      xfac( 4, 1) =     3.447327d0 !   Beryllium -     Hydrogen
      alpb( 4, 2) =     3.367214d0 !   Beryllium -       Helium
      xfac( 4, 2) =    12.563185d0 !   Beryllium -       Helium
      alpb( 4, 3) =     2.432991d0 !   Beryllium -      Lithium
      xfac( 4, 3) =    11.306301d0 !   Beryllium -      Lithium
      alpb( 4, 4) =     2.042783d0 !   Beryllium -    Beryllium
      xfac( 4, 4) =     1.680607d0 !   Beryllium -    Beryllium
 !
      alpb( 5, 1) =     2.314226d0 !       Boron -     Hydrogen
      xfac( 5, 1) =     1.135897d0 !       Boron -     Hydrogen
      alpb( 5, 2) =     3.163140d0 !       Boron -       Helium
      xfac( 5, 2) =     1.974170d0 !       Boron -       Helium
      alpb( 5, 3) =     3.000118d0 !       Boron -      Lithium
      xfac( 5, 3) =     5.549420d0 !       Boron -      Lithium
      alpb( 5, 4) =     1.991820d0 !       Boron -    Beryllium
      xfac( 5, 4) =     1.213171d0 !       Boron -    Beryllium
      alpb( 5, 5) =     2.181999d0 !       Boron -        Boron
      xfac( 5, 5) =     0.964011d0 !       Boron -        Boron
 !
      alpb( 6, 1) =     1.038716d0 !      Carbon -     Hydrogen
      xfac( 6, 1) =     0.204582d0 !      Carbon -     Hydrogen
      alpb( 6, 2) =     3.042705d0 !      Carbon -       Helium
      xfac( 6, 2) =     3.213971d0 !      Carbon -       Helium
      alpb( 6, 3) =     3.061752d0 !      Carbon -      Lithium
      xfac( 6, 3) =     7.351225d0 !      Carbon -      Lithium
      alpb( 6, 4) =     2.798358d0 !      Carbon -    Beryllium
      xfac( 6, 4) =     2.803943d0 !      Carbon -    Beryllium
      alpb( 6, 5) =     2.650092d0 !      Carbon -        Boron
      xfac( 6, 5) =     1.456708d0 !      Carbon -        Boron
      alpb( 6, 6) =     2.655746d0 !      Carbon -       Carbon
      xfac( 6, 6) =     0.937818d0 !      Carbon -       Carbon
 !
      alpb( 7, 1) =     1.049958d0 !    Nitrogen -     Hydrogen
      xfac( 7, 1) =     0.173966d0 !    Nitrogen -     Hydrogen
      alpb( 7, 2) =     2.814339d0 !    Nitrogen -       Helium
      xfac( 7, 2) =     1.077861d0 !    Nitrogen -       Helium
      alpb( 7, 3) =     2.205062d0 !    Nitrogen -      Lithium
      xfac( 7, 3) =     1.070477d0 !    Nitrogen -      Lithium
      alpb( 7, 4) =     2.391056d0 !    Nitrogen -    Beryllium
      xfac( 7, 4) =     1.588438d0 !    Nitrogen -    Beryllium
      alpb( 7, 5) =     2.264882d0 !    Nitrogen -        Boron
      xfac( 7, 5) =     0.794623d0 !    Nitrogen -        Boron
      alpb( 7, 6) =     2.734636d0 !    Nitrogen -       Carbon
      xfac( 7, 6) =     0.959163d0 !    Nitrogen -       Carbon
      alpb( 7, 7) =     2.786676d0 !    Nitrogen -     Nitrogen
      xfac( 7, 7) =     0.863430d0 !    Nitrogen -     Nitrogen
 !
      alpb( 8, 1) =     1.508030d0 !      Oxygen -     Hydrogen
      xfac( 8, 1) =     0.160602d0 !      Oxygen -     Hydrogen
      alpb( 8, 2) =     3.595772d0 !      Oxygen -       Helium
      xfac( 8, 2) =     6.688705d0 !      Oxygen -       Helium
      alpb( 8, 3) =     2.090902d0 !      Oxygen -      Lithium
      xfac( 8, 3) =     0.499065d0 !      Oxygen -      Lithium
      alpb( 8, 4) =     3.081366d0 !      Oxygen -    Beryllium
      xfac( 8, 4) =     2.735812d0 !      Oxygen -    Beryllium
      alpb( 8, 5) =     2.694696d0 !      Oxygen -        Boron
      xfac( 8, 5) =     1.106827d0 !      Oxygen -        Boron
      alpb( 8, 6) =     2.850931d0 !      Oxygen -       Carbon
      xfac( 8, 6) =     0.826848d0 !      Oxygen -       Carbon
      alpb( 8, 7) =     2.841313d0 !      Oxygen -     Nitrogen
      xfac( 8, 7) =     0.699190d0 !      Oxygen -     Nitrogen
      alpb( 8, 8) =     2.457135d0 !      Oxygen -       Oxygen
      xfac( 8, 8) =     0.356439d0 !      Oxygen -       Oxygen
 !
      alpb( 9, 1) =     3.104012d0 !    Fluorine -     Hydrogen
      xfac( 9, 1) =     0.587775d0 !    Fluorine -     Hydrogen
      alpb( 9, 2) =     2.856543d0 !    Fluorine -       Helium
      xfac( 9, 2) =     0.745107d0 !    Fluorine -       Helium
      alpb( 9, 3) =     3.036354d0 !    Fluorine -      Lithium
      xfac( 9, 3) =     0.756703d0 !    Fluorine -      Lithium
      alpb( 9, 4) =     3.328175d0 !    Fluorine -    Beryllium
      xfac( 9, 4) =     2.273878d0 !    Fluorine -    Beryllium
      alpb( 9, 5) =     2.841215d0 !    Fluorine -        Boron
      xfac( 9, 5) =     0.892467d0 !    Fluorine -        Boron
      alpb( 9, 6) =     3.230134d0 !    Fluorine -       Carbon
      xfac( 9, 6) =     0.946026d0 !    Fluorine -       Carbon
      alpb( 9, 7) =     3.207271d0 !    Fluorine -     Nitrogen
      xfac( 9, 7) =     0.894011d0 !    Fluorine -     Nitrogen
      alpb( 9, 8) =     3.140885d0 !    Fluorine -       Oxygen
      xfac( 9, 8) =     0.698033d0 !    Fluorine -       Oxygen
      alpb( 9, 9) =     4.492940d0 !    Fluorine -     Fluorine
      xfac( 9, 9) =     3.111004d0 !    Fluorine -     Fluorine
 !
      alpb(10, 1) =     5.999680d0 !        Neon -     Hydrogen
      xfac(10, 1) =     5.535021d0 !        Neon -     Hydrogen
      alpb(10, 2) =     3.677758d0 !        Neon -       Helium
      xfac(10, 2) =     1.960924d0 !        Neon -       Helium
      alpb(10, 3) =     2.242969d0 !        Neon -      Lithium
      xfac(10, 3) =     0.642933d0 !        Neon -      Lithium
      alpb(10, 4) =     0.832530d0 !        Neon -    Beryllium
      xfac(10, 4) =     0.140208d0 !        Neon -    Beryllium
      alpb(10, 5) =     2.756190d0 !        Neon -        Boron
      xfac(10, 5) =     2.764140d0 !        Neon -        Boron
      alpb(10, 6) =     3.441188d0 !        Neon -       Carbon
      xfac(10, 6) =     5.468780d0 !        Neon -       Carbon
      alpb(10, 7) =     4.426370d0 !        Neon -     Nitrogen
      xfac(10, 7) =    29.999609d0 !        Neon -     Nitrogen
      alpb(10, 8) =     2.906840d0 !        Neon -       Oxygen
      xfac(10, 8) =     0.753518d0 !        Neon -       Oxygen
      alpb(10, 9) =     3.675611d0 !        Neon -     Fluorine
      xfac(10, 9) =     2.706754d0 !        Neon -     Fluorine
      alpb(10,10) =     5.180440d0 !        Neon -         Neon
      xfac(10,10) =     0.500000d0 !        Neon -         Neon
 !
      alpb(11, 1) =     1.619287d0 !      Sodium -     Hydrogen
      xfac(11, 1) =     1.966963d0 !      Sodium -     Hydrogen
      alpb(11, 2) =     2.171840d0 !      Sodium -       Helium
      xfac(11, 2) =     4.590369d0 !      Sodium -       Helium
      alpb(11, 3) =     0.898897d0 !      Sodium -      Lithium
      xfac(11, 3) =     0.655446d0 !      Sodium -      Lithium
      alpb(11, 4) =     1.255480d0 !      Sodium -    Beryllium
      xfac(11, 4) =     3.121620d0 !      Sodium -    Beryllium
      alpb(11, 5) =     2.476698d0 !      Sodium -        Boron
      xfac(11, 5) =     8.164739d0 !      Sodium -        Boron
      alpb(11, 6) =     2.394648d0 !      Sodium -       Carbon
      xfac(11, 6) =     6.318544d0 !      Sodium -       Carbon
      alpb(11, 7) =     2.482865d0 !      Sodium -     Nitrogen
      xfac(11, 7) =     5.755473d0 !      Sodium -     Nitrogen
      alpb(11, 8) =     2.699675d0 !      Sodium -       Oxygen
      xfac(11, 8) =     8.556302d0 !      Sodium -       Oxygen
      alpb(11, 9) =     3.036873d0 !      Sodium -     Fluorine
      xfac(11, 9) =     9.250233d0 !      Sodium -     Fluorine
      alpb(11,10) =     1.469568d0 !      Sodium -         Neon
      xfac(11,10) =     0.697745d0 !      Sodium -         Neon
      alpb(11,11) =     1.994455d0 !      Sodium -       Sodium
      xfac(11,11) =     9.335783d0 !      Sodium -       Sodium
 !
      alpb(12, 1) =     2.423259d0 !   Magnesium -     Hydrogen
      xfac(12, 1) =     6.170068d0 !   Magnesium -     Hydrogen
      alpb(12, 2) =     2.289485d0 !   Magnesium -       Helium
      xfac(12, 2) =     3.779366d0 !   Magnesium -       Helium
      alpb(12, 3) =     1.374791d0 !   Magnesium -      Lithium
      xfac(12, 3) =     2.510632d0 !   Magnesium -      Lithium
      alpb(12, 4) =     1.593445d0 !   Magnesium -    Beryllium
      xfac(12, 4) =     2.960809d0 !   Magnesium -    Beryllium
      alpb(12, 5) =     2.466919d0 !   Magnesium -        Boron
      xfac(12, 5) =     6.072802d0 !   Magnesium -        Boron
      alpb(12, 6) =     2.321772d0 !   Magnesium -       Carbon
      xfac(12, 6) =     3.390341d0 !   Magnesium -       Carbon
      alpb(12, 7) =     2.025732d0 !   Magnesium -     Nitrogen
      xfac(12, 7) =     2.115961d0 !   Magnesium -     Nitrogen
      alpb(12, 8) =     2.730174d0 !   Magnesium -       Oxygen
      xfac(12, 8) =     2.888295d0 !   Magnesium -       Oxygen
      alpb(12, 9) =     3.378507d0 !   Magnesium -     Fluorine
      xfac(12, 9) =     5.439497d0 !   Magnesium -     Fluorine
      alpb(12,10) =     0.922342d0 !   Magnesium -         Neon
      xfac(12,10) =     0.452636d0 !   Magnesium -         Neon
      alpb(12,11) =     1.682212d0 !   Magnesium -       Sodium
      xfac(12,11) =     8.429332d0 !   Magnesium -       Sodium
      alpb(12,12) =     2.055257d0 !   Magnesium -    Magnesium
      xfac(12,12) =    20.557591d0 !   Magnesium -    Magnesium
 !
      alpb(13, 1) =     1.610842d0 !    Aluminum -     Hydrogen
      xfac(13, 1) =     1.183718d0 !    Aluminum -     Hydrogen
      alpb(13, 2) =     2.255830d0 !    Aluminum -       Helium
      xfac(13, 2) =     2.701400d0 !    Aluminum -       Helium
      alpb(13, 3) =     1.248327d0 !    Aluminum -      Lithium
      xfac(13, 3) =     0.929842d0 !    Aluminum -      Lithium
      alpb(13, 4) =     1.916502d0 !    Aluminum -    Beryllium
      xfac(13, 4) =     4.229824d0 !    Aluminum -    Beryllium
      alpb(13, 5) =     1.990198d0 !    Aluminum -        Boron
      xfac(13, 5) =     2.676137d0 !    Aluminum -        Boron
      alpb(13, 6) =     2.058949d0 !    Aluminum -       Carbon
      xfac(13, 6) =     3.161882d0 !    Aluminum -       Carbon
      alpb(13, 7) =     1.477524d0 !    Aluminum -     Nitrogen
      xfac(13, 7) =     0.883919d0 !    Aluminum -     Nitrogen
      alpb(13, 8) =     2.054038d0 !    Aluminum -       Oxygen
      xfac(13, 8) =     1.619036d0 !    Aluminum -       Oxygen
      alpb(13, 9) =     2.253927d0 !    Aluminum -     Fluorine
      xfac(13, 9) =     1.368035d0 !    Aluminum -     Fluorine
      alpb(13,10) =     2.528574d0 !    Aluminum -         Neon
      xfac(13,10) =     1.702157d0 !    Aluminum -         Neon
      alpb(13,11) =     1.141388d0 !    Aluminum -       Sodium
      xfac(13,11) =     1.128163d0 !    Aluminum -       Sodium
      alpb(13,12) =     1.455074d0 !    Aluminum -    Magnesium
      xfac(13,12) =     1.821180d0 !    Aluminum -    Magnesium
      alpb(13,13) =     1.224852d0 !    Aluminum -     Aluminum
      xfac(13,13) =     1.669052d0 !    Aluminum -     Aluminum
 !
      alpb(14, 1) =     1.542308d0 !     Silicon -     Hydrogen
      xfac(14, 1) =     0.688945d0 !     Silicon -     Hydrogen
      alpb(14, 2) =     2.028628d0 !     Silicon -       Helium
      xfac(14, 2) =     1.976149d0 !     Silicon -       Helium
      alpb(14, 3) =     1.911808d0 !     Silicon -      Lithium
      xfac(14, 3) =     2.989391d0 !     Silicon -      Lithium
      alpb(14, 4) =     2.162457d0 !     Silicon -    Beryllium
      xfac(14, 4) =     4.322374d0 !     Silicon -    Beryllium
      alpb(14, 5) =     1.915795d0 !     Silicon -        Boron
      xfac(14, 5) =     1.162577d0 !     Silicon -        Boron
      alpb(14, 6) =     1.673306d0 !     Silicon -       Carbon
      xfac(14, 6) =     0.501779d0 !     Silicon -       Carbon
      alpb(14, 7) =     1.854197d0 !     Silicon -     Nitrogen
      xfac(14, 7) =     0.671576d0 !     Silicon -     Nitrogen
      alpb(14, 8) =     1.824047d0 !     Silicon -       Oxygen
      xfac(14, 8) =     0.502254d0 !     Silicon -       Oxygen
      alpb(14, 9) =     2.160762d0 !     Silicon -     Fluorine
      xfac(14, 9) =     0.564372d0 !     Silicon -     Fluorine
      alpb(14,10) =     2.655346d0 !     Silicon -         Neon
      xfac(14,10) =    12.754805d0 !     Silicon -         Neon
      alpb(14,11) =     1.842304d0 !     Silicon -       Sodium
      xfac(14,11) =     9.125996d0 !     Silicon -       Sodium
      alpb(14,12) =     1.157990d0 !     Silicon -    Magnesium
      xfac(14,12) =     0.527802d0 !     Silicon -    Magnesium
      alpb(14,13) =     1.300963d0 !     Silicon -     Aluminum
      xfac(14,13) =     1.290056d0 !     Silicon -     Aluminum
      alpb(14,14) =     1.109923d0 !     Silicon -      Silicon
      xfac(14,14) =     0.369696d0 !     Silicon -      Silicon
 !
      alpb(15, 1) =     1.694074d0 !  Phosphorus -     Hydrogen
      xfac(15, 1) =     1.829199d0 !  Phosphorus -     Hydrogen
      alpb(15, 2) =     2.076667d0 !  Phosphorus -       Helium
      xfac(15, 2) =     1.493985d0 !  Phosphorus -       Helium
      alpb(15, 3) =     1.727121d0 !  Phosphorus -      Lithium
      xfac(15, 3) =     5.864987d0 !  Phosphorus -      Lithium
      alpb(15, 4) =     1.872176d0 !  Phosphorus -    Beryllium
      xfac(15, 4) =     2.310001d0 !  Phosphorus -    Beryllium
      alpb(15, 5) =     1.742693d0 !  Phosphorus -        Boron
      xfac(15, 5) =     2.541187d0 !  Phosphorus -        Boron
      alpb(15, 6) =     1.770260d0 !  Phosphorus -       Carbon
      xfac(15, 6) =     1.381969d0 !  Phosphorus -       Carbon
      alpb(15, 7) =     1.958580d0 !  Phosphorus -     Nitrogen
      xfac(15, 7) =     1.312039d0 !  Phosphorus -     Nitrogen
      alpb(15, 8) =     2.293080d0 !  Phosphorus -       Oxygen
      xfac(15, 8) =     1.314436d0 !  Phosphorus -       Oxygen
      alpb(15, 9) =     2.513083d0 !  Phosphorus -     Fluorine
      xfac(15, 9) =     1.086188d0 !  Phosphorus -     Fluorine
      alpb(15,10) =     2.243688d0 !  Phosphorus -         Neon
      xfac(15,10) =     0.762937d0 !  Phosphorus -         Neon
      alpb(15,11) =     1.518961d0 !  Phosphorus -       Sodium
      xfac(15,11) =     3.713750d0 !  Phosphorus -       Sodium
      alpb(15,12) =     1.297069d0 !  Phosphorus -    Magnesium
      xfac(15,12) =     1.585367d0 !  Phosphorus -    Magnesium
      alpb(15,13) =     1.375504d0 !  Phosphorus -     Aluminum
      xfac(15,13) =     3.249399d0 !  Phosphorus -     Aluminum
      alpb(15,14) =     0.895674d0 !  Phosphorus -      Silicon
      xfac(15,14) =     0.616954d0 !  Phosphorus -      Silicon
      alpb(15,15) =     1.329209d0 !  Phosphorus -   Phosphorus
      xfac(15,15) =     2.558579d0 !  Phosphorus -   Phosphorus
 !
      alpb(16, 1) =     2.182464d0 !      Sulfur -     Hydrogen
      xfac(16, 1) =     0.703252d0 !      Sulfur -     Hydrogen
      alpb(16, 2) =     1.959149d0 !      Sulfur -       Helium
      xfac(16, 2) =     0.437618d0 !      Sulfur -       Helium
      alpb(16, 3) =     1.737806d0 !      Sulfur -      Lithium
      xfac(16, 3) =     0.566769d0 !      Sulfur -      Lithium
      alpb(16, 4) =     2.575836d0 !      Sulfur -    Beryllium
      xfac(16, 4) =     3.179465d0 !      Sulfur -    Beryllium
      alpb(16, 5) =     2.363313d0 !      Sulfur -        Boron
      xfac(16, 5) =     1.177082d0 !      Sulfur -        Boron
      alpb(16, 6) =     2.429136d0 !      Sulfur -       Carbon
      xfac(16, 6) =     0.843145d0 !      Sulfur -       Carbon
      alpb(16, 7) =     2.653791d0 !      Sulfur -     Nitrogen
      xfac(16, 7) =     1.197307d0 !      Sulfur -     Nitrogen
      alpb(16, 8) =     2.508022d0 !      Sulfur -       Oxygen
      xfac(16, 8) =     0.729340d0 !      Sulfur -       Oxygen
      alpb(16, 9) =     2.533157d0 !      Sulfur -     Fluorine
      xfac(16, 9) =     0.534080d0 !      Sulfur -     Fluorine
      alpb(16,10) =     2.787058d0 !      Sulfur -         Neon
      xfac(16,10) =     3.296160d0 !      Sulfur -         Neon
      alpb(16,11) =     2.614090d0 !      Sulfur -       Sodium
      xfac(16,11) =     6.263298d0 !      Sulfur -       Sodium
      alpb(16,12) =     1.442313d0 !      Sulfur -    Magnesium
      xfac(16,12) =     0.578881d0 !      Sulfur -    Magnesium
      alpb(16,13) =     1.706655d0 !      Sulfur -     Aluminum
      xfac(16,13) =     1.677290d0 !      Sulfur -     Aluminum
      alpb(16,14) =     1.647931d0 !      Sulfur -      Silicon
      xfac(16,14) =     0.553963d0 !      Sulfur -      Silicon
      alpb(16,15) =     1.596824d0 !      Sulfur -   Phosphorus
      xfac(16,15) =     1.189185d0 !      Sulfur -   Phosphorus
      alpb(16,16) =     1.985120d0 !      Sulfur -       Sulfur
      xfac(16,16) =     0.509363d0 !      Sulfur -       Sulfur
 !
      alpb(17, 1) =     2.548456d0 !    Chlorine -     Hydrogen
      xfac(17, 1) =     0.721003d0 !    Chlorine -     Hydrogen
      alpb(17, 2) =     1.671634d0 !    Chlorine -       Helium
      xfac(17, 2) =     0.500002d0 !    Chlorine -       Helium
      alpb(17, 3) =     2.838217d0 !    Chlorine -      Lithium
      xfac(17, 3) =     2.531354d0 !    Chlorine -      Lithium
      alpb(17, 4) =     2.716560d0 !    Chlorine -    Beryllium
      xfac(17, 4) =     2.638266d0 !    Chlorine -    Beryllium
      alpb(17, 5) =     2.228737d0 !    Chlorine -        Boron
      xfac(17, 5) =     0.742613d0 !    Chlorine -        Boron
      alpb(17, 6) =     2.450283d0 !    Chlorine -       Carbon
      xfac(17, 6) =     0.690352d0 !    Chlorine -       Carbon
      alpb(17, 7) =     2.385624d0 !    Chlorine -     Nitrogen
      xfac(17, 7) =     0.659831d0 !    Chlorine -     Nitrogen
      alpb(17, 8) =     2.288635d0 !    Chlorine -       Oxygen
      xfac(17, 8) =     0.437776d0 !    Chlorine -       Oxygen
      alpb(17, 9) =     2.575402d0 !    Chlorine -     Fluorine
      xfac(17, 9) =     0.503186d0 !    Chlorine -     Fluorine
      alpb(17,10) =     1.732740d0 !    Chlorine -         Neon
      xfac(17,10) =     0.499482d0 !    Chlorine -         Neon
      alpb(17,11) =     2.536945d0 !    Chlorine -       Sodium
      xfac(17,11) =    10.364642d0 !    Chlorine -       Sodium
      alpb(17,12) =     2.292455d0 !    Chlorine -    Magnesium
      xfac(17,12) =     2.207847d0 !    Chlorine -    Magnesium
      alpb(17,13) =     1.678498d0 !    Chlorine -     Aluminum
      xfac(17,13) =     1.079875d0 !    Chlorine -     Aluminum
      alpb(17,14) =     1.818389d0 !    Chlorine -      Silicon
      xfac(17,14) =     0.590060d0 !    Chlorine -      Silicon
      alpb(17,15) =     1.297513d0 !    Chlorine -   Phosphorus
      xfac(17,15) =     0.496813d0 !    Chlorine -   Phosphorus
      alpb(17,16) =     2.167945d0 !    Chlorine -       Sulfur
      xfac(17,16) =     0.624384d0 !    Chlorine -       Sulfur
      alpb(17,17) =     2.469990d0 !    Chlorine -     Chlorine
      xfac(17,17) =     0.900377d0 !    Chlorine -     Chlorine
 !
      alpb(18, 1) =     4.056167d0 !       Argon -     Hydrogen
      xfac(18, 1) =     3.933445d0 !       Argon -     Hydrogen
      alpb(18, 2) =     2.716562d0 !       Argon -       Helium
      xfac(18, 2) =     1.177211d0 !       Argon -       Helium
      alpb(18, 3) =     3.001334d0 !       Argon -      Lithium
      xfac(18, 3) =     2.193788d0 !       Argon -      Lithium
      alpb(18, 4) =     3.227598d0 !       Argon -    Beryllium
      xfac(18, 4) =     2.700296d0 !       Argon -    Beryllium
      alpb(18, 5) =     2.674207d0 !       Argon -        Boron
      xfac(18, 5) =     2.017996d0 !       Argon -        Boron
      alpb(18, 6) =     1.471309d0 !       Argon -       Carbon
      xfac(18, 6) =     0.122309d0 !       Argon -       Carbon
      alpb(18, 7) =     2.326805d0 !       Argon -     Nitrogen
      xfac(18, 7) =     0.562581d0 !       Argon -     Nitrogen
      alpb(18, 8) =     2.228209d0 !       Argon -       Oxygen
      xfac(18, 8) =     0.367713d0 !       Argon -       Oxygen
      alpb(18, 9) =     3.920658d0 !       Argon -     Fluorine
      xfac(18, 9) =     9.269715d0 !       Argon -     Fluorine
      alpb(18,10) =     2.963747d0 !       Argon -         Neon
      xfac(18,10) =     1.304697d0 !       Argon -         Neon
      alpb(18,11) =     2.667734d0 !       Argon -       Sodium
      xfac(18,11) =     5.946915d0 !       Argon -       Sodium
      alpb(18,12) =     1.996514d0 !       Argon -    Magnesium
      xfac(18,12) =     2.030224d0 !       Argon -    Magnesium
      alpb(18,13) =     2.716128d0 !       Argon -     Aluminum
      xfac(18,13) =     1.838228d0 !       Argon -     Aluminum
      alpb(18,14) =     1.935869d0 !       Argon -      Silicon
      xfac(18,14) =     1.288907d0 !       Argon -      Silicon
      alpb(18,15) =     3.998905d0 !       Argon -   Phosphorus
      xfac(18,15) =     0.173766d0 !       Argon -   Phosphorus
      alpb(18,16) =     2.049398d0 !       Argon -       Sulfur
      xfac(18,16) =     0.653769d0 !       Argon -       Sulfur
      alpb(18,17) =     2.554449d0 !       Argon -     Chlorine
      xfac(18,17) =     2.256094d0 !       Argon -     Chlorine
      alpb(18,18) =     2.306432d0 !       Argon -        Argon
      xfac(18,18) =     0.972699d0 !       Argon -        Argon
 !
      alpb(19, 1) =     2.304518d0 !   Potassium -     Hydrogen
      xfac(19, 1) =    29.964954d0 !   Potassium -     Hydrogen
      alpb(19, 2) =     2.140614d0 !   Potassium -       Helium
      xfac(19, 2) =     6.673621d0 !   Potassium -       Helium
      alpb(19, 3) =     1.108062d0 !   Potassium -      Lithium
      xfac(19, 3) =     4.364297d0 !   Potassium -      Lithium
      alpb(19, 4) =     3.000365d0 !   Potassium -    Beryllium
      xfac(19, 4) =     6.514383d0 !   Potassium -    Beryllium
      alpb(19, 5) =     2.507524d0 !   Potassium -        Boron
      xfac(19, 5) =    28.190857d0 !   Potassium -        Boron
      alpb(19, 6) =     1.769643d0 !   Potassium -       Carbon
      xfac(19, 6) =     2.489951d0 !   Potassium -       Carbon
      alpb(19, 7) =     1.907394d0 !   Potassium -     Nitrogen
      xfac(19, 7) =     3.943077d0 !   Potassium -     Nitrogen
      alpb(19, 8) =     2.151119d0 !   Potassium -       Oxygen
      xfac(19, 8) =     4.281570d0 !   Potassium -       Oxygen
      alpb(19, 9) =     3.065393d0 !   Potassium -     Fluorine
      xfac(19, 9) =    17.321092d0 !   Potassium -     Fluorine
      alpb(19,10) =     1.653125d0 !   Potassium -         Neon
      xfac(19,10) =     1.093188d0 !   Potassium -         Neon
      alpb(19,11) =     0.944935d0 !   Potassium -       Sodium
      xfac(19,11) =     6.450008d0 !   Potassium -       Sodium
      alpb(19,12) =     1.272102d0 !   Potassium -    Magnesium
      xfac(19,12) =     2.832505d0 !   Potassium -    Magnesium
      alpb(19,13) =     1.849469d0 !   Potassium -     Aluminum
      xfac(19,13) =    27.774025d0 !   Potassium -     Aluminum
      alpb(19,14) =     1.674691d0 !   Potassium -      Silicon
      xfac(19,14) =     8.047633d0 !   Potassium -      Silicon
      alpb(19,15) =     1.415563d0 !   Potassium -   Phosphorus
      xfac(19,15) =     4.258021d0 !   Potassium -   Phosphorus
      alpb(19,16) =     2.428403d0 !   Potassium -       Sulfur
      xfac(19,16) =    30.000181d0 !   Potassium -       Sulfur
      alpb(19,17) =     2.346443d0 !   Potassium -     Chlorine
      xfac(19,17) =    12.630753d0 !   Potassium -     Chlorine
      alpb(19,18) =     2.436124d0 !   Potassium -        Argon
      xfac(19,18) =     8.318024d0 !   Potassium -        Argon
      alpb(19,19) =     1.492751d0 !   Potassium -    Potassium
      xfac(19,19) =     6.173527d0 !   Potassium -    Potassium
 !
      alpb(20, 1) =     1.997037d0 !     Calcium -     Hydrogen
      xfac(20, 1) =     5.125996d0 !     Calcium -     Hydrogen
      alpb(20, 2) =     2.150217d0 !     Calcium -       Helium
      xfac(20, 2) =     5.381385d0 !     Calcium -       Helium
      alpb(20, 5) =     1.700010d0 !     Calcium -        Boron
      xfac(20, 5) =     1.700010d0 !     Calcium -        Boron
      alpb(20, 6) =     3.376881d0 !     Calcium -       Carbon
      xfac(20, 6) =    45.518388d0 !     Calcium -       Carbon
      alpb(20, 7) =     2.335548d0 !     Calcium -     Nitrogen
      xfac(20, 7) =     3.063067d0 !     Calcium -     Nitrogen
      alpb(20, 8) =     3.347983d0 !     Calcium -       Oxygen
      xfac(20, 8) =     8.353090d0 !     Calcium -       Oxygen
      alpb(20, 9) =     3.871263d0 !     Calcium -     Fluorine
      xfac(20, 9) =    14.692101d0 !     Calcium -     Fluorine
      alpb(20,10) =     1.247453d0 !     Calcium -         Neon
      xfac(20,10) =     0.493997d0 !     Calcium -         Neon
      alpb(20,11) =     2.172223d0 !     Calcium -       Sodium
      xfac(20,11) =    10.049083d0 !     Calcium -       Sodium
      alpb(20,12) =     1.612133d0 !     Calcium -    Magnesium
      xfac(20,12) =     5.062878d0 !     Calcium -    Magnesium
      alpb(20,13) =     1.612565d0 !     Calcium -     Aluminum
      xfac(20,13) =     4.188555d0 !     Calcium -     Aluminum
      alpb(20,14) =     1.730018d0 !     Calcium -      Silicon
      xfac(20,14) =     4.282139d0 !     Calcium -      Silicon
      alpb(20,15) =     1.922605d0 !     Calcium -   Phosphorus
      xfac(20,15) =    15.033250d0 !     Calcium -   Phosphorus
      alpb(20,16) =     1.481189d0 !     Calcium -       Sulfur
      xfac(20,16) =     0.561550d0 !     Calcium -       Sulfur
      alpb(20,17) =     2.785624d0 !     Calcium -     Chlorine
      xfac(20,17) =     8.996518d0 !     Calcium -     Chlorine
      alpb(20,18) =     1.544903d0 !     Calcium -        Argon
      xfac(20,18) =     0.699868d0 !     Calcium -        Argon
      alpb(20,19) =     1.210391d0 !     Calcium -    Potassium
      xfac(20,19) =     1.755307d0 !     Calcium -    Potassium
      alpb(20,20) =     1.477787d0 !     Calcium -      Calcium
      xfac(20,20) =     5.134189d0 !     Calcium -      Calcium
 !
      alpb(21, 1) =     2.630734d0 !    Scandium -     Hydrogen
      xfac(21, 1) =     5.354101d0 !    Scandium -     Hydrogen
      alpb(21, 6) =     2.774943d0 !    Scandium -       Carbon
      xfac(21, 6) =    13.452840d0 !    Scandium -       Carbon
      alpb(21, 7) =     2.081124d0 !    Scandium -     Nitrogen
      xfac(21, 7) =     1.980291d0 !    Scandium -     Nitrogen
      alpb(21, 8) =     2.238586d0 !    Scandium -       Oxygen
      xfac(21, 8) =     1.567669d0 !    Scandium -       Oxygen
      alpb(21, 9) =     3.226175d0 !    Scandium -     Fluorine
      xfac(21, 9) =     7.919620d0 !    Scandium -     Fluorine
      alpb(21,13) =     1.003550d0 !    Scandium -     Aluminum
      xfac(21,13) =     0.500620d0 !    Scandium -     Aluminum
      alpb(21,14) =     1.849600d0 !    Scandium -      Silicon
      xfac(21,14) =     2.767826d0 !    Scandium -      Silicon
      alpb(21,15) =     1.919608d0 !    Scandium -   Phosphorus
      xfac(21,15) =     4.663061d0 !    Scandium -   Phosphorus
      alpb(21,16) =     1.111949d0 !    Scandium -       Sulfur
      xfac(21,16) =     0.498540d0 !    Scandium -       Sulfur
      alpb(21,17) =     2.094163d0 !    Scandium -     Chlorine
      xfac(21,17) =     2.355302d0 !    Scandium -     Chlorine
      alpb(21,21) =     2.106571d0 !    Scandium -     Scandium
      xfac(21,21) =    30.002441d0 !    Scandium -     Scandium
 !
      alpb(22, 1) =     1.447725d0 !    Titanium -     Hydrogen
      xfac(22, 1) =     0.603333d0 !    Titanium -     Hydrogen
      alpb(22, 3) =     1.514050d0 !    Titanium -      Lithium
      xfac(22, 3) =     0.502488d0 !    Titanium -      Lithium
      alpb(22, 5) =     1.628710d0 !    Titanium -        Boron
      xfac(22, 5) =     0.649360d0 !    Titanium -        Boron
      alpb(22, 6) =     1.798067d0 !    Titanium -       Carbon
      xfac(22, 6) =     0.562296d0 !    Titanium -       Carbon
      alpb(22, 7) =     1.638936d0 !    Titanium -     Nitrogen
      xfac(22, 7) =     0.543706d0 !    Titanium -     Nitrogen
      alpb(22, 8) =     1.962314d0 !    Titanium -       Oxygen
      xfac(22, 8) =     0.872204d0 !    Titanium -       Oxygen
      alpb(22, 9) =     2.186657d0 !    Titanium -     Fluorine
      xfac(22, 9) =     0.836131d0 !    Titanium -     Fluorine
      alpb(22,11) =     1.124786d0 !    Titanium -       Sodium
      xfac(22,11) =     1.987793d0 !    Titanium -       Sodium
      alpb(22,12) =     1.900606d0 !    Titanium -    Magnesium
      xfac(22,12) =     6.889073d0 !    Titanium -    Magnesium
      alpb(22,13) =     1.833384d0 !    Titanium -     Aluminum
      xfac(22,13) =     8.952566d0 !    Titanium -     Aluminum
      alpb(22,14) =     1.373954d0 !    Titanium -      Silicon
      xfac(22,14) =     0.561089d0 !    Titanium -      Silicon
      alpb(22,15) =     1.610003d0 !    Titanium -   Phosphorus
      xfac(22,15) =     3.074680d0 !    Titanium -   Phosphorus
      alpb(22,16) =     2.309450d0 !    Titanium -       Sulfur
      xfac(22,16) =     1.781817d0 !    Titanium -       Sulfur
      alpb(22,17) =     1.953656d0 !    Titanium -     Chlorine
      xfac(22,17) =     0.831301d0 !    Titanium -     Chlorine
      alpb(22,20) =     1.268314d0 !    Titanium -      Calcium
      xfac(22,20) =     0.513504d0 !    Titanium -      Calcium
      alpb(22,22) =     2.445684d0 !    Titanium -     Titanium
      xfac(22,22) =    29.795082d0 !    Titanium -     Titanium
 !
      alpb(23, 1) =     1.454900d0 !    Vanadium -     Hydrogen
      xfac(23, 1) =     0.350807d0 !    Vanadium -     Hydrogen
      alpb(23, 6) =     1.904429d0 !    Vanadium -       Carbon
      xfac(23, 6) =     0.489034d0 !    Vanadium -       Carbon
      alpb(23, 7) =     2.139547d0 !    Vanadium -     Nitrogen
      xfac(23, 7) =     0.964593d0 !    Vanadium -     Nitrogen
      alpb(23, 8) =     2.076717d0 !    Vanadium -       Oxygen
      xfac(23, 8) =     0.789091d0 !    Vanadium -       Oxygen
      alpb(23, 9) =     2.483525d0 !    Vanadium -     Fluorine
      xfac(23, 9) =     1.056377d0 !    Vanadium -     Fluorine
      alpb(23,11) =     2.548904d0 !    Vanadium -       Sodium
      xfac(23,11) =     8.346697d0 !    Vanadium -       Sodium
      alpb(23,15) =     2.205190d0 !    Vanadium -   Phosphorus
      xfac(23,15) =     6.763663d0 !    Vanadium -   Phosphorus
      alpb(23,16) =     2.407934d0 !    Vanadium -       Sulfur
      xfac(23,16) =     1.374332d0 !    Vanadium -       Sulfur
      alpb(23,17) =     2.395745d0 !    Vanadium -     Chlorine
      xfac(23,17) =     1.590959d0 !    Vanadium -     Chlorine
      alpb(23,19) =     1.361275d0 !    Vanadium -    Potassium
      xfac(23,19) =     1.893631d0 !    Vanadium -    Potassium
      alpb(23,23) =     1.859935d0 !    Vanadium -     Vanadium
      xfac(23,23) =     0.953942d0 !    Vanadium -     Vanadium
 !
      alpb(24, 1) =     1.710489d0 !    Chromium -     Hydrogen
      xfac(24, 1) =     0.451845d0 !    Chromium -     Hydrogen
      alpb(24, 3) =     1.554282d0 !    Chromium -      Lithium
      xfac(24, 3) =     1.523425d0 !    Chromium -      Lithium
      alpb(24, 6) =     2.200250d0 !    Chromium -       Carbon
      xfac(24, 6) =     0.723497d0 !    Chromium -       Carbon
      alpb(24, 7) =     1.978476d0 !    Chromium -     Nitrogen
      xfac(24, 7) =     0.431966d0 !    Chromium -     Nitrogen
      alpb(24, 8) =     2.226688d0 !    Chromium -       Oxygen
      xfac(24, 8) =     0.603066d0 !    Chromium -       Oxygen
      alpb(24, 9) =     2.545695d0 !    Chromium -     Fluorine
      xfac(24, 9) =     0.581501d0 !    Chromium -     Fluorine
      alpb(24,11) =     1.742438d0 !    Chromium -       Sodium
      xfac(24,11) =     7.141413d0 !    Chromium -       Sodium
      alpb(24,12) =     1.949255d0 !    Chromium -    Magnesium
      xfac(24,12) =     9.004042d0 !    Chromium -    Magnesium
      alpb(24,14) =     1.632536d0 !    Chromium -      Silicon
      xfac(24,14) =     1.831750d0 !    Chromium -      Silicon
      alpb(24,15) =     0.965663d0 !    Chromium -   Phosphorus
      xfac(24,15) =     0.488071d0 !    Chromium -   Phosphorus
      alpb(24,16) =     2.022399d0 !    Chromium -       Sulfur
      xfac(24,16) =     0.610052d0 !    Chromium -       Sulfur
      alpb(24,17) =     2.494604d0 !    Chromium -     Chlorine
      xfac(24,17) =     0.987014d0 !    Chromium -     Chlorine
      alpb(24,19) =     1.827441d0 !    Chromium -    Potassium
      xfac(24,19) =    14.122878d0 !    Chromium -    Potassium
      alpb(24,20) =     1.748419d0 !    Chromium -      Calcium
      xfac(24,20) =     3.971766d0 !    Chromium -      Calcium
      alpb(24,24) =     2.859778d0 !    Chromium -     Chromium
      xfac(24,24) =    21.294482d0 !    Chromium -     Chromium
 !
      alpb(25, 1) =     1.815287d0 !   Manganese -     Hydrogen
      xfac(25, 1) =     1.334984d0 !   Manganese -     Hydrogen
      alpb(25, 6) =     2.122570d0 !   Manganese -       Carbon
      xfac(25, 6) =     1.646822d0 !   Manganese -       Carbon
      alpb(25, 7) =     2.625097d0 !   Manganese -     Nitrogen
      xfac(25, 7) =     2.366982d0 !   Manganese -     Nitrogen
      alpb(25, 8) =     3.225970d0 !   Manganese -       Oxygen
      xfac(25, 8) =     3.636943d0 !   Manganese -       Oxygen
      alpb(25, 9) =     3.508953d0 !   Manganese -     Fluorine
      xfac(25, 9) =     2.404476d0 !   Manganese -     Fluorine
      alpb(25,13) =     1.231200d0 !   Manganese -     Aluminum
      xfac(25,13) =     1.130368d0 !   Manganese -     Aluminum
      alpb(25,14) =     1.881580d0 !   Manganese -      Silicon
      xfac(25,14) =     3.934609d0 !   Manganese -      Silicon
      alpb(25,15) =     1.879268d0 !   Manganese -   Phosphorus
      xfac(25,15) =     5.259289d0 !   Manganese -   Phosphorus
      alpb(25,16) =     2.205580d0 !   Manganese -       Sulfur
      xfac(25,16) =     2.583375d0 !   Manganese -       Sulfur
      alpb(25,17) =     2.275167d0 !   Manganese -     Chlorine
      xfac(25,17) =     2.025304d0 !   Manganese -     Chlorine
      alpb(25,19) =     1.328545d0 !   Manganese -    Potassium
      xfac(25,19) =     1.921563d0 !   Manganese -    Potassium
      alpb(25,20) =     1.298445d0 !   Manganese -      Calcium
      xfac(25,20) =     0.520488d0 !   Manganese -      Calcium
      alpb(25,22) =     1.633575d0 !   Manganese -     Titanium
      xfac(25,22) =     4.212201d0 !   Manganese -     Titanium
      alpb(25,25) =     2.502150d0 !   Manganese -    Manganese
      xfac(25,25) =    23.014869d0 !   Manganese -    Manganese
 !
      alpb(26, 1) =     2.325000d0 !        Iron -     Hydrogen
      xfac(26, 1) =     0.797044d0 !        Iron -     Hydrogen
      alpb(26, 6) =     2.439391d0 !        Iron -       Carbon
      xfac(26, 6) =     0.840113d0 !        Iron -       Carbon
      alpb(26, 7) =     2.710121d0 !        Iron -     Nitrogen
      xfac(26, 7) =     1.307687d0 !        Iron -     Nitrogen
      alpb(26, 8) =     2.977229d0 !        Iron -       Oxygen
      xfac(26, 8) =     1.669098d0 !        Iron -       Oxygen
      alpb(26, 9) =     3.266034d0 !        Iron -     Fluorine
      xfac(26, 9) =     1.572783d0 !        Iron -     Fluorine
      alpb(26,14) =     2.261269d0 !        Iron -      Silicon
      xfac(26,14) =     1.302779d0 !        Iron -      Silicon
      alpb(26,15) =     1.425836d0 !        Iron -   Phosphorus
      xfac(26,15) =     0.597968d0 !        Iron -   Phosphorus
      alpb(26,16) =     2.922342d0 !        Iron -       Sulfur
      xfac(26,16) =     3.055008d0 !        Iron -       Sulfur
      alpb(26,17) =     2.803764d0 !        Iron -     Chlorine
      xfac(26,17) =     1.475990d0 !        Iron -     Chlorine
      alpb(26,19) =     0.914983d0 !        Iron -    Potassium
      xfac(26,19) =     0.471163d0 !        Iron -    Potassium
      alpb(26,22) =     2.152071d0 !        Iron -     Titanium
      xfac(26,22) =     1.718797d0 !        Iron -     Titanium
      alpb(26,24) =     2.320197d0 !        Iron -     Chromium
      xfac(26,24) =     1.605266d0 !        Iron -     Chromium
      alpb(26,26) =     3.253806d0 !        Iron -         Iron
      xfac(26,26) =    25.101048d0 !        Iron -         Iron
 !
      alpb(27, 1) =     2.212022d0 !      Cobalt -     Hydrogen
      xfac(27, 1) =     0.781287d0 !      Cobalt -     Hydrogen
      alpb(27, 3) =     1.930303d0 !      Cobalt -      Lithium
      xfac(27, 3) =     0.523612d0 !      Cobalt -      Lithium
      alpb(27, 5) =     3.200000d0 !      Cobalt -        Boron
      xfac(27, 5) =     1.000000d0 !      Cobalt -        Boron
      alpb(27, 6) =     1.369735d0 !      Cobalt -       Carbon
      xfac(27, 6) =     0.101941d0 !      Cobalt -       Carbon
      alpb(27, 7) =     2.018692d0 !      Cobalt -     Nitrogen
      xfac(27, 7) =     0.371117d0 !      Cobalt -     Nitrogen
      alpb(27, 8) =     2.512985d0 !      Cobalt -       Oxygen
      xfac(27, 8) =     0.617937d0 !      Cobalt -       Oxygen
      alpb(27, 9) =     3.169014d0 !      Cobalt -     Fluorine
      xfac(27, 9) =     1.042929d0 !      Cobalt -     Fluorine
      alpb(27,11) =     1.130004d0 !      Cobalt -       Sodium
      xfac(27,11) =     0.525429d0 !      Cobalt -       Sodium
      alpb(27,14) =     2.247195d0 !      Cobalt -      Silicon
      xfac(27,14) =     1.130253d0 !      Cobalt -      Silicon
      alpb(27,15) =     2.298868d0 !      Cobalt -   Phosphorus
      xfac(27,15) =     3.189088d0 !      Cobalt -   Phosphorus
      alpb(27,16) =     2.144853d0 !      Cobalt -       Sulfur
      xfac(27,16) =     0.522339d0 !      Cobalt -       Sulfur
      alpb(27,17) =     2.604673d0 !      Cobalt -     Chlorine
      xfac(27,17) =     0.979572d0 !      Cobalt -     Chlorine
      alpb(27,19) =     1.347379d0 !      Cobalt -    Potassium
      xfac(27,19) =     1.363649d0 !      Cobalt -    Potassium
      alpb(27,24) =     1.965685d0 !      Cobalt -     Chromium
      xfac(27,24) =     0.907585d0 !      Cobalt -     Chromium
      alpb(27,27) =     1.072023d0 !      Cobalt -       Cobalt
      xfac(27,27) =     0.082968d0 !      Cobalt -       Cobalt
 !
      alpb(28, 1) =     1.921141d0 !      Nickel -     Hydrogen
      xfac(28, 1) =     0.694497d0 !      Nickel -     Hydrogen
      alpb(28, 5) =     2.332207d0 !      Nickel -        Boron
      xfac(28, 5) =     0.529685d0 !      Nickel -        Boron
      alpb(28, 6) =     2.135123d0 !      Nickel -       Carbon
      xfac(28, 6) =     0.429059d0 !      Nickel -       Carbon
      alpb(28, 7) =     2.259589d0 !      Nickel -     Nitrogen
      xfac(28, 7) =     0.403691d0 !      Nickel -     Nitrogen
      alpb(28, 8) =     2.452312d0 !      Nickel -       Oxygen
      xfac(28, 8) =     0.284888d0 !      Nickel -       Oxygen
      alpb(28, 9) =     3.145389d0 !      Nickel -     Fluorine
      xfac(28, 9) =     0.559407d0 !      Nickel -     Fluorine
      alpb(28,14) =     2.260625d0 !      Nickel -      Silicon
      xfac(28,14) =     3.024544d0 !      Nickel -      Silicon
      alpb(28,15) =     1.646184d0 !      Nickel -   Phosphorus
      xfac(28,15) =     0.793563d0 !      Nickel -   Phosphorus
      alpb(28,16) =     2.360866d0 !      Nickel -       Sulfur
      xfac(28,16) =     0.923582d0 !      Nickel -       Sulfur
      alpb(28,17) =     2.771621d0 !      Nickel -     Chlorine
      xfac(28,17) =     1.509842d0 !      Nickel -     Chlorine
      alpb(28,19) =     1.110139d0 !      Nickel -    Potassium
      xfac(28,19) =     0.642360d0 !      Nickel -    Potassium
      alpb(28,24) =     2.774356d0 !      Nickel -     Chromium
      xfac(28,24) =    29.999969d0 !      Nickel -     Chromium
      alpb(28,28) =     1.626235d0 !      Nickel -       Nickel
      xfac(28,28) =     0.339558d0 !      Nickel -       Nickel
 !
      alpb(29, 1) =     2.941555d0 !      Copper -     Hydrogen
      xfac(29, 1) =     1.781622d0 !      Copper -     Hydrogen
      alpb(29, 6) =     3.018944d0 !      Copper -       Carbon
      xfac(29, 6) =     1.413488d0 !      Copper -       Carbon
      alpb(29, 7) =     2.566300d0 !      Copper -     Nitrogen
      xfac(29, 7) =     0.429906d0 !      Copper -     Nitrogen
      alpb(29, 8) =     1.911057d0 !      Copper -       Oxygen
      xfac(29, 8) =     0.098068d0 !      Copper -       Oxygen
      alpb(29, 9) =     3.176529d0 !      Copper -     Fluorine
      xfac(29, 9) =     0.411293d0 !      Copper -     Fluorine
      alpb(29,11) =     1.306695d0 !      Copper -       Sodium
      xfac(29,11) =     0.785487d0 !      Copper -       Sodium
      alpb(29,13) =     2.320517d0 !      Copper -     Aluminum
      xfac(29,13) =    12.995965d0 !      Copper -     Aluminum
      alpb(29,15) =     0.858794d0 !      Copper -   Phosphorus
      xfac(29,15) =     5.035151d0 !      Copper -   Phosphorus
      alpb(29,16) =     2.053844d0 !      Copper -       Sulfur
      xfac(29,16) =     0.296518d0 !      Copper -       Sulfur
      alpb(29,17) =     2.475894d0 !      Copper -     Chlorine
      xfac(29,17) =     0.372668d0 !      Copper -     Chlorine
      alpb(29,19) =     2.087357d0 !      Copper -    Potassium
      xfac(29,19) =     7.795310d0 !      Copper -    Potassium
      alpb(29,29) =     3.103277d0 !      Copper -       Copper
      xfac(29,29) =     3.391704d0 !      Copper -       Copper
 !
      alpb(30, 1) =     1.874800d0 !        Zinc -     Hydrogen
      xfac(30, 1) =     1.696831d0 !        Zinc -     Hydrogen
      alpb(30, 6) =     2.171605d0 !        Zinc -       Carbon
      xfac(30, 6) =     2.386580d0 !        Zinc -       Carbon
      alpb(30, 7) =     1.805998d0 !        Zinc -     Nitrogen
      xfac(30, 7) =     0.900539d0 !        Zinc -     Nitrogen
      alpb(30, 8) =     2.079887d0 !        Zinc -       Oxygen
      xfac(30, 8) =     1.116990d0 !        Zinc -       Oxygen
      alpb(30, 9) =     1.859561d0 !        Zinc -     Fluorine
      xfac(30, 9) =     0.499581d0 !        Zinc -     Fluorine
      alpb(30,11) =     1.588584d0 !        Zinc -       Sodium
      xfac(30,11) =     5.694720d0 !        Zinc -       Sodium
      alpb(30,14) =     1.890360d0 !        Zinc -      Silicon
      xfac(30,14) =     6.865738d0 !        Zinc -      Silicon
      alpb(30,15) =     1.398572d0 !        Zinc -   Phosphorus
      xfac(30,15) =     1.863594d0 !        Zinc -   Phosphorus
      alpb(30,16) =     1.379514d0 !        Zinc -       Sulfur
      xfac(30,16) =     0.533478d0 !        Zinc -       Sulfur
      alpb(30,17) =     1.588143d0 !        Zinc -     Chlorine
      xfac(30,17) =     0.547720d0 !        Zinc -     Chlorine
      alpb(30,20) =     0.974041d0 !        Zinc -      Calcium
      xfac(30,20) =     1.296565d0 !        Zinc -      Calcium
      alpb(30,23) =     1.513777d0 !        Zinc -     Vanadium
      xfac(30,23) =     1.442270d0 !        Zinc -     Vanadium
      alpb(30,24) =     2.071878d0 !        Zinc -     Chromium
      xfac(30,24) =     2.312424d0 !        Zinc -     Chromium
      alpb(30,30) =     1.998115d0 !        Zinc -         Zinc
      xfac(30,30) =    19.124599d0 !        Zinc -         Zinc
 !
      alpb(31, 1) =     2.170771d0 !     Gallium -     Hydrogen
      xfac(31, 1) =     2.091955d0 !     Gallium -     Hydrogen
      alpb(31, 6) =     2.188866d0 !     Gallium -       Carbon
      xfac(31, 6) =     1.617568d0 !     Gallium -       Carbon
      alpb(31, 7) =     1.949999d0 !     Gallium -     Nitrogen
      xfac(31, 7) =     0.867734d0 !     Gallium -     Nitrogen
      alpb(31, 8) =     2.408216d0 !     Gallium -       Oxygen
      xfac(31, 8) =     1.379976d0 !     Gallium -       Oxygen
      alpb(31, 9) =     3.055971d0 !     Gallium -     Fluorine
      xfac(31, 9) =     2.319957d0 !     Gallium -     Fluorine
      alpb(31,14) =     2.169690d0 !     Gallium -      Silicon
      xfac(31,14) =     5.031330d0 !     Gallium -      Silicon
      alpb(31,15) =     1.600000d0 !     Gallium -   Phosphorus
      xfac(31,15) =     4.000000d0 !     Gallium -   Phosphorus
      alpb(31,16) =     2.514000d0 !     Gallium -       Sulfur
      xfac(31,16) =     4.204343d0 !     Gallium -       Sulfur
      alpb(31,17) =     2.104228d0 !     Gallium -     Chlorine
      xfac(31,17) =     1.129276d0 !     Gallium -     Chlorine
      alpb(31,31) =     2.390223d0 !     Gallium -      Gallium
      xfac(31,31) =    11.941483d0 !     Gallium -      Gallium
 !
      alpb(32, 1) =     2.470301d0 !   Germanium -     Hydrogen
      xfac(32, 1) =     2.398259d0 !   Germanium -     Hydrogen
      alpb(32, 6) =     2.351577d0 !   Germanium -       Carbon
      xfac(32, 6) =     1.605487d0 !   Germanium -       Carbon
      alpb(32, 7) =     2.239698d0 !   Germanium -     Nitrogen
      xfac(32, 7) =     1.028521d0 !   Germanium -     Nitrogen
      alpb(32, 8) =     2.217395d0 !   Germanium -       Oxygen
      xfac(32, 8) =     0.690557d0 !   Germanium -       Oxygen
      alpb(32, 9) =     1.727325d0 !   Germanium -     Fluorine
      xfac(32, 9) =     0.165644d0 !   Germanium -     Fluorine
      alpb(32,14) =     2.053934d0 !   Germanium -      Silicon
      xfac(32,14) =     3.121907d0 !   Germanium -      Silicon
      alpb(32,15) =     1.831652d0 !   Germanium -   Phosphorus
      xfac(32,15) =     4.212771d0 !   Germanium -   Phosphorus
      alpb(32,16) =     2.358433d0 !   Germanium -       Sulfur
      xfac(32,16) =     1.947726d0 !   Germanium -       Sulfur
      alpb(32,17) =     2.506796d0 !   Germanium -     Chlorine
      xfac(32,17) =     1.783333d0 !   Germanium -     Chlorine
      alpb(32,25) =     1.937769d0 !   Germanium -    Manganese
      xfac(32,25) =     2.470135d0 !   Germanium -    Manganese
      alpb(32,27) =     2.852610d0 !   Germanium -       Cobalt
      xfac(32,27) =     2.151850d0 !   Germanium -       Cobalt
      alpb(32,32) =     2.215455d0 !   Germanium -    Germanium
      xfac(32,32) =     5.884206d0 !   Germanium -    Germanium
 !
      alpb(33, 1) =     1.749762d0 !     Arsenic -     Hydrogen
      xfac(33, 1) =     0.763924d0 !     Arsenic -     Hydrogen
      alpb(33, 6) =     1.805305d0 !     Arsenic -       Carbon
      xfac(33, 6) =     0.604465d0 !     Arsenic -       Carbon
      alpb(33, 7) =     2.035339d0 !     Arsenic -     Nitrogen
      xfac(33, 7) =     0.784041d0 !     Arsenic -     Nitrogen
      alpb(33, 8) =     2.387990d0 !     Arsenic -       Oxygen
      xfac(33, 8) =     1.076670d0 !     Arsenic -       Oxygen
      alpb(33, 9) =     2.783517d0 !     Arsenic -     Fluorine
      xfac(33, 9) =     1.196884d0 !     Arsenic -     Fluorine
      alpb(33,11) =     1.763497d0 !     Arsenic -       Sodium
      xfac(33,11) =     2.673425d0 !     Arsenic -       Sodium
      alpb(33,13) =     1.332670d0 !     Arsenic -     Aluminum
      xfac(33,13) =     1.322056d0 !     Arsenic -     Aluminum
      alpb(33,14) =     1.771030d0 !     Arsenic -      Silicon
      xfac(33,14) =     1.384298d0 !     Arsenic -      Silicon
      alpb(33,16) =     1.826372d0 !     Arsenic -       Sulfur
      xfac(33,16) =     0.732648d0 !     Arsenic -       Sulfur
      alpb(33,17) =     1.947927d0 !     Arsenic -     Chlorine
      xfac(33,17) =     0.747436d0 !     Arsenic -     Chlorine
      alpb(33,19) =     1.267957d0 !     Arsenic -    Potassium
      xfac(33,19) =     2.276204d0 !     Arsenic -    Potassium
      alpb(33,22) =     1.711955d0 !     Arsenic -     Titanium
      xfac(33,22) =     1.371503d0 !     Arsenic -     Titanium
      alpb(33,27) =     1.514923d0 !     Arsenic -       Cobalt
      xfac(33,27) =     2.030232d0 !     Arsenic -       Cobalt
      alpb(33,30) =     1.618734d0 !     Arsenic -         Zinc
      xfac(33,30) =     2.700385d0 !     Arsenic -         Zinc
      alpb(33,31) =     1.534812d0 !     Arsenic -      Gallium
      xfac(33,31) =     1.196640d0 !     Arsenic -      Gallium
      alpb(33,33) =     1.707277d0 !     Arsenic -      Arsenic
      xfac(33,33) =     1.325873d0 !     Arsenic -      Arsenic
 !
      alpb(34, 1) =     2.547234d0 !    Selenium -     Hydrogen
      xfac(34, 1) =     1.229099d0 !    Selenium -     Hydrogen
      alpb(34, 6) =     2.186857d0 !    Selenium -       Carbon
      xfac(34, 6) =     0.654796d0 !    Selenium -       Carbon
      alpb(34, 7) =     1.980885d0 !    Selenium -     Nitrogen
      xfac(34, 7) =     0.448537d0 !    Selenium -     Nitrogen
      alpb(34, 8) =     2.612643d0 !    Selenium -       Oxygen
      xfac(34, 8) =     0.860233d0 !    Selenium -       Oxygen
      alpb(34, 9) =     2.463196d0 !    Selenium -     Fluorine
      xfac(34, 9) =     0.473969d0 !    Selenium -     Fluorine
      alpb(34,11) =     1.115555d0 !    Selenium -       Sodium
      xfac(34,11) =     0.902628d0 !    Selenium -       Sodium
      alpb(34,14) =     2.318601d0 !    Selenium -      Silicon
      xfac(34,14) =     2.051717d0 !    Selenium -      Silicon
      alpb(34,15) =     1.865719d0 !    Selenium -   Phosphorus
      xfac(34,15) =     2.359419d0 !    Selenium -   Phosphorus
      alpb(34,16) =     1.492756d0 !    Selenium -       Sulfur
      xfac(34,16) =     0.530796d0 !    Selenium -       Sulfur
      alpb(34,17) =     2.170000d0 !    Selenium -     Chlorine
      xfac(34,17) =     0.869163d0 !    Selenium -     Chlorine
      alpb(34,19) =     1.680151d0 !    Selenium -    Potassium
      xfac(34,19) =     3.871380d0 !    Selenium -    Potassium
      alpb(34,25) =     1.981410d0 !    Selenium -    Manganese
      xfac(34,25) =     2.170787d0 !    Selenium -    Manganese
      alpb(34,27) =     2.523450d0 !    Selenium -       Cobalt
      xfac(34,27) =     2.202410d0 !    Selenium -       Cobalt
      alpb(34,30) =     1.163289d0 !    Selenium -         Zinc
      xfac(34,30) =     0.367711d0 !    Selenium -         Zinc
      alpb(34,32) =     1.604107d0 !    Selenium -    Germanium
      xfac(34,32) =     0.556002d0 !    Selenium -    Germanium
      alpb(34,33) =     1.514823d0 !    Selenium -      Arsenic
      xfac(34,33) =     0.541956d0 !    Selenium -      Arsenic
      alpb(34,34) =     1.524158d0 !    Selenium -     Selenium
      xfac(34,34) =     0.334506d0 !    Selenium -     Selenium
 !
      alpb(35, 1) =     2.339252d0 !     Bromine -     Hydrogen
      xfac(35, 1) =     1.270390d0 !     Bromine -     Hydrogen
      alpb(35, 2) =     2.127598d0 !     Bromine -       Helium
      xfac(35, 2) =     1.062013d0 !     Bromine -       Helium
      alpb(35, 3) =     2.143819d0 !     Bromine -      Lithium
      xfac(35, 3) =     2.241404d0 !     Bromine -      Lithium
      alpb(35, 4) =     2.283569d0 !     Bromine -    Beryllium
      xfac(35, 4) =     2.659130d0 !     Bromine -    Beryllium
      alpb(35, 5) =     2.307098d0 !     Bromine -        Boron
      xfac(35, 5) =     1.849590d0 !     Bromine -        Boron
      alpb(35, 6) =     2.252349d0 !     Bromine -       Carbon
      xfac(35, 6) =     0.968921d0 !     Bromine -       Carbon
      alpb(35, 7) =     3.015469d0 !     Bromine -     Nitrogen
      xfac(35, 7) =     4.148435d0 !     Bromine -     Nitrogen
      alpb(35, 8) =     2.739280d0 !     Bromine -       Oxygen
      xfac(35, 8) =     1.425004d0 !     Bromine -       Oxygen
      alpb(35, 9) =     2.756179d0 !     Bromine -     Fluorine
      xfac(35, 9) =     0.915922d0 !     Bromine -     Fluorine
      alpb(35,10) =     2.483203d0 !     Bromine -         Neon
      xfac(35,10) =     1.001506d0 !     Bromine -         Neon
      alpb(35,11) =     2.327183d0 !     Bromine -       Sodium
      xfac(35,11) =    11.511433d0 !     Bromine -       Sodium
      alpb(35,12) =     2.350023d0 !     Bromine -    Magnesium
      xfac(35,12) =     7.183008d0 !     Bromine -    Magnesium
      alpb(35,13) =     1.514344d0 !     Bromine -     Aluminum
      xfac(35,13) =     1.504277d0 !     Bromine -     Aluminum
      alpb(35,14) =     1.715815d0 !     Bromine -      Silicon
      xfac(35,14) =     1.126731d0 !     Bromine -      Silicon
      alpb(35,15) =     1.664343d0 !     Bromine -   Phosphorus
      xfac(35,15) =     1.746224d0 !     Bromine -   Phosphorus
      alpb(35,16) =     2.099922d0 !     Bromine -       Sulfur
      xfac(35,16) =     1.004759d0 !     Bromine -       Sulfur
      alpb(35,17) =     1.906403d0 !     Bromine -     Chlorine
      xfac(35,17) =     0.581542d0 !     Bromine -     Chlorine
      alpb(35,18) =     2.454724d0 !     Bromine -        Argon
      xfac(35,18) =     3.261699d0 !     Bromine -        Argon
      alpb(35,19) =     1.887799d0 !     Bromine -    Potassium
      xfac(35,19) =     7.523969d0 !     Bromine -    Potassium
      alpb(35,20) =     2.558257d0 !     Bromine -      Calcium
      xfac(35,20) =    12.875179d0 !     Bromine -      Calcium
      alpb(35,21) =     1.531278d0 !     Bromine -     Scandium
      xfac(35,21) =     1.063920d0 !     Bromine -     Scandium
      alpb(35,22) =     1.760015d0 !     Bromine -     Titanium
      xfac(35,22) =     1.534465d0 !     Bromine -     Titanium
      alpb(35,23) =     1.909502d0 !     Bromine -     Vanadium
      xfac(35,23) =     1.394543d0 !     Bromine -     Vanadium
      alpb(35,24) =     1.781866d0 !     Bromine -     Chromium
      xfac(35,24) =     0.746857d0 !     Bromine -     Chromium
      alpb(35,25) =     2.183298d0 !     Bromine -    Manganese
      xfac(35,25) =     2.530861d0 !     Bromine -    Manganese
      alpb(35,26) =     2.388196d0 !     Bromine -         Iron
      xfac(35,26) =     1.415413d0 !     Bromine -         Iron
      alpb(35,27) =     2.124301d0 !     Bromine -       Cobalt
      xfac(35,27) =     0.759690d0 !     Bromine -       Cobalt
      alpb(35,28) =     2.543159d0 !     Bromine -       Nickel
      xfac(35,28) =     1.656698d0 !     Bromine -       Nickel
      alpb(35,29) =     3.040037d0 !     Bromine -       Copper
      xfac(35,29) =     2.655465d0 !     Bromine -       Copper
      alpb(35,30) =     1.594962d0 !     Bromine -         Zinc
      xfac(35,30) =     1.241996d0 !     Bromine -         Zinc
      alpb(35,31) =     1.934418d0 !     Bromine -      Gallium
      xfac(35,31) =     1.722754d0 !     Bromine -      Gallium
      alpb(35,32) =     2.062366d0 !     Bromine -    Germanium
      xfac(35,32) =     2.031652d0 !     Bromine -    Germanium
      alpb(35,33) =     1.750449d0 !     Bromine -      Arsenic
      xfac(35,33) =     0.949379d0 !     Bromine -      Arsenic
      alpb(35,34) =     1.788806d0 !     Bromine -     Selenium
      xfac(35,34) =     0.682982d0 !     Bromine -     Selenium
      alpb(35,35) =     2.147378d0 !     Bromine -      Bromine
      xfac(35,35) =     1.599562d0 !     Bromine -      Bromine
 !
      alpb(36, 1) =     3.770453d0 !     Krypton -     Hydrogen
      xfac(36, 1) =     5.125897d0 !     Krypton -     Hydrogen
      alpb(36, 2) =     1.996943d0 !     Krypton -       Helium
      xfac(36, 2) =     0.627701d0 !     Krypton -       Helium
      alpb(36, 3) =     3.004783d0 !     Krypton -      Lithium
      xfac(36, 3) =     8.377143d0 !     Krypton -      Lithium
      alpb(36, 4) =     3.289764d0 !     Krypton -    Beryllium
      xfac(36, 4) =    10.264026d0 !     Krypton -    Beryllium
      alpb(36, 5) =     2.559201d0 !     Krypton -        Boron
      xfac(36, 5) =     2.931148d0 !     Krypton -        Boron
      alpb(36, 6) =     2.076738d0 !     Krypton -       Carbon
      xfac(36, 6) =     0.652623d0 !     Krypton -       Carbon
      alpb(36, 7) =     1.644052d0 !     Krypton -     Nitrogen
      xfac(36, 7) =     0.199606d0 !     Krypton -     Nitrogen
      alpb(36, 8) =     0.297001d0 !     Krypton -       Oxygen
      xfac(36, 8) =     0.004690d0 !     Krypton -       Oxygen
      alpb(36, 9) =     3.452321d0 !     Krypton -     Fluorine
      xfac(36, 9) =     4.134407d0 !     Krypton -     Fluorine
      alpb(36,10) =     2.813679d0 !     Krypton -         Neon
      xfac(36,10) =     1.433722d0 !     Krypton -         Neon
      alpb(36,11) =     2.562062d0 !     Krypton -       Sodium
      xfac(36,11) =     9.817859d0 !     Krypton -       Sodium
      alpb(36,12) =     1.296221d0 !     Krypton -    Magnesium
      xfac(36,12) =     1.119449d0 !     Krypton -    Magnesium
      alpb(36,13) =     2.493834d0 !     Krypton -     Aluminum
      xfac(36,13) =     5.076857d0 !     Krypton -     Aluminum
      alpb(36,14) =     1.545354d0 !     Krypton -      Silicon
      xfac(36,14) =     0.639030d0 !     Krypton -      Silicon
      alpb(36,17) =     1.884662d0 !     Krypton -     Chlorine
      xfac(36,17) =     0.520353d0 !     Krypton -     Chlorine
      alpb(36,18) =     1.995125d0 !     Krypton -        Argon
      xfac(36,18) =     0.554874d0 !     Krypton -        Argon
      alpb(36,19) =     2.296640d0 !     Krypton -    Potassium
      xfac(36,19) =     8.532309d0 !     Krypton -    Potassium
      alpb(36,20) =     1.559229d0 !     Krypton -      Calcium
      xfac(36,20) =     1.305808d0 !     Krypton -      Calcium
      alpb(36,35) =     1.608300d0 !     Krypton -      Bromine
      xfac(36,35) =     0.499653d0 !     Krypton -      Bromine
      alpb(36,36) =     1.913342d0 !     Krypton -      Krypton
      xfac(36,36) =     0.252431d0 !     Krypton -      Krypton
 !
      alpb(37, 1) =     1.890495d0 !    Rubidium -     Hydrogen
      xfac(37, 1) =     4.316836d0 !    Rubidium -     Hydrogen
      alpb(37, 2) =     1.543436d0 !    Rubidium -       Helium
      xfac(37, 2) =     1.804024d0 !    Rubidium -       Helium
      alpb(37, 5) =     2.989999d0 !    Rubidium -        Boron
      xfac(37, 5) =    10.280532d0 !    Rubidium -        Boron
      alpb(37, 6) =     2.287377d0 !    Rubidium -       Carbon
      xfac(37, 6) =    24.603216d0 !    Rubidium -       Carbon
      alpb(37, 7) =     2.205212d0 !    Rubidium -     Nitrogen
      xfac(37, 7) =    19.919815d0 !    Rubidium -     Nitrogen
      alpb(37, 8) =     1.572166d0 !    Rubidium -       Oxygen
      xfac(37, 8) =     0.791546d0 !    Rubidium -       Oxygen
      alpb(37, 9) =     3.131045d0 !    Rubidium -     Fluorine
      xfac(37, 9) =     2.629683d0 !    Rubidium -     Fluorine
      alpb(37,10) =     2.429710d0 !    Rubidium -         Neon
      xfac(37,10) =     7.683406d0 !    Rubidium -         Neon
      alpb(37,13) =     0.931060d0 !    Rubidium -     Aluminum
      xfac(37,13) =    19.138062d0 !    Rubidium -     Aluminum
      alpb(37,15) =     0.922429d0 !    Rubidium -   Phosphorus
      xfac(37,15) =     0.526907d0 !    Rubidium -   Phosphorus
      alpb(37,16) =     1.285680d0 !    Rubidium -       Sulfur
      xfac(37,16) =     1.380751d0 !    Rubidium -       Sulfur
      alpb(37,17) =     1.349244d0 !    Rubidium -     Chlorine
      xfac(37,17) =     0.714916d0 !    Rubidium -     Chlorine
      alpb(37,18) =     2.581073d0 !    Rubidium -        Argon
      xfac(37,18) =    18.431817d0 !    Rubidium -        Argon
      alpb(37,19) =     1.719577d0 !    Rubidium -    Potassium
      xfac(37,19) =     1.174003d0 !    Rubidium -    Potassium
      alpb(37,23) =     2.024277d0 !    Rubidium -     Vanadium
      xfac(37,23) =    12.360809d0 !    Rubidium -     Vanadium
      alpb(37,24) =     1.569973d0 !    Rubidium -     Chromium
      xfac(37,24) =     5.846691d0 !    Rubidium -     Chromium
      alpb(37,35) =     1.737339d0 !    Rubidium -      Bromine
      xfac(37,35) =     7.887658d0 !    Rubidium -      Bromine
      alpb(37,36) =     2.413547d0 !    Rubidium -      Krypton
      xfac(37,36) =    15.315026d0 !    Rubidium -      Krypton
      alpb(37,37) =     0.684616d0 !    Rubidium -     Rubidium
      xfac(37,37) =     4.280624d0 !    Rubidium -     Rubidium
 !
      alpb(38, 1) =     2.918082d0 !   Strontium -     Hydrogen
      xfac(38, 1) =    18.457472d0 !   Strontium -     Hydrogen
      alpb(38, 6) =     2.802150d0 !   Strontium -       Carbon
      xfac(38, 6) =    11.414713d0 !   Strontium -       Carbon
      alpb(38, 7) =     3.025117d0 !   Strontium -     Nitrogen
      xfac(38, 7) =    14.103810d0 !   Strontium -     Nitrogen
      alpb(38, 8) =     3.921499d0 !   Strontium -       Oxygen
      xfac(38, 8) =    17.067149d0 !   Strontium -       Oxygen
      alpb(38, 9) =     3.309755d0 !   Strontium -     Fluorine
      xfac(38, 9) =     3.758667d0 !   Strontium -     Fluorine
      alpb(38,14) =     2.416117d0 !   Strontium -      Silicon
      xfac(38,14) =    29.996030d0 !   Strontium -      Silicon
      alpb(38,15) =     2.273841d0 !   Strontium -   Phosphorus
      xfac(38,15) =    23.946650d0 !   Strontium -   Phosphorus
      alpb(38,16) =     3.325295d0 !   Strontium -       Sulfur
      xfac(38,16) =    41.563327d0 !   Strontium -       Sulfur
      alpb(38,17) =     3.501974d0 !   Strontium -     Chlorine
      xfac(38,17) =    39.960719d0 !   Strontium -     Chlorine
      alpb(38,22) =     2.880030d0 !   Strontium -     Titanium
      xfac(38,22) =     2.817250d0 !   Strontium -     Titanium
      alpb(38,31) =     1.489463d0 !   Strontium -      Gallium
      xfac(38,31) =     2.800419d0 !   Strontium -      Gallium
      alpb(38,35) =     3.086374d0 !   Strontium -      Bromine
      xfac(38,35) =    19.218824d0 !   Strontium -      Bromine
      alpb(38,38) =     2.194036d0 !   Strontium -    Strontium
      xfac(38,38) =    31.817350d0 !   Strontium -    Strontium
 !
      alpb(39, 1) =     2.322175d0 !     Yttrium -     Hydrogen
      xfac(39, 1) =     6.935667d0 !     Yttrium -     Hydrogen
      alpb(39, 3) =     1.212009d0 !     Yttrium -      Lithium
      xfac(39, 3) =     0.577598d0 !     Yttrium -      Lithium
      alpb(39, 6) =     2.541211d0 !     Yttrium -       Carbon
      xfac(39, 6) =    19.957240d0 !     Yttrium -       Carbon
      alpb(39, 7) =     2.084245d0 !     Yttrium -     Nitrogen
      xfac(39, 7) =     3.253368d0 !     Yttrium -     Nitrogen
      alpb(39, 8) =     2.086475d0 !     Yttrium -       Oxygen
      xfac(39, 8) =     1.424444d0 !     Yttrium -       Oxygen
      alpb(39, 9) =     3.245964d0 !     Yttrium -     Fluorine
      xfac(39, 9) =     9.257528d0 !     Yttrium -     Fluorine
      alpb(39,13) =     1.003500d0 !     Yttrium -     Aluminum
      xfac(39,13) =     0.500670d0 !     Yttrium -     Aluminum
      alpb(39,14) =     2.016820d0 !     Yttrium -      Silicon
      xfac(39,14) =     3.219030d0 !     Yttrium -      Silicon
      alpb(39,15) =     1.172165d0 !     Yttrium -   Phosphorus
      xfac(39,15) =     1.726458d0 !     Yttrium -   Phosphorus
      alpb(39,16) =     1.345475d0 !     Yttrium -       Sulfur
      xfac(39,16) =     0.961448d0 !     Yttrium -       Sulfur
      alpb(39,17) =     1.882700d0 !     Yttrium -     Chlorine
      xfac(39,17) =     2.186706d0 !     Yttrium -     Chlorine
      alpb(39,19) =     0.947193d0 !     Yttrium -    Potassium
      xfac(39,19) =     1.143281d0 !     Yttrium -    Potassium
      alpb(39,35) =     1.359600d0 !     Yttrium -      Bromine
      xfac(39,35) =     1.090173d0 !     Yttrium -      Bromine
      alpb(39,39) =     1.533049d0 !     Yttrium -      Yttrium
      xfac(39,39) =    15.620872d0 !     Yttrium -      Yttrium
 !
      alpb(40, 1) =     1.536594d0 !   Zirconium -     Hydrogen
      xfac(40, 1) =     0.414278d0 !   Zirconium -     Hydrogen
      alpb(40, 6) =     1.738320d0 !   Zirconium -       Carbon
      xfac(40, 6) =     0.715392d0 !   Zirconium -       Carbon
      alpb(40, 7) =     1.986255d0 !   Zirconium -     Nitrogen
      xfac(40, 7) =     1.386137d0 !   Zirconium -     Nitrogen
      alpb(40, 8) =     2.093741d0 !   Zirconium -       Oxygen
      xfac(40, 8) =     1.293822d0 !   Zirconium -       Oxygen
      alpb(40, 9) =     2.406399d0 !   Zirconium -     Fluorine
      xfac(40, 9) =     1.606098d0 !   Zirconium -     Fluorine
      alpb(40,13) =     1.270620d0 !   Zirconium -     Aluminum
      xfac(40,13) =     0.874060d0 !   Zirconium -     Aluminum
      alpb(40,14) =     1.605795d0 !   Zirconium -      Silicon
      xfac(40,14) =     1.566491d0 !   Zirconium -      Silicon
      alpb(40,15) =     0.963910d0 !   Zirconium -   Phosphorus
      xfac(40,15) =     1.020015d0 !   Zirconium -   Phosphorus
      alpb(40,16) =     0.957666d0 !   Zirconium -       Sulfur
      xfac(40,16) =     0.200522d0 !   Zirconium -       Sulfur
      alpb(40,17) =     2.352409d0 !   Zirconium -     Chlorine
      xfac(40,17) =     2.273630d0 !   Zirconium -     Chlorine
      alpb(40,35) =     1.617591d0 !   Zirconium -      Bromine
      xfac(40,35) =     1.243772d0 !   Zirconium -      Bromine
      alpb(40,40) =     2.714671d0 !   Zirconium -    Zirconium
      xfac(40,40) =    29.768192d0 !   Zirconium -    Zirconium
 !
      alpb(41, 1) =     2.321651d0 !     Niobium -     Hydrogen
      xfac(41, 1) =     6.958727d0 !     Niobium -     Hydrogen
      alpb(41, 6) =     2.277928d0 !     Niobium -       Carbon
      xfac(41, 6) =     1.991488d0 !     Niobium -       Carbon
      alpb(41, 7) =     2.810017d0 !     Niobium -     Nitrogen
      xfac(41, 7) =     4.374260d0 !     Niobium -     Nitrogen
      alpb(41, 8) =     2.715670d0 !     Niobium -       Oxygen
      xfac(41, 8) =     2.793681d0 !     Niobium -       Oxygen
      alpb(41, 9) =     3.115376d0 !     Niobium -     Fluorine
      xfac(41, 9) =     3.297982d0 !     Niobium -     Fluorine
      alpb(41,11) =     2.551010d0 !     Niobium -       Sodium
      xfac(41,11) =     8.276020d0 !     Niobium -       Sodium
      alpb(41,15) =     1.922968d0 !     Niobium -   Phosphorus
      xfac(41,15) =     6.219347d0 !     Niobium -   Phosphorus
      alpb(41,16) =     2.279550d0 !     Niobium -       Sulfur
      xfac(41,16) =     3.225637d0 !     Niobium -       Sulfur
      alpb(41,17) =     2.757523d0 !     Niobium -     Chlorine
      xfac(41,17) =     6.452483d0 !     Niobium -     Chlorine
      alpb(41,19) =     4.521360d0 !     Niobium -    Potassium
      xfac(41,19) =     2.026590d0 !     Niobium -    Potassium
      alpb(41,35) =     2.531918d0 !     Niobium -      Bromine
      xfac(41,35) =     8.316457d0 !     Niobium -      Bromine
      alpb(41,41) =     2.030464d0 !     Niobium -      Niobium
      xfac(41,41) =    10.027153d0 !     Niobium -      Niobium
 !
      alpb(42, 1) =     2.139004d0 !  Molybdenum -     Hydrogen
      xfac(42, 1) =     1.177934d0 !  Molybdenum -     Hydrogen
      alpb(42, 3) =     2.201335d0 !  Molybdenum -      Lithium
      xfac(42, 3) =     5.209247d0 !  Molybdenum -      Lithium
      alpb(42, 6) =     2.140063d0 !  Molybdenum -       Carbon
      xfac(42, 6) =     1.042667d0 !  Molybdenum -       Carbon
      alpb(42, 7) =     2.293955d0 !  Molybdenum -     Nitrogen
      xfac(42, 7) =     1.330858d0 !  Molybdenum -     Nitrogen
      alpb(42, 8) =     2.197353d0 !  Molybdenum -       Oxygen
      xfac(42, 8) =     0.864597d0 !  Molybdenum -       Oxygen
      alpb(42, 9) =     2.593518d0 !  Molybdenum -     Fluorine
      xfac(42, 9) =     1.107779d0 !  Molybdenum -     Fluorine
      alpb(42,11) =     2.440770d0 !  Molybdenum -       Sodium
      xfac(42,11) =     8.286550d0 !  Molybdenum -       Sodium
      alpb(42,15) =     1.850441d0 !  Molybdenum -   Phosphorus
      xfac(42,15) =     1.522846d0 !  Molybdenum -   Phosphorus
      alpb(42,16) =     2.343350d0 !  Molybdenum -       Sulfur
      xfac(42,16) =     1.822187d0 !  Molybdenum -       Sulfur
      alpb(42,17) =     2.358706d0 !  Molybdenum -     Chlorine
      xfac(42,17) =     1.570824d0 !  Molybdenum -     Chlorine
      alpb(42,19) =     1.594941d0 !  Molybdenum -    Potassium
      xfac(42,19) =    10.522232d0 !  Molybdenum -    Potassium
      alpb(42,24) =     1.873232d0 !  Molybdenum -     Chromium
      xfac(42,24) =     0.484096d0 !  Molybdenum -     Chromium
      alpb(42,26) =     2.239581d0 !  Molybdenum -         Iron
      xfac(42,26) =     3.257115d0 !  Molybdenum -         Iron
      alpb(42,35) =     1.934589d0 !  Molybdenum -      Bromine
      xfac(42,35) =     1.291525d0 !  Molybdenum -      Bromine
      alpb(42,37) =     2.971399d0 !  Molybdenum -     Rubidium
      xfac(42,37) =     0.874676d0 !  Molybdenum -     Rubidium
      alpb(42,42) =     1.078447d0 !  Molybdenum -   Molybdenum
      xfac(42,42) =     0.223956d0 !  Molybdenum -   Molybdenum
 !
      alpb(43, 1) =     2.576199d0 !  Technetium -     Hydrogen
      xfac(43, 1) =     5.418951d0 !  Technetium -     Hydrogen
      alpb(43, 6) =     2.815972d0 !  Technetium -       Carbon
      xfac(43, 6) =     3.999428d0 !  Technetium -       Carbon
      alpb(43, 7) =     2.177956d0 !  Technetium -     Nitrogen
      xfac(43, 7) =     0.980071d0 !  Technetium -     Nitrogen
      alpb(43, 8) =     2.535619d0 !  Technetium -       Oxygen
      xfac(43, 8) =     1.303538d0 !  Technetium -       Oxygen
      alpb(43, 9) =     3.385092d0 !  Technetium -     Fluorine
      xfac(43, 9) =     3.880884d0 !  Technetium -     Fluorine
      alpb(43,15) =     0.930051d0 !  Technetium -   Phosphorus
      xfac(43,15) =     0.470758d0 !  Technetium -   Phosphorus
      alpb(43,16) =     2.141702d0 !  Technetium -       Sulfur
      xfac(43,16) =     1.449910d0 !  Technetium -       Sulfur
      alpb(43,17) =     2.360242d0 !  Technetium -     Chlorine
      xfac(43,17) =     1.744657d0 !  Technetium -     Chlorine
      alpb(43,32) =     2.852820d0 !  Technetium -    Germanium
      xfac(43,32) =     2.152060d0 !  Technetium -    Germanium
      alpb(43,34) =     2.523660d0 !  Technetium -     Selenium
      xfac(43,34) =     2.202620d0 !  Technetium -     Selenium
      alpb(43,35) =     2.688330d0 !  Technetium -      Bromine
      xfac(43,35) =     6.426037d0 !  Technetium -      Bromine
      alpb(43,43) =     2.153000d0 !  Technetium -   Technetium
      xfac(43,43) =     2.572063d0 !  Technetium -   Technetium
 !
      alpb(44, 1) =     3.031003d0 !   Ruthenium -     Hydrogen
      xfac(44, 1) =     0.462490d0 !   Ruthenium -     Hydrogen
      alpb(44, 6) =     2.661734d0 !   Ruthenium -       Carbon
      xfac(44, 6) =     0.434352d0 !   Ruthenium -       Carbon
      alpb(44, 7) =     1.951233d0 !   Ruthenium -     Nitrogen
      xfac(44, 7) =     0.271221d0 !   Ruthenium -     Nitrogen
      alpb(44, 8) =     1.928484d0 !   Ruthenium -       Oxygen
      xfac(44, 8) =     0.339590d0 !   Ruthenium -       Oxygen
      alpb(44, 9) =     2.719488d0 !   Ruthenium -     Fluorine
      xfac(44, 9) =     0.680978d0 !   Ruthenium -     Fluorine
      alpb(44,14) =     2.775910d0 !   Ruthenium -      Silicon
      xfac(44,14) =     0.849430d0 !   Ruthenium -      Silicon
      alpb(44,15) =     1.440298d0 !   Ruthenium -   Phosphorus
      xfac(44,15) =     0.482587d0 !   Ruthenium -   Phosphorus
      alpb(44,16) =     3.002139d0 !   Ruthenium -       Sulfur
      xfac(44,16) =     0.788319d0 !   Ruthenium -       Sulfur
      alpb(44,17) =     3.340740d0 !   Ruthenium -     Chlorine
      xfac(44,17) =     1.986295d0 !   Ruthenium -     Chlorine
      alpb(44,32) =     2.852320d0 !   Ruthenium -    Germanium
      xfac(44,32) =     2.151560d0 !   Ruthenium -    Germanium
      alpb(44,34) =     2.523160d0 !   Ruthenium -     Selenium
      xfac(44,34) =     2.202120d0 !   Ruthenium -     Selenium
      alpb(44,35) =     2.611647d0 !   Ruthenium -      Bromine
      xfac(44,35) =     3.893512d0 !   Ruthenium -      Bromine
      alpb(44,44) =     2.341541d0 !   Ruthenium -    Ruthenium
      xfac(44,44) =     0.984874d0 !   Ruthenium -    Ruthenium
 !
      alpb(45, 1) =     2.716287d0 !     Rhodium -     Hydrogen
      xfac(45, 1) =     1.728302d0 !     Rhodium -     Hydrogen
      alpb(45, 5) =     2.400000d0 !     Rhodium -        Boron
      xfac(45, 5) =     2.000000d0 !     Rhodium -        Boron
      alpb(45, 6) =     3.007700d0 !     Rhodium -       Carbon
      xfac(45, 6) =     0.562962d0 !     Rhodium -       Carbon
      alpb(45, 7) =     3.028135d0 !     Rhodium -     Nitrogen
      xfac(45, 7) =     1.013618d0 !     Rhodium -     Nitrogen
      alpb(45, 8) =     3.452408d0 !     Rhodium -       Oxygen
      xfac(45, 8) =     1.534037d0 !     Rhodium -       Oxygen
      alpb(45, 9) =     3.083507d0 !     Rhodium -     Fluorine
      xfac(45, 9) =     0.772245d0 !     Rhodium -     Fluorine
      alpb(45,14) =     2.776490d0 !     Rhodium -      Silicon
      xfac(45,14) =     0.850010d0 !     Rhodium -      Silicon
      alpb(45,15) =     2.236601d0 !     Rhodium -   Phosphorus
      xfac(45,15) =     0.738916d0 !     Rhodium -   Phosphorus
      alpb(45,16) =     3.005420d0 !     Rhodium -       Sulfur
      xfac(45,16) =     0.970563d0 !     Rhodium -       Sulfur
      alpb(45,17) =     3.542676d0 !     Rhodium -     Chlorine
      xfac(45,17) =     0.628186d0 !     Rhodium -     Chlorine
      alpb(45,32) =     2.852900d0 !     Rhodium -    Germanium
      xfac(45,32) =     2.152140d0 !     Rhodium -    Germanium
      alpb(45,34) =     2.523740d0 !     Rhodium -     Selenium
      xfac(45,34) =     2.202700d0 !     Rhodium -     Selenium
      alpb(45,35) =     2.893677d0 !     Rhodium -      Bromine
      xfac(45,35) =     1.441509d0 !     Rhodium -      Bromine
      alpb(45,45) =     3.281577d0 !     Rhodium -      Rhodium
      xfac(45,45) =    17.154616d0 !     Rhodium -      Rhodium
 !
      alpb(46, 1) =     3.052992d0 !   Palladium -     Hydrogen
      xfac(46, 1) =     0.675244d0 !   Palladium -     Hydrogen
      alpb(46, 6) =     1.449994d0 !   Palladium -       Carbon
      xfac(46, 6) =     0.040769d0 !   Palladium -       Carbon
      alpb(46, 7) =     2.319285d0 !   Palladium -     Nitrogen
      xfac(46, 7) =     0.327063d0 !   Palladium -     Nitrogen
      alpb(46, 8) =     2.362481d0 !   Palladium -       Oxygen
      xfac(46, 8) =     0.394849d0 !   Palladium -       Oxygen
      alpb(46, 9) =     3.117188d0 !   Palladium -     Fluorine
      xfac(46, 9) =     0.610235d0 !   Palladium -     Fluorine
      alpb(46,13) =     1.572720d0 !   Palladium -     Aluminum
      xfac(46,13) =     1.057290d0 !   Palladium -     Aluminum
      alpb(46,14) =     2.714212d0 !   Palladium -      Silicon
      xfac(46,14) =     1.381243d0 !   Palladium -      Silicon
      alpb(46,15) =     0.876896d0 !   Palladium -   Phosphorus
      xfac(46,15) =     0.223289d0 !   Palladium -   Phosphorus
      alpb(46,16) =     3.134436d0 !   Palladium -       Sulfur
      xfac(46,16) =     0.568359d0 !   Palladium -       Sulfur
      alpb(46,17) =     2.966363d0 !   Palladium -     Chlorine
      xfac(46,17) =     0.764165d0 !   Palladium -     Chlorine
      alpb(46,35) =     2.087790d0 !   Palladium -      Bromine
      xfac(46,35) =     0.491172d0 !   Palladium -      Bromine
      alpb(46,46) =     1.712149d0 !   Palladium -    Palladium
      xfac(46,46) =     0.297913d0 !   Palladium -    Palladium
 !
      alpb(47, 1) =     1.866268d0 !      Silver -     Hydrogen
      xfac(47, 1) =     0.669745d0 !      Silver -     Hydrogen
      alpb(47, 5) =     1.454270d0 !      Silver -        Boron
      xfac(47, 5) =     2.733745d0 !      Silver -        Boron
      alpb(47, 6) =     2.401775d0 !      Silver -       Carbon
      xfac(47, 6) =     1.108319d0 !      Silver -       Carbon
      alpb(47, 7) =     2.835438d0 !      Silver -     Nitrogen
      xfac(47, 7) =     1.090232d0 !      Silver -     Nitrogen
      alpb(47, 8) =     2.453629d0 !      Silver -       Oxygen
      xfac(47, 8) =     0.372450d0 !      Silver -       Oxygen
      alpb(47, 9) =     3.119532d0 !      Silver -     Fluorine
      xfac(47, 9) =     0.897783d0 !      Silver -     Fluorine
      alpb(47,13) =     1.683750d0 !      Silver -     Aluminum
      xfac(47,13) =     1.093559d0 !      Silver -     Aluminum
      alpb(47,15) =     1.305572d0 !      Silver -   Phosphorus
      xfac(47,15) =     0.482631d0 !      Silver -   Phosphorus
      alpb(47,16) =     2.575670d0 !      Silver -       Sulfur
      xfac(47,16) =     1.817701d0 !      Silver -       Sulfur
      alpb(47,17) =     3.198107d0 !      Silver -     Chlorine
      xfac(47,17) =     3.386746d0 !      Silver -     Chlorine
      alpb(47,19) =     2.092259d0 !      Silver -    Potassium
      xfac(47,19) =     5.619211d0 !      Silver -    Potassium
      alpb(47,24) =     2.700428d0 !      Silver -     Chromium
      xfac(47,24) =    21.083639d0 !      Silver -     Chromium
      alpb(47,35) =     3.259287d0 !      Silver -      Bromine
      xfac(47,35) =     6.111850d0 !      Silver -      Bromine
      alpb(47,46) =     4.000000d0 !      Silver -    Palladium
      xfac(47,46) =     2.000000d0 !      Silver -    Palladium
      alpb(47,47) =     1.489404d0 !      Silver -       Silver
      xfac(47,47) =     0.178879d0 !      Silver -       Silver
 !
      alpb(48, 1) =     1.875490d0 !     Cadmium -     Hydrogen
      xfac(48, 1) =     3.377913d0 !     Cadmium -     Hydrogen
      alpb(48, 6) =     1.940388d0 !     Cadmium -       Carbon
      xfac(48, 6) =     3.855359d0 !     Cadmium -       Carbon
      alpb(48, 7) =     1.769441d0 !     Cadmium -     Nitrogen
      xfac(48, 7) =     1.481460d0 !     Cadmium -     Nitrogen
      alpb(48, 8) =     2.668165d0 !     Cadmium -       Oxygen
      xfac(48, 8) =     5.349517d0 !     Cadmium -       Oxygen
      alpb(48, 9) =     3.174783d0 !     Cadmium -     Fluorine
      xfac(48, 9) =     8.351869d0 !     Cadmium -     Fluorine
      alpb(48,11) =     2.000000d0 !     Cadmium -       Sodium
      xfac(48,11) =     6.000000d0 !     Cadmium -       Sodium
      alpb(48,14) =     1.286882d0 !     Cadmium -      Silicon
      xfac(48,14) =     2.345912d0 !     Cadmium -      Silicon
      alpb(48,16) =     1.735391d0 !     Cadmium -       Sulfur
      xfac(48,16) =     2.929257d0 !     Cadmium -       Sulfur
      alpb(48,17) =     1.870170d0 !     Cadmium -     Chlorine
      xfac(48,17) =     2.254752d0 !     Cadmium -     Chlorine
      alpb(48,19) =     1.033580d0 !     Cadmium -    Potassium
      xfac(48,19) =     2.093242d0 !     Cadmium -    Potassium
      alpb(48,34) =     1.881368d0 !     Cadmium -     Selenium
      xfac(48,34) =     6.139995d0 !     Cadmium -     Selenium
      alpb(48,35) =     1.918455d0 !     Cadmium -      Bromine
      xfac(48,35) =     5.550415d0 !     Cadmium -      Bromine
      alpb(48,48) =     1.428097d0 !     Cadmium -      Cadmium
      xfac(48,48) =    10.662907d0 !     Cadmium -      Cadmium
 !
      alpb(49, 1) =     1.852461d0 !      Indium -     Hydrogen
      xfac(49, 1) =     1.773147d0 !      Indium -     Hydrogen
      alpb(49, 5) =     1.735480d0 !      Indium -        Boron
      xfac(49, 5) =     1.951651d0 !      Indium -        Boron
      alpb(49, 6) =     1.810115d0 !      Indium -       Carbon
      xfac(49, 6) =     1.041540d0 !      Indium -       Carbon
      alpb(49, 7) =     2.052217d0 !      Indium -     Nitrogen
      xfac(49, 7) =     1.529722d0 !      Indium -     Nitrogen
      alpb(49, 8) =     2.178110d0 !      Indium -       Oxygen
      xfac(49, 8) =     1.467957d0 !      Indium -       Oxygen
      alpb(49, 9) =     2.319418d0 !      Indium -     Fluorine
      xfac(49, 9) =     1.018315d0 !      Indium -     Fluorine
      alpb(49,16) =     2.430104d0 !      Indium -       Sulfur
      xfac(49,16) =     4.796933d0 !      Indium -       Sulfur
      alpb(49,17) =     2.211880d0 !      Indium -     Chlorine
      xfac(49,17) =     2.224354d0 !      Indium -     Chlorine
      alpb(49,31) =     1.596053d0 !      Indium -      Gallium
      xfac(49,31) =     2.473577d0 !      Indium -      Gallium
      alpb(49,33) =     1.520977d0 !      Indium -      Arsenic
      xfac(49,33) =     1.375570d0 !      Indium -      Arsenic
      alpb(49,34) =     1.362364d0 !      Indium -     Selenium
      xfac(49,34) =     0.598029d0 !      Indium -     Selenium
      alpb(49,35) =     1.862313d0 !      Indium -      Bromine
      xfac(49,35) =     2.138634d0 !      Indium -      Bromine
      alpb(49,37) =     0.859259d0 !      Indium -     Rubidium
      xfac(49,37) =     4.688357d0 !      Indium -     Rubidium
      alpb(49,49) =     2.601789d0 !      Indium -       Indium
      xfac(49,49) =    24.204383d0 !      Indium -       Indium
 !
      alpb(50, 1) =     1.855042d0 !         Tin -     Hydrogen
      xfac(50, 1) =     1.459969d0 !         Tin -     Hydrogen
      alpb(50, 6) =     1.818782d0 !         Tin -       Carbon
      xfac(50, 6) =     0.961947d0 !         Tin -       Carbon
      alpb(50, 7) =     1.783560d0 !         Tin -     Nitrogen
      xfac(50, 7) =     0.731228d0 !         Tin -     Nitrogen
      alpb(50, 8) =     1.959102d0 !         Tin -       Oxygen
      xfac(50, 8) =     0.723272d0 !         Tin -       Oxygen
      alpb(50, 9) =     2.593459d0 !         Tin -     Fluorine
      xfac(50, 9) =     1.477352d0 !         Tin -     Fluorine
      alpb(50,13) =     1.597939d0 !         Tin -     Aluminum
      xfac(50,13) =     2.367990d0 !         Tin -     Aluminum
      alpb(50,16) =     2.065722d0 !         Tin -       Sulfur
      xfac(50,16) =     1.909070d0 !         Tin -       Sulfur
      alpb(50,17) =     1.887044d0 !         Tin -     Chlorine
      xfac(50,17) =     0.944374d0 !         Tin -     Chlorine
      alpb(50,19) =     2.238329d0 !         Tin -    Potassium
      xfac(50,19) =     3.440153d0 !         Tin -    Potassium
      alpb(50,32) =     2.016055d0 !         Tin -    Germanium
      xfac(50,32) =     3.500376d0 !         Tin -    Germanium
      alpb(50,34) =     1.393411d0 !         Tin -     Selenium
      xfac(50,34) =     0.413851d0 !         Tin -     Selenium
      alpb(50,35) =     1.594297d0 !         Tin -      Bromine
      xfac(50,35) =     0.954605d0 !         Tin -      Bromine
      alpb(50,50) =     1.045406d0 !         Tin -          Tin
      xfac(50,50) =     0.300460d0 !         Tin -          Tin
 !
      alpb(51, 1) =     1.091035d0 !    Antimony -     Hydrogen
      xfac(51, 1) =     0.408876d0 !    Antimony -     Hydrogen
      alpb(51, 6) =     1.240714d0 !    Antimony -       Carbon
      xfac(51, 6) =     0.327493d0 !    Antimony -       Carbon
      alpb(51, 7) =     0.846645d0 !    Antimony -     Nitrogen
      xfac(51, 7) =     0.137604d0 !    Antimony -     Nitrogen
      alpb(51, 8) =     1.462059d0 !    Antimony -       Oxygen
      xfac(51, 8) =     0.346536d0 !    Antimony -       Oxygen
      alpb(51, 9) =     1.622505d0 !    Antimony -     Fluorine
      xfac(51, 9) =     0.283768d0 !    Antimony -     Fluorine
      alpb(51,11) =     1.106800d0 !    Antimony -       Sodium
      xfac(51,11) =     0.547287d0 !    Antimony -       Sodium
      alpb(51,13) =     1.085906d0 !    Antimony -     Aluminum
      xfac(51,13) =     1.291895d0 !    Antimony -     Aluminum
      alpb(51,14) =     2.519702d0 !    Antimony -      Silicon
      xfac(51,14) =     8.707039d0 !    Antimony -      Silicon
      alpb(51,16) =     1.016407d0 !    Antimony -       Sulfur
      xfac(51,16) =     0.211102d0 !    Antimony -       Sulfur
      alpb(51,17) =     1.170710d0 !    Antimony -     Chlorine
      xfac(51,17) =     0.217072d0 !    Antimony -     Chlorine
      alpb(51,25) =     1.698753d0 !    Antimony -    Manganese
      xfac(51,25) =     2.384408d0 !    Antimony -    Manganese
      alpb(51,27) =     2.204630d0 !    Antimony -       Cobalt
      xfac(51,27) =     2.276050d0 !    Antimony -       Cobalt
      alpb(51,35) =     1.227775d0 !    Antimony -      Bromine
      xfac(51,35) =     0.567204d0 !    Antimony -      Bromine
      alpb(51,43) =     2.204850d0 !    Antimony -   Technetium
      xfac(51,43) =     2.276260d0 !    Antimony -   Technetium
      alpb(51,44) =     2.968084d0 !    Antimony -    Ruthenium
      xfac(51,44) =     2.509269d0 !    Antimony -    Ruthenium
      alpb(51,45) =     2.204930d0 !    Antimony -      Rhodium
      xfac(51,45) =     2.276340d0 !    Antimony -      Rhodium
      alpb(51,49) =     1.011173d0 !    Antimony -       Indium
      xfac(51,49) =     0.470521d0 !    Antimony -       Indium
      alpb(51,51) =     0.657753d0 !    Antimony -     Antimony
      xfac(51,51) =     0.219843d0 !    Antimony -     Antimony
 !
      alpb(52, 1) =     2.879705d0 !   Tellurium -     Hydrogen
      xfac(52, 1) =     7.645321d0 !   Tellurium -     Hydrogen
      alpb(52, 5) =     2.443355d0 !   Tellurium -        Boron
      xfac(52, 5) =     2.926026d0 !   Tellurium -        Boron
      alpb(52, 6) =     2.858205d0 !   Tellurium -       Carbon
      xfac(52, 6) =     7.513380d0 !   Tellurium -       Carbon
      alpb(52, 7) =     2.548060d0 !   Tellurium -     Nitrogen
      xfac(52, 7) =     2.356842d0 !   Tellurium -     Nitrogen
      alpb(52, 8) =     2.359294d0 !   Tellurium -       Oxygen
      xfac(52, 8) =     1.147602d0 !   Tellurium -       Oxygen
      alpb(52, 9) =     3.109030d0 !   Tellurium -     Fluorine
      xfac(52, 9) =     2.199214d0 !   Tellurium -     Fluorine
      alpb(52,13) =     1.783994d0 !   Tellurium -     Aluminum
      xfac(52,13) =     9.305330d0 !   Tellurium -     Aluminum
      alpb(52,15) =     1.482343d0 !   Tellurium -   Phosphorus
      xfac(52,15) =     1.459960d0 !   Tellurium -   Phosphorus
      alpb(52,16) =     2.969323d0 !   Tellurium -       Sulfur
      xfac(52,16) =    14.279019d0 !   Tellurium -       Sulfur
      alpb(52,17) =     1.475730d0 !   Tellurium -     Chlorine
      xfac(52,17) =     0.514830d0 !   Tellurium -     Chlorine
      alpb(52,19) =     1.257635d0 !   Tellurium -    Potassium
      xfac(52,19) =     2.073198d0 !   Tellurium -    Potassium
      alpb(52,30) =     1.704782d0 !   Tellurium -         Zinc
      xfac(52,30) =     4.125260d0 !   Tellurium -         Zinc
      alpb(52,32) =     2.049526d0 !   Tellurium -    Germanium
      xfac(52,32) =     7.601044d0 !   Tellurium -    Germanium
      alpb(52,33) =     1.275249d0 !   Tellurium -      Arsenic
      xfac(52,33) =     0.866529d0 !   Tellurium -      Arsenic
      alpb(52,34) =     1.585819d0 !   Tellurium -     Selenium
      xfac(52,34) =     1.322800d0 !   Tellurium -     Selenium
      alpb(52,35) =     2.316655d0 !   Tellurium -      Bromine
      xfac(52,35) =     4.158560d0 !   Tellurium -      Bromine
      alpb(52,48) =     1.759718d0 !   Tellurium -      Cadmium
      xfac(52,48) =     8.405812d0 !   Tellurium -      Cadmium
      alpb(52,49) =     1.913212d0 !   Tellurium -       Indium
      xfac(52,49) =     9.943252d0 !   Tellurium -       Indium
      alpb(52,50) =     2.265433d0 !   Tellurium -          Tin
      xfac(52,50) =    11.004064d0 !   Tellurium -          Tin
      alpb(52,51) =     1.634994d0 !   Tellurium -     Antimony
      xfac(52,51) =     0.575666d0 !   Tellurium -     Antimony
      alpb(52,52) =     3.032862d0 !   Tellurium -    Tellurium
      xfac(52,52) =    29.604279d0 !   Tellurium -    Tellurium
 !
      alpb(53, 1) =     2.303917d0 !      Iodine -     Hydrogen
      xfac(53, 1) =     2.456367d0 !      Iodine -     Hydrogen
      alpb(53, 2) =     2.264096d0 !      Iodine -       Helium
      xfac(53, 2) =     2.613098d0 !      Iodine -       Helium
      alpb(53, 3) =     1.392191d0 !      Iodine -      Lithium
      xfac(53, 3) =     1.220335d0 !      Iodine -      Lithium
      alpb(53, 4) =     2.137694d0 !      Iodine -    Beryllium
      xfac(53, 4) =     4.012926d0 !      Iodine -    Beryllium
      alpb(53, 5) =     1.949150d0 !      Iodine -        Boron
      xfac(53, 5) =     1.926808d0 !      Iodine -        Boron
      alpb(53, 6) =     2.148076d0 !      Iodine -       Carbon
      xfac(53, 6) =     1.536813d0 !      Iodine -       Carbon
      alpb(53, 7) =     2.204300d0 !      Iodine -     Nitrogen
      xfac(53, 7) =     1.197247d0 !      Iodine -     Nitrogen
      alpb(53, 8) =     2.031236d0 !      Iodine -       Oxygen
      xfac(53, 8) =     0.673908d0 !      Iodine -       Oxygen
      alpb(53, 9) =     2.168508d0 !      Iodine -     Fluorine
      xfac(53, 9) =     0.518622d0 !      Iodine -     Fluorine
      alpb(53,10) =     2.572520d0 !      Iodine -         Neon
      xfac(53,10) =     1.449278d0 !      Iodine -         Neon
      alpb(53,11) =     1.999781d0 !      Iodine -       Sodium
      xfac(53,11) =    12.909796d0 !      Iodine -       Sodium
      alpb(53,12) =     1.832289d0 !      Iodine -    Magnesium
      xfac(53,12) =     4.415343d0 !      Iodine -    Magnesium
      alpb(53,13) =     1.515624d0 !      Iodine -     Aluminum
      xfac(53,13) =     2.691541d0 !      Iodine -     Aluminum
      alpb(53,14) =     1.472015d0 !      Iodine -      Silicon
      xfac(53,14) =     1.272495d0 !      Iodine -      Silicon
      alpb(53,15) =     1.560276d0 !      Iodine -   Phosphorus
      xfac(53,15) =     2.308251d0 !      Iodine -   Phosphorus
      alpb(53,16) =     2.108468d0 !      Iodine -       Sulfur
      xfac(53,16) =     1.287638d0 !      Iodine -       Sulfur
      alpb(53,17) =     1.674480d0 !      Iodine -     Chlorine
      xfac(53,17) =     0.582734d0 !      Iodine -     Chlorine
      alpb(53,18) =     1.583967d0 !      Iodine -        Argon
      xfac(53,18) =     0.297828d0 !      Iodine -        Argon
      alpb(53,19) =     1.527318d0 !      Iodine -    Potassium
      xfac(53,19) =     6.255639d0 !      Iodine -    Potassium
      alpb(53,20) =     1.931292d0 !      Iodine -      Calcium
      xfac(53,20) =     5.485613d0 !      Iodine -      Calcium
      alpb(53,21) =     1.888645d0 !      Iodine -     Scandium
      xfac(53,21) =     4.305507d0 !      Iodine -     Scandium
      alpb(53,22) =     1.569430d0 !      Iodine -     Titanium
      xfac(53,22) =     2.273746d0 !      Iodine -     Titanium
      alpb(53,23) =     1.204771d0 !      Iodine -     Vanadium
      xfac(53,23) =     0.566891d0 !      Iodine -     Vanadium
      alpb(53,24) =     1.505878d0 !      Iodine -     Chromium
      xfac(53,24) =     0.754833d0 !      Iodine -     Chromium
      alpb(53,25) =     1.920970d0 !      Iodine -    Manganese
      xfac(53,25) =     2.239969d0 !      Iodine -    Manganese
      alpb(53,26) =     1.995455d0 !      Iodine -         Iron
      xfac(53,26) =     1.244120d0 !      Iodine -         Iron
      alpb(53,27) =     2.394155d0 !      Iodine -       Cobalt
      xfac(53,27) =     3.145732d0 !      Iodine -       Cobalt
      alpb(53,28) =     2.491283d0 !      Iodine -       Nickel
      xfac(53,28) =     3.452112d0 !      Iodine -       Nickel
      alpb(53,29) =     3.049738d0 !      Iodine -       Copper
      xfac(53,29) =     5.342329d0 !      Iodine -       Copper
      alpb(53,30) =     1.785943d0 !      Iodine -         Zinc
      xfac(53,30) =     4.270507d0 !      Iodine -         Zinc
      alpb(53,31) =     1.903558d0 !      Iodine -      Gallium
      xfac(53,31) =     3.519264d0 !      Iodine -      Gallium
      alpb(53,32) =     1.431330d0 !      Iodine -    Germanium
      xfac(53,32) =     0.946363d0 !      Iodine -    Germanium
      alpb(53,33) =     1.454624d0 !      Iodine -      Arsenic
      xfac(53,33) =     0.863506d0 !      Iodine -      Arsenic
      alpb(53,34) =     1.464103d0 !      Iodine -     Selenium
      xfac(53,34) =     0.509254d0 !      Iodine -     Selenium
      alpb(53,35) =     1.793757d0 !      Iodine -      Bromine
      xfac(53,35) =     1.192163d0 !      Iodine -      Bromine
      alpb(53,36) =     1.242469d0 !      Iodine -      Krypton
      xfac(53,36) =     0.195416d0 !      Iodine -      Krypton
      alpb(53,37) =     0.893509d0 !      Iodine -     Rubidium
      xfac(53,37) =     0.753057d0 !      Iodine -     Rubidium
      alpb(53,38) =     2.702289d0 !      Iodine -    Strontium
      xfac(53,38) =    32.561240d0 !      Iodine -    Strontium
      alpb(53,39) =     1.443236d0 !      Iodine -      Yttrium
      xfac(53,39) =     2.307839d0 !      Iodine -      Yttrium
      alpb(53,40) =     1.402802d0 !      Iodine -    Zirconium
      xfac(53,40) =     1.833851d0 !      Iodine -    Zirconium
      alpb(53,41) =     2.001333d0 !      Iodine -      Niobium
      xfac(53,41) =     4.678302d0 !      Iodine -      Niobium
      alpb(53,42) =     2.042051d0 !      Iodine -   Molybdenum
      xfac(53,42) =     3.618240d0 !      Iodine -   Molybdenum
      alpb(53,43) =     2.576693d0 !      Iodine -   Technetium
      xfac(53,43) =     9.860653d0 !      Iodine -   Technetium
      alpb(53,44) =     1.432008d0 !      Iodine -    Ruthenium
      xfac(53,44) =     0.552218d0 !      Iodine -    Ruthenium
      alpb(53,45) =     2.347687d0 !      Iodine -      Rhodium
      xfac(53,45) =     1.588054d0 !      Iodine -      Rhodium
      alpb(53,46) =     1.720521d0 !      Iodine -    Palladium
      xfac(53,46) =     0.587118d0 !      Iodine -    Palladium
      alpb(53,47) =     2.959757d0 !      Iodine -       Silver
      xfac(53,47) =     9.538157d0 !      Iodine -       Silver
      alpb(53,48) =     1.751947d0 !      Iodine -      Cadmium
      xfac(53,48) =     6.820820d0 !      Iodine -      Cadmium
      alpb(53,49) =     1.830626d0 !      Iodine -       Indium
      xfac(53,49) =     4.302750d0 !      Iodine -       Indium
      alpb(53,50) =     2.479003d0 !      Iodine -          Tin
      xfac(53,50) =    24.450811d0 !      Iodine -          Tin
      alpb(53,51) =     1.114193d0 !      Iodine -     Antimony
      xfac(53,51) =     0.767547d0 !      Iodine -     Antimony
      alpb(53,52) =     2.102109d0 !      Iodine -    Tellurium
      xfac(53,52) =     4.751442d0 !      Iodine -    Tellurium
      alpb(53,53) =     1.619225d0 !      Iodine -       Iodine
      xfac(53,53) =     1.278518d0 !      Iodine -       Iodine
 !
      alpb(54, 1) =     1.356861d0 !       Xenon -     Hydrogen
      xfac(54, 1) =     0.701016d0 !       Xenon -     Hydrogen
      alpb(54, 2) =     2.497832d0 !       Xenon -       Helium
      xfac(54, 2) =     2.599471d0 !       Xenon -       Helium
      alpb(54, 3) =     1.697716d0 !       Xenon -      Lithium
      xfac(54, 3) =     4.467048d0 !       Xenon -      Lithium
      alpb(54, 4) =     6.000011d0 !       Xenon -    Beryllium
      xfac(54, 4) =     0.654334d0 !       Xenon -    Beryllium
      alpb(54, 5) =     3.233962d0 !       Xenon -        Boron
      xfac(54, 5) =     1.995594d0 !       Xenon -        Boron
      alpb(54, 6) =     1.704440d0 !       Xenon -       Carbon
      xfac(54, 6) =     0.826727d0 !       Xenon -       Carbon
      alpb(54, 7) =     1.932952d0 !       Xenon -     Nitrogen
      xfac(54, 7) =     0.925624d0 !       Xenon -     Nitrogen
      alpb(54, 8) =     2.566313d0 !       Xenon -       Oxygen
      xfac(54, 8) =     1.623526d0 !       Xenon -       Oxygen
      alpb(54, 9) =     2.837749d0 !       Xenon -     Fluorine
      xfac(54, 9) =     2.086480d0 !       Xenon -     Fluorine
      alpb(54,10) =     1.330202d0 !       Xenon -         Neon
      xfac(54,10) =     0.293862d0 !       Xenon -         Neon
      alpb(54,11) =     1.291138d0 !       Xenon -       Sodium
      xfac(54,11) =     5.076100d0 !       Xenon -       Sodium
      alpb(54,12) =     2.756089d0 !       Xenon -    Magnesium
      xfac(54,12) =     9.774960d0 !       Xenon -    Magnesium
      alpb(54,13) =     2.420691d0 !       Xenon -     Aluminum
      xfac(54,13) =     7.358944d0 !       Xenon -     Aluminum
      alpb(54,14) =     2.796986d0 !       Xenon -      Silicon
      xfac(54,14) =    16.526889d0 !       Xenon -      Silicon
      alpb(54,17) =     1.389615d0 !       Xenon -     Chlorine
      xfac(54,17) =     0.593028d0 !       Xenon -     Chlorine
      alpb(54,18) =     0.591520d0 !       Xenon -        Argon
      xfac(54,18) =     0.049266d0 !       Xenon -        Argon
      alpb(54,19) =     0.886811d0 !       Xenon -    Potassium
      xfac(54,19) =     1.526138d0 !       Xenon -    Potassium
      alpb(54,20) =     1.698890d0 !       Xenon -      Calcium
      xfac(54,20) =     2.050654d0 !       Xenon -      Calcium
      alpb(54,35) =     1.400900d0 !       Xenon -      Bromine
      xfac(54,35) =     0.711370d0 !       Xenon -      Bromine
      alpb(54,36) =     0.551561d0 !       Xenon -      Krypton
      xfac(54,36) =     0.049793d0 !       Xenon -      Krypton
      alpb(54,37) =     1.345397d0 !       Xenon -     Rubidium
      xfac(54,37) =     1.856289d0 !       Xenon -     Rubidium
      alpb(54,53) =     1.187975d0 !       Xenon -       Iodine
      xfac(54,53) =     0.555791d0 !       Xenon -       Iodine
      alpb(54,54) =     1.912510d0 !       Xenon -        Xenon
      xfac(54,54) =     9.565337d0 !       Xenon -        Xenon
 !
      alpb(55, 1) =     1.719572d0 !      Cesium -     Hydrogen
      xfac(55, 1) =     2.711386d0 !      Cesium -     Hydrogen
      alpb(55, 5) =     3.000034d0 !      Cesium -        Boron
      xfac(55, 5) =    10.289233d0 !      Cesium -        Boron
      alpb(55, 6) =     2.251416d0 !      Cesium -       Carbon
      xfac(55, 6) =    17.858749d0 !      Cesium -       Carbon
      alpb(55, 7) =     2.465681d0 !      Cesium -     Nitrogen
      xfac(55, 7) =    28.270100d0 !      Cesium -     Nitrogen
      alpb(55, 8) =     1.517551d0 !      Cesium -       Oxygen
      xfac(55, 8) =     0.871027d0 !      Cesium -       Oxygen
      alpb(55, 9) =     1.636155d0 !      Cesium -     Fluorine
      xfac(55, 9) =     0.551707d0 !      Cesium -     Fluorine
      alpb(55,15) =     0.917812d0 !      Cesium -   Phosphorus
      xfac(55,15) =     0.499881d0 !      Cesium -   Phosphorus
      alpb(55,16) =     1.348833d0 !      Cesium -       Sulfur
      xfac(55,16) =     1.767711d0 !      Cesium -       Sulfur
      alpb(55,17) =     1.241351d0 !      Cesium -     Chlorine
      xfac(55,17) =     0.942491d0 !      Cesium -     Chlorine
      alpb(55,19) =     1.722882d0 !      Cesium -    Potassium
      xfac(55,19) =     1.188658d0 !      Cesium -    Potassium
      alpb(55,23) =     2.002665d0 !      Cesium -     Vanadium
      xfac(55,23) =    11.159719d0 !      Cesium -     Vanadium
      alpb(55,35) =     1.949820d0 !      Cesium -      Bromine
      xfac(55,35) =    13.999636d0 !      Cesium -      Bromine
      alpb(55,39) =     0.929803d0 !      Cesium -      Yttrium
      xfac(55,39) =     1.057035d0 !      Cesium -      Yttrium
      alpb(55,49) =     0.852581d0 !      Cesium -       Indium
      xfac(55,49) =     4.457697d0 !      Cesium -       Indium
      alpb(55,53) =     1.239277d0 !      Cesium -       Iodine
      xfac(55,53) =     3.226708d0 !      Cesium -       Iodine
      alpb(55,55) =     1.267283d0 !      Cesium -       Cesium
      xfac(55,55) =    29.382256d0 !      Cesium -       Cesium
 !
      alpb(56, 1) =     3.120384d0 !      Barium -     Hydrogen
      xfac(56, 1) =    27.058819d0 !      Barium -     Hydrogen
      alpb(56, 6) =     1.318794d0 !      Barium -       Carbon
      xfac(56, 6) =     0.549254d0 !      Barium -       Carbon
      alpb(56, 7) =     2.188957d0 !      Barium -     Nitrogen
      xfac(56, 7) =     4.679835d0 !      Barium -     Nitrogen
      alpb(56, 8) =     2.337452d0 !      Barium -       Oxygen
      xfac(56, 8) =     4.174798d0 !      Barium -       Oxygen
      alpb(56, 9) =     2.539909d0 !      Barium -     Fluorine
      xfac(56, 9) =     3.008132d0 !      Barium -     Fluorine
      alpb(56,12) =     1.432600d0 !      Barium -    Magnesium
      xfac(56,12) =    10.342497d0 !      Barium -    Magnesium
      alpb(56,13) =     2.891358d0 !      Barium -     Aluminum
      xfac(56,13) =    15.460538d0 !      Barium -     Aluminum
      alpb(56,14) =     0.996995d0 !      Barium -      Silicon
      xfac(56,14) =     0.887820d0 !      Barium -      Silicon
      alpb(56,15) =     1.646819d0 !      Barium -   Phosphorus
      xfac(56,15) =     8.719637d0 !      Barium -   Phosphorus
      alpb(56,16) =     1.637082d0 !      Barium -       Sulfur
      xfac(56,16) =     1.742576d0 !      Barium -       Sulfur
      alpb(56,17) =     1.987384d0 !      Barium -     Chlorine
      xfac(56,17) =     2.636334d0 !      Barium -     Chlorine
      alpb(56,20) =     1.342035d0 !      Barium -      Calcium
      xfac(56,20) =     2.833561d0 !      Barium -      Calcium
      alpb(56,22) =     1.702345d0 !      Barium -     Titanium
      xfac(56,22) =     4.943061d0 !      Barium -     Titanium
      alpb(56,29) =     1.699850d0 !      Barium -       Copper
      xfac(56,29) =     1.896329d0 !      Barium -       Copper
      alpb(56,35) =     1.806723d0 !      Barium -      Bromine
      xfac(56,35) =     2.830984d0 !      Barium -      Bromine
      alpb(56,51) =     1.329425d0 !      Barium -     Antimony
      xfac(56,51) =    12.262981d0 !      Barium -     Antimony
      alpb(56,53) =     1.370665d0 !      Barium -       Iodine
      xfac(56,53) =     2.112756d0 !      Barium -       Iodine
      alpb(56,56) =     1.860576d0 !      Barium -       Barium
      xfac(56,56) =    57.199345d0 !      Barium -       Barium
 !
      alpb(57, 1) =     1.073406d0 !   Lanthanum -     Hydrogen
      xfac(57, 1) =     0.399521d0 !   Lanthanum -     Hydrogen
      alpb(57, 6) =     2.129683d0 !   Lanthanum -       Carbon
      xfac(57, 6) =     4.650201d0 !   Lanthanum -       Carbon
      alpb(57, 7) =     2.329214d0 !   Lanthanum -     Nitrogen
      xfac(57, 7) =     2.192625d0 !   Lanthanum -     Nitrogen
      alpb(57, 8) =     1.940554d0 !   Lanthanum -       Oxygen
      xfac(57, 8) =     1.648001d0 !   Lanthanum -       Oxygen
      alpb(57, 9) =     2.228378d0 !   Lanthanum -     Fluorine
      xfac(57, 9) =     1.892928d0 !   Lanthanum -     Fluorine
      alpb(57,13) =     1.003510d0 !   Lanthanum -     Aluminum
      xfac(57,13) =     0.500540d0 !   Lanthanum -     Aluminum
      alpb(57,14) =     2.872867d0 !   Lanthanum -      Silicon
      xfac(57,14) =     1.218295d0 !   Lanthanum -      Silicon
      alpb(57,15) =     1.991054d0 !   Lanthanum -   Phosphorus
      xfac(57,15) =    18.284518d0 !   Lanthanum -   Phosphorus
      alpb(57,16) =     1.158196d0 !   Lanthanum -       Sulfur
      xfac(57,16) =     0.486832d0 !   Lanthanum -       Sulfur
      alpb(57,17) =     1.835651d0 !   Lanthanum -     Chlorine
      xfac(57,17) =     1.631876d0 !   Lanthanum -     Chlorine
      alpb(57,35) =     1.253581d0 !   Lanthanum -      Bromine
      xfac(57,35) =     0.731795d0 !   Lanthanum -      Bromine
      alpb(57,53) =     1.612519d0 !   Lanthanum -       Iodine
      xfac(57,53) =     3.278712d0 !   Lanthanum -       Iodine
      alpb(57,57) =     2.066209d0 !   Lanthanum -    Lanthanum
      xfac(57,57) =    29.272376d0 !   Lanthanum -    Lanthanum
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
      alpb(64,14) =     2.112525d0 !  Gadolinium -      Silicon
      xfac(64,14) =     3.203995d0 !  Gadolinium -      Silicon
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
      alpb(71, 1) =     2.089118d0 !    Lutetium -     Hydrogen
      xfac(71, 1) =     7.421490d0 !    Lutetium -     Hydrogen
      alpb(71, 6) =     2.386830d0 !    Lutetium -       Carbon
      xfac(71, 6) =     6.432873d0 !    Lutetium -       Carbon
      alpb(71, 7) =     1.647895d0 !    Lutetium -     Nitrogen
      xfac(71, 7) =     0.783123d0 !    Lutetium -     Nitrogen
      alpb(71, 8) =     1.889190d0 !    Lutetium -       Oxygen
      xfac(71, 8) =     0.868896d0 !    Lutetium -       Oxygen
      alpb(71, 9) =     1.888274d0 !    Lutetium -     Fluorine
      xfac(71, 9) =     1.730185d0 !    Lutetium -     Fluorine
      alpb(71,15) =     1.345992d0 !    Lutetium -   Phosphorus
      xfac(71,15) =     8.048165d0 !    Lutetium -   Phosphorus
      alpb(71,17) =     2.558367d0 !    Lutetium -     Chlorine
      xfac(71,17) =     8.330639d0 !    Lutetium -     Chlorine
      alpb(71,35) =     1.381701d0 !    Lutetium -      Bromine
      xfac(71,35) =     0.992835d0 !    Lutetium -      Bromine
      alpb(71,53) =     1.436788d0 !    Lutetium -       Iodine
      xfac(71,53) =     4.313665d0 !    Lutetium -       Iodine
      alpb(71,71) =     1.403653d0 !    Lutetium -     Lutetium
      xfac(71,71) =    21.889048d0 !    Lutetium -     Lutetium
 !
      alpb(72, 1) =     2.088799d0 !     Hafnium -     Hydrogen
      xfac(72, 1) =     3.833288d0 !     Hafnium -     Hydrogen
      alpb(72, 5) =     1.617370d0 !     Hafnium -        Boron
      xfac(72, 5) =     0.588837d0 !     Hafnium -        Boron
      alpb(72, 6) =     2.294622d0 !     Hafnium -       Carbon
      xfac(72, 6) =     4.159075d0 !     Hafnium -       Carbon
      alpb(72, 7) =     2.521801d0 !     Hafnium -     Nitrogen
      xfac(72, 7) =     5.468404d0 !     Hafnium -     Nitrogen
      alpb(72, 8) =     2.446232d0 !     Hafnium -       Oxygen
      xfac(72, 8) =     2.857484d0 !     Hafnium -       Oxygen
      alpb(72, 9) =     2.979096d0 !     Hafnium -     Fluorine
      xfac(72, 9) =     4.736067d0 !     Hafnium -     Fluorine
      alpb(72,11) =     1.840619d0 !     Hafnium -       Sodium
      xfac(72,11) =     8.832085d0 !     Hafnium -       Sodium
      alpb(72,12) =     1.911350d0 !     Hafnium -    Magnesium
      xfac(72,12) =     4.330250d0 !     Hafnium -    Magnesium
      alpb(72,13) =     0.949150d0 !     Hafnium -     Aluminum
      xfac(72,13) =     0.622520d0 !     Hafnium -     Aluminum
      alpb(72,14) =     2.189300d0 !     Hafnium -      Silicon
      xfac(72,14) =     3.382300d0 !     Hafnium -      Silicon
      alpb(72,15) =     2.099591d0 !     Hafnium -   Phosphorus
      xfac(72,15) =     5.936976d0 !     Hafnium -   Phosphorus
      alpb(72,16) =     2.327110d0 !     Hafnium -       Sulfur
      xfac(72,16) =     1.666760d0 !     Hafnium -       Sulfur
      alpb(72,17) =     1.953166d0 !     Hafnium -     Chlorine
      xfac(72,17) =     1.685929d0 !     Hafnium -     Chlorine
      alpb(72,20) =     2.054500d0 !     Hafnium -      Calcium
      xfac(72,20) =     4.319510d0 !     Hafnium -      Calcium
      alpb(72,33) =     1.799500d0 !     Hafnium -      Arsenic
      xfac(72,33) =     1.280820d0 !     Hafnium -      Arsenic
      alpb(72,35) =     2.237896d0 !     Hafnium -      Bromine
      xfac(72,35) =     6.312154d0 !     Hafnium -      Bromine
      alpb(72,53) =     2.354639d0 !     Hafnium -       Iodine
      xfac(72,53) =    18.443532d0 !     Hafnium -       Iodine
      alpb(72,56) =     2.264830d0 !     Hafnium -       Barium
      xfac(72,56) =     9.022520d0 !     Hafnium -       Barium
      alpb(72,72) =     2.216588d0 !     Hafnium -      Hafnium
      xfac(72,72) =    29.394192d0 !     Hafnium -      Hafnium
 !
      alpb(73, 1) =     1.786631d0 !    Tantalum -     Hydrogen
      xfac(73, 1) =     1.893110d0 !    Tantalum -     Hydrogen
      alpb(73, 6) =     1.450720d0 !    Tantalum -       Carbon
      xfac(73, 6) =     0.581370d0 !    Tantalum -       Carbon
      alpb(73, 7) =     2.013737d0 !    Tantalum -     Nitrogen
      xfac(73, 7) =     1.152896d0 !    Tantalum -     Nitrogen
      alpb(73, 8) =     2.494885d0 !    Tantalum -       Oxygen
      xfac(73, 8) =     2.316225d0 !    Tantalum -       Oxygen
      alpb(73, 9) =     2.732769d0 !    Tantalum -     Fluorine
      xfac(73, 9) =     2.163084d0 !    Tantalum -     Fluorine
      alpb(73,11) =     2.551120d0 !    Tantalum -       Sodium
      xfac(73,11) =     8.276130d0 !    Tantalum -       Sodium
      alpb(73,15) =     2.513800d0 !    Tantalum -   Phosphorus
      xfac(73,15) =     6.261880d0 !    Tantalum -   Phosphorus
      alpb(73,16) =     2.091335d0 !    Tantalum -       Sulfur
      xfac(73,16) =     3.201126d0 !    Tantalum -       Sulfur
      alpb(73,17) =     2.003584d0 !    Tantalum -     Chlorine
      xfac(73,17) =     1.661719d0 !    Tantalum -     Chlorine
      alpb(73,19) =     4.521470d0 !    Tantalum -    Potassium
      xfac(73,19) =     2.026700d0 !    Tantalum -    Potassium
      alpb(73,35) =     1.962327d0 !    Tantalum -      Bromine
      xfac(73,35) =     3.310460d0 !    Tantalum -      Bromine
      alpb(73,53) =     1.500797d0 !    Tantalum -       Iodine
      xfac(73,53) =     1.995370d0 !    Tantalum -       Iodine
      alpb(73,73) =     0.982767d0 !    Tantalum -     Tantalum
      xfac(73,73) =     0.831956d0 !    Tantalum -     Tantalum
 !
      alpb(74, 1) =     2.665390d0 !    Tungsten -     Hydrogen
      xfac(74, 1) =     5.441909d0 !    Tungsten -     Hydrogen
      alpb(74, 6) =     2.600118d0 !    Tungsten -       Carbon
      xfac(74, 6) =     4.729842d0 !    Tungsten -       Carbon
      alpb(74, 7) =     2.505903d0 !    Tungsten -     Nitrogen
      xfac(74, 7) =     3.877751d0 !    Tungsten -     Nitrogen
      alpb(74, 8) =     2.343168d0 !    Tungsten -       Oxygen
      xfac(74, 8) =     1.878859d0 !    Tungsten -       Oxygen
      alpb(74, 9) =     2.411756d0 !    Tungsten -     Fluorine
      xfac(74, 9) =     1.368205d0 !    Tungsten -     Fluorine
      alpb(74,11) =     1.090156d0 !    Tungsten -       Sodium
      xfac(74,11) =     0.686226d0 !    Tungsten -       Sodium
      alpb(74,12) =     1.434249d0 !    Tungsten -    Magnesium
      xfac(74,12) =     1.904971d0 !    Tungsten -    Magnesium
      alpb(74,15) =     1.715627d0 !    Tungsten -   Phosphorus
      xfac(74,15) =     4.472129d0 !    Tungsten -   Phosphorus
      alpb(74,16) =     2.045564d0 !    Tungsten -       Sulfur
      xfac(74,16) =     2.401567d0 !    Tungsten -       Sulfur
      alpb(74,17) =     1.907817d0 !    Tungsten -     Chlorine
      xfac(74,17) =     1.349560d0 !    Tungsten -     Chlorine
      alpb(74,19) =     1.521243d0 !    Tungsten -    Potassium
      xfac(74,19) =     2.096182d0 !    Tungsten -    Potassium
      alpb(74,20) =     1.870733d0 !    Tungsten -      Calcium
      xfac(74,20) =     8.590544d0 !    Tungsten -      Calcium
      alpb(74,26) =     1.787925d0 !    Tungsten -         Iron
      xfac(74,26) =     1.977394d0 !    Tungsten -         Iron
      alpb(74,28) =     1.775099d0 !    Tungsten -       Nickel
      xfac(74,28) =     1.430746d0 !    Tungsten -       Nickel
      alpb(74,30) =     1.928464d0 !    Tungsten -         Zinc
      xfac(74,30) =     5.376323d0 !    Tungsten -         Zinc
      alpb(74,35) =     2.143627d0 !    Tungsten -      Bromine
      xfac(74,35) =     3.993357d0 !    Tungsten -      Bromine
      alpb(74,37) =     0.900113d0 !    Tungsten -     Rubidium
      xfac(74,37) =     4.075269d0 !    Tungsten -     Rubidium
      alpb(74,40) =     2.023641d0 !    Tungsten -    Zirconium
      xfac(74,40) =    19.994079d0 !    Tungsten -    Zirconium
      alpb(74,53) =     1.997307d0 !    Tungsten -       Iodine
      xfac(74,53) =     5.825642d0 !    Tungsten -       Iodine
      alpb(74,55) =     0.899625d0 !    Tungsten -       Cesium
      xfac(74,55) =     4.061044d0 !    Tungsten -       Cesium
      alpb(74,56) =     1.566159d0 !    Tungsten -       Barium
      xfac(74,56) =     1.861828d0 !    Tungsten -       Barium
      alpb(74,74) =     2.141401d0 !    Tungsten -     Tungsten
      xfac(74,74) =    13.807246d0 !    Tungsten -     Tungsten
 !
      alpb(75, 1) =     1.748317d0 !     Rhenium -     Hydrogen
      xfac(75, 1) =     0.497281d0 !     Rhenium -     Hydrogen
      alpb(75, 6) =     2.109510d0 !     Rhenium -       Carbon
      xfac(75, 6) =     0.646616d0 !     Rhenium -       Carbon
      alpb(75, 7) =     2.474230d0 !     Rhenium -     Nitrogen
      xfac(75, 7) =     1.438570d0 !     Rhenium -     Nitrogen
      alpb(75, 8) =     2.403640d0 !     Rhenium -       Oxygen
      xfac(75, 8) =     1.075351d0 !     Rhenium -       Oxygen
      alpb(75, 9) =     2.790322d0 !     Rhenium -     Fluorine
      xfac(75, 9) =     1.384171d0 !     Rhenium -     Fluorine
      alpb(75,14) =     2.775930d0 !     Rhenium -      Silicon
      xfac(75,14) =     0.849450d0 !     Rhenium -      Silicon
      alpb(75,15) =     1.316878d0 !     Rhenium -   Phosphorus
      xfac(75,15) =     0.761808d0 !     Rhenium -   Phosphorus
      alpb(75,16) =     2.637193d0 !     Rhenium -       Sulfur
      xfac(75,16) =     3.055234d0 !     Rhenium -       Sulfur
      alpb(75,17) =     2.857608d0 !     Rhenium -     Chlorine
      xfac(75,17) =     3.265852d0 !     Rhenium -     Chlorine
      alpb(75,32) =     2.852340d0 !     Rhenium -    Germanium
      xfac(75,32) =     2.151580d0 !     Rhenium -    Germanium
      alpb(75,34) =     2.523170d0 !     Rhenium -     Selenium
      xfac(75,34) =     2.202140d0 !     Rhenium -     Selenium
      alpb(75,35) =     2.195052d0 !     Rhenium -      Bromine
      xfac(75,35) =     1.575571d0 !     Rhenium -      Bromine
      alpb(75,51) =     2.204360d0 !     Rhenium -     Antimony
      xfac(75,51) =     2.275780d0 !     Rhenium -     Antimony
      alpb(75,53) =     2.239594d0 !     Rhenium -       Iodine
      xfac(75,53) =     3.240592d0 !     Rhenium -       Iodine
      alpb(75,75) =     2.195649d0 !     Rhenium -      Rhenium
      xfac(75,75) =     1.776660d0 !     Rhenium -      Rhenium
 !
      alpb(76, 1) =     2.399448d0 !      Osmium -     Hydrogen
      xfac(76, 1) =     3.609773d0 !      Osmium -     Hydrogen
      alpb(76, 6) =     1.938959d0 !      Osmium -       Carbon
      xfac(76, 6) =     0.616916d0 !      Osmium -       Carbon
      alpb(76, 7) =     2.139750d0 !      Osmium -     Nitrogen
      xfac(76, 7) =     0.730399d0 !      Osmium -     Nitrogen
      alpb(76, 8) =     2.539022d0 !      Osmium -       Oxygen
      xfac(76, 8) =     1.230187d0 !      Osmium -       Oxygen
      alpb(76, 9) =     2.210417d0 !      Osmium -     Fluorine
      xfac(76, 9) =     0.562952d0 !      Osmium -     Fluorine
      alpb(76,11) =     2.550740d0 !      Osmium -       Sodium
      xfac(76,11) =     8.275750d0 !      Osmium -       Sodium
      alpb(76,15) =     2.060122d0 !      Osmium -   Phosphorus
      xfac(76,15) =     4.267629d0 !      Osmium -   Phosphorus
      alpb(76,16) =     2.809500d0 !      Osmium -       Sulfur
      xfac(76,16) =     4.186050d0 !      Osmium -       Sulfur
      alpb(76,17) =     2.080978d0 !      Osmium -     Chlorine
      xfac(76,17) =     1.177666d0 !      Osmium -     Chlorine
      alpb(76,19) =     1.351484d0 !      Osmium -    Potassium
      xfac(76,19) =     0.875486d0 !      Osmium -    Potassium
      alpb(76,35) =     2.225810d0 !      Osmium -      Bromine
      xfac(76,35) =     2.709104d0 !      Osmium -      Bromine
      alpb(76,53) =     2.189487d0 !      Osmium -       Iodine
      xfac(76,53) =     4.869377d0 !      Osmium -       Iodine
      alpb(76,76) =     1.661052d0 !      Osmium -       Osmium
      xfac(76,76) =     0.928334d0 !      Osmium -       Osmium
 !
      alpb(77, 1) =     1.634365d0 !     Iridium -     Hydrogen
      xfac(77, 1) =     0.406470d0 !     Iridium -     Hydrogen
      alpb(77, 6) =     1.604977d0 !     Iridium -       Carbon
      xfac(77, 6) =     0.185955d0 !     Iridium -       Carbon
      alpb(77, 7) =     2.997358d0 !     Iridium -     Nitrogen
      xfac(77, 7) =     1.790021d0 !     Iridium -     Nitrogen
      alpb(77, 8) =     3.116069d0 !     Iridium -       Oxygen
      xfac(77, 8) =     2.303902d0 !     Iridium -       Oxygen
      alpb(77, 9) =     2.612609d0 !     Iridium -     Fluorine
      xfac(77, 9) =     0.714245d0 !     Iridium -     Fluorine
      alpb(77,11) =     2.550820d0 !     Iridium -       Sodium
      xfac(77,11) =     8.275830d0 !     Iridium -       Sodium
      alpb(77,15) =     2.714060d0 !     Iridium -   Phosphorus
      xfac(77,15) =     6.284670d0 !     Iridium -   Phosphorus
      alpb(77,16) =     3.009199d0 !     Iridium -       Sulfur
      xfac(77,16) =     2.680449d0 !     Iridium -       Sulfur
      alpb(77,17) =     2.575683d0 !     Iridium -     Chlorine
      xfac(77,17) =     0.858848d0 !     Iridium -     Chlorine
      alpb(77,19) =     4.521170d0 !     Iridium -    Potassium
      xfac(77,19) =     2.026400d0 !     Iridium -    Potassium
      alpb(77,35) =     2.058351d0 !     Iridium -      Bromine
      xfac(77,35) =     0.804901d0 !     Iridium -      Bromine
      alpb(77,53) =     2.031222d0 !     Iridium -       Iodine
      xfac(77,53) =     1.787121d0 !     Iridium -       Iodine
      alpb(77,55) =     1.559526d0 !     Iridium -       Cesium
      xfac(77,55) =     1.027369d0 !     Iridium -       Cesium
      alpb(77,77) =     1.465795d0 !     Iridium -      Iridium
      xfac(77,77) =     0.190914d0 !     Iridium -      Iridium
 !
      alpb(78, 1) =     3.062604d0 !    Platinum -     Hydrogen
      xfac(78, 1) =     2.051954d0 !    Platinum -     Hydrogen
      alpb(78, 6) =     2.296772d0 !    Platinum -       Carbon
      xfac(78, 6) =     0.370388d0 !    Platinum -       Carbon
      alpb(78, 7) =     2.347134d0 !    Platinum -     Nitrogen
      xfac(78, 7) =     0.447775d0 !    Platinum -     Nitrogen
      alpb(78, 8) =     2.680367d0 !    Platinum -       Oxygen
      xfac(78, 8) =     0.827827d0 !    Platinum -       Oxygen
      alpb(78, 9) =     3.157007d0 !    Platinum -     Fluorine
      xfac(78, 9) =     1.031240d0 !    Platinum -     Fluorine
      alpb(78,13) =     1.572360d0 !    Platinum -     Aluminum
      xfac(78,13) =     1.056930d0 !    Platinum -     Aluminum
      alpb(78,14) =     0.999990d0 !    Platinum -      Silicon
      xfac(78,14) =     0.099990d0 !    Platinum -      Silicon
      alpb(78,15) =     1.307810d0 !    Platinum -   Phosphorus
      xfac(78,15) =     0.485582d0 !    Platinum -   Phosphorus
      alpb(78,16) =     2.919597d0 !    Platinum -       Sulfur
      xfac(78,16) =     2.008326d0 !    Platinum -       Sulfur
      alpb(78,17) =     3.034813d0 !    Platinum -     Chlorine
      xfac(78,17) =     1.610994d0 !    Platinum -     Chlorine
      alpb(78,19) =     1.495407d0 !    Platinum -    Potassium
      xfac(78,19) =     2.058817d0 !    Platinum -    Potassium
      alpb(78,35) =     2.596546d0 !    Platinum -      Bromine
      xfac(78,35) =     1.409311d0 !    Platinum -      Bromine
      alpb(78,47) =     1.387422d0 !    Platinum -       Silver
      xfac(78,47) =     5.456551d0 !    Platinum -       Silver
      alpb(78,53) =     2.228284d0 !    Platinum -       Iodine
      xfac(78,53) =     1.174520d0 !    Platinum -       Iodine
      alpb(78,78) =     3.276872d0 !    Platinum -     Platinum
      xfac(78,78) =     8.178033d0 !    Platinum -     Platinum
 !
      alpb(79, 1) =     2.006469d0 !        Gold -     Hydrogen
      xfac(79, 1) =     0.748516d0 !        Gold -     Hydrogen
      alpb(79, 6) =     2.119485d0 !        Gold -       Carbon
      xfac(79, 6) =     0.603200d0 !        Gold -       Carbon
      alpb(79, 7) =     2.395362d0 !        Gold -     Nitrogen
      xfac(79, 7) =     0.620935d0 !        Gold -     Nitrogen
      alpb(79, 8) =     2.323131d0 !        Gold -       Oxygen
      xfac(79, 8) =     0.355344d0 !        Gold -       Oxygen
      alpb(79, 9) =     3.153884d0 !        Gold -     Fluorine
      xfac(79, 9) =     0.880186d0 !        Gold -     Fluorine
      alpb(79,13) =     1.572570d0 !        Gold -     Aluminum
      xfac(79,13) =     1.057140d0 !        Gold -     Aluminum
      alpb(79,15) =     1.360881d0 !        Gold -   Phosphorus
      xfac(79,15) =     0.477023d0 !        Gold -   Phosphorus
      alpb(79,16) =     1.908644d0 !        Gold -       Sulfur
      xfac(79,16) =     0.423252d0 !        Gold -       Sulfur
      alpb(79,17) =     2.495913d0 !        Gold -     Chlorine
      xfac(79,17) =     1.153029d0 !        Gold -     Chlorine
      alpb(79,19) =     1.098797d0 !        Gold -    Potassium
      xfac(79,19) =     0.777289d0 !        Gold -    Potassium
      alpb(79,34) =     1.840962d0 !        Gold -     Selenium
      xfac(79,34) =     1.396436d0 !        Gold -     Selenium
      alpb(79,35) =     1.633736d0 !        Gold -      Bromine
      xfac(79,35) =     0.259683d0 !        Gold -      Bromine
      alpb(79,46) =     1.311827d0 !        Gold -    Palladium
      xfac(79,46) =     0.663311d0 !        Gold -    Palladium
      alpb(79,53) =     2.017017d0 !        Gold -       Iodine
      xfac(79,53) =     1.635086d0 !        Gold -       Iodine
      alpb(79,79) =     1.539843d0 !        Gold -         Gold
      xfac(79,79) =     0.352184d0 !        Gold -         Gold
 !
      alpb(80, 1) =     1.953060d0 !     Mercury -     Hydrogen
      xfac(80, 1) =     3.306359d0 !     Mercury -     Hydrogen
      alpb(80, 6) =     1.702331d0 !     Mercury -       Carbon
      xfac(80, 6) =     0.911944d0 !     Mercury -       Carbon
      alpb(80, 7) =     1.715039d0 !     Mercury -     Nitrogen
      xfac(80, 7) =     1.016140d0 !     Mercury -     Nitrogen
      alpb(80, 8) =     2.151298d0 !     Mercury -       Oxygen
      xfac(80, 8) =     2.032727d0 !     Mercury -       Oxygen
      alpb(80, 9) =     1.836494d0 !     Mercury -     Fluorine
      xfac(80, 9) =     0.631905d0 !     Mercury -     Fluorine
      alpb(80,11) =     1.459803d0 !     Mercury -       Sodium
      xfac(80,11) =     2.437893d0 !     Mercury -       Sodium
      alpb(80,14) =     2.770860d0 !     Mercury -      Silicon
      xfac(80,14) =     3.680740d0 !     Mercury -      Silicon
      alpb(80,15) =     0.891179d0 !     Mercury -   Phosphorus
      xfac(80,15) =     1.351633d0 !     Mercury -   Phosphorus
      alpb(80,16) =     1.900145d0 !     Mercury -       Sulfur
      xfac(80,16) =     1.772158d0 !     Mercury -       Sulfur
      alpb(80,17) =     1.838378d0 !     Mercury -     Chlorine
      xfac(80,17) =     1.142381d0 !     Mercury -     Chlorine
      alpb(80,22) =     3.414630d0 !     Mercury -     Titanium
      xfac(80,22) =     2.957200d0 !     Mercury -     Titanium
      alpb(80,34) =     1.607270d0 !     Mercury -     Selenium
      xfac(80,34) =     1.015811d0 !     Mercury -     Selenium
      alpb(80,35) =     1.705395d0 !     Mercury -      Bromine
      xfac(80,35) =     1.997568d0 !     Mercury -      Bromine
      alpb(80,52) =     1.536568d0 !     Mercury -    Tellurium
      xfac(80,52) =     4.486299d0 !     Mercury -    Tellurium
      alpb(80,53) =     1.476731d0 !     Mercury -       Iodine
      xfac(80,53) =     2.489683d0 !     Mercury -       Iodine
      alpb(80,80) =     2.288223d0 !     Mercury -      Mercury
      xfac(80,80) =    29.334203d0 !     Mercury -      Mercury
 !
      alpb(81, 1) =     2.098110d0 !    Thallium -     Hydrogen
      xfac(81, 1) =     2.104104d0 !    Thallium -     Hydrogen
      alpb(81, 5) =     1.558857d0 !    Thallium -        Boron
      xfac(81, 5) =     8.505888d0 !    Thallium -        Boron
      alpb(81, 6) =     2.721075d0 !    Thallium -       Carbon
      xfac(81, 6) =     5.320930d0 !    Thallium -       Carbon
      alpb(81, 7) =     1.975560d0 !    Thallium -     Nitrogen
      xfac(81, 7) =     1.168533d0 !    Thallium -     Nitrogen
      alpb(81, 8) =     3.327926d0 !    Thallium -       Oxygen
      xfac(81, 8) =    14.162059d0 !    Thallium -       Oxygen
      alpb(81, 9) =     3.188782d0 !    Thallium -     Fluorine
      xfac(81, 9) =     5.857103d0 !    Thallium -     Fluorine
      alpb(81,13) =     1.458567d0 !    Thallium -     Aluminum
      xfac(81,13) =     7.820754d0 !    Thallium -     Aluminum
      alpb(81,16) =     2.648560d0 !    Thallium -       Sulfur
      xfac(81,16) =     7.355971d0 !    Thallium -       Sulfur
      alpb(81,17) =     3.127130d0 !    Thallium -     Chlorine
      xfac(81,17) =    15.106797d0 !    Thallium -     Chlorine
      alpb(81,35) =     2.596395d0 !    Thallium -      Bromine
      xfac(81,35) =     9.792163d0 !    Thallium -      Bromine
      alpb(81,37) =     1.512698d0 !    Thallium -     Rubidium
      xfac(81,37) =    19.748653d0 !    Thallium -     Rubidium
      alpb(81,53) =     2.616072d0 !    Thallium -       Iodine
      xfac(81,53) =    30.000873d0 !    Thallium -       Iodine
      alpb(81,81) =     2.597707d0 !    Thallium -     Thallium
      xfac(81,81) =    32.531404d0 !    Thallium -     Thallium
 !
      alpb(82, 1) =     2.827636d0 !        Lead -     Hydrogen
      xfac(82, 1) =    11.387111d0 !        Lead -     Hydrogen
      alpb(82, 3) =     0.947660d0 !        Lead -      Lithium
      xfac(82, 3) =     1.191773d0 !        Lead -      Lithium
      alpb(82, 5) =     1.718658d0 !        Lead -        Boron
      xfac(82, 5) =     2.379671d0 !        Lead -        Boron
      alpb(82, 6) =     2.608618d0 !        Lead -       Carbon
      xfac(82, 6) =     6.379969d0 !        Lead -       Carbon
      alpb(82, 7) =     1.830414d0 !        Lead -     Nitrogen
      xfac(82, 7) =     1.025862d0 !        Lead -     Nitrogen
      alpb(82, 8) =     2.980453d0 !        Lead -       Oxygen
      xfac(82, 8) =     4.952015d0 !        Lead -       Oxygen
      alpb(82, 9) =     3.676637d0 !        Lead -     Fluorine
      xfac(82, 9) =     9.548494d0 !        Lead -     Fluorine
      alpb(82,14) =     1.284024d0 !        Lead -      Silicon
      xfac(82,14) =     1.429987d0 !        Lead -      Silicon
      alpb(82,15) =     1.955648d0 !        Lead -   Phosphorus
      xfac(82,15) =    13.812691d0 !        Lead -   Phosphorus
      alpb(82,16) =     2.362854d0 !        Lead -       Sulfur
      xfac(82,16) =     5.436090d0 !        Lead -       Sulfur
      alpb(82,17) =     1.499678d0 !        Lead -     Chlorine
      xfac(82,17) =     0.736101d0 !        Lead -     Chlorine
      alpb(82,20) =     1.691921d0 !        Lead -      Calcium
      xfac(82,20) =     6.709141d0 !        Lead -      Calcium
      alpb(82,23) =     1.501708d0 !        Lead -     Vanadium
      xfac(82,23) =     3.819119d0 !        Lead -     Vanadium
      alpb(82,24) =     1.305185d0 !        Lead -     Chromium
      xfac(82,24) =     0.974428d0 !        Lead -     Chromium
      alpb(82,30) =     1.434109d0 !        Lead -         Zinc
      xfac(82,30) =     2.865002d0 !        Lead -         Zinc
      alpb(82,33) =     1.792215d0 !        Lead -      Arsenic
      xfac(82,33) =     4.943875d0 !        Lead -      Arsenic
      alpb(82,34) =     2.893161d0 !        Lead -     Selenium
      xfac(82,34) =    29.986812d0 !        Lead -     Selenium
      alpb(82,35) =     2.364003d0 !        Lead -      Bromine
      xfac(82,35) =     6.777519d0 !        Lead -      Bromine
      alpb(82,41) =     1.500000d0 !        Lead -      Niobium
      xfac(82,41) =     1.000000d0 !        Lead -      Niobium
      alpb(82,42) =     1.759074d0 !        Lead -   Molybdenum
      xfac(82,42) =     5.265939d0 !        Lead -   Molybdenum
      alpb(82,52) =     3.242448d0 !        Lead -    Tellurium
      xfac(82,52) =   176.768383d0 !        Lead -    Tellurium
      alpb(82,53) =     2.179090d0 !        Lead -       Iodine
      xfac(82,53) =     8.112077d0 !        Lead -       Iodine
      alpb(82,74) =     1.517042d0 !        Lead -     Tungsten
      xfac(82,74) =     1.512242d0 !        Lead -     Tungsten
      alpb(82,82) =     2.529682d0 !        Lead -         Lead
      xfac(82,82) =    38.479040d0 !        Lead -         Lead
 !
      alpb(83, 1) =     1.727556d0 !     Bismuth -     Hydrogen
      xfac(83, 1) =     1.225129d0 !     Bismuth -     Hydrogen
      alpb(83, 3) =     1.221685d0 !     Bismuth -      Lithium
      xfac(83, 3) =     2.187383d0 !     Bismuth -      Lithium
      alpb(83, 6) =     1.970985d0 !     Bismuth -       Carbon
      xfac(83, 6) =     1.397988d0 !     Bismuth -       Carbon
      alpb(83, 7) =     1.976984d0 !     Bismuth -     Nitrogen
      xfac(83, 7) =     1.315182d0 !     Bismuth -     Nitrogen
      alpb(83, 8) =     2.337898d0 !     Bismuth -       Oxygen
      xfac(83, 8) =     1.621567d0 !     Bismuth -       Oxygen
      alpb(83, 9) =     2.029420d0 !     Bismuth -     Fluorine
      xfac(83, 9) =     0.490733d0 !     Bismuth -     Fluorine
      alpb(83,11) =     1.532800d0 !     Bismuth -       Sodium
      xfac(83,11) =     2.410886d0 !     Bismuth -       Sodium
      alpb(83,16) =     1.866193d0 !     Bismuth -       Sulfur
      xfac(83,16) =     1.624988d0 !     Bismuth -       Sulfur
      alpb(83,17) =     1.405944d0 !     Bismuth -     Chlorine
      xfac(83,17) =     0.496440d0 !     Bismuth -     Chlorine
      alpb(83,19) =     1.417970d0 !     Bismuth -    Potassium
      xfac(83,19) =     2.123183d0 !     Bismuth -    Potassium
      alpb(83,34) =     1.609528d0 !     Bismuth -     Selenium
      xfac(83,34) =     1.139985d0 !     Bismuth -     Selenium
      alpb(83,35) =     1.750597d0 !     Bismuth -      Bromine
      xfac(83,35) =     1.792203d0 !     Bismuth -      Bromine
      alpb(83,37) =     1.528441d0 !     Bismuth -     Rubidium
      xfac(83,37) =     2.435372d0 !     Bismuth -     Rubidium
      alpb(83,53) =     1.592333d0 !     Bismuth -       Iodine
      xfac(83,53) =     2.364966d0 !     Bismuth -       Iodine
      alpb(83,55) =     1.567880d0 !     Bismuth -       Cesium
      xfac(83,55) =     2.314878d0 !     Bismuth -       Cesium
      alpb(83,83) =     1.756620d0 !     Bismuth -      Bismuth
      xfac(83,83) =     7.710187d0 !     Bismuth -      Bismuth
 !
      alpb(87, 7) =     2.218810d0 !    Francium -     Nitrogen
      xfac(87, 7) =     1.012630d0 !    Francium -     Nitrogen
      alpb(87, 9) =     2.218810d0 !    Francium -     Fluorine
      xfac(87, 9) =     1.012630d0 !    Francium -     Fluorine
      alpb(87,17) =     1.579660d0 !    Francium -     Chlorine
      xfac(87,17) =     0.761560d0 !    Francium -     Chlorine
      alpb(87,87) =     1.579660d0 !    Francium -     Francium
      xfac(87,87) =     0.761560d0 !    Francium -     Francium
    end subroutine alpb_and_xfac_pm7
  end module Parameters_for_PM7_C
