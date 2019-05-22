  module Parameters_for_PM8_C
    use vast_kind_param, ONLY:  double  
    real(double), dimension(107) :: uss8, upp8, udd8, zs8, zp8, zd8, betas8, &
    betap8, betad8, gss8, gsp8, gpp8, gp28, hsp8, polvo8, poc_8, &
    zsn8, zpn8, zdn8, f0sd8, g2sd8, alp8, &
    CPE_Zet8, CPE_Z08, CPE_B8, CPE_Xlo8, CPE_Xhi8
    real(double) :: v_par8(60) 
    real(double), dimension(107,4) :: gues81, gues82, gues83
!
!                    Data for Element   1         Hydrogen
!
      data     uss8(  1)/       -11.396427D0/
      data   betas8(  1)/        -6.173787D0/
      data      zs8(  1)/         1.188078D0/
      data     alp8(  1)/         2.882324D0/
      data     gss8(  1)/        12.848000D0/
      data   polvo8(  1)/         0.178581D0/
      data CPE_Zet8(  1)/         1.342400D0/
      data  CPE_Z08(  1)/         2.255100D0/
      data   CPE_B8(  1)/         0.856600D0/
      data CPE_Xlo8(  1)/         0.379600D0/
      data CPE_Xhi8(  1)/         1.379600D0/
      data gues81(  1,1)/         0.122796D0/
      data gues82(  1,1)/         5.000000D0/
      data gues83(  1,1)/         1.200000D0/
      data gues81(  1,2)/         0.005090D0/
      data gues82(  1,2)/         5.000000D0/
      data gues83(  1,2)/         1.800000D0/
      data gues81(  1,3)/        -0.018336D0/
      data gues82(  1,3)/         2.000000D0/
      data gues83(  1,3)/         2.100000D0/
!
!                    Data for Element   2           Helium
!
      data     uss8(  2)/       -35.493117D0/
      data     upp8(  2)/         9.999843D0/
      data   betas8(  2)/       -17.554502D0/
      data   betap8(  2)/       -26.052606D0/
      data      zs8(  2)/         1.771076D0/
      data      zp8(  2)/         6.901826D0/
      data     alp8(  2)/         5.781393D0/
      data     gss8(  2)/         8.059014D0/
      data     gsp8(  2)/        11.015223D0/
      data     gpp8(  2)/        19.384334D0/
      data     gp28(  2)/        11.351722D0/
      data     hsp8(  2)/         0.532120D0/
      data gues81(  2,1)/         0.198598D0/
      data gues82(  2,1)/         1.000959D0/
      data gues83(  2,1)/         0.503151D0/
!
!                    Data for Element   3          Lithium
!
      data     uss8(  3)/        -5.300000D0/
      data     upp8(  3)/        -3.400000D0/
      data   betas8(  3)/        -0.550000D0/
      data   betap8(  3)/        -1.500000D0/
      data      zs8(  3)/         0.650000D0/
      data      zp8(  3)/         0.750000D0/
      data     alp8(  3)/         1.255000D0/
      data     gss8(  3)/         4.500000D0/
      data     gsp8(  3)/         3.000000D0/
      data     gpp8(  3)/         5.250000D0/
      data     gp28(  3)/         4.500000D0/
      data     hsp8(  3)/         0.150000D0/
      data gues81(  3,1)/        -0.450000D0/
      data gues82(  3,1)/         5.000000D0/
      data gues83(  3,1)/         1.000000D0/
      data gues81(  3,2)/         0.800000D0/
      data gues82(  3,2)/         6.500000D0/
      data gues83(  3,2)/         1.000000D0/
!
!                    Data for Element   4        Beryllium
!
      data     uss8(  4)/       -17.264752D0/
      data     upp8(  4)/       -11.304243D0/
      data   betas8(  4)/        -3.962053D0/
      data   betap8(  4)/        -2.780684D0/
      data      zs8(  4)/         0.877439D0/
      data      zp8(  4)/         1.508755D0/
      data     alp8(  4)/         1.593536D0/
      data     gss8(  4)/         9.012851D0/
      data     gsp8(  4)/         6.576199D0/
      data     gpp8(  4)/         6.057182D0/
      data     gp28(  4)/         9.005219D0/
      data     hsp8(  4)/         0.544679D0/
      data gues81(  4,1)/         1.631572D0/
      data gues82(  4,1)/         2.672962D0/
      data gues83(  4,1)/         1.791686D0/
      data gues81(  4,2)/        -2.110959D0/
      data gues82(  4,2)/         1.968594D0/
      data gues83(  4,2)/         1.755871D0/
!
!                    Data for Element   5            Boron
!
      data     uss8(  5)/       -50.477683D0/
      data     upp8(  5)/       -37.411983D0/
      data   betas8(  5)/       -10.549726D0/
      data   betap8(  5)/        -3.999595D0/
      data      zs8(  5)/         1.531260D0/
      data      zp8(  5)/         1.143460D0/
      data     alp8(  5)/         2.210416D0/
      data     gss8(  5)/        18.278280D0/
      data     gsp8(  5)/        15.333067D0/
      data     gpp8(  5)/        12.315858D0/
      data     gp28(  5)/        11.178535D0/
      data     hsp8(  5)/         0.599788D0/
      data gues81(  5,1)/        -0.351841D0/
      data gues82(  5,1)/         3.000862D0/
      data gues83(  5,1)/         0.824118D0/
!
!                    Data for Element   6           Carbon
!
      data     uss8(  6)/       -52.028658D0/
      data     upp8(  6)/       -39.614239D0/
      data   betas8(  6)/       -15.715783D0/
      data   betap8(  6)/        -7.719283D0/
      data      zs8(  6)/         1.808665D0/
      data      zp8(  6)/         1.685116D0/
      data     alp8(  6)/         2.648274D0/
      data     gss8(  6)/        12.230000D0/
      data     gsp8(  6)/        11.470000D0/
      data     gpp8(  6)/        11.080000D0/
      data     gp28(  6)/         9.840000D0/
      data     hsp8(  6)/         2.430000D0/
      data   polvo8(  6)/         1.242730D0/
      data CPE_Zet8(  6)/         1.167040D0/
      data  CPE_Z08(  6)/         1.178300D0/
      data   CPE_B8(  6)/         0.004800D0/
      data CPE_Xlo8(  6)/         1.086200D0/
      data CPE_Xhi8(  6)/         2.353000D0/
      data gues81(  6,1)/         0.011355D0/
      data gues82(  6,1)/         5.000000D0/
      data gues83(  6,1)/         1.600000D0/
      data gues81(  6,2)/         0.045924D0/
      data gues82(  6,2)/         5.000000D0/
      data gues83(  6,2)/         1.850000D0/
      data gues81(  6,3)/        -0.020061D0/
      data gues82(  6,3)/         5.000000D0/
      data gues83(  6,3)/         2.050000D0/
      data gues81(  6,4)/        -0.001260D0/
      data gues82(  6,4)/         5.000000D0/
      data gues83(  6,4)/         2.650000D0/
!
!                    Data for Element   7         Nitrogen
!
      data     uss8(  7)/       -49.335672D0/
      data     upp8(  7)/       -47.509736D0/
      data   betas8(  7)/       -14.062521D0/
      data   betap8(  7)/       -20.043848D0/
      data      zs8(  7)/         2.028094D0/
      data      zp8(  7)/         2.313728D0/
      data     alp8(  7)/         2.830545D0/
      data     gss8(  7)/        11.904787D0/
      data     gsp8(  7)/         7.348565D0/
      data     gpp8(  7)/        11.754672D0/
      data     gp28(  7)/        10.807277D0/
      data     hsp8(  7)/         1.136713D0/
      data   polvo8(  7)/         0.920190D0/
      data CPE_Zet8(  7)/         1.378880D0/
      data  CPE_Z08(  7)/         2.029200D0/
      data   CPE_B8(  7)/         0.323800D0/
      data CPE_Xlo8(  7)/         1.651100D0/
      data CPE_Xhi8(  7)/         2.292100D0/
      data gues81(  7,1)/         1.501674D0/
      data gues82(  7,1)/         5.901148D0/
      data gues83(  7,1)/         1.710740D0/
      data gues81(  7,2)/        -1.505772D0/
      data gues82(  7,2)/         6.004658D0/
      data gues83(  7,2)/         1.716149D0/
!
!                    Data for Element   8           Oxygen
!
      data     uss8(  8)/       -97.830000D0/
      data     upp8(  8)/       -78.262380D0/
      data   betas8(  8)/       -29.272773D0/
      data   betap8(  8)/       -29.272773D0/
      data      zs8(  8)/         3.108032D0/
      data      zp8(  8)/         2.524039D0/
      data     alp8(  8)/         4.455371D0/
      data     gss8(  8)/        15.420000D0/
      data     gsp8(  8)/        14.480000D0/
      data     gpp8(  8)/        14.520000D0/
      data     gp28(  8)/        12.980000D0/
      data     hsp8(  8)/         3.940000D0/
      data   polvo8(  8)/         0.434605D0/
      data CPE_Zet8(  8)/         1.585280D0/
      data  CPE_Z08(  8)/         4.322700D0/
      data   CPE_B8(  8)/         0.045100D0/
      data CPE_Xlo8(  8)/         3.483200D0/
      data CPE_Xhi8(  8)/         3.605000D0/
      data gues81(  8,1)/         0.280962D0/
      data gues82(  8,1)/         5.000000D0/
      data gues83(  8,1)/         0.847918D0/
      data gues81(  8,2)/         0.081430D0/
      data gues82(  8,2)/         7.000000D0/
      data gues83(  8,2)/         1.445071D0/
!
!                    Data for Element   9         Fluorine
!
      data     uss8(  9)/      -110.435303D0/
      data     upp8(  9)/      -105.685047D0/
      data   betas8(  9)/       -48.405939D0/
      data   betap8(  9)/       -27.744660D0/
      data      zs8(  9)/         4.708555D0/
      data      zp8(  9)/         2.491178D0/
      data     alp8(  9)/         3.358921D0/
      data     gss8(  9)/        10.496667D0/
      data     gsp8(  9)/        16.073689D0/
      data     gpp8(  9)/        14.817256D0/
      data     gp28(  9)/        14.418393D0/
      data     hsp8(  9)/         0.727763D0/
      data   polvo8(  9)/         0.235893D0/
      data gues81(  9,1)/        -0.012166D0/
      data gues82(  9,1)/         6.023574D0/
      data gues83(  9,1)/         1.856859D0/
      data gues81(  9,2)/        -0.002852D0/
      data gues82(  9,2)/         6.003717D0/
      data gues83(  9,2)/         2.636158D0/
!
!                    Data for Element  10             Neon
!
      data     uss8( 10)/         9.655474D0/
      data     upp8( 10)/       -71.146647D0/
      data   betas8( 10)/        -0.151266D0/
      data   betap8( 10)/       -23.847019D0/
      data      zs8( 10)/         5.999874D0/
      data      zp8( 10)/         4.175260D0/
      data     alp8( 10)/         2.537726D0/
      data     gss8( 10)/         0.498937D0/
      data     gsp8( 10)/        10.159204D0/
      data     gpp8( 10)/        18.944998D0/
      data     gp28( 10)/         8.564859D0/
      data     hsp8( 10)/         0.299999D0/
      data gues81( 10,1)/         0.235211D0/
      data gues82( 10,1)/         8.811661D0/
      data gues83( 10,1)/         1.092788D0/
!
!                    Data for Element  11           Sodium
!
      data     uss8( 11)/        -5.014059D0/
      data     upp8( 11)/        -2.925230D0/
      data   betas8( 11)/        -3.948634D0/
      data   betap8( 11)/        -4.240181D0/
      data      zs8( 11)/         2.661894D0/
      data      zp8( 11)/         0.883742D0/
      data     alp8( 11)/         0.910918D0/
      data     gss8( 11)/         5.797525D0/
      data     gsp8( 11)/        10.000006D0/
      data     gpp8( 11)/         1.224790D0/
      data     gp28( 11)/         0.999650D0/
      data     hsp8( 11)/         0.399936D0/
!
!                    Data for Element  12        Magnesium
!
      data     uss8( 12)/       -14.623688D0/
      data     upp8( 12)/       -14.173460D0/
      data   betas8( 12)/        -2.071691D0/
      data   betap8( 12)/        -0.569581D0/
      data      zs8( 12)/         0.698552D0/
      data      zp8( 12)/         1.483453D0/
      data     alp8( 12)/         1.329147D0/
      data     gss8( 12)/         6.694300D0/
      data     gsp8( 12)/         6.793995D0/
      data     gpp8( 12)/         6.910446D0/
      data     gp28( 12)/         7.090823D0/
      data     hsp8( 12)/         0.543300D0/
      data gues81( 12,1)/         2.117050D0/
      data gues82( 12,1)/         6.009477D0/
      data gues83( 12,1)/         2.084406D0/
      data gues81( 12,2)/        -2.547767D0/
      data gues82( 12,2)/         4.395370D0/
      data gues83( 12,2)/         2.063674D0/
!
!                    Data for Element  13         Aluminum
!
      data     uss8( 13)/       -24.845404D0/
      data     upp8( 13)/       -22.264159D0/
      data   betas8( 13)/        -0.594301D0/
      data   betap8( 13)/        -0.956550D0/
      data      zs8( 13)/         1.702888D0/
      data      zp8( 13)/         1.073629D0/
      data     alp8( 13)/         1.521703D0/
      data     gss8( 13)/         5.776737D0/
      data     gsp8( 13)/        11.659856D0/
      data     gpp8( 13)/         6.347790D0/
      data     gp28( 13)/         6.121077D0/
      data     hsp8( 13)/         4.006245D0/
      data gues81( 13,1)/        -0.473090D0/
      data gues82( 13,1)/         1.915825D0/
      data gues83( 13,1)/         1.451728D0/
      data gues81( 13,2)/        -0.154051D0/
      data gues82( 13,2)/         6.005086D0/
      data gues83( 13,2)/         2.519997D0/
!
!                    Data for Element  14          Silicon
!
      data     uss8( 14)/       -26.763483D0/
      data     upp8( 14)/       -22.813635D0/
      data   betas8( 14)/        -2.862145D0/
      data   betap8( 14)/        -3.933148D0/
      data      zs8( 14)/         1.635075D0/
      data      zp8( 14)/         1.313088D0/
      data     alp8( 14)/         2.135809D0/
      data     gss8( 14)/         5.047196D0/
      data     gsp8( 14)/         5.949057D0/
      data     gpp8( 14)/         6.759367D0/
      data     gp28( 14)/         5.161297D0/
      data     hsp8( 14)/         0.919832D0/
      data gues81( 14,1)/        -0.390600D0/
      data gues82( 14,1)/         6.000054D0/
      data gues83( 14,1)/         0.632262D0/
      data gues81( 14,2)/         0.057259D0/
      data gues82( 14,2)/         6.007183D0/
      data gues83( 14,2)/         2.019987D0/
!
!                    Data for Element  15       Phosphorus
!
      data     uss8( 15)/       -40.413096D0/
      data     upp8( 15)/       -29.593052D0/
      data   betas8( 15)/       -12.615879D0/
      data   betap8( 15)/        -4.160040D0/
      data      zs8( 15)/         2.017563D0/
      data      zp8( 15)/         1.504732D0/
      data     alp8( 15)/         1.940534D0/
      data     gss8( 15)/         7.801615D0/
      data     gsp8( 15)/         5.186949D0/
      data     gpp8( 15)/         6.618478D0/
      data     gp28( 15)/         6.062002D0/
      data     hsp8( 15)/         1.542809D0/
      data gues81( 15,1)/        -0.611421D0/
      data gues82( 15,1)/         1.997272D0/
      data gues83( 15,1)/         0.794624D0/
      data gues81( 15,2)/        -0.093935D0/
      data gues82( 15,2)/         1.998360D0/
      data gues83( 15,2)/         1.910677D0/
!
!                    Data for Element  16           Sulfur
!
      data     uss8( 16)/       -49.895371D0/
      data     upp8( 16)/       -44.392583D0/
      data   betas8( 16)/        -8.827465D0/
      data   betap8( 16)/        -8.091415D0/
      data      zs8( 16)/         1.891185D0/
      data      zp8( 16)/         1.658972D0/
      data     alp8( 16)/         2.269706D0/
      data     gss8( 16)/         8.964667D0/
      data     gsp8( 16)/         6.785936D0/
      data     gpp8( 16)/         9.968164D0/
      data     gp28( 16)/         7.970247D0/
      data     hsp8( 16)/         4.041836D0/
      data gues81( 16,1)/        -0.399191D0/
      data gues82( 16,1)/         6.000669D0/
      data gues83( 16,1)/         0.962123D0/
      data gues81( 16,2)/        -0.054899D0/
      data gues82( 16,2)/         6.001845D0/
      data gues83( 16,2)/         1.579944D0/
!
!                    Data for Element  17         Chlorine
!
      data     uss8( 17)/      -100.626747D0/
      data     upp8( 17)/       -53.614396D0/
      data   betas8( 17)/       -27.528560D0/
      data   betap8( 17)/       -11.593922D0/
      data      zs8( 17)/         2.246210D0/
      data      zp8( 17)/         2.151010D0/
      data     alp8( 17)/         2.517296D0/
      data     gss8( 17)/        16.013601D0/
      data     gsp8( 17)/         8.048115D0/
      data     gpp8( 17)/         7.522215D0/
      data     gp28( 17)/         7.504154D0/
      data     hsp8( 17)/         3.481153D0/
      data   polvo8( 17)/         1.855600D0/
      data gues81( 17,1)/        -0.171591D0/
      data gues82( 17,1)/         6.000802D0/
      data gues83( 17,1)/         1.087502D0/
      data gues81( 17,2)/        -0.013458D0/
      data gues82( 17,2)/         1.966618D0/
      data gues83( 17,2)/         2.292891D0/
!
!                    Data for Element  18            Argon
!
      data     uss8( 18)/         4.275991D0/
      data     upp8( 18)/       -75.100410D0/
      data   betas8( 18)/        -2.432005D0/
      data   betap8( 18)/       -22.985694D0/
      data      zs8( 18)/         0.982170D0/
      data      zp8( 18)/         5.999715D0/
      data     alp8( 18)/         2.211332D0/
      data     gss8( 18)/         3.237741D0/
      data     gsp8( 18)/        16.357397D0/
      data     gpp8( 18)/        19.998465D0/
      data     gp28( 18)/        10.934324D0/
      data     hsp8( 18)/         0.816915D0/
      data gues81( 18,1)/        -0.737650D0/
      data gues82( 18,1)/         3.888254D0/
      data gues83( 18,1)/         0.715067D0/
!
!                    Data for Element  19        Potassium
!
      data     uss8( 19)/        -4.259647D0/
      data     upp8( 19)/        -2.642509D0/
      data   betas8( 19)/        -0.425039D0/
      data   betap8( 19)/        -3.199814D0/
      data      zs8( 19)/         0.810169D0/
      data      zp8( 19)/         0.957834D0/
      data     alp8( 19)/         0.725216D0/
      data     gss8( 19)/         6.778831D0/
      data     gsp8( 19)/         9.347280D0/
      data     gpp8( 19)/         3.496351D0/
      data     gp28( 19)/         2.741680D0/
      data     hsp8( 19)/         1.659246D0/
!
!                    Data for Element  20          Calcium
!
      data     uss8( 20)/       -11.960816D0/
      data     upp8( 20)/        -8.635086D0/
      data   betas8( 20)/        -0.986827D0/
      data   betap8( 20)/        -3.968427D0/
      data      zs8( 20)/         1.208741D0/
      data      zp8( 20)/         0.940937D0/
      data     alp8( 20)/         0.499998D0/
      data     gss8( 20)/         5.971991D0/
      data     gsp8( 20)/         4.960778D0/
      data     gpp8( 20)/         3.721490D0/
      data     gp28( 20)/         3.711650D0/
      data     hsp8( 20)/         0.792815D0/
!
!                    Data for Element  21         Scandium
!
      data     udd8( 21)/       -22.018012D0/
      data   betad8( 21)/        -0.999956D0/
      data     zsn8( 21)/         0.501877D0/
      data     zpn8( 21)/         1.105305D0/
      data     zdn8( 21)/         1.687497D0/
      data     gss8( 21)/         2.743711D0/
      data     gsp8( 21)/         3.351775D0/
      data     gpp8( 21)/         6.584268D0/
      data     gp28( 21)/         5.771753D0/
      data     hsp8( 21)/         0.271247D0/
      data    f0sd8( 21)/         6.960452D0/
      data    g2sd8( 21)/         1.164321D0/
!
!                    Data for Element  22         Titanium
!
      data     udd8( 22)/       -43.136457D0/
      data   betad8( 22)/       -12.069692D0/
      data     zsn8( 22)/         1.778401D0/
      data     zpn8( 22)/         1.564743D0/
      data     zdn8( 22)/         1.957840D0/
      data     gss8( 22)/         9.722339D0/
      data     gsp8( 22)/         9.065692D0/
      data     gpp8( 22)/         9.321123D0/
      data     gp28( 22)/         8.170874D0/
      data     hsp8( 22)/         2.134882D0/
      data    f0sd8( 22)/        11.039032D0/
      data    g2sd8( 22)/         4.883814D0/
!
!                    Data for Element  23         Vanadium
!
      data     udd8( 23)/       -49.364215D0/
      data   betad8( 23)/        -4.862075D0/
      data     zsn8( 23)/         1.151717D0/
      data     zpn8( 23)/         0.649694D0/
      data     zdn8( 23)/         1.706071D0/
      data     gss8( 23)/         6.296318D0/
      data     gsp8( 23)/         4.238285D0/
      data     gpp8( 23)/         3.870206D0/
      data     gp28( 23)/         3.392613D0/
      data     hsp8( 23)/         0.576541D0/
      data    f0sd8( 23)/         7.631125D0/
      data    g2sd8( 23)/         1.387749D0/
!
!                    Data for Element  24         Chromium
!
      data     udd8( 24)/       -55.626233D0/
      data   betad8( 24)/       -15.691040D0/
      data     zsn8( 24)/         0.918980D0/
      data     zpn8( 24)/         0.665550D0/
      data     zdn8( 24)/         1.504446D0/
      data     gss8( 24)/         5.023973D0/
      data     gsp8( 24)/         4.120900D0/
      data     gpp8( 24)/         3.964660D0/
      data     gp28( 24)/         3.475411D0/
      data     hsp8( 24)/         0.831435D0/
      data    f0sd8( 24)/         7.933107D0/
      data    g2sd8( 24)/         1.011301D0/
!
!                    Data for Element  26             Iron
!
      data     udd8( 26)/       -95.674336D0/
      data   betad8( 26)/       -11.951200D0/
      data     zsn8( 26)/         1.245692D0/
      data     zpn8( 26)/         1.821407D0/
      data     zdn8( 26)/         1.935644D0/
      data     gss8( 26)/         6.810073D0/
      data     gsp8( 26)/         7.829047D0/
      data     gpp8( 26)/        10.850066D0/
      data     gp28( 26)/         9.511142D0/
      data     hsp8( 26)/         1.471440D0/
      data    f0sd8( 26)/         8.917910D0/
      data    g2sd8( 26)/         1.187159D0/
!
!                    Data for Element  27           Cobalt
!
      data     udd8( 27)/       -88.956780D0/
      data   betad8( 27)/       -17.205910D0/
      data     zsn8( 27)/         0.817217D0/
      data     zpn8( 27)/         2.025146D0/
      data     zdn8( 27)/         1.600832D0/
      data     gss8( 27)/         4.467643D0/
      data     gsp8( 27)/         5.497106D0/
      data     gpp8( 27)/        12.063731D0/
      data     gp28( 27)/        10.575037D0/
      data     hsp8( 27)/         0.313579D0/
      data    f0sd8( 27)/         6.789595D0/
      data    g2sd8( 27)/         1.349020D0/
!
!                    Data for Element  28           Nickel
!
      data     udd8( 28)/       -68.078673D0/
      data   betad8( 28)/       -26.880782D0/
      data     zsn8( 28)/         2.015549D0/
      data     zpn8( 28)/         4.000002D0/
      data     zdn8( 28)/         0.929788D0/
      data     gss8( 28)/        11.018802D0/
      data     gsp8( 28)/        13.336339D0/
      data     gpp8( 28)/        23.827887D0/
      data     gp28( 28)/        20.887467D0/
      data     hsp8( 28)/         1.413045D0/
      data    f0sd8( 28)/         8.028073D0/
      data    g2sd8( 28)/         2.468167D0/
!
!                    Data for Element  29           Copper
!
      data     udd8( 29)/      -173.906811D0/
      data   betad8( 29)/       -25.652027D0/
      data     zsn8( 29)/         3.945734D0/
      data     zpn8( 29)/         1.105085D0/
      data     zdn8( 29)/         2.426250D0/
      data     gss8( 29)/        21.570928D0/
      data     gsp8( 29)/         7.502713D0/
      data     gpp8( 29)/         6.582958D0/
      data     gp28( 29)/         5.770605D0/
      data     hsp8( 29)/         0.109224D0/
      data    f0sd8( 29)/        18.892745D0/
      data    g2sd8( 29)/         9.509410D0/
!
!                    Data for Element  30             Zinc
!
      data     uss8( 30)/       -18.532198D0/
      data     upp8( 30)/       -11.047409D0/
      data   betas8( 30)/        -0.715578D0/
      data   betap8( 30)/        -6.351864D0/
      data      zs8( 30)/         1.819989D0/
      data      zp8( 30)/         1.506922D0/
      data     alp8( 30)/         1.350126D0/
      data     gss8( 30)/         9.677196D0/
      data     gsp8( 30)/         7.736204D0/
      data     gpp8( 30)/         4.980174D0/
      data     gp28( 30)/         4.669656D0/
      data     hsp8( 30)/         0.600413D0/
      data gues81( 30,1)/        -0.111234D0/
      data gues82( 30,1)/         6.001478D0/
      data gues83( 30,1)/         1.516032D0/
      data gues81( 30,2)/        -0.132370D0/
      data gues82( 30,2)/         1.995839D0/
      data gues83( 30,2)/         2.519642D0/
!
!                    Data for Element  31          Gallium
!
      data     uss8( 31)/       -29.855593D0/
      data     upp8( 31)/       -21.875371D0/
      data   betas8( 31)/        -4.945618D0/
      data   betap8( 31)/        -0.407053D0/
      data      zs8( 31)/         1.847040D0/
      data      zp8( 31)/         0.839411D0/
      data     alp8( 31)/         1.605115D0/
      data     gss8( 31)/         8.458554D0/
      data     gsp8( 31)/         8.925619D0/
      data     gpp8( 31)/         5.086855D0/
      data     gp28( 31)/         4.983045D0/
      data     hsp8( 31)/         2.051260D0/
      data gues81( 31,1)/        -0.560179D0/
      data gues82( 31,1)/         5.623273D0/
      data gues83( 31,1)/         1.531780D0/
      data gues81( 31,2)/        -0.272731D0/
      data gues82( 31,2)/         1.991843D0/
      data gues83( 31,2)/         2.183864D0/
!
!                    Data for Element  32        Germanium
!
      data     uss8( 32)/       -35.467196D0/
      data     upp8( 32)/       -31.586358D0/
      data   betas8( 32)/        -5.325002D0/
      data   betap8( 32)/        -2.250157D0/
      data      zs8( 32)/         2.237353D0/
      data      zp8( 32)/         1.592432D0/
      data     alp8( 32)/         1.972337D0/
      data     gss8( 32)/         5.376963D0/
      data     gsp8( 32)/        10.209529D0/
      data     gpp8( 32)/         7.671865D0/
      data     gp28( 32)/         6.924266D0/
      data     hsp8( 32)/         1.337020D0/
      data gues81( 32,1)/         0.963173D0/
      data gues82( 32,1)/         6.012013D0/
      data gues83( 32,1)/         2.163365D0/
      data gues81( 32,2)/        -0.959389D0/
      data gues82( 32,2)/         5.749180D0/
      data gues83( 32,2)/         2.169372D0/
!
!                    Data for Element  33          Arsenic
!
      data     uss8( 33)/       -38.507424D0/
      data     upp8( 33)/       -35.152415D0/
      data   betas8( 33)/        -8.232165D0/
      data   betap8( 33)/        -5.017386D0/
      data      zs8( 33)/         2.636177D0/
      data      zp8( 33)/         1.703889D0/
      data     alp8( 33)/         1.794477D0/
      data     gss8( 33)/         8.789001D0/
      data     gsp8( 33)/         5.397983D0/
      data     gpp8( 33)/         8.287250D0/
      data     gp28( 33)/         8.210346D0/
      data     hsp8( 33)/         1.951034D0/
      data gues81( 33,1)/        -0.460095D0/
      data gues82( 33,1)/         1.983115D0/
      data gues83( 33,1)/         1.086793D0/
      data gues81( 33,2)/        -0.088996D0/
      data gues82( 33,2)/         1.992944D0/
      data gues83( 33,2)/         2.140058D0/
!
!                    Data for Element  34         Selenium
!
      data     uss8( 34)/       -55.378135D0/
      data     upp8( 34)/       -49.823076D0/
      data   betas8( 34)/        -6.157822D0/
      data   betap8( 34)/        -5.493039D0/
      data      zs8( 34)/         2.828051D0/
      data      zp8( 34)/         1.732536D0/
      data     alp8( 34)/         3.043957D0/
      data     gss8( 34)/         7.432591D0/
      data     gsp8( 34)/        10.060461D0/
      data     gpp8( 34)/         9.568326D0/
      data     gp28( 34)/         7.724289D0/
      data     hsp8( 34)/         4.016558D0/
      data gues81( 34,1)/         0.047873D0/
      data gues82( 34,1)/         6.007400D0/
      data gues83( 34,1)/         2.081717D0/
      data gues81( 34,2)/         0.114720D0/
      data gues82( 34,2)/         6.008672D0/
      data gues83( 34,2)/         1.516423D0/
!
!                    Data for Element  35          Bromine
!
      data     uss8( 35)/      -116.619311D0/
      data     upp8( 35)/       -74.227129D0/
      data   betas8( 35)/       -31.171342D0/
      data   betap8( 35)/        -6.814013D0/
      data      zs8( 35)/         5.348457D0/
      data      zp8( 35)/         2.127590D0/
      data     alp8( 35)/         2.511842D0/
      data     gss8( 35)/        15.943425D0/
      data     gsp8( 35)/        16.061680D0/
      data     gpp8( 35)/         8.282763D0/
      data     gp28( 35)/         7.816849D0/
      data     hsp8( 35)/         0.578869D0/
      data   polvo8( 35)/         2.845420D0/
      data gues81( 35,1)/         0.960458D0/
      data gues82( 35,1)/         5.976508D0/
      data gues83( 35,1)/         2.321654D0/
      data gues81( 35,2)/        -0.954916D0/
      data gues82( 35,2)/         5.944703D0/
      data gues83( 35,2)/         2.328142D0/
!
!                    Data for Element  36          Krypton
!
      data     uss8( 36)/         9.798694D0/
      data     upp8( 36)/       -73.859566D0/
      data   betas8( 36)/        -1.209471D0/
      data   betap8( 36)/        -8.432827D0/
      data      zs8( 36)/         3.560862D0/
      data      zp8( 36)/         1.983206D0/
      data     alp8( 36)/         1.702543D0/
      data     gss8( 36)/         0.499417D0/
      data     gsp8( 36)/         9.526126D0/
      data     gpp8( 36)/        19.998729D0/
      data     gp28( 36)/        11.198605D0/
      data     hsp8( 36)/         2.181834D0/
      data gues81( 36,1)/         0.780439D0/
      data gues82( 36,1)/         8.837140D0/
      data gues83( 36,1)/         0.500042D0/
!
!                    Data for Element  37         Rubidium
!
      data     uss8( 37)/        -4.592054D0/
      data     upp8( 37)/        -3.011921D0/
      data   betas8( 37)/       -12.089466D0/
      data   betap8( 37)/        -1.999927D0/
      data      zs8( 37)/         4.000042D0/
      data      zp8( 37)/         1.013459D0/
      data     alp8( 37)/         0.998532D0/
      data     gss8( 37)/         9.275721D0/
      data     gsp8( 37)/        20.000004D0/
      data     gpp8( 37)/        13.369414D0/
      data     gp28( 37)/        19.000003D0/
      data     hsp8( 37)/         4.998715D0/
!
!                    Data for Element  38        Strontium
!
      data     uss8( 38)/       -10.918344D0/
      data     upp8( 38)/        -7.983522D0/
      data   betas8( 38)/       -10.000007D0/
      data   betap8( 38)/        -5.675539D0/
      data      zs8( 38)/         1.279453D0/
      data      zp8( 38)/         1.391250D0/
      data     alp8( 38)/         1.646755D0/
      data     gss8( 38)/         5.036147D0/
      data     gsp8( 38)/         3.928401D0/
      data     gpp8( 38)/         3.153314D0/
      data     gp28( 38)/         3.245731D0/
      data     hsp8( 38)/         0.759621D0/
!
!                    Data for Element  40        Zirconium
!
      data     udd8( 40)/       -42.629405D0/
      data   betad8( 40)/       -11.194855D0/
      data     zsn8( 40)/         1.893199D0/
      data     zpn8( 40)/         1.418504D0/
      data     zdn8( 40)/         2.687782D0/
      data     gss8( 40)/         8.487900D0/
      data     gsp8( 40)/         7.107285D0/
      data     gpp8( 40)/         6.955633D0/
      data     gp28( 40)/         6.061690D0/
      data     hsp8( 40)/         1.470983D0/
      data    f0sd8( 40)/         9.961974D0/
      data    g2sd8( 40)/         3.022841D0/
!
!                    Data for Element  42       Molybdenum
!
      data     udd8( 42)/       -55.480489D0/
      data   betad8( 42)/       -15.234860D0/
      data     zsn8( 42)/         1.799704D0/
      data     zpn8( 42)/         1.780167D0/
      data     zdn8( 42)/         1.923161D0/
      data     gss8( 42)/         8.068729D0/
      data     gsp8( 42)/         8.024421D0/
      data     gpp8( 42)/         8.729045D0/
      data     gp28( 42)/         7.607181D0/
      data     hsp8( 42)/         1.997360D0/
      data    f0sd8( 42)/         7.669693D0/
      data    g2sd8( 42)/         0.999982D0/
!
!                    Data for Element  46        Palladium
!
      data     udd8( 46)/      -114.704547D0/
      data   betad8( 46)/       -14.176361D0/
      data     zsn8( 46)/         2.341668D0/
      data     zpn8( 46)/         3.991025D0/
      data     zdn8( 46)/         2.241264D0/
      data     gss8( 46)/        10.498550D0/
      data     gsp8( 46)/        12.324887D0/
      data     gpp8( 46)/        19.569984D0/
      data     gp28( 46)/        17.054833D0/
      data     hsp8( 46)/         1.632081D0/
      data    f0sd8( 46)/         8.962630D0/
      data    g2sd8( 46)/         2.818435D0/
!
!                    Data for Element  47           Silver
!
      data     udd8( 47)/      -137.085891D0/
      data   betad8( 47)/       -15.089393D0/
      data     zsn8( 47)/         1.488088D0/
      data     zpn8( 47)/         1.105479D0/
      data     zdn8( 47)/         2.442462D0/
      data     gss8( 47)/         6.671642D0/
      data     gsp8( 47)/         5.551865D0/
      data     gpp8( 47)/         5.420713D0/
      data     gp28( 47)/         4.724038D0/
      data     hsp8( 47)/         1.136355D0/
      data    f0sd8( 47)/         7.950899D0/
      data    g2sd8( 47)/         1.658340D0/
!
!                    Data for Element  48          Cadmium
!
      data     uss8( 48)/       -15.828584D0/
      data     upp8( 48)/         8.749795D0/
      data   betas8( 48)/        -8.581944D0/
      data   betap8( 48)/        -0.601034D0/
      data      zs8( 48)/         1.679351D0/
      data      zp8( 48)/         2.066412D0/
      data     alp8( 48)/         1.525382D0/
      data     gss8( 48)/         9.206960D0/
      data     gsp8( 48)/         8.231539D0/
      data     gpp8( 48)/         4.948104D0/
      data     gp28( 48)/         4.669656D0/
      data     hsp8( 48)/         1.656234D0/
!
!                    Data for Element  49           Indium
!
      data     uss8( 49)/       -26.176205D0/
      data     upp8( 49)/       -20.005822D0/
      data   betas8( 49)/        -2.993319D0/
      data   betap8( 49)/        -1.828908D0/
      data      zs8( 49)/         2.016116D0/
      data      zp8( 49)/         1.445350D0/
      data     alp8( 49)/         1.418385D0/
      data     gss8( 49)/         6.554900D0/
      data     gsp8( 49)/         8.229873D0/
      data     gpp8( 49)/         6.299269D0/
      data     gp28( 49)/         4.984211D0/
      data     hsp8( 49)/         2.631461D0/
      data gues81( 49,1)/        -0.343138D0/
      data gues82( 49,1)/         1.994034D0/
      data gues83( 49,1)/         1.625516D0/
      data gues81( 49,2)/        -0.109532D0/
      data gues82( 49,2)/         5.683217D0/
      data gues83( 49,2)/         2.867009D0/
!
!                    Data for Element  50              Tin
!
      data     uss8( 50)/       -34.550192D0/
      data     upp8( 50)/       -25.894419D0/
      data   betas8( 50)/        -2.785802D0/
      data   betap8( 50)/        -2.005999D0/
      data      zs8( 50)/         2.373328D0/
      data      zp8( 50)/         1.638233D0/
      data     alp8( 50)/         1.699650D0/
      data     gss8( 50)/        10.190033D0/
      data     gsp8( 50)/         7.235327D0/
      data     gpp8( 50)/         5.673810D0/
      data     gp28( 50)/         5.182214D0/
      data     hsp8( 50)/         1.033157D0/
      data gues81( 50,1)/        -0.150353D0/
      data gues82( 50,1)/         6.005694D0/
      data gues83( 50,1)/         1.704642D0/
      data gues81( 50,2)/        -0.044417D0/
      data gues82( 50,2)/         2.257381D0/
      data gues83( 50,2)/         2.469869D0/
!
!                    Data for Element  51         Antimony
!
      data     uss8( 51)/       -56.432196D0/
      data     upp8( 51)/       -29.434954D0/
      data   betas8( 51)/       -14.794217D0/
      data   betap8( 51)/        -2.817948D0/
      data      zs8( 51)/         2.343039D0/
      data      zp8( 51)/         1.899992D0/
      data     alp8( 51)/         2.034301D0/
      data     gss8( 51)/         9.238277D0/
      data     gsp8( 51)/         5.277680D0/
      data     gpp8( 51)/         6.350000D0/
      data     gp28( 51)/         6.250000D0/
      data     hsp8( 51)/         2.424464D0/
      data gues81( 51,1)/         3.002028D0/
      data gues82( 51,1)/         6.005342D0/
      data gues83( 51,1)/         0.853060D0/
      data gues81( 51,2)/        -0.018892D0/
      data gues82( 51,2)/         6.011478D0/
      data gues83( 51,2)/         2.793311D0/
!
!                    Data for Element  52        Tellurium
!
      data     uss8( 52)/       -44.938036D0/
      data     upp8( 52)/       -46.314099D0/
      data   betas8( 52)/        -2.665146D0/
      data   betap8( 52)/        -3.895430D0/
      data      zs8( 52)/         4.165492D0/
      data      zp8( 52)/         1.647555D0/
      data     alp8( 52)/         2.485019D0/
      data     gss8( 52)/        10.255073D0/
      data     gsp8( 52)/         8.169145D0/
      data     gpp8( 52)/         7.777592D0/
      data     gp28( 52)/         7.755121D0/
      data     hsp8( 52)/         3.772462D0/
      data gues81( 52,1)/         0.033391D0/
      data gues82( 52,1)/         5.956379D0/
      data gues83( 52,1)/         2.277575D0/
      data gues81( 52,2)/        -1.921867D0/
      data gues82( 52,2)/         4.973219D0/
      data gues83( 52,2)/         0.524243D0/
!
!                    Data for Element  53           Iodine
!
      data     uss8( 53)/       -96.454037D0/
      data     upp8( 53)/       -61.091582D0/
      data   betas8( 53)/       -14.494234D0/
      data   betap8( 53)/        -5.894703D0/
      data      zs8( 53)/         7.001013D0/
      data      zp8( 53)/         2.454354D0/
      data     alp8( 53)/         1.990185D0/
      data     gss8( 53)/        13.631943D0/
      data     gsp8( 53)/        14.990406D0/
      data     gpp8( 53)/         7.288330D0/
      data     gp28( 53)/         5.966407D0/
      data     hsp8( 53)/         2.630035D0/
      data   polvo8( 53)/         4.799730D0/
      data gues81( 53,1)/        -0.131481D0/
      data gues82( 53,1)/         5.206417D0/
      data gues83( 53,1)/         1.748824D0/
      data gues81( 53,2)/        -0.036897D0/
      data gues82( 53,2)/         6.010117D0/
      data gues83( 53,2)/         2.710373D0/
!
!                    Data for Element  54            Xenon
!
      data     uss8( 54)/         5.920589D0/
      data     upp8( 54)/       -86.982832D0/
      data   betas8( 54)/        -3.604845D0/
      data   betap8( 54)/        -4.967374D0/
      data      zs8( 54)/         4.990079D0/
      data      zp8( 54)/         2.692925D0/
      data     alp8( 54)/         1.794855D0/
      data     gss8( 54)/         2.187476D0/
      data     gsp8( 54)/         4.903868D0/
      data     gpp8( 54)/        12.493918D0/
      data     gp28( 54)/        14.952790D0/
      data     hsp8( 54)/         2.476853D0/
      data gues81( 54,1)/        -0.245535D0/
      data gues82( 54,1)/         2.058008D0/
      data gues83( 54,1)/         1.717330D0/
!
!                    Data for Element  55           Cesium
!
      data     uss8( 55)/        -3.203656D0/
      data     upp8( 55)/        -1.745197D0/
      data   betas8( 55)/        -0.602876D0/
      data   betap8( 55)/        -5.938609D0/
      data      zs8( 55)/         3.596030D0/
      data      zp8( 55)/         0.925517D0/
      data     alp8( 55)/         0.523834D0/
      data     gss8( 55)/         2.160566D0/
      data     gsp8( 55)/         4.166750D0/
      data     gpp8( 55)/         5.414048D0/
      data     gp28( 55)/         6.290494D0/
      data     hsp8( 55)/         0.399556D0/
!
!                    Data for Element  56           Barium
!
      data     uss8( 56)/       -10.102810D0/
      data     upp8( 56)/        -6.674384D0/
      data   betas8( 56)/       -10.000005D0/
      data   betap8( 56)/       -10.000010D0/
      data      zs8( 56)/         1.925822D0/
      data      zp8( 56)/         1.451991D0/
      data     alp8( 56)/         0.499997D0/
      data     gss8( 56)/         4.837241D0/
      data     gsp8( 56)/         3.194222D0/
      data     gpp8( 56)/         2.124328D0/
      data     gp28( 56)/         2.215394D0/
      data     hsp8( 56)/         0.399924D0/
!
!                    Data for Element  78         Platinum
!
      data     udd8( 78)/      -111.371876D0/
      data   betad8( 78)/       -10.325035D0/
      data     zsn8( 78)/         1.896554D0/
      data     zpn8( 78)/         1.668468D0/
      data     zdn8( 78)/         2.736124D0/
      data     gss8( 78)/         7.214948D0/
      data     gsp8( 78)/         6.717575D0/
      data     gpp8( 78)/         6.963140D0/
      data     gp28( 78)/         6.039313D0/
      data     hsp8( 78)/         1.634205D0/
      data    f0sd8( 78)/         8.380502D0/
      data    g2sd8( 78)/         2.420265D0/
!
!                    Data for Element  80          Mercury
!
      data     uss8( 80)/       -17.762229D0/
      data     upp8( 80)/       -18.330751D0/
      data   betas8( 80)/        -3.101365D0/
      data   betap8( 80)/        -3.464031D0/
      data      zs8( 80)/         1.476885D0/
      data      zp8( 80)/         2.479951D0/
      data     alp8( 80)/         1.529377D0/
      data     gss8( 80)/         6.624720D0/
      data     gsp8( 80)/        10.639297D0/
      data     gpp8( 80)/        14.709283D0/
      data     gp28( 80)/        16.000740D0/
      data     hsp8( 80)/         2.036311D0/
      data gues81( 80,1)/         1.082720D0/
      data gues82( 80,1)/         6.496598D0/
      data gues83( 80,1)/         1.195146D0/
      data gues81( 80,2)/        -0.096553D0/
      data gues82( 80,2)/         3.926281D0/
      data gues83( 80,2)/         2.627160D0/
!
!                    Data for Element  81         Thallium
!
      data     uss8( 81)/       -30.053170D0/
      data     upp8( 81)/       -26.920637D0/
      data   betas8( 81)/        -1.084495D0/
      data   betap8( 81)/        -7.946799D0/
      data      zs8( 81)/         6.867921D0/
      data      zp8( 81)/         1.969445D0/
      data     alp8( 81)/         1.340951D0/
      data     gss8( 81)/        10.460412D0/
      data     gsp8( 81)/        11.223883D0/
      data     gpp8( 81)/         4.992785D0/
      data     gp28( 81)/         8.962727D0/
      data     hsp8( 81)/         2.530406D0/
      data gues81( 81,1)/        -1.361399D0/
      data gues82( 81,1)/         3.557226D0/
      data gues83( 81,1)/         1.092802D0/
      data gues81( 81,2)/        -0.045401D0/
      data gues82( 81,2)/         2.306995D0/
      data gues83( 81,2)/         2.965029D0/
!
!                    Data for Element  82             Lead
!
      data     uss8( 82)/       -30.322756D0/
      data     upp8( 82)/       -24.425834D0/
      data   betas8( 82)/        -6.126024D0/
      data   betap8( 82)/        -1.395430D0/
      data      zs8( 82)/         3.141289D0/
      data      zp8( 82)/         1.892418D0/
      data     alp8( 82)/         1.620045D0/
      data     gss8( 82)/         7.011992D0/
      data     gsp8( 82)/         6.793782D0/
      data     gpp8( 82)/         5.183780D0/
      data     gp28( 82)/         5.045651D0/
      data     hsp8( 82)/         1.566302D0/
      data gues81( 82,1)/        -0.122576D0/
      data gues82( 82,1)/         6.003062D0/
      data gues83( 82,1)/         1.901597D0/
      data gues81( 82,2)/        -0.056648D0/
      data gues82( 82,2)/         4.743705D0/
      data gues83( 82,2)/         2.861879D0/
!
!                    Data for Element  83          Bismuth
!
      data     uss8( 83)/       -33.495938D0/
      data     upp8( 83)/       -35.521026D0/
      data   betas8( 83)/        -5.607283D0/
      data   betap8( 83)/        -5.800152D0/
      data      zs8( 83)/         4.916451D0/
      data      zp8( 83)/         1.934935D0/
      data     alp8( 83)/         1.857431D0/
      data     gss8( 83)/         4.989480D0/
      data     gsp8( 83)/         6.103308D0/
      data     gpp8( 83)/         8.696007D0/
      data     gp28( 83)/         8.335447D0/
      data     hsp8( 83)/         0.599122D0/
      data gues81( 83,1)/         2.581693D0/
      data gues82( 83,1)/         5.094022D0/
      data gues83( 83,1)/         0.499787D0/
      data gues81( 83,2)/         0.060320D0/
      data gues82( 83,2)/         6.001538D0/
      data gues83( 83,2)/         2.427970D0/
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
!                    Data for Element 102      Capped bond
!
      data     uss8(102)/       -11.906276D0/
      data   betas8(102)/  -9999999.000000D0/
      data      zs8(102)/         4.000000D0/
      data      zp8(102)/         0.300000D0/
      data     alp8(102)/         2.544134D0/
      data     gss8(102)/        12.848000D0/
      data     hsp8(102)/         0.100000D0/
!
!                    Data for Element 104        + Sparkle
!
      data     alp8(104)/         1.500000D0/
!
!                    Data for Element 106        - Sparkle
!
      data     alp8(106)/         1.500000D0/
!
!
!                     Global parameters
!
!
      data   v_par8(1)/    8.947612d0/  ! Used in ccrep for scalar correction of C-C triple bonds.
      data   v_par8(2)/    6.024265d0/  ! Used in ccrep for exponent correction of C-C triple bonds.
      data   v_par8(3)/   -0.012037d0/  ! Used in ccrep for scalar correction of O-H term.
      data   v_par8(4)/    0.701333d0/  ! Used in ccrep for exponent correction of C-C triple bonds.
      data   v_par8(7)/    0.900000d0/  ! Used in dftd3 to set "s6"  in D3H4
      data   v_par8(8)/   14.000000d0/  ! Used in dftd3 to set "alp" in D3H4
      data   v_par8(9)/    1.561000d0/  ! Used in dftd3 to set "rs6" in D3H4
  contains
  subroutine alpb_and_xfac_pm8
    use parameters_C, only : xfac, alpb
 !
      alpb(11, 1) =     1.800472d0 !      Sodium -     Hydrogen
      xfac(11, 1) =     3.171946d0 !      Sodium -     Hydrogen
      alpb(11, 6) =     1.321600d0 !      Sodium -       Carbon
      xfac(11, 6) =     0.876072d0 !      Sodium -       Carbon
      alpb(11, 7) =     0.999895d0 !      Sodium -     Nitrogen
      xfac(11, 7) =     0.295812d0 !      Sodium -     Nitrogen
      alpb(11, 8) =     2.116028d0 !      Sodium -       Oxygen
      xfac(11, 8) =     3.392634d0 !      Sodium -       Oxygen
      alpb(11, 9) =     1.860719d0 !      Sodium -     Fluorine
      xfac(11, 9) =     1.605474d0 !      Sodium -     Fluorine
      alpb(11,11) =     1.024150d0 !      Sodium -       Sodium
      xfac(11,11) =     1.126080d0 !      Sodium -       Sodium
 !
      alpb(16,11) =     0.999990d0 !      Sulfur -       Sodium
      xfac(16,11) =     0.241710d0 !      Sulfur -       Sodium
 !
      alpb(17,11) =     1.420756d0 !    Chlorine -       Sodium
      xfac(17,11) =     1.168589d0 !    Chlorine -       Sodium
 !
      alpb(19, 1) =     1.027701d0 !   Potassium -     Hydrogen
      xfac(19, 1) =     0.893642d0 !   Potassium -     Hydrogen
      alpb(19, 6) =     1.453975d0 !   Potassium -       Carbon
      xfac(19, 6) =     2.009100d0 !   Potassium -       Carbon
      alpb(19, 7) =     0.999984d0 !   Potassium -     Nitrogen
      xfac(19, 7) =     0.425244d0 !   Potassium -     Nitrogen
      alpb(19, 8) =     1.470651d0 !   Potassium -       Oxygen
      xfac(19, 8) =     0.999945d0 !   Potassium -       Oxygen
      alpb(19, 9) =     0.999898d0 !   Potassium -     Fluorine
      xfac(19, 9) =     0.342063d0 !   Potassium -     Fluorine
      alpb(19,16) =     1.000576d0 !   Potassium -       Sulfur
      xfac(19,16) =     0.584513d0 !   Potassium -       Sulfur
      alpb(19,17) =     0.999989d0 !   Potassium -     Chlorine
      xfac(19,17) =     0.584823d0 !   Potassium -     Chlorine
      alpb(19,19) =     1.532875d0 !   Potassium -    Potassium
      xfac(19,19) =     8.250024d0 !   Potassium -    Potassium
 !
      alpb(20, 1) =     1.427846d0 !     Calcium -     Hydrogen
      xfac(20, 1) =     1.805489d0 !     Calcium -     Hydrogen
      alpb(20, 6) =     0.999993d0 !     Calcium -       Carbon
      xfac(20, 6) =     0.320183d0 !     Calcium -       Carbon
      alpb(20, 7) =     0.999993d0 !     Calcium -     Nitrogen
      xfac(20, 7) =     0.275261d0 !     Calcium -     Nitrogen
      alpb(20, 8) =     1.934512d0 !     Calcium -       Oxygen
      xfac(20, 8) =     0.480203d0 !     Calcium -       Oxygen
      alpb(20, 9) =     1.964914d0 !     Calcium -     Fluorine
      xfac(20, 9) =     0.594612d0 !     Calcium -     Fluorine
      alpb(20,16) =     0.999960d0 !     Calcium -       Sulfur
      xfac(20,16) =     0.212682d0 !     Calcium -       Sulfur
      alpb(20,17) =     1.669363d0 !     Calcium -     Chlorine
      xfac(20,17) =     0.928683d0 !     Calcium -     Chlorine
      alpb(20,20) =     0.999174d0 !     Calcium -      Calcium
      xfac(20,20) =     3.160373d0 !     Calcium -      Calcium
 !
      alpb(35,11) =     1.517965d0 !     Bromine -       Sodium
      xfac(35,11) =     1.741658d0 !     Bromine -       Sodium
      alpb(35,19) =     1.401604d0 !     Bromine -    Potassium
      xfac(35,19) =     2.440790d0 !     Bromine -    Potassium
      alpb(35,20) =     1.509690d0 !     Bromine -      Calcium
      xfac(35,20) =     0.874420d0 !     Bromine -      Calcium
 !
      alpb(37, 1) =     2.066936d0 !    Rubidium -     Hydrogen
      xfac(37, 1) =     9.999919d0 !    Rubidium -     Hydrogen
      alpb(37, 5) =     1.996668d0 !    Rubidium -        Boron
      xfac(37, 5) =    10.000119d0 !    Rubidium -        Boron
      alpb(37, 8) =     1.999948d0 !    Rubidium -       Oxygen
      xfac(37, 8) =     1.817233d0 !    Rubidium -       Oxygen
      alpb(37, 9) =     3.083855d0 !    Rubidium -     Fluorine
      xfac(37, 9) =    10.000006d0 !    Rubidium -     Fluorine
      alpb(37,17) =     2.371423d0 !    Rubidium -     Chlorine
      xfac(37,17) =     9.970607d0 !    Rubidium -     Chlorine
      alpb(37,35) =     2.071332d0 !    Rubidium -      Bromine
      xfac(37,35) =     7.407687d0 !    Rubidium -      Bromine
      alpb(37,37) =     0.539344d0 !    Rubidium -     Rubidium
      xfac(37,37) =     2.654922d0 !    Rubidium -     Rubidium
 !
      alpb(38, 1) =     1.385332d0 !   Strontium -     Hydrogen
      xfac(38, 1) =     7.195639d0 !   Strontium -     Hydrogen
      alpb(38, 6) =     1.392807d0 !   Strontium -       Carbon
      xfac(38, 6) =     5.455084d0 !   Strontium -       Carbon
      alpb(38, 7) =     1.392514d0 !   Strontium -     Nitrogen
      xfac(38, 7) =     5.293460d0 !   Strontium -     Nitrogen
      alpb(38, 8) =     2.471371d0 !   Strontium -       Oxygen
      xfac(38, 8) =     2.147711d0 !   Strontium -       Oxygen
      alpb(38, 9) =     2.525416d0 !   Strontium -     Fluorine
      xfac(38, 9) =     5.359983d0 !   Strontium -     Fluorine
      alpb(38,16) =     1.045944d0 !   Strontium -       Sulfur
      xfac(38,16) =     0.594774d0 !   Strontium -       Sulfur
      alpb(38,17) =     1.473507d0 !   Strontium -     Chlorine
      xfac(38,17) =     1.176805d0 !   Strontium -     Chlorine
      alpb(38,35) =     1.338861d0 !   Strontium -      Bromine
      xfac(38,35) =     1.082828d0 !   Strontium -      Bromine
      alpb(38,38) =     1.393785d0 !   Strontium -    Strontium
      xfac(38,38) =     5.138251d0 !   Strontium -    Strontium
 !
      alpb(53,11) =     1.148784d0 !      Iodine -       Sodium
      xfac(53,11) =     0.747959d0 !      Iodine -       Sodium
      alpb(53,19) =     1.402113d0 !      Iodine -    Potassium
      xfac(53,19) =     2.568324d0 !      Iodine -    Potassium
      alpb(53,20) =     1.513600d0 !      Iodine -      Calcium
      xfac(53,20) =     0.856865d0 !      Iodine -      Calcium
      alpb(53,38) =     1.730393d0 !      Iodine -    Strontium
      xfac(53,38) =     3.855491d0 !      Iodine -    Strontium
 !
      alpb(55, 1) =     1.448617d0 !      Cesium -     Hydrogen
      xfac(55, 1) =     6.721688d0 !      Cesium -     Hydrogen
      alpb(55, 6) =     0.999871d0 !      Cesium -       Carbon
      xfac(55, 6) =     5.579858d0 !      Cesium -       Carbon
      alpb(55, 7) =     0.999559d0 !      Cesium -     Nitrogen
      xfac(55, 7) =     0.242498d0 !      Cesium -     Nitrogen
      alpb(55, 8) =     3.222538d0 !      Cesium -       Oxygen
      xfac(55, 8) =     7.547802d0 !      Cesium -       Oxygen
      alpb(55, 9) =     2.748652d0 !      Cesium -     Fluorine
      xfac(55, 9) =     3.360596d0 !      Cesium -     Fluorine
      alpb(55,15) =     0.999863d0 !      Cesium -   Phosphorus
      xfac(55,15) =     0.099894d0 !      Cesium -   Phosphorus
      alpb(55,16) =     0.999607d0 !      Cesium -       Sulfur
      xfac(55,16) =     0.869429d0 !      Cesium -       Sulfur
      alpb(55,17) =     0.999578d0 !      Cesium -     Chlorine
      xfac(55,17) =     0.288955d0 !      Cesium -     Chlorine
      alpb(55,35) =     0.999651d0 !      Cesium -      Bromine
      xfac(55,35) =     0.350614d0 !      Cesium -      Bromine
      alpb(55,53) =     0.998419d0 !      Cesium -       Iodine
      xfac(55,53) =     0.377054d0 !      Cesium -       Iodine
      alpb(55,55) =     0.997434d0 !      Cesium -       Cesium
      xfac(55,55) =     8.549732d0 !      Cesium -       Cesium
 !
      alpb(56, 1) =     0.999997d0 !      Barium -     Hydrogen
      xfac(56, 1) =     0.099997d0 !      Barium -     Hydrogen
      alpb(56, 6) =     0.999997d0 !      Barium -       Carbon
      xfac(56, 6) =     0.099997d0 !      Barium -       Carbon
      alpb(56, 7) =     0.999997d0 !      Barium -     Nitrogen
      xfac(56, 7) =     0.099997d0 !      Barium -     Nitrogen
      alpb(56, 8) =     1.249857d0 !      Barium -       Oxygen
      xfac(56, 8) =     0.352542d0 !      Barium -       Oxygen
      alpb(56, 9) =     1.886689d0 !      Barium -     Fluorine
      xfac(56, 9) =     1.528347d0 !      Barium -     Fluorine
      alpb(56,16) =     0.999862d0 !      Barium -       Sulfur
      xfac(56,16) =     0.904098d0 !      Barium -       Sulfur
      alpb(56,17) =     1.490059d0 !      Barium -     Chlorine
      xfac(56,17) =     1.846584d0 !      Barium -     Chlorine
      alpb(56,35) =     1.628325d0 !      Barium -      Bromine
      xfac(56,35) =     3.652542d0 !      Barium -      Bromine
      alpb(56,53) =     1.461169d0 !      Barium -       Iodine
      xfac(56,53) =     2.514512d0 !      Barium -       Iodine
      alpb(56,56) =     0.999997d0 !      Barium -       Barium
      xfac(56,56) =     0.099997d0 !      Barium -       Barium
    end subroutine alpb_and_xfac_pm8
  end module Parameters_for_PM8_C
