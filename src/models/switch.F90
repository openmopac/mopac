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

      subroutine switch
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE parameters_C, only : alp, guess1, guess2, guess3, &
      betas, betap, betad, uss, upp, udd, zs, zp, zd, zsn, zpn, zdn, &
      gss, gsp, gpp, gp2, hsp, f0sd, g2sd, f0sd_store, g2sd_store, v_par, &
      polvol, pocord, xfac, alpb, CPE_Zeta, CPE_Z0, CPE_B, CPE_Xlo, CPE_Xhi, &
      ams, npq, natorb
!
!
      USE molkst_C, only : keywrd, keywrd_quoted, method_mndo, method_pm3, method_indo, &
      & method_mndod, method_pm6, method_rm1, method_pm7, method_PM7_ts, method_pm6_org, method_pm8
!
!
      USE parameters_for_INDO_C, only: isoki, nbfai, zetai, zetadi, &
        zetawti, zcoreai, betaai, fgi

      USE reimers_C, only: isok, nprin, nbfa, zcorea, weight,&
        zeta, zetad, zetawt, betaa, fg, mg1sp, ifgfac
!
!
      USE parameters_for_mndod_C, only : zsnd, zpnd, zdnd, ussd, uppd, uddd, &
      & zsd, zpd, zdd, betasd, betapd, betadd, alpd, gssd, gspd, &
      & gp2d, hspd, gppd, poc_d
!
!
      USE parameters_for_mndo_C, only :  ussm, uppm, uddm, zsm, zpm, zdm, betasm, &
        betapm, betadm, alpm, gssm, gspm, gppm, gp2m, hspm, polvolm, guesm1, &
        guesm2, guesm3, f0sdm, g2sdm, pocm, zsnm, zpnm, zdnm
!
!
      USE parameters_for_PM6_C, only :  uss6, upp6, udd6, zs6, zp6, zd6, betas6, &
        betap6, betad6, gss6, gsp6, gpp6, gp26, hsp6, polvo6, gues61, f0sd6, g2sd6, &
        gues62, gues63, poc_6, zsn6, zpn6, zdn6, alpb_and_xfac_pm6, alp6, v_par6, &
        CPE_Zet6, CPE_Z06, CPE_B6, CPE_Xlo6, CPE_Xhi6
!
!
      USE parameters_for_PM7_C, only :  uss7, upp7, udd7, zs7, zp7, zd7, betas7, &
        betap7, betad7, gss7, gsp7, gpp7, gp27, hsp7, polvo7, gues71, f0sd7, g2sd7, &
        gues72, gues73, poc_7, zsn7, zpn7, zdn7, alpb_and_xfac_pm7, alp7, v_par7
!
!
        USE parameters_for_PM6_ORG_C, only :  uss_org, upp_org, udd_org, zs_org, zp_org, zd_org, betas_org, &
        betap_org, betad_org, gss_org, gsp_org, gpp_org, gp2_org, hsp_org, polvo_org, gues_org1, f0sd_org, g2sd_org, &
        gues_org2, gues_org3, poc__org, zsn_org, zpn_org, zdn_org, alpb_and_xfac_pm6_org, alp_org, v_par_org
!
!
      USE parameters_for_PM8_C, only :  uss8, upp8, udd8, zs8, zp8, zd8, betas8, &
        betap8, betad8, gss8, gsp8, gpp8, gp28, hsp8, polvo8, gues81, f0sd8, g2sd8, &
        gues82, gues83, poc_8, zsn8, zpn8, zdn8, alpb_and_xfac_pm8, alp8, v_par8
!
!
       USE parameters_for_PM7_TS_C, only : uss7_TS, upp7_TS, udd7_TS, zs7_TS, zp7_TS, zd7_TS, betas7_TS, &
         betap7_TS, betad7_TS, gss7_TS, gsp7_TS, gpp7_TS, gp27_TS, hsp7_TS, polvo7_TS, poc_7_TS, &
         zsn7_TS, zpn7_TS, zdn7_TS, f0sd7_TS, g2sd7_TS, alp7_TS, v_par7_TS, alpb_and_xfac_pm7_ts, &
         gues7_ts1, gues7_ts2, gues7_ts3
!
!
      USE parameters_for_PM3_C, only : usspm3, upppm3, zspm3, zppm3, &
        betasp, betapp, alppm3, gsspm3, gsppm3, gpppm3, polvolpm3, &
        gp2pm3, hsppm3, guesp1, guesp2, guesp3
!
!
      use parameters_for_PM3_Sparkles_C, only : gssPM3sp, alpPM3sp, guesPM3sp1, guesPM3sp2, guesPM3sp3
!
!
      USE parameters_for_AM1_C, only : zsam1, zpam1, zdam1, ussam1, uppam1, uddam1, alpam1, &
      & gssam1, gppam1, gspam1, gp2am1, betasa, betapa, hspam1, guesa1, &
      & guesa2, guesa3, betada, zsnam1, zpnam1, zdnam1, f0sdam1, g2sdam1, polvolam1
!
!
      use parameters_for_AM1_Sparkles_C, only : gssam1sp, alpam1sp, guesam1sp1, guesam1sp2, guesam1sp3
!
!
      use parameters_for_PM6_Sparkles_C, only : gss6sp, alp6sp, gues6sp1, gues6sp2, gues6sp3
!
!
      Use Parameters_for_PM7_Sparkles_C, only : gss7sp, alp7sp, gues7sp1, gues7sp2, gues7sp3
!
!
      Use Parameters_for_RM1_Sparkles_C, only : gssrm1sp, alprm1sp, guesrm1sp1, guesrm1sp2, guesrm1sp3
!
!
      use parameters_for_RM1_C, only : ussRM1, uppRM1, uddRM1, zsRM1, zpRM1, zdRM1, betasRM1, &
      & betapRM1, betadRM1, gssRM1, gspRM1, gppRM1, gp2RM1, hspRM1, alpRM1, zsnRM1, zpnRM1, zdnRM1, &
      & poc_RM1, f0sdRM1, g2sdRM1, guess1RM1, guess2RM1, guess3RM1
!
!
!***********************************************************************
!
!   SWITCH copies data from the reference modules into the
!          arrays used by MOPAC.  The choice of reference data used
!          is given by the keywords MNDO, AM1, PM3, or MNDOD
!
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j
        pocord = 0.d0
        if (method_mndo) then
!
!    SWITCH IN MNDO PARAMETERS
!
          guess1 = guesm1
          guess2 = guesm2
          guess3 = guesm3
          polvol = polvolm
          zs = zsm
          zp = zpm
          zd = zdm
          zsn = zsnm
          zpn = zpnm
          zdn = zdnm
          uss = ussm
          upp = uppm
          udd = uddm
          betas = betasm
          betap = betapm
          betad = betadm
          alp = alpm
          gss = gssm
          gpp = gppm
          gsp = gspm
          gp2 = gp2m
          hsp = hspm
          f0sd = f0sdm
          g2sd = g2sdm
          pocord = pocm
          call alpb_and_xfac_mndo
        else if (method_indo) then
!
!    SWITCH IN INDO PARAMETERS
!
          isok   = isoki
          nbfa   = nbfai
          zeta   = zetai
          zetad  = zetadi
          zetawt = zetawti
          zcorea = zcoreai
          betaa  = betaai
          fg     = fgi
          do i= 1,80
            natorb(i) = nbfa(i)
            nprin(i)  = npq(i,1)
            weight(i) = ams(i)
            do j= 1,11
              fg(mg1sp-1+j,i)= fg(mg1sp-1+j,i) / ifgfac(j)
            end do
            zs(i) = zeta(i)
            zp(i) = zeta(i)
          end do
        else if (method_pm3) then
!
!    SWITCH IN MNDO-PM3 PARAMETERS
!
          guess1 = guesp1
          guess2 = guesp2
          guess3 = guesp3
          polvol = polvolpm3
          zs = zspm3
          zp = zppm3
          zd = 0.d0
          uss = usspm3
          upp = upppm3
          betas = betasp
          betap = betapp
          alp = alppm3
          gss = gsspm3
          gpp = gpppm3
          gsp = gsppm3
          gp2 = gp2pm3
          hsp = hsppm3
          call alpb_and_xfac_pm3
          if (index(keywrd, " SPARK") /= 0) then
            do i = 1, 102
              if (alpPM3sp(i) > 0.1d0) then
                zd(i) = 0.d0
                zp(i) = 0.d0
                zs(i) = 0.d0
                zsn(i) = 0.d0
                zpn(i) = 0.d0
                zdn(i) = 0.d0
                uss(i) = 0.d0
                upp(i) = 0.d0
                udd(i) = 0.d0
                betas(i) = 0.d0
                betap(i) = 0.d0
                betad(i) = 0.d0
                alp(i) = alpPM3sp(i)
                gss(i) = gssPM3sp(i)
                gpp(i) = 0.d0
                gp2(i) = 0.d0
                hsp(i) = 0.d0
                gsp(i) = 0.d0
                f0sd(i) = 0.d0
                g2sd(i) = 0.d0
                pocord(i) = 0.d0
                guess1(i,1) = guesPM3sp1(i,1)
                guess2(i,1) = guesPM3sp2(i,1)
                guess3(i,1) = guesPM3sp3(i,1)
                guess1(i,2) = guesPM3sp1(i,2)
                guess2(i,2) = guesPM3sp2(i,2)
                guess3(i,2) = guesPM3sp3(i,2)
                guess1(i,3:4) = 0.d0
                guess2(i,3:4) = 0.d0
                guess3(i,3:4) = 0.d0
              end if
            end do
            alpb(:100,57:71) = 0.d0
            alpb(57:71,:100) = 0.d0
            xfac(:100,57:71) = 0.d0
            xfac(57:71,:100) = 0.d0
          end if
        else if (method_pm6_org) then
!
!    SWITCH IN PM6-ORG PARAMETERS
!
          guess1 = gues_org1
          guess2 = gues_org2
          guess3 = gues_org3
          zs = zs_org
          zp = zp_org
          zd = zd_org
          zsn = zsn_org
          zpn = zpn_org
          zdn = zdn_org
          uss = uss_org
          upp = upp_org
          udd = udd_org
          betas = betas_org
          betap = betap_org
          betad = betad_org
          gss = gss_org
          gpp = gpp_org
          gsp = gsp_org
          gp2 = gp2_org
          hsp = hsp_org
          f0sd = f0sd_org
          g2sd = g2sd_org
          alp = alp_org
          pocord = poc__org
          polvol = polvo_org
          v_par = v_par_org
          call alpb_and_xfac_pm6_org
        else if (method_pm8) then
!
!    SWITCH IN PM8 PARAMETERS
!
          guess1 = gues81
          guess2 = gues82
          guess3 = gues83
          zs = zs8
          zp = zp8
          zd = zd8
          zsn = zsn8
          zpn = zpn8
          zdn = zdn8
          uss = uss8
          upp = upp8
          udd = udd8
          betas = betas8
          betap = betap8
          betad = betad8
          gss = gss8
          gpp = gpp8
          gsp = gsp8
          gp2 = gp28
          hsp = hsp8
          f0sd = f0sd8
          g2sd = g2sd8
          alp = alp8
          pocord = poc_8
          polvol = polvo8
          v_par = v_par8
          call alpb_and_xfac_pm8
        else if (method_pm6) then
!
!    SWITCH IN PM6 PARAMETERS
!
          guess1 = gues61
          guess2 = gues62
          guess3 = gues63
          zs = zs6
          zp = zp6
          zd = zd6
          zsn = zsn6
          zpn = zpn6
          zdn = zdn6
          uss = uss6
          upp = upp6
          udd = udd6
          betas = betas6
          betap = betap6
          betad = betad6
          gss = gss6
          gpp = gpp6
          gsp = gsp6
          gp2 = gp26
          hsp = hsp6
          f0sd = f0sd6
          g2sd = g2sd6
          alp = alp6
          pocord = poc_6
          polvol = polvo6
          CPE_Zeta = CPE_Zet6
          CPE_Z0   = CPE_Z06
          CPE_B    = CPE_B6
          CPE_Xlo  = CPE_Xlo6
          CPE_Xhi  = CPE_Xhi6
          v_par = v_par6
          call alpb_and_xfac_pm6
          if (index(keywrd, " SPARK") /= 0) then
            do i = 1, 107
              if (alp6sp(i) > 0.1d0) then
                zd(i) = 0.d0
                zp(i) = 0.d0
                zs(i) = 0.d0
                zsn(i) = 0.d0
                zpn(i) = 0.d0
                zdn(i) = 0.d0
                uss(i) = 0.d0
                upp(i) = 0.d0
                udd(i) = 0.d0
                betas(i) = 0.d0
                betap(i) = 0.d0
                betad(i) = 0.d0
                alp(i) = alp6sp(i)
                gss(i) = gss6sp(i)
                gpp(i) = 0.d0
                gp2(i) = 0.d0
                hsp(i) = 0.d0
                gsp(i) = 0.d0
                f0sd(i) = 0.d0
                g2sd(i) = 0.d0
                pocord(i) = 0.d0
                guess1(i,1) = gues6sp1(i,1)
                guess2(i,1) = gues6sp2(i,1)
                guess3(i,1) = gues6sp3(i,1)
                guess1(i,2) = gues6sp1(i,2)
                guess2(i,2) = gues6sp2(i,2)
                guess3(i,2) = gues6sp3(i,2)
                guess1(i,3:4) = 0.d0
                guess2(i,3:4) = 0.d0
                guess3(i,3:4) = 0.d0
              end if
            end do
            alpb(:100,57:71) = 0.d0
            alpb(57:71,:100) = 0.d0
            xfac(:100,57:71) = 0.d0
            xfac(57:71,:100) = 0.d0
          end if
        else if (method_pm7_ts) then
!
!    SWITCH IN PM7_TS PARAMETERS
!
          guess1 = gues7_TS1
          guess2 = gues7_TS2
          guess3 = gues7_TS3
          zs = zs7_TS
          zp = zp7_TS
          zd = zd7_TS
          zsn = zsn7_TS
          zpn = zpn7_TS
          zdn = zdn7_TS
          uss = uss7_TS
          upp = upp7_TS
          udd = udd7_TS
          betas = betas7_TS
          betap = betap7_TS
          betad = betad7_TS
          gss = gss7_TS
          gpp = gpp7_TS
          gsp = gsp7_TS
          gp2 = gp27_TS
          hsp = hsp7_TS
          f0sd = f0sd7_TS
          g2sd = g2sd7_TS
          alp = alp7_TS
          pocord = poc_7_TS
          polvol = polvo7_TS
          v_par = v_par7_TS
          call alpb_and_xfac_pm7_TS
        else if (method_pm7) then
!
!    SWITCH IN PM7 PARAMETERS
!
          guess1 = gues71
          guess2 = gues72
          guess3 = gues73
          zs = zs7
          zp = zp7
          zd = zd7
          zsn = zsn7
          zpn = zpn7
          zdn = zdn7
          uss = uss7
          upp = upp7
          udd = udd7
          betas = betas7
          betap = betap7
          betad = betad7
          gss = gss7
          gpp = gpp7
          gsp = gsp7
          gp2 = gp27
          hsp = hsp7
          f0sd = f0sd7
          g2sd = g2sd7
          alp = alp7
          pocord = poc_7
          polvol = polvo7
          v_par = v_par7
          call alpb_and_xfac_pm7
          if (index(keywrd, " SPARK") /= 0) then
            do i = 1, 102
              if (alp7sp(i) > 0.1d0) then
                zd(i) = 0.d0
                zp(i) = 0.d0
                zs(i) = 0.d0
                zsn(i) = 0.d0
                zpn(i) = 0.d0
                zdn(i) = 0.d0
                uss(i) = 0.d0
                upp(i) = 0.d0
                udd(i) = 0.d0
                betas(i) = 0.d0
                betap(i) = 0.d0
                betad(i) = 0.d0
                alp(i) = alp7sp(i)
                gss(i) = gss7sp(i)
                gpp(i) = 0.d0
                gp2(i) = 0.d0
                hsp(i) = 0.d0
                gsp(i) = 0.d0
                f0sd(i) = 0.d0
                g2sd(i) = 0.d0
                pocord(i) = 0.d0
                guess1(i,1) = gues7sp1(i,1)
                guess2(i,1) = gues7sp2(i,1)
                guess3(i,1) = gues7sp3(i,1)
                guess1(i,2) = gues7sp1(i,2)
                guess2(i,2) = gues7sp2(i,2)
                guess3(i,2) = gues7sp3(i,2)
                guess1(i,3:4) = 0.d0
                guess2(i,3:4) = 0.d0
                guess3(i,3:4) = 0.d0
              end if
            end do
            alpb(:100,57:71) = 0.d0
            alpb(57:71,:100) = 0.d0
            xfac(:100,57:71) = 0.d0
            xfac(57:71,:100) = 0.d0
          end if
        else if (method_mndod) then
!
!    SWITCH IN MNDOD PARAMETERS
!
          uss = ussd
          upp = uppd
          udd = uddd
          zs  = zsd
          zp  = zpd
          zd  = zdd
          zsn = zsnd
          zpn = zpnd
          zdn = zdnd
          betas = betasd
          betap = betapd
          betad = betadd
          gss = gssd
          gsp = gspd
          gpp = gppd
          gp2 = gp2d
          hsp = hspd
          alp = alpd
          pocord = poc_d
          guess1 = 0.d0
          guess2 = 0.d0
          guess3 = 0.d0
          f0sd = 0.d0
          g2sd = 0.d0
          call alpb_and_xfac_mndod
        else if (method_rm1) then
!
!    SWITCH IN RM1 PARAMETERS
!
          guess1 = guess1RM1
          guess2 = guess2RM1
          guess3 = guess3RM1
          uss = ussRM1
          upp = uppRM1
          udd = uddRM1
          zs  = zsRM1
          zp  = zpRM1
          zd  = zdRM1
          zsn = zsnRM1
          zpn = zpnRM1
          zdn = zdnRM1
          betas = betasRM1
          betap = betapRM1
          betad = betadRM1
          gss = gssRM1
          gsp = gspRM1
          gpp = gppRM1
          gp2 = gp2RM1
          hsp = hspRM1
          alp = alpRM1
          f0sd = f0sdRM1
          g2sd = g2sdRM1
          alp = alpRM1
          pocord = poc_RM1
          if (index(keywrd, " SPARK") /= 0) then
            do i = 1, 102
              if (alp7sp(i) > 0.1d0) then
                zd(i) = 0.d0
                zp(i) = 0.d0
                zs(i) = 0.d0
                zsn(i) = 0.d0
                zpn(i) = 0.d0
                zdn(i) = 0.d0
                uss(i) = 0.d0
                upp(i) = 0.d0
                udd(i) = 0.d0
                betas(i) = 0.d0
                betap(i) = 0.d0
                betad(i) = 0.d0
                alp(i) = alprm1sp(i)
                gss(i) = gssrm1sp(i)
                gpp(i) = 0.d0
                gp2(i) = 0.d0
                hsp(i) = 0.d0
                gsp(i) = 0.d0
                f0sd(i) = 0.d0
                g2sd(i) = 0.d0
                pocord(i) = 0.d0
                guess1(i,1) = guesrm1sp1(i,1)
                guess2(i,1) = guesrm1sp2(i,1)
                guess3(i,1) = guesrm1sp3(i,1)
                guess1(i,2) = guesrm1sp1(i,2)
                guess2(i,2) = guesrm1sp2(i,2)
                guess3(i,2) = guesrm1sp3(i,2)
                guess1(i,3:4) = 0.d0
                guess2(i,3:4) = 0.d0
                guess3(i,3:4) = 0.d0
              end if
            end do
          end if
        else
!
!    SWITCH IN AM1 PARAMETERS
!
          guess1 = guesa1
          guess2 = guesa2
          guess3 = guesa3
          polvol = polvolam1
          zs = zsam1
          zp = zpam1
          zd = zdam1
          zsn = zsnam1
          zpn = zpnam1
          zdn = zdnam1
          f0sd = f0sdam1
          g2sd = g2sdam1
          uss = ussam1
          upp = uppam1
          udd = uddam1
          betas = betasa
          betap = betapa
          betad = betada
          alp = alpam1
          gss = gssam1
          gpp = gppam1
          gsp = gspam1
          gp2 = gp2am1
          hsp = hspam1
          call alpb_and_xfac_am1
!
!   Unique parameter for Voityuk's Molybdenum
!
          pocord(42) = 1.334d0
          if (index(keywrd, " SPARK") /= 0) then
            do i = 1, 102
              if (alpam1sp(i) > 0.1d0) then
                zd(i) = 0.d0
                zp(i) = 0.d0
                zs(i) = 0.d0
                zsn(i) = 0.d0
                zpn(i) = 0.d0
                zdn(i) = 0.d0
                uss(i) = 0.d0
                upp(i) = 0.d0
                udd(i) = 0.d0
                betas(i) = 0.d0
                betap(i) = 0.d0
                betad(i) = 0.d0
                alp(i) = alpam1sp(i)
                gss(i) = gssam1sp(i)
                gpp(i) = 0.d0
                gp2(i) = 0.d0
                hsp(i) = 0.d0
                gsp(i) = 0.d0
                f0sd(i) = 0.d0
                g2sd(i) = 0.d0
                pocord(i) = 0.d0
                guess1(i,1) = guesam1sp1(i,1)
                guess2(i,1) = guesam1sp2(i,1)
                guess3(i,1) = guesam1sp3(i,1)
                guess1(i,2) = guesam1sp1(i,2)
                guess2(i,2) = guesam1sp2(i,2)
                guess3(i,2) = guesam1sp3(i,2)
                guess1(i,3:4) = 0.d0
                guess2(i,3:4) = 0.d0
                guess3(i,3:4) = 0.d0
              end if
            end do
            alpb(:100,57:71) = 0.d0
            alpb(57:71,:100) = 0.d0
            xfac(:100,57:71) = 0.d0
            xfac(57:71,:100) = 0.d0
          end if
        end if
!
!
        f0sd_store = f0sd
        g2sd_store = g2sd

!
!  Symmetrize the alpb and xfac arrays
!
      do i = 1, 100
        alpb(:i,i) = alpb(i,:i)
        xfac(:i,i) = xfac(i,:i)
      end do
      call fractional_metal_ion
      if (index(keywrd_quoted,'EXTERNAL') /= 0) return
      if (.not. method_indo .and. uss(1) > (-1.D0)) then
        call mopend (&
          'THE HAMILTONIAN REQUESTED IS NOT AVAILABLE IN THIS PROGRAM')
        return
      end if
      return
      end subroutine switch
      subroutine fractional_metal_ion
      use parameters_C, only : tore, zs, zp, gss, gpp, gp2, gsp, hsp, &
      betas, betap, upp, alp
      integer :: i
!
!  Pseudo - halide metal ion with core charge of -1/2 - for use with
!  transition metal ions only
    tore(85) = -0.5d0
!
!  Pseudo - alkali metal ion with core charge of 1/2 - for use with
!  transition metal ions only
      tore(87) = 0.5d0
      do i = 85, 87, 2
        upp(i) = 0.d0
        alp(i) = 3.0d0
        zs(i) = 0.d0
        zp(i) = 0.d0
        betas(i) = 0.d0
        betap(i) = 0.d0
        gss(i) = 10.d0
        gsp(i) =0.d0
        gpp(i) = 0.d0
        gp2(i) = 0.d0
        hsp(i) = 0.d0
      end do
      end subroutine fractional_metal_ion
