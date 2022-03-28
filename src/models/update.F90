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

  subroutine update(iparam, ielmnt, param, c1)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!
    use parameters_C, only : natorb, guess1, guess2, guess3, zs, zp, zd, &
    betas, betap, betad, alp, zsn, zpn, zdn, uss, upp, udd, gss, gpp, &
    gsp, gp2, hsp, pocord, alpb, xfac, f0sd_store, g2sd_store, dorbs, v_par, &
    f0sd, g2sd, CPE_Zeta, CPE_Z0, CPE_B, CPE_Xlo, CPE_Xhi, &
    n_partyp_fn, n_partyp_alpb
      use chanel_C, only : iw
      use molkst_C, only : method_indo
      use reimers_C, only : isok, nbfa, zcorea, zeta, zetad, zetawt, betaa, fg
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: iparam
      integer , intent(in) :: ielmnt
      double precision , intent(in) :: c1, param
!-----------------------------------------------
!***********************************************************************
!
!  UPDATE UPDATES THE MODULES WHICH HOLD ALL THE PARAMETERS FOR
!         RUNNING MOPAC.
!         IPARAM REFERS TO THE TYPE OF PARAMETER,
!         IELMNT REFERS TO THE ELEMENT,
!         PARAM IS THE VALUE OF THE PARAMETER.
!         C1 is zero if the parameter is to be re-set to PARAM
!         C1 is one if the parameter is to be perturbed
!
!***********************************************************************
!------------------------------------------------------------
    integer :: i, kfn, ni, nj, jparam
    intrinsic Nint
!------------------------------------------------------------
    kfn = 0
    jparam = iparam
    if (jparam > n_partyp_fn - 1 .and. jparam < n_partyp_alpb) then
      kfn = (jparam - n_partyp_fn) / 3
      jparam = jparam - kfn * 3
      kfn = kfn + 1
    end if
    select case (jparam)
    case (2)
      upp(ielmnt) = upp(ielmnt) * c1 + param
    case (3)
      udd(ielmnt) = udd(ielmnt) * c1 + param
    case (4)
      zs(ielmnt) = zs(ielmnt) * c1 + param
    case (5)
      zp(ielmnt) = zp(ielmnt) * c1 + param
    case (6)
      zd(ielmnt) = zd(ielmnt) * c1 + param
    case (7)
!
!  Add INDO switch for beta parameters
!
      if (method_indo) then
        betaa(1, ielmnt) =  betaa(1, ielmnt) * c1 + param
      else
        betas(ielmnt) = betas(ielmnt) * c1 + param
      end if
    case (8)
      if (method_indo) then
        betaa(2, ielmnt) =  betaa(2, ielmnt) * c1 + param
      else
        betap(ielmnt) = betap(ielmnt) * c1 + param
      end if
    case (9)
      if (method_indo) then
        betaa(3, ielmnt) =  betaa(3, ielmnt) * c1 + param
      else
        betad(ielmnt) = betad(ielmnt) * c1 + param
      end if
    case (10)
      gss(ielmnt) = gss(ielmnt) * c1 + param
    case (11)
      gsp(ielmnt) = gsp(ielmnt) * c1 + param
    case (12)
      gpp(ielmnt) = gpp(ielmnt) * c1 + param
    case (13)
      gp2(ielmnt) = gp2(ielmnt) * c1 + param
    case (14)
      hsp(ielmnt) = hsp(ielmnt) * c1 + param
    case (15)
      f0sd_store(ielmnt) = f0sd_store(ielmnt) * c1 + param
      if (c1 < 1.d-20) f0sd(ielmnt) = param
    case (16)
      g2sd_store(ielmnt) = g2sd_store(ielmnt) * c1 + param
      if (c1 < 1.d-20) g2sd(ielmnt) = param
    case (17)
      pocord(ielmnt) = pocord(ielmnt) * c1 + param
    case (18)
      alp(ielmnt) = alp(ielmnt) * c1 + param
    case (19)
      zsn(ielmnt) = zsn(ielmnt) * c1 + param
    case (20)
      zpn(ielmnt) = zpn(ielmnt) * c1 + param
    case (21)
      zdn(ielmnt) = zdn(ielmnt) * c1 + param
    case (22)
      CPE_Zeta(ielmnt) = CPE_Zeta(ielmnt) * c1 + param
    case (23)
      CPE_Z0(ielmnt) = CPE_Z0(ielmnt) * c1 + param
    case (24)
      CPE_B(ielmnt) = CPE_B(ielmnt) * c1 + param
    case (25)
      CPE_Xlo(ielmnt) = CPE_Xlo(ielmnt) * c1 + param
    case (26)
      CPE_Xhi(ielmnt) = CPE_Xhi(ielmnt) * c1 + param
    case (27)
      guess1(ielmnt, kfn) = guess1(ielmnt, kfn) * c1 + param
    case (28)
      guess2(ielmnt, kfn) = guess2(ielmnt, kfn) * c1 + param
    case (29)
      guess3(ielmnt, kfn) = guess3(ielmnt, kfn) * c1 + param
    case (42)
      if (method_indo) then
        nbfa(ielmnt) = Nint (param)
      else
        natorb(ielmnt) = Nint (param)
      end if
      dorbs(ielmnt) = (natorb(ielmnt) == 9)
      i = Nint (param)
      if (i /= 9 .and. i /= 4 .and. i /= 1) then
        write (iw, "(///10x,' UNACCEPTABLE VALUE FOR NO. OF ORBITALS ON ATOM ')")
        stop
      end if
    case (n_partyp_alpb)
      nj = ielmnt/200
      ni = ielmnt - nj*200
      alpb(ni,nj) = alpb(ni,nj)*c1 + param
      alpb(nj,ni) = alpb(ni,nj)
    case (n_partyp_alpb + 1)
      nj = ielmnt/200
      ni = ielmnt - nj*200
      xfac(ni,nj) = xfac(ni,nj)*c1 + param
      xfac(nj,ni) = xfac(ni,nj)
    case (41)
      v_par(ielmnt) = v_par(ielmnt) * c1 + param
!
!  Add INDO parameters
!
    case (43)
      zcorea(ielmnt) = zcorea(ielmnt) * c1 + param
    case (44)
      zeta(ielmnt)  = zeta(ielmnt) * c1 + param
    case (45)
      zetad(1,ielmnt) = zetad(1,ielmnt) * c1 + param
    case (46)
      zetad(2,ielmnt) = zetad(2,ielmnt) * c1 + param
    case (47)
      zetawt(1,ielmnt) = zetawt(1,ielmnt) * c1 + param
    case (48)
      zetawt(2,ielmnt) = zetawt(2,ielmnt) * c1 + param
    case (49)
      fg( 1,ielmnt) = fg( 1,ielmnt) * c1 + param
    case (50)
      fg( 2,ielmnt) = fg( 2,ielmnt) * c1 + param
    case (51)
      fg( 3,ielmnt) = fg( 3,ielmnt) * c1 + param
    case (52)
      fg( 4,ielmnt) = fg( 4,ielmnt) * c1 + param
    case (53)
      fg( 5,ielmnt) = fg( 5,ielmnt) * c1 + param
    case (54)
      fg( 6,ielmnt) = fg( 6,ielmnt) * c1 + param
    case (55)
      fg( 7,ielmnt) = fg( 7,ielmnt) * c1 + param
    case (56)
      fg( 8,ielmnt) = fg( 8,ielmnt) * c1 + param
    case (57)
      fg( 9,ielmnt) = fg( 9,ielmnt) * c1 + param
    case (58)
      fg(10,ielmnt) = fg(10,ielmnt) * c1 + param
    case (59)
      fg(11,ielmnt) = fg(11,ielmnt) * c1 + param
    case (60)
      fg(12,ielmnt) = fg(12,ielmnt) * c1 + param
    case (61)
      fg(13,ielmnt) = fg(13,ielmnt) * c1 + param
    case (62)
      fg(14,ielmnt) = fg(14,ielmnt) * c1 + param
    case (63)
      fg(15,ielmnt) = fg(15,ielmnt) * c1 + param
    case (64)
      fg(16,ielmnt) = fg(16,ielmnt) * c1 + param
    case (65)
      fg(17,ielmnt) = fg(17,ielmnt) * c1 + param
    case (66)
      fg(18,ielmnt) = fg(18,ielmnt) * c1 + param
    case (67)
      fg(19,ielmnt) = fg(19,ielmnt) * c1 + param
    case (68)
      fg(20,ielmnt) = fg(20,ielmnt) * c1 + param
    case (69)
      fg(21,ielmnt) = fg(21,ielmnt) * c1 + param
    case (70)
      fg(22,ielmnt) = fg(22,ielmnt) * c1 + param
    case (71)
      fg(23,ielmnt) = fg(23,ielmnt) * c1 + param
    case (72)
      fg(24,ielmnt) = fg(24,ielmnt) * c1 + param
    case default
      uss(ielmnt) = uss(ielmnt) * c1 + param
    end select
!
! If parameters are read in for an element, assume element is OK to use.
!
    if (method_indo) isok(ielmnt) = 1
  return
  end subroutine update
