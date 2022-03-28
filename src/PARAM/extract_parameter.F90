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

 subroutine extract_parameter (iparam, ielmnt, param)
!
    use parameters_C, only : guess1, guess2, guess3, zs, zp, zd, &
    betas, betap, betad, alp, zsn, zpn, zdn, uss, upp, udd, gss, gpp, &
    gsp, gp2, hsp, pocord, xfac, alpb, f0sd_store, g2sd_store, v_par, &
    n_partyp_fn, n_partyp_alpb
    implicit none
    integer, intent (in) :: iparam
    integer, intent (in) :: ielmnt
    double precision, intent (inout) :: param
!------------------------------------------------------------
    integer :: kfn, ni, nj, jparam
!------------------------------------------------------------
    jparam = iparam
    kfn = 0
    if (jparam > n_partyp_fn - 1 .and. jparam < n_partyp_alpb) then
      kfn = (jparam - n_partyp_fn) / 3
      jparam = jparam - kfn * 3
      kfn = kfn + 1
    end if
    select case (jparam)
    case (2)
      param = upp(ielmnt)
    case (3)
      param = udd(ielmnt)
    case (4)
      param = zs(ielmnt)
    case (5)
      param = zp(ielmnt)
    case (6)
      param = zd(ielmnt)
    case (7)
      param = betas(ielmnt)
    case (8)
      param = betap(ielmnt)
    case (9)
      param = betad(ielmnt)
    case (10)
      param = gss(ielmnt)
    case (11)
      param = gsp(ielmnt)
    case (12)
      param = gpp(ielmnt)
    case (13)
      param = gp2(ielmnt)
    case (14)
      param = hsp(ielmnt)
    case (15)
      param = f0sd_store(ielmnt)
    case (16)
      param = g2sd_store(ielmnt)
    case (17)
      param = pocord(ielmnt)
    case (18)
      param = alp(ielmnt)
    case (19)
      param = zsn(ielmnt)
    case (20)
      param = zpn(ielmnt)
    case (21)
      param = zdn(ielmnt)
!
!  Add in C_ZET to C_XHI here
!
    case (27)
      param = guess1(ielmnt, kfn)
    case (28)
      param = guess2(ielmnt, kfn)
    case (29)
      param = guess3(ielmnt, kfn)
    case (25)
      write (6, "(' YOU ARE NOT ALLOWED TO OPTIMIZE THIS PARAMETER!')")
      stop
    case (39)
      nj = ielmnt/200
      ni = ielmnt - nj*200
      param = alpb(ni, nj)
    case (40)
      nj = ielmnt/200
      ni = ielmnt - nj*200
      param = xfac(ni, nj)
    case (41)
      param = v_par(ielmnt)
    case default
      param = uss(ielmnt)
    end select
    return
end subroutine extract_parameter
