  subroutine update(iparam, ielmnt, param, c1)
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!
    use parameters_C, only : natorb, guess1, guess2, guess3, zs, zp, zd, &
    betas, betap, betad, alp, zsn, zpn, zdn, uss, upp, udd, gss, gpp, &
    gsp, gp2, hsp, pocord, alpb, xfac, f0sd_store, g2sd_store, dorbs, v_par, &
    f0sd, g2sd, CPE_Zeta, CPE_Z0, CPE_B, CPE_Xlo, CPE_Xhi, &
    n_partyp_fn, n_partyp_alpb, n_partyp
      use chanel_C, only : iw
!...Translated by Pacific-Sierra Research 77to90  4.4G  11:05:04  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: iparam 
      integer , intent(in) :: ielmnt 
      real(double) , intent(in) :: c1, param 
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
      betas(ielmnt) = betas(ielmnt) * c1 + param
    case (8)
      betap(ielmnt) = betap(ielmnt) * c1 + param
    case (9)
      betad(ielmnt) = betad(ielmnt) * c1 + param
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
    case (n_partyp)
      natorb(ielmnt) = Nint (param)
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
    case default
      uss(ielmnt) = uss(ielmnt) * c1 + param
    end select
    return
  end subroutine update 
