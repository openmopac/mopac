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

    subroutine calpar
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE parameters_C, only : ios, iop, iod, qq, am, ad, aq, dd, &
      gpp, gp2, hsp, gss, gsp, &
      zs, zp, uss, upp, udd, eisol, &
      f0sd, g2sd, f0sd_store, g2sd_store, &
         dsd, dpd, ddd, zd
      USE funcon_C, only : ev
      USE molkst_C, only : keywrd, method_indo
      USE reimers_C, only: zetad, zetawt, nbfa
#if MOPAC_F2003
      USE, INTRINSIC :: IEEE_ARITHMETIC
#endif
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      integer , dimension(107) :: nspqn
      integer :: i, k, l, jmax, j
      double precision, dimension(107) :: gssc, gspc, hspc, gp2c, gppc
      double precision :: p, p4,  hpp, qn, gdd1, d1, d2, df, hsp1, hsp2, d3, &
        gqq, q1, q2, qf, hpp1, hpp2, q3
      double precision :: zda, zdb, zwa, zwb, dsda, dsdb, dpda, dpdb, ddda, dddb
      save nspqn
!-----------------------------------------------
      data nspqn/ 2*1, 8*2, 8*3, 18*4, 18*5, 32*6, 21*0/
!
      if (method_indo) then
        do i = 1, 107
#ifdef MOPAC_F2003
          if (ieee_is_nan(zs(i))) zs(i) = 0.d0
          if (ieee_is_nan(zp(i))) zp(i) = 0.d0
          if (ieee_is_nan(zd(i))) zd(i) = 1.d-6
#else
          if (isnan(zs(i))) zs(i) = 0.d0
          if (isnan(zp(i))) zp(i) = 0.d0
          if (isnan(zd(i))) zd(i) = 1.d-6
#endif
          if(zd(i) == 0.d0) zd(i) = 1.d-6
        end do
      end if
!
!
!     SET SCALING PARAMETER.
      p = 2.D0
      p4 = p**4
      call sp_two_electron
      am = 0.d0
      ad = 0.d0
      aq = 0.d0
      dd = 0.d0
      qq = 0.d0
      do i = 2, 97
!
!  GSSC is the number of two-electron terms of type <SS|SS>
!
        gssc(i) = max(ios(i)-1,0)
        k = iop(i)
!
!  GSPC is the number of two-electron terms of type <SS|PP>
!
        gspc(i) = ios(i)*k
        l = min(k,6 - k)
!
!  GP2C is the number of two-electron terms of type <PP|PP>
!       plus 0.5 of the number of HPP integrals.
!  (HPP is not used; instead it is replaced by 0.5(GPP-GP2))
!
        gp2c(i) = (k*(k - 1))/2 + 0.5D0*(l*(l - 1))/2
!
!  GPPC is minus 0.5 times the number of HPP integrals.
!
        gppc(i) = -0.5D0*(l*(l - 1))/2
!
!  HSPC is the number of two-electron terms of type <SP|SP>.
!       (S and P must have the same spin.  In all cases, if
!  P is non-zero, there are two S electrons)
!
        hspc(i) = -k*ios(i)*0.5d0
!
!
        if (zp(i)<1.D-4 .and. zs(i)<1.D-4) cycle
!*********************************************************************
!
!   CONSTRAINTS ON THE POSSIBLE VALUES OF PARAMETERS
!
!*********************************************************************
        zp(i) = dmax1(0.3D0,zp(i))
!  PUT IN ANY CONSTRAINTS AT THIS POINT
!*********************************************************************
        hpp = 0.5D0*(gpp(i)-gp2(i))
        hpp = max(0.1D0,hpp)
        hsp(i) = max(1.D-7,hsp(i))
        eisol(i) = uss(i)*ios(i) + upp(i)*iop(i) + udd(i)*iod(i) + gss(i)*gssc(i) + &
          gpp(i)*gppc(i) + gsp(i)*gspc(i) + gp2(i)*gp2c(i) + hsp(i)*hspc(i)
        qn = nspqn(i)
        dd(i) = (2.D0*qn + 1)*(4.D0*zs(i)*zp(i))**(qn + 0.5D0)/(zs(i)+zp(i))**(2.D0*qn + 2)/sqrt(3.D0)
        qq(i) = sqrt((4.D0*qn*qn + 6.D0*qn + 2.D0)/20.D0)/zp(i)
        if (method_indo .and. i<=80) then
          if (nbfa(i) > 4) then
!           INDO d orbitals - since two Slater basis functions, calculate two
!           terms and weights based on the coefficients
            zda = zetad(1,i)
            zdb = max(zetad(2,i), 1.d-8)
            zwa = zetawt(1,i)
            zwb = zetawt(2,i)

            dsda = (2.D0*zs(i))**(qn + 0.5D0) * (2.D0*zda)**(qn - 0.5D0) / (zs(i) + zda)**(2.D0*qn + 2)
            dsda = dsda * (2.D0*qn + 1) * sqrt(2.D0*qn*(2.D0*qn - 1))

            dsdb = (2.D0*zs(i))**(qn + 0.5D0) * (2.D0*zdb)**(qn - 0.5D0) / (zs(i) + zdb)**(2.D0*qn + 2)
            dsdb = dsdb * (2.D0*qn + 1) * sqrt(2.D0*qn*(2.D0*qn - 1))

            dsd(i) = sqrt((dsda*zwa + dsdb*zwb) / sqrt(15.D0))


            dpda = (2.D0*zp(i))**(qn + 0.5D0) * (2.D0*zda)**(qn - 0.5D0) / (zp(i) + zda)**(2.D0*qn + 1)
            dpda = dpda * sqrt(2.D0*qn*(2.D0*qn-1))

            dpdb = (2.D0*zp(i))**(qn + 0.5D0) * (2.D0*zdb)**(qn - 0.5D0) / (zp(i) + zdb)**(2.D0*qn + 1)
            dpdb = dpdb * sqrt(2.D0*qn*(2.D0*qn-1))

            dpd(i)= (dpda*zwa + dpdb*zwb) / sqrt(5.D0)


            ddda = (4.D0*qn*qn - 2.D0*qn) / (2*zda)**2
            dddb = (4.D0*qn*qn - 2.D0*qn) / (2*zdb)**2

            ddd(i)= sqrt((ddda*zwa + dddb*zwb) / 7.D0)
          else
! Single-Slater d orbitals
            dsd(i) = sqrt( (2.D0*zs(i))**(qn + 0.5D0) * (2.D0*zd(i))**(qn - 0.5D0) &
                     / (zs(i) + zd(i))**(2.D0*qn + 2) &
                     * (2.D0*qn + 1) * sqrt(2.D0*qn*(2.D0*qn - 1)) / sqrt(15.0))
            dpd(i) = (2.D0*zp(i))**(qn + 0.5D0) * (2.D0*zd(i))**(qn - 0.5D0) &
                     / (zp(i) + zd(i))**(2.D0*qn + 1) * sqrt(2.D0*qn*(2.D0*qn-1)) / sqrt(5.D0)
            ddd(i) = sqrt( (4.D0*qn*qn - 2.D0*qn) / (2*zd(i))**2 / 7.D0)
          end if
        end if
!     CALCULATE ADDITIVE TERMS, IN ATOMIC UNITS.
        jmax = 5
        gdd1 = (hsp(i)/(ev*dd(i)**2))**(1.D0/3.D0)
        d1 = gdd1
        d2 = gdd1 + 0.04D0
        do j = 1, jmax
          df = d2 - d1
          hsp1 = 0.5D0*d1 - 0.5D0/sqrt(4.D0*dd(i)**2+1.D0/d1**2)
          hsp2 = 0.5D0*d2 - 0.5D0/sqrt(4.D0*dd(i)**2+1.D0/d2**2)
          if (abs(hsp2 - hsp1) < 1.D-25) exit
          d3 = d1 + df*(hsp(i)/ev-hsp1)/(hsp2 - hsp1)
          d1 = d2
          d2 = d3
        end do
        gqq = (p4*hpp/(ev*48.D0*qq(i)**4))**0.2D0
        q1 = gqq
        q2 = gqq + 0.04D0
        do j = 1, jmax
          qf = q2 - q1
          hpp1 = 0.25D0*q1 - 0.5D0/sqrt(4.D0*qq(i)**2+1.D0/q1**2) + 0.25D0/&
            sqrt(8.D0*qq(i)**2+1.D0/q1**2)
          hpp2 = 0.25D0*q2 - 0.5D0/sqrt(4.D0*qq(i)**2+1.D0/q2**2) + 0.25D0/&
            sqrt(8.D0*qq(i)**2+1.D0/q2**2)
          if (abs(hpp2 - hpp1) < 1.D-25) exit
          q3 = q1 + qf*(hpp/ev - hpp1)/(hpp2 - hpp1)
          q1 = q2
          q2 = q3
        end do
        am(i) = gss(i)/ev
        ad(i) = d2
        aq(i) = q2
      end do
      do i = 1, 107
        if (am(i) < 1.d-20) then
          if (gss(i) > 1.d-20) then
            am(i) = gss(i)/ev
          else
            am(i) = 1.d0
          end if
        end if
      end do
      eisol(1) = uss(1)
      am(1) = gss(1)/ev
      ad(1) = am(1)
      aq(1) = am(1)
      do i = 1, 100
        if (f0sd_store(i) < 1.d-20) f0sd(i) = 0.d0 ! Force f0sd and g2sd to zero if not already defined
        if (g2sd_store(i) < 1.d-20) g2sd(i) = 0.d0
      end do

      call inid        ! Calculate derived parameters for "d" orbital work
!
!   Atomic number 102 is the capped bond.  It should have a very
!   small AM to prevent division by zero in REPPD, and to avoid
!   spurious effects due to normal AM values.
!
      am(102) = 1.D-10
!
!     DEBUG PRINTING.
!     THIS IS FORMATTED FOR DIRECT INSERTION INTO 'PARAM'
!
      if (index(keywrd,' DEP ') == 0) return
      call create_parameters_for_PMx_C
      return
    end subroutine calpar
!
!
!
!
  subroutine create_parameters_for_PMx_C
!
    USE molkst_C, only : line, keywrd
    USE parameters_C, only : guess1, guess2, guess3, gpp, gp2, hsp, gss, gsp, betas, betap, betad, &
    zs, zp, zd, uss, upp, udd, zsn, zpn, zdn, pocord, alpb, xfac, &
    f0sd, g2sd, main_group, alp, polvol, CPE_Zeta, CPE_Z0, CPE_B, CPE_Xlo, CPE_Xhi, v_par, t_par
    USE chanel_C, only : iw, param_out
    USE elemts_C, only : atom_names
!
    implicit none
!
!  Local
!
    integer :: i, j, l1, l2, l3
    logical :: lnew
    character :: num1*2, num2*2, name_1*10, name_2*10, name_3*10, file_name*120
    double precision :: abond, fff
    double precision, external :: reada
!
!   Define all text variables here
!
! Module names
!
    name_1 = "6_ORG"  ! PM6_ORG
    l1 = len_trim(name_1)
!
! Name of parameter set
!
    name_2 = "_org"
    l2 = len_trim(name_2)
!
! Name of parameter set for global parameters
!
    name_3 = "_org("
    l3 = len_trim(name_3)
!
!  End if definitions
!
    file_name = "parameters_for_PM"//name_1(:l1)//"_C.F90"
    i = index(keywrd, "ifiles_8")
    if (i > 0) then
      i = nint(reada(keywrd, i))
    else
      i = iw
    end if
    write (i, '(//10x, a, /)')" A new file called '"//trim(file_name)//"' will be written"
    call add_path(file_name)
    j = 0
97  open (unit=param_out, file=trim(file_name), status="UNKNOWN", iostat=i)
    if ( i /= 0 .and. j < 10)then
!
!  The file can exist, but is not currently accessible
!
      write(iw,*) "Problem writing '"//trim(file_name)//"'"
      call sleep(3)
      j = j + 1
      goto 97
    end if
    if (j > 9) return
    rewind (param_out)
    write (param_out,"(a)") "! Molecular Orbital PACkage (MOPAC)"
    write (param_out,"(a)") "! Copyright (C) 2021, Virginia Polytechnic Institute and State University"
    write (param_out,"(a)") "!"
    write (param_out,"(a)") "! MOPAC is free software: you can redistribute it and/or modify it under"
    write (param_out,"(a)") "! the terms of the GNU Lesser General Public License as published by"
    write (param_out,"(a)") "! the Free Software Foundation, either version 3 of the License, or"
    write (param_out,"(a)") "! (at your option) any later version."
    write (param_out,"(a)") "!"
    write (param_out,"(a)") "! MOPAC is distributed in the hope that it will be useful,"
    write (param_out,"(a)") "! but WITHOUT ANY WARRANTY; without even the implied warranty of"
    write (param_out,"(a)") "! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the"
    write (param_out,"(a)") "! GNU Lesser General Public License for more details."
    write (param_out,"(a)") "!"
    write (param_out,"(a)") "! You should have received a copy of the GNU Lesser General Public License"
    write (param_out,"(a)") "! along with this program.  If not, see <https://www.gnu.org/licenses/>."
    write (param_out,"(a)") " "
    write (param_out,"(a)") "  module Parameters_for_PM"//name_1(:l1)//"_C"
    write (param_out,"(a)") "    double precision, dimension(107) :: uss"//name_2(:l2)//", upp"//name_2(:l2)//", udd" &
      //name_2(:l2)//", zs"//name_2(:l2)//", zp"//name_2(:l2)//", zd"//name_2(:l2)//", betas"//name_2(:l2)//", &"
    write (param_out,"(a)") "    betap"//name_2(:l2)//", betad"//name_2(:l2)//", gss"//name_2(:l2)//", gsp" &
      //name_2(:l2)//", gpp"//name_2(:l2)//", gp2"//name_2(:l2)//", hsp"//name_2(:l2)//", polvo" &
      //name_2(:l2)//", poc_"//name_2(:l2)//", &"
    write (param_out,"(a)") "    zsn"//name_2(:l2)//", zpn"//name_2(:l2)//", zdn"//name_2(:l2)//", f0sd" &
      //name_2(:l2)//", g2sd"//name_2(:l2)//", alp"//name_2(:l2)//", &"
    write (param_out,"(a)") "    CPE_Zet"//name_2(:l2)//", CPE_Z0"//name_2(:l2)//", CPE_B"//name_2(:l2)//", CPE_Xlo" &
      //name_2(:l2)//", CPE_Xhi"//name_2(:l2)//""
    write (param_out,"(a)") "    double precision :: v_par"//name_2(:l2)//"(60) "
    write (param_out,"(a)") "    double precision, dimension(107,4) :: gues"//name_2(:l2)//"1, gues" &
      //name_2(:l2)//"2, gues"//name_2(:l2)//"3"
    do i = 1, 107
      if (zs(i) < 1.d-20 .and. gss(i) < 1.d-20 .and. Abs(guess1(i,1)) &
      + Abs(alp(i)) < 1.d-20) cycle
      write (param_out,"('!')")
      write (param_out, '("!",20X,"Data for Element ",I3,5x,a)') i, atom_names(i)
      write (param_out,"('!')")
      if (uss(i) /= 0.D0) write (param_out, &
        '(6X,"data     uss'//name_2(:l2)//'(",I3,")/",f17.8,"D0/")') i, uss(i)
      if (upp(i) /= 0.D0) write (param_out, &
        '(6X,"data     upp'//name_2(:l2)//'(",I3,")/",f17.8,"D0/")') i, upp(i)
      if (udd(i) /= 0.D0) write (param_out, &
        '(6X,"data     udd'//name_2(:l2)//'(",I3,")/",f17.8,"D0/")') i, udd(i)
      if (betas(i) /= 0.D0) write (param_out, &
        '(6X,"data   betas'//name_2(:l2)//'(",I3,")/",f17.8,"D0/")') i, betas(i)
      if (betap(i) /= 0.D0) write (param_out, &
        '(6X,"data   betap'//name_2(:l2)//'(",I3,")/",f17.8,"D0/")') i, betap(i)
      if (betad(i) /= 0.D0) write (param_out, &
        '(6X,"data   betad'//name_2(:l2)//'(",I3,")/",f17.8,"D0/")') i, betad(i)
      if (zs(i) /= 0.D0) write (param_out, &
        '(6X,"data      zs'//name_2(:l2)//'(",I3,")/",f17.8,"D0/")') i, zs(i)
      if (zp(i) /= 0.D0) write (param_out, &
        '(6X,"data      zp'//name_2(:l2)//'(",I3,")/",f17.8,"D0/")') i, zp(i)
      if (zd(i) /= 0.D0) write (param_out, &
        '(6X,"data      zd'//name_2(:l2)//'(",I3,")/",f17.8,"D0/")') i, zd(i)
      if (zsn(i) /= 0.D0) write (param_out, &
        '(6X,"data     zsn'//name_2(:l2)//'(",I3,")/",f17.8,"D0/")') i, zsn(i)
      if (zpn(i) /= 0.D0) write (param_out, &
        '(6X,"data     zpn'//name_2(:l2)//'(",I3,")/",f17.8,"D0/")') i, zpn(i)
      if (zdn(i) /= 0.D0) write (param_out, &
        '(6X,"data     zdn'//name_2(:l2)//'(",I3,")/",f17.8,"D0/")') i, zdn(i)
      if (alp(i) /= 0.D0) write (param_out, &
        '(6X,"data     alp'//name_2(:l2)//'(",I3,")/",f17.8,"D0/")') i, alp(i)
      if (gss(i) /= 0.D0) write (param_out, &
        '(6X,"data     gss'//name_2(:l2)//'(",I3,")/",f17.8,"D0/")') i, gss(i)
      if (gsp(i) /= 0.D0) write (param_out, &
        '(6X,"data     gsp'//name_2(:l2)//'(",I3,")/",f17.8,"D0/")') i, gsp(i)
      if (gpp(i) /= 0.D0) write (param_out, &
        '(6X,"data     gpp'//name_2(:l2)//'(",I3,")/",f17.8,"D0/")') i, gpp(i)
      if (gp2(i) /= 0.D0) write (param_out, &
        '(6X,"data     gp2'//name_2(:l2)//'(",I3,")/",f17.8,"D0/")') i, gp2(i)
      if (hsp(i) /= 0.D0) write (param_out, &
        '(6X,"data     hsp'//name_2(:l2)//'(",I3,")/",f17.8,"D0/")') i, hsp(i)
      if (pocord(i) /= 0.D0) write (param_out, &
        '(6X,"data    poc_'//name_2(:l2)//'(",I3,")/",f17.8,"D0/")') i, pocord(i)
      if (polvol(i) /= 0.D0) write (param_out, &
        '(6X,"data   polvo'//name_2(:l2)//'(",I3,")/",f17.8,"D0/")') i, polvol(i)
      if (CPE_Zeta(i) /= 0.D0) write (param_out, &
        '(6X,"data CPE_Zet'//name_2(:l2)//'(",I3,")/",f17.8,"D0/")') i, CPE_Zeta(i)
      if (CPE_Z0(i) /= 0.D0) write (param_out, &
        '(6X,"data  CPE_Z0'//name_2(:l2)//'(",I3,")/",f17.8,"D0/")') i, CPE_Z0(i)
      if (CPE_B(i) /= 0.D0) write (param_out, &
        '(6X,"data   CPE_B'//name_2(:l2)//'(",I3,")/",f17.8,"D0/")') i, CPE_B(i)
      if (CPE_Xlo(i) /= 0.D0) write (param_out, &
        '(6X,"data CPE_Xlo'//name_2(:l2)//'(",I3,")/",f17.8,"D0/")') i, CPE_Xlo(i)
      if (CPE_Xhi(i) /= 0.D0) write (param_out, &
        '(6X,"data CPE_Xhi'//name_2(:l2)//'(",I3,")/",f17.8,"D0/")') i, CPE_Xhi(i)
      if (.not. main_group(i) .and. f0sd(i) /= 0.D0) write (param_out, &
        '(6X,"data    f0sd'//name_2(:l2)//'(",I3,")/",f17.8,"D0/")') i, f0sd(i)
      if (.not. main_group(i) .and. g2sd(i) /= 0.D0) write (param_out, &
        '(6X,"data    g2sd'//name_2(:l2)//'(",I3,")/",f17.8,"D0/")') i, g2sd(i)
      do j = 1, 4
        if (guess1(i,j) /= 0.D0) &
        write (param_out, &
        '(6X,"data gues'//name_2(:l2)//'1(",I3,",",I1,")/",          f17.8,"D0/")')&
          i, j, guess1(i,j)
        if (guess2(i,j) /= 0.D0) &
        write (param_out, &
        '(6X,"data gues'//name_2(:l2)//'2(",I3,",",I1,")/",          f17.8,"D0/")')&
          i, j, guess2(i,j)
        if (guess3(i,j) /= 0.D0) &
        write (param_out, &
        '(6X,"data gues'//name_2(:l2)//'3(",I3,",",I1,")/",          f17.8,"D0/")')&
          i, j, guess3(i,j)
      end do
    end do
!
! Write out the global parameters
!
    write (param_out,'(2("!",/),"!",21x,a,2(/,"!"))')"Global parameters"
    do i = 1, 60
      if (abs(v_par(i)) > 1.d-10) then
          num1 = "2"
          num2 = "1"
        line = " "
        if (t_par(i) /= " ") line = "  ! "//trim(t_par(i))
        if (i > 9) then
          num1 = "1"
          num2 = "2"
        end if
        write(param_out,'(6x, "data  ", a, i'//num2//', a2, f1'//num1//'.8,"d0/", a)') &
          "v_par"//name_3(:l3), i, ")/", v_par(i), trim(line)
      end if
    end do
    write (param_out,"(a)") '  contains'

    write (param_out,"(a)") &
    & '  subroutine alpb_and_xfac_pm'//name_1(:l1)//'', &
    & '    use parameters_C, only : xfac, alpb'
!
!  Write out all the diatomic parameters
!
    do i = 1, 100
      lnew = .true.
      do j = 1, i
        abond = alpb(i,j)
        if (abond > 1.d-4) then
          fff = xfac(i,j)
          if ( lnew ) then
          lnew = .false.
          write (param_out,*)"!"
          end if
            write(param_out,"(6x,a,i2,a,i2,a,f14.8,a)") &
            "alpb(",i,",",j,") = ", abond,"d0 !"//atom_names(i)//" - "//atom_names(j)
            write(param_out,"(6x,a,i2,a,i2,a,f14.8,a)") &
            "xfac(",i,",",j,") = ", fff,"d0 !"//atom_names(i)//" - "//atom_names(j)
        end if
      end do
    end do
    write (param_out,"(a)") '    end subroutine alpb_and_xfac_pm'//name_1(:l1)
    write (param_out,"(a)") "  end module Parameters_for_PM"//name_1(:l1)//"_C"
    close (param_out)
    return
  end subroutine create_parameters_for_PMx_C
