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

      subroutine inid
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE parameters_C, only : am, ad, aq, dd, qq, dorbs, po, ddp, pocord, &
      zdn, natorb, main_group
      use mndod_C, only : repd, aij
!     *
!     DEFINE SEVERAL PARAMETERS FOR D-ORBITAL CALCULATIONS.
!     *
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ni, i
!-----------------------------------------------
      repd = 0.d0
      aij = 0.d0
      do ni = 1, 107
        if (.not.dorbs(ni)) cycle
        call aijm (ni)
        if (zdn(ni) > 1.d-4) call inighd (ni)
        call ddpo (ni)
      end do
      do i = 1, 106
        if( natorb(i) < 6 .or. main_group(i)) then
          if (am(i) < 1.D-4) am(i) = 1.D0
          po(1,i) = 0.5D0/am(i)
          if (ad(i) > 1.D-5) po(2,i) = 0.5D0/ad(i)
          if (aq(i) > 1.D-5) po(3,i) = 0.5D0/aq(i)
          po(7,i) = po(1,i)
          ddp(2,i) = dd(i)
          ddp(3,i) = qq(i)*sqrt(2.D0)
        end if
        po(9,i) = po(1,i)
        if (pocord(i) > 1.D-5) po(9,i) = pocord(i)
      end do
      po(2,1) = 0.D0
      po(3,1) = 0.D0
      return
      end subroutine inid


      subroutine inighd(ni)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use parameters_C, only : dorbs, f0dd, f2dd, f4dd, f0sd, g2sd, f0pd, f2pd, g1pd, &
      & g3pd
      use mndod_C, only : repd
!     *
!     ONE-CENTER TWO-ELECTRON INTEGRALS FOR SPD-BASIS.
!     *
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: ni
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      double precision :: s3, s5, s15, r066, r266, r466, r016, r244, r036, r236, &
        r155, r355, r125, r234, r246
!-----------------------------------------------
      data s3/ 1.7320508D0/
      data s5/ 2.23606797D0/
      data s15/ 3.87298334D0/
!
      if (.not.dorbs(ni)) return
!
!     SLATER-CONDON PARAMETERS (RLIJ).
!     FIRST  DIGIT (L)  L QUANTUM NUMBER OF SLATER-CONDON PARAMETER.
!     SECOND DIGIT (I)  SS 1, SP 2, PP 3, SD 4, PD 5, DD 6 - ELECTRON 1.
!     SECOND DIGIT (J)  SS 1, SP 2, PP 3, SD 4, PD 5, DD 6 - ELECTRON 2.
      call scprm (ni, r066, r266, r466, r016, r244, r036, r236, r155, r355, &
        r125, r234, r246)
      if (f0sd(ni) > 0.001D0) r016 = f0sd(ni)
      if (g2sd(ni) > 0.001D0) r244 = g2sd(ni)
      call eiscor (r016, r066, r244, r266, r466, ni)
      repd(1,ni) = r016
      repd(2,ni) = 2.d0/(3.D0*s5)*r125
      repd(3,ni) = 1.d0/s15*r125
      repd(4,ni) = 2.d0/(5.D0*s5)*r234
      repd(5,ni) = r036 + 4.D0/35.D0*r236
      repd(6,ni) = r036 + 2.D0/35.D0*r236
      repd(7,ni) = r036 - 4.D0/35.D0*r236
      repd(8,ni) = -1.d0/(3.D0*s5)*r125
      repd(9,ni) = sqrt(3.D0/125.D0)*r234
      repd(10,ni) = s3/35.D0*r236
      repd(11,ni) = 3.D0/35.D0*r236
      repd(12,ni) = -1.d0/(5.D0*s5)*r234
      repd(13,ni) = r036 - 2.D0/35.D0*r236
      repd(14,ni) = -2.D0*s3/35.D0*r236
      repd(15,ni) = -repd(3,ni)
      repd(16,ni) = -repd(11,ni)
      repd(17,ni) = -repd(9,ni)
      repd(18,ni) = -repd(14,ni)
      repd(19,ni) = 1.d0/5.D0*r244
      repd(20,ni) = 2.D0/(7.D0*s5)*r246
      repd(21,ni) = repd(20,ni)/2.d0
      repd(22,ni) = -repd(20,ni)
      repd(23,ni) = 4.D0/15.D0*r155 + 27.D0/245.D0*r355
      repd(24,ni) = 2.D0*s3/15.D0*r155 - 9.D0*s3/245.D0*r355
      repd(25,ni) = 1.d0/15.D0*r155 + 18.D0/245.D0*r355
      repd(26,ni) = (-s3/15.D0*r155) + 12.D0*s3/245.D0*r355
      repd(27,ni) = (-s3/15.D0*r155) - 3.D0*s3/245.D0*r355
      repd(28,ni) = -repd(27,ni)
      repd(29,ni) = r066 + 4.D0/49.D0*r266 + 4.D0/49.D0*r466
      repd(30,ni) = r066 + 2.D0/49.D0*r266 - 24.D0/441.D0*r466
      repd(31,ni) = r066 - 4.D0/49.D0*r266 + 6.D0/441.D0*r466
      repd(32,ni) = sqrt(3.D0/245.D0)*r246
      repd(33,ni) = 1.d0/5.D0*r155 + 24.D0/245.D0*r355
      repd(34,ni) = 1.d0/5.D0*r155 - 6.D0/245.D0*r355
      repd(35,ni) = 3.D0/49.D0*r355
      repd(36,ni) = 1.d0/49.D0*r266 + 30.D0/441.D0*r466
      repd(37,ni) = s3/49.D0*r266 - 5.D0*s3/441.D0*r466
      repd(38,ni) = r066 - 2.D0/49.D0*r266 - 4.D0/441.D0*r466
      repd(39,ni) = (-2.D0*s3/49.D0*r266) + 10.D0*s3/441.D0*r466
      repd(40,ni) = -repd(32,ni)
      repd(41,ni) = -repd(34,ni)
      repd(42,ni) = -repd(35,ni)
      repd(43,ni) = -repd(37,ni)
      repd(44,ni) = 3.D0/49.D0*r266 + 20.D0/441.D0*r466
      repd(45,ni) = -repd(39,ni)
      repd(46,ni) = 1.D0/5.D0*r155 - 3.D0/35.D0*r355
      repd(47,ni) = -repd(46,ni)
      repd(48,ni) = 4.D0/49.D0*r266 + 15.D0/441.D0*r466
      repd(49,ni) = 3.D0/49.D0*r266 - 5.D0/147.D0*r466
      repd(50,ni) = -repd(49,ni)
      repd(51,ni) = r066 + 4.D0/49.D0*r266 - 34.D0/441.D0*r466
      repd(52,ni) = 35.D0/441.D0*r466
      f0dd(ni) = r066
      f2dd(ni) = r266
      f4dd(ni) = r466
      f0sd(ni) = r016
      g2sd(ni) = r244
      f0pd(ni) = r036
      f2pd(ni) = r236
      g1pd(ni) = r155
      g3pd(ni) = r355
      return
      end subroutine inighd

      double precision function poij (l, d, fg)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use funcon_C, only : ev
!     *
!     DETERMINE ADDITIVE TERMS RHO=POIJ FOR TWO-CENTER TWO-ELECTRON
!     INTEGRALS FROM THE REQUIREMENT THAT THE APPROPRIATE ONE-CENTER
!     TWO-ELECTRON INTEGRALS ARE REPRODUCED.
!     *
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: l
      double precision , intent(in) :: d
      double precision , intent(in) :: fg
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      double precision, parameter :: epsil = 1.0D-08
      double precision, parameter :: g1 = 0.382D0
      double precision, parameter :: g2 = 0.618D0
      integer, parameter :: niter = 100
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      double precision :: dsq, ev4, ev8, a1, a2, delta, y1, y2, f1, f2
!-----------------------------------------------
      f1 = 0.d0
      f2 = 0.d0
      if (l == 0) then
        poij = 0.5d0*ev/fg
        return
      end if
! *** HIGHER TERMS.
      dsq = d*d
      ev4 = ev*0.25d0
      ev8 = ev/8.0D0
      a1 = 0.1D0
      a2 = 5.0D0
      if (l == 1) then
        do i = 1, niter
          delta = a2 - a1
          if (delta < epsil) exit
          y1 = a1 + delta*g1
          y2 = a1 + delta*g2
          f1 = (ev4*(1.d0/y1 - 1.d0/sqrt(y1**2 + dsq)) - fg)**2
          f2 = (ev4*(1.d0/y2 - 1.d0/sqrt(y2**2 + dsq)) - fg)**2
          if (f1 < f2) then
            a2 = y2
          else
            a1 = y1
          end if
        end do
      else
        if (l == 2) then
          do i = 1, niter
            delta = a2 - a1
            if (delta < epsil) exit
            y1 = a1 + delta*g1
            y2 = a1 + delta*g2
            f1 = (ev8*(1.d0/y1 - 2.d0/sqrt(y1**2 + dsq*0.5d0) + 1.d0/sqrt(y1**2 + &
              dsq)) - fg)**2
            f2 = (ev8*(1.d0/y2 - 2.d0/sqrt(y2**2 + dsq*0.5d0) + 1.d0/sqrt(y2**2 + &
              dsq)) - fg)**2
            if (f1 < f2) then
              a2 = y2
            else
              a1 = y1
            end if
          end do
        end if
      end if
!     DEFINE ADDITIVE TERM AFTER CONVERGENCE OF ITERATIONS.
      if (f1 >= f2) then
        poij = a2
      else
        poij = a1
      end if
      return
      end function poij


      subroutine printp(i, para, value, txt)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE chanel_C, only : iw
#if MOPAC_F2003
      USE, INTRINSIC :: IEEE_ARITHMETIC
#endif
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: i
      double precision , intent(in) :: value
      character , intent(in) :: para*(*)
      character , intent(in) :: txt*(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      double precision :: loc_value
#ifdef MOPAC_F2003
      if (ieee_is_nan(value)) then
#else
      if (isnan(value)) then
#endif
        loc_value = 0.d0
      else
        loc_value = value
      end if
      if (abs(loc_value) > 1.D-5) write (iw, '(I4,A7,2X,F13.8,2X,A)') i, para, &
        loc_value, txt
      return
      end subroutine printp

      subroutine prtpar
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use common_arrays_C, only : nat
      use molkst_C, only : numat
      use chanel_C, only : iw
      USE parameters_C, only : alp, tore, uss, upp, udd, zs, zp, zd, zsn, &
      zpn, zdn, gpp, gp2, hsp, gss, gsp, eisol, eheat, betas, betap, betad, &
      po, ddp, f0dd, f2dd, f4dd, f0sd, g2sd, f0pd, f2pd, alpb, xfac, &
      & g1pd, g3pd, guess1, guess2, guess3
#if MOPAC_F2003
      USE, INTRINSIC :: IEEE_ARITHMETIC
#endif
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      logical , dimension(107) :: used
      integer ::  i, j
!-----------------------------------------------
      used = .false.
      used(nat(:numat)) = .true.
      write (iw, *)
      write (iw, *) 'PARAMETER VALUES USED IN THE CALCULATION'
      write (iw, *)
      write (iw, *) ' NI    TYPE        VALUE     UNIT'
      write (iw, *)
      do i = 1, 100
        if ( .not. used(i)) cycle
        write (iw, *)
        call printp (i, 'USS  ', uss(i), 'EV        ONE-CENTER ENERGY FOR S')
        call printp (i, 'UPP  ', upp(i), 'EV        ONE-CENTER ENERGY FOR P')
        call printp (i, 'UDD  ', udd(i), 'EV        ONE-CENTER ENERGY FOR D')
        call printp (i, 'ZS   ', zs(i), 'AU        ORBITAL EXPONENT  FOR S')
        call printp (i, 'ZP   ', zp(i), 'AU        ORBITAL EXPONENT  FOR P')
        call printp (i, 'ZD   ', zd(i), 'AU        ORBITAL EXPONENT  FOR D')
        call printp (i, 'BETAS', betas(i), 'EV        BETA PARAMETER    FOR S')
        call printp (i, 'BETAP', betap(i), 'EV        BETA PARAMETER    FOR P')
        call printp (i, 'BETAD', betad(i), 'EV        BETA PARAMETER    FOR D')
        call printp (i, 'ALP  ', alp(i), '(1/A)     ALPHA PARAMETER   FOR CORE'&
          )
        call printp (i, 'GSS  ', gss(i), &
          'EV        ONE-CENTER INTEGRAL (SS,SS)')
        call printp (i, 'GPP  ', gpp(i), &
          'EV        ONE-CENTER INTEGRAL (PP,PP)')
        call printp (i, 'GSP  ', gsp(i), &
          'EV        ONE-CENTER INTEGRAL (SS,PP)')
        call printp (i, 'GP2  ', gp2(i), &
          'EV        ONE-CENTER INTEGRAL (PP*,PP*)')
        call printp (i, 'HSP  ', hsp(i), &
          'EV        ONE-CENTER INTEGRAL (SP,SP)')
        call printp (i, 'ZSN  ', zsn(i), &
          'AU        INTERNAL EXPONENT FOR S - (IJ,KL)')
        call printp (i, 'ZPN  ', zpn(i), &
          'AU        INTERNAL EXPONENT FOR P - (IJ,KL)')
        call printp (i, 'ZDN  ', zdn(i), &
          'AU        INTERNAL EXPONENT FOR D - (IJ,KL)')
        call printp (i, 'F0DD ', f0dd(i), &
          'EV        SLATER-CONDON PARAMETER F0DD')
        call printp (i, 'F2DD ', f2dd(i), &
          'EV        SLATER-CONDON PARAMETER F2DD')
        call printp (i, 'F4DD ', f4dd(i), &
          'EV        SLATER-CONDON PARAMETER F4DD')
        call printp (i, 'F0SD ', f0sd(i), &
          'EV        SLATER-CONDON PARAMETER F0SD')
        call printp (i, 'G2SD ', g2sd(i), &
          'EV        SLATER-CONDON PARAMETER G2SD')
        call printp (i, 'F0PD ', f0pd(i), &
          'EV        SLATER-CONDON PARAMETER F0PD')
        call printp (i, 'F2PD ', f2pd(i), &
          'EV        SLATER-CONDON PARAMETER F2PD')
        call printp (i, 'G1PD ', g1pd(i), &
          'EV        SLATER-CONDON PARAMETER G1PD')
        call printp (i, 'G3PD ', g3pd(i), &
          'EV        SLATER-CONDON PARAMETER G3PD')
        call printp (i, 'DD2  ', ddp(2,i), &
          'BOHR      CHARGE SEPARATION, SP, L=1')
        call printp (i, 'DD3  ', ddp(3,i), &
          'BOHR      CHARGE SEPARATION, PP, L=2')
        call printp (i, ' =   ', ddp(3,i)/sqrt(2.d0), &
          'BOHR      USING ORIGINAL MNDO PAPER FORMULA')
        call printp (i, 'DD4  ', ddp(4,i), &
          'BOHR      CHARGE SEPARATION, SD, L=2')
        call printp (i, 'DD5  ', ddp(5,i), &
          'BOHR      CHARGE SEPARATION, PD, L=1')
        call printp (i, 'DD6  ', ddp(6,i), &
          'BOHR      CHARGE SEPARATION, DD, L=2')
        call printp (i, 'PO1  ', po(1,i), &
          'BOHR      KLOPMAN-OHNO TERM, SS, L=0')
        call printp (i, 'PO2  ', po(2,i), &
          'BOHR      KLOPMAN-OHNO TERM, SP, L=1')
        call printp (i, 'PO3  ', po(3,i), &
          'BOHR      KLOPMAN-OHNO TERM, PP, L=2')
        call printp (i, 'PO4  ', po(4,i), &
          'BOHR      KLOPMAN-OHNO TERM, SD, L=2')
        call printp (i, 'PO5  ', po(5,i), &
          'BOHR      KLOPMAN-OHNO TERM, PD, L=1')
        call printp (i, 'PO6  ', po(6,i), &
          'BOHR      KLOPMAN-OHNO TERM, DD, L=2')
        call printp (i, 'PO7  ', po(7,i), &
          'BOHR      KLOPMAN-OHNO TERM, PP, L=0')
        call printp (i, 'PO8  ', po(8,i), &
          'BOHR      KLOPMAN-OHNO TERM, DD, L=0')
        call printp (i, 'PO9  ', po(9,i), 'BOHR      KLOPMAN-OHNO TERM, CORE')
        call printp (i, 'CORE ', tore(i), 'E         CORE CHARGE')
        call printp (i, 'EHEAT', eheat(i), &
          'KCAL/MOL  HEAT OF FORMATION OF THE ATOM (EXP)')
        call printp (i, 'EISOL', eisol(i), &
          'EV        TOTAL ENERGY OF THE ATOM (CALC)')
        call printp(i, 'FN11 ', guess1(i,1), 'CORE-CORE VDW MULTIPLIER 1')
        call printp(i, 'FN21 ', guess2(i,1), 'CORE-CORE VDW EXPONENT 1')
        call printp(i, 'FN31 ', guess3(i,1), 'CORE-CORE VDW POSITION 1')
        call printp(i, 'FN12 ', guess1(i,2), 'CORE-CORE VDW MULTIPLIER 2')
        call printp(i, 'FN22 ', guess2(i,2), 'CORE-CORE VDW EXPONENT 2')
        call printp(i, 'FN32 ', guess3(i,2), 'CORE-CORE VDW POSITION 2')
        call printp(i, 'FN13 ', guess1(i,3), 'CORE-CORE VDW MULTIPLIER 3')
        call printp(i, 'FN23 ', guess2(i,3), 'CORE-CORE VDW EXPONENT 3')
        call printp(i, 'FN33 ', guess3(i,3), 'CORE-CORE VDW POSITION 3')
        call printp(i, 'FN14 ', guess1(i,4), 'CORE-CORE VDW MULTIPLIER 4')
        call printp(i, 'FN24 ', guess2(i,4), 'CORE-CORE VDW EXPONENT 4')
        call printp(i, 'FN34 ', guess3(i,4), 'CORE-CORE VDW POSITION 4')
        do j = 1, 100
#ifdef MOPAC_F2003
          if (ieee_is_nan(alpb(i,j))) alpb(i,j)= 0.d0
#else
          if (isnan(alpb(i,j))) alpb(i,j)= 0.d0
#endif
        if (Abs (alpb(i,j)) > 1.d-5 .and. used(j)) then
            write (iw, "(I4,A6,i2,F13.8,2X,A)") i, "ALPB_", j,alpb(i,j), "ALPB factor"
            write (iw, "(I4,A6,i2,F13.8,2X,A)") i, "XFAC_", j,xfac(i,j), "XFAC factor"
          end if
       end do
         end do
      return
      end subroutine prtpar


      subroutine reppd(ni, nj, rij, ri, gab)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE parameters_C, only : natorb, dd, qq, am, ad, aq, po
      use molkst_C, only : l_feather
      use funcon_C, only : ev, a0
    !  USE chanel_C, only : iw
!***********************************************************************
!
!..VECTOR VERSION WRITTEN BY ERNEST R. DAVIDSON, INDIANA UNIVERSITY
!
!
!  REPPD CALCULATES THE TWO-ELECTRON REPULSION INTEGRALS AND THE
!       NUCLEAR ATTRACTION INTEGRALS.
!
!     ON INPUT RIJ     = INTERATOMIC DISTANCE
!              NI      = ATOM NUMBER OF FIRST ATOM
!              NJ      = ATOM NUMBER OF SECOND ATOM
!
!
!    ON OUTPUT RI      = ARRAY OF TWO-ELECTRON REPULSION INTEGRALS
!
!
! *** THIS ROUTINE COMPUTES THE TWO-CENTRE REPULSION INTEGRALS AND THE
! *** NUCLEAR ATTRACTION INTEGRALS.
! *** THE TWO-CENTER REPULSION INTEGRALS (OVER LOCAL COORDINATES) ARE
! *** STORED AS FOLLOWS (WHERE P-SIGMA = O,  AND P-PI = P AND P* )
!     (SS/SS)=1,   (SO/SS)=2,   (OO/SS)=3,   (PP/SS)=4,   (SS/OS)=5,
!     (SO/SO)=6,   (SP/SP)=7,   (OO/SO)=8,   (PP/SO)=9,   (PO/SP)=10,
!     (SS/OO)=11,  (SS/PP)=12,  (SO/OO)=13,  (SO/PP)=14,  (SP/OP)=15,
!     (OO/OO)=16,  (PP/OO)=17,  (OO/PP)=18,  (PP/PP)=19,  (PO/PO)=20,
!     (PP/P*P*)=21,   (P*P/P*P)=22.
! *** NI AND NJ ARE THE ATOMIC NUMBERS OF THE TWO ELEMENTS.
!
!***********************************************************************
!***********************************************************************
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: ni
      integer , intent(in) :: nj
      double precision , intent(in) :: rij
      double precision , intent(out) :: ri(22), gab
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(22) :: nri
      integer :: i
      double precision, dimension(72) :: arg, sqr
      double precision :: td, pp, ev1, ev2, ev3, ev4, r, aee, da, qa, ade, aqe, &
        rsq, xxx, ee, db, qb, aed, aeq, axx, adq, aqd, aqq, yyy, zzz, www, &
        dze, qzze, qxxe, edz, eqzz, eqxx, dxdx, dzdz, dzqxx, qxxdz, dzqzz, &
        qzzdz, qxxqxx, qxxqyy, qxxqzz, qzzqxx, qzzqzz, dxqxz, qxzdx, qxzqxz, const, &
        point
      logical :: si, sj

      save td, pp
!-----------------------------------------------

      data td/ 2.D00/
      data pp/ 0.5D00/
      data nri/ 1, -1, 1, 1, -1, 1, 1, -1, -1, -1, 1, 1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1/
        ev1 = ev/2
        ev2 = ev1/2
        ev3 = ev2/2
        ev4 = ev3/2
        ri = 0.D0
        r = rij/a0
        si = natorb(ni) >= 3
        sj = natorb(nj) >= 3
!
!
!  Calculate the <ss|ss> integral and G_(AB) here - so that any changes to the
!  two-electron two-center integrals can also be made to the core-core term as well.
!
        aee = po(9,ni)+po(9,nj)
        aee = aee*aee
        gab = ev/sqrt(r*r + aee)  !  Used by the core-core term only
!
        aee = pp/am(ni) + pp/am(nj) !  Used by the two-electron, two-center term only
        aee = aee*aee
!
        if (.not.si .and. .not.sj) then
!
!     HYDROGEN - HYDROGEN  (SS/SS)
!
          ri(1) = ev/sqrt(r*r + aee)
!
        else if (si .and. .not.sj) then
!
!     HEAVY ATOM - HYDROGEN
!

          da = dd(ni)
          qa = qq(ni)*td
          ade = pp/ad(ni) + pp/am(nj)
          ade = ade*ade
          aqe = pp/aq(ni) + pp/am(nj)
          aqe = aqe*aqe
          rsq = r*r
          arg(1) = rsq + aee
          xxx = r + da
          arg(2) = xxx*xxx + ade
          xxx = r - da
          arg(3) = xxx*xxx + ade
          xxx = r + qa
          arg(4) = xxx*xxx + aqe
          xxx = r - qa
          arg(5) = xxx*xxx + aqe
          arg(6) = rsq + aqe
          arg(7) = arg(6) + qa*qa
          do i = 1, 7
            sqr(i) = sqrt(arg(i))
          end do
          ee = ev/sqr(1)
          ri(1) = ee
          ri(2) = ev1/sqr(2) - ev1/sqr(3)
          ri(3) = ee + ev2/sqr(4) + ev2/sqr(5) - ev1/sqr(6)
          ri(4) = ee + ev1/sqr(7) - ev1/sqr(6)
!
        else if (.not.si .and. sj) then
!
!     HYDROGEN - HEAVY ATOM
!
          db = dd(nj)
          qb = qq(nj)*td
          aed = pp/am(ni) + pp/ad(nj)
          aed = aed*aed
          aeq = pp/am(ni) + pp/aq(nj)
          aeq = aeq*aeq
          rsq = r*r
          arg(1) = rsq + aee
          xxx = r - db
          arg(2) = xxx*xxx + aed
          xxx = r + db
          arg(3) = xxx*xxx + aed
          xxx = r - qb
          arg(4) = xxx*xxx + aeq
          xxx = r + qb
          arg(5) = xxx*xxx + aeq
          arg(6) = rsq + aeq
          arg(7) = arg(6) + qb*qb
          do i = 1, 7
            sqr(i) = sqrt(arg(i))
          end do
          ee = ev/sqr(1)
          ri(1) = ee
          ri(5) = ev1/sqr(2) - ev1/sqr(3)
          ri(11) = ee + ev2/sqr(4) + ev2/sqr(5) - ev1/sqr(6)
          ri(12) = ee + ev1/sqr(7) - ev1/sqr(6)
!
        else
!
!     HEAVY ATOM - HEAVY ATOM
!
!     DEFINE CHARGE SEPARATIONS.
          da = dd(ni)
          db = dd(nj)
          qa = qq(ni)*td
          qb = qq(nj)*td
!
          ade = pp/ad(ni) + pp/am(nj)
          ade = ade*ade
          aqe = pp/aq(ni) + pp/am(nj)
          aqe = aqe*aqe
          aed = pp/am(ni) + pp/ad(nj)
          aed = aed*aed
          aeq = pp/am(ni) + pp/aq(nj)
          aeq = aeq*aeq
          axx = pp/ad(ni) + pp/ad(nj)
          axx = axx*axx
          adq = pp/ad(ni) + pp/aq(nj)
          adq = adq*adq
          aqd = pp/aq(ni) + pp/ad(nj)
          aqd = aqd*aqd
          aqq = pp/aq(ni) + pp/aq(nj)
          aqq = aqq*aqq
          rsq = r*r
          arg(1) = rsq + aee
          xxx = r + da
          arg(2) = xxx*xxx + ade
          xxx = r - da
          arg(3) = xxx*xxx + ade
          xxx = r - qa
          arg(4) = xxx*xxx + aqe
          xxx = r + qa
          arg(5) = xxx*xxx + aqe
          arg(6) = rsq + aqe
          arg(7) = arg(6) + qa*qa
          xxx = r - db
          arg(8) = xxx*xxx + aed
          xxx = r + db
          arg(9) = xxx*xxx + aed
          xxx = r - qb
          arg(10) = xxx*xxx + aeq
          xxx = r + qb
          arg(11) = xxx*xxx + aeq
          arg(12) = rsq + aeq
          arg(13) = arg(12) + qb*qb
          xxx = da - db
          arg(14) = rsq + axx + xxx*xxx
          xxx = da + db
          arg(15) = rsq + axx + xxx*xxx
          xxx = r + da - db
          arg(16) = xxx*xxx + axx
          xxx = r - da + db
          arg(17) = xxx*xxx + axx
          xxx = r - da - db
          arg(18) = xxx*xxx + axx
          xxx = r + da + db
          arg(19) = xxx*xxx + axx
          xxx = r + da
          arg(20) = xxx*xxx + adq
          arg(21) = arg(20) + qb*qb
          xxx = r - da
          arg(22) = xxx*xxx + adq
          arg(23) = arg(22) + qb*qb
          xxx = r - db
          arg(24) = xxx*xxx + aqd
          arg(25) = arg(24) + qa*qa
          xxx = r + db
          arg(26) = xxx*xxx + aqd
          arg(27) = arg(26) + qa*qa
          xxx = r + da - qb
          arg(28) = xxx*xxx + adq
          xxx = r - da - qb
          arg(29) = xxx*xxx + adq
          xxx = r + da + qb
          arg(30) = xxx*xxx + adq
          xxx = r - da + qb
          arg(31) = xxx*xxx + adq
          xxx = r + qa - db
          arg(32) = xxx*xxx + aqd
          xxx = r + qa + db
          arg(33) = xxx*xxx + aqd
          xxx = r - qa - db
          arg(34) = xxx*xxx + aqd
          xxx = r - qa + db
          arg(35) = xxx*xxx + aqd
          arg(36) = rsq + aqq
          xxx = qa - qb
          arg(37) = arg(36) + xxx*xxx
          xxx = qa + qb
          arg(38) = arg(36) + xxx*xxx
          arg(39) = arg(36) + qa*qa
          arg(40) = arg(36) + qb*qb
          arg(41) = arg(39) + qb*qb
          xxx = r - qb
          arg(42) = xxx*xxx + aqq
          arg(43) = arg(42) + qa*qa
          xxx = r + qb
          arg(44) = xxx*xxx + aqq
          arg(45) = arg(44) + qa*qa
          xxx = r + qa
          arg(46) = xxx*xxx + aqq
          arg(47) = arg(46) + qb*qb
          xxx = r - qa
          arg(48) = xxx*xxx + aqq
          arg(49) = arg(48) + qb*qb
          xxx = r + qa - qb
          arg(50) = xxx*xxx + aqq
          xxx = r + qa + qb
          arg(51) = xxx*xxx + aqq
          xxx = r - qa - qb
          arg(52) = xxx*xxx + aqq
          xxx = r - qa + qb
          arg(53) = xxx*xxx + aqq
          qa = qq(ni)
          qb = qq(nj)
          xxx = da - qb
          xxx = xxx*xxx
          yyy = r - qb
          yyy = yyy*yyy
          zzz = da + qb
          zzz = zzz*zzz
          www = r + qb
          www = www*www
          arg(54) = xxx + yyy + adq
          arg(55) = xxx + www + adq
          arg(56) = zzz + yyy + adq
          arg(57) = zzz + www + adq
          xxx = qa - db
          xxx = xxx*xxx
          yyy = qa + db
          yyy = yyy*yyy
          zzz = r + qa
          zzz = zzz*zzz
          www = r - qa
          www = www*www
          arg(58) = zzz + xxx + aqd
          arg(59) = www + xxx + aqd
          arg(60) = zzz + yyy + aqd
          arg(61) = www + yyy + aqd
          xxx = qa - qb
          xxx = xxx*xxx
          arg(62) = arg(36) + td*xxx
          yyy = qa + qb
          yyy = yyy*yyy
          arg(63) = arg(36) + td*yyy
          arg(64) = arg(36) + td*(qa*qa + qb*qb)
          zzz = r + qa - qb
          zzz = zzz*zzz
          arg(65) = zzz + xxx + aqq
          arg(66) = zzz + yyy + aqq
          zzz = r + qa + qb
          zzz = zzz*zzz
          arg(67) = zzz + xxx + aqq
          arg(68) = zzz + yyy + aqq
          zzz = r - qa - qb
          zzz = zzz*zzz
          arg(69) = zzz + xxx + aqq
          arg(70) = zzz + yyy + aqq
          zzz = r - qa + qb
          zzz = zzz*zzz
          arg(71) = zzz + xxx + aqq
          arg(72) = zzz + yyy + aqq
          do i = 1, 72
            sqr(i) = sqrt(arg(i))
          end do
          ee = ev/sqr(1)
          dze = (-ev1/sqr(2)) + ev1/sqr(3)
          qzze = ev2/sqr(4) + ev2/sqr(5) - ev1/sqr(6)
          qxxe = ev1/sqr(7) - ev1/sqr(6)
          edz = (-ev1/sqr(8)) + ev1/sqr(9)
          eqzz = ev2/sqr(10) + ev2/sqr(11) - ev1/sqr(12)
          eqxx = ev1/sqr(13) - ev1/sqr(12)
          dxdx = ev1/sqr(14) - ev1/sqr(15)
          dzdz = ev2/sqr(16) + ev2/sqr(17) - ev2/sqr(18) - ev2/sqr(19)
          dzqxx = ev2/sqr(20) - ev2/sqr(21) - ev2/sqr(22) + ev2/sqr(23)
          qxxdz = ev2/sqr(24) - ev2/sqr(25) - ev2/sqr(26) + ev2/sqr(27)
          dzqzz = (-ev3/sqr(28)) + ev3/sqr(29) - ev3/sqr(30) + ev3/sqr(31) - &
            ev2/sqr(22) + ev2/sqr(20)
          qzzdz = (-ev3/sqr(32)) + ev3/sqr(33) - ev3/sqr(34) + ev3/sqr(35) + &
            ev2/sqr(24) - ev2/sqr(26)
          qxxqxx = ev3/sqr(37) + ev3/sqr(38) - ev2/sqr(39) - ev2/sqr(40) + ev2/&
            sqr(36)
          qxxqyy = ev2/sqr(41) - ev2/sqr(39) - ev2/sqr(40) + ev2/sqr(36)
          qxxqzz = ev3/sqr(43) + ev3/sqr(45) - ev3/sqr(42) - ev3/sqr(44) - ev2/&
            sqr(39) + ev2/sqr(36)
          qzzqxx = ev3/sqr(47) + ev3/sqr(49) - ev3/sqr(46) - ev3/sqr(48) - ev2/&
            sqr(40) + ev2/sqr(36)
          qzzqzz = ev4/sqr(50) + ev4/sqr(51) + ev4/sqr(52) + ev4/sqr(53) - ev3/&
            sqr(48) - ev3/sqr(46) - ev3/sqr(42) - ev3/sqr(44) + ev2/sqr(36)
          dxqxz = (-ev2/sqr(54)) + ev2/sqr(55) + ev2/sqr(56) - ev2/sqr(57)
          qxzdx = (-ev2/sqr(58)) + ev2/sqr(59) + ev2/sqr(60) - ev2/sqr(61)
          qxzqxz = ev3/sqr(65) - ev3/sqr(67) - ev3/sqr(69) + ev3/sqr(71) - ev3/&
            sqr(66) + ev3/sqr(68) + ev3/sqr(70) - ev3/sqr(72)
          ri(1) = ee
          ri(2) = -dze
          ri(3) = ee + qzze
          ri(4) = ee + qxxe
          ri(5) = -edz
          ri(6) = dzdz
          ri(7) = dxdx
          ri(8) = (-edz) - qzzdz
          ri(9) = (-edz) - qxxdz
          ri(10) = -qxzdx
          ri(11) = ee + eqzz
          ri(12) = ee + eqxx
          ri(13) = (-dze) - dzqzz
          ri(14) = (-dze) - dzqxx
          ri(15) = -dxqxz
          ri(16) = ee + eqzz + qzze + qzzqzz
          ri(17) = ee + eqzz + qxxe + qxxqzz
          ri(18) = ee + eqxx + qzze + qzzqxx
          ri(19) = ee + eqxx + qxxe + qxxqxx
          ri(20) = qxzqxz
          ri(21) = ee + eqxx + qxxe + qxxqyy
          ri(22) = pp*(qxxqxx - qxxqyy)
!
        end if
!
        if (l_feather) then
          call to_point(rij, point, const)
          ri(1) = ri(1)*const + (1.d0 - const)*point
          ri(2) = ri(2)*const
          ri(3) = ri(3)*const + (1.d0 - const)*point
          ri(4) = ri(4)*const + (1.d0 - const)*point
          ri(5) = ri(5)*const
          ri(6) = ri(6)*const
          ri(7) = ri(7)*const
          ri(8) = ri(8)*const
          ri(9) = ri(9)*const
          ri(10) = ri(10)*const
          ri(11) = ri(11)*const + (1.d0 - const)*point
          ri(12) = ri(12)*const + (1.d0 - const)*point
          ri(13) = ri(13)*const
          ri(14) = ri(14)*const
          ri(15) = ri(15)*const
          ri(16) = ri(16)*const + (1.d0 - const)*point
          ri(17) = ri(17)*const + (1.d0 - const)*point
          ri(18) = ri(18)*const + (1.d0 - const)*point
          ri(19) = ri(19)*const + (1.d0 - const)*point
          ri(20) = ri(20)*const
          ri(21) = ri(21)*const + (1.d0 - const)*point
          ri(22) = ri(22)*const
          gab = gab*const + (1.d0 - const)*point
!
!  Use the following block for checking the behavior of the electrostatic function
!
    !      do i = 10,80
    !        qa = i*0.1d0
    !        r = qa/a0
    !        gab = ev/sqrt(r*r + aee)
    !        call to_point(qa, point, const)
    !        da = gab*const + (1.d0 - const)*point
    !        write(iw,'(f8.1, 7f10.4)') qa, gab, da, da - gab, point
    !      end do
    !      stop
!
!  End of check
!
        end if
        ri = ri*nri
        return
      end subroutine reppd


      subroutine reppd2(ni, nj, r, ri, rep, core)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use parameters_C, only : tore, dorbs
      use mndod_C, only : indexd, ind2, isym
      use funcon_C, only : ev, a0
      use molkst_C, only : l_feather
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: ni
      integer  :: nj
      double precision  :: r
      double precision , intent(in) :: ri(22)
      double precision , intent(out) :: rep(491)
      double precision , intent(out) :: core(10,2)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(9) :: lorb
      integer , dimension(34) :: ipos
      integer :: i, lasti, lastk, li, j, lj, ij, k, lk, l, ll, kl, numb, &
        nold
      logical :: coul, coulomb
      double precision :: const, point
      double precision, external :: rijkl
!-----------------------------------------------
!     *
!     LOCAL TWO-CENTER TWO-ELECTRON INTEGRALS (SPD)
!     *
      data ipos/ 1, 5, 11, 12, 12, 2, 6, 13, 14, 14, 3, 8, 16, 18, 18, 7, 15, &
        10, 20, 4, 9, 17, 19, 21, 7, 15, 10, 20, 22, 4, 9, 17, 21, 19/
      data lorb/ 0, 3*1, 5*2/
!
      rep(:34) = ri(ipos)
      if (dorbs(ni) .or. dorbs(nj)) then
        if (dorbs(ni)) then
          lasti = 9
        else if (ni < 3) then
          lasti = 1
        else
          lasti = 4
        end if
!
        if (dorbs(nj)) then
          lastk = 9
        else if (nj < 3) then
          lastk = 1
        else
          lastk = 4
        end if
!
        do i = 1, lasti
          li = lorb(i)
          do j = 1, i
            coul = (i == j)
            lj = lorb(j)
            ij = indexd(i,j)
!
            do k = 1, lastk
              lk = lorb(k)
              do l = 1, k
                coulomb = (coul .and. k == l)
                ll = lorb(l)
                kl = indexd(k,l)
!
                numb = ind2(ij,kl)
                if (numb <= 34) cycle
                nold = isym(numb)
                select case (nold)
                case (35:)
                  rep(numb) = rep(nold)
                case (:(-35))
                  rep(numb) = -rep((-nold))
                case (0)
                  rep(numb) = rijkl(ni,nj,ij,kl,li,lj,lk,ll,0,r)*ev
                  if (l_feather) then
                    call to_point(r*a0, point, const)
                    if (coulomb) then
                      rep(numb) =rep(numb)*const + (1.d0 - const)*point
                    else
                      rep(numb) =rep(numb)*const
                    end if
                  end if
                end select
!
!      WRITE(11,'(I3,6X,''INT2C('',I2,'','',I2,'') ='',I5)')
!    -    NUMB,IJ,KL,NSYM
!       WRITE(11,'(6X,''INT2C('',I2,'','',I2,'') ='',I5)')IJ,KL,NSYM
!            IND2(IJ,KL) = NUMB
!            WRITE(6,'(2X,I4,'' <'',I2,I2,''|'',I2,I2,''>'',F12.4)')
!    -                  NUMB ,       I, J,       K, L,      REP(NUMB)
!
              end do
            end do
          end do
        end do
!
        core(5:10,1) = 0.D0
        core(5:10,2) = 0.D0
!
        if (dorbs(nj)) then
! --- <S S | D S>
          kl = indexd(5,1)
          core(5,2) = -rijkl(ni,nj,ij,kl,0,0,2,0,1,r)*ev*tore(ni)
! --- <S S | D P >
          kl = indexd(5,2)
          core(6,2) = -rijkl(ni,nj,ij,kl,0,0,2,1,1,r)*ev*tore(ni)
! --- <S S | D D >
          kl = indexd(5,5)
          core(7,2) = -rijkl(ni,nj,ij,kl,0,0,2,2,1,r)*ev*tore(ni)
! --- <S S | D+P+>
          kl = indexd(6,3)
          core(8,2) = -rijkl(ni,nj,ij,kl,0,0,2,1,1,r)*ev*tore(ni)
! --- <S S | D+D+>
          kl = indexd(6,6)
          core(9,2) = -rijkl(ni,nj,ij,kl,0,0,2,2,1,r)*ev*tore(ni)
! --- <S S | D#D#>
          kl = indexd(8,8)
          core(10,2) = -rijkl(ni,nj,ij,kl,0,0,2,2,1,r)*ev*tore(ni)
        end if
!*
        if (dorbs(ni)) then
! --- <D S | S S>
          kl = indexd(5,1)
          core(5,1) = -rijkl(ni,nj,kl,ij,2,0,0,0,2,r)*ev*tore(nj)
! --- <D P | S S >
          kl = indexd(5,2)
          core(6,1) = -rijkl(ni,nj,kl,ij,2,1,0,0,2,r)*ev*tore(nj)
! --- <D D | S S >
          kl = indexd(5,5)
          core(7,1) = -rijkl(ni,nj,kl,ij,2,2,0,0,2,r)*ev*tore(nj)
! --- <D+P+| S S >
          kl = indexd(6,3)
          core(8,1) = -rijkl(ni,nj,kl,ij,2,1,0,0,2,r)*ev*tore(nj)
! --- <D+D+| S S >
          kl = indexd(6,6)
          core(9,1) = -rijkl(ni,nj,kl,ij,2,2,0,0,2,r)*ev*tore(nj)
! --- <D#D#| S S >
          kl = indexd(8,8)
          core(10,1) = -rijkl(ni,nj,kl,ij,2,2,0,0,2,r)*ev*tore(nj)
!
        end if
      end if
!*    WRITE(6,'('' DCORE:'',/(2X,6F12.4))') CORE
      return
      end subroutine reppd2


      double precision function rijkl (ni, nj, ij, kl, li, lj, lk, ll, ic, r)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use parameters_C, only : ddp, po
      use mndod_C, only : indx, ch
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: ni
      integer , intent(in) :: nj
      integer , intent(in) :: ij
      integer , intent(in) :: kl
      integer , intent(in) :: li
      integer , intent(in) :: lj
      integer , intent(in) :: lk
      integer , intent(in) :: ll
      integer , intent(in) :: ic
      double precision  :: r
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: l1min, l1max, lij, l2min, l2max, lkl, l1, l2, lmin, m, mm
      double precision :: sum, pij, dij, pkl, dkl, add, s1, ccc
      double precision, external :: charg
!-----------------------------------------------
!     *
!
      pij = 0.d0
      pkl = 0.d0
      l1min = iabs(li - lj)
      l1max = li + lj
      lij = indx(li+1,lj+1)
      l2min = iabs(lk - ll)
      l2max = lk + ll
      lkl = indx(lk+1,ll+1)
      l1max = min(l1max,2)
      l1min = min(l1min,2)
      l2max = min(l2max,2)
      l2min = min(l2min,2)
      sum = 0.D00
!
      do l1 = l1min, l1max
        if (l1 == 0) then
          select case (lij)
          case (1)
            pij = po(1,ni)
            if (ic == 1) pij = po(9,ni)
          case (3)
            pij = po(7,ni)
          case (6)
            pij = po(8,ni)
          end select
        else
          dij = ddp(lij,ni)
          pij = po(lij,ni)
        end if
!
        do l2 = l2min, l2max
          if (l2 == 0) then
            select case (lkl)
            case (1)
              pkl = po(1,nj)
              if (ic == 2) pkl = po(9,nj)
            case (3)
              pkl = po(7,nj)
            case (6)
              pkl = po(8,nj)
            end select
          else
            dkl = ddp(lkl,nj)
            pkl = po(lkl,nj)
          end if
!
          add = (pij + pkl)**2
          lmin = min(l1,l2)
          s1 = 0.D00
          do m = -lmin, lmin
            ccc = ch(ij,l1,m)*ch(kl,l2,m)
            if (ccc == 0.000) cycle
            mm = iabs(m)
            s1 = s1 + charg(r,l1,l2,mm,dij,dkl,add)*ccc
          end do
          sum = sum + s1
        end do
      end do
!
      rijkl = sum
      return
      end function rijkl


      subroutine rotatd(ni, nj, ci, cj, w, kr, enuc)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use funcon_C, only : a0, ev
      use mndod_C, only : indx, indexd, sp, sd, pp, dp, d_d, cored, inddd
      use parameters_C, only : natorb, iod, tore
      use molkst_C, only : numcal, l_feather, method_PM7
!     *
!     CALCULATION OF TWO-CENTER TWO-ELECTRON INTEGRALS
!     IN THE MOLECULAR COODINATE SYSTEM BY 2.d0-STEP PROCEDURE
!     *
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
      integer, intent(in) :: ni, nj
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: kr
      double precision  :: ci(3)
      double precision  :: cj(3)
      double precision, dimension(*)  :: w
      double precision  :: enuc
      double precision, dimension(171) :: en
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(45) :: met
      integer :: limkl, kl, ii, kk, limij, istep, i, i1, j1, ij, jj, mm, k, l, &
        iw, iminus, j, ij1, jw, icalcn = 0, li, lj
      double precision, dimension(45,45) :: v
      double precision, dimension(22) :: ri
      double precision, dimension(491) :: rep
      double precision, dimension(2025) :: ww
      double precision :: r, wrepp, cc, gab, sum, rij, const, point
      logical, dimension(45,45) :: logv
!-----------------------------------------------
!
      data met/ 1, 2, 3, 2, 3, 3, 2, 3, 3, 3, 4, 5, 5, 5, 6, 4, 5, 5, 5, 6, 6, &
        4, 5, 5, 5, 6, 6, 6, 4, 5, 5, 5, 6, 6, 6, 6, 4, 5, 5, 5, 5*6/
!
      if (icalcn /= numcal) then
        icalcn = numcal
      end if
!
      call rotmat (nj, ni, ci, cj, r)
!
      call reppd (ni, nj, r, ri, gab) ! Work out the 22 s-p terms
      rij = r
      r = r/a0
      call spcore (ni, nj, r, cored)
      if (l_feather) then
        call to_point(rij, point, const)
      else
        const = 1.d0
        point = 0.d0
      end if
      call reppd2 (ni, nj, r, ri, rep, cored) ! add in "d" terms
      point = -(eV/r)*tore(nj)
      cored(1,1) = cored(1,1)*const + (1.d0 - const)*point
      cored(2,1) = cored(2,1)*const
      cored(3,1) = cored(3,1)*const + (1.d0 - const)*point
      cored(4,1) = cored(4,1)*const + (1.d0 - const)*point
      cored(5,1) = cored(5,1)*const
      cored(6,1) = cored(6,1)*const
      cored(7,1) = cored(7,1)*const + (1.d0 - const)*point
      cored(8,1) = cored(8,1)*const
      cored(9,1) = cored(9,1)*const + (1.d0 - const)*point
      cored(10,1) = cored(10,1)*const + (1.d0 - const)*point
      point = -(eV/r)*tore(ni)
      cored(1,2) = cored(1,2)*const + (1.d0 - const)*point
      cored(2,2) = cored(2,2)*const
      cored(3,2) = cored(3,2)*const + (1.d0 - const)*point
      cored(4,2) = cored(4,2)*const + (1.d0 - const)*point
      cored(5,2) = cored(5,2)*const
      cored(6,2) = cored(6,2)*const
      cored(7,2) = cored(7,2)*const + (1.d0 - const)*point
      cored(8,2) = cored(8,2)*const
      cored(9,2) = cored(9,2)*const + (1.d0 - const)*point
      cored(10,2) = cored(10,2)*const + (1.d0 - const)*point

!
      ii = natorb(ni)
      kk = natorb(nj)
      ww = 0
      if (ii*kk > 0) then
        limij = indx(ii,ii)
        limkl = indx(kk,kk)
        istep = limkl*limij
        ww(:istep) = 0.0D0
  !
        call tx (ii, kk, rep, logv, v)
  !
        do i1 = 1, ii
          do j1 = 1, i1
  !
            ij = indexd(i1,j1)
            jj = indx(i1,j1)
            mm = met(jj)
  !
            do k = 1, kk
              do l = 1, k
                kl = indx(k,l)
                if (.not.logv(ij,kl)) cycle
                wrepp = v(ij,kl)
  !     GO TO (1,2,3,4,5,6),MM
  !
                select case (mm)
                case (1)
                  iw = indw(1,1)
                  ww(iw) = wrepp
                case (2)
                  do i = 1, 3
                    iw = indw(i + 1,1)
                    ww(iw) = ww(iw) + sp(i1-1,i)*wrepp
                  end do
                case (3)
                  do i = 1, 3
                    cc = pp(i,i1-1,j1-1)
                    iw = indw(i + 1,i + 1)
                    ww(iw) = ww(iw) + cc*wrepp
                    iminus = i - 1
                    if (iminus == 0) cycle
                    do j = 1, iminus
                      cc = pp(1+i+j,i1-1,j1-1)
                      iw = indw(i + 1,j + 1)
                      ww(iw) = ww(iw) + cc*wrepp
                    end do
                  end do
                case (4)
                  do i = 1, 5
                    iw = indw(i + 4,1)
                    ww(iw) = ww(iw) + sd(i1-4,i)*wrepp
                  end do
                case (5)
                  do i = 1, 5
                    do j = 1, 3
                      iw = indw(i + 4,j + 1)
                      ij1 = 3*(i - 1) + j
                      ww(iw) = ww(iw) + dp(ij1,i1-4,j1-1)*wrepp
                    end do
                  end do
                case (6)
                  do i = 1, 5
                    cc = d_d(i,i1-4,j1-4)
                    iw = indw(i + 4,i + 4)
                    ww(iw) = ww(iw) + cc*wrepp
                    iminus = i - 1
                    if (iminus == 0) cycle
                    do j = 1, iminus
                      ij1 = inddd(i,j)
                      cc = d_d(ij1,i1-4,j1-4)
                      iw = indw(i + 4,j + 4)
                      ww(iw) = ww(iw) + cc*wrepp
                    end do
                  end do
                end select
               end do
            end do
          end do
        end do
      end if
      if (method_PM7) then
!
!  Elements with "d"-electrons have an imbalance between the nuclear-nuclear term and the other
!  two: electron-nuclear and electron-electron. To correct this specific error, the average value
!  of the two-electron two-center "d" coulomb integrals is set equal to the <ss|ss> integral.
! (The <ss|ss> integral is used in the nuclear-nuclear term)
!
        if (iod(ni) > 0) then
          sum = 0
          if (natorb(nj) == 9) then
              k = 45
            else if (natorb(nj) == 4) then
              k = 10
            else
              k = 1
            end if
          if (natorb(nj) > 1) then
            do i = 5,9                      !  "p" on nj with "d" on ni
              j = k*((i*(i + 1))/2 - 1)
              sum = sum + ww(j + 3) + ww(j + 6) + ww(j + 10)
            end do
            sum = (ww(1) - sum/15.d0)
            do i = 5,9
              j = k*((i*(i + 1))/2 - 1)
              do l = 2,4
                ww((l*(l+1))/2 + j) = ww((l*(l+1))/2 + j) + sum
              end do
            end do
            sum = cored(1,1) - (cored(3,1) + 2.d0*cored(4,1))/3.d0
            cored(3,1) = cored(3,1) + sum
            cored(4,1) = cored(4,1) + sum
          end if
          sum = 0.d0
          do i = 5,9
              sum = sum +ww(k*((i*(i + 1))/2 - 1) + 1)
          end do
          sum = (ww(1) - sum/5.d0)
          do i = 5,9                        !  "s" on nj with "d" on ni
            ww(k*((i*(i + 1))/2 - 1) + 1) = ww(k*((i*(i + 1))/2 - 1) + 1) + sum
          end do
          sum = cored(1,1) - (cored(7,1) + 2.d0*cored(9,1) + 2.d0*cored(10,1))/5.d0
          cored( 7,1) = cored( 7,1) + sum
          cored( 9,1) = cored( 9,1) + sum
          cored(10,1) = cored(10,1) + sum
        end if
        if (iod(nj) > 0) then
          sum = 0.d0
          if (iod(ni) > 0) then
            do i = 5, 9                     !  "d" with "d"  -  common to both atoms
              j = 45*((i*(i + 1))/2 - 1)
              sum = sum + (ww(j + 15) + ww(j + 21) + ww(j + 28) + ww(j + 36) + ww(j + 45))
            end do
            sum = (ww(1) - sum/25.d0)
            do i = 5, 9
              j = 45*((i*(i + 1))/2 - 1)
              do k = 5,9
                ww((k*(k+1))/2 + j) = ww((k*(k+1))/2 + j) + sum
              end do
            end do
          end if
          if (natorb(ni) > 1) then
            sum = 0.d0
            do i = 2,4                      !  "p" on ni with "d" on nj
              j = 45*((i*(i + 1))/2 - 1)
              sum = sum + (ww(j + 15) + ww(j + 21) + ww(j + 28) + ww(j + 36) + ww(j + 45))
            end do
            sum = (ww(1) - sum/15.d0)

            do i = 2,4
              j = 45*((i*(i + 1))/2 - 1)
              do k = 5,9
                ww((k*(k+1))/2 + j) = ww((k*(k+1))/2 + j) + sum
              end do
            end do
            sum = cored(1,2) - (cored(3,2) + 2.d0*cored(4,2))/3.d0
            cored(3,2) = cored(3,2) + sum
            cored(4,2) = cored(4,2) + sum
          end if
          sum = (ww(1) - (ww(15) + ww(21) + ww(28) + ww(36) + ww(45))/5.d0)
          do k = 5,9                        !  "s" on ni with "d" on nj
            ww((k*(k+1))/2) = ww((k*(k+1))/2) + sum
          end do
          sum = (cored(1,2) - (cored(7,2) + 2.d0*cored(9,2) + 2.d0*cored(10,2))/5.d0)
          cored( 7,2) = cored( 7,2) + sum
          cored( 9,2) = cored( 9,2) + sum
          cored(10,2) = cored(10,2) + sum
        end if
      end if
      iw = (ii*(ii + 1))/2
      jw = (kk*(kk + 1))/2
      call w2mat (ww, w, kr, iw, jw)
      sum = 0.d0
      do i = 1,9
        do j = 1,9
          ii = 45*((i*(i+1))/2 - 1)
          jj = (j*(j+1))/2
          sum = sum + ww(ii+jj)
        end do
      end do
      en(:) = 0.d0
      li = natorb(ni)
      lj = natorb(nj)
      call elenuc (1, li, li + 1, li + lj, en)
! *** CORE-CORE REPULSIONS FOR MNDO.
      call ccrep (ni, nj, r, enuc, gab)
      return
      contains


      integer function indw (i, j)
      integer, intent(in) :: i
      integer, intent(in) :: j
      indw = (indx(i,j)-1)*limkl + kl
      return
      end function indw
      end subroutine rotatd


      subroutine rotmat(nj, ni, coordi, coordj, r)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use parameters_C, only : dorbs
      use mndod_C, only : sp, pp, sd, dp, d_d
!     *
!     ROTATION MATRIX FOR A GIVEN ATOM PAIR I-J (I.GT.J).
!     *
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nj
      integer , intent(in) :: ni
      double precision , intent(out) :: r
      double precision , intent(in) :: coordi(3)
      double precision , intent(in) :: coordj(3)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      double precision, parameter :: small = 1.0D-07
      double precision, parameter :: pt5sq3 = 0.8660254037841D0
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k
      double precision, dimension(3,3) :: p
      double precision, dimension(5,5) :: d
      double precision :: x11, x22, x33, b, sqb, sb, ca, sa, cb, c2a, c2b, s2a, s2b
!-----------------------------------------------
! *** CALCULATE GEOMETRIC DATA AND INTERATOMIC DISTANCE.
!     CA  = COS(PHI)    , SA  = SIN(PHI)
!     CB  = COS(THETA)  , SB  = SIN(THETA)
!     C2A = COS(2*PHI)  , S2A = SIN(2*PHI)
!     C2B = COS(2*THETA), S2B = SIN(2*PHI)
      x11 = coordj(1) - coordi(1)
      x22 = coordj(2) - coordi(2)
      x33 = coordj(3) - coordi(3)
      b = x11*x11 + x22*x22
      r = sqrt(b + x33*x33)
      sqb = sqrt(b)
      sb = sqb/r
!     CHECK FOR SPECIAL CASE (BOTH ATOMS ON Z AXIS).
      if (sb > small) then
        ca = x11/sqb
        sa = x22/sqb
        cb = x33/r
      else
        sa = 0.d0
        sb = 0.d0
        if (x33 < 0.d0) then
          ca = -1.d0
          cb = -1.d0
        else if (x33 > 0.d0) then
          ca = 1.d0
          cb = 1.d0
        else
          ca = 0.d0
          cb = 0.d0
        end if
      end if
!     CONVERT DISTANCE TO ATOMIC UNITS.
!#      R      = R/A0
! *** CALCULATE ROTATION MATRIX ELEMENTS.
      p(1,1) = ca*sb
      p(2,1) = ca*cb
      p(3,1) = -sa
      p(1,2) = sa*sb
      p(2,2) = sa*cb
      p(3,2) = ca
      p(1,3) = cb
      p(2,3) = -sb
      p(3,3) = 0.d0
      if (dorbs(ni) .or. dorbs(nj)) then
        c2a = 2.d0*ca*ca - 1.d0
        c2b = 2.d0*cb*cb - 1.d0
        s2a = 2.d0*sa*ca
        s2b = 2.d0*sb*cb
        d(1,1) = pt5sq3*c2a*sb*sb
        d(2,1) = 0.5d0*c2a*s2b
        d(3,1) = -s2a*sb
        d(4,1) = c2a*(cb*cb + 0.5d0*sb*sb)
        d(5,1) = -s2a*cb
        d(1,2) = pt5sq3*ca*s2b
        d(2,2) = ca*c2b
        d(3,2) = -sa*cb
        d(4,2) = -0.5d0*ca*s2b
        d(5,2) = sa*sb
        d(1,3) = cb*cb - 0.5d0*sb*sb
        d(2,3) = -pt5sq3*s2b
        d(3,3) = 0.d0
        d(4,3) = pt5sq3*sb*sb
        d(5,3) = 0.d0
        d(1,4) = pt5sq3*sa*s2b
        d(2,4) = sa*c2b
        d(3,4) = ca*cb
        d(4,4) = -0.5d0*sa*s2b
        d(5,4) = -ca*sb
        d(1,5) = pt5sq3*s2a*sb*sb
        d(2,5) = 0.5d0*s2a*s2b
        d(3,5) = c2a*sb
        d(4,5) = s2a*(cb*cb + 0.5d0*sb*sb)
        d(5,5) = c2a*cb
      end if
!     *
!
!  S-P
      sp = p
!  P-P
!     DATA INDPP /1,4,5,4,2,6,5,6,3/
!     DATA INDDP /1,4,7,10,13,2,5,8,11,14,3,6,9,12,15/
      do k = 1, 3
        pp(1,k,k) = p(k,1)*p(k,1)
        pp(2,k,k) = p(k,2)*p(k,2)
        pp(3,k,k) = p(k,3)*p(k,3)
        pp(4,k,k) = p(k,1)*p(k,2)
        pp(5,k,k) = p(k,1)*p(k,3)
        pp(6,k,k) = p(k,2)*p(k,3)
        if (k == 1) cycle
        pp(1,k,:k-1) = 2.D0*p(k,1)*p(:k-1,1)
        pp(2,k,:k-1) = 2.D0*p(k,2)*p(:k-1,2)
        pp(3,k,:k-1) = 2.D0*p(k,3)*p(:k-1,3)
        pp(4,k,:k-1) = p(k,1)*p(:k-1,2) + p(k,2)*p(:k-1,1)
        pp(5,k,:k-1) = p(k,1)*p(:k-1,3) + p(k,3)*p(:k-1,1)
        pp(6,k,:k-1) = p(k,2)*p(:k-1,3) + p(k,3)*p(:k-1,2)
      end do
!
      if (dorbs(ni) .or. dorbs(nj)) then
!  S-D
        sd = d
!  D-P
        do k = 1, 5
          dp(1,k,:) = d(k,1)*p(:,1)
          dp(2,k,:) = d(k,1)*p(:,2)
          dp(3,k,:) = d(k,1)*p(:,3)
          dp(4,k,:) = d(k,2)*p(:,1)
          dp(5,k,:) = d(k,2)*p(:,2)
          dp(6,k,:) = d(k,2)*p(:,3)
          dp(7,k,:) = d(k,3)*p(:,1)
          dp(8,k,:) = d(k,3)*p(:,2)
          dp(9,k,:) = d(k,3)*p(:,3)
          dp(10,k,:) = d(k,4)*p(:,1)
          dp(11,k,:) = d(k,4)*p(:,2)
          dp(12,k,:) = d(k,4)*p(:,3)
          dp(13,k,:) = d(k,5)*p(:,1)
          dp(14,k,:) = d(k,5)*p(:,2)
          dp(15,k,:) = d(k,5)*p(:,3)
        end do
!  D-D
        do k = 1, 5
          d_d(1,k,k) = d(k,1)*d(k,1)
          d_d(2,k,k) = d(k,2)*d(k,2)
          d_d(3,k,k) = d(k,3)*d(k,3)
          d_d(4,k,k) = d(k,4)*d(k,4)
          d_d(5,k,k) = d(k,5)*d(k,5)
          d_d(6,k,k) = d(k,1)*d(k,2)
          d_d(7,k,k) = d(k,1)*d(k,3)
          d_d(8,k,k) = d(k,2)*d(k,3)
          d_d(9,k,k) = d(k,1)*d(k,4)
          d_d(10,k,k) = d(k,2)*d(k,4)
          d_d(11,k,k) = d(k,3)*d(k,4)
          d_d(12,k,k) = d(k,1)*d(k,5)
          d_d(13,k,k) = d(k,2)*d(k,5)
          d_d(14,k,k) = d(k,3)*d(k,5)
          d_d(15,k,k) = d(k,4)*d(k,5)
          if (k == 1) cycle
          d_d(1,k,:k-1) = 2.D0*d(k,1)*d(:k-1,1)
          d_d(2,k,:k-1) = 2.D0*d(k,2)*d(:k-1,2)
          d_d(3,k,:k-1) = 2.D0*d(k,3)*d(:k-1,3)
          d_d(4,k,:k-1) = 2.D0*d(k,4)*d(:k-1,4)
          d_d(5,k,:k-1) = 2.D0*d(k,5)*d(:k-1,5)
          d_d(6,k,:k-1) = d(k,1)*d(:k-1,2) + d(k,2)*d(:k-1,1)
          d_d(7,k,:k-1) = d(k,1)*d(:k-1,3) + d(k,3)*d(:k-1,1)
          d_d(8,k,:k-1) = d(k,2)*d(:k-1,3) + d(k,3)*d(:k-1,2)
          d_d(9,k,:k-1) = d(k,1)*d(:k-1,4) + d(k,4)*d(:k-1,1)
          d_d(10,k,:k-1) = d(k,2)*d(:k-1,4) + d(k,4)*d(:k-1,2)
          d_d(11,k,:k-1) = d(k,3)*d(:k-1,4) + d(k,4)*d(:k-1,3)
          d_d(12,k,:k-1) = d(k,1)*d(:k-1,5) + d(k,5)*d(:k-1,1)
          d_d(13,k,:k-1) = d(k,2)*d(:k-1,5) + d(k,5)*d(:k-1,2)
          d_d(14,k,:k-1) = d(k,3)*d(:k-1,5) + d(k,5)*d(:k-1,3)
          d_d(15,k,:k-1) = d(k,4)*d(:k-1,5) + d(k,5)*d(:k-1,4)
        end do
      end if
      return
      end subroutine rotmat


      double precision function rsc (k, na, ea, nb, eb, nc, ec, nd, ed)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE mndod_C, only : fx, b
      use funcon_C, only : ev
!     *
!     CALCULATE THE RADIAL PART OF ONE-CENTER TWO-ELECTRON INTEGRALS
!     (SLATER-CONDON PARAMETER)
!     K- TYPE OF INTEGRAL ,   CAN BE EQUAL TO 0,1,2,3,4 IN SPD-BASIS
!     NA,NB -PRINCIPLE QUANTUM NUMBER OF AO,CORRESPONDING ELECTRON 1
!     EA,EB -EXPONENTS OF AO,CORRESPONDING ELECTRON 1
!     NC,ND -PRINCIPLE QUANTUM NUMBER OF AO,CORRESPONDING ELECTRON 2
!     EC,ED -EXPONENTS OF AO,CORRESPONDING ELECTRON 2
!     *
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: k
      integer , intent(in) :: na
      integer , intent(in) :: nb
      integer , intent(in) :: nc
      integer , intent(in) :: nd
      double precision , intent(in) :: ea
      double precision , intent(in) :: eb
      double precision , intent(in) :: ec
      double precision , intent(in) :: ed
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nab, ncd, n, m, i, m1, m2
      double precision :: aea, aeb, aec, aed, ecd, eab, e, ae, a2, acd, aab, ff, c, s0, s1, s2, s3
!-----------------------------------------------
      aea = log(ea)
      aeb = log(eb)
      aec = log(ec)
      aed = log(ed)
      nab = na + nb
      ncd = nc + nd
      ecd = ec + ed
      eab = ea + eb
      e = ecd + eab
      n = nab + ncd
      ae = log(e)
      a2 = log(2.d0)
      acd = log(ecd)
      aab = log(eab)
      ff = fx(n)/sqrt(fx(2*na+1)*fx(2*nb+1)*fx(2*nc+1)*fx(2*nd+1))
      c = ev*ff*exp(na*aea + nb*aeb + nc*aec + nd*aed + 0.5d0*(aea + aeb + aec + &
        aed) + a2*(n + 2) - ae*n)
      s0 = 1.d0/e
      s1 = 0.d0
      s2 = 0.d0
      m = ncd - k
      do i = 1, m
        s0 = s0*e/ecd
        s1 = s1 + s0*(b(ncd-k,i)-b(ncd+k+1,i))/b(n,i)
      end do
      m1 = m + 1
      m2 = ncd + k + 1
      do i = m1, m2
        s0 = s0*e/ecd
        s2 = s2 + s0*b(m2,i)/b(n,i)
      end do
      s3 = exp(ae*n - acd*m2 - aab*(nab - k))/b(n,m2)
      rsc = c*(s1 - s2 + s3)
      return
      end function rsc


      subroutine scprm(ni, r066, r266, r466, r016, r244, r036, r236, r155, r355, r125, r234, r246)
      USE mndod_C, only : iii, iiid
      use parameters_C, only : zsn, zpn, zdn
      implicit none
      integer , intent(in) :: ni
      double precision , intent(out) :: r066
      double precision , intent(out) :: r266
      double precision , intent(out) :: r466
      double precision , intent(out) :: r016
      double precision , intent(out) :: r244
      double precision , intent(out) :: r036
      double precision , intent(out) :: r236
      double precision , intent(out) :: r155
      double precision , intent(out) :: r355
      double precision , intent(out) :: r125
      double precision , intent(out) :: r234
      double precision , intent(out) :: r246
      integer :: ns, nd
      double precision :: es, ep, ed
      double precision, external :: rsc
!-----------------------------------------------
      ns = iii(ni)
      nd = iiid(ni)
      es = zsn(ni)
      ep = zpn(ni)
      ed = zdn(ni)
      r016 = rsc(0,ns,es,ns,es,nd,ed,nd,ed)
      r036 = rsc(0,ns,ep,ns,ep,nd,ed,nd,ed)
      r066 = rsc(0,nd,ed,nd,ed,nd,ed,nd,ed)
      r155 = rsc(1,ns,ep,nd,ed,ns,ep,nd,ed)
      r125 = rsc(1,ns,es,ns,ep,ns,ep,nd,ed)
      r244 = rsc(2,ns,es,nd,ed,ns,es,nd,ed)
      r236 = rsc(2,ns,ep,ns,ep,nd,ed,nd,ed)
      r266 = rsc(2,nd,ed,nd,ed,nd,ed,nd,ed)
      r234 = rsc(2,ns,ep,ns,ep,ns,es,nd,ed)
      r246 = rsc(2,ns,es,nd,ed,nd,ed,nd,ed)
      r355 = rsc(3,ns,ep,nd,ed,ns,ep,nd,ed)
      r466 = rsc(4,nd,ed,nd,ed,nd,ed,nd,ed)
      return
      end subroutine scprm


      subroutine spcore(ni, nj, r, core)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use parameters_C, only : tore, po, ddp
      use funcon_C, only : ev
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: ni
      integer , intent(in) :: nj
      double precision , intent(in) :: r
      double precision , intent(out) :: core(10,2)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      double precision, dimension(7) :: pxy
      double precision, dimension(4) :: ai, aj
      double precision, dimension(7) :: xi, xj
      double precision :: r2, aci, acj, ssi, ssj, ppj, da, qa, twoqa, adj, aqj, ppi&
        , db, qb, adi, aqi, twoqb
!-----------------------------------------------
!     *
!     CALCULATE THE NUCLEAR ATTRACTION INTEGRALS IN LOCAL COORDINATES.
!     *
      data pxy/ 1.0d0, -0.5d0, -0.5d0, 0.5d0, 0.25d0, 0.25d0, 0.5d0/
      core = 0.d0
! *** INITIALIZATION.
      r2 = r*r
      aci = po(9,ni)
      acj = po(9,nj)
! *** SS -CORE INTERACTION
      ssi = (aci + po(1,nj))**2
      ssj = (acj + po(1,ni))**2
      core(1,1) = -tore(nj)*ev/sqrt(r2 + ssj)
      core(1,2) = -tore(ni)*ev/sqrt(r2 + ssi)
      if (ni>=3 .or. nj>=3) then
! *** NI -  HEAVY ATOM
        if (ni >= 3) then
          ppj = (acj + po(7,ni))**2
          da = ddp(2,ni)
          qa = ddp(3,ni)/sqrt(2.d0)
          twoqa = qa + qa
          adj = (po(2,ni)+acj)**2
          aqj = (po(3,ni)+acj)**2
          xj(1) = r2 + ppj
          xj(2) = r2 + aqj
          xj(3) = (r + da)**2 + adj
          xj(4) = (r - da)**2 + adj
          xj(5) = (r - twoqa)**2 + aqj
          xj(6) = (r + twoqa)**2 + aqj
          xj(7) = r2 + twoqa*twoqa + aqj
          do i = 1, 7
            xj(i) = pxy(i)/sqrt(xj(i))
          end do
          aj(2) = (xj(3)+xj(4))*ev
          aj(3) = (xj(1)+xj(2)+xj(5)+xj(6))*ev
          aj(4) = (xj(1)+xj(2)+xj(7))*ev
          core(2,1) = -tore(nj)*aj(2)
          core(3,1) = -tore(nj)*aj(3)
          core(4,1) = -tore(nj)*aj(4)
        end if
! *** NJ- HEAVY ATOM
        if (nj >= 3) then
          ppi = (aci + po(7,nj))**2
          db = ddp(2,nj)
          qb = ddp(3,nj)/sqrt(2.d0)
          adi = (po(2,nj)+aci)**2
          aqi = (po(3,nj)+aci)**2
          twoqb = qb + qb
          xi(1) = r2 + ppi
          xi(2) = r2 + aqi
          xi(3) = (r + db)**2 + adi
          xi(4) = (r - db)**2 + adi
          xi(5) = (r - twoqb)**2 + aqi
          xi(6) = (r + twoqb)**2 + aqi
          xi(7) = r2 + twoqb*twoqb + aqi
          do i = 1, 7
            xi(i) = pxy(i)/sqrt(xi(i))
          end do
          ai(2) = -(xi(3)+xi(4))*ev
          ai(3) = (xi(1)+xi(2)+xi(5)+xi(6))*ev
          ai(4) = (xi(1)+xi(2)+xi(7))*ev
          core(2,2) = -tore(ni)*ai(2)
          core(3,2) = -tore(ni)*ai(3)
          core(4,2) = -tore(ni)*ai(4)
        end if
      end if
      return
      end subroutine spcore


      subroutine eiscor(r016, r066, r244, r266, r466, ni)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use parameters_C, only : eisol
!***********************************************************************
!
!    EISCOR adds in the ONE-CENTER terms for the atomic energy of
!           those atoms that have partly filled "d" shells
!
!***********************************************************************
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: ni
      double precision , intent(in) :: r016
      double precision , intent(in) :: r066
      double precision , intent(in) :: r244
      double precision , intent(in) :: r266
      double precision , intent(in) :: r466
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(100) :: ir016, ir066, ir244, ir266, ir466
      integer :: i
!-----------------------------------------------
      data (ir016(i),ir066(i),ir244(i),ir266(i),ir466(i),i=1,20)/ 100*0/
!
!                                Sc  Ti   V  Cr  Mn  Fe  Co  Ni  Cu
! Atomic orbital population   4s  2   2   2   1   2   2   2   2   1
! of gaseous atom             3d  1   2   3   5   5   6   7   8  10
!
! State term:                    2D  3F  4F  7S  6S  5D  4F  3F  2S
!
      data (ir016(i),i=21,29)/    2,  4,  6,  5, 10, 12, 14, 16, 10/
      data (ir066(i),i=21,29)/    0,  1,  3, 10, 10, 15, 21, 28, 45/
      data (ir244(i),i=21,29)/    1,  2,  3,  5,  5,  6,  7,  8,  5/
      data (ir266(i),i=21,29)/    0,  8, 15, 35, 35, 35, 43, 50, 70/
      data (ir466(i),i=21,29)/    0,  1,  8, 35, 35, 35, 36, 43, 70/
!
      data (ir016(i),ir066(i),ir244(i),ir266(i),ir466(i),i=30,38)/ 45*0/
!
!                                 Y  Zr  Nb  Mo  Tc  Ru  Rh  Pd  Ag
! Atomic orbital population   5s  2   2   1   1   2   1   1   0   1
! of gaseous atom             4d  1   2   4   5   5   7   8  10  10
!
! State term:                    2D  3F  6D  7S  6D  5F  4F  1S  2S
!
      data (ir016(i),i=39,47)/    2,  4,  4,  5, 10,  7,  8,  0, 10/
      data (ir066(i),i=39,47)/    0,  1,  6, 10, 10, 21, 28, 45, 45/
      data (ir244(i),i=39,47)/    1,  2,  4,  5,  5,  5,  5,  0,  5/
      data (ir266(i),i=39,47)/    0,  8, 21, 35, 35, 43, 50, 70, 70/
      data (ir466(i),i=39,47)/    0,  1, 21, 35, 35, 36, 43, 70, 70/
!
      data(ir016(i), ir066(i), ir244(i), ir266(i), ir466(i), i=48,56) / 45 * 0 /
      data ir016(57)  / 2 /
      data ir066(57)  / 0 /  !  Lanthanum
      data ir244(57)  / 1 /  !  For the methods used in MOPAC, use the
      data ir266(57)  / 0 /  !  configuration 6s(2)5d(1), not 6s(2)4f(1).
      data ir466(57)  / 0 /  !  This is a low-lying excited state.
      data(ir016(i), ir066(i), ir244(i), ir266(i), ir466(i), i=58, 70) / 65 * 0 /
      data ir016(71)  / 2 /
      data ir066(71)  / 0 /  !  Lutetium
      data ir244(71)  / 1 /  !
      data ir266(71)  / 0 /  !
      data ir466(71)  / 0 /  !
!                                    Hf  Ta   W  Re  Os  Ir  Pt  Au Hg
! Atomic orbital population   6s      2   2   1   2   2   2   1   1  2
! of gaseous atom             5d      2   3   5   5   6   7   9  10  0
!
! State term:                        3F  4F  7S  6S  5D  4F  3D  2S 1S
!
      data (ir016(i),i=72,80)/        4,  6,  5, 10, 12, 14,  9, 10, 0/
      data (ir066(i),i=72,80)/        1,  3, 10, 10, 15, 21, 36, 45, 0/
      data (ir244(i),i=72,80)/        2,  3,  5,  5,  6,  7,  5,  5, 0/
      data (ir266(i),i=72,80)/        8, 15, 35, 35, 35, 43, 56, 70, 0/
      data (ir466(i),i=72,80)/        1,  8, 35, 35, 35, 36, 56, 70, 0/
!
!
!
!     R016:  <SS|DD>
!     R066:  <DD|DD> "0" term
!     R244:  <SD|SD>
!     R266:  <DD|DD> "2" term
!     R466:  <DD|DD> "4" term
!
      eisol(ni) = eisol(ni) + ir016(ni)*r016 + ir066(ni)*r066 - ir244(ni)*r244/&
        5 - ir266(ni)*r266/49 - ir466(ni)*r466/49
      return
      end subroutine eiscor


      subroutine tx(ii, kk, rep, logv, v)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use mndod_C, only : indexd, indx, ind2, sp, sd, pp, dp, d_d
!     *
!     ROTATION OF TWO-ELECTRON TWO-CENTER INTEGRALS IN SPD BASIS
!     FIRST STEP
!     *
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: ii
      integer , intent(in) :: kk
      double precision , intent(in) :: rep(491)
      double precision , intent(out) :: v(45,45)
      logical , intent(out) :: logv(45,45)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(45) :: met
      integer :: limkl, k, i1, j1, ij, k1, l1, kl, nd, ll, mm, l
      double precision :: wrepp
!-----------------------------------------------
      data met/ 1, 2, 3, 2, 3, 3, 2, 3, 3, 3, 4, 5, 5, 5, 6, 4, 5, 5, 5, 6, 6, &
        4, 5, 5, 5, 6, 6, 6, 4, 5, 5, 5, 6, 6, 6, 6, 4, 5, 5, 5, 5*6/
!
      limkl = indx(kk,kk)
      logv(:,:limkl) = .FALSE.
      v(:,:limkl) = 0.D0
!
      do i1 = 1, ii
        do j1 = 1, i1
          ij = indexd(i1,j1)
!
          do k1 = 1, kk
!
            do l1 = 1, k1
              kl = indexd(k1,l1)
              nd = ind2(ij,kl)
              if (nd == 0) cycle
!
              wrepp = rep(nd)
              ll = indx(k1,l1)
              mm = met(ll)
!
              select case (mm)
              case (1)
                v(ij,1) = wrepp
              case (2)
                k = k1 - 1
                v(ij,2) = v(ij,2) + sp(k,1)*wrepp
                v(ij,4) = v(ij,4) + sp(k,2)*wrepp
                v(ij,7) = v(ij,7) + sp(k,3)*wrepp
              case (3)
                k = k1 - 1
                l = l1 - 1
                v(ij,3) = v(ij,3) + pp(1,k,l)*wrepp
                v(ij,6) = v(ij,6) + pp(2,k,l)*wrepp
                v(ij,10) = v(ij,10) + pp(3,k,l)*wrepp
                v(ij,5) = v(ij,5) + pp(4,k,l)*wrepp
                v(ij,8) = v(ij,8) + pp(5,k,l)*wrepp
                v(ij,9) = v(ij,9) + pp(6,k,l)*wrepp
                cycle
              case (4)
                k = k1 - 4
                v(ij,11) = v(ij,11) + sd(k,1)*wrepp
                v(ij,16) = v(ij,16) + sd(k,2)*wrepp
                v(ij,22) = v(ij,22) + sd(k,3)*wrepp
                v(ij,29) = v(ij,29) + sd(k,4)*wrepp
                v(ij,37) = v(ij,37) + sd(k,5)*wrepp
              case (5)
                k = k1 - 4
                l = l1 - 1
                v(ij,12) = v(ij,12) + dp(1,k,l)*wrepp
                v(ij,13) = v(ij,13) + dp(2,k,l)*wrepp
                v(ij,14) = v(ij,14) + dp(3,k,l)*wrepp
                v(ij,17) = v(ij,17) + dp(4,k,l)*wrepp
                v(ij,18) = v(ij,18) + dp(5,k,l)*wrepp
                v(ij,19) = v(ij,19) + dp(6,k,l)*wrepp
                v(ij,23) = v(ij,23) + dp(7,k,l)*wrepp
                v(ij,24) = v(ij,24) + dp(8,k,l)*wrepp
                v(ij,25) = v(ij,25) + dp(9,k,l)*wrepp
                v(ij,30) = v(ij,30) + dp(10,k,l)*wrepp
                v(ij,31) = v(ij,31) + dp(11,k,l)*wrepp
                v(ij,32) = v(ij,32) + dp(12,k,l)*wrepp
                v(ij,38) = v(ij,38) + dp(13,k,l)*wrepp
                v(ij,39) = v(ij,39) + dp(14,k,l)*wrepp
                v(ij,40) = v(ij,40) + dp(15,k,l)*wrepp
              case (6)
                k = k1 - 4
                l = l1 - 4
                v(ij,15) = v(ij,15) + d_d(1,k,l)*wrepp
                v(ij,21) = v(ij,21) + d_d(2,k,l)*wrepp
                v(ij,28) = v(ij,28) + d_d(3,k,l)*wrepp
                v(ij,36) = v(ij,36) + d_d(4,k,l)*wrepp
                v(ij,45) = v(ij,45) + d_d(5,k,l)*wrepp
                v(ij,20) = v(ij,20) + d_d(6,k,l)*wrepp
                v(ij,26) = v(ij,26) + d_d(7,k,l)*wrepp
                v(ij,27) = v(ij,27) + d_d(8,k,l)*wrepp
                v(ij,33) = v(ij,33) + d_d(9,k,l)*wrepp
                v(ij,34) = v(ij,34) + d_d(10,k,l)*wrepp
                v(ij,35) = v(ij,35) + d_d(11,k,l)*wrepp
                v(ij,41) = v(ij,41) + d_d(12,k,l)*wrepp
                v(ij,42) = v(ij,42) + d_d(13,k,l)*wrepp
                v(ij,43) = v(ij,43) + d_d(14,k,l)*wrepp
                v(ij,44) = v(ij,44) + d_d(15,k,l)*wrepp
!
              end select
            end do
          end do
          where (v(ij,:limkl) /= 0.00D00)
            logv(ij,:limkl) = .TRUE.
          end where
        end do
      end do
!
      return
      end subroutine tx


      subroutine w2mat(ww, w, kr, limij, limkl)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!     *
!     STORE TWO-CENTER TWO-ELECTRON INTEGRALS IN A SQUARE MATRIX.
!     *
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(inout) :: kr
      integer , intent(in) :: limij
      integer , intent(in) :: limkl
      double precision , intent(in) :: ww(limkl,limij)
      double precision , intent(out) :: w(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    integer :: ij, kl, l
   !
   ! ... Executable Statements ...
   !
    l = 0
    do ij = 1, limij
      do kl = 1, limkl
        l = l + 1
        w(l) = ww(kl, ij)
      end do
    end do
    kr = kr + l
      return
      end subroutine w2mat


   subroutine wstore (w, kr, ni, ilim)
   !     *
   !     COMPLETE DEFINITION OF SQUARE MATRIX OF MNDO TWO-ELECTRON
   !     INTEGRALS BY INCLUDING THE ONE-CENTER TERMS AND THE TERMS
   !     WITH TRANSPOSED INDICES.
   !     *
   !     *
    use mndod_C, only: intij, intkl, intrep, repd
    use parameters_C, only: gss, gsp, gpp, gp2, hsp, natorb
   !
   !.. Implicit Declarations ..
    implicit none
   !
   !.. Formal Arguments ..
    integer, intent (in) :: ilim
    double precision, dimension (ilim, ilim), intent (out) :: w
    integer, intent (inout) :: kr
    integer, intent (in) :: ni
   !
   !.. Local Scalars ..
    integer :: i, ij, ij0, Int, ip, ipx, ipy, ipz, j, kl
   !
   ! ... Executable Statements ...
   !
   ! *** INCLUDE NONZERO ONE-CENTER TERMS.
    do i = 1, ilim
      do j = 1, ilim
        w(i, j) = 0.d0
      end do
    end do
    ip = 1
    w(ip, ip) = gss(ni)
    if (natorb(ni) > 2) then
      ipx = ip + 2
      ipy = ip + 5
      ipz = ip + 9
      w(ipx, ip) = gsp(ni)
      w(ipy, ip) = gsp(ni)
      w(ipz, ip) = gsp(ni)
      w(ip, ipx) = gsp(ni)
      w(ip, ipy) = gsp(ni)
      w(ip, ipz) = gsp(ni)
      w(ipx, ipx) = gpp(ni)
      w(ipy, ipy) = gpp(ni)
      w(ipz, ipz) = gpp(ni)
      w(ipy, ipx) = gp2(ni)
      w(ipz, ipx) = gp2(ni)
      w(ipz, ipy) = gp2(ni)
      w(ipx, ipy) = gp2(ni)
      w(ipx, ipz) = gp2(ni)
      w(ipy, ipz) = gp2(ni)
      w(ip+1, ip+1) = hsp(ni)
      w(ip+3, ip+3) = hsp(ni)
      w(ip+6, ip+6) = hsp(ni)
      w(ip+4, ip+4) = 0.5d0*(gpp(ni)-gp2(ni))
      w(ip+7, ip+7) = 0.5d0*(gpp(ni)-gp2(ni))
      w(ip+8, ip+8) = 0.5d0*(gpp(ni)-gp2(ni))
      if (ilim > 10) then
        ij0 = ip - 1
        do i = 1, 243
          ij = intij(i)
          kl = intkl(i)
          Int = intrep(i)
          w(ij+ij0, kl+ij0) = repd(Int, ni)
        end do
      end if
    end if
    kr = kr + ilim ** 2
end subroutine wstore

      subroutine aijm(ni)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE mndod_C, only : aij,  iii, iiid, fx
      use parameters_C, only : zs, zp, zd, dorbs
!     *
!     AIJ-VALUES FOR EVALUATION OF TWO-CENTER TWO-ELECTRON INTEGRALS
!     AND OF HYBRID CONTRIBUTION TO DIPOLE MOMENT IN MNDO-D.
!     DEFINITION SEE EQUATION (7) OF TCA PAPER.
!     *
!     RESULTS ARE STORED IN ARRAY AIJ(6,5,107). CONVENTIONS:
!     FIRST  INDEX      1 SS, 2 SP, 3 PP, 4 SD, 5 PD, 6 DD.
!     SECOND INDEX      L+1 FROM DEFINITION OF MULTIPOLE.
!     THIRD  INDEX      ATOMIC NUMBER.
!     *
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: ni
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nsp, nd
      double precision :: z1, z2, z3, zz
!-----------------------------------------------
!     *
!
!
      z1 = zs(ni)
      z2 = zp(ni)
      z3 = zd(ni)
      nsp = iii(ni)
      if (ni < 3) return
      zz = z1*z2
      if (zz < 0.01d0) return
      aij(2,ni) = aijl(z1,z2,nsp,nsp,1)
      aij(3,ni) = aijl(z2,z2,nsp,nsp,2)
      if (dorbs(ni)) then
        nd = iiid(ni)
        aij(4,ni) = aijl(z1,z3,nsp,nd,2)
        aij(5,ni) = aijl(z2,z3,nsp,nd,1)
        aij(6,ni) = aijl(z3,z3,nd,nd,2)
      end if
      return
      contains


      double precision function aijl (z1, z2, n1, n2, l)
      double precision, intent(in) :: z1
      double precision, intent(in) :: z2
      integer, intent(in) :: n1
      integer, intent(in) :: n2
      integer, intent(in) :: l
      double precision :: zz
      zz = z1 + z2 + 1.d-20
      aijl = fx(n1+n2+l+1)/sqrt(fx(2*n1+1)*fx(2*n2+1))*(2*z1/zz)**n1*sqrt(2&
        *z1/zz)*(2*z2/zz)**n2*sqrt(2*z2/zz)*2**l/zz**l
      return
      end function aijl
      end subroutine aijm

      double precision function charg (r, l1, l2, m, da, db, add)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!     *
!     INTERACTION BETWEEN 2 POINT-CHARGE CONFIGURATIONS (MNDO-D).
!     *
!     R      DISTANCE IN ATOMIC UNITS.
!     L1,M   QUANTUM NUMBERS FOR MULTIPOLE OF CONFIGURATION 1.
!     L2,M   QUANTUM NUMBERS FOR MULTIPOLE OF CONFIGURATION 2.
!     DA     CHARGE SEPARATION OF CONFIGURATION 1.
!     DB     CHARGE SEPARATION OF CONFIGURATION 2.
!     ADD    ADDITIVE TERM
!     *
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: l1
      integer , intent(in) :: l2
      integer , intent(in) :: m
      double precision , intent(in) :: r
      double precision , intent(in) :: da
      double precision , intent(in) :: db
      double precision , intent(in) :: add
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      double precision :: dzdz, dxdx, qqzz, qzzq, dzqzz, qzzdz, zzzz, xyxy, ab, &
        dxqxz, aa, qxzdx, qxzqxz
!-----------------------------------------------
!
      charg = 0.0D00
!     Q - Q.
      if (l1==0 .and. l2==0) then
        charg = 1.D00/sqrt(r**2 + add)
!     Z - Q.
      else if (l1==1 .and. l2==0) then
        charg = (-1.D00/sqrt((r + da)**2 + add)) + 1.D00/sqrt((r - da)**2 + add&
          )
        charg = charg/2.0D00
!     Q - Z.
      else if (l1==0 .and. l2==1) then
        charg = 1.D00/sqrt((r + db)**2 + add) - 1.D00/sqrt((r - db)**2 + add)
        charg = charg/2.D00
!     Z - Z.
      else if (l1==1 .and. l2==1 .and. m==0) then
        dzdz = 1.D00/sqrt((r + da - db)**2 + add) + 1.D00/sqrt((r - da + db)**2&
           + add) - 1.D00/sqrt((r - da - db)**2 + add) - 1.D00/sqrt((r + da + &
          db)**2 + add)
        charg = dzdz/4.D00
!     X - X
      else if (l1==1 .and. l2==1 .and. m==1) then
        dxdx = 2.D00/sqrt(r**2 + (da - db)**2 + add) - 2.D00/sqrt(r**2 + (da + &
          db)**2 + add)
        charg = dxdx*0.25D00
!     Q - ZZ
      else if (l1==0 .and. l2==2) then
        qqzz = 1.D00/sqrt((r - db)**2 + add) - 2.D00/sqrt(r**2 + db**2 + add)&
           + 1.D00/sqrt((r + db)**2 + add)
        charg = qqzz/4.D00
!     ZZ -Q
      else if (l1==2 .and. l2==0) then
        qzzq = 1.D00/sqrt((r - da)**2 + add) - 2.D00/sqrt(r**2 + da**2 + add)&
           + 1.D00/sqrt((r + da)**2 + add)
        charg = qzzq/4.D00
!     Z - ZZ
      else if (l1==1 .and. l2==2 .and. m==0) then
        dzqzz = 1.D00/sqrt((r - da - db)**2 + add) - 2.D00/sqrt((r - da)**2 + &
          db**2 + add) + 1.D00/sqrt((r + db - da)**2 + add) - 1.D00/sqrt((r - &
          db + da)**2 + add) + 2.D00/sqrt((r + da)**2 + db**2 + add) - 1.D00/&
          sqrt((r + da + db)**2 + add)
        charg = dzqzz/8.D00
!     ZZ - Z
      else if (l1==2 .and. l2==1 .and. m==0) then
        qzzdz = (-1.D00/sqrt((r - da - db)**2 + add)) + 2.D00/sqrt((r - db)**2&
           + da**2 + add) - 1.D00/sqrt((r + da - db)**2 + add) + 1.D00/sqrt((r&
           - da + db)**2 + add) - 2.D00/sqrt((r + db)**2 + da**2 + add) + 1.D00&
          /sqrt((r + da + db)**2 + add)
        charg = qzzdz/8.D00
!     ZZ - ZZ
      else if (l1==2 .and. l2==2 .and. m==0) then
        zzzz = 1.D00/sqrt((r - da - db)**2 + add) + 1.D00/sqrt((r + da + db)**2&
           + add) + 1.D00/sqrt((r - da + db)**2 + add) + 1.D00/sqrt((r + da - &
          db)**2 + add) - 2.D00/sqrt((r - da)**2 + db**2 + add) - 2.D00/sqrt((r&
           - db)**2 + da**2 + add) - 2.D00/sqrt((r + da)**2 + db**2 + add) - &
          2.D00/sqrt((r + db)**2 + da**2 + add) + 2.D00/sqrt(r**2 + (da - db)**&
          2 + add) + 2.D00/sqrt(r**2 + (da + db)**2 + add)
        xyxy = 4.D00/sqrt(r**2 + (da - db)**2 + add) + 4.D00/sqrt(r**2 + (da + &
          db)**2 + add) - 8.D00/sqrt(r**2 + da**2 + db**2 + add)
        charg = zzzz/16.D00 - xyxy/64.D00
!     X - ZX
      else if (l1==1 .and. l2==2 .and. m==1) then
        ab = db/sqrt(2.D0)
        dxqxz = (-2.D00/sqrt((r - ab)**2 + (da - ab)**2 + add)) + 2.D00/sqrt((r&
           + ab)**2 + (da - ab)**2 + add) + 2.D00/sqrt((r - ab)**2 + (da + ab)&
          **2 + add) - 2.D00/sqrt((r + ab)**2 + (da + ab)**2 + add)
        charg = dxqxz/8.D00
!     ZX - X
      else if (l1==2 .and. l2==1 .and. m==1) then
        aa = da/sqrt(2.D0)
        qxzdx = (-2.D00/sqrt((r + aa)**2 + (aa - db)**2 + add)) + 2.D00/sqrt((r&
           - aa)**2 + (aa - db)**2 + add) + 2.D00/sqrt((r + aa)**2 + (aa + db)&
          **2 + add) - 2.D00/sqrt((r - aa)**2 + (aa + db)**2 + add)
        charg = qxzdx/8.D00
!     ZX - ZX
      else if (l1==2 .and. l2==2 .and. m==1) then
        aa = da/sqrt(2.D0)
        ab = db/sqrt(2.D0)
        qxzqxz = 2.D00/sqrt((r + aa - ab)**2 + (aa - ab)**2 + add) - 2.D00/&
          sqrt((r + aa + ab)**2 + (aa - ab)**2 + add) - 2.D00/sqrt((r - aa - ab&
          )**2 + (aa - ab)**2 + add) + 2.D00/sqrt((r - aa + ab)**2 + (aa - ab)&
          **2 + add) - 2.D00/sqrt((r + aa - ab)**2 + (aa + ab)**2 + add) + &
          2.D00/sqrt((r + aa + ab)**2 + (aa + ab)**2 + add) + 2.D00/sqrt((r - &
          aa - ab)**2 + (aa + ab)**2 + add) - 2.D00/sqrt((r - aa + ab)**2 + (aa&
           + ab)**2 + add)
        charg = qxzqxz/16.D00
!     XX - XX
      else if (l1==2 .and. l2==2 .and. m==2) then
        xyxy = 4.D00/sqrt(r**2 + (da - db)**2 + add) + 4.D00/sqrt(r**2 + (da + &
          db)**2 + add) - 8.D00/sqrt(r**2 + da**2 + db**2 + add)
        charg = xyxy/16.D00
      end if
      return
      end function charg


      subroutine ddpo(ni)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use parameters_C, only : dorbs, gss, hsp, gpp, gp2, ddp, po
      use mndod_C, only :  aij, repd
!     *
!     CALCULATION OF CHARGE SEPARATIONS AND ADDITIVE TERMS USED
!     TO COMPUTE THE TWO-CENTER TWO-ELECTRON INTEGRALS IN MNDO/D.
!     *
!     CHARGE SEPARATIONS DD(6,107) FROM ARRAY AIJ COMPUTED IN AIJM.
!     ADDITIVE TERMS     PO(9,107) FROM FUNCTION POIJ.
!     SECOND INDEX OF DD AND PO    SS 1, SP 2,PP 8, PP 3, SD 4, PD 5, DD
!     SEE EQUATIONS (12)-(16) OF TCA PAPER FOR DD.
!     SEE EQUATIONS (19)-(26) OF TCA PAPER FOR PO.
!     SPECIAL CONVENTION FOR ATOMIC CORE: ADDITIVE TERM PO(9,NI)
!     USED IN THE EVALUATION OF THE CORE-ELECTRON ATTRACTIONS AND
!     CORE-CORE REPULSIONS.
!     *
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: ni
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      double precision :: fg, d, da
      double precision, external :: poij
!-----------------------------------------------
! *** ADDITIVE TERM FOR SS.
      fg = gss(ni)
      if (fg > 0.1d0) po(1,ni) = poij(0,1.d0,fg)
      if (ni >= 3) then
! *** OTHER TERMS FOR SP BASIS.
!     SP
        d = aij(2,ni)/sqrt(12.0D0)
        fg = hsp(ni)
        ddp(2,ni) = d
        po(2,ni) = poij(1,d,fg)
!     PP
        po(7,ni) = po(1,ni)
        d = sqrt(aij(3,ni)*0.1D0)
        fg = 0.5D0*(gpp(ni)-gp2(ni))
        ddp(3,ni) = d
        po(3,ni) = poij(2,d,fg)
! *** TERMS INVOLVING D ORBITALS.
        if (dorbs(ni)) then
!     SD
          da = sqrt(1.d0/60.0D0)
          d = sqrt(aij(4,ni)*da)
          fg = repd(19,ni)
          ddp(4,ni) = d
          po(4,ni) = poij(2,d,fg)
!     PD
          d = aij(5,ni)/sqrt(20.0D0)
          fg = repd(23,ni) - 1.8D0*repd(35,ni)
          ddp(5,ni) = d
          po(5,ni) = poij(1,d,fg)
!     DD
          fg = 0.2D0*(repd(29,ni)+2.d0*repd(30,ni)+2.d0*repd(31,ni))
          if (fg > 1.d-5) then
            po(8,ni) = poij(0,1.d0,fg)
          else
            po(8,ni) = 1.d5
          end if
          d = sqrt(aij(6,ni)/14.0D0)
          fg = repd(44,ni) - (20.0D0/35.0D0)*repd(52,ni)
          ddp(6,ni) = d
          po(6,ni) = poij(2,d,fg)
        end if
      end if
      return
      end subroutine ddpo



      subroutine elenuc(ia, ib, ja, jb, h)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use molkst_C, only : mpack
      use mndod_C, only : sp, sd, pp, dp, d_d, cored, inddd, inddp, indpp
!***********************************************************************
!
!   ELENUC - Nuclear stabilization terms added to one-electron matrix
!
!***********************************************************************
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: ia
      integer , intent(in) :: ib
      integer , intent(in) :: ja
      integer , intent(in) :: jb
      double precision , intent(inout) :: h(mpack)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k, l, n, i, ind1, j, ind2, m, ipp, idp, idd
!-----------------------------------------------
!
!
!     (SS/)=1,   (SO/)=2,   (OO/)=3,   (PP/)=4
!
!     (D S/)=5,  (D P /)=6, (D D /)=7, (D+P+/)=8, (D+D+/)=9, (D#D#/)=10
      k = ia
      l = ib
      n = 1
!
   10 continue
      do i = k, l
        ind1 = i - k
        if (ind1 == 0) then
          do j = k, i
            ind2 = j - k
            m = (i*(i - 1))/2 + j
            if (ind2 == 0) then
! -- (SS/)
              h(m) = h(m) + cored(1,n)
! -- (SD/)
            else
              if (ind2 < 4) then
! -- (PP/)
                ipp = indpp(ind1,ind2)
                h(m) = h(m) + cored(3,n)*pp(ipp,1,1) + cored(4,n)*(pp(ipp,2,2)+pp&
                  (ipp,3,3))
! -- (PD/)
              else
! -- (d_d/)
                idd = inddd(ind1-3,ind2-3)
                h(m) = h(m) + cored(7,n)*d_d(idd,1,1) + cored(9,n)*(d_d(idd,2,2)+d_d&
                  (idd,3,3)) + cored(10,n)*(d_d(idd,4,4)+d_d(idd,5,5))
              end if
            end if
!
          end do
        else
          if (ind1 < 4) then
            do j = k, i
              ind2 = j - k
              m = (i*(i - 1))/2 + j
              if (ind2 == 0) then
! -- (SP/)
                h(m) = h(m) + sp(1,ind1)*cored(2,n)
! -- (SD/)
              else
                if (ind2 < 4) then
! -- (PP/)
                  ipp = indpp(ind1,ind2)
                  h(m) = h(m) + cored(3,n)*pp(ipp,1,1) + cored(4,n)*(pp(ipp,2,2)+&
                    pp(ipp,3,3))
! -- (PD/)
                else
! -- (d_d/)
                  idd = inddd(ind1-3,ind2-3)
                  h(m) = h(m) + cored(7,n)*d_d(idd,1,1) + cored(9,n)*(d_d(idd,2,2)+&
                    d_d(idd,3,3)) + cored(10,n)*(d_d(idd,4,4)+d_d(idd,5,5))
                end if
              end if
!
            end do
          else
            do j = k, i
              ind2 = j - k
              m = (i*(i - 1))/2 + j
              if (ind2 == 0) then
! -- (SD/)
                h(m) = h(m) + sd(1,ind1-3)*cored(5,n)
              else
                if (ind2 < 4) then
! -- (PD/)
                  idp = inddp(ind1-3,ind2)
                  h(m) = h(m) + cored(6,n)*dp(idp,1,1) + cored(8,n)*(dp(idp,2,2)+&
                    dp(idp,3,3))
                else
! -- (d_d/)
                  idd = inddd(ind1-3,ind2-3)
                  h(m) = h(m) + cored(7,n)*d_d(idd,1,1) + cored(9,n)*(d_d(idd,2,2)+&
                    d_d(idd,3,3)) + cored(10,n)*(d_d(idd,4,4)+d_d(idd,5,5))
                end if
              end if
!
            end do
          end if
        end if
      end do
!
      if (n == 2) go to 30
!
      k = ja
      l = jb
      n = 2
!
      go to 10
   30 continue
      return
      end subroutine elenuc


      subroutine fbx
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use mndod_C, only : fx, b
!     *
!     DEFINE FACTORIALS AND BINOMIAL COEFFICIENTS.
!     fx(30)     FACTORIALS.
!     B(30,30)      BINOMIAL COEFFICIENTS.
!     *
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
!-----------------------------------------------
      fx(1) = 1.D0
      do i = 2, 30
        fx(i) = fx(i-1)*dble(i - 1)
      end do
      b(:,1) = 1.D0
      b(:,2:30) = 0.D0
      do i = 2, 30
        b(i,2:i) = b(i-1,:i-1) + b(i-1,2:i)
      end do
      return
      end subroutine fbx

      subroutine fordd
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE mndod_C, only : indx, indexd, ch, ind2, isym, inddd, inddp, indpp
!     *
!     DEFINITION OF INDICES AND LOGICAL VARIABLES.
!     *
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, l
!-----------------------------------------------
!     *
!     INDEX(I,J) AND INDX(I,J) DEFINE ADDRESSES IN LOWER TRIANGLE
!     BY DIFFERENT CONVENTIONS.
!     INDEX(I,J) = 1,2,3 ..  FOR (I,J)=(1,1),(2,1),(3,1),..,(2,2),(3,2)
!     INDX(I,J)  = 1,2,3 ..  FOR (I,J)=(1,1),(2,1),(2,2),(3,1),(3,2) ..
!     INDX(I,J)    CORRESPONDS TO THE USUAL PAIR INDEX.
!     INDEX(J,I) = INDEX(I,J)
!     INDX(J,I)  = INDX(I,J)
!     *
!     *
!     COEFFICIENTS RELATING ANALYTICAL AND POINT-CHARGE MULTIPOLE
!     MOMENTS, SEE EQUATION (17) AND TABLE 2 OF TCA PAPER.
!     *
!     FIRST  INDEX          STANDARD PAIR INDEX (SPD BASIS)
!     SECOND INDEX          L QUANTUM NUMBER FOR MULTIPOLE MOMENT
!     THIRD  INDEX          M QUANTUM NUMBER FOR MULTIPOLE MOMENT
!     *
!
      do i = 1, 9
        do j = 1, i
          indexd(i,j) = (-(j*(j - 1))/2) + i + 9*(j - 1)
          indx(i,j) = (i*(i - 1))/2 + j
          indexd(j,i) = indexd(i,j)
          indx(j,i) = indx(i,j)
        end do
      end do
!
      ind2 = 0
!
!   SP-SP
      ind2(1,1) = 1
      ind2(1,2) = 2
      ind2(1,10) = 3
      ind2(1,18) = 4
      ind2(1,25) = 5
      ind2(2,1) = 6
      ind2(2,2) = 7
      ind2(2,10) = 8
      ind2(2,18) = 9
      ind2(2,25) = 10
      ind2(10,1) = 11
      ind2(10,2) = 12
      ind2(10,10) = 13
      ind2(10,18) = 14
      ind2(10,25) = 15
      ind2(3,3) = 16
      ind2(3,11) = 17
      ind2(11,3) = 18
      ind2(11,11) = 19
      ind2(18,1) = 20
      ind2(18,2) = 21
      ind2(18,10) = 22
      ind2(18,18) = 23
      ind2(18,25) = 24
      ind2(4,4) = 25
      ind2(4,12) = 26
      ind2(12,4) = 27
      ind2(12,12) = 28
      ind2(19,19) = 29
      ind2(25,1) = 30
      ind2(25,2) = 31
      ind2(25,10) = 32
      ind2(25,18) = 33
      ind2(25,25) = 34
!   SPD-SPD
      ind2(1,5) = 35
      ind2(1,13) = 36
      ind2(1,31) = 37
      ind2(1,21) = 38
      ind2(1,36) = 39
      ind2(1,28) = 40
      ind2(1,40) = 41
      ind2(1,43) = 42
      ind2(1,45) = 43
      ind2(2,5) = 44
      ind2(2,13) = 45
      ind2(2,31) = 46
      ind2(2,21) = 47
      ind2(2,36) = 48
      ind2(2,28) = 49
      ind2(2,40) = 50
      ind2(2,43) = 51
      ind2(2,45) = 52
      ind2(10,5) = 53
      ind2(10,13) = 54
      ind2(10,31) = 55
      ind2(10,21) = 56
      ind2(10,36) = 57
      ind2(10,28) = 58
      ind2(10,40) = 59
      ind2(10,43) = 60
      ind2(10,45) = 61
      ind2(3,20) = 62
      ind2(3,6) = 63
      ind2(3,14) = 64
      ind2(3,32) = 65
      ind2(3,23) = 66
      ind2(3,38) = 67
      ind2(3,30) = 68
      ind2(3,42) = 69
      ind2(11,20) = 70
      ind2(11,6) = 71
      ind2(11,14) = 72
      ind2(11,32) = 73
      ind2(11,23) = 74
      ind2(11,38) = 75
      ind2(11,30) = 76
      ind2(11,42) = 77
      ind2(18,5) = 78
      ind2(18,13) = 79
      ind2(18,31) = 80
      ind2(18,21) = 81
      ind2(18,36) = 82
      ind2(18,28) = 83
      ind2(18,40) = 84
      ind2(18,8) = 85
      ind2(18,16) = 86
      ind2(18,34) = 87
      ind2(18,43) = 88
      ind2(18,45) = 89
      ind2(4,26) = 90
      ind2(4,7) = 91
      ind2(4,15) = 92
      ind2(4,33) = 93
      ind2(4,29) = 94
      ind2(4,41) = 95
      ind2(4,24) = 96
      ind2(4,39) = 97
      ind2(12,26) = 98
      ind2(12,7) = 99
      ind2(12,15) = 100
      ind2(12,33) = 101
      ind2(12,29) = 102
      ind2(12,41) = 103
      ind2(12,24) = 104
      ind2(12,39) = 105
      ind2(19,27) = 106
      ind2(19,22) = 107
      ind2(19,37) = 108
      ind2(19,9) = 109
      ind2(19,17) = 110
      ind2(19,35) = 111
      ind2(25,5) = 112
      ind2(25,13) = 113
      ind2(25,31) = 114
      ind2(25,21) = 115
      ind2(25,36) = 116
      ind2(25,28) = 117
      ind2(25,40) = 118
      ind2(25,8) = 119
      ind2(25,16) = 120
      ind2(25,34) = 121
      ind2(25,43) = 122
      ind2(25,45) = 123
      ind2(5,1) = 124
      ind2(5,2) = 125
      ind2(5,10) = 126
      ind2(5,18) = 127
      ind2(5,25) = 128
      ind2(5,5) = 129
      ind2(5,13) = 130
      ind2(5,31) = 131
      ind2(5,21) = 132
      ind2(5,36) = 133
      ind2(5,28) = 134
      ind2(5,40) = 135
      ind2(5,43) = 136
      ind2(5,45) = 137
      ind2(13,1) = 138
      ind2(13,2) = 139
      ind2(13,10) = 140
      ind2(13,18) = 141
      ind2(13,25) = 142
      ind2(13,5) = 143
      ind2(13,13) = 144
      ind2(13,31) = 145
      ind2(13,21) = 146
      ind2(13,36) = 147
      ind2(13,28) = 148
      ind2(13,40) = 149
      ind2(13,43) = 150
      ind2(13,45) = 151
      ind2(20,3) = 152
      ind2(20,11) = 153
      ind2(20,20) = 154
      ind2(20,6) = 155
      ind2(20,14) = 156
      ind2(20,32) = 157
      ind2(20,23) = 158
      ind2(20,38) = 159
      ind2(20,30) = 160
      ind2(20,42) = 161
      ind2(26,4) = 162
      ind2(26,12) = 163
      ind2(26,26) = 164
      ind2(26,7) = 165
      ind2(26,15) = 166
      ind2(26,33) = 167
      ind2(26,29) = 168
      ind2(26,41) = 169
      ind2(26,24) = 170
      ind2(26,39) = 171
      ind2(31,1) = 172
      ind2(31,2) = 173
      ind2(31,10) = 174
      ind2(31,18) = 175
      ind2(31,25) = 176
      ind2(31,5) = 177
      ind2(31,13) = 178
      ind2(31,31) = 179
      ind2(31,21) = 180
      ind2(31,36) = 181
      ind2(31,28) = 182
      ind2(31,40) = 183
      ind2(31,43) = 184
      ind2(31,45) = 185
      ind2(6,3) = 186
      ind2(6,11) = 187
      ind2(6,20) = 188
      ind2(6,6) = 189
      ind2(6,14) = 190
      ind2(6,32) = 191
      ind2(6,23) = 192
      ind2(6,38) = 193
      ind2(6,30) = 194
      ind2(6,42) = 195
      ind2(14,3) = 196
      ind2(14,11) = 197
      ind2(14,20) = 198
      ind2(14,6) = 199
      ind2(14,14) = 200
      ind2(14,32) = 201
      ind2(14,23) = 202
      ind2(14,38) = 203
      ind2(14,30) = 204
      ind2(14,42) = 205
      ind2(21,1) = 206
      ind2(21,2) = 207
      ind2(21,10) = 208
      ind2(21,18) = 209
      ind2(21,25) = 210
      ind2(21,5) = 211
      ind2(21,13) = 212
      ind2(21,31) = 213
      ind2(21,21) = 214
      ind2(21,36) = 215
      ind2(21,28) = 216
      ind2(21,40) = 217
      ind2(21,8) = 218
      ind2(21,16) = 219
      ind2(21,34) = 220
      ind2(21,43) = 221
      ind2(21,45) = 222
      ind2(27,19) = 223
      ind2(27,27) = 224
      ind2(27,22) = 225
      ind2(27,37) = 226
      ind2(27,9) = 227
      ind2(27,17) = 228
      ind2(27,35) = 229
      ind2(32,3) = 230
      ind2(32,11) = 231
      ind2(32,20) = 232
      ind2(32,6) = 233
      ind2(32,14) = 234
      ind2(32,32) = 235
      ind2(32,23) = 236
      ind2(32,38) = 237
      ind2(32,30) = 238
      ind2(32,42) = 239
      ind2(36,1) = 240
      ind2(36,2) = 241
      ind2(36,10) = 242
      ind2(36,18) = 243
      ind2(36,25) = 244
      ind2(36,5) = 245
      ind2(36,13) = 246
      ind2(36,31) = 247
      ind2(36,21) = 248
      ind2(36,36) = 249
      ind2(36,28) = 250
      ind2(36,40) = 251
      ind2(36,8) = 252
      ind2(36,16) = 253
      ind2(36,34) = 254
      ind2(36,43) = 255
      ind2(36,45) = 256
      ind2(7,4) = 257
      ind2(7,12) = 258
      ind2(7,26) = 259
      ind2(7,7) = 260
      ind2(7,15) = 261
      ind2(7,33) = 262
      ind2(7,29) = 263
      ind2(7,41) = 264
      ind2(7,24) = 265
      ind2(7,39) = 266
      ind2(15,4) = 267
      ind2(15,12) = 268
      ind2(15,26) = 269
      ind2(15,7) = 270
      ind2(15,15) = 271
      ind2(15,33) = 272
      ind2(15,29) = 273
      ind2(15,41) = 274
      ind2(15,24) = 275
      ind2(15,39) = 276
      ind2(22,19) = 277
      ind2(22,27) = 278
      ind2(22,22) = 279
      ind2(22,37) = 280
      ind2(22,9) = 281
      ind2(22,17) = 282
      ind2(22,35) = 283
      ind2(28,1) = 284
      ind2(28,2) = 285
      ind2(28,10) = 286
      ind2(28,18) = 287
      ind2(28,25) = 288
      ind2(28,5) = 289
      ind2(28,13) = 290
      ind2(28,31) = 291
      ind2(28,21) = 292
      ind2(28,36) = 293
      ind2(28,28) = 294
      ind2(28,40) = 295
      ind2(28,8) = 296
      ind2(28,16) = 297
      ind2(28,34) = 298
      ind2(28,43) = 299
      ind2(28,45) = 300
      ind2(33,4) = 301
      ind2(33,12) = 302
      ind2(33,26) = 303
      ind2(33,7) = 304
      ind2(33,15) = 305
      ind2(33,33) = 306
      ind2(33,29) = 307
      ind2(33,41) = 308
      ind2(33,24) = 309
      ind2(33,39) = 310
      ind2(37,19) = 311
      ind2(37,27) = 312
      ind2(37,22) = 313
      ind2(37,37) = 314
      ind2(37,9) = 315
      ind2(37,17) = 316
      ind2(37,35) = 317
      ind2(40,1) = 318
      ind2(40,2) = 319
      ind2(40,10) = 320
      ind2(40,18) = 321
      ind2(40,25) = 322
      ind2(40,5) = 323
      ind2(40,13) = 324
      ind2(40,31) = 325
      ind2(40,21) = 326
      ind2(40,36) = 327
      ind2(40,28) = 328
      ind2(40,40) = 329
      ind2(40,8) = 330
      ind2(40,16) = 331
      ind2(40,34) = 332
      ind2(40,43) = 333
      ind2(40,45) = 334
      ind2(8,18) = 335
      ind2(8,25) = 336
      ind2(8,21) = 337
      ind2(8,36) = 338
      ind2(8,28) = 339
      ind2(8,40) = 340
      ind2(8,8) = 341
      ind2(8,16) = 342
      ind2(8,34) = 343
      ind2(16,18) = 344
      ind2(16,25) = 345
      ind2(16,21) = 346
      ind2(16,36) = 347
      ind2(16,28) = 348
      ind2(16,40) = 349
      ind2(16,8) = 350
      ind2(16,16) = 351
      ind2(16,34) = 352
      ind2(23,3) = 353
      ind2(23,11) = 354
      ind2(23,20) = 355
      ind2(23,6) = 356
      ind2(23,14) = 357
      ind2(23,32) = 358
      ind2(23,23) = 359
      ind2(23,38) = 360
      ind2(23,30) = 361
      ind2(23,42) = 362
      ind2(29,4) = 363
      ind2(29,12) = 364
      ind2(29,26) = 365
      ind2(29,7) = 366
      ind2(29,15) = 367
      ind2(29,33) = 368
      ind2(29,29) = 369
      ind2(29,41) = 370
      ind2(29,24) = 371
      ind2(29,39) = 372
      ind2(34,18) = 373
      ind2(34,25) = 374
      ind2(34,21) = 375
      ind2(34,36) = 376
      ind2(34,28) = 377
      ind2(34,40) = 378
      ind2(34,8) = 379
      ind2(34,16) = 380
      ind2(34,34) = 381
      ind2(38,3) = 382
      ind2(38,11) = 383
      ind2(38,20) = 384
      ind2(38,6) = 385
      ind2(38,14) = 386
      ind2(38,32) = 387
      ind2(38,23) = 388
      ind2(38,38) = 389
      ind2(38,30) = 390
      ind2(38,42) = 391
      ind2(41,4) = 392
      ind2(41,12) = 393
      ind2(41,26) = 394
      ind2(41,7) = 395
      ind2(41,15) = 396
      ind2(41,33) = 397
      ind2(41,29) = 398
      ind2(41,41) = 399
      ind2(41,24) = 400
      ind2(41,39) = 401
      ind2(43,1) = 402
      ind2(43,2) = 403
      ind2(43,10) = 404
      ind2(43,18) = 405
      ind2(43,25) = 406
      ind2(43,5) = 407
      ind2(43,13) = 408
      ind2(43,31) = 409
      ind2(43,21) = 410
      ind2(43,36) = 411
      ind2(43,28) = 412
      ind2(43,40) = 413
      ind2(43,43) = 414
      ind2(43,45) = 415
      ind2(9,19) = 416
      ind2(9,27) = 417
      ind2(9,22) = 418
      ind2(9,37) = 419
      ind2(9,9) = 420
      ind2(9,17) = 421
      ind2(9,35) = 422
      ind2(17,19) = 423
      ind2(17,27) = 424
      ind2(17,22) = 425
      ind2(17,37) = 426
      ind2(17,9) = 427
      ind2(17,17) = 428
      ind2(17,35) = 429
      ind2(24,4) = 430
      ind2(24,12) = 431
      ind2(24,26) = 432
      ind2(24,7) = 433
      ind2(24,15) = 434
      ind2(24,33) = 435
      ind2(24,29) = 436
      ind2(24,41) = 437
      ind2(24,24) = 438
      ind2(24,39) = 439
      ind2(30,3) = 440
      ind2(30,11) = 441
      ind2(30,20) = 442
      ind2(30,6) = 443
      ind2(30,14) = 444
      ind2(30,32) = 445
      ind2(30,23) = 446
      ind2(30,38) = 447
      ind2(30,30) = 448
      ind2(30,42) = 449
      ind2(35,19) = 450
      ind2(35,27) = 451
      ind2(35,22) = 452
      ind2(35,37) = 453
      ind2(35,9) = 454
      ind2(35,17) = 455
      ind2(35,35) = 456
      ind2(39,4) = 457
      ind2(39,12) = 458
      ind2(39,26) = 459
      ind2(39,7) = 460
      ind2(39,15) = 461
      ind2(39,33) = 462
      ind2(39,29) = 463
      ind2(39,41) = 464
      ind2(39,24) = 465
      ind2(39,39) = 466
      ind2(42,3) = 467
      ind2(42,11) = 468
      ind2(42,20) = 469
      ind2(42,6) = 470
      ind2(42,14) = 471
      ind2(42,32) = 472
      ind2(42,23) = 473
      ind2(42,38) = 474
      ind2(42,30) = 475
      ind2(42,42) = 476
      ind2(44,44) = 477
      ind2(45,1) = 478
      ind2(45,2) = 479
      ind2(45,10) = 480
      ind2(45,18) = 481
      ind2(45,25) = 482
      ind2(45,5) = 483
      ind2(45,13) = 484
      ind2(45,31) = 485
      ind2(45,21) = 486
      ind2(45,36) = 487
      ind2(45,28) = 488
      ind2(45,40) = 489
      ind2(45,43) = 490
      ind2(45,45) = 491
      isym = 0
!
      isym(40) = 38
      isym(41) = 39
      isym(43) = 42
      isym(49) = 47
      isym(50) = 48
      isym(52) = 51
      isym(58) = 56
      isym(59) = 57
      isym(61) = 60
      isym(68) = 66
      isym(69) = 67
      isym(76) = 74
      isym(77) = 75
      isym(89) = 88
      isym(90) = 62
      isym(91) = 63
      isym(92) = 64
      isym(93) = 65
      isym(94) = -66
      isym(95) = -67
      isym(96) = 66
      isym(97) = 67
      isym(98) = 70
      isym(99) = 71
      isym(100) = 72
      isym(101) = 73
      isym(102) = -74
      isym(103) = -75
      isym(104) = 74
      isym(105) = 75
      isym(106) = 86
      isym(107) = 86
      isym(109) = 85
      isym(110) = 86
      isym(111) = 87
      isym(112) = 78
      isym(113) = 79
      isym(114) = 80
      isym(115) = 83
      isym(116) = 84
      isym(117) = 81
      isym(118) = 82
      isym(119) = -85
      isym(120) = -86
      isym(121) = -87
      isym(122) = 88
      isym(123) = 88
      isym(128) = 127
      isym(134) = 132
      isym(135) = 133
      isym(137) = 136
      isym(142) = 141
      isym(148) = 146
      isym(149) = 147
      isym(151) = 150
      isym(160) = 158
      isym(161) = 159
      isym(162) = 152
      isym(163) = 153
      isym(164) = 154
      isym(165) = 155
      isym(166) = 156
      isym(167) = 157
      isym(168) = -158
      isym(169) = -159
      isym(170) = 158
      isym(171) = 159
      isym(176) = 175
      isym(182) = 180
      isym(183) = 181
      isym(185) = 184
      isym(194) = 192
      isym(195) = 193
      isym(204) = 202
      isym(205) = 203
      isym(222) = 221
      isym(224) = 219
      isym(225) = 219
      isym(227) = 218
      isym(228) = 219
      isym(229) = 220
      isym(238) = 236
      isym(239) = 237
      isym(256) = 255
      isym(257) = 186
      isym(258) = 187
      isym(259) = 188
      isym(260) = 189
      isym(261) = 190
      isym(262) = 191
      isym(263) = -192
      isym(264) = -193
      isym(265) = 192
      isym(266) = 193
      isym(267) = 196
      isym(268) = 197
      isym(269) = 198
      isym(270) = 199
      isym(271) = 200
      isym(272) = 201
      isym(273) = -202
      isym(274) = -203
      isym(275) = 202
      isym(276) = 203
      isym(277) = 223
      isym(278) = 219
      isym(279) = 219
      isym(280) = 226
      isym(281) = 218
      isym(282) = 219
      isym(283) = 220
      isym(284) = 206
      isym(285) = 207
      isym(286) = 208
      isym(287) = 210
      isym(288) = 209
      isym(289) = 211
      isym(290) = 212
      isym(291) = 213
      isym(292) = 216
      isym(293) = 217
      isym(294) = 214
      isym(295) = 215
      isym(296) = -218
      isym(297) = -219
      isym(298) = -220
      isym(299) = 221
      isym(300) = 221
      isym(301) = 230
      isym(302) = 231
      isym(303) = 232
      isym(304) = 233
      isym(305) = 234
      isym(306) = 235
      isym(307) = -236
      isym(308) = -237
      isym(309) = 236
      isym(310) = 237
      isym(312) = 253
      isym(313) = 253
      isym(315) = 252
      isym(316) = 253
      isym(317) = 254
      isym(318) = 240
      isym(319) = 241
      isym(320) = 242
      isym(321) = 244
      isym(322) = 243
      isym(323) = 245
      isym(324) = 246
      isym(325) = 247
      isym(326) = 250
      isym(327) = 251
      isym(328) = 248
      isym(329) = 249
      isym(330) = -252
      isym(331) = -253
      isym(332) = -254
      isym(333) = 255
      isym(334) = 255
      isym(336) = -335
      isym(339) = -337
      isym(340) = -338
      isym(342) = 337
      isym(344) = 223
      isym(345) = -223
      isym(346) = 219
      isym(347) = 226
      isym(348) = -219
      isym(349) = -226
      isym(350) = 218
      isym(351) = 219
      isym(352) = 220
      isym(363) = -353
      isym(364) = -354
      isym(365) = -355
      isym(366) = -356
      isym(367) = -357
      isym(368) = -358
      isym(369) = 359
      isym(370) = 360
      isym(371) = -361
      isym(372) = -362
      isym(374) = -373
      isym(377) = -375
      isym(378) = -376
      isym(380) = 375
      isym(392) = -382
      isym(393) = -383
      isym(394) = -384
      isym(395) = -385
      isym(396) = -386
      isym(397) = -387
      isym(398) = 388
      isym(399) = 389
      isym(400) = -390
      isym(401) = -391
      isym(406) = 405
      isym(412) = 410
      isym(413) = 411
      isym(416) = 335
      isym(417) = 337
      isym(418) = 337
      isym(419) = 338
      isym(420) = 341
      isym(421) = 337
      isym(422) = 343
      isym(423) = 223
      isym(424) = 219
      isym(425) = 219
      isym(426) = 226
      isym(427) = 218
      isym(428) = 219
      isym(429) = 220
      isym(430) = 353
      isym(431) = 354
      isym(432) = 355
      isym(433) = 356
      isym(434) = 357
      isym(435) = 358
      isym(436) = -361
      isym(437) = -362
      isym(438) = 359
      isym(439) = 360
      isym(440) = 353
      isym(441) = 354
      isym(442) = 355
      isym(443) = 356
      isym(444) = 357
      isym(445) = 358
      isym(446) = 361
      isym(447) = 362
      isym(448) = 359
      isym(449) = 360
      isym(450) = 373
      isym(451) = 375
      isym(452) = 375
      isym(453) = 376
      isym(454) = 379
      isym(455) = 375
      isym(456) = 381
      isym(457) = 382
      isym(458) = 383
      isym(459) = 384
      isym(460) = 385
      isym(461) = 386
      isym(462) = 387
      isym(463) = -390
      isym(464) = -391
      isym(465) = 388
      isym(466) = 389
      isym(467) = 382
      isym(468) = 383
      isym(469) = 384
      isym(470) = 385
      isym(471) = 386
      isym(472) = 387
      isym(473) = 390
      isym(474) = 391
      isym(475) = 388
      isym(476) = 389
      isym(478) = 402
      isym(479) = 403
      isym(480) = 404
      isym(481) = 405
      isym(482) = 405
      isym(483) = 407
      isym(484) = 408
      isym(485) = 409
      isym(486) = 410
      isym(487) = 411
      isym(488) = 410
      isym(489) = 411
      isym(490) = 415
      isym(491) = 414
!   *
      do i = 1, 45
        do l = 0, 2
          ch(i,l,(-l):l) = 0.d0
        end do
      end do
!   *
! *** THE STANDARD MNDO93 CODE DEFINES THE FOLLOWING CONSTANTS
! *** MORE PRECISELY (BUT WITH A DIFFERENT NUMBERING SCHEME).
!     PARAMETER (CLM3  = 0.13333333333333D+01)
!     PARAMETER (CLM6  =-0.66666666666667D+00)
!     PARAMETER (CLM10 =-0.66666666666667D+00)
!     PARAMETER (CLM11 = 0.11547005383793D+01)
!     PARAMETER (CLM12 = 0.11547005383793D+01)
!     PARAMETER (CLM13 =-0.57735026918963D+00)
!     PARAMETER (CLM15 = 0.13333333333333D+01)
!     PARAMETER (CLM20 = 0.57735026918963D+00)
!     PARAMETER (CLM21 = 0.66666666666667D+00)
!     PARAMETER (CLM28 = 0.66666666666667D+00)
!     PARAMETER (CLM33 =-0.11547005383793D+01)
!     PARAMETER (CLM36 =-0.13333333333333D+01)
!     PARAMETER (CLM45 =-0.13333333333333D+01)
! *** PLEASE MAKE THE OBVIOUS CORRECTIONS.
      ch(1,0,0) = 1.d0
      ch(2,1,0) = 1.d0
      ch(3,1,1) = 1.d0
      ch(4,1,-1) = 1.d0
      ch(5,2,0) = 1.15470054D0
      ch(6,2,1) = 1.d0
      ch(7,2,-1) = 1.d0
      ch(8,2,2) = 1.d0
      ch(9,2,-2) = 1.d0
      ch(10,0,0) = 1.d0
      ch(10,2,0) = 1.33333333D0
      ch(11,2,1) = 1.d0
      ch(12,2,-1) = 1.d0
      ch(13,1,0) = 1.15470054D0
      ch(14,1,1) = 1.d0
      ch(15,1,-1) = 1.d0
      ch(18,0,0) = 1.d0
      ch(18,2,0) = -.66666667D0
      ch(18,2,2) = 1.d0
      ch(19,2,-2) = 1.d0
      ch(20,1,1) = -.57735027D0
      ch(21,1,0) = 1.d0
      ch(23,1,1) = 1.d0
      ch(24,1,-1) = 1.d0
      ch(25,0,0) = 1.d0
      ch(25,2,0) = -.66666667D0
      ch(25,2,2) = -1.d0
      ch(26,1,-1) = -.57735027D0
      ch(28,1,0) = 1.d0
      ch(29,1,-1) = -1.d0
      ch(30,1,1) = 1.d0
      ch(31,0,0) = 1.d0
      ch(31,2,0) = 1.33333333D0
      ch(32,2,1) = .57735027D0
      ch(33,2,-1) = .57735027D0
      ch(34,2,2) = -1.15470054D0
      ch(35,2,-2) = -1.15470054D0
      ch(36,0,0) = 1.d0
      ch(36,2,0) = .66666667D0
      ch(36,2,2) = 1.d0
      ch(37,2,-2) = 1.d0
      ch(38,2,1) = 1.d0
      ch(39,2,-1) = 1.d0
      ch(40,0,0) = 1.d0
      ch(40,2,0) = .66666667D0
      ch(40,2,2) = -1.d0
      ch(41,2,-1) = -1.d0
      ch(42,2,1) = 1.d0
      ch(43,0,0) = 1.d0
      ch(43,2,0) = -1.33333333D0
      ch(45,0,0) = 1.d0
      ch(45,2,0) = -1.33333333D0
!
!   INDPP
      indpp(1,1) = 1
      indpp(2,1) = 4
      indpp(3,1) = 5
      indpp(1,2) = 4
      indpp(2,2) = 2
      indpp(3,2) = 6
      indpp(1,3) = 5
      indpp(2,3) = 6
      indpp(3,3) = 3
!   INDDP
      inddp(1,1) = 1
      inddp(2,1) = 4
      inddp(3,1) = 7
      inddp(4,1) = 10
      inddp(5,1) = 13
      inddp(1,2) = 2
      inddp(2,2) = 5
      inddp(3,2) = 8
      inddp(4,2) = 11
      inddp(5,2) = 14
      inddp(1,3) = 3
      inddp(2,3) = 6
      inddp(3,3) = 9
      inddp(4,3) = 12
      inddp(5,3) = 15
!   INDDD
      inddd(1,1) = 1
      inddd(2,1) = 6
      inddd(3,1) = 7
      inddd(4,1) = 9
      inddd(5,1) = 12
      inddd(1,2) = 6
      inddd(2,2) = 2
      inddd(3,2) = 8
      inddd(4,2) = 10
      inddd(5,2) = 13
      inddd(1,3) = 7
      inddd(2,3) = 8
      inddd(3,3) = 3
      inddd(4,3) = 11
      inddd(5,3) = 14
      inddd(1,4) = 9
      inddd(2,4) = 10
      inddd(3,4) = 11
      inddd(4,4) = 4
      inddd(5,4) = 15
      inddd(1,5) = 12
      inddd(2,5) = 13
      inddd(3,5) = 14
      inddd(4,5) = 15
      inddd(5,5) = 5
      return
      end subroutine fordd
      subroutine to_point(r, point, const)
        use molkst_C, only : trunc_1, trunc_2
        use funcon_C, only : ev, a0
        double precision, intent(in) :: r
        double precision, intent(out) :: point, const
!
! smooth the transition from NDDO-type two-electron integrals to exact point-charge
! values.  The function used here is the reverse of the intuitive function:
! it starts off with the Gaussian value being small and builds up to 1.0 at the
! distance "trunc_1".  This behavior is designed to have the smallest effect where
! the difference between point change and NDDO values is largest.
!
! On input, r = interatomic distance in Angstroms
!
! On output, point = e-e repulsion point-charge energy in eV
!            const = fraction of NDDO term to be used
!
        point = ev*a0/r
        if (r < trunc_1) then
          const = 1.d0 - exp(-(r - trunc_1)**2*trunc_2)
        else
          const = 0.d0
        end if
        return
      end subroutine to_point
