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

      subroutine thermo(a, b, c, linear, sym, vibs, nvibs, escf)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use funcon_C, only : fpc_5, fpc_6, fpc_7, fpc_8, fpc_10, pi
      use molkst_C, only : keywrd, title, koment, mol_weight, &
      cp298 => temp_1, s298 => temp_2, ilim
      use common_arrays_C, only : T_range, HOF_tot, H_tot, Cp_tot, S_tot
      use chanel_C, only : iw
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n T s
!-----------------------------------------------
      integer , intent(in) :: nvibs
      double precision , intent(in) :: a
      double precision , intent(in) :: b
      double precision , intent(in) :: c
      double precision , intent(in) :: sym
      double precision , intent(in) :: escf
      logical , intent(in) :: linear
      double precision , intent(inout) :: vibs(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: it1, it2, istep, i, itemp, ir
      double precision :: R, h, ak, ac, qtr2, T, c1, qv, hv, cpv, sv1, wi, ewj, &
        ewjr, e0, sv, qr, cpr, hr, sr, qint, hint, cpint, sint, qtr, cptr, htr&
        , str, cptot, stot, htot, h298
      character :: tmpkey*241
      double precision, external :: reada
!-----------------------------------------------
      if (allocated(T_range)) deallocate(T_range, HOF_tot, H_tot, Cp_tot, S_tot)
      allocate(T_range(300), HOF_tot(300), H_tot(300), Cp_tot(300), S_tot(300))
      R = fpc_5
      h = fpc_6
      ak = fpc_7
      ac = fpc_8
      it1 = 200
      it2 = 400
      istep = 10
      tmpkey = trim(keywrd)
      i = index(tmpkey,'THERMO(') + index(tmpkey,'THERMO=(')
      if (i /= 0) then
!
!   ERASE ALL TEXT FROM TMPKEY EXCEPT THERMO DATA
!
        tmpkey = tmpkey(i + index(tmpkey(i:), "("):)
        tmpkey(index(tmpkey,')'):) = ' '
        it1 = nint(reada(tmpkey,1))
        if (it1 < 100) then
          write (iw, &
      '(2/10X,''TEMPERATURE RANGE STARTS TOO LOW, LOWER BOUND IS RESET TO 100K'')')
          it1 = 100
        end if
        istep = 0
        do i = 1, len_trim(tmpkey)
          if (tmpkey(i:i) == ",") istep = istep + 1
        end do
        select case (istep)
          case (0)
!
!  Only starting temperatures provided
!
            it2 = it1 + 200
            istep = 10
          case (1)
!
!  Starting and ending temperatures provided
!
            i = index(tmpkey,',')
            it2 = nint(reada(tmpkey,i))
            istep = (it2 - it1)/20
            if (istep == 0) then
              istep = 1
            else if (istep <   5) then
              istep = 2
            else if (istep <  10) then
              istep = 5
            else if (istep <  20) then
              istep = 10
            else if (istep <  50) then
              istep = 20
            else if (istep < 100) then
              istep = 50
            else
              istep = 100
            end if
          case (2)
!
!  Starting, step-size, and ending temperatures provided
!
          i = index(tmpkey,',')
          tmpkey(i:i) = ' '
          it2 = nint(reada(tmpkey,i))
          istep = max0(1,istep)
          i = index(tmpkey,',')
          istep = nint(reada(tmpkey,i))
        end select
      end if
      if (nvibs > 0) then
        write (iw, '(2/,A)') trim(title)
        write (iw, '(A)') trim(koment)
        if (linear) then
          write (iw, '(2/10X,''MOLECULE IS LINEAR'')')
        else
          write (iw, '(2/10X,''MOLECULE IS NOT LINEAR'')')
        end if
        if (nvibs > 99) then
          write (iw, &
          '(/10X,''THERE ARE'',I5,'' GENUINE VIBRATIONS IN THIS '',''SYSTEM'')') &
          nvibs
        else
          write (iw, &
          '(/10X,''THERE ARE'',I3,'' GENUINE VIBRATIONS IN THIS '',''SYSTEM'')') &
          nvibs
        end if
        write (iw, 20)
  20    format(10x,'THIS THERMODYNAMICS CALCULATION IS LIMITED TO',/,10x,&
          'MOLECULES WHICH HAVE NO INTERNAL ROTATIONS'/,/)
      end if
      write (iw, '(2/20X,''CALCULATED THERMODYNAMIC PROPERTIES'')')
      write (iw, '(42X,''*'')')
      write (iw, &
      '(''   TEMP. (K)   PARTITION FUNCTION   H.O.F.    ENTHALPY  HEAT CAPACITY ENTROPY'')')
      write (iw, &
      '(''                                   KCAL/MOL   CAL/MOLE    CAL/K/MOL  CAL/K/MOL'',/)')
      do i = 1, nvibs
        vibs(i) = abs(vibs(i))
      end do
      ilim = 1
      if (istep > it2) then
        i = istep
        istep = it2
        it2 = i
      end if
      do itemp = it1, it2, istep
        ilim = ilim + 1
        T_range(ilim) = itemp
      end do
      T_range(1) = 298.D0
      qtr2 = 2.D0*pi*ac*mol_weight/fpc_10
      if (ilim == 2 .and. Abs(T_range(1) - T_range(2)) < 1.d-4) ilim = 1
      do ir = 1, ilim
        itemp = int(T_range(ir))
        T = itemp
!   ***   INITIALISE SOME VARIABLES   ***
        c1 = h*ac/ak/T
        qv = 1.0D0
        hv = 0.0D0
        cpv = 0.0D0
        sv1 = 0.0D0
!   ***   CONSTRUCT THE FREQUENCY DEPENDENT PARTS OF PARTITION FUNCTION
        do i = 1, nvibs
          wi = vibs(i)
!     INSERTED BY PROF. HIRANO
!        To exclude imaginary Wi (here expressed in negative value) and
!        small Wi arising from the imperfect geometry optimization,
!        Wi less than 100 cm-1 (an arbitrary threshold) are neglected
!        in the calculations of vibrational contributions.
          if (wi < 1.d-4) cycle
!
          ewj = dexp((-wi*c1))
          ewjr = 1.D0 - ewj
          qv = qv/ewjr
          e0 = wi*ewj/ewjr
          hv = hv + e0
          cpv = cpv + e0*e0/ewj
          sv1 = sv1 + dlog(ewjr)
        end do
!   ***   FINISH CALCULATION OF VIBRATIONAL PARTS   ***
        e0 = R*c1
        sv = hv*e0 - R*sv1
        hv = hv*e0*T
        cpv = cpv*e0*c1
        if (nvibs == 0) then
          qr = 0.d0 !  This is an atom, therefore no vibrations
          cpr = 0.d0
          hr = cpr*T
          sr = 0.d0
!   ***   NOW CALCULATE THE ROTATIONAL PARTS  (FIRST NON-LINEAR MOLECULES)
        else if (.not.linear) then
          e0 = pi/(a*b*c)
          qr = dsqrt(e0/c1)/c1/sym
          cpr = 1.5D0*R
          hr = cpr*T
          sr = R*(1.5D0*dlog(1.D0/c1) - dlog(sym) + dlog(e0)/2.D0 + 1.5D0)
        else
          qr = 1.D0/(c1*a*sym)
          cpr = R
          hr = cpr*T
          sr = R*dlog(qr) + R
        end if
!   ***   CALCULATE INTERNAL CONTRIBUTIONS   ***
        qint = qv*qr
        hint = hv + hr
        cpint = cpv + cpr
        sint = sv + sr
!   ***   CONSTRUCT TRANSLATION CONTRIBUTIONS   ***
        qtr = dsqrt(qtr2/c1/h)**3
        cptr = 2.5D0*R
        htr = cptr*T
!     UPDATED BY PROF.HIRANO AND MODIFIED BY YI
        str = 4.96804D0*(log(T) + 0.6D0*log(mol_weight)) - 2.31482D0
!
!  Which is best?
!
!   STR=2.2868D0*(5.D0*LOG10(T)+3.D0*LOG10(mol_weight))-2.3135D0
!   STR=0.993608D0*(5.D0*LOG(T)+3.D0*LOG(mol_weight))-2.31482D0
!   STR=4.96804D0*(LOG(T)+0.6D0*LOG(mol_weight))-2.31482D0
!
!   ***   CONSTRUCT TOTALS   ***
        cptot = cptr + cpint
        stot = str + sint
        htot = htr + hint
!   ***   OUTPUT SECTION   ***
!     UPDATED BY PROF. HIRANO
        if (ir == 1) h298 = htot
        write (iw, '(/,F7.2,''  VIB.'',D15.4,10X,f17.4,2F11.4)') T, qv, hv, cpv, sv
        write (iw, '(9X,''ROT.'',D15.4,10X,f17.4,2F11.4)') qr, hr, cpr, sr
        write (iw, '(9X,''INT.'',D15.4,10X,f17.4,2F11.4)') qint, hint, cpint, sint
        write (iw, '(9X,''TRA.'',D15.4,10X,f17.4,2F11.4)') qtr, htr, cptr, str
        write (iw, '(9X,''TOT.'',12X,F17.3,F13.4,2F11.4)') escf + (htot - h298)/1000.D0, htot, cptot, stot
        HOF_tot(ir) = escf + (htot - h298)/1000.D0
        H_tot(ir)   = htot
        Cp_tot(ir)  = cptot
        S_tot(ir)   = stot
      end do
      if (ilim == 1 .and. Abs(T_range(1) - 298) < 1.d1) then
        cp298 = cptot
        s298  = stot
      end if
      write (iw, &
      '(/3X,''*: NOTE: HEATS OF FORMATION ARE RELATIVE TO THE'',/12X,&
        &''ELEMENTS IN THEIR STANDARD STATE AT 298K'')')
!     INSERTED BY PROF.HIRANO
      write (iw, '(12X,''  (=  Standard Enthalpy of Formation).'')')  !,        /12X,&
   !     &''Heat of Formation at T is calculated as [Escf + (HT - H298.15)-.'')')
      write (iw, &
      & '(/12X,''Hvib:  Zero-point energy is not included.'',    /12X,&
      &''       frequencies of less than zero cm-1 are not included.'', &
      & /12X,''Hrot = (3/2)RT''/12X,''Htra = (3/2)RT + pV = (5/2)RT'')')
      write (iw, '(/12X,''Heat capacity is Cp, not Cv'')')
      call web_message(iw,"thermochemistry.html")
!  END OF INSERT
      return
      end subroutine thermo
