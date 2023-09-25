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

      subroutine deriv(geo, gradnt)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE molkst_C, ONLY: numat, norbs, nclose, nopen, fract, natoms, numcal, &
      & ndep, nvar, keywrd, cosine, moperr, mpack, isok, id, l123, line, nscf, &
      pressure, l1u, l2u, l3u, method_PM7, method_pm8, method_pm6_org
      use common_arrays_C, only : dxyz, loc, errfn, aicorr, tvec
      USE symmetry_C, ONLY: locpar, idepfn
      USE chanel_C, only : iw, ir, job_fn
      use funcon_C, only : pi
      use derivs_C, only : aidref, work2
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
!
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      double precision, intent(in)  :: geo(3,natoms)
      double precision, intent(inout)  :: gradnt(nvar)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: idelta, icalcn, i, j, nstep, nw2, ncol
      double precision, dimension(3) :: change
      double precision :: coord(3*natoms), gold(3*natoms), xparam(3*natoms)


      double precision :: grlim, sum, gnorm, step, press
      double precision, dimension(3,3) :: tderiv
      double precision, external :: dot, volume
      logical :: scf1, halfe, slow, aifrst, debug, precis, intn, geochk, ci, &
        aic, noanci, field, saddle, DH_correction, l_redo_bonds

      save change, scf1, halfe, idelta, slow, icalcn, aifrst, debug, l_redo_bonds, &
        precis, intn, geochk, ci, aic, grlim, nw2, field, DH_correction
!-----------------------------------------------
!***********************************************************************
!
!    DERIV CALCULATES THE DERIVATIVES OF THE ENERGY WITH RESPECT TO THE
!          INTERNAL COORDINATES. THIS IS DONE BY FINITE DIFFERENCES.
!
!    THE MAIN ARRAYS IN DERIV ARE:
!        LOC    INTEGER ARRAY, LOC(1,I) CONTAINS THE ADDRESS OF THE ATOM
!               INTERNAL COORDINATE LOC(2,I) IS TO BE USED IN THE
!               DERIVATIVE CALCULATION.
!        GEO    ARRAY \GEO\ HOLDS THE INTERNAL COORDINATES.
!        gradnt   ON EXIT, CONTAINS THE DERIVATIVES
!
!***********************************************************************
      data icalcn/ 0/
      if (icalcn /= numcal) then
        aifrst = index(keywrd,' RESTART') == 0
        saddle = index(keywrd, " SADDLE") /= 0
        debug = index(keywrd,' DERIV') /= 0
        field = index(keywrd,' FIELD') /= 0
        precis = index(keywrd,' PREC') /= 0
        DH_correction = (index(keywrd,' PM6-D') + index(keywrd,' PM6-H') /= 0 .or. &
          method_PM7 .or. method_pm6_org .or. method_pm8 )
        intn = index(keywrd,'  XYZ') == 0
        if (saddle) then
          nw2 = max(mpack, 9*natoms**2)
        else
          nw2 = max(mpack, 6*natoms*l123)
        end if
        errfn = 0.d0
        if (allocated(aidref)) deallocate(aidref)
        if (allocated(work2))  deallocate(work2)
        allocate(aidref(nvar), work2(nw2))
        aidref = 0.d0
!
!   GEOCHK is true if the system is a transition state in internal
!          coordinates with all coordinates marked for optimization.
!          It can be over-ridden by GEO-OK
!
        geochk = index(keywrd,' TS') + index(keywrd,' NLLSQ') + index(keywrd,&
          ' SIGMA') /= 0
        geochk = geochk .and. intn .and. nvar>=numat*3-6 .and. id==0 .and. &
          index(keywrd,'GEO-OK')==0
        geochk = geochk .and. index(keywrd,' XYZ')==0
        ci = index(keywrd,' C.I.') /= 0
        scf1 = index(keywrd,' 1SCF') /= 0
        aic = index(keywrd,'AIDER') /= 0
        if (aic .and. aifrst) then
          open(unit=ir, file=job_fn, status='OLD', blank='ZERO',position='asis')
          rewind ir
!
!  ISOK IS SET FALSE: ONLY ONE SYSTEM ALLOWED
!
          isok = .FALSE.
          do i = 1, 1000
            read (ir, '(A)') line
            call upcase (line, 80)
            if (index(line,'AIDER') == 0) cycle
            exit
          end do
          do j = 1, 1000
            read (ir, '(A)', end=40, err=40) line
            call upcase (line, 80)
            if (index(line,'AIDER') /= 0) go to 60
          end do
   40     continue
          write (iw, '(2/,A)') ' KEYWORD "AIDER" SPECIFIED, BUT NOT'
          write (iw, '(A)') ' PRESENT AFTER Z-MATRIX.  JOB STOPPED'
          call mopend (&
         'KEYWORD "AIDER" SPECIFIED, BUT NOT PRESENT AFTER Z-MATRIX.  JOB STOPPED')
          return
   50     continue
          write (iw, '(2/,A)') '  FAULT IN READ OF AB INITIO DERIVATIVES'
          write (iw, '(A)') '  DERIVATIVES READ IN ARE AS FOLLOWS'
          write (iw, '(6F12.6)') (aidref(j),j=1,i)
          call mopend ('FAULT IN READ OF AB INITIO DERIVATIVES')
          return
   60     continue
          if (natoms > 2) then
            j = 3*natoms - 6
          else
            j = 1
          end if
          read (ir, *, end=50, err=50) (aidref(i),i=1,j)
          write (iw, '(/,A,/)') &
            ' AB-INITIO DERIVATIVES IN KCAL/MOL/(ANGSTROM OR RADIAN)'
          write (iw, '(5F12.6)') (aidref(i),i=1,j)
          do i = 1, nvar
            if (loc(1,i) > 3) then
              j = 3*loc(1,i) + loc(2,i) - 9
            else if (loc(1,i) == 3) then
              j = loc(2,i) + 1
            else
              j = 1
            end if
            aidref(i) = aidref(j)
          end do
          write (iw, '(/,A,/)') ' AB-INITIO DERIVATIVES FOR VARIABLES'
          write (iw, '(5F12.6)') (aidref(i),i=1,nvar)
          if (ndep /= 0) then
            do i = 1, nvar
              sum = aidref(i)
              aidref(i) = aidref(i) + count(loc(1,i)==locpar(:ndep) .and. &
              (loc(2,i)==idepfn(:ndep) .or. loc(2,i)==3 .and. &
              idepfn(:ndep)==14))*sum
            end do
            write (iw, '(/,A,/)') &
              ' AB-INITIO DERIVATIVES AFTER SYMMETRY WEIGHTING'
            write (iw, '(5F12.6)') (aidref(j),j=1,nvar)
          end if
          close(ir, status='KEEP')
        end if
        l_redo_bonds = (index(keywrd,' FORCE') + index(keywrd,' IRC=') + &
          index(keywrd,' THERM') + index(keywrd,' DFORCE') == 0)
        grlim = 0.01D0
        if (precis) grlim = 0.0001D0
        halfe = nopen>nclose .and. Abs(fract - 2.d0) > 1.d-20 .and. Abs(fract) > 1.d-20 .or. ci
        idelta = -7
!
!   IDELTA IS A MACHINE-PRECISION DEPENDANT INTEGER
!
        change(1) = 10.D0**idelta
        change(2) = 10.D0**idelta
        change(3) = 10.D0**idelta
!
!    CHANGE(I) IS THE STEP SIZE USED IN CALCULATING THE DERIVATIVES.
!    FOR "CARTESIAN" DERIVATIVES, CALCULATED USING DCART,AN
!    INFINITESIMAL STEP, HERE 0.000001, IS ACCEPTABLE. IN THE
!    HALF-ELECTRON METHOD A QUITE LARGE STEP IS NEEDED AS FULL SCF
!    CALCULATIONS ARE NEEDED, AND THE DIFFERENCE BETWEEN THE TOTAL
!    ENERGIES IS USED. THE STEP CANNOT BE VERY LARGE, AS THE SECOND
!    DERIVITIVE IN FLEPO IS CALCULATED FROM THE DIFFERENCES OF TWO
!    FIRST DERIVATIVES. CHANGE(1) IS FOR CHANGE IN BOND LENGTH,
!    (2) FOR ANGLE, AND (3) FOR DIHEDRAL.
!
      end if
      if (nvar == 0) return
      if (debug) then
        write (iw, '(10X, "GEOMETRY AT THE START OF DERIV")')
        call geout(-iw)
      end if
      gnorm = 0.D0
      do i = 1, nvar
        gold(i) = gradnt(i)
        xparam(i) = geo(loc(2,i),loc(1,i))
        gnorm = gnorm + gradnt(i)**2
      end do
      gnorm = sqrt(gnorm)
      slow = .FALSE.
      noanci = .FALSE.
      if (halfe) then
        noanci = (index(keywrd,'NOANCI') /= 0 .or. nopen == norbs)
        slow = (noanci .and. (gnorm<grlim .or. scf1))
      else
        slow = (index(keywrd,'NOANCI') /= 0)
      end if
      if (ndep /= 0) call symtry
      call gmetry (geo, coord)
!
!  COORD NOW HOLDS THE CARTESIAN COORDINATES
!
      if (halfe .and. .not.noanci .and. numat > 1) then
        if (debug) write (iw, '(10x,a)') 'DOING ANALYTICAL C.I. DERIVATIVES'
        if (debug) write (iw, '(" NUMBER  ATOM  ",5X,"X",12X,"Y",12X,"Z",/)')
        call dernvo ()
        if (moperr) return
      else
        if (debug) write (iw, '(10x,a)') 'DOING VARIATIONALLY OPTIMIZED DERIVATIVES'
        call dcart (coord, dxyz)
      end if
      if (l_redo_bonds .and. mod(nscf,10) == 4 .and. nscf /= 0 .and. id == 0) then
!
!  There is a possibility that a bond might be made or broken during a geometry optimization
!  or a gradient minimization.  To allow for this, the bonds array should be updated every
!  few SCF calculations. The values in the mod test are "intelligent guesses"
!
        call lewis (.true.)
        if (moperr) then
           write (iw, '(/10x,A,/)') ' Geometry at the point this error was detected'
          call geout(iw)
        end if
      end if
      if (DH_correction) call post_scf_corrections(sum, .true.)
!
!   THE CARTESIAN DERIVATIVES ARE IN DXYZ
!
      if (field) call dfield ()
      if (Abs (pressure) > 1.d-4) then
        if (id == 1) then
      !
      !  Add in gradient tension contribution
      !
          i = 3 * l123 / 2 - 1
        !
        !  For polymers, pressure is pull in Newtons for 1 mole
        !  1N = J/M = 10**(-3)/4.184 kcal/M = 4.184*10**(-3)*10**(-10) kcal/Angstrom
          press = pressure / Sqrt (dot(tvec(1, 1), tvec(1, 1), 3))
          do j = 1, 3
            dxyz(j + i) = dxyz(j + i) - tvec(j, 1) * press
            dxyz(j + i + 3) = dxyz(j + i + 3) + tvec(j, 1) * press
          end do
        else if (id == 3) then
  ! Transition vector derivatives of pressure times volume
  !
  ! Derivatives are calculated in a redundant coordinate system where the
  ! atomic coordinates in the central cell plus translation vectors are
  ! replaced by atomic coordinates in the central cell plus atomic coordinates
  ! in all coupled non-central cells. The formula for volume that is being
  ! differentiated uses the coordinate difference between a central atom and
  ! its neighbors as proxies for the translation vectors in the volume formula.
  ! The choice of atom in this formula and the resulting derivatives are completely
  ! arbitrary and do not change the derivatives in the original coordinate system.
          press = ((tvec(2, 1)*tvec(3, 2)-tvec(3, 1)*tvec(2, 2))*tvec(1, 3) + &
                  (tvec(3, 1)*tvec(1, 2)-tvec(1, 1)*tvec(3, 2))*tvec(2, 3) + &
                  (tvec(1, 1)*tvec(2, 2)-tvec(2, 1)*tvec(1, 2))*tvec(3, 3))
          tderiv(1, 1) = tvec(2, 2)*tvec(3, 3) - tvec(3, 2)*tvec(2, 3)
          tderiv(2, 1) = tvec(3, 2)*tvec(1, 3) - tvec(1, 2)*tvec(3, 3)
          tderiv(3, 1) = tvec(1, 2)*tvec(2, 3) - tvec(2, 2)*tvec(1, 3)
          tderiv(1, 2) = tvec(3, 1)*tvec(2, 3) - tvec(2, 1)*tvec(3, 3)
          tderiv(2, 2) = tvec(1, 1)*tvec(3, 3) - tvec(3, 1)*tvec(1, 3)
          tderiv(3, 2) = tvec(2, 1)*tvec(1, 3) - tvec(1, 1)*tvec(2, 3)
          tderiv(1, 3) = tvec(2, 1)*tvec(3, 2) - tvec(3, 1)*tvec(2, 2)
          tderiv(2, 3) = tvec(3, 1)*tvec(1, 2) - tvec(1, 1)*tvec(3, 2)
          tderiv(3, 3) = tvec(1, 1)*tvec(2, 2) - tvec(2, 1)*tvec(1, 2)
          tderiv = tderiv*pressure*sign(1.d0, press)
  !
  ! Add in gradient pressure term
  !
  !  First, the central unit cell
  !
          i = 3*(l1u * (2*l2u+1) * (2*l3u+1) + l2u * (2*l3u+1) + l3u - 1)
          do j = 1, 3
            dxyz(j + i) = dxyz(j + i) - tderiv(j, 1)
            dxyz(j + i) = dxyz(j + i) - tderiv(j, 2)
            dxyz(j + i) = dxyz(j + i) - tderiv(j, 3)
          end do
  !
  !  Cell in 0,0,1 position
  !
          i = 3*(l1u * (2*l2u+1) * (2*l3u+1) + l2u * (2*l3u+1) + l3u)
          do j = 1, 3
            dxyz(j + i) = dxyz(j + i) + tderiv(j, 3)
          end do
  !
  !  Cell in 0,1,0 position
  !
          i = 3*(l1u * (2*l2u+1) * (2*l3u+1) + (l2u+1) * (2*l3u+1) + l3u - 1)
          do j = 1, 3
            dxyz(j + i) = dxyz(j + i) + tderiv(j, 2)
          end do
  !
  !  Cell in 1,0,0 position
  !
          i = 3*((l1u+1) * (2*l2u+1) * (2*l3u+1) + l2u * (2*l3u+1) + l3u - 1)
          do j = 1, 3
            dxyz(j + i) = dxyz(j + i) + tderiv(j, 1)
          end do
        end if
      end if
      if ((DH_correction .or. Abs (pressure) > 1.d-4) .and. debug) then
         call print_dxyz("  (Includes post-SCF corrections)")
      end if

      step = change(1)
      nstep = nw2/(3*numat*l123)
      do i = 1, nvar, nstep
        j = min(i + nstep - 1,nvar)
        call jcarin (xparam, step, precis, work2, ncol, i, j)
        call mxm (work2, j - i + 1, dxyz, ncol, gradnt(i), 1)
!
      end do
      if (precis) then
        step = 0.5D0/step
      else
        step = 1.0D0/step
      end if
      gradnt(:nvar) = gradnt(:nvar)*step
!
!  NOW TO ENSURE THAT INTERNAL DERIVATIVES ACCURATELY REFLECT CARTESIAN
!  DERIVATIVES
!
      if (geochk) then
        sum = dot(gradnt,gradnt,nvar)
        if (sum<2.D0 .and. dot(dxyz,dxyz,3*numat)>max(4.D0,sum*4.D0)) then
!
! OOPS, LOOKS LIKE AN ERROR.
!
          do i = 1, nvar
            j = int(xparam(i)/3.141D0)
            if (.not.(loc(2,i)==2 .and. loc(1,i)>3 .and. abs(xparam(i)-j*pi)<0.005D0)) cycle
!
!  ERROR LOCATED, BUT CANNOT CORRECT IN THIS RUN
!
            write (iw, '(2/,3(A,/),I3,A)') &
              ' INTERNAL COORDINATE DERIVATIVES DO NOT REFLECT', &
              ' CARTESIAN COORDINATE DERIVATIVES', &
              ' TO CORRECT ERROR, INCREASE DIHEDRAL OF ATOM', loc(1,i), &
              ' BY 90 DEGREES'
            write (iw, '(2/,A)') '     CURRENT GEOMETRY'
            call geout (iw)
            call mopend (&
       ' INTERNAL COORDINATE DERIVATIVES DO NOT REFLECT CARTESIAN COORDINATE DERIVATIVES')
            return
          end do
        end if
      end if
!
!  THIS CODE IS ONLY USED IF THE KEYWORD NOANCI IS SPECIFIED
      if (slow) then
        if (debug) write (iw, *) 'DOING FULL SCF DERIVATIVES'
        call deritr ()
        icalcn = numcal
!
! THE ARRAY ERRFN HOLDS THE EXACT DERIVATIVES MINUS THE APPROXIMATE
! DERIVATIVES
        errfn(:nvar) = errfn(:nvar) - gradnt(:nvar)
      end if
!
!  AT THIS POINT, THE INTERNAL DERIVATIVES ARE KNOWN, AND ARE IN GRADNT.
!
      cosine = dot(gradnt,gold,nvar)/sqrt(dot(gradnt,gradnt,nvar)*dot(gold,gold,nvar)&
         + 1.D-20)
      if (slow) gradnt(:nvar) = gradnt(:nvar) + errfn(:nvar)
      if (aic) then
        if (aifrst) then
          aifrst = .FALSE.
          aicorr(:nvar) = (-aidref(:nvar)) - gradnt(:nvar)
        end if
        gradnt(:nvar) = gradnt(:nvar) + aicorr(:nvar)
      end if
      if (debug) then
        write (iw, '('' GRADIENTS'')')
        write (iw, '(8F10.3)') (gradnt(i),i=1,nvar)
        write (iw, '('' XPARAM'')')
        write (iw, '(8F10.3)') (xparam(i),i=1,nvar)
        if (slow) then
          write (iw, '('' ERROR FUNCTION'')')
          write (iw, '(8F10.3)') (errfn(i),i=1,nvar)
        end if
      end if
      if (debug) write (iw, '('' COSINE OF SEARCH DIRECTION ='',F30.6)') cosine
      return
      end subroutine deriv
