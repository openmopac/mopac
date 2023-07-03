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

subroutine parfg (errors, ttype, nerr, loop, do_scf)
  !**********************************************************************
  !
  ! PARFG computes the values of the errors for molecule number LOOP
  !       using the current set of parameters.
  !
  !  On exit  ERRORS = Array of errors for reference functions
  !           NERR   = Number of errors in ERRORS
  !           TTYPE  = (text): name of the type of reference function.
  !
  !**********************************************************************
    use molkst_C, only : last, emin, escf, nclose, nopen, fract, &
    & nscf,  nvar, numat, rjkab1
    use param_global_C, only : ifiles_8, ihrefs, hofcal, refgeo, refher, &
     refder, refcer, refger, reftot, wtgeo, wthof, refhof, wtips, refips, &
    & wtdip, refdip, wtz, refpKa, wtpKa, atom_pKa, atom_pKa_O, &
    spin_state, i_r, qn
    use parameters_C, only : uss, upp, udd, zs, tore
    use meci_C, only : eig, lab, ispin, nelec, labsiz
    use funcon_C, only : fpc_9
    use symmetry_C, only : namo, jndex
    use cosmo_C, only: useps, area
    use common_arrays_C, only :  nat, coord, nfirst, nlast,  uspd, &
    & xparam, grad, eigs, loc, p, c
    implicit none
    double precision, dimension (300), intent (out) :: errors
    character (len=20), dimension (300), intent (out) :: ttype
    integer, intent (out) :: nerr
    integer, intent (in) :: loop
    logical :: do_scf
    integer :: i, ndif, nopn, loop_ref, j, qn_temp(20,20), irep(20), m, l, k, &
      iloop = -10
    double precision :: caldip, calhof, calips, cerr, derr, gerr, herr, toterr, &
    & sum
    double precision, dimension (300) :: q
    character , dimension(10) :: tspin*8
    double precision, external :: meci, parips, pardip, pargeo
    character :: srep(1000)*4
    save :: ndif, nopn, loop_ref, iloop
    data tspin/ 'SINGLET ', 'DOUBLET ', 'TRIPLET ', 'QUARTET ', 'QUINTET ', &
      'SEXTET  ', 'SEPTET  ', 'OCTET   ', 'NONET   ', '??????? '/
!-----------------------------------------------------------------------
  !
  !  Need to fill USPD - this would normally be filled in MOLDAT...
  !
     call filusp (nat, nfirst, nlast, uspd)
  !***********************************************************************
  !
  !  Compute the Function
  !
  !  set last to 3 to force MECI to calculate all roots. (Dataset must have
  !  keyword MECI in order for this to work)
  !
    last = 3
    emin = escf
    nscf=0
    calhof = 0.d0
    if (do_scf) then
      grad(1:nvar) = 0.d0
      call compfg (xparam, .true., calhof, .true., grad, wtgeo > 1.d-4)
      escf = calhof
    else
      calhof = escf
    end if
    toterr = 0.d0
    nerr = 0
    gerr = 0.d0
!
!  wtz = 1.0 if not a solid.
!  wtz = 1.0/N where N = number of stoichiometric units in cluster, if a solid.
!
!
    calhof = calhof*wtz
    if(loop /= 0) loop_ref = loop
  !***********************************************************************
  !
  !                     Heat of Formation Error
  !
  !***********************************************************************
    if (wthof > 1.d-8) then
      if(qn(loop_ref) /= 0)then
      if (labsiz == 11) sum = meci()
        call symtrz (c(1, nelec + 1), eig, 3, .TRUE.)
        qn_temp = 0
        j = 0
        l = 0
        do i = 1, lab
          if (i > 1) then
    !
    !  If part of a degenerate manifold, use previous quantum number
    !
            if (ispin(i) == ispin(i-1) .and. namo(i) == namo(i-1) .and. abs(eig(i) - eig(i - 1)) < 0.01d0) then
              jndex(i) = jndex(i-1)
              cycle
            end if
          end if
    !
    !  Find Irreducible Representation
    !
          do k = 1, j
            if (srep(k) == namo(i)) exit
          end do
          if (k > j) then
            j = j + 1
            srep(j) = namo(i)
          end if
  !
          do m = 1, l
            if (irep(m) == ispin(i)) exit
          end do
          if (m > l) then
            l = l + 1
            irep(l) = ispin(i)
          end if
          qn_temp(k, m) = qn_temp(k, m) + 1
          jndex(i) = qn_temp(k, m)
        end do
!
!  Use excited state relative energy.  First, identify desired state.
!
        do i = 1,lab
          if( qn(loop_ref) /= jndex(i) ) cycle
          if( spin_state(loop_ref) /= ispin(i) ) cycle
          if( i_r(loop_ref) /= namo(i) ) cycle
!
!  State "i" is the one needed.
!
           calhof = (eig(i) - eig(1))*fpc_9 + 1.d-15
!
!  Extra "1.d-15" to ensure that calhof is not zero
!
          exit
        end do
        if(i > lab) then
          if (iloop /= loop .and. loop > 0) then
            write(ifiles_8,"(a,i1,a,i1,2a)") &
      & " Root requested not found: ", &
      & qn(loop_ref),",",spin_state(loop_ref),",",i_r(loop_ref)
            write(ifiles_8,"(/,a,/)")" Possible roots"
            write(ifiles_8,"(i5,2x,a,2x,a,f12.4)")(jndex(i), tspin(ispin(i)), namo(i), eig(i), i = 1, lab)
            iloop = loop
          end if
!
!  Instead of stopping the run, set the calculated state energy to the observed
!  state energy.  That will give rise to a zero error.  On the next cycle,
!  the failure should correct itself, or the user might see the error
!  message and fix the problem.
!
          calhof = refhof
        end if
      end if
      herr = calhof - refhof + 1.d-15 !(Ensure that the H.o.F. error is not exactly zero)
      nerr = nerr + 1
      errors(nerr) = herr * wthof
      ttype (nerr) = "HOF "
    else
      herr = 0.d0
    end if
    if (wtpKa > 1.d-8) then
!
!
!   Use dummy atom parameters for parameters used in unconventional ways.
!
!  Here, pKa = (Charge on ionizable hydrogen)*(constant_1) + (constant_2)
      call chrge (p, q)
      if (useps) then
        if (atom_pKa_O == 0) then
          herr = 1.d6
          do i = 1, numat
            if (atom_pKa == i) cycle
            sum = (coord(1,i) - coord(1,atom_pKa))**2 + &
         & (coord(2,i) - coord(2,atom_pKa))**2 + &
         & (coord(3,i) - coord(3,atom_pKa))**2
            if (sum < herr) then
              herr = sum
              atom_pKa_O = i
            end if
          end do
        end if
        sum = Sqrt((coord(1,atom_pKa_O) - coord(1,atom_pKa))**2 + &
          &      (coord(2,atom_pKa_O) - coord(2,atom_pKa))**2 + &
          &      (coord(3,atom_pKa_O) - coord(3,atom_pKa))**2 )
    !    herr = (tore(nat(atom_pKa)) - q(atom_pKa))*uss(99) +  sum*upp(99) + zs(99) - refpKa
        herr = (tore(nat(atom_pKa)) - q(atom_pKa))*uss(99) + area*udd(99) + sum*upp(99) + zs(99) - refpKa
    !  if (loop /= 0) write(ifiles_8,'(a,f8.4)')" Charge:", (tore(nat(atom_pKa)) - q(atom_pKa))
      else
        herr = (tore(nat(atom_pKa)) - q(atom_pKa))*uss(99) + upp(99)  - refpKa
      end if
      nerr = nerr + 1
      errors(nerr) = herr * wtpKa
      ttype (nerr) = "pKa "
    end if
  !***********************************************************************
  !
  !                     Ionization Potential Error
  !
  !***********************************************************************
    if (wtips > 1.d-8) then
      i = nclose
      if (fract > 1.99d0) then
        i = nopen
      end if
      nopn = nopen - i
    !   CORRECTION TO I.P. OF DOUBLETS
      calips = parips(eigs, loop_ref)
      if (nopn == 1) then
        calips = calips + 0.5d0 * rjkab1
      end if
      if (Abs(calips) > 1.d-5) then
        cerr = calips - refips
        nerr = nerr + 1
        errors(nerr) = cerr * wtips
        ttype (nerr) = "I.P."
      else
        cerr = 0.d0
      end if
    else
      cerr = 0.d0
    end if
  !***********************************************************************
  !
  !                       Dipole Moment Error
  !
  !***********************************************************************
    if (wtdip > 1.d-8) then
      caldip = pardip(coord, nat)
      derr = caldip - refdip
      nerr = nerr + 1
      errors(nerr) = derr * wtdip
      ttype (nerr) = "DIPO"
    else
      derr = 0.d0
    end if
  !***********************************************************************
  !
  !                 Geometry errors (grad)
  !
  !***********************************************************************
    if (wtgeo > 1.d-8) then
      gerr = pargeo(grad, wtgeo / wtz**0.666d0, refgeo, loc, errors, ndif)
      nerr = nerr + ndif
      j = 0
      do i = 1, Min(nvar,30)
        if (refgeo(i) /= " ") then
          j = j + 1
          ttype (j) = "GEOM "//refgeo(i)
        end if
      end do
    else
      gerr = 0.d0
    end if
  !***********************************************************************
      if (loop > 0) then
      hofcal(loop) = calhof
      refher = herr
      do i = 1, ihrefs(1, loop)
        refher = refher - hofcal(ihrefs(i+1, loop))
      end do
      refder = derr
      refcer = cerr
      refger = gerr
      reftot = toterr
    end if
end subroutine parfg
