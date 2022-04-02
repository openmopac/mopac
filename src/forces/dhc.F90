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

subroutine dhc (p, pa, pb, xi, nat, if, il, jf, jl, dener, mode)
    use molkst_C, only : numcal, cutofp, id, uhf, clower
    implicit none
    integer, intent (in) :: if, il, jf, jl, mode
    double precision, intent (out) :: dener
    integer, dimension (2), intent (in) :: nat
    double precision, dimension (*), intent (in) :: p, pa, pb
    double precision, dimension (3, 2), intent (in) :: xi
    logical, save :: ignor1, ignor2, lpoint, cutoff
    integer, save :: icalcn = 0
    double precision, save :: cutof
    double precision, save :: wlim
    integer :: ic
    integer :: i, i1, i2, ia, ii, j, j1, ja, jc, jj, k, kr, linear, ni, nj
    double precision :: ee, enuclr_loc, rij
    integer, dimension (2) :: nfirst_loc, nlast_loc
    double precision, dimension (45) :: e1b = 0.d0, e2a = 0.d0
    double precision, dimension (171) :: f, h
    double precision, dimension (9, 9) :: shmat
    double precision, dimension (2026), save :: w = 0.d0, wk = 0.d0
    double precision, external :: helect
    if (icalcn /= numcal) then
      !
      !  Switch to point charges for separations of more than 14 Angstroms
      !  (This is related to FACT in subroutine REPPD)
      !
      cutof = clower**2
      lpoint = .false.
      icalcn = numcal
      ignor1 = .false.
      ignor2 = .false.
      wlim = 4.d0
      if (id == 0) then
        wlim = 0.d0
      end if
    end if
    rij = (xi(1, 1)-xi(1, 2)) ** 2 + (xi(2, 1)-xi(2, 2)) ** 2 &
         & + (xi(3, 1)-xi(3, 2)) ** 2
    if (mode == 1 .and. id /= 0) then
      ignor1 = (rij > 225.d0)  ! Ignore overlaps involving atoms separated by more than 15 A.
      ignor2 = (rij >(4.d0/3.d0*cutofp)**2)
      if (ignor1 .and. ignor2) then
        dener = 0.d0
        return
      end if
    end if
    if (mode == 2 .and. id /= 0 .and. ignor2 .and. ignor1) then
      dener = 0.d0
      return
    end if
    nfirst_loc(1) = 1
    nlast_loc(1) = il - if + 1
    nfirst_loc(2) = nlast_loc(1) + 1
    nlast_loc(2) = nfirst_loc(2) + jl - jf
    linear = (nlast_loc(2)*(nlast_loc(2)+1)) / 2
    do i = 1, linear
      f(i) = 0.d0
      h(i) = 0.0d00
    end do
    kr = 0
    ia = nfirst_loc(2)
    ic = nlast_loc(2)
    ja = nfirst_loc(1)
    jc = nlast_loc(1)
    j = 1
    i = 2
    nj = nat(1)
    ni = nat(2)
    if ( .not. ignor1) then
      !
      !   Evaluate the one-electron term
      !
      call h1elec (nj, ni, xi(1, 1), xi(1, 2), shmat)
      if (nat(1) == 102 .or. nat(2) == 102) then
        k = (ic*(ic+1)) / 2
        h(1:k) = 0.d0
      else
        j1 = 0
        do j = ia, ic
          jj = j * (j-1) / 2
          j1 = j1 + 1
          i1 = 0
          do i = ja, jc
            jj = jj + 1
            i1 = i1 + 1
            h(jj) = shmat(i1, j1)
            f(jj) = shmat(i1, j1)
          end do
        end do
      end if
    end if
    if ( .not. ignor2) then
      !
      !   Evaluate the two-electron term
      !
      if (id > 0) then
        if (mode == 1) then
          lpoint = (rij > cutof)
        end if
      end if
      if (lpoint) then
        rij = sqrt(rij)
        call point (rij, ni, nj, w, kr, e2a, e1b, enuclr_loc)
      else
        kr = 1
        call rotate (ni, nj, xi(1, 2), xi(1, 1), w, kr, e2a, e1b, &
       & enuclr_loc)
      end if
      if (id /= 0) then
        do i = 1, kr - 1
          wk(i) = w(i)
        end do
        if (mode == 1) then
          cutoff = (w(1) < wlim)
        end if
        if (cutoff) then
          do i = 1, kr - 1
            wk(i) = 0.d0
          end do
        end if
      end if
      !
      !    * ENUCLR IS SUMMED OVER CORE-CORE REPULSION INTEGRALS.
      !
      i2 = 0
      do i1 = ja, jc
        ii = i1 * (i1-1) / 2 + ja - 1
        do j1 = ja, i1
          ii = ii + 1
          i2 = i2 + 1
          h(ii) = h(ii) + e1b(i2)
          f(ii) = f(ii) + e1b(i2)
        end do
      end do
      i2 = 0
      do i1 = ia, ic
        ii = i1 * (i1-1) / 2 + ia - 1
        do j1 = ia, i1
          ii = ii + 1
          i2 = i2 + 1
          h(ii) = h(ii) + e2a(i2)
          f(ii) = f(ii) + e2a(i2)
        end do
      end do
      i = -2
      call fock2 (f, p, pa, w, w, wk, i, nfirst_loc, nlast_loc, 1)
    else
      enuclr_loc = 0.d0
    end if
    ee = helect (nlast_loc(2), pa, h, f)
    if (uhf) then
      do i = 1, linear
        f(i) = h(i)
      end do
      i = -2
      call fock2 (f, p, pb, w, w, wk, i, nfirst_loc, nlast_loc, 1)
      ee = ee + helect (nlast_loc(2), pb, h, f)
    else
      ee = ee * 2.d0
    end if
    dener = ee + enuclr_loc
end subroutine dhc
