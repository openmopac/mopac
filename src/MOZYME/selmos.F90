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

subroutine selmos (nmos_loc, nc, ic, n01, c, n02, nnc_loc, ncmo, ln, mn, iws, &
& jopt, ncnew, ncmnew, nncnew, mode)
    use molkst_C, only: numat
    use MOZYME_C, only : numred, nelred, norred
    use chanel_C, only: iw
    implicit none
    integer, intent (in) :: ln, mn, mode, n01, n02, nmos_loc
    integer, dimension (n01), intent (inout) :: ic
    integer, dimension (nmos_loc), intent (inout) :: iws, nc, ncmnew, ncmo, &
   & ncnew, nnc_loc, nncnew
    integer, dimension (numat), intent (in) :: jopt
    double precision, dimension (n02), intent (inout) :: c
    integer :: i, ibot, itop, j, jbot, jtop, l, l1, l2, m, n, nbot, ncoefs, &
   & nreal, ntop
   !
   !   Move LMOs for atoms in JOPT to start of storage space
   !
   !    MN, LN   = END OF STORAGE ALREADY USED
   !    N02, N01 = END OF STORAGE
   !
   !   FIRST, MOVE ALL VECTORS TO END OF STORAGE
   !
    l = n02 - mn
    do i = mn, 1, -1
      c(i+l) = c(i)
    end do
   !
   !  Modify pointers
   !
    do i = 1, nmos_loc
      ncmo(i) = ncmo(i) + l
    end do
   !
   !  Move atom addresses
   !
    l = n01 - ln
    do i = ln, 1, -1
      ic(i+l) = ic(i)
    end do
   !
   !  Modify atom address pointers
   !
    do i = 1, nmos_loc
      nnc_loc(i) = nnc_loc(i) + l
    end do
   !
   !        NBOT       refers to IC                -> IC
   !        JBOT, JTOP refer to C                  -> C
   !        IBOT, ITOP refer to NC, NNC, and NCMO -> NCNEW, NNCNEW, NCMNEW
   !
   !   Select those LMOs that are to be used in the SCF.  Move these to the
   !   start of the storage
   !
    jbot = 0
    nbot = 0
    jtop = n02 - mn - 1
    ntop = n01 - ln - 1
    ibot = 0
    itop = nmos_loc + 1
    nreal = 0
   !
   !   Put number of coefficients into IWS.
   !
    do i = nmos_loc, 2, -1
      iws(i) = iws(i) - iws(i-1)
    end do
    do i = 1, nmos_loc
      l = nnc_loc(i)
      l1 = ic(1+l)
      l2 = ic(2+l)
      do m = 1, numred
        if (l1 == jopt(m) .or. l2 == jopt(m)) go to 1000
      end do
      !
      ! PUT VECTOR AT END OF IC AND C.  This is the default: the LMO is
      ! NOT selected.
      !
      itop = itop - 1
      nncnew(itop) = nnc_loc(i)
      ncnew(itop) = nc(i)
      ncmnew(itop) = ncmo(i)
      cycle
1000  if (nbot+nc(i) > ntop .or. jbot+iws(i) > jtop) then
        !
        !   There is no space for the remaining vectors.  Remove unused space
        !   from between the vectors that will not be used in the SCF.
        !
        call compct (nncnew, ncnew, ncmnew, itop, nc, ic, iws, n01, c, n02, &
       & nmos_loc, i, ntop, jtop, l, ncmo(i))
      end if
      !
      ! PUT VECTOR AT START OF IC AND C
      !
      ibot = ibot + 1
      nncnew(ibot) = nbot
      ncnew(ibot) = nc(i)
      ncmnew(ibot) = jbot
      do n = 1, nc(i)
        ic(nbot+n) = ic(l+n)
      end do
      nbot = nbot + nc(i)
      nc(i) = 0
      j = ncmo(i)
      ncoefs = iws(i)
      do n = 1, ncoefs
        c(jbot+n) = c(j+n)
      end do
      jbot = jbot + ncoefs
      nreal = nreal + 1
    end do
    do i = 1, ibot
      nc(i) = ncnew(i)
      nnc_loc(i) = nncnew(i)
      ncmo(i) = ncmnew(i)
    end do
    j = nmos_loc + 1
    do i = itop, nmos_loc
      j = j - 1
      nc(i) = ncnew(j)
      nnc_loc(i) = nncnew(j)
      ncmo(i) = ncmnew(j)
    end do
   !
   !  HAVE ARRAY BOUNDS BEEN DAMAGED?
   !
    if (jbot > jtop) then
      if (mode == 1) then
        write (iw,*) " OCCUPIED C VECTOR DAMAGED IN SELMOS"
        write (iw,*) " (THIS IS A BUG IN THE PROGRAM)"
      else
        write (iw,*) " VIRTUAL C VECTOR DAMAGED IN SELMOS"
        write (iw,*) "  (THIS IS A BUG IN THE PROGRAM)"
      end if
      call mopend ("BUG IN SELMOS")
    end if
    if (nbot > ntop) then
      if (mode == 1) then
        write (iw,*) " OCCUPIED IC VECTOR DAMAGED IN SELMOS"
        write (iw,*) "  (THIS IS A BUG IN THE PROGRAM)"
      else
        write (iw,*) " VIRTUAL IC VECTOR DAMAGED IN SELMOS"
        write (iw,*) "  (THIS IS A BUG IN THE PROGRAM)"
      end if
      call mopend ("BUG IN SELMOS")
    end if
    if (mode == 1) then
      !
      !   Number of electrons to be used in the SCF is NELRED.
      !
      nelred = 2 * nreal
    else
      !
      !   Number of LMOs to be used in the SCF is NORRED.
      !
      norred = nreal + nelred / 2
    end if
    norred = norred
    nelred = nelred
end subroutine selmos
