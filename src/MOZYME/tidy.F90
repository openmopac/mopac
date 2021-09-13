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

subroutine tidy (nmos_loc, nc, ic, n01, c, n02, nnc_loc, ncmo, ln, mn, mode)
    use MOZYME_C, only: iorbs, jopt, thresh, numred
    use molkst_C, only: numat, step_num, norbs, moperr, keywrd, numcal
    use chanel_C, only: iw
    implicit none
    integer, intent (in) :: mode,  n02, nmos_loc
    integer, intent (inout) :: n01
    integer, intent (out) :: ln, mn
    integer, dimension (n01), intent (inout) :: ic
    integer, dimension (nmos_loc), intent (inout) :: nc, ncmo, nnc_loc
    double precision, dimension (n02), intent (inout) :: c
!
    integer :: alloc_stat
    logical, save :: debug, large
    integer :: icalcn = 0
    integer, save :: ireset
    integer :: i, in, isnew, ispace, j, j1, jdash, jsav, jspace, jtop, &
         & k, l, li, ll, mdash, mm, momax, momin, msav, mtop, n
    integer :: nmol = 0
    double precision :: sum
    integer, dimension(:), allocatable :: iused, jused, kused, lused
    integer, dimension (2) :: imode
    data imode / 2 * 0 /
    allocate (iused(norbs), jused(norbs), kused(norbs), lused(norbs), &
         & stat=alloc_stat)
    if (alloc_stat /= 0) then
      call memory_error ("tidy")
      goto 100
    end if
   !
    if (numcal /= icalcn) then
      icalcn = numcal
      debug = (Index (keywrd, " TIDY") /= 0)
      large = (Index (keywrd, " LARGE") /= 0)
    end if
    if (nmos_loc == 0) return
    isnew = 0
    if (nmol /= numcal) then
      !
      !  A new molecule.  Therefore, set IMODE to step_num to prevent
      !  SELMOS from being called.
      !
      imode(1) = step_num
      imode(2) = step_num
      ireset = 0
    end if
    do
      !
      !  COMPRESS  IC AND C
      !
      momax = 0
      ireset = ireset - 1
      if (ireset == 0) then
        thresh = thresh * 100.d0
      end if
      momin = 10000
      ln = 0
      mn = 0
      in = 0
      li = 0
      do i = 1, nmos_loc
        ll = nnc_loc(i)
        li = li - ll
        mm = ncmo(i)
        ncmo(i) = mn
        nnc_loc(i) = 0
        if (i /= 1) then
          nnc_loc(i) = nnc_loc(i-1) + in
        end if
        in = 0
        momax = Max (momax, nc(i))
        if (nc(i) < momin) then
          momin = nc(i)
        end if
        do j1 = 1, nc(i)
          ll = ll + 1
          sum = 0.d0
          j = ic(ll)
          do k = 1, iorbs(j)
            mm = mm + 1
            sum = sum + c(mm) ** 2
          end do
          if (sum > thresh) then
            mm = mm - iorbs(j)
            do k = 1, iorbs(j)
              mn = mn + 1
              mm = mm + 1
              c(mn) = c(mm)
            end do
            in = in + 1
            ln = ln + 1
            ic(ln) = j
          end if
        end do
        nc(i) = in
        iused(i) = mn
        li = li + ll
      end do
      !
      !   Do the LMOs of the SCF need to be put at the start of the storage?
      !
      if (isnew /= 0 .or. step_num <= 1 .or. step_num == imode(mode)) exit
      isnew = 2
      if (numred >= numat-1) then
        isnew = 1
      end if
      imode(mode) = step_num
      if (isnew == 2) then
         !
         !   Move LMO's for atoms in JOPT to start of storage space
         !
        call selmos (nmos_loc, nc, ic, n01, c, n02, nnc_loc, ncmo, ln, mn, &
       & iused, jopt, jused, kused, lused, mode)
      else
        exit
      end if
    end do
   !
   !   DIVIDE AVAILABLE SPACE EQUALLY AMONG THE LMO'S
   !
    ispace = (n01-ln) / nmos_loc
    jspace = (n02-mn) / nmos_loc
   !
   !   LN     = Highest address in IC used
   !            (the end of the compressed LMOs)
   !   LN/NMOS= Average number of atoms in a LMO.
   !   ISPACE = Average amount of free space that LMOs can expand into.
   !
   !   Is there enough space for LMOs to expand?  If not, stop the run
   !
   !  In order to ensure that the next annihilation step can be run,
   !  the amount of storage for each LMO, (LN/NMOS+ISPACE), should be
   !  more than NUMAT, or ISPACE must be greater than the larger of
   !  20% of the available space and 20. That is, there must be at least
   !  20% unused space for a LMO to expand into.
   !
    i = n01 / nmos_loc
    if (i + ispace < numat .and. ispace < Max (i/5, 20)) then
      call pinout (1, .false.)
  !    write (iw, "(/,A,/,A,I8)") &
  !   & " There is not enough unused space left to ensure ", &
  !   & " that the next iteration can be done.  New memory being allocated:", &
  !   & Nint (n01*0.2)
      moperr = .true.
      return
    end if
    if (ispace > numat .and. jspace > nmos_loc) then
      !
      !   There is more than enough space to guarantee that every LMO can
      !   be filled completely, therefore:
      ispace = 0
      jspace = 0
      if (debug) then
        write (iw, "(6(A,I5))") " MAX.:", momax, " MIN.:", momin, " AVE. IC:", &
       & ln / nmos_loc, " AVE. C:", mn / nmos_loc
      end if
    else if (debug) then
      write (iw, "(6(A,I5))") " MAX.:", momax, " MIN.:", momin, " AVE. IC:", &
           & ln / nmos_loc, " AVE. C:", mn / nmos_loc, " FREE IC:", &
           & ispace, " FREE C:", jspace
    end if
    jtop = n01
    mtop = n02
    jsav = n01 - ln - ispace * nmos_loc
    msav = n02 - mn - jspace * nmos_loc
    do i = nmos_loc, 2, -1
      j = Min (nc(i)+ispace, numat)
      if (j < nc(i)+ispace) then
         !
         !   Some space is unused - save it for later
         !
        jsav = jsav + nc(i) + ispace - j
        jtop = jtop - j
      else
         !
         !  USE UP SOME OF THE SAVED SPACE
         !
        jdash = Min (numat, Max (jsav, 0)+j)
        jtop = jtop - jdash
        jsav = jsav - jdash + j
      end if
      l = nnc_loc(i)
      nnc_loc(i) = jtop
      do k = nc(i), 1, -1
        ic(jtop+k) = ic(l+k)
      end do
      n = iused(i) - iused(i-1)
      j = Min (n+jspace, norbs)
      if (j < n+jspace) then
         !
         !   Some space is unused - save it for later
         !
        msav = msav + n + jspace - j
        mtop = mtop - j
      else
         !
         !  USE UP SOME OF THE SAVED SPACE
         !
        mdash = Min (norbs, Max (msav, 0)+j)
        mtop = mtop - mdash
        msav = msav - mdash + j
      end if
      ncmo(i) = mtop
      do k = n, 1, -1
        c(k+mtop) = c(k+iused(i-1))
      end do
    end do
    if (debug .and. large) then
      !
      !   PRINT ADDRESSES SET UP BY TIDY
      !
      write (iw,*) " ADDRESSES IN TIDY"
      write (iw, "(A)") "               NNC       NC     NCMO "
      do i = 1, nmos_loc
        write (iw, "(7I9)") i, nnc_loc (i), nc (i), ncmo (i)
      end do
      write (iw,*) " Contents of IC"
      do i = 1, nmos_loc
        write (iw, "(I4,I7,18I4,10(/,I11,18I4))") i, (ic(j), j=nnc_loc(i)+1, &
       & nnc_loc(i)+nc(i))
      end do
    end if
    if (mode == 2) then
      nmol = numcal
    end if
   !
    deallocate (iused, jused, kused, lused)
100 continue
end subroutine tidy
