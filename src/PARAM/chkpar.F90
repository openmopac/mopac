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

subroutine chkpar
    use param_global_C, only : ifiles_8, numvar, nfns, diffns, locvar
    use elemts_C, only : elemnt
    use parameters_C, only : partyp
    implicit none
    integer :: i, il, iu, j, j2, k, l, lim, ll, lu
    double precision :: sum
    double precision, dimension (:), allocatable :: par
    double precision, dimension (5) :: eiglim
    double precision, dimension (:), allocatable :: parmod
    double precision, dimension (:), allocatable :: vects
    intrinsic Min, Sqrt
  !
  !.. Data Declarations ..
    data eiglim / 1.d-5, 1.d-0, 1.d1, 1.d2, 1.d3 /
  !
  ! ... Executable Statements ...
  !
    write (ifiles_8, "(//20X,A)") " CHECK OF PARAMETER INDEPENDENCE"
    write (ifiles_8, "(//15X,A,I5)") " Number of reference data:", nfns
    write (ifiles_8, "(  15X,A,I5)") " Number of parameters:    ", numvar
     i = (numvar*(numvar+1))/2
     allocate (par(i), parmod(numvar), vects(numvar**2))
    if (nfns < 30 .and. numvar < 9) then
      write (ifiles_8, "(/3X,2A,/)") "Matrix of Differentials ", "(Rows: Refer&
     &ence data, Columns: Parameters)"
      do i = 1, nfns
        write (ifiles_8, "(8F10.4)") (diffns(j, i), j=1, numvar)
      end do
    end if
    l = 0
    do i = 1, numvar
      sum = 1.d-20
      do k = 1, nfns
        sum = sum + diffns(i, k) * diffns(i, k)
      end do
      parmod(i) = 1.d0 / Sqrt (sum)
      parmod(i) = 1.d0
      do j = 1, i
        sum = 0.d0
        do k = 1, nfns
          sum = sum + diffns(i, k) * diffns(j, k)
        end do
        l = l + 1
        par(l) = sum * parmod(i) * parmod(j)
      end do
    end do
    if (numvar < 9) then
      write (ifiles_8, "(/15X,A,/)") "Parameter Hessian"
      l = 1
      do i = 1, numvar
        ll = l
        l = l + i
        lu = l - 1
        write (ifiles_8, "(8F10.1)") (par(j), j=ll, lu)
      end do
    end if
    call rsp (par, numvar, parmod, vects)
    write (ifiles_8, "(//20X,A,/)") " Eigenvalues of Parameter Hessian"
    do j = 1, numvar
      if (parmod(j) > 1.d7) then
        j2 = j - 1
        go to 1000
      end if
    end do
    j2 = numvar
1000 write (ifiles_8, "(8F10.1)") (parmod(i), i=1, j2)
    if (j2 /= numvar) then
      write (ifiles_8, "(/20X,A)") "(Eigenvalues over 10,000,000 omitted)"
    end if
    write (ifiles_8, "(/)")
    if (parmod(1) < eiglim(1)) then
      lim = 1
      write (ifiles_8, "(20X,A)") " Parameters are ill-defined"
    else if (parmod(1) < eiglim(2)) then
      lim = 2
      write (ifiles_8, "(20X,A)") " Parameters are very poorly defined"
    else if (parmod(1) < eiglim(3)) then
      lim = 3
      write (ifiles_8, "(20X,A)") " Parameters are poorly defined"
    else
      lim = 4
      write (ifiles_8, "(20X,A)") " Parameters are defined by reference data"
    end if
    if (lim == 4) return
    do i = 1, numvar
      if (parmod(i) > eiglim(lim)) go to 1100
    end do
    i = numvar
1100 write (ifiles_8, "(/15X,A,/)") " Parameter Eigenvectors of Low Eigenvalue&
   &s"
    i = Min (i, 8)
    write (ifiles_8, "(12X,8F8.1)") (parmod(j), j=1, i)
    locvar(1,numvar+1) = 0
    locvar(2,numvar+1) = 0
    do j = 1, numvar
      il = (j-1) * numvar + 1
      iu = il + i - 1
      write (ifiles_8, "(3X,A,2X,A,2X,8F8.3)") partyp (locvar(1, j)), elemnt &
     & (locvar(2, j)), (vects(k), k=il, iu)
    end do
end subroutine chkpar
