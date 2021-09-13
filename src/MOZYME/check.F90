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

subroutine check (nvec, nnc, nc, icvec, ic_dim, iorbs, ncvec, cvec, c_dim)
   !***********************************************************************
   !   CHECK corrects very small normalization errors in the LMOs.
   !   If the errors are large, greater than 0.1, then a diagnostic
   !   message is printed, and the run stopped.  In practice, the errors
   !   are either (a) very small, in which case correcting them is
   !   appropriate, or (b) very large, in which case there is a bug in the
   !   program.
   !
   !***********************************************************************
    use molkst_C, only: numat
    use chanel_C, only: iw
    use MOZYME_C, only : ws
   !
   !.. Implicit Declarations ..
    implicit none
   !
   !.. Formal Arguments ..
    integer, intent (in) :: ic_dim, c_dim, nvec
    integer, dimension (ic_dim), intent (in) :: icvec
    integer, dimension (nvec+1), intent (in) :: nc, ncvec, nnc
    integer, dimension (numat), intent (in) :: iorbs
    double precision, dimension (c_dim), intent (inout) :: cvec
   !
   !.. Local Scalars ..
    integer :: i, j, k, l, m, mm, n
    double precision :: error, sum
   !
   !.. Intrinsic Functions ..
    intrinsic Abs, Sqrt
   !
   ! ... Executable Statements ...
   !
    error = 0.d0
   !
   !   RENORMALIZE L.M.O.s
   !
    do i = 1, nvec
      m = nnc(i)
      sum = 0.d0
      n = 0
      do j = 1, nc(i)
        m = m + 1
        k = icvec(m)
        do mm = 1, iorbs(k)
          n = n + 1
          sum = sum + cvec(ncvec(i)+n) ** 2
        end do
      end do
      error = error + Abs (1.d0-sum)
      ws(i) = sum
    end do
    do i = 1, nvec
      m = nnc(i)
      sum = 1.d0 / Sqrt (ws(i))
      n = 0
      do j = 1, nc(i)
        m = m + 1
        k = icvec(m)
        do mm = 1, iorbs(k)
          n = n + 1
          cvec(ncvec(i)+n) = cvec(ncvec(i)+n) * sum
        end do
      end do
    end do
    if (error <= 0.1d0) return
   !
   !   SEVERE ERROR IN LMO's.  QUIT THE JOB
   !
    call mopend(" ERROR DETECTED IN SUBROUTINE CHECK")
    write (iw, "(/,A,F12.6,/)") " MAGNITUDE OF ERROR:", error
    write (iw,*) " NUMBER OF VECTORS:", nvec
    write (iw, "(/,A,/)") " L.M.O.s WERE NORMALIZED TO"
    write (iw, "(5F16.12)") (ws(i), i=1, nvec)
    write (iw,*)
    do i = 1, nvec
      if (Abs (ws(i)-1.d0) > 0.1d0) then
        write (iw, "(/,A,I6,A,F12.6/)") "  DETAILS OF LMO", i, " N:", ws (i)
        write (iw,*) " NUMBER OF ATOMS IN LMO:         ", nc (i)
        write (iw,*) " STARTING ADDRESS OF LMO ATOMS:  ", nnc (i)
        write (iw,*) " STARTING ADDRESS OF LMO COEFFS.:", ncvec (i)
        write (iw,*) " ATOM NUMBERS OF ATOMS IN LMO"
        write (iw, "(10I7)") (icvec(j), j=nnc(i)+1, nnc(i)+nc(i))
        write (iw,*) " NNC(I+1):", nnc (i+1)
        write (iw,*) " NCVEC(I+1):", ncvec (i+1)
        l = 0
        write (iw,*) " ORBITAL COEFFICIENTS"
        sum = Sqrt (ws(i))
        do j = 1, nc(i)
          k = icvec(nnc(i)+j)
          write (iw, "(2I6,9F12.6)") j, k, (cvec(m)*sum, m=ncvec(i)+l+1, &
         & ncvec(i)+l+iorbs(k))
          l = l + iorbs(k)
        end do
        write (iw,*) " NUMBER OF ORBITALS IN LMO:", l
      end if
    end do
end subroutine check
