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

      subroutine vecprt(a, numm)
      use molkst_C, only : numat, mozyme
      use common_arrays_C, only : nfirst, nlast, nat
      use chanel_C, only : iw
      use elemts_C, only : elemnt
      implicit none
      integer , intent(in) :: numm
      double precision , intent(inout) :: a(*)
!
! Local
!
      integer , dimension(:), allocatable :: natom
      integer :: i, jlo, jhi, l, k, numb, limit, kk, na, ll, m, ma, n
      double precision :: sumax, fact
      character :: line(21)*6, atorbs(9)*2, fmt*31
      character, dimension(:), allocatable :: itext*2, jtext*2
      save atorbs
!**********************************************************************
!
!  VECPRT PRINTS A LOWER-HALF TRIANGLE OF A SQUARE MATRIX, THE
!         LOWER-HALF TRIANGLE BEING STORED IN PACKED FORM IN THE
!         ARRAY "A"
!
! ON INPUT:
!      A      = ARRAY TO BE PRINTED
!      NUMM   = SIZE OF ARRAY TO BE PRINTED
!(REF) NUMAT  = NUMBER OF ATOMS IN THE MOLECULE (THIS IS NEEDED TO
!               DECIDE IF AN ATOMIC ARRAY OR ATOMIC ORBITAL ARRAY IS
!               TO BE PRINTED
!(REF) NAT    = LIST OF ATOMIC NUMBERS
!(REF) NFIRST = LIST OF ORBITAL COUNTERS
!(REF) NLAST  = LIST OF ORBITAL COUNTERS
!
!  NONE OF THE ARGUMENTS ARE ALTERED BY THE CALL OF VECPRT
!
!*********************************************************************
      data atorbs/ ' S', 'PX', 'PY', 'PZ', 'X2', 'XZ', 'Z2', 'YZ', 'XY'/
      l = 0
      if (mozyme) then
        call vecprt_for_MOZYME(a,numm)
        return
      end if
      numb = abs(numm)
      allocate (natom(numb), itext(numb), jtext(numb))
      sumax = 1.D0
      do i = 1, numm
        sumax = max(abs(a((i*(i+1))/2)), sumax)
      end do
      i = int(log10(sumax))
      if (i == 1 .or. i == 2) i = 0
      fact = 10.D0**(-i)
      if (abs(fact - 1.D0) > 0.001D0) then
        write (iw, '(A,F12.1)') 'Diagonal Terms should be Multiplied by', 1.D0/fact
        do i = 1, numm
          a((i*(i+1))/2) = a((i*(i+1))/2)*fact
        end do
      end if
      if (numat /= 0 .and. numat == numm) then
!
!    PRINT OVER ATOM COUNT
!
        do i = 1, numat
          itext(i) = '  '
          jtext(i) = elemnt(nat(i))
          natom(i) = i
        end do
      else
        if (numat /= 0 .and. nlast(numat) == numm) then
          do i = 1, numat
            jlo = nfirst(i)
            jhi = nlast(i)
            l = nat(i)
            k = 0
            itext(jlo:jhi) = atorbs(:jhi-jlo+1)
            jtext(jlo:jhi) = elemnt(l)
            natom(jlo:jhi) = i
          end do
        else
          do i = 1, numb
            itext(i) = '  '
            jtext(i) = '  '
            natom(i) = i
          end do
        end if
      end if
      line = '------'
      limit = (numb*(numb + 1))/2
      kk = 8
      na = 1
      if (numb > 999) then
        fmt = "(/,/,12x,10(1x,a2,1x,a2,i4,1x))"
      else
        fmt = "(/,/,13x,10(1x,a2,1x,a2,i3,2x))"
      end if
   80 continue
      ll = 0
      m = min0(numb + 1 - na, 6)
      ma = 2*m + 1
      m = na + m - 1
      write (iw, fmt) (itext(i), jtext(i), natom(i), i = na, m)
      write (iw, "(' ',21a6)") (line(k),k=1,ma)
      do i = na, numb
        ll = ll + 1
        k = (i*(i - 1))/2
        l = min0(k + m,k + i)
        k = k + na
        if (kk + ll > 50) then
          write (iw, fmt) (itext(n), jtext(n), natom(n) ,n = na, m)
          write (iw, "(' ',21a6)") (line(n),n=1,ma)
          kk = 4
          ll = 0
        end if
        write (iw, "(' ', a2, 1x, a2, i5, 10f11.6)") itext(i), jtext(i), natom(i), (a(n), n = k, l)
      end do
      if (l >= limit) go to 110
      kk = kk + ll + 4
      na = m + 1
      if (kk + numb + 1 - na <= 50) go to 80
      kk = 4
      go to 80
  110 continue
      if (abs(fact - 1.D0) > 0.001D0) then
        do i = 1, numm
          a((i*(i+1))/2) = a((i*(i+1))/2)/fact
        end do
      end if
      return
      end subroutine vecprt
