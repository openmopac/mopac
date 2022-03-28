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

      subroutine matou1(a, b, ncx, nr, ndim, iflag)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use molkst_C, only : numat, norbs, nalpha, nclose, keywrd
      use common_arrays_C, only : nfirst, nlast, nat
      use chanel_C, only : iw
      use elemts_C, only : elemnt
      use to_screen_C, only : to_a, to_b
      USE symmetry_C, only :  jndex, namo
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: ncx
      integer , intent(inout) :: nr
      integer , intent(in) :: ndim
      integer , intent(in) :: iflag
      double precision , intent(in) :: a(nr,nr)
      double precision , intent(in) :: b(ndim)

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(nr) :: natom
      integer :: nc, nsave, nfix, i, jlo, jhi, l, j, ka, kc, kb, la, lc, lb, &
        nup, ndown
      logical :: allprt
      character , dimension(9) :: atorbs*2
      character , dimension(nr) :: itext*2, jtext*2
      character, dimension(3) :: xyz*2

      save atorbs, xyz
!-----------------------------------------------
!**********************************************************************
!
!      MATOUT PRINTS A SQUARE MATRIX OF EIGENVECTORS AND EIGENVALUES
!
!    ON INPUT A CONTAINS THE MATRIX TO BE PRINTED.
!             B CONTAINS THE EIGENVALUES.
!             NC NUMBER OF MOLECULAR ORBITALS TO BE PRINTED.
!             NR IS THE SIZE OF THE SQUARE ARRAY TO BE PRINTED.
!             NDIM IS THE ACTUAL SIZE OF THE SQUARE ARRAY "A".
!             NFIRST AND NLAST CONTAIN ATOM ORBITAL COUNTERS.
!             NAT = ARRAY OF ATOMIC NUMBERS OF ATOMS.
!
!
!     OUTPUT TYPE (ROW LABELING)
!       IFLAG=1 : ORBITALS
!       IFLAG=2 : ORBITALS + SYMMETRY-DESIGNATORS
!       IFLAG=3 : ATOMS
!       IFLAG=4 : NUMBERS ONLY
!       IFLAG=5 : VIBRATIONS + SYMMETRY-DESIGNATIONS
!
!
!***********************************************************************
      data xyz/ ' x', ' y', ' z'/
      data atorbs/ 'S ', 'Px', 'Py', 'Pz', 'x2', 'xz', 'z2', 'yz', 'xy'/
!      -------------------------------------------------
      nfix = 0
      allprt = index(keywrd,'ALLVEC') /= 0
      i = index(keywrd," VECTORS")
      nup = 7
      ndown = 9
      if (i > 0) then
        do j = i + 1, i + 20
          if (keywrd(j:j) == " ") exit
        end do
        if (index(keywrd(i:j), ")") > 0) then
!
!  User has specified number of occupied and virtual M.O.s to be printed
!  Format:  VECTORS=(12,34)
!
        do i = i + 1, j
            if (keywrd(i:i) >= "0" .and. keywrd(i:i) <= "9") exit
          end do
          ndown = ichar(keywrd(i:i)) - ichar("0")
          do i = i + 1, j
            if (keywrd(i:i) < "0" .or. keywrd(i:i) > "9") exit
            ndown = ndown*10 + ichar(keywrd(i:i)) - ichar("0")
          end do
          do i = i, j
            if (keywrd(i:i) >= "0" .and. keywrd(i:i) <= "9") exit
          end do
          nup = ichar(keywrd(i:i)) - ichar("0")
          do i = i + 1, j
            if (keywrd(i:i) < "0" .or. keywrd(i:i) > "9") exit
            nup = nup*10 + ichar(keywrd(i:i)) - ichar("0")
          end do
        end if
        ndown = ndown - 1
        if (nup == -16) then
          nup = (ndown + 1)/2
          ndown = ndown - nup
        end if
      end if
      nc = ncx
      if (iflag > 2 .and. iflag /= 5) go to 50
      if ( .not. allprt) then
        nsave = ncx
        nfix = max(nalpha,nclose)
        if (iflag == 2 .and. nc > 16) nc = nfix + nup
        nc = min0(nsave,nc)
      end if
      if (numat == 0) go to 50
      if (nlast(numat) /= nr) go to 50
      do i = 1, numat
        jlo = nfirst(i)
        jhi = nlast(i)
        l = nat(i)
        if (iflag <= 2) then
          itext(jlo:jhi) = atorbs(:jhi-jlo+1)
          jtext(jlo:jhi) = elemnt(l)
          natom(jlo:jhi) = i
        else
          jhi = 3*(i - 1)
          itext(1+jhi:3+jhi) = xyz
          jtext(1+jhi:3+jhi) = elemnt(l)
          natom(1+jhi:3+jhi) = i
        end if
      end do
      go to 70
   50 continue
      nr = abs(nr)
      if (iflag == 3) then
        do i = 1, nr
          itext(i) = '  '
          jtext(i) = elemnt(nat(i))
          natom(i) = i
        end do
      else
        do i = 1, nr
          itext(i) = '  '
          jtext(i) = '  '
          natom(i) = i
        end do
      end if
   70 continue
      ka = 1
      kc = 8
      if (.not.allprt) then
        if (iflag == 2 .and. norbs > 16) ka = nfix - ndown
        ka = max0(1,ka)
        if (iflag == 2 .and. norbs > 16) kc = ka + 7
      end if
   90 continue
      kb = min0(kc,nc)
      write (iw, 130) (i,i=ka,kb)
      if (iflag == 2 .or. iflag == 5) write (iw, 180) (jndex(i),namo(i),i=ka,kb)
      if (b(1) /= 0.D0) then
        if (iflag == 5) then
          write (iw, 150) (b(i),i=ka,kb)
        else
          write (iw, 140) (b(i),i=ka,kb)
        end if
      end if
      write (iw, 160)
      la = 1
      lc = 40
  100 continue
      lb = min0(lc,nr)
      do i = la, lb
        if (itext(i) == ' S') write (iw, 160)
        write (iw, 170) itext(i), jtext(i), natom(i), (a(i,j),j=ka,kb)
      end do
      if (lb == nr) go to 120
      la = lc + 1
      lc = lc + 40
      go to 100
  120 continue
      if (kb == nc) then
        if (allocated (to_a)) deallocate(to_a)
        if (allocated (to_b)) deallocate(to_b)
        allocate(to_a(nr,nr), to_b(ndim))
        to_a = a
        to_b = b
        deallocate(to_a, to_b)
        return
      end if
      ka = kc + 1
      kc = kc + 8
      go to 90

  130 format(/,/,2x,' Root No.',i8,11i10)
  140 format(/,12x,10f10.3)
  150 format(/,12x,10f10.1)
  160 format('  ')
  170 format(' ',2(1x,a2),i5,f10.4,10f10.4)
  180 format(/,13x,10(i5,1x,a4))
      end subroutine matou1
