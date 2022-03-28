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

  subroutine fillij (count)
      use molkst_C, only: numat, natoms, cutofp, id, n2elec, l1u, l2u, l3u, keywrd, mpack, &
        ispd, line
      use common_arrays_C, only : tvec, coord
      use chanel_C, only: iw
      use MOZYME_C, only : cutofs, direct, semidr, &
        nijbo, lijbo, iijj, iij, ij_dim, ijall, numij, morb, iorbs
      use overlaps_C, only : cutof1, cutof2
!
      implicit none
      !
      !***********************************************************************
      !
      ! fillij creates either the array ijbo, if there is enough memory, or
      ! the arrays iijj, ijall, and iij, if there is not enough memory.
      !
      !***********************************************************************
      logical :: count
!
      integer :: i, ii, iloop, io, ip, j, jloop, jo, jp, &
           & kp, lp, ix
      double precision :: r, rmin, rr, x1, x2, x3
      save :: ix
      logical :: first
      double precision, dimension (3) :: xj
      double precision, external :: reada
!
      if (.not. allocated(nijbo)) then
        ix = 0
        n2elec = 0
        mpack = 0
        first = .true.
        morb = 4
        if (ispd > 0) morb = 9
      else
!
!   Array nijbo was already allocated, therefore it was filled.  Do not reset mpack1, instead add any new interactions
!   to array nijbo
!
        first = .false.
      end if
!
      if (.not. count) then
        ! Allocate either nijbo or iijj and ijall
        if ( .not. allocated(nijbo)) then
          i = 1
          if (lijbo) then
!
!  Try to allocate memory for a simple square array
!
            allocate (nijbo(numat, numat), stat = i)
            if (i == 0) nijbo(:, :) = 0
          end if
          if (i /= 0) then

!
!  Array nijbo could not be created.  Therefore set lijbo false and
!  create smaller, but more CPU intensive, arrays
!
            lijbo = .false.
            allocate (iijj(ij_dim), ijall(ij_dim), iij(natoms), numij(natoms), &
               & stat = i)
            if (i /= 0) then
              call mopend("There is not enough memory to create the ijbo-type arrays")
              return
            end if
            iijj(:) = 0
            ijall(:) = 0
            iij(:) = 0
            numij(:) = 0
          end if
        end if
      else
        lijbo = .true.
      end if
      !
      !   Set CUTOF values  CUTOF2 = First cutoff (NDDO-dipolar)
      !                     CUTOF1 = Second cutoff (dipolar-point charge)
      !

      cutof1 = 10.d0 ** 2
      if (cutofp < 100.d0) cutof1 = cutofp ** 2 + 1.d-4
      if (numat < 30)      cutof1 = 1.d6 + 10
      cutof2 = 9.9d0 ** 2
      if (cutofp < 100.d0) cutof2 = cutofp ** 2
      if (numat < 30)      cutof2 = 1.d6
      line = trim(keywrd)
      if (index(line," GEO_DAT") /= 0) then
        i = index(line," GEO_DAT") + 9
        j = index(line(i + 10:),'" ') + i + 9
        line(i:j) = " "
      end if

      i = Index (line, " CUTOFF=")
      if (i /= 0) then
        cutof1 = reada (line, i+8)
        cutof2 = (cutof1 - 1.d-1)**2
        cutof1 = cutof1**2
      else
        i = Index (line, " CUTOF1=")
        if (i /= 0) cutof1 = reada (keywrd, i+8) ** 2
        i = Index (line, " CUTOF2=")
        if (i /= 0) cutof2 = reada (keywrd, i+8) ** 2
      end if
      !
      i = Index (line, " CUTOFS=")
      if (i /= 0) then
        cutofs = reada (line, i+8) ** 2
      else
        cutofs = 7.d0 ** 2
      end if
      !
      !   Check that CUTOF1 is greater than CUTOF2
      !
      if (cutof1 < cutof2) then
        cutof1 = cutof2 + 1.d-4
        write (iw, "(/,A,F12.4,/)") " CUTOF1 WAS SET SMALLER THAN CUTOF2," // &
             & " RESET TO", Sqrt (cutof1)
      end if
      !#aab
      if ((Index (keywrd, " NODIRECT") /= 0) .or. (id /= 0) ) then
        direct = .false.
        semidr = .false.
      else if (Index (keywrd, " SEMIDIRECT") /= 0) then
        direct = .false.
        semidr = .true.
      else
        direct = .true.
        semidr = .true.
      end if
      !#aab - end
      !
      rmin = 100.d0
      !
      !     Identify atom pairs involved in calculation.  This is done based
      !     on interatomic distance.
      !
      i = 0
      do iloop = 1, numat
          io = iorbs(iloop)
          i = i + 1
          x1 = coord(1, i)
          x2 = coord(2, i)
          x3 = coord(3, i)
          ii = (io*(io+1)) / 2
          if (first) n2elec = n2elec + ii * ii ! for one-centre two-electron integrals
          if (.not. lijbo .and. .not. count) then
            iij(i) = ix + 1
          end if
          j = 0
          do jloop = 1, iloop - 1
              jo = iorbs(jloop)
              j = j + 1
              if (id == 0) then
                r = (x1-coord(1, j)) ** 2 + (x2-coord(2, j)) ** 2 &
                     & + (x3-coord(3, j)) ** 2
              else
                r = 1.d10
                do ip = -l1u, l1u
                  do jp = -l2u, l2u
                    do kp = -l3u, l3u
                      do lp = 1, 3
                        xj(lp) = coord(lp, j) + tvec(lp, 1) * ip &
                             & + tvec(lp, 2) * jp + tvec(lp, 3) * kp
                      end do
                      rr = (x1-xj(1)) ** 2 + (x2-xj(2)) ** 2 + (x3-xj(3)) ** 2
                      r = Min (r, rr)
                    end do
                  end do
                end do
              end if

              if (r < cutof2) then
                ix = ix + 1
                if (.not. count) then
                  if (lijbo) then
                    if (first .or. nijbo(i,j) < 0) then
                      nijbo(i, j) = mpack
                      nijbo(j, i) = mpack
                    else
                      mpack = mpack - io * jo  !  prevent mpack from being incremented - the element of
                                                 !  nijbo was already set.
                    end if
                  else
                    iijj(ix) = mpack
                    ijall(ix) = j
                  end if
                end if
                mpack = mpack + io * jo
                !
                if ( .not. direct) then
                  n2elec = n2elec + (jo*(jo+1)) / 2 * ii
                end if
                if (r < rmin) then
                  rmin = r
                end if
              else if (r < cutof1) then
                !
                !   USE CHARGES AND DIPOLES
                !
                if ( .not. semidr) then
                  !
                  if (io > 1) then
                    if (jo > 1) then
                      n2elec = n2elec + 7
                    else
                      n2elec = n2elec + 4
                    end if
                  else if (jo > 1) then
                    n2elec = n2elec + 4
                  else
                    n2elec = n2elec + 1
                  end if
                end if !#aab - end
                !
                if (lijbo .and. (.not. count)) then
                  if (first .or. nijbo(i,j) == -1) then
                    nijbo(i, j) = -2
                    nijbo(j, i) = -2
                  end if
                end if
                !#aab
              else
                !
                !  USE POINT CHARGE APPROXIMATION
                !
                if ( .not. semidr) then
                  n2elec = n2elec + 1
                end if
                !
                if (lijbo .and. (.not. count)) then
                  if (first) then
                    nijbo(i, j) = -1
                    nijbo(j, i) = -1
                  end if
                end if
                !#aab
            end if
          end do
          !
          ix = ix + 1
          if (.not. count) then
            if (lijbo) then
              if (first) then
                nijbo(i, i) = mpack
              else
                mpack = mpack - (io*(io+1)) / 2  !  prevent mpack from being incremented - the element of
                                                   !  nijbo was already set.
              end if
            else
              iijj(ix) = mpack
              ijall(ix) = i
              numij(i) = ix
            end if
          end if
          if (id /= 0) then
            n2elec = n2elec + ((io*(io+1))/2) ** 2
          end if
          mpack = mpack + (io*(io+1)) / 2
      end do
      !
      if (count) then
        ij_dim = ix
      end if
      !
      !  Leave room for dummy atom integrals and for gradients
      !
      if (n2elec < 2025) then
        n2elec = 2025
      end if
      if (first) then
        n2elec = n2elec + 10
        if (direct .and. ispd == 0) n2elec = n2elec + 100
        if (direct .and. ispd /= 0) n2elec = n2elec + 2025
      end if
      if (id /= 0) then
        n2elec = n2elec * 2
      end if
    end subroutine fillij
