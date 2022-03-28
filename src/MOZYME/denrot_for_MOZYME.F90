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

subroutine denrot_for_MOZYME ()
    use molkst_C, only: numat, maxtxt, id, l1u, l2u, &
       & l3u
!
    use MOZYME_C, only : iorbs
    use chanel_C, only: iw
    use elemts_C, only: elemnt
    use common_arrays_C, only: tvec, coord, p, txtatm, nat
    implicit none
    character (len=40) :: blank
    logical :: first
    integer :: i, i1, ii, ij, imax, io, ipq, is_loc, j, j1, jj, jo, jpq, js, &
   & ks, l, l1, l2, lsp
    double precision :: delx, dely, delz, r, rmin, sum
    integer, dimension (4), save :: isp
    integer, dimension (5, 10), save :: irot
    double precision, dimension (3) :: xjuc
    double precision, dimension (4) :: sigpi
    double precision, dimension (4, 4) :: arot, pab, vect
    double precision, dimension (3, 5, 5) :: c
    integer, external :: ijbo
    data irot / 1, 1, 1, 3, 3, 2, 2, 2, 4, 3, 3, 2, 2, 2, 3, 4, 2, 2, 3, 3, 2, &
   & 3, 2, 4, 2, 3, 3, 2, 2, 2, 4, 3, 2, 3, 2, 2, 4, 2, 4, 4, 3, 4, 2, 2, 4, &
   & 4, 4, 2, 3, 4 /
    data isp / 1, 2, 3, 3 /
    blank = " "
    if (maxtxt /= 0) then
      write (iw, "(2A)") blank(:maxtxt/2 + 3)//"ATOM"//blank(:maxtxt + 4)//"ATOM" //&
      & blank(:maxtxt/2 + 7)//"s-s Sig  s-p Sig  p-p Sig  p-p Pi"
    else
      write (iw, "(2A)") "   ATOM    ATOM    S-S (Sig) " // &
                       & "  S-P (Sig)   P-P (Sig)    P-P (PI)"
    end if
    write (iw,*)
    imax = maxtxt + 10
    pab(:,:) = 0.d0
    do i = 1, numat
      io = min(iorbs(i),4)  !  Sometime this subroutine should be expanded to include "d" orbitals
      first = .true.
      ipq = Max (1, io-2)
      do j = 1, numat
        if (i /= j .and. ijbo (i, j) >= 0) then
          ij = ijbo (i, j)
          jo = min(iorbs(j),9)
          jpq = Max (1, jo-2)
          if (id == 0) then
            delx = coord(1, j) - coord(1, i)
            dely = coord(2, j) - coord(2, i)
            delz = coord(3, j) - coord(3, i)
          else
            rmin = 1.d6
            do is_loc = -l1u, l1u
              do js = -l2u, l2u
                do ks = -l3u, l3u
                  do l = 1, 3
                    xjuc(l) = coord(l, j) + tvec(l, 1) * is_loc + tvec(l, &
                   & 2) * js + tvec(l, 3) * ks - coord(l, i)
                  end do
                  r = xjuc(1) ** 2 + xjuc(2) ** 2 + xjuc(3) ** 2
                  if (r < rmin) then
                    rmin = r
                    delx = xjuc(1)
                    dely = xjuc(2)
                    delz = xjuc(3)
                  end if
                end do
              end do
            end do
          end if
          sum = 0.d0
          do ii = 1, io
            do jj = 1, jo
              ij = ij + 1
              sum = sum + p(ij) ** 2
              pab(ii, jj) = p(ij)
            end do
          end do
          if (sum >= 0.1d0) then
            call coe (delx, dely, delz, ipq, jpq, c, r)
            do i1 = 1, 4
              do j1 = 1, 4
                arot(i1, j1) = 0.d0
              end do
            end do
            do i1 = 1, 10
              arot(irot(1, i1), irot(2, i1)) = c(irot(3, i1), irot(4, i1), &
             & irot(5, i1))
            end do
            do i1 = 1, io
              do j1 = 1, jo
                vect(i1, j1) = 0.d0
              end do
            end do
            if (i /= j) then
              ij = Max (io, jo)
              do i1 = 1, io
                do j1 = 1, jo
                  sum = 0.d0
                  do l1 = 1, ij
                    do l2 = 1, ij
                      sum = sum + arot(l1, i1) * pab(l1, l2) * arot(l2, j1)
                    end do
                  end do
                  vect(isp(i1), isp(j1)) = vect(isp(i1), isp(j1)) + sum ** 2
                end do
              end do
            end if
            sigpi(2) = 0.d0
            lsp = 1
               !
               !    S-S SIGMA
               !
            sigpi(1) = vect(1, 1)
            if (jo > 1) then
              lsp = 2
                  !
                  !    S-P SIGMA
                  !
              sigpi(2) = vect(1, 2)
              if (io > 1) then
                lsp = 4
                     !
                     !    P-S SIGMA (CONT)
                     !
                sigpi(2) = sigpi(2) + vect(2, 1)
                     !
                     !    P-P SIGMA
                     !
                sigpi(3) = vect(2, 2)
                     !
                     !    P-P PI
                     !
                sigpi(4) = vect(3, 3)
              end if
            else if (io > 1) then
              lsp = 2
                  !
                  !    P-S SIGMA
                  !
              sigpi(2) = vect(2, 1)
            end if
            if (first) then
              if (maxtxt > 0) then
                write (iw, "(I6,A,I6,A,4F9.5)") i, elemnt (nat(i)) // "(" // &
               & txtatm (i) (:maxtxt) // ")", j, elemnt (nat(j)) // "(" // &
               & txtatm (j) (:maxtxt) // ")", (sigpi(ii), ii=1, lsp)
              else
                write (iw, "(I6,A,I6,A,4F12.5)") i, elemnt (nat(i)), j, elemnt &
               & (nat(j)), (sigpi(ii), ii=1, lsp)
              end if
              first = .false.
            else if (maxtxt > 0) then
              write (iw, "(A,I6,A,4F9.5)") blank (:imax), j, elemnt (nat(j)) &
             & // "(" // txtatm (j) (:maxtxt) // ")", (sigpi(ii), ii=1, lsp)
            else
              write (iw, "(8X,I6,A,4F12.5)") j, elemnt (nat(j)), (sigpi(ii), &
             & ii=1, lsp)
            end if
          end if
        end if
      end do
    end do
end subroutine denrot_for_MOZYME
