! Molecular Orbital PACkage (MOPAC)
! Copyright 2021 Virginia Polytechnic Institute and State University
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!    http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

  subroutine lyse
!
!  Before a residue sequence can be worked out, all inter-protein bonds
!  must first be broken.  These bonds are of type N-S, O-S, and S-S
!
    use common_arrays_C, only : nat, nbonds, ibonds
    use molkst_C, only : numat
    implicit none
    integer :: i, j, k, l, ll, m, n
    logical :: special
    double precision :: r_min, r
    double precision, external :: distance
 !
 !   FIRST, GET RID OF H ATOMS AND ATOMS ALREADY DEFINED
 !
    do i = 1, numat
      if (nat(i) /= 1) then
        if (nat(i) == 6) then
          l = 0
          do k = 1, nbonds(i)
            l = l + 1
            ibonds(l, i) = ibonds(k, i)
          end do
        else
          special = .false.
          if (nat(i) == 8) then
!
!  Atom is oxygen - check to see if it's attached to a sulfate or phosphate.
!  If so, do NOT lyse the bond.
!
            do k = 1, nbonds(i)
              if (nat(ibonds(k, i)) == 16 .or. nat(ibonds(k, i)) == 15) then
                if (nbonds(ibonds(k,i)) == 4) special = .true.
              end if
            end do
          end if
          if (special) then
            l = nbonds(i)
            cycle
          end if
!
! EXCLUDE ALL N-S, O-S, and S-S BONDS
!
          l = 0
          do k = 1, nbonds(i)
            if (nat(ibonds(k, i)) /= 16) then
              l = l + 1
              ibonds(l, i) = ibonds(k, i)
            end if
          end do
          ll = 0
          do k = 1, l
            j = ibonds(k,i)
            if(nat(j) == 7 .or. nat(j) == 6) ll = 1
          end do
          if (ll == 0) cycle
!
! EXCLUDE ALL BONDS TO NON-PROTEIN ATOMS
!
          ll = l
          l = 0
          do k = 1, ll
            j = ibonds(k, i)
            if (nat(j) == 1 .or. nat(j) == 6 .or. nat(j) == 7 .or. nat(j) == 8 &
           & .or. nat(j) == 16) then
              l = l + 1
              ibonds(l, i) = ibonds(k, i)
            else
!
!  Also remove bond from non-protein atom
!
              m = 0
              do n = 1, nbonds(j)
                if (ibonds(n,j) /= i) then
                  m = m + 1
                  ibonds(m,j) = ibonds(n,j)
                end if
              end do
              nbonds(j) = nbonds(j) - 1
            end if
          end do
        end if
        nbonds(i) = l
      else
        if (nbonds(i) > 1) then
!
! Break all bonds except the shortest
!
          r_min = 1000.d0
          k = 0
          do j = 1, nbonds(i)
            r = distance(ibonds(j,i), i)
!
!  Chose the shortest bond, but not shorter than a normal X-H
!
            if (r < r_min .and. r > 0.95d0) then
              r_min = r
              k = ibonds(j,i)
            end if
          end do
          if ( k > 0) then
!
! Shortest bond is i to k
! Break all bonds to atom i
!
            do j = 1, nbonds(i)
              l = ibonds(j,i)
              do m = 1, nbonds(l)
                if (ibonds(m,l) == i) exit
              end do
              nbonds(l) = nbonds(l) - 1
              do n = m, nbonds(l)
                ibonds(n,l) = ibonds(n + 1,l)
              end do
            end do
!
!  Create bond between atoms i and k
!
            nbonds(i) = 1
            ibonds(1,i) = k
            nbonds(k) = nbonds(k) + 1
            ibonds(nbonds(k),k) = i
          end if
        end if
      end if
    end do
  end subroutine lyse
