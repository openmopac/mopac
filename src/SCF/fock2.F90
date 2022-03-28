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

    subroutine fock2(f, ptot, p, w, wj, wk, numat, nfirst, nlast, mode)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use molkst_C, only : numcal, norbs, mpack, n2elec, id, numat_ref => numat
      use cosmo_C, only : useps

!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: numat, mode
      integer  :: nfirst(numat_ref)
      integer  :: nlast(numat_ref)
      double precision  :: f(mpack)
      double precision , intent(in) :: ptot(mpack)
      double precision  :: p(mpack)
      double precision  :: w(n2elec)
      double precision , intent(in) :: wj(n2elec)
      double precision , intent(in) :: wk(n2elec)


!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(:), allocatable :: ifact, i1fact
      double precision, dimension(:,:), allocatable  :: ptot2
      integer :: ione, icalcn
      integer , dimension(256) :: jindex
      integer :: i, m, j, ij, ji, k, ik, l, kl, lk, jl, i1, ia, ib, jk, kj, kk&
        , ii, iminus, jj, ja, jb, ll, j1, ka, kb, kc, il
      double precision, dimension(81) :: pk
      double precision, dimension(16) :: pja, pjb
      double precision :: sumdia, sumoff, sum, aa, bb, aj, ak, a
      logical :: lid, deriv

      save ifact, i1fact, ione, lid, icalcn, jindex, ptot2
!-----------------------------------------------
!***********************************************************************
!
! FOCK2 FORMS THE TWO-ELECTRON TWO-CENTER REPULSION PART OF THE FOCK
! MATRIX
! ON INPUT  PTOT = TOTAL DENSITY MATRIX.
!           P    = ALPHA OR BETA DENSITY MATRIX.
!           W    = TWO-ELECTRON INTEGRAL MATRIX.
!
!  ON OUTPUT F   = PARTIAL FOCK MATRIX
!***********************************************************************
      data icalcn/ 0/
      deriv = (numat < 0)
      numat = abs(numat)
      if (numat == 0) then
        if (allocated(ptot2))  deallocate(ptot2)
        if (allocated(ifact))  deallocate(ifact)
        if (allocated(i1fact)) deallocate(i1fact)
         return
      end if
      if (icalcn /= numcal) then
        if (allocated(ptot2))  deallocate(ptot2)
        if (allocated(ifact))  deallocate(ifact)
        if (allocated(i1fact)) deallocate(i1fact)
        allocate(ptot2(numat,81), ifact(3 + norbs), i1fact(3 + norbs))
        icalcn = numcal
!
!   SET UP ARRAY OF LOWER HALF TRIANGLE INDICES (PASCAL'S TRIANGLE)
!
        do i = 1, norbs
          ifact(i) = (i*(i - 1))/2
          i1fact(i) = ifact(i) + i
        end do
!
!   SET UP GATHER-SCATTER TYPE ARRAYS FOR USE WITH TWO-ELECTRON
!   INTEGRALS.  JINDEX ARE THE INDICES OF THE J-INTEGRALS FOR ATOM I
!   INTEGRALS.  JJNDEX ARE THE INDICES OF THE J-INTEGRALS FOR ATOM J
!               KINDEX ARE THE INDICES OF THE K-INTEGRALS
!
        m = 0
        do i = 1, 4
          do j = 1, 4
            ij = min(i,j)
            ji = i + j - ij
            do k = 1, 4
              ik = min(i,k)
              do l = 1, 4
                m = m + 1
                kl = min(k,l)
                lk = k + l - kl
                jl = min(j,l)
                jindex(m) = (ifact(ji)+ij)*10 + ifact(lk) + kl - 10
              end do
            end do
          end do
        end do
        l = 0
        do i = 1, 4
          i1 = (i - 1)*4
          if (i > 0) then
            l = i + l
          end if
        end do
        lid = id == 0
        ione = 1
        if (id /= 0) ione = 0
!
!      END OF INITIALIZATION
!
      end if
!
!     START OF MNDO, AM1, OR PM3 OPTION
!
      l = 0
      do i = 1, numat
        ia = nfirst(i)
        ib = nlast(i)
        m = 0
        do j = ia, ib
          do k = ia, ib
            m = m + 1
            jk = min(j,k)
            kj = k + j - jk
            jk = jk + (kj*(kj - 1))/2
            ptot2(i,m) = ptot(jk)
          end do
        end do
      end do
      kk = 0
      do ii = 1, numat
        ia = nfirst(ii)
        ib = nlast(ii)
!
!  IF NUMAT=2 THEN WE ARE IN A DERIVATIVE IN A SOLID STATE OR IN A MOLECULE CALCULATION
!
        if (deriv) then
          iminus = ii - 1
        else
          iminus = ii - ione
        end if
        do jj = 1, iminus
          ja = nfirst(jj)
          jb = nlast(jj)
          if (lid) then    !  System does not use periodic boundary conditions
            if (ib - ia >= 6 .or. jb - ja >= 6) then
              call fockdorbs(ia, ib, ja, jb, f, p, ptot, w, kk, ifact)

            else if (ib - ia >= 3 .and. jb - ja >= 3) then
!
!                         HEAVY-ATOM  - HEAVY-ATOM
!
!   EXTRACT COULOMB TERMS
!
              pja = ptot2(ii,:16)
              pjb = ptot2(jj,:16)
!
!  COULOMB TERMS
!
              call jab (ia, ja, pja, pjb, w(kk+1), f)
!
!  EXCHANGE TERMS
!
!
!  EXTRACT INTERSECTION OF ATOMS II AND JJ IN THE SPIN DENSITY MATRIX
!
!  The following loop had been written in a more compact style, but
!  this had caused problems with the INTEL compiler in RELEASE mode when QuickWin was used
!
              l = 1
              do i = ia, ib
                i1 = ifact(i) + ja
                do j = l, l + 3
                  pk(j) = p(i1)
                  i1 = i1 + 1
                end do
                l = 4 + l
              end do
              call kab (ia, ja, pk, w(kk+1), f)
              kk = kk + 100
            else if (ib - ia >= 3 .and. ja == jb) then
!
!                         LIGHT-ATOM  - HEAVY-ATOM
!
!
!   COULOMB TERMS
!
              sumdia = 0.D0
              sumoff = 0.D0
              ll = i1fact(ja)
              k = 0
              do i = 0, 3
                j1 = ifact(ia+i) + ia - 1
                if (i > 0) then
                  do j = 1, i
                    f(j+j1) = f(j+j1) + ptot(ll)*w(j+kk+k)
                    sumoff = sumoff + ptot(j+j1)*w(j+kk+k)
                  end do
                  k = i + k
                  j1 = i + j1
                end if
                j1 = j1 + 1
                k = k + 1
                f(j1) = f(j1) + ptot(ll)*w(kk+k)
                sumdia = sumdia + ptot(j1)*w(kk+k)
              end do
              f(ll) = f(ll) + sumoff*2.D0 + sumdia
!
!  EXCHANGE TERMS
!
!
!  EXTRACT INTERSECTION OF ATOMS II AND JJ IN THE SPIN DENSITY MATRIX
!
              k = 0
              do i = ia, ib
                i1 = ifact(i) + ja
                sum = 0.D0
                do j = 1, ib - ia + 1
                  sum = sum + p(ifact(j-1+ia)+ja)*w(kk+jindex(j+k))
                end do
                k = ib - ia + 1 + k
                f(i1) = f(i1) - sum
              end do
              kk = kk + 10
            else if (jb - ja >= 3 .and. ia == ib) then
!
!                         HEAVY-ATOM - LIGHT-ATOM
!
!
!   COULOMB TERMS
!
              sumdia = 0.D0
              sumoff = 0.D0
              ll = i1fact(ia)
              k = 0
              do i = 0, 3
                j1 = ifact(ja+i) + ja - 1
                if (i > 0) then
                  do j = 1, i
                    f(j+j1) = f(j+j1) + ptot(ll)*w(j+kk+k)
                    sumoff = sumoff + ptot(j+j1)*w(j+kk+k)
                  end do
                  k = i + k
                  j1 = i + j1
                end if
                j1 = j1 + 1
                k = k + 1
                f(j1) = f(j1) + ptot(ll)*w(kk+k)
                sumdia = sumdia + ptot(j1)*w(kk+k)
              end do
              f(ll) = f(ll) + sumoff*2.D0 + sumdia
!
!  EXCHANGE TERMS
!
!
!  EXTRACT INTERSECTION OF ATOMS II AND JJ IN THE SPIN DENSITY MATRIX
!
              k = ifact(ia) + ja
              j = 0
              do i = k, k + 3
                sum = 0.D0
                do l = 1, 4
                  sum = sum + p(l-1+k)*w(kk+jindex(l+j))
                end do
                j = 4 + j
                f(i) = f(i) - sum
              end do
              kk = kk + 10
            else if (jb == ja .and. ia == ib) then
!
!                         LIGHT-ATOM - LIGHT-ATOM
!
              i1 = i1fact(ia)
              j1 = i1fact(ja)
              ij = i1 + ja - ia
              f(i1) = f(i1) + ptot(j1)*w(kk+1)
              f(j1) = f(j1) + ptot(i1)*w(kk+1)
              f(ij) = f(ij) - p(ij)*w(kk+1)
              kk = kk + 1
            end if
          else
            do i = ia, ib
              ka = ifact(i)
              do j = ia, i
                kb = ifact(j)
                ij = ka + j
                aa = 2.0D00
                if (i == j) aa = 1.0D00
                do k = ja, jb
                  kc = ifact(k)
                  if (i >= k) then
                    ik = ka + k
                  else
                    ik = 0
                  end if
                  if (j >= k) then
                    jk = kb + k
                  else
                    jk = 0
                  end if
                  do l = ja, k
                    if (i >= l) then
                      il = ka + l
                    else
                      il = 0
                    end if
                    if (j >= l) then
                      jl = kb + l
                    else
                      jl = 0
                    end if
                    kl = kc + l
                    bb = 2.0D00
                    if (k == l) bb = 1.0D00
                    kk = kk + 1
                    aj = wj(kk)
                    ak = wk(kk)
!
!     A  IS THE REPULSION INTEGRAL (I,J/K,L) WHERE ORBITALS I AND J ARE
!     ON ATOM II, AND ORBITALS K AND L ARE ON ATOM JJ.
!     AA AND BB ARE CORRECTION FACTORS SINCE
!     (I,J/K,L)=(J,I/K,L)=(I,J/L,K)=(J,I/L,K)
!     IJ IS THE LOCATION OF THE MATRIX ELEMENTS BETWEEN ATOMIB ORBITALS
!     I AND J.  SIMILARLY FOR IK ETC.
!
! THIS FORMS THE TWO-ELECTRON TWO-CENTER REPULSION PART OF THE FOCK
! MATRIX.  THE CODE HERE IS HARD TO FOLLOW, AND IMPOSSIBLE TO MODIFY!,
! BUT IT WORKS,
                    if (kl > ij) cycle
                    if (i == k .and. aa + bb < 2.1D0) then
                      f(ij) = f(ij) + aj*ptot(kl)
                    else
                      f(ij) = f(ij) + bb*aj*ptot(kl)
                      f(kl) = f(kl) + aa*aj*ptot(ij)
                      a = ak*aa*bb*0.25D0
                      if (jl > 0) f(ik) = f(ik) - a*p(jl)
                      if (jk > 0) f(il) = f(il) - a*p(jk)
                      if (jk > 0) f(jk) = f(jk) - a*p(il)
                      if (jl > 0) f(jl) = f(jl) - a*p(ik)
                    end if
                  end do
                end do
              end do
            end do
          end if
        end do
        if (mode == 2) then
          i = ((ib - ia + 1)*(ib - ia + 2))/2
          call fock1(f, ptot, p, mpack, w(kk+1), kk, ia, ib, i)
        end if
      end do
      if (useps) then
        call addfck (f,ptot)
      end if
      return
    end subroutine fock2



    subroutine fockdorbs (ia, ib, ja, jb, f, p, ptot, w, kr, ifact)
      use molkst_C, only: norbs, mpack
      implicit none
      integer, intent (in) :: ia, ib, ja, jb
      integer, intent (inout) :: kr
      integer, dimension (norbs), intent (in) :: ifact
      double precision, dimension (mpack), intent (in) :: p, ptot
      double precision, dimension (*), intent (inout) :: f
      double precision, dimension (*), intent (in) :: w
      integer :: i, ij, ik, il, j, jk, jl, k, ka, kb, kc, kl, l, n1, n2, nn, kref
      double precision :: a, aa, bb
      if (ia > ja) then
        do i = ia, ib
          ka = ifact(i)
          aa = 2.0d00
          do j = ia, i
            if (i == j) then
              aa = 1.0d00
            end if
            kb = ifact(j)
            ij = ka + j
            do k = ja, jb
              kc = ifact(k)
              ik = ka + k
              jk = kb + k
              bb = 2.0d00
              do l = ja, k
                if (k == l) then
                  bb = 1.0d00
                end if
                il = ka + l
                jl = kb + l
                kl = kc + l
                kr = kr + 1
                a = w(kr)
                f(ij) = f(ij) + bb * a * ptot(kl)
                f(kl) = f(kl) + aa * a * ptot(ij)
                a = a * aa * bb * 0.25d0
                f(ik) = f(ik) - a * p(jl)
                f(il) = f(il) - a * p(jk)
                f(jk) = f(jk) - a * p(il)
                f(jl) = f(jl) - a * p(ik)
              end do
            end do
          end do
        end do
      else
        kref = kr
        nn = jb - ja + 1
        nn = (nn*(nn+1)) / 2
        n1 = 0
        do i = ja, jb
          ka = ifact(i)
          aa = 2.0d00
          do j = ja, i
            n1 = n1 + 1
            if (i == j) then
              aa = 1.0d00
            end if
            kb = ifact(j)
            ij = ka + j
            n2 = 0
            do k = ia, ib
              kc = ifact(k)
              ik = ka + k
              jk = kb + k
              bb = 2.0d00
              do l = ia, k
                n2 = n2 + 1
                if (k == l) then
                  bb = 1.0d00
                end if
                il = ka + l
                jl = kb + l
                kl = kc + l
                kr = kr + 1
                a = w(kref+(n2-1)*nn+n1)
                f(ij) = f(ij) + bb * a * ptot(kl)
                f(kl) = f(kl) + aa * a * ptot(ij)
                a = a * aa * bb * 0.25d0
                f(ik) = f(ik) - a * p(jl)
                f(il) = f(il) - a * p(jk)
                f(jk) = f(jk) - a * p(il)
                f(jl) = f(jl) - a * p(ik)
              end do
            end do
          end do
        end do
      end if
    end subroutine fockdorbs
