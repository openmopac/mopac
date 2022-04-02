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

      subroutine ciosci(vects, lroot, oscil, conf)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use molkst_C, only: numat,  norbs, keywrd
      use common_arrays_C, only : coord, nfirst, nlast, nat
      USE parameters_C, only : dd, zp, zd
      use funcon_C, only : a0
      use meci_C, only : nstate, lab, nmos, microa, microb
      use chanel_C, only: iw
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      double precision , intent(in) :: vects(norbs,norbs)
      double precision , intent(inout) :: oscil(3,*)
      double precision, dimension(:,:), allocatable :: t4
      double precision , intent(in) :: conf(lab*lab)
      integer :: lroot
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(14) :: ll
      integer :: loop, iloop, iatom, i, j, k, l, m, ii, ij, jj, istate, &
      nd, npp, ns
      double precision, dimension(norbs,nmos) :: vect1, vect2
      double precision, dimension(nmos,nmos) :: t2
      double precision, dimension(lab) :: work
      integer, dimension (107), save :: nspqn
      double precision, dimension (numat) :: sppol
      double precision, dimension (:), allocatable :: pdpol
      character (len=1) :: type (3)
      double precision :: sum
      logical :: debug
!***********************************************************************
!
!    CIOSCI evaluates the expectation value <STATE|x or y or z|STATE>
!
!   Info:   VECTS  = Molecular Orbitals
!           CONF   = State Eigenvectors
!           NSTATE = State to be excited from
!
!
!***********************************************************************
      data nspqn / 2 * 1, 8 * 2, 8 * 3, 18 * 4, 18 * 5, 32 * 6, 21 * 0 /
      data type / "x", "y", "z" /
      debug = (Index (keywrd, "CIOSCI") /= 0)
      allocate(t4(lab, lab))
   !
   !   a0 = 0.529
   !
   ! dd(i) = (2*n+1)*2**(2n+1)*(exp(s)*exp(p))**(n+1/2)
   !           ----------------------------------------  See calpar
   !          (sqrt(3)*(exp(s)+exp(p))**(2*n+2))
   !
   ! dd is the transition length <ns|r|np> in bohr
   !
      sppol(:numat) = a0 * dd(nat(:numat))
      ll(1) = 1
      do i = 2, 14
        ll(i) = ll(i-1) * i
      end do
  !
  !  For p-d transitions,
  !  pdpop(i) = (np+nd+1)!*2**(np+nd+1)*exp(p)**(np+1/2)*exp(d)**(nd+1/2)
  !              --------------------------------------------------------
  !             sqrt(5)*(exp(p)+exp(d)**(np+nd+2)*sqrt( (2*np)! * (2*nd)! )
  !
      allocate (pdpol(numat))
      do i = 1, numat
        if (nlast(i)-nfirst(i) > 5) then
          j = nat(i)
          nd = nspqn(j)
          npp = nd
          if (j > 20) npp = npp - 1
          pdpol(i) = a0 * ll(npp+nd+1) * 2 ** (npp+nd+1) * &
               & zp(j) ** (npp+ 0.5d0) * zd(j) ** (nd+0.5d0) / &
               & (Sqrt(5.d0)*(zp(j)+zd(j))**(npp+nd+2) * &
               & Sqrt(Dble(ll((2*npp))*ll(2*nd))))
        end if
      end do
      oscil(:3, :lab) = 0.d0
      do loop = 1, 3
        !
        !    SET UP THE M.O. UNITARY MATRIX, <PSI | IOPER | PSI>
        !
        do iloop = 1, nmos
          do iatom = 1, numat
            do i = nfirst(iatom), nlast(iatom)
              vect1(i, iloop) = vects(i, iloop)
              vect2(i, iloop) = vects(i, iloop) * coord(loop, iatom)
              if (nat(iatom) /= 1) then
  !
  !   Include the s-p transition on the same atom
  !
                if (i == nfirst(iatom)) then
                  vect2(i, iloop) = vect2(i, iloop) + sppol(iatom) * &
                       & vects(i+loop, iloop)
                else if (i-nfirst(iatom) == loop) then
                  vect2(i, iloop) = vect2(i, iloop) + sppol(iatom) * &
                       & vects(i-loop, iloop)
                end if
                if (i == nfirst(iatom)+8) then
  !
  !  k = start of "p" shell - 1;  l = start of "d" shell - 1.
  !
                  k = nfirst(iatom)
                  l = k + 3
  !
  ! Add in p - d transition.  Do all three (x, y, z) at one time
  !
                  if (loop == 1) then
                    vect2(k+1, iloop) = vect2(k+1, iloop) + pdpol(iatom) * &
                         & (vects(l+1, iloop)-1/Sqrt(3.d0)*vects(l+3, iloop))
                 !
                    vect2(k+2, iloop) = vect2(k+2, iloop) + pdpol(iatom) * &
                         & vects(l+5, iloop)
                 !
                    vect2(k+3, iloop) = vect2(k+3, iloop) + pdpol(iatom) * &
                         & vects(l+2, iloop)
  !
  !  "d" set
  !
                    vect2(l+1, iloop) = vect2(l+1, iloop) + pdpol(iatom) * &
                         & vects(k+1, iloop)
                    vect2(l+2, iloop) = vect2(l+2, iloop) + pdpol(iatom) * &
                         & vects(k+3, iloop)
                    vect2(l+3, iloop) = vect2(l+3, iloop) - 1 / Sqrt (3.d0) * &
                         & pdpol(iatom) * vects(k+1, iloop)
                    vect2(l+5, iloop) = vect2(l+5, iloop) + pdpol(iatom) * &
                         & vects(k+2, iloop)
                  else if (loop == 2) then
                    vect2(k+1, iloop) = vect2(k+1, iloop) + pdpol(iatom) * &
                         & vects(l+5, iloop)
                    vect2(k+2, iloop) = vect2(k+2, iloop) + pdpol(iatom) * &
                         & (-vects(l+1, iloop)-1/Sqrt(3.d0)*vects(l+3, iloop))
                    vect2(k+3, iloop) = vect2(k+3, iloop) + pdpol(iatom) * &
                         & vects(l+4, iloop)
                    vect2(l+1, iloop) = vect2(l+1, iloop) - pdpol(iatom) * &
                         & vects(k+2, iloop)
                    vect2(l+3, iloop) = vect2(l+3, iloop) - 1 / Sqrt (3.d0) * &
                   & pdpol(iatom) * vects(k+2, iloop)
                    vect2(l+4, iloop) = vect2(l+4, iloop) + pdpol(iatom) * &
                         & vects(k+3, iloop)
                    vect2(l+5, iloop) = vect2(l+5, iloop) + pdpol(iatom) * &
                         & vects(k+1, iloop)
                  else
                    vect2(k+1, iloop) = vect2(k+1, iloop) + pdpol(iatom) * &
                         & vects(l+2, iloop)
                    vect2(k+2, iloop) = vect2(k+2, iloop) + pdpol(iatom) * &
                         & vects(l+4, iloop)
                    vect2(k+3, iloop) = vect2(k+3, iloop) + pdpol(iatom) * &
                         & 2.d0 / Sqrt (3.d0) * vects(l+3, iloop)
                 !
                    vect2(l+2, iloop) = vect2(l+2, iloop) + pdpol(iatom) * &
                         & vects(k+1, iloop)
                    vect2(l+3, iloop) = vect2(l+3, iloop) + 2 / Sqrt (3.d0) * &
                         & pdpol(iatom) * vects(k+3, iloop)
                    vect2(l+4, iloop) = vect2(l+4, iloop) + pdpol(iatom) * &
                         & vects(k+2, iloop)
                  end if
                end if
              end if
            end do
          end do
        end do
        !
        !    VECT1 HOLDS THE MOLECULAR ORBITALS IN THE COORDINATE SYSTEM OF
        !          SYMTRZ.
        !    VECT2 HOLDS THE SAME ORBITALS, AFTER BEING OPERATED ON BY THE
        !          DIPOLE OPERATORS 'X', 'Y', AND 'Z'.
        !
        !     T2 HOLD THE UNITARY TRANSFORM FOR THE SET OF M.O.S
        !        <VECT1|IOPER|VECT2>
        !
        do i = 1, nmos
          t2(i, i) = 0.d0
          do j = 1, i - 1
            if (i /= j) then
              sum = 0.d0
              do k = 1, norbs
                sum = sum + vect1(k, i) * vect2(k, j)
              end do
              t2(i, j) = sum
              t2(j, i) = -sum
            end if
          end do
        end do
        if (debug) then
          write (iw, "(/,A)") "              EFFECT OF DIPOLE OPERATOR ON M.O.S"
          do j = 1, nmos, 2
            jj = Min (j+1, nmos)
            if (jj /= j) then
              write (iw, "(4(a,i1,a,i1),a)") "    A.O.           |M.O.", j, &
                   & ">         " // type (loop) // "|M.O.", j, &
                   & ">          |M.O.", j + 1, ">         " // type (loop) // &
                   & "|M.O.", j + 1, ">"
            else
              write (iw, "(2(a,i1,a,i1),a)") "    A.O.           |M.O.", j, &
                   & ">        " // type (loop) // "|M.O.", j, ">"
            end if
            do i = 1, norbs
              write (iw, "(i6,6x,4f17.12)") i, (vect1(i, k), vect2(i, k), &
                   & k=j, jj)
            end do
          end do
          write (iw, "(/6x,a)") "         Integrals <M.O.|" // type (loop) // &
               & "|M.O.>"
          write (iw, "(i12,8i17)") (i, i=1, nmos)
          do i = 1, nmos
            write (iw, "(I6,10(10F17.12,/6x))") i, (t2(j, i), j=1, nmos)
          end do
        end if
        !
        !   THE BIG LOOP TO FILL T4
        !
        do i = 1, lab
          do j = 1, i
            l = 0
            m = 0
            do k = 1, nmos
              if (microa(k, i) /= microa(k, j)) then
                l = l + 1
                ll(l) = k
              end if
              if (microb(k, i) /= microb(k, j)) then
                m = m + 1
                ll(m) = k
              end if
            end do
            if (l == 2 .and. m == 0 .or. m == 2 .and. l == 0) then
              t4(i, j) = t2(ll(1), ll(2))
              if (l == 2) then
                !
                !   Two microstates differing in two ALPHA electrons
                !
                ii = ll(1)
                !
                !   Add in MICROA(II,I) to allow for the fact that M.O. LL(1)
                !   is ALWAYS less than M.O. LL(2).
                !
                !   Step over MICROB(II,I)
                !
                ij = microb(ii, i) + microa(ii, i)
                do jj = ii + 1, ll(2) - 1
                  ij = ij + microa(jj, i) + microb(jj, i)
                end do
              else
                !
                !   Two microstates differing in two BETA electrons
                !
                ii = ll(1)
                !
                !   Add in MICROB(II,I) to allow for the fact that M.O. LL(1)
                !   is ALWAYS less than M.O. LL(2).
                !
                ij = microb(ii, i)
                do jj = ii + 1, ll(2) - 1
                  ij = ij + microa(jj, i) + microb(jj, i)
                end do
                    !
                    !  Must jump over MICROA(LL(2),I) to get to the beta M.O.
                    !
                ij = ij + microa(ll(2), i)
              end if
              if (Mod(ij, 2) == 1) then
                t4(i, j) = -t4(i, j)
              end if
                 !
                 !  Matrix is antisymmetric, therefore:
                 !
              t4(j, i) = -t4(i, j)
            else
              t4(i, j) = 0.d0
              t4(j, i) = 0.d0
            end if
          end do
        end do
        if (debug) then
          write (iw, "(//10x, ' EFFECT OF DIPOLE OPERATOR ON MICROSTATES')")
          if (lab < 6) then
            write (iw, "(i12,5i17)") (i, i=1, lab)
            do i = 1, lab
              write (iw, "(I6,3(5F17.12,2x))") i, (t4(j, i), j=1, lab)
            end do
          else
            i = -lab
            call matout (t4, t4, lab,i, lab)
          end if
        end if
        !
        !    NOW TO PERFORM <STATE(LROOT)| X OR Y OR Z  | STATE(ISTATE)>
        !
        do ns = 1, nstate
          ii = (lroot+ns-2) * lab
          do istate = 1, lab
            do j = 1, lab
              sum = 0.d0
              do k = 1, lab
                sum = sum + conf(k+ii) * t4(j, k)
              end do
              work(j) = sum
            end do
            sum = 0.d0
            do k = 1, lab
              sum = sum + work(k) * conf(k+ (istate-1)*lab)
            end do
            oscil(loop, istate) = oscil(loop, istate) + sum ** 2
          end do
        end do
      end do
      if (Index (keywrd, " LARGE") /= 0 .or. debug) then
        write (iw, "(//,A)") "           EFFECT OF DIPOLE OPERATOR ON STATES"
        write (iw, "(' STATE         X                Y                Z')")
        do i = 1, lab
          write (iw, "(I4,3F17.8)") i, (oscil(j, i), j=1, 3)
        end do
      end if
      deallocate (pdpol)
end subroutine ciosci
