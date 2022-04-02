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

      subroutine getocc(aoc, boc, ao1, bo1, io, iv, nio)
      use reimers_C, only : nov, ncore
      implicit none
      logical*1 ::     aoc(nov), boc(nov), &
                      ao1(nov), bo1(nov)
      integer ::       io(4), iv(4), nio
      integer ::      i, j, niv

!***********************************************************************
!   RMG - get the occupied/virtual orbitals involved in a configuration
!   aoc, boc = occupations of excited configuration
!   ao1, bo1 = occupations of ground configuration
!   io, iv = lists of occupied/virtural orbitals
!   nio     = number of excitations in the configuration
!***********************************************************************
      do i = 1, 4
        io(i) = 0
        iv(i) = 0
      end do

      nio = 0
      niv = 0

!      write (6, *) aoc, boc
!      write (6, *) ao1, bo1

!     **** loop over orbitals, check diffs from first ref det ****
      do j = 1, nov
        if (ao1(j) .neqv. aoc(j)) then
!         **** alpha differs from first ref det ****
          if (ao1(j)) then
            nio = nio + 1
            io(nio) = j + ncore
          else
            niv = niv + 1
            iv(niv) = j + ncore
          end if
        end if

        if (bo1(j) .neqv. boc(j)) then
!         **** beta differs from first ref det ****
          if (bo1(j)) then
            nio = nio + 1
            io(nio) = j + ncore
          else
            niv = niv + 1
            iv(niv) = j + ncore
          end if
        end if
      end do
      if (nio /= niv) write (6, *) 'Error: Different number of occ ', &
         'and vir orbitals found', nio, niv, io, iv
      return
      end subroutine getocc

      subroutine diagci(cimat, aocc, bocc, diagsl)
      use molkst_C, only : mpack
      use reimers_C, only : nconf, nov, nex, matind, mspn, &
                           aor1, bor1

      implicit none
      double precision :: cimat(*), deltap(mpack), &
                      diagsl(nconf)
      logical*1 :: aocc(nov, mspn, nex), bocc(nov, mspn, nex)
      integer :: iconf, io(4), iv(4), nm, nio
      double precision :: ener
!***********************************************************************
!   RMG - compute the correction to the diagonal terms of the CI matrix
!   due to solvation, as in Klamt, J. Phys. Chem. 1996, 100, 3349.
!
!   For each excitation, compute the density change, then compute the
!   correction to the diagonal energy term. Correct the CI matrix, but
!   save all corrections in dcie for later use in unsolv.
!***********************************************************************
      diagsl(1) = 0.d0
      do iconf = 2, nconf
        call getocc(aocc(1:nov, 1, iconf), bocc(1:nov, 1, iconf), &
                    aor1(1:nov), bor1(1:nov), io, iv, nio)
! Compute the density change
        call exdeltap(io, iv, nio, deltap)

! Compute the energy change
        call solenr(deltap, ener)
        diagsl(iconf) = ener
        nm = matind(iconf) + iconf
        cimat(nm) = cimat(nm) + ener
      end do

      return
      end subroutine diagci

      subroutine corrci(cimat, cist, aocc, bocc, diagsl, eval)
      use molkst_C, only : mpack, keywrd
      use reimers_C, only : nconf, nov, nex, mspn
      USE chanel_C, only : iw
      implicit none
      double precision :: cimat(nconf*(nconf + 1)/2), deltap(mpack), &
                      diagsl(nconf), eval(nconf), cist(nconf, nconf), &
!                      dpall(mpack, nconf), ucener(nconf), eorig(nconf)
                      ucener(nconf), eorig(nconf)
      logical*1 :: aocc(nov, mspn, nex), bocc(nov, mspn, nex)

      integer :: ist
      double precision :: ener
!***********************************************************************
!   RMG - compute the correction to the final CI states. First compute
!   the energies without the diagonal corrections. Then compute the
!   state densities and updated energy corrections for the proper
!   densities.
!***********************************************************************

      do ist = 1, nconf
        eorig(ist) = eval(ist)
      end do

! Compute uncorrected state energies if necessary
      if (index(keywrd, ' KL2') == 0) then
        call unsolv(cimat, cist, diagsl, ucener)
      else
        ucener = eorig
      end if

! For each state, compute the proper energetic correction
      do ist = 1, nconf
        call stdeltap(cist, aocc, bocc, deltap, ist)
        call solenr(deltap, ener)
        eval(ist) = ucener(ist) + ener
      end do

!      write (6, *) 'Energy corrections'
      write (iw, *) ''
      write (iw, *) 'COSMO corrections to state energies'
      write (iw, *) 'State    Diag E (eV)   Uncorr E (eV)   Final E (eV)'
      do ist = 1, nconf
        !write (6, '(i6, 3f15.5)') ist, eorig(ist), ucener(ist), eval(ist)
        write (iw, '(i6, 3f15.5)') ist, eorig(ist), ucener(ist), eval(ist)
      end do

      return
      end subroutine corrci


      subroutine unsolv(cimat, cist, diagsl, ucener)
      use reimers_C, only : nconf, matind
      implicit none
      double precision :: cimat(nconf*(nconf + 1)/2), cist(nconf, nconf), &
                      diagsl(nconf), ucener(nconf)
      integer :: iconf, ist, nm, jconf
      double precision :: ener, etmp
!***********************************************************************
!   RMG - Compute the energies of the CI states without the diagonal
!   energy corrections. Un - correct the CI matrix, then multiply with the
!   earlier eigenvectors to get the new 'eigenvalues'.
!***********************************************************************
! Uncorrect the CI matrix
      do iconf = 1, nconf
        nm = matind(iconf) + iconf
        cimat(nm) = cimat(nm) - diagsl(iconf)
      end do
!
! Compute expectation values for the energies of the old eigenvectors
! with the new CI matrix to get uncorrected energies
      do ist = 1, nconf
        ener = 0.d0
        do iconf = 1, nconf
          etmp = 0.d0
          do jconf = 1, nconf
            nm = matind(max(iconf, jconf)) + min(iconf, jconf)
            etmp = etmp + cimat(nm)*cist(jconf, ist)
          end do
          ener = ener + cist(iconf, ist)*etmp
        end do
        ucener(ist) = ener
      end do

      return
      end subroutine unsolv


      subroutine stdeltap(cist, aocc, bocc, dpstat, ist)
      use molkst_c, only : mpack
      use reimers_C, only : nconf, nov, nex, &
          mspn, aor1, bor1
      implicit none
      double precision :: dpstat(mpack), cist(nconf, nconf), &
                      deltap(mpack)
      logical*1 :: aocc(nov, mspn, nex), bocc(nov, mspn, nex)
      integer :: iorb, ist, iconf
      integer :: io(4), iv(4), nio
      double precision :: prob(nconf)
!***********************************************************************
!   RMG - Compute the AO density for the full CI states. Loop over all
!   involved excitations to get the overall change in density.
!***********************************************************************

! Set up an empty matrix to hold all changes in density
      do iorb = 1, mpack
        dpstat(iorb) = 0.d0
      end do

! For each excitation, compute the change in density
      do iconf = 1, nconf
! Determine occupied & unoccupied orbitals in this configuration
        call getocc(aocc(1:nov, 1, iconf), bocc(1:nov, 1, iconf), &
                    aor1(1:nov), bor1(1:nov), io, iv, nio)
! Compute the density change
        call exdeltap(io, iv, nio, deltap)

! For each CI state, add the change in density with the appropriate
! coefficient
        prob(ist) = cist(iconf, ist)**2
        if (prob(ist) > 0.0001d0) then
          do iorb = 1, mpack
            dpstat(iorb) = dpstat(iorb) + prob(ist)*deltap(iorb)
          end do
        end if
      end do

      return
      end subroutine stdeltap


      subroutine exdeltap(io, iv, nio, deltap)
      use molkst_c, only : mpack
      use reimers_C, only : n, matind, cc0
      implicit none
      double precision  :: deltap(mpack)
      integer :: i, j, ia, occ, vir, &
                 io(4), iv(4), nio
!***********************************************************************
!   RMG - compute the change in density for a particular excitation.
!***********************************************************************
      do i = 1, mpack
        deltap(i) = 0.d0
      end do

      do j = 1, nio
        occ = io(j)
        vir = iv(j)
        do i = 1, n
          ia = matind(i) + i
          deltap(ia) = deltap(ia) + cc0(vir, i)**2 - cc0(occ, i)**2
        end do
      end do

      return
      end subroutine exdeltap


      subroutine solenr(deltap, dcie)
      use cosmo_C, only : amat, bmat, nden, gden, ipiden, nps, &
           nsetf, fnsq
      use molkst_c, only : mpack
      use funcon_C, only : a0, ev
      implicit none
      double precision  :: deltap(mpack), qsc(nps), phi(nps), qden(nden)
      integer :: i, j
      double precision :: dcie, fcon
!***********************************************************************
!   RMG - adapted from dmecip subroutine in MOPAC
!
!   Given the density deltap (change in AO occupation upon excitation),
!   computes the energetic stabilization due to electronic screening of
!   the density difference (dcie).
!
!   See Klamt, J. Phys. Chem. 1996, 100, 3349. Correction computed here
!   corresponds to the final term in Equation 10.
!***********************************************************************

      dcie = 0.d0
      fcon = a0 * ev

      do i = 1, nden
        qden(i) = gden(i) * deltap(ipiden(i))
      end do
      do i = 1, nps
        phi(i) = 0.d0
        do j = 1, nden
          phi(i) = phi(i) + bmat(j, i) * qden(j)
        end do
      end do

      call coscl2 (amat, nsetf, qsc, phi, nps)
      do i = 1, nps
        dcie = dcie + qsc(i) * phi(i)
      end do

      dcie = -dcie * fnsq * fcon / 2

      return
      end subroutine solenr

      subroutine staticsolv(ao, bo, ener)
      use cosmo_C, only : bmat, nden, gden, ipiden, nps, qscnet
      use molkst_c, only : mpack
      use funcon_C, only : a0, ev
      use reimers_C, only : n, nov, cc0, ncore

      implicit none
      double precision  :: p(mpack), phi(nps), qden(nden)
      logical*1 :: ao(nov), bo(nov)
      integer :: a, i, j, k, k0, kocc
      double precision :: ener, fcon
!***********************************************************************
!   RMG - Compute the interaction energy of an excited-state electron
!   density (p) with the ground-state surface charges (qdenet(2)).
!***********************************************************************
      ener = 0.d0
      fcon = a0 * ev

      do i = 1, mpack
        p(i) = 0.d0
      end do

      a = 0
      do i = 1, n
        do j = 1, i
           a = a + 1
           p(a) = p(a) + sum(cc0(1:ncore, i)*cc0(1:ncore, j)) * 2.d0
           do k = 1, nov
              k0 = k + ncore
              kocc = 0
              if (ao(k)) kocc = kocc + 1
              if (bo(k)) kocc = kocc + 1
              p(a) = p(a) + cc0(k0, i) * cc0(k0, j) * kocc

           end do
        end do
      end do

! Set up charge density for solvation
      do i = 1, nden
        qden(i) = gden(i) * p(ipiden(i))
      end do

! Compute electrostatic potential at each surface point
      do i = 1, nps
        phi(i) = 0.d0
        do j = 1, nden
          phi(i) = phi(i) + bmat(j, i) * qden(j)
        end do
      end do

! Calculate interaction energy of charge density with surface points
      do i = 1, nps
        ener = ener + qscnet(i, 2) * phi(i)
      end do
      ener = ener * fcon
      return
      end subroutine staticsolv
