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

      subroutine scfmat (v, shft)

!     ************************************************************************
!     **** generates the Fock matrices F for each shell plus determins the ***
!     **** variational energy using trace formula, V.			   ***
!     **** then it constracts the total fock operator in A, 		   ***
!     **** using SHFT as a shift of the virtual energy levels in AII.	   ***
!     **** D = each shells density matrix (d, ff in symmetric storage mode)  ***
!     ****   D is destroyed if open shell or damped closed shell	   ***
!     **** P = atomic charges from each shell; BETA, U = 1 e matr; GAMMA = 2e ***
!     **** C = previous set of MO coefficients if first element  >  - 1.     ***
!     ************************************************************************
!     **** References:							   ***
!     ****   ROHF: Edwards and Zerner Theor Chim Acta 72 (1987) 347	   ***
!     ****	   eqn 27 for v, eqn 33 for shell Fock matrices		   ***
!     ****	   eqn 24 is solved in previous MO basis for the total	    **
!     ****	   Fock operator by putting				   ***
!     ****	     alam = n(nu)/n(mu)*LAMBDA(mu, nu) - LAMBDA(nu, mu)	   ***
!     ****	     so that cross term is alam*(F(nu) - F(mu))*n(mu)/n(nu)  ***
!     ****	     alam is arbitrary non - zero, is set to 1 herein	   ***
!     ****	   eqn 24 is extended by a VIRTUAL - VIRTUAL term designed   ***
!     ****	     to produce better virtual eigenvalues		   ***
!     ****	   Errors: eqn 24, no 1/n(mu); eqn 33 F(hat) = not F = 	   ***
!     ****   INDO (Rotationally invariant for Transition Metals):	   ***
!     ****	   Bacon and Zerner Theor Chim Acta 53 (1979) 21	   ***
!     ****	   with additional R integrals included			   ***
!     ************************************************************************

      use molkst_C, only : mpack, norbs
      use reimers_C, only : n, na, nb2, matind, occfr, nel, nshell, norbl, &
          norbh, nbf, ibf, iat, ig1, ig2, ig3, ig4, g, nirreg, natt, &
          avec1, bvec1, avec2, bvec2, ppg, pg, dia, dd, ff, aa, cc0, &
          beta, gamma

      implicit none
      double precision ::   q, v, avec10, &
                         ga, gk, hg, off, shft, wt
      integer ::            i, j, k, irreg, &
                         ib, is, ishell, jb, jshell, &
                         k1, k2, k3, k4, mu, nm, nu

      double precision :: eltot
      double precision :: dd1(mpack, 2)

      avec10 = 0.d0
      hg = 0.d0
      v = 0.d0

!     **** diagonal term, sum Coulomb integrals from all other atoms ****
      if (allocated(pg)) deallocate(pg)
      if (allocated(ppg)) deallocate(ppg)
      allocate(pg(na, 0:2))
      allocate(ppg(norbs, 0:2))

      do ishell = 0, nshell
        do ib = 1, n
          ppg(ib, ishell) = 0.d0
        end do
        do i = 1, na
          pg(i, ishell) = 0.d0
        end do
      end do

      do ib = 1, n
        i = iat(ib)

        do jb = 1, n
          j = iat(jb)
          ga = gamma(ib, jb)
          if (ib < jb) ga = gamma(jb, ib)
          nm = matind(jb) + jb
          call veccou (nm, avec1, bvec1)
          do ishell = 0, nshell
            if (ib == jb) then
!	      **** weighted charges ****
              pg(i, ishell) = pg(i, ishell) + avec1(ishell)
            else if (i /= j) then
!	      **** weighted sums over all other atoms ****
              ppg(ib, ishell) = ppg(ib, ishell) + avec1(ishell) * ga
            end if
          end do
        end do

      end do

!     ************ Shell Fock matrices for different Hamiltonians *************

!	**** INDO Hamiltonian ****
      do i = 1, na
        do mu = ibf(i), ibf(i) + nbf(i) - 1
!	  **** terms in diagonal atom - atom block ****

!	  **** zero counters for the diagonal energy ****
          do is = 0, nshell
            dia(is) = 0.d0
          end do

!	  **** loop over all other densities on this atom ****
        do nu = ibf(i), ibf(i) + nbf(i) - 1

!	    **** load the (nu, nu) elements of density matrix ****
            nm = matind(nu) + nu
            call veccou (nm, avec1, bvec1)
            if (mu > nu) then
!	      **** load the (mu, nu) elements of density matrix ****
              nm = matind(mu) + nu
              call veccou (nm, avec2, bvec2)
            end if

            do is = 0, nshell
              if (mu == nu) then
                dia(is) = dia(is) + gamma(mu, mu) * (avec1(is) - bvec1(is))
                avec10 = avec1(0)
              else if (mu > nu) then
!	        **** gamma(mu, nu) = (mu, mu|nu, nu); gamma(nu, mu) = (mu, nu|mu, nu) ***
                dia(is) = dia(is) + gamma(mu, nu)*avec1(is)&
     &                          - gamma(nu, mu)*bvec1(is)
!		**** same atom but different basis functions ****
                off = gamma(nu, mu) * (2.d0*avec2(is) - bvec2(is))&
     &            - gamma(mu, nu) * bvec2(is)
                if (is == 0) then
                  hg = off
                  v = v + off * avec2(0)
                else
                  ff(nm, is) = hg - off
                  v = v - off * dd(nm, is)
                end if
              else
                dia(is) = dia(is) + gamma(nu, mu)*avec1(is)&
     &                          - gamma(mu, nu)*bvec1(is)
              end if
            end do
          end do

!	  **** store fully diagonal terms ****
          nm = matind(mu) + mu
          hg = beta(nm) + ppg(mu, 0) + dia(0)
          v = v + (beta(nm) + hg) * avec10 * 0.5d0
          do is = 1, nshell
            q = ppg(mu, is) + dia(is)
            ff(nm, is) = hg - q
            v = v - q * dd(nm, is) * 0.5d0
          end do

!	  **** off - diagonal atom - atom blocks and energy contribution ****
          do j = 1, i - 1
            do nu = ibf(j), ibf(j) + nbf(j) - 1
              nm = matind(mu) + nu
              call veccou (nm, avec1, bvec1)
              hg = beta(nm) - bvec1(0) * gamma(mu, nu)
              v = v + (beta(nm) + hg) * avec1(0)
              do is = 1, nshell
                q = -bvec1(is) * gamma(mu, nu)
                ff(nm, is) = hg - q
                v = v - q * dd(nm, is)
              end do
            end do
          end do
        end do

!	**** Zerners irregular integrals ****
        k = natt(i)
        if (k > 0) then
          mu = ibf(i)
          do irreg = 1, nirreg
            k1 = ig1(irreg, k)
            k2 = ig2(irreg, k)
            k3 = ig3(irreg, k)
            k4 = ig4(irreg, k)
            gk = g(irreg, k)
            call scfirr (v, mu, gk, k1, k2, k3, k4)
            call scfirr (v, mu, gk, k1, k2, k4, k3)
            call scfirr (v, mu, gk, k3, k4, k1, k2)
            call scfirr (v, mu, gk, k4, k3, k1, k2)
            if (k1 /= k2) then
              call scfirr (v, mu, gk, k2, k1, k3, k4)
              call scfirr (v, mu, gk, k2, k1, k4, k3)
              call scfirr (v, mu, gk, k3, k4, k2, k1)
              call scfirr (v, mu, gk, k4, k3, k2, k1)
            end if
          end do
        end if

      end do

!     ********************* construct total fock matrix **********************

      if (nshell == 1 .and. shft == 0.d0) then
!	**** Closed shell, no damping so just copy single Fock matrix  ****
        do nm = 1, nb2
          aa(nm) = ff(nm, 1)
        end do

      else
       dd1 = dd

!	***** first, zero all matrix elements of total fock matrix ****
        do nm = 1, nb2
          aa(nm) = 0.d0
        end do
!	**** transform Fock matrices into previous MO basis set ****
!	**** doing projections of E and Z eqn 24		****
!	**** uses first two density matrices as work space	****

!	**** these coupling coeffs are almost arbitrary, choose irrational ****
        eltot = 0.d0
        do ishell = 1, nshell
          eltot = eltot + nel(ishell)
        end do

        do ishell = 1, nshell
          wt = 1.d0
!	  **** diagonal terms ISHELL - ISHELL ****
          call ao2mo1 (aa, ff(:, ishell), cc0, dd1, norbl(ishell), norbh(ishell), &
     &          norbl(ishell), norbh(ishell), wt)
!	  **** terms ISHELL - VIRTUAL ****
          call ao2mo1 (aa, ff(:, ishell), cc0, dd1, norbl(ishell), norbh(ishell), &
     &          norbh(nshell) + 1, n, wt)
!	  **** JRRs xtra VIRTUAL - VIRTUAL term (improves order of virtuals) ****
          wt = nel(ishell) / eltot
          call ao2mo1 (aa, ff(:, ishell), cc0, dd1, norbh(nshell) + 1, n, &
     &          norbh(nshell) + 1, n, wt)

!	  **** off - diagonal occupied terms ISHELL - JSHELL ****
          do jshell = ishell + 1, nshell
            wt = 1.d0
            call ao2mo1 (aa, ff(1, ishell), cc0, dd1, norbl(ishell), norbh(ishell), &
     &            norbl(jshell), norbh(jshell), wt)
            wt = -occfr(jshell) / occfr(ishell)
            call ao2mo1 (aa, ff(1, jshell), cc0, dd1, norbl(ishell), norbh(ishell), &
     &            norbl(jshell), norbh(jshell), wt)

          end do
        end do

!	**** shift virtual energies to dampen SCF convergence problems ****

        do is = 1, max(1, nshell - 1)
          do i = norbh(is) + 1, n
            nm = matind(i) + i
            aa(nm) = aa(nm) + shft
          end do
        end do

!	**** back - transform Fock matrix into AO basis ****

        call mo2ao (aa, cc0, dd1, n)

        if (cc0(1, 1) == 0.d0) then
          do nm = 1, nb2
            aa(nm) = ff(nm, 1)
          end do
        end if
      end if
      return
      end subroutine scfmat

!     *************************************************************************

      subroutine scfirr (v, mu, g, k1, k2, k3, k4)

      use reimers_C, only : nshell, avec1, bvec1, avec2, bvec2, matind, ff, dd

      implicit none
      double precision ::   g, q, v, hg, v1
      integer ::            is, k1, k2, k3, k4, mu, nm

!     **** inserts Zerners irregular ints into atomic Fock matrix block ****
!     **** add and contribution to total energy v			 ****
!     **** mu = nber of first basis function of this atom		 ****

      if (k2 >= k1) then
!	**** Coulomb integral, use A vec ****

!	**** mixing of dens(3, 4) ****
        nm = matind(mu + max(k3, k4)) + mu + min(k3, k4)
        call veccou (nm, avec1, bvec1)

!	**** mixing of dens(1, 2) ****
        nm = matind(mu + k2) + mu + k1
        call veccou (nm, avec2, bvec2)

        hg = avec1(0) * g
        v1 = avec2(0) * hg
        do is = 1, nshell
          q = avec1(is) * g
          ff(nm, is) = ff(nm, is) + hg - q
          v1 = v1 - q * dd(nm, is)
        end do
        if (k2 == k1) v1 = 0.5d0 * v1
        v = v + v1
      end if

      if (k3 >= k1) then
!	**** Exchange integral, use - B vec ****

!	**** mixing of dens(2, 4) ****
        nm = matind(mu + max(k2, k4)) + mu + min(k2, k4)
        call veccou (nm, avec1, bvec1)

!	**** mixing of dens(3, 1) ****
        nm = matind(mu + k3) + mu + k1
        call veccou (nm, avec2, bvec2)

        hg = -bvec1(0) * g
        v1 = avec2(0) * hg
        nm = matind(mu + k3) + mu + k1
        do is = 1, nshell
          q = -bvec1(is) * g
          ff(nm, is) = ff(nm, is) + hg - q
          v1 = v1 - q * dd(nm, is)
        end do
        if (k3 == k1) v1 = 0.5d0 * v1
        v = v + v1
      end if

      return
      end subroutine scfirr

!     *************************************************************************

      subroutine veccou (nm, avec1, bvec1)

!     **** determins A and B coeffs for shells of Fock matrix using	****
!     **** E&Z eqn 30c; closed shell terms are set to zero.		****
!     **** nm = matrix element index, shell 0 is the total density	****
!     **** the B vec has the half factor of eqn 30b included		****
      use reimers_C, only : vca, vcb, nshell, dd

      implicit none
      double precision ::   avec1(0:4), bvec1(0:4)
      integer ::            ishell, jshell, nm

      avec1(0) = dd(nm, 1)
      avec1(1) = 0.d0
      bvec1(1) = 0.d0
      do ishell = 2, nshell
        avec1(0) = avec1(0) + dd(nm, ishell)
        avec1(ishell) = 0.d0
        bvec1(ishell) = 0.d0
        do jshell = 2, nshell
          avec1(ishell) = avec1(ishell) + (1.d0 - vca(ishell, jshell)) *&
     &                                dd(nm, jshell)
          bvec1(ishell) = bvec1(ishell) + (1.d0 - vcb(ishell, jshell)) *&
     &                                dd(nm, jshell)
        end do
        bvec1(ishell) = bvec1(ishell) / 2.d0
      end do
      bvec1(0) = avec1(0) / 2.d0

      return
      end subroutine veccou

!     *************************************************************************

! ____________________________________________________________________________
! ________________________ SCF OUTPUT ROUTINES _______________________________
! ____________________________________________________________________________

      subroutine output (c, aii)

      USE chanel_C, only : iw
      use reimers_C, only : n, nshell, norbl, norbh, &
          nmrep, nptg, nr, nsym

      implicit none
      double precision ::  aii(n), c(n, n)
      character*1       ll(5)
      integer ::           notype(8, 5), notot(8)
      integer ::           i, ir, is, nmoou1, nmoou2

!     **** write number of each orb of each symmetry in each shell ****
!     **** get shell operator name ****

      nmoou1 = 1
      nmoou2 = n
      nr = 1
      nptg = 1

      if (nshell == 1) then
        ll(1) = '1'
        ll(2) = 'v'
      else
        ll(1) = '1'
        ll(2) = '2'
        ll(3) = 'v'
      end if

      if (allocated(nsym)) deallocate(nsym)
      allocate(nsym(n))
      do i = 1, n
        nsym(i) = 1
      end do

      do is = 1, nshell + 1
        do ir = 1, nr
          notype(ir, is) = 0
        end do
      end do

      do is = 1, nshell
        do i = norbl(is), norbh(is)
          notype(1, is) = notype(1, is) + 1
        end do
      end do
      is = nshell + 1
      do i = norbh(nshell) + 1, n
          notype(1, is) = notype(1, is) + 1
      end do
      do ir = 1, nr
        notot(ir) = 0
        do is = 1, nshell + 1
          notot(ir) = notot(ir) + notype(ir, is)
        end do
      end do

      write (iw, "(/' Summary of symmetry occupancy in each shell:'//' SHELL  ', 9(a3, 1x))") &
        (nmrep(ir, nptg), ir = 1, nr)
      do is = 1, nshell + 1
        write (iw, '(4x, a1, 1x, 8i4)') ll(is), (notype(ir, is), ir = 1, nr)
      end do
      write (iw, '(1x, a5, 8i4)') 'TOTAL', (notot(ir), ir = 1, nr)

!     ***** page eject for MO dump ****


      if (nmoou1 > nmoou2) return

!     **** first eigenvalues before main dump ****

!     **** dump of molecular orbitals ****
      write (iw, "(/,/,20x,' MOLECULAR ORBITALS ',/,/)")
      call matout(c, aii, n, n, n)
      return
      end subroutine output

!     *************************************************************************

      subroutine gsdip (xz, d, zcore, dm)

!     **** calculates ground-state dipole moment and Mulliken charges ****
      use reimers_C, only : na, nb2, matind, debye, nbf, ibf, nbt, &
          qgs, dipgs, ndtype
      implicit none
      double precision ::  d(nb2), dm(nb2, 3), xz(na, 3), zcore(na), cm(3), &
                        pt(4), pm(3)
      integer ::           ind(0:8), i, j, k, ib, mu, nm, nu

!     ****		identify AO as s, p, d(eg) or d(t2g) ****
      data ind          /1, 2, 2, 2, 3, 3, 4, 4, 4/


      ndtype = 1

      if (allocated(qgs))  deallocate(qgs)
      if (allocated(dipgs))  deallocate(dipgs)
      allocate(qgs(na))
      allocate(dipgs(na, 3))

!     **** zero counters for atomic and hybridization contributions ****
      do i = 1, 3
        cm(i) = 0.d0
        pm(i) = 0.d0
      end do

      do i = 1, na
        qgs(i) = zcore(i)
        dipgs(i, 1) = 0.d0
        dipgs(i, 2) = 0.d0
        dipgs(i, 3) = 0.d0

!	**** break charge up into s, p, d(eg) and d(t2g) contributions ****
        do j = 1, 4
          pt(j) = 0.d0
        end do
        do ib = ibf(i), ibf(i) + nbf(i) - 1
          nm = matind(ib) + ib
          j = ind(nbt(ib))
          pt(j) = pt(j) + d(nm)
        end do

!	**** atomic charge and its contribution to dipole moment ****
        qgs(i) = qgs(i) - pt(1) - pt(2) - pt(3) - pt(4)
        do j = 1, 3
          cm(j) = cm(j) + qgs(i) * xz(i, j)
        end do

        if (ndtype /= 1 .or. nbf(i) <= 1) then
!	  **** no hybridization contribution from this atom ****
        else
!	  **** determine hybridization contribution ****
          do mu = ibf(i), ibf(i) + nbf(i) - 1
            do nu = ibf(i), mu - 1
              nm = matind(mu) + nu
              do k = 1, 3
                dipgs(i, k) = dipgs(i, k) + d(nm) * dm(nm, k)
              end do
            end do
          end do
          do j = 1, 3
            dipgs(i, j) = dipgs(i, j) * debye * 2.d0
            pm(j) = pm(j) + dipgs(i, j)
          end do
        end if

      end do

      do j = 1, 3
        cm(j) = cm(j) * debye
      end do
      return
      end subroutine gsdip
