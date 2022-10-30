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

  subroutine gdisp(r0ab, rs6, alp6, c6ab, s6, s8, mxc, r2r4, rcov, rs8, alp8, dxyz_temp)
!
! Calculates the derivative (gradient) of the D3 dispersion term,
! based on material provided by Stefan Grimme, University of Muenster, Germany
! v3.1 of the DFTD3 library
!
  use common_arrays_C, only: nat, cell_ijk, Vab
  use molkst_C, only : numat, l123, l1u, l2u, l3u
  implicit none
  double precision, parameter :: k1 = 16.d0
  integer, parameter :: max_elem = 94, maxc = 5 ! maximum coordination number references per element
  double precision :: &
  c6ab(max_elem, max_elem, maxc, maxc, 3), & ! C6 for all element pairs
  dxyz_temp(3,numat),                      & ! Contribution to gradient
  rcov(max_elem),                          & ! covalent radii
  cn(numat),                               & ! coordination numbers of the atoms
  r0ab(max_elem, max_elem),                & ! cut  -  off radii for all element pairs
  rs6, alp6, s6, s8, r2r4(max_elem), rs8, alp8
  integer :: mxc(max_elem)
!
!  Local variables
!
  integer :: i, j, linij, iii, jjj, i_cell, j_cell, kkkk
  double precision :: R0, r2, damp6, damp8, c6, tmp1, tmp2, r, dc6_rest, rij(3), dc6iji, dc6ijj, r6, r7, t6, t8, &
  rcovij, expterm, dcn,x1, r42, r8, r9
  double precision, allocatable :: drij(:), dc6i(:)
  double precision, external :: distance
  integer, external :: lin
!
! In this subroutine, array "coord" is in atomic units, not Angstroms.
!
  i = (numat*(numat + 1))/2
  allocate(drij(i), dc6i(numat))
  call ncoord(numat, rcov, nat, cn)
  drij = 0.d0
  dc6i = 0.d0
  dc6_rest = 0.0d0
  do i = 2, numat
    do j = 1, i - 1
      r = distance(j,i)
      if (r < 0.1d0) cycle
      r2 = r**2
      if (r2 > 10000.d0) cycle
      linij = lin(i,j)
      r0 = r0ab(nat(j),nat(i))
      r42=r2r4(nat(i))*r2r4(nat(j))
      call get_dC6_dCNij(maxc, max_elem, c6ab, mxc(nat(i)), &
            mxc(nat(j)), cn(i), cn(j), nat(i), nat(j), &
            c6, dc6iji, dc6ijj)
      r6 = r2*r2*r2
      r7 = r6*r
      r8=r6*r2
      r9=r8*r
!
!  Calculate damping functions
!
      t6  =  (r/(rs6*R0))**( - alp6)
      damp6  = 1.d0/( 1.d0 + 6.d0*t6 )
      t8 = (r/(rs8*R0))**(-alp8)
          damp8 =1.d0/( 1.d0+6.d0*t8 )
      tmp1 = s6*6.d0*damp6*C6/r7
      tmp2=s8*6.d0*C6*r42*damp8/r9
      drij(linij) = drij(linij) - tmp1 - 4.d0*tmp2
      drij(linij) = drij(linij)  + tmp1*alp6*t6*damp6 +3.d0*tmp2*alp8*t8*damp8
      dc6_rest = s6/r6*damp6 +3.d0*s8*r42/r8*damp8
!
!     saving all f_dmp/r6*dC6(ij)/dCN(i) for each atom for later
!
      dc6i(i) = dc6i(i) + dc6_rest*dc6iji
      dc6i(j) = dc6i(j) + dc6_rest*dc6ijj
    end do
  end do
!
! After calculating all derivatives dE/dr_ij w.r.t. distances,
! the grad w.r.t. the coordinates is calculated dE/dr_ij * dr_ij/dxyz_i
!
  do i = 2, numat
    do j = 1, i - 1
      linij = lin(i,j)
      r = distance(j,i)
      if (r < 0.1d0) cycle
      r2 = r**2
      rij = Vab
      r2 = sum(rij*rij)
      r = dsqrt(r2)
      if (r2 < 100.d0) then
        rcovij = rcov(nat(i)) + rcov(nat(j))
        expterm = exp( - k1*(rcovij/r - 1.d0))
        dcn =  - k1*rcovij*expterm/(r*r*(expterm + 1.d0)*(expterm + 1.d0))
      else
        dcn = 0.d0
      end if
      x1 = drij(linij) + dcn*(dc6i(i) + dc6i(j))
      iii = l123*(i - 1)
      jjj = l123*(j - 1)
      kkkk = (l3u - cell_ijk(3)) + (2*l3u + 1)*(l2u - cell_ijk(2) + (2*l2u + 1)*(l1u - cell_ijk(1))) + 1
      i_cell = iii + kkkk 
      j_cell = jjj + kkkk       
      dxyz_temp(:,i_cell) = dxyz_temp(:,i_cell) + x1*rij/r
      dxyz_temp(:,j_cell) = dxyz_temp(:,j_cell) - x1*rij/r
    end do
  end do
  return
  end subroutine gdisp


  subroutine ncoord(natoms, rcov, nat, cn)
  implicit none
  integer nat(*),natoms, i, j
  double precision cn(*), rcov(94)
  double precision r, damp, xn, rr, rco
  double precision, external :: distance
  do i = 1, natoms
    xn = 0.0d0
    do j = 1, natoms
        if (j /= i) then
          r = distance(i, j)
! covalent distance in Bohr
          rco = rcov(nat(i)) + rcov(nat(j))
          rr = rco/r
! counting function exponential has a better long - range behavior than MHGs inverse damping
          damp = 1.d0/(1.d0 + exp( -16.d0*(rr - 1.0d0)))
          xn = xn + damp
        end if
    end do
    cn(i) = xn
  end do
  return
  end subroutine ncoord

  subroutine get_dC6_dCNij(maxc ,max_elem, c6ab, mxci, mxcj,  &
    cni, cnj, izi, izj, c6check, dc6i, dc6j)
    implicit none
    integer, intent (in) :: maxc ,max_elem, mxci, mxcj
    double precision, intent (in) :: c6ab(max_elem,max_elem,maxc,maxc,3), cni, cnj
    integer :: izi, izj, a, b
    double precision :: dc6i, dc6j, c6check, term
    double precision :: zaehler, nenner, dzaehler_i, dnenner_i, dzaehler_j, dnenner_j
    double precision :: expterm, cn_refi, cn_refj, c6ref, r
    double precision :: c6mem, r_save
    c6mem =  -1.d99
    r_save = 9999.d0
    zaehler = 0.0d0
    nenner = 0.0d0
    dzaehler_i = 0.d0
    dnenner_i = 0.d0
    dzaehler_j = 0.d0
    dnenner_j = 0.d0
    do a = 1,mxci
      do b = 1,mxcj
        c6ref = c6ab(izi, izj, a, b, 1)
        if (c6ref > 0.d0) then
          cn_refi = c6ab(izi, izj, a, b, 2)
          cn_refj = c6ab(izi, izj, a, b, 3)
          r = (cn_refi - cni)*(cn_refi - cni) + (cn_refj - cnj)*(cn_refj - cnj)
          if (r < r_save) then
            r_save = r
            c6mem = c6ref
          end if
          expterm = exp(-4.d0*r)
          zaehler = zaehler + c6ref*expterm
          nenner = nenner + expterm
          expterm = -4.d0*expterm*2.d0
          term = expterm*(cni - cn_refi)
          dzaehler_i = dzaehler_i + c6ref*term
          dnenner_i  = dnenner_i + term
          term = expterm*(cnj - cn_refj)
          dzaehler_j = dzaehler_j + c6ref*term
          dnenner_j  = dnenner_j +      term
        end if
      end do
    end do

    if (nenner > 1.0d-99) then
      c6check = zaehler/nenner
      dc6i = ((dzaehler_i*nenner) - (dnenner_i*zaehler))/(nenner*nenner)
      dc6j = ((dzaehler_j*nenner) - (dnenner_j*zaehler))/(nenner*nenner)
    else
      c6check = c6mem
      dc6i = 0.0d0
      dc6j = 0.0d0
    end if
end subroutine get_dC6_dCNij
