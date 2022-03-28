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

module cosmo_C
  implicit none
  logical :: iseps, noeps, useps, lpka
  integer :: nspa, nps, nps2, nden, lenabc, nppa = 1082, &
  amat_dim, isude_dim, nipc, ioldcv
  integer, dimension (2) :: n0
  integer, dimension(:), allocatable :: &
  & iatsp,   & !
  & nar_csm, & !
  & nsetf,i, & !  Watch out for the "i"
  & ipiden,  & !
  & idenat,  & !
  & nset       !
  integer, dimension(:,:), allocatable :: &
  & isude,   & !
  & nn
  double precision ::&
    fepsi,       & ! Dielectric factor =  (e-1)/(e+0.5), e = dielectric constant
                   !
    rds, disex2, &
                   !
    ediel,       & ! Dielectric energy, in eV
                   !
    solv_energy, & ! Solvation energy, in eV
                   !
    area,        & ! Surface area, in square Angstroms
                   !
    cosvol,      & ! Molecular volume, in cubic Angstroms
                   !
    cif1, cif2,  & !
    cdiagi,      &
   fnsq, rsolv
  double precision :: ffact = 0.d0
  double precision, dimension(:), allocatable :: diagsl
  double precision, dimension(3,3,1082) :: tm
  double precision, dimension(4,1082) :: dirsm, dirvec
  double precision, dimension(:), allocatable :: &
  & amat,    & !
  & cmat,    & !
  & gden,    & !
  & qscat,   & !
  & arat,    & !
  & srad,    & !
  & abcmat,  & !
  & qden,    & !
  & bh,      & !
  & cdiag
  double precision, dimension(:,:), allocatable :: &
  & bmat,    & !
  & phinet,  & !
  & qscnet,  & !
  & qdenet,  & !
  & cosurf,  & !
  & sude,    & !
  & xsp,     & !
  & cxy
end module cosmo_C
