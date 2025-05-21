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
   fnsq, rsolv
  double precision :: ffact = 0.d0
  double precision, dimension(:), allocatable :: diagsl
  double precision, dimension(3,3,1082) :: tm
  double precision, dimension(4,1082) :: dirsm, dirvec
  double precision, dimension(:), allocatable :: &
  & amat,    & !
  & gden,    & !
  & qscat,   & !
  & arat,    & !
  & srad,    & !
  & qden
  double precision, dimension(:,:), allocatable :: &
  & bmat,    & !
  & phinet,  & !
  & qscnet,  & !
  & qdenet,  & !
  & cosurf,  & !
  & sude
end module cosmo_C
