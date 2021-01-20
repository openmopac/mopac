module cosmo_C
  USE vast_kind_param, ONLY:  double
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
  real(double) ::&
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
