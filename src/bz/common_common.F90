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

module common_common
!
!  Scalars
!
  implicit none
  double precision ::   &
    tvec(3,3),          & !
    trans(3,3),         & !
    xop(11,149),        & !
    allt1(9,9,130),     & !
    r1,                 & !
    top_l(3),           & !
    top_r(3),           & !
    bottom_l(3),        & !
    bottom_r(3),        & !
    
    dummy_r
  real ::               & !
    xscale,             & !
    yscale,             & !
    xoffset,            & !
    yoffset,            & !
    size_x,             & !
    size_y
    
  double precision, allocatable :: &
    grid(:,:,:)
    
  integer ::            &
    id,                 & !   Dimensionality of the system.  id = 1 = polymer; id = 2 = layer; id = 3 = solid
    numat,              & !   Number of atoms in the unit cell
    nvecs,              & !   Number of vectors = number of atomic orbitals = number of M.O.s in the unit cell,
                        !   or three times the number of atoms (in phonon work)
    mr1,                & !
    mr2,                & !
    mr3,                & ! 
    iop,                & !   Number of symmetry operations
    ir,                 & !   Input channel number
    ir_new,             & !   
    iw,                 & !   Output channel number
    iw_new=19,          & !
    nijk(3,1000),       & !   Indices of all atoms in all unit cells
    nijk_cuc(3,500),    & !   Indices of all atoms in the Central Unit Cell
    ixyz,               & !   Index of a specific unit cell
    mop(4,300,130),     & !   
    ncell,              & !
    im,                 & !
    jm,                 & !
    isurf,              & !
    dummy_i             !   Not used: dummy_i is present to terminate this block
  integer*2 ::          &
    top_left_x,         &
    top_left_y
  character ::          &
    jobnam*2000,        & !  Name of the job = name of the brz file less ".brz"
    data_set_name*2000, & !  Name of input file containing the specific work BZ has to do
    keywrd*2000,        & !
    units*9,            & !  Name of units (eV or cm**(-1))
    title(149)*12,      & !   
    Sym_Oper*12,        & !  Name of Symmetry Operation
    line*2000           !  Temporary storage
  logical ::            &
    bcc,                & !
    lijk(600),          & !  for lijk(x) .true. if there exists an nijk(-i,-j,-k,x) for  nijk(i,j,k,x)
    l_read_dat,         & ! .true. if <file>.dat is NOT used, in which case create and save a new <file>.dat
    phonon
  integer, dimension(:), allocatable :: nfirst, nlast, per_atom
  double precision, allocatable, dimension(:, :) :: coord, ref_coord, sec_det
end module common_common
