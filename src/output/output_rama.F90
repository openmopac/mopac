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

subroutine output_rama
!
! Write out Ramachandran angles, if keyword "RAMA" is present.
!
! These angles are calculated in get_angles each time they are printed.
!
  use common_arrays_C, only : txtatm
  use molkst_C, only : keywrd
  use MOZYME_C, only : angles, res_start, uni_res
  use chanel_C, only: iw
  implicit none
  integer :: i, j
  if (index(keywrd, " RAMA") == 0) return
  call get_angles()
  if (uni_res == 0) return
  write(iw,"(/22x,a)")"Ramachandran Angles"
  write(iw,"(/15x, a, 8x, a/)")"Residue","Phi    Psi  Omega"
  do i = 1, uni_res
    if (Abs(angles(1,i)) + Abs(angles(3,i))> 1.d-20 .and. res_start(i) > 0 .and. &
      txtatm(res_start(i))(1:4) == "ATOM") then
      if (Abs(angles(1,i)) > 1.d-20 .and. Abs(angles(2,i)) > 1.d-20 ) then
        write(iw,"(14x,a, 3x, 3f7.1, a)")txtatm(res_start(i))(18:26), (angles(j,i),j = 1,3)
        else if  (Abs(angles(1,i)) > 1.d-20) then
        write(iw,"(14x,a, 3x,f7.1, 2a)")txtatm(res_start(i))(18:26), angles(1,i), &
        "    -  ","    -  "
        else
        write(iw,"(14x,a, 3x,a, 3f7.1)")txtatm(res_start(i))(18:26), &
        "    -  ", (angles(j,i),j = 2,3)
      end if
    end if
  end do
  write(iw,*)" "
  return
end subroutine output_rama


subroutine get_angles ()
!
! Calculate the Ramachandran angles and the number of unique residues.
!
! Nothing else should be altered.
!
! Array "ions" is needed by subroutine "names"
!
  use common_arrays_C, only : nat, txtatm
  use MOZYME_C, only : ions, angles, allres, ib, allr, iopt, &
    start_res, lstart_res, uni_res
  use molkst_C, only: natoms, numat, moperr, maxatoms
  use atomradii_C, only: atom_radius_covalent
  implicit none
!
  double precision, dimension(:), allocatable :: radius
  logical, save :: l_protein
  logical, dimension (:), allocatable :: ioptl
  integer :: i, ires, nfrag, io, j, n1, alloc_stat, mres, &
    nn1,  delta_res, max_frag
!
  if (allocated(ions))    deallocate (ions)
  if (allocated(iopt))    deallocate(iopt)
  if (Allocated (ib))     deallocate (ib)
  allocate (ions(maxatoms),  ib(maxatoms),  &
          & ioptl(maxatoms), iopt(maxatoms), radius(maxatoms), &
          stat=alloc_stat)
  if (alloc_stat /= 0) return
  l_protein = .false.
  ib(:) = 0
  mres = 0
  ions = 0
  call extvdw_for_MOZYME (radius, atom_radius_covalent)
!
!   WORK OUT WHAT ATOMS ARE BONDED TO EACH OTHER.
!
  call lewis (.true.)
  if (moperr) return
!
!  FIND THE NITROGEN ATOM OF THE N END OF THE PROTEIN.
!
  ioptl(:numat) = .false.
  call findn1 (n1, ioptl, io, delta_res)
  txtatm(:numat) = " "
  angles = 0.d0
  allres = " "
  l_protein = (n1 /= 0)
  ib(:numat) = -100000
  ires = 0
  uni_res = 0
!
!  Break all intra-chain bonds, so that the residues can easily be
!  identified.
!
  call lyse !
  allr = " "
  do max_frag = 1, 10000
    if (start_res(max_frag) == -200) exit
  end do
  max_frag = max_frag - 1
  lstart_res = .false.
  nfrag = 0
  do
    if (max_frag > 0) then
      do nfrag = 1, max_frag
        if (.not. lstart_res(start_res(nfrag) + 1)) exit
      end do
    else
      nfrag = nfrag + 1
    end if
    if (max_frag > 0) then
      if (nfrag > max_frag) exit
      ires = start_res(nfrag) - delta_res
    end if
!
    call names (ioptl, ib, n1, ires, nfrag, io, uni_res, mres)
    if (moperr) return
!
    nn1 = n1
    call findn1 (n1, ioptl, io, delta_res)
!
! n1 is the start of the next fragment
! delta_res is the distance back along the chain to the start of the next fragment.
!
    if (.not. l_protein) nfrag = 0
    if (n1 == 0) exit
    if (n1 == nn1) ioptl(n1) = .true.
  end do
!
!  Re-evaluate all residues
!
  j = 1
  allres(j) = txtatm(1)(18:20)
  do i = 2, natoms
    if (txtatm(i) == " ")exit
    if (nat(i) /= 1 .and. txtatm(i)(23:27) /= txtatm(i - 1)(23:27)) then
      j = j + 1
      allres(j) = txtatm(i)(18:20)
    end if
  end do
  ires = j
  iopt(:natoms) = ib(:natoms)
!
!   LABEL THE ATOMS IN ANY NON-PROTEIN MOLECULES IN THE SYSTEM
!
  nfrag = nfrag + 1
  if (start_res(max(1,nfrag)) == -200) then
    if (.not. l_protein) ires = 0
  else
    ires = start_res(nfrag)
  end if
  call ligand (ires, start_res, nfrag)
!
!  Add chain letters
!
  call update_txtatm(.true., .true.)
  return
end subroutine get_angles
