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

subroutine set_up_rapid (txt)
!
!  set_up_rapid creates the environment where SCF and gradient calculations can be run
!  using a small number of atoms instead of the whole set.  When first called (option ON)
!  a full SCF is done, then the atoms of interest are selected and moved to the start.
!
!  Until OFF is specified, the SCF and gradient calculations use the results of the sub-set plus
!  "depleted" arrays - arrays representing the atoms not used in the SCF.
!
!  txt has three setting:
!      ON    - turn on the RAPID option
!      OFF   - turn off the RAPID option
!      RESET - reverse the OFF option, i.e., turn RAPID on but don't do work such as SCF calculations
!
  use common_arrays_C, only : xparam, grad, coord, dxyz
  use molkst_C, only : escf, step_num, iflepo, moperr, nelecs, numat, use_ref_geo, keywrd
  use MOZYME_C, only : mode, nelred, numred
  implicit none
  character :: txt*2
  logical :: store_use_ref_geo
  integer, save :: store_numred, store_nelred, store_mode
  if (txt == "OF") then
!
!  Turn off the RAPID method
!
      store_numred = numred
      store_mode = mode
      store_nelred = nelred
      numred = numat
      mode = 0
      nelred = nelecs
      return
  else if (txt == "RE") then
!
!  Turn off the RAPID method
!
      numred = store_numred
      mode = store_mode
      nelred = store_nelred
      return
  end if
  if (txt /= "ON") then
    continue ! panic
  end if
  store_use_ref_geo = use_ref_geo
  use_ref_geo = .false.
  numred = numat
  mode = 0
  nelred = nelecs
!
!  When RAPID is used, do a full SCF + derivatives.  Then
!  identify those atoms that move.  Subsequent SCFs and derivatives
!  only involve those atoms that move.
!
!  First, run a normal, single-point SCF calculation.  This fills all
!  the arrays that will be used in future SCF calculations.
!  In particular, the arrays P, H, F, and dxyz are filled.
!
  call picopt(-1)
  mode = 0
  grad = 0.d0
  call compfg (xparam, .true., escf, .true., grad, (index(keywrd," RAPID") /= 0))
!  From here on, any call to COMPFG will use the subset of atoms.
!
  call pinout (1, .true.)
!
!  Identify all the atoms that move.  The number of moving atoms is
!  NUMRED
!
  call picopt (1)
!
  iflepo = 1
!
!  Increment step_num: this indicates that startup in the SCF should
!  be re-done.
!
  step_num = step_num + 1
!
!  The RAPID technique uses "depleted arrays"  - normal arrays that
!  would be used in a full SCF, but with all terms that refer to the
!  moving atoms deleted.  This is done in the following steps:
! 1.  The array is multipled by -1.
! 2.  The terms arising from the moving atoms are added to the array.
! 3.  The resulting array is again multiplied by -1.
!
!  That is, each array = the normal array - terms due to the moving atoms.
!  Each "depleted array" has the same "part<x>" where <x> refers to the
!  name of the normal array.
!
!  Second, recalculate the one-electron matrix, but this time build the array parth
!
  mode = -1
  call hcore_for_MOZYME ()
  if (moperr) return
!  Build the depleted array part_dxyz if RAPID is used
  if (index(keywrd," RAPID") /= 0)  call dcart(coord, dxyz)
!  Finally, set mode to +1.  The first time COMPFG is run, the following arrays will be built:
!    partp
!    partf
!
  mode = 1
  use_ref_geo = store_use_ref_geo
  return
end subroutine set_up_rapid
