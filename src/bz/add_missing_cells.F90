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

  subroutine add_missing_cells()
    use common_common, only : ncell, nijk
!
!  For every cell (i,j,k) there must exist a cell (-i,j,k) and a cell (i,-j,k) etc.
!
!  If any cells are missing, add them now.
!
    implicit none
    integer :: iloop, jloop, new_cells(4,500), ia, ib, ic, ja, jb, jc, &
      to_add
    logical :: cell_ok(500)
!
!  Make sure that for every cell (i,j,k) there exists a partner cell (-i,-j,-k)
!
    cell_ok = .false.
    to_add = 0
    do iloop = 1, ncell
      if (cell_ok(iloop)) cycle
      ia = nijk(1,iloop)
      ib = nijk(2,iloop)
      ic = nijk(3,iloop)
      j_loop: do jloop = iloop + 1, ncell
        if (cell_ok(jloop)) cycle
        ja = nijk(1,jloop)
        jb = nijk(2,jloop)
        jc = nijk(3,jloop)
!
!  Check cell
!
        if ((ia == 0 .or. ia == -ja) .and. (ib == 0 .or. ib == -jb) .and. (ic == 0 .or. ic == -jc)) then
!
!   The cell (ja, jb, jc) already exists.  Mark it as existing and move on.
!
            if (ia /= 0 .or. ib /= 0 .or. ic /= 0) cell_ok(jloop) = .true.
            exit j_loop
        end if
      end do j_loop
!
!    Cells needed are in need_cell
!
      if (jloop > ncell) then
        to_add = to_add + 1
        new_cells(1,to_add) = -ia
        new_cells(2,to_add) = -ib
        new_cells(3,to_add) = -ic
        new_cells(4,to_add) = iloop
      end if
    end do
  !  write(iw_new,'(a)')"  Existing cells"
  !  write(iw_new,'(3i3)')nijk(:,:ncell)
  !  write(iw_new,'(a)')"  Missing cells"
  !  write(iw_new,'(4i3)')new_cells(:,:to_add)
    call expand(new_cells, to_add)
  end subroutine add_missing_cells
  subroutine expand(new_cells, to_add)
    use common_common, only : ncell, nijk, nvecs, sec_det
    implicit none
    integer :: to_add, new_cells(4, to_add)
!
    integer :: iloop, jloop, kloop
!
!  Expand the set of cells
!
    do kloop =  1, to_add
      nijk(:,kloop + ncell) = new_cells(:3,kloop)
    end do
!
!  Expand the secular determinant
!
    do kloop = 1, to_add
      do iloop = 1, nvecs
        do jloop = 1, nvecs
          sec_det(iloop, jloop + (kloop + ncell - 1)*nvecs) = sec_det(jloop, iloop + (new_cells(4, kloop) - 1)*nvecs)
        end do
      end do
    end do
    ncell = ncell + to_add
  end subroutine expand