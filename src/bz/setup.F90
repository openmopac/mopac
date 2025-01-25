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

subroutine setup (MOPAC_sec_det, nvecs, sec_det)
  !
  !
  !   Convert the secular determinant MOPAC_sec_det (which is a lower-half
  !   triangle over the Large Unit Cell) into a set of small 
  !   interaction matrices, sec_det.  sec_det has the structure:
  !
  ! sec_det(1:nvecs,(i*nr2*nr3 + j*nr3 + k )*nvecs+1:(i*nr2*nr3 + j*nr3 + k +1)*nvecs) 
  ! =  unit cell (0,0,0) interacting with cell (i,j,k)
  !
  use common_common, only : mr1, mr2, mr3, bcc, keywrd, iw_new
  implicit none
  integer, intent(in) :: nvecs
  double precision, dimension(*), intent(in) :: MOPAC_sec_det
  double precision, dimension(nvecs, nvecs*mr1*mr2*mr3), intent(out) :: sec_det
!
  logical :: prt
  integer :: i, ii, j, jj, k, l, lvecs, m, mprt, n, ncell, nprt
  double precision, external :: reada
  ! 
  !
  !   Set up individual Fock matrices.
  !
  prt = (Index (keywrd, " SEC ") /= 0)
  if (prt) then
    write(iw_new, "(a)") " Secular Determinant, as supplied by MOPAC"
    i = Index (keywrd, " SEC=")
    mprt = 2000
    if (i /= 0) then
      mprt = Nint (reada (keywrd, i+3))
    end if
  end if
  nprt = Min (8, nvecs)
  l = 0
  do i = 1, nvecs
    do j = 1, i
      l = l + 1
      sec_det(i, j) = MOPAC_sec_det(l)
      sec_det(j, i) = MOPAC_sec_det(l)
    end do
  end do
  if (prt) then
    write(iw_new, "(A,3I3)") " Unit Cell:", 0, 0, 0
    do jj = 1, nvecs
      write(iw_new, "(10f8.4)") (sec_det(ii, jj), ii = 1, nprt)
    end do
  end if
  ncell = 1
  do i = 1, mr1
    do j = 1, mr2
      do k = 1, mr3
        if (i + j + k /= 3) then
          if (.not.bcc .or. Mod(i + j + k, 2) == 1) then
            lvecs = nvecs * ncell
            ncell = ncell + 1
            do m = 1, nvecs
              lvecs = lvecs + 1
              do n = 1, nvecs
                l = (lvecs*(lvecs-1))/2 + n
                sec_det(n, lvecs) = MOPAC_sec_det(l)
              end do
            end do
            if (prt .and. ncell<=mprt) then
              write(iw_new, "(A,3I3)") " Unit Cell:", i-1, j-1, k-1
              do jj = lvecs-nvecs+1, lvecs
                write(iw_new, "(10F8.4)") (sec_det(ii, jj), ii = 1, nprt)
              end do
            end if
          end if
        end if
      end do
    end do
  end do
end subroutine setup
!
