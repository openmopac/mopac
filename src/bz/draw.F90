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

subroutine draw (iz, jz, jdir, a, b, c, pz, i1, i2, i3, mdim)
  use common_common, only : im, jm, r1
  implicit none
  integer, intent(in) :: i1, i2, i3, iz, jdir, jz, mdim
  double precision, intent(in) :: pz
  logical, dimension(40000), intent(inout) :: b
  logical, dimension(20000, 2), intent(out) :: c
  double precision, dimension(mdim, mdim), intent(in) :: a
!
  logical :: ai, aj, ak
  integer :: i, icount, iddx, iddy, idir, idx, idy, ipen, j, k, m, number
  double precision :: aa, ab, ad, cx = 0.d0, cy = 0.d0, ac, fa, factor, px, py
  integer, dimension(6) :: idirx, idiry
  data idiry/0, 1, 0, -1, 0, 1/
  data idirx/1, 0, -1, 0, 1, 0/
  ! 
  ! ... Executable Statements ...
  ! 
  icount = -5
  idir = jdir
  ipen = 2
  i = iz
  j = jz
  do
    aa = a(i, j)
    ai = (aa<0.)
    idx = idirx(idir)
    idy = idiry(idir)
    ab = a(i+idx, j+idy)
    factor = aa / (aa-ab)
    py = (im-i-factor*idx) / (im-1)
    px = (jm-j-factor*idy) / (jm-1)
    call euler (px*2.0, py*2.0, pz*2.0, i1, i2, i3, ipen)
    if (ipen==3 .or. r1>0.d0) then
      r1 = r1 + Sqrt ((cx-px)**2 + (cy-py)**2)
    end if
    cy = py
    cx = px
    icount = icount + 1
    if (idir/2*2 /= idir) then
      if (j==1 .or. j==jm) then
        number = idir * j
        if (number==3 .or. number==jm) exit
      end if
    else
      if (idir /= 2) then
        number = j - 1
        m = im
      else
        number = j
        m = 1
      end if
      k = (number-1)*im + i
      if (.not. b(k)) return
      b(k) = .false.
      if (i == m) return
    end if
    ipen = 3
    iddx = i + idirx(idir+1)
    iddy = j + idiry(idir+1)
    ac = a(iddx+idx, iddy+idy)
    aj = (ac<0.)
    ad = a(iddx, iddy)
    ak = (ad<0.)
    fa = 1
    if (aj .and. ak) then
      if (ai) goto 1000
    elseif (.not. (aj.or.ak)) then
      if (.not. ai) goto 1000
    else
      fa = 0.d0
      if ((ai.or..not.ak) .and. (.not.ai.or.ak)) goto 1000
      fa = aa*ac - ab*ad
      if (fa >= 0.) goto 1000
    end if
    idir = idir + 1
    if (idir == 5) then
      idir = 1
    end if
    cycle
    1000 i = iddx
    j = iddy
    if (fa /= 0.) then
      i = i + idx
      j = j + idy
      idir = idir - 1
      if (idir == 0) then
        idir = 4
      end if
    end if
  end do
  if (number == 3) then
    c(i-1, 1) = .false.
  else
    c(i, 2) = .false.
  end if
end subroutine draw
