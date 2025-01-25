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

subroutine orient (nvecs, sec_det, natot, coord)
!
!   Rotate system (atoms and interaction matrices) so that the first
!   translation vector is along the "x" axis, and the second vector,
!   if present, is in the "x-y" plane.
!
  use common_common, only : numat, id, tvec, xop, iop, phonon, jobnam, &
    nfirst, nlast, per_atom, ncell, nijk, title, mr1, mr2, mr3, bcc, trans, &
      keywrd, iw_new
  implicit none
  integer, intent(in) :: nvecs
  integer, intent(inout) :: natot
  double precision, dimension(3, natot+3), intent(inout) :: coord
  double precision, dimension(nvecs, *), intent(inout) :: sec_det
!
  logical :: prtok, prtops
  integer :: i, i1, i2, i3, i4, i99, ii, il, iofset, iu, j, j1, j2, j3, &
       & jj, jl, ju, k, k1, k2, l1, l2, mprt, nprt
  double precision :: sum, sum1, summax, summin, x, y, z, &
  a, b, c, d, e, f, g, h, o, two, six, ca, sa
  character :: line*120
  double precision, dimension(3) :: storen
  double precision, dimension(3, 3) :: xtoc, p_rot
  double precision, dimension(9, 9) :: t1, t2, e1, e2
  double precision, external :: reada
  equivalence (p_rot(1, 1), a)
  equivalence (p_rot(2, 1), b)
  equivalence (p_rot(3, 1), c)
  equivalence (p_rot(1, 2), d)
  equivalence (p_rot(2, 2), e)
  equivalence (p_rot(3, 2), f)
  equivalence (p_rot(1, 3), g)
  equivalence (p_rot(2, 3), h)
  equivalence (p_rot(3, 3), o)
! 
! ... Executable Statements ...
! 
  trans = 0.d0
  do i = 1, 3
    trans(i, i) = 1.d0
  end do
!
!  Make normalized translation vectors
!
  do i = 1, id
    sum = 0.d0
    do j = 1, 3
      sum = sum + tvec(j, i)**2
    end do
    storen(i) = Sqrt (sum)
    sum = 1.d0 / storen(i)
    do j = 1, 3
      xtoc(j, i) = tvec(j, i) * sum
    end do
  end do
!
!   Construct vectors perpendicular to known vectors,
!   if system is of dimension 1 or 2
!
  if (id == 1) then
!
!    Make a vector orthogonal to the first vector
!
    summax = 0.d0
    summin = 1.d0
    do j = 1, 3
      if (Abs (summax) < Abs (xtoc(j, 1))) then
        j1 = j
        summax = xtoc(j, 1)
      end if
      if (Abs (summin) >= Abs (xtoc(j, 1))) then
        j2 = j
        summin = xtoc(j, 1)
      end if
    end do
!
!  If possible, make XTOC(2,2) large.
!
    if (Abs (summax) > 0.80d0 .and. j2 /= 2) then
      j2 = 2
      j1 = 1
      j3 = 3
    else
      j3 = 6 - j1 - j2
    end if
    xtoc(j3, 2) = 0.d0
    xtoc(j2, 2) = summax
    xtoc(j1, 2) = -summin
  end if
  if (id /= 3) then
!
!  Make a vector orthogonal to the first and second vectors.
!
    summin = 1.d0
    do j = 1, 3
      sum = xtoc(j, 1)**2 + xtoc(j, 2)**2
      if (summin > sum) then
        j2 = j
        summin = sum
      end if
    end do
    xtoc(j2, 3) = 1.d0
    do j = 1, 3
      if (j /= j2) then
        sum = 0.d0
        do i = 1, 2
          sum = sum + xtoc(j, i)*xtoc(j2, i)
        end do
        xtoc(j, 3) = -sum
      end if
    end do
    sum = 1.d0 / Sqrt (xtoc(1, 3)**2+xtoc(2, 3)**2+xtoc(3, 3)**2)
    xtoc(:, 3) = xtoc(:, 3) * sum
  end if
!
!   Expand normalized vectors to their full size
!
  do i = 1, id
    do j = 1, 3
      xtoc(j, i) = xtoc(j, i) * storen(i)
    end do
  end do
  do i = id + 1, 3
    do j = 1, 3
      tvec(j, i) = xtoc(j, i)
    end do
  end do
!
!   Convert Cartesian coordinates into fractional unit cell coordinates.
!
!   Add in unit vectors at the end of coord, to be used in making trans or p_rot
  natot = natot + 3
  k = 0
  do i = natot-2, natot
    k = k + 1
    do j = 1, 3
      coord(j, i) = 0.d0
    end do
    coord(k, i) = 1.d0
  end do
  sum = Sqrt (tvec(2, 1)**2 + tvec(3, 1)**2)
  if (sum > 1.d-6) then
!
!    Rotate to eliminate TVEC(3,1)
!
    ca = tvec(2, 1) / sum
    sa = tvec(3, 1) / sum
    do i = 1, 3
      sum1 = tvec(2, i)*ca + tvec(3, i)*sa
      tvec(3, i) = -tvec(2, i)*sa + tvec(3, i)*ca
      tvec(2, i) = sum1
    end do
    do i = 1, natot
      sum1 = coord(2, i)*ca + coord(3, i)*sa
      coord(3, i) = -coord(2, i)*sa + coord(3, i)*ca
      coord(2, i) = sum1
    end do
!
!    Rotate to eliminate TVEC(2,1)
!
    sum = Sqrt (tvec(1, 1)**2 + tvec(2, 1)**2)
    ca = tvec(1, 1) / sum
    sa = tvec(2, 1) / sum
    do i = 1, 3
      sum1 = tvec(1, i)*ca + tvec(2, i)*sa
      tvec(2, i) = -tvec(1, i)*sa + tvec(2, i)*ca
      tvec(1, i) = sum1
    end do
    do i = 1, natot
      sum1 = coord(1, i)*ca + coord(2, i)*sa
      coord(2, i) = -coord(1, i)*sa + coord(2, i)*ca
      coord(1, i) = sum1
    end do
  end if
!
!    Rotate to eliminate TVEC(3,2)
!
  sum = Sqrt (tvec(2, 2)**2 + tvec(3, 2)**2)
  if (sum > 1.d-6) then
    ca = tvec(2, 2) / sum
    sa = tvec(3, 2) / sum
    do i = 2, 3
      sum1 = tvec(2, i)*ca + tvec(3, i)*sa
      tvec(3, i) = -tvec(2, i)*sa + tvec(3, i)*ca
      tvec(2, i) = sum1
    end do
    do i = 1, natot
      sum1 = coord(2, i)*ca + coord(3, i)*sa
      coord(3, i) = -coord(2, i)*sa + coord(3, i)*ca
      coord(2, i) = sum1
    end do
  end if
  natot = natot - 3
  do i = 1, 3
    do j = 1, 3
      trans(i, j) = coord(i, j+natot)
    end do
  end do
  natot = natot + 3
  do i = 1, natot
    coord(3, i) = coord(3, i) / tvec(3, 3)
    coord(2, i) = coord(2, i) - coord(3, i)*tvec(2, 3)
    coord(1, i) = coord(1, i) - coord(3, i)*tvec(1, 3)
    coord(2, i) = coord(2, i) / tvec(2, 2)
    coord(1, i) = coord(1, i) - coord(2, i)*tvec(1, 2)
    coord(1, i) = coord(1, i) / tvec(1, 1)
  end do
  natot = natot - 3
  if (Index (keywrd, " FRACT") /= 0) then
    write(iw_new, "(a)") " Fractional Unit Cell Coordinates"
    write(iw_new, "(i4,3f12.5)") (i, coord(:, i), i = 1, natot)
  end if
  do i = 1, 3
    if (tvec(i, i) < 0) then
      do j = 1, 3
        trans(i, j) = -trans(i, j)
      end do
    end if
  end do
!
!   TRANS is the unitary transform to rotate from the Cartesian frame
!   into the crystallographic frame.
!
  if (Index (keywrd, " TRANS") /= 0) then
    write(iw_new, "(a)") " Unitary Transform: Cartesian to Crystal"
    write(iw_new, "(3f12.5)") (trans(:, i), i = 1, 3)
  end if
!***********************************************************************
!
!   Rotate Secular Determinant into crystal coordinate frame.
!
!***********************************************************************
  do i = 1, 4
    do j = 1, 4
      t1(i, j) = 0.d0
    end do
  end do
  if(phonon)then
    t1(1:3,1:3) = trans(1:3,1:3)
  else
    do i = 1, 3
      do j = 1, 3
        t1(i + 1, j + 1) = trans(i, j)
      end do
    end do
    t1(1, 1) = 1.d0
    do i = 1,3
      do j = 1,3
        p_rot(i,j) = trans(j,i)
      end do
    end do
    p_rot = trans
!
!  "d"-orbital transform
!
    two = Sqrt (0.5d0)  
    six = Sqrt (1/6.d0)  
!
!   d(x^2 - y^2)
!
    t1(5, 5) = (a*a - b*b - d*d + e*e)/2.d0
    t1(5, 6) = a*c - d*f  
    t1(5, 7) = (2.d0*c*c - b*b - a*a - 2.d0*f*f + e*e + d*d) * two * six  
    t1(5, 8) = b*c - e*f  
    t1(5, 9) = a*b - d*e 
!
!  d(xz)
! 
    t1(6, 5) = a*g - b*h  
    t1(6, 6) = a*o + c*g  
    t1(6, 7) = (2.d0*c*o - b*h - a*g) * six/two  
    t1(6, 8) = b*o + c*h  
    t1(6, 9) = a*h + b*g  
!
!  d(2z^2 - x^2 - y^2)
!
    t1(7, 5) = (2.d0*g*g - d*d - a*a - 2.d0*h*h + e*e + b*b) * two * six  
    t1(7, 6) = (2.d0*g*o - d*f - a*c) * six/two
    t1(7, 7) = (4.d0*o*o - 2.d0*(h*h + g*g + f*f + c*c) + e*e + d*d + b*b + a*a)/6.d0 
    t1(7, 8) = (2.d0*h*o - e*f - c*b) * six/two  
    t1(7, 9) = (2.d0*g*h - e*d - a*b) * six/two 
!
!  d(yz)
!
    t1(8, 5) = d*g - e*h  
    t1(8, 6) = d*o + f*g  
    t1(8, 7) = (2.d0*o*f - e*h - d*g) * six/two  
    t1(8, 8) = e*o + f*h  
    t1(8, 9) = d*h + e*g  
!
!  d(xy)
!
    t1(9, 5) = a*d - b*e  
    t1(9, 6) = a*f + c*d  
    t1(9, 7) = (2.d0*c*f - e*b - a*d) * six/two  
    t1(9, 8) = b*f + c*e  
    t1(9, 9) = a*e + d*b
!
    do i = 5,9
      do j = 5,i
        sum = t1(i,j)
        t1(i,j) = t1(j,i)
        t1(j,i) = sum
      end do
    end do
    if (.false.) then
      write(*,*)"  p_rot"
      write(*,'(3f16.10)')((p_rot(i,j), i = 1,3),j = 1,3)
      write(*,*)
      do i = 1, 3
        do j = 1, 3
          sum = 0.d0
          sum1 = 0.d0
          do k = 1, 3
            sum = sum + p_rot(i,k)*p_rot(j,k)
            sum1 = sum1 + p_rot(k,i)*p_rot(k,j)
          end do
          e1(i,j) = sum
          e2(i,j) = sum1
        end do
      end do
      write(*,*) " p_rot squared"
      write(*,'(3f16.10)')((e1(i,j), i = 1,3),j = 1,3)
      write(*,*)
      write(*,'(3f16.10)')((e2(i,j), i = 1,3),j = 1,3)
      write(*,*)
      do i = 5, 9
        do j = 5, 9
          sum = 0.d0
          sum1 = 0.d0
          do k = 5,9
            sum = sum + t1(i,k)*t1(j,k)
            sum1 = sum1 + t1(k,i)*t1(k,j)
          end do
          e1(i,j) = sum
          e2(i,j) = sum1
        end do
      end do
      write(*,*) " d_rot "
      write(*,'(5f16.10)')((t1(i,j), i = 5,9),j = 5,9)
      write(*,*)
      write(*,*) " d_rot squared"
      write(*,'(5f16.10)')((e1(i,j), i = 5,9),j = 5,9)
      write(*,*)
      write(*,'(5f16.10)')((e2(i,j), i = 5,9),j = 5,9)
    end if
  end if
  nprt = Min (8, nvecs)
  prtok = (Index (keywrd, " ROTSEC")/=0)
  if (prtok) then
    write(iw_new, "(/2X,'SECULAR DETERMINANT IN CRYSTAL ORIENTATION',/)")
    i = Index (keywrd, " ROTSEC=")
    mprt = 2000
    if (i /= 0) then
      mprt = Nint (reada (keywrd, i+3))
    end if
  end if
!
!  Rotate entire secular determinant into crystal orientation
!
  ncell = 0
  do i = 1, mr1
    do j = 1, mr2
      do k = 1, mr3
        if (.not.bcc .or. Mod(i + j + k, 2) == 1) then
          iofset = nvecs * ncell
          ncell = ncell + 1
          do ii = 1, numat
            il = nfirst(ii)
            iu = nlast(ii)
            do jj = 1, numat
              jl = nfirst(jj)
              ju = nlast(jj)
              do i1 = 1, per_atom(ii)
                do j1 = 1, per_atom(jj)
                  sum = 0.d0
                  k2 = 0
                  do k1 = il, iu
                    k2 = k2 + 1
                    l2 = 0
                    do l1 = jl, ju
                      l2 = l2 + 1 
                      sum = sum + t1(i1, k2)*sec_det(k1, l1 + iofset)*t1(j1, l2)
                    end do
                  end do
                  t2(i1, j1) = sum
                end do
              end do
              i2 = 0
              do i1 = il, iu
                i2 = i2 + 1
                j2 = 0
                do j1 = jl, ju
                  j2 = j2 + 1
                  sec_det(i1, j1 + iofset) = t2(i2, j2)
                end do
              end do
            end do
          end do
          if (prtok .and. ncell <= mprt) then
            write(iw_new, "(A,3I3)") " Unit Cell:", i-1, j-1, k-1
            do jj = iofset + 1, iofset + nvecs
              write(iw_new, "(10f8.4)") (sec_det(ii, jj), ii = 1, nprt)
            end do
          end if
        end if
      end do
    end do
  end do
  prtops = (Index (keywrd, " OPS")/=0)
!  i = Index (jobnam, " ") - 1
!
!
!                                   Read in symmetry operations
!
!
  open (unit = 4, file = trim(jobnam)//".ops", status = "OLD", iostat = i99)
  i = 1
  if (i99 == 0) then
    if (prtops) then
      write(iw_new, "(//11X,'SYMMETRY OPERATIONS OF SPACE GROUP',//)")
      write(iw_new, '( "    INV.  NON-PRIMITIVE ROTATION   AXIS OF               CENTER          NAME", &
        & / "    1=YES TRANSLATIONS   ANGLE    ROTATION           OF OPERATION", &
        & /, "             X   Y   Z", 12x, "X     Y     Z       X     Y     Z "/)')
    end if
    do i = 1, 149
      do
        read(4,'(a)') line
        if (line(1:1) /= "*") exit
      end do
      read (line, *, end = 1000, err = 1000) i1, i2, i3, i4, &
           & (xop(j, i), j = 5, 8), xop(9, i), xop(10, i), xop(11, i), title(i)
      xop(1, i) = i1
      if (Abs (xop(6, i)) + Abs (xop(7, i)) < 1.d-5) then
        xop(8, i) = 1.d0
      end if
      if (i1 == 2) exit
      x = 0.d0
      if (i2 /= 0) x = 1.d0 / i2
      y = 0.d0
      if (i3 /= 0) y = 1.d0 / i3
      z = 0.d0
      if (i4 /= 0) z = 1.d0 / i4
      xop(2, i) = x
      xop(3, i) = y
      xop(4, i) = z
      j = Nint (xop(5, i)*360.d0)
      xop(5, i) = j
      if (prtops) write(iw_new, '(i2, i4, f9.1, 2f4.1, i6, f8.2, 2f6.2, f8.2, 2f6.2, 3x, a12)') &
          i, i1, (xop(j, i), j = 2, 4), Int (xop(5, i)+0.5d0), (xop(j, i), j = 6, 11), title(i)
      xop(5, i) = xop(5, i) / 360.d0
    end do
  end if
  1000 iop = i - 1
  if (iop == 0) then
!
!  No symmetry operations were supplied, so force in the Identity
!
         iop = 1
         do i = 1, 11
           xop(i, 1) = 0.d0
         end do
         xop(8, 1) = 1.d0
         title(1) = "Identity"
       end if
!***********************************************************************
!
!     Symmetrize Secular Determinant
!
!***********************************************************************
  prtok = (Index (keywrd, " SYMSEC") /= 0)
  if (iop /= 1 .and. prtok) then
    write(iw_new, *) " Symmetrizing Secular Determinant."
  end if
  call symtrz (sec_det, nvecs, coord)
  if (.not. bcc) call add_missing_cells()
  prtok = (Index (keywrd, " SYMSEC") /= 0)
  if (.not. prtok) return
  write(iw_new, "(/2X,'SYMMETRIZED SECULAR DETERMINANT',/)")
  i = Index (keywrd, " SYMSEC=")
  mprt = 2000
  if (i /= 0) mprt = Nint (reada (keywrd, i+3))
  mprt = Min (ncell, mprt)
  nprt = Min (8, nvecs)
  do ii = 1, mprt
    write(iw_new, "(/3X,'Unit Cell:',3I3)") (nijk(j, ii), j = 1, 3)
    do i = (ii-1)*nvecs + 1, ii*nvecs
      write(iw_new, "(10F8.4)") (sec_det(jj, i), jj = 1, nprt)
    end do
  end do
end subroutine orient
!
