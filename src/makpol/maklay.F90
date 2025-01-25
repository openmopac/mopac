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

subroutine maklay (labels, coord, geo, na, nb, nc, mers, loc, nmax)
    use common_systm, only : maxtxt, nalpha, numat, id, iw, l_int, natoms, nvar, &
      line, ndep, l_use_tv
    use common_keywrd, only : keywrd
    use common_ctvec, only : tvec
    use common_jobnam, only : jobnam
    use common_sizes, only : numatm
    
    implicit none
    integer, intent (in) :: nmax
    integer, dimension (3), intent (in) :: mers
    integer, dimension (nmax), intent (inout) :: labels, na, nb, nc
    integer, dimension (2, 3*nmax), intent (inout) :: loc
    double precision, dimension (3, nmax), intent (inout) :: coord, geo
    logical :: debug, lbcc, sym, exists
    integer :: i, ioff = 0, j, k, l, loop1, loop2, loop3, ncell1
    double precision :: degree,  r, xparam(3*numatm)
    double precision, external :: snapth
    debug = (Index (keywrd, " DEBUG") /= 0)
    lbcc = (Index (keywrd, " BCC") /= 0)
    sym = (Index (keywrd, " SYMM") /= 0) 
    if (debug) then
      write(iw,"(/5x,A)")"At entry to MAKLAY"
      write(iw,'(/5x,a, i4)')"Number of atoms:",numat
      write(iw,'(5x,a, i4)')"Number of Tv:   ",id
      write(iw,'(5x,a,3i4)')"MERS:           ",mers(:id)
      write(iw,'(/15x,a,9F9.4)')"TVEC           "
      write(iw,'(2x,3f9.4)')tvec(:,:id)
      write(iw,'(/10x,a,/)')"Geometry (atom number, label, coordinates"      
      do i = 1, numat
        write (iw, "(2I6,3F18.12,3I4)") i, labels (i), (coord(k, i), k = 1, 3)
      end do
    end if
    nalpha = numat
    do loop1 = 0, mers(3) - 1
      do loop2 = 0, mers(2) - 1
        do loop3 = 0, mers(1) - 1
!
!   Find unit cell for relating current unit cell to
!
          if (loop3 /= 0 .or. loop2 /= 0 .or. loop1 /= 0) then
            if ( .not. lbcc .or. Mod (loop1+loop2+loop3, 2) /= 1) then
!
!  All the atoms in the unit cell LOOP1, LOOP2, LOOP3
!
              ioff = ioff + numat
              do j = 1, numat
                labels(ioff+j) = labels(j)
                do k = 1, 3
                  coord(k, ioff+j) = coord(k, j) + tvec(k, 3) * loop1 + &
                 & tvec(k, 2) * loop2 + tvec(k, 1) * loop3
                end do
                k = j + ioff - numat
                if (na(k) /= 0 .and. nb(k) /= 0 .and. nc(k) /= 0) then
                  na(ioff+j) = na(k) + numat
                  nb(ioff+j) = nb(k) + numat
                  nc(ioff+j) = nc(k) + numat
                end if
              end do
            end if
          end if
        end do
      end do
    end do
!
!  Convert from
!
!  All atoms in Unit cell 1
!  All atoms in Unit cell 2
!  All atoms in Unit cell
!
!  to
!
!  Atom 1 (all unit cells)
!  Atom 2 (all unit cells)
!  Atom 3 (all unit cells)
!
    i = Index (keywrd, " SORT")
    if (i /= 0) then
!
!  Delete keywords SORT and MERS
!
      keywrd (i:i+4) = " "
      i = Index (keywrd, " MERS")
      j = Index (keywrd(i+1:), " ") + i
      keywrd (i:j) = " "
      ncell1 = ioff / numat + 1
      do i = 1, numat * ncell1
        na(i) = labels(i)
        do j = 1, 3
          geo(j, i) = coord(j, i)
        end do
      end do
      ncell1 = ncell1 - 1
      l = 0
      do i = 1, numat
        do j = 0, ncell1
          l = l + 1
          labels(l) = na(i)
          do k = 1, 3
            coord(k, l) = geo(k, j*numat+i)
          end do
        end do
      end do
      do i = 1, l
        na(i) = 0
        nb(i) = 0
        nc(i) = 0
      end do
    end if
    numat = ioff + numat     
!
!  Increase size of Tv by factor MERS
!
    do i = 1, id
      do j = 1, 3
        tvec(j, i) = tvec(j, i) * mers(i) 
      end do
    end do
    do i = 2, numat
      do j = 1, i - 1
        r = (coord(1, i)-coord(1, j))**2 + (coord(2, i)-coord(2, j))**2 + &
      & (coord(3, i)-coord(3, j))**2
        if (r < 0.2d0 .and. labels(i) < 99 .and. labels(j) < 99) labels(i) = 200
      end do
    end do
    j = 1
    do i = 2, numat
      if (labels(i) /= 200) then
        j = j + 1
        coord(1,j) = coord(1,i)
        coord(2,j) = coord(2,i)
        coord(3,j) = coord(3,i)
        labels(j) = labels(i)
      end if
    end do
    numat = j
    if (l_int) then
!
!  Add in dummy atoms for translation vectors to use
!
      do i = 1, id
        numat = numat + 1
        labels(numat) = 99
        do j = 1, 3
          coord(j, numat) = tvec(j, i)
        end do
      end do  
      do i = 1, numat
        if (labels(i) == 98) then
          labels(i) = 99
        end if
      end do

      do i = 1, numat
        if (labels(i) == 99) then
          labels(i) = 98
        end if
      end do
!
!  Convert back to internal coordinates
!
      degree = 1.d0
      na = 0; nb = 0; nc = 0;
      call xyzint (coord, na, nb, nc, degree, geo)
!
!   If any angles are near to important angles (such as 109.47...)
!   snap the angle to the exact angle
!
      do i = 1, numat
        if (na(i) /= 0) then
          geo(2, i) = snapth (geo(2, i))
          geo(3, i) = snapth (geo(3, i))
        end if
      end do
!
!   All real atoms plus "ID" dummy atoms are now in GEO, in
!   internal coordinates.
!
!
!  Now to add the translation vectors
!
      do i = 1, id
        numat = numat + 1
        labels(numat) = 107
        geo(1, numat) = Sqrt (tvec(1, i)**2+tvec(2, i)**2+tvec(3, i)**2)
        geo(2, numat) = 0.d0
        geo(3, numat) = 0.d0
        na(numat) = 1
        nb(numat) = numat - id
        nc(numat) = 2
      end do
      go to 1100
    else
      do i = 1, numat
        do j = 1, 3
          geo(j, i) = coord(j, i)
        end do
      end do
!
!  Add in Cartesian translation vectors
!
      do i = 1, id
        numat = numat + 1
        labels(numat) = 107
        do j = 1, 3
          geo(j, numat) = tvec(j, i)
        end do
      end do
    end if
1100 natoms = numat
!
!   Make parameters from all possible geometric variables.
!
    nvar = 0
    do i = 1, numat
      if (l_int) then
        k = Min (3, i-1)
!
!  For Tv, only use "bond length" - do NOT optimize angles or dihedrals
!
        if (labels(i) == 107) then
          k = 1
        end if
      else
        k = 3
      end if
      do j = 1, k
        nvar = nvar + 1
        xparam(nvar) = geo(j, i)
        loc(1, nvar) = i
        loc(2, nvar) = j
      end do
    end do
    if (sym) then
      call maksym (loc, xparam)
    end if
    if (l_int) then
      i = index(keywrd, "     ")
      keywrd (i:i + 3) = " INT"
    end if
!
!   Unprotect any dummy atoms
!
    do i = 1, numat
      if (labels(i) == 98) then
        labels(i) = 99
      end if
    end do
    j = 0
    do i = 1, numat
      if (na(i) /= 0) j = 1
    end do
    if (.false. .and. j == 0) then
!
!  Put atom 1 at the origin
!
      do i = 2, numat - 3
        geo(:,i) = geo(:,i) - geo(:,1)
      end do
      geo(:,1) = 0.d0
    end if
    maxtxt = 0
    if (debug) then
      write(iw,'(///10x,a,///)')"MOPAC data-set"
    end if
    i = index(keywrd, " INT")
    if (i > 0) keywrd(i:i + 3) = " "
    call geout (iw, geo, na, nb, nc, labels, loc)
    line = jobnam
 !   do i = 1, len_trim(line)
 !     if (line(i:i) == "_") line(i:i) = " "
 !   end do
    line = "Make_"//trim(line)//".dat"
    call add_path(line)
    inquire (file = trim(line), exist = exists)
    if ( .not. exists .and. l_use_tv) then
      open (unit=15, file=trim(line), status="NEW")
      ndep = 0
      k = index(keywrd, "BCC")
      if (k /= 0) then
        keywrd(k:k + 2) = " "
        i = index(keywrd, "MERS")
        if (i /= 0) then
          j = index(keywrd(i:), " ") + i
          write(keywrd(i:j),'(a,i1,a,i1,a,i1,a)')"MERS=(", mers(1)/2, ",", mers(2)/2, ",", mers(3)/2, ")"
        end if           
      end if
      call geout (15, geo, na, nb, nc, labels, loc)
    end if
    stop
end subroutine maklay
