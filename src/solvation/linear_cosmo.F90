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

!**************************************************************************
!
!  Based on the work of:
!
!  Dr Andrey Bliznuk, Australian National University, Canberra, Australia.
!
!  See also:
!  Jana Khandogin, Anguang Hu, Darrin M. York, "Electronic structure properties
!  of solvated biomolecules: A quantum approach for macromolecular characterization
!  Volume 21, Issue 16 , Pages 1562 - 1571 (2000)
!
!**************************************************************************
  module cosmo_mini
!
!  This structure, the use of a module inside a module containing private data,
!  is necessary for debugging only.  The data here can be monitored by a debugger.
!  The other data in linear_cosmo cannot be monitored.
!
    logical :: new_surface, new_iteration
    double precision, dimension(:, :), allocatable :: a_block
    double precision, dimension(:), allocatable :: r_vec, p_vec, &
       & q_vec, z_vec, a_diag, m_vec, a_part
    double precision, dimension(:, :, :), allocatable :: tm
    integer, dimension(:), allocatable :: iblock_pos
  end module cosmo_mini
  module linear_cosmo
  use afmm_C
  use cosmo_mini, only : new_surface, new_iteration, a_block, tm, &
  r_vec, p_vec, q_vec, z_vec, a_diag, m_vec, a_part, iblock_pos
  implicit none

  public :: ini_linear_cosmo, coscavz, addnucz, addfckz, am1dft_solve
  double precision, public :: c_proc


  private

  integer, parameter :: na1max = 6000     ! when use O(n*log(n))
  integer, parameter :: na2max = 6000     ! when use O(n) algorithm
                                          ! for A*q vector

  integer, parameter :: nb1max = 8000     ! when use O(n*log(n))
  integer, parameter :: nb2max = 8000     ! when use O(n) algorithm
                                          ! for B*q vector

  logical, parameter :: compute_a_part = .true. ! precompute A matrix
  logical, parameter :: use_a_blocks = .true.   ! use block-diagonal
                                                ! preconditioner

  !!!!
  !!!! Note. use_a_blocks assumes that compute_a_part is .true.
  !!!!

  double precision, dimension(:, :), allocatable :: rsc

  integer, dimension(:), allocatable :: nipsrs
  integer, dimension(:), allocatable :: nset
  integer, dimension(:), allocatable :: npoints
  integer, dimension(:), allocatable :: iatom_pos

  integer, dimension(:), allocatable :: ijbo_diag



  integer :: max_block_size
  integer :: atom_handle, surface_handle
  integer, parameter :: dcmplx_kind = Kind ((1.d0, 1.d0))

contains

  subroutine ini_linear_cosmo

    use molkst_C, only: numat
    use cosmo_C, only : lenabc, solv_energy
    use common_arrays_C, only: nfirst, nlast
    implicit none
    integer :: maxrs, i, nnn
    intrinsic Allocated
    integer, external :: ijbo
    !
    !       Allocate Linear-Cosmo specific arrays
    !
    maxrs = 70 * numat
    if (Allocated (rsc)) then
      deallocate (rsc, tm, iatom_pos, nipsrs, nset, npoints, &
                & ijbo_diag, stat = i)

      if (i /= 0) then
        call mopend ("CosmoZ (1):  Deallocate error")
      end if
    end if

    allocate (rsc (4, maxrs), tm(3, 3, numat), iatom_pos(numat), &
         & nipsrs(lenabc), nset(1082*numat), npoints(numat+1), &
         ijbo_diag(numat), stat=i)

    if (i /= 0) then
      call mopend ("CosmoZ (2):  Allocate error ")
    end if

    iatom_pos(1) = 0
    ijbo_diag(1) = ijbo(1, 1)
    do i = 2, numat
      nnn = nlast(i-1) - nfirst(i-1) + 1
      nnn = ((nnn+1)*nnn)/2
      iatom_pos(i) = iatom_pos(i-1) + nnn
      ijbo_diag(i) = ijbo(i, i)
    end do

    solv_energy = 0.d0

    call afmm_ini

  end subroutine ini_linear_cosmo

  subroutine bpnew_vec (v)
    !
    ! Interaction of electronic densities with induced charges
    !

    use molkst_C, only : numat
    use cosmo_C, only : nps, cosurf
    use common_arrays_C, only: coord, nfirst, nlast, nat, p
    use parameters_C, only : dd, qq
    use funcon_C, only : a0
    implicit none
    double precision, dimension(nps), intent (out) :: v
    integer :: i, j, k, l, m, jj, nao, ni
    double precision :: dx, dy, dz, r, t, dip, quad
    double precision :: w(45)
    v(1:nps) = 0.d0
    dip = 0.d0
    quad = 0.d0

    do j=1, numat
      nao = nlast(j)-nfirst(j)
      jj = ijbo_diag(j )

      if (nao > 0) then
        ni = nat(j)
        dip = dd(ni) * a0
        quad = (a0 *qq(ni))** 2
      end if

      do i=1, nps

        dx= cosurf(1, i) - coord(1, j)
        dy= cosurf(2, i) - coord(2, j)
        dz= cosurf(3, i) - coord(3, j)

        r = 1.d0 / Sqrt (dx**2 + dy**2 + dz**2)

        if (nao == 0) then ! --- S-element
          v(i) = v(i) - p(jj+1) * r
          cycle

        else ! --- sp element

           w(1) = r

           dx = dx * r
           dy = dy * r
           dz = dz * r

           t = dip * r*r

           w(2) =  dx * t
           w(4) =  dy * t
           w(7) =  dz * t

           t = quad * r*r*r
           w(3)  = r + (3.d0*dx**2 - 1.d0) * t
           w(6)  = r + (3.d0*dy**2 - 1.d0) * t
           w(10) = r + (3.d0*dz**2 - 1.d0) * t

           w(5) = 3.d0*dx*dy * t
           w(8) = 3.d0*dx*dz * t
           w(9) = 3.d0*dy*dz * t

           if (nao > 3) then ! --- d - element
              w(11:45) = 0
              w(15) = r
              w(21) = r
              w(28) = r
              w(36) = r
              w(45) = r
           end if
        end if

        m = 0
        t = 0.d0
        do k = 1, nao+1
          do l = 1, k-1
            m = m + 1
            t = t - 2*p(jj+m) * w(m)
          end do
          m = m + 1
          t = t - p(jj+m) * w(m)
        end do

        !!write(26,*) ' Pair :', i, j
        !!write(26,*) ( w(k),k=1,((nao+1)*(nao+2))/2 )

        v(i) = v(i) + t
      end do
    end do

  end subroutine bpnew_vec

  subroutine bz_vec (v)
    !
    ! Interaction of atom cores with induced charges
    !

    use molkst_C, only : numat
    use common_arrays_C, only: coord, nat
    use cosmo_C, only : cosurf, nps
    use parameters_C, only: tore
    implicit none
    double precision, dimension(nps), intent (out) :: v
    integer :: i, j
    double precision :: t
    intrinsic Sqrt


    do i = 1, nps
      t = 0
      do j = 1, numat

        t = t + tore(nat(j)) / &
             & Sqrt ((coord(1, j) - cosurf(1, i))**2 + &
             & (coord(2, j) - cosurf(2, i))**2 + &
             & (coord(3, j) - cosurf(3, i))**2)
      end do

      v(i) = t

    end do
  end subroutine bz_vec

  subroutine get_bvec (x1, x2, nao, ni, w)
    use funcon_C, only : a0
    use parameters_C, only : dd, qq
    implicit none
    double precision, dimension(3), intent (in) :: x1, x2
    integer, intent (in) :: nao, ni
    double precision, dimension(*), intent (out) :: w
    double precision, dimension(3) :: dx
    double precision :: r, dip, quad
    integer :: i
    intrinsic Sqrt
    !
    dx(1) = x1(1) - x2(1)
    dx(2) = x1(2) - x2(2)
    dx(3) = x1(3) - x2(3)

    r = 1.d0 / Sqrt (dx(1)**2 + dx(2)**2 + dx(3)**2)

    w(1) = r

    if (nao == 0) then ! --- H
      return
    else ! --- sp - element
      dip = dd(ni) * a0 * r * r
      quad = (a0 *qq(ni))** 2 * r * r * r
      dx(1) = dx(1) * r
      dx(2) = dx(2) * r
      dx(3) = dx(3) * r

      w(2) =  dx(1) * dip
      w(4) =  dx(2) * dip
      w(7) =  dx(3) * dip

      w(3)  = r + (3.d0*dx(1)**2 - 1.d0) * quad
      w(6)  = r + (3.d0*dx(2)**2 - 1.d0) * quad
      w(10) = r + (3.d0*dx(3)**2 - 1.d0) * quad

      w(5) = 3.d0*dx(1)*dx(2) * quad
      w(8) = 3.d0*dx(1)*dx(3) * quad
      w(9) = 3.d0*dx(2)*dx(3) * quad

      if (nao < 4) then
        return
      else ! --- d - element
        do i = 11, 45
          w(i) = 0
        end do
        w(15) = r
        w(21) = r
        w(28) = r
        w(36) = r
        w(45) = r
      end if
    end if
  end subroutine get_bvec

  subroutine coscavz (coord, nat)
    !
    ! Prepares molecular surface & tesselations
    !

    use molkst_C, only: numat
    use overlaps_C, only : cutof2
    use cosmo_C, only : nps, disex2, ioldcv, lenabc, iatsp, nar_csm, nsetf, phinet, &
         & cosurf, srad, arat
    use chanel_C, only: iw
    implicit none
    double precision, dimension(3, numat), intent (inout) :: coord
    integer, dimension(numat), intent (in) :: nat
    integer :: i, j, ndim, maxrs
    integer, dimension(2) :: ind
    new_surface = .true.

    call coscanz (srad, coord, nat, cosurf, iatsp, nar_csm, nsetf, &
         & phinet, arat)

 !   write(26,*) ' Cosmo radii:'
 !   write(26,'(F12.5)') (srad(i), i=1, numat)

    call prepare_tesselations (numat, nps, ind, 2, i)

    if (i /= 0) then
      call mopend ("CosmoZ (3):  Tesselation ini error ")
      return
    end if

    atom_handle = ind(1)
    surface_handle = ind(2)

    if (numat > nb2max .or. numat > nb1max) then
 !     write (iw, *) " NUMAT Doing divide_box by", Sqrt(cutof2)
      call divide_box (coord, 3, numat, Sqrt(cutof2), 10, 0, atom_handle, i)
    end if

    if (i /= 0) then
      call mopend("Error in divide box by atoms ")
      return
    end if

 !   write (iw, *) " Number of surface elements = ", nps

    if (nps > na2max .or. nps > na1max) then
  !    write (iw, *) " NPS Doing divide_box by", Sqrt(disex2)
      call divide_box (cosurf, 4, nps, Sqrt(disex2), 20, 0, &
           & surface_handle, i)
    end if

    if (i /= 0) then
      call mopend("Error in divide box by surface ")
      return
    end if

    c_proc = 0

    if (compute_a_part) then

      if (allocated (a_part)) then
        deallocate (a_part)
      end if
      if (nps > na2max .or. nps > na1max) then

        call set_tesselation (surface_handle, i)
        if (i /= 0) then
          call mopend ("CosmoZ:  Tesselation ini error ")
          return
        end if
        ndim = count_short_ints  (cosurf, 4, simulate_aq_dir_int, .false.)

      else

        maxrs = 60*numat

        call simulate_aq_vec (coord, srad, numat, cosurf, nps, nar_csm, &
       & nsetf, nset, rsc, nipsrs, iatsp, tm, ioldcv, maxrs, lenabc, &
       & .false., ndim)

      end if

  !    write(iw,*) ' Number of close-range interactions: ', ndim

      allocate (a_part(ndim) , stat=i)
      if (i /= 0) then
         write(iw,*) ' Allocation failed for a_part'
         call mopend (' Allocation failed for a_part')
         return
      end if

    end if

    if (allocated (m_vec)) then
      deallocate (m_vec)
    end if

    if (use_a_blocks) then
      if (allocated (iblock_pos)) then
        deallocate (iblock_pos)
        deallocate (a_block)
      end if

      allocate (iblock_pos(numat) , stat=i)
      if (i /= 0) then
         write(iw, *) ' Allocation failed for iblock_pos'
         call mopend (' Allocation failed for iblock_pos')
         return
      end if

      ndim = 0
      max_block_size = 0
      do i=1, numat
        iblock_pos(i) = ndim +1
        j = npoints(i+1) - npoints(i)
        if (j > max_block_size) then
          max_block_size = j
        end if
        ndim = ndim + (j*(j+1)) / 2
      end do

  !    write(iw, *) ' Size of the preconditioner = ', ndim
  !    write(iw, *) ' Size of max_block = ', max_block_size

      allocate (m_vec(ndim) , stat=i)
      if (i /= 0) then
         write(iw, *) ' Allocation failed for m_vec'
         call mopend (' Allocation failed for m_vec')
         return
      end if
      allocate (a_block(max_block_size, max_block_size) , stat=i)
      if (i /= 0) then
         write(iw, *) ' Allocation failed for a_block'
         call mopend (' Allocation failed for a_block')
         return
      end if

      a_block = 0

    else

      write(iw, *) ' Size of the preconditioner = ', nps

      allocate (m_vec(nps) , stat=i)
      if (i /= 0) then
         write(iw, *) ' Allocation failed for m_vec'
         call mopend (' Allocation failed for m_vec')
         return
      end if

    end if

  end subroutine coscavz

  subroutine coscanz (srad, coord, nat, cosurf, iatsp, nar_loc, nsetf, &
       & phinet, arat)
    !***********************************************************************
    !
    ! THIS ROUTINE CONSTRUCTS OR UPDATES THE SOLVENT-ACCESSIBLE
    ! SURFACE (SAS)
    ! (Linear-scaling version of the coscan routine from Cosmo)
    !
    !***********************************************************************
    use molkst_C, only: numat
    use cosmo_C, only : lenabc, nps, area, cosvol, rsolv, ioldcv, dirvec, dirsm, &
    n0, isude, sude
    use chanel_C, only: iw
    implicit none
    integer, dimension(lenabc+1), intent (inout) :: iatsp, nar_loc, nsetf
    integer, dimension(numat), intent (in) :: nat
    double precision, dimension(numat), intent (in) :: srad
    double precision, dimension(numat), intent (out) :: arat
    double precision, dimension(lenabc+1, 3), intent (inout) :: phinet
    double precision, dimension(3, numat), intent (inout) :: coord
    double precision, dimension(4, lenabc), intent (inout) :: cosurf
    integer :: i, i0, ik, ilipa, inset, ipm, ips, j, jmax, jps, k, &
         & l, nara, narea, niter, nps0, maxrs, alloc_stat, indi
    double precision :: dist, dist1, dist2, dist3, r, ri, &
         & ri2, rr, sininv, sp, spm, x1, x2, x3, dists

    double precision :: xmax, ymax, zmax, xmin, ymin, zmin, r_max, o_2r
    integer :: nx, ny, nz, ix, iy, iz, i1, i2, l1, k1, ii, jj, &
         & ix1, iy1, iz1, ixfrom, ixto, iyfrom, iyto, izfrom, izto
    integer, dimension(:), allocatable :: atoms_in_box, s_box, box_number, &
         & temp
    logical, allocatable, dimension(:) :: din
    integer, dimension(1082) :: iseg
    double precision, dimension(3) :: xa, xx
    double precision, dimension(3, 1082) :: dirtm
    integer, dimension(:), allocatable :: isort, ipsrs, nipa, lipa
    integer, dimension(:, :), allocatable :: nn
    intrinsic Cos, Sqrt, Min, Max, Int, Allocated
    !
    ! MAKE COORDINATES A BIT ASYMMETRIC IN ORDER TO AVOID
    ! SYMMETRY PROBLEMS WITH CAVITY CONSTRUCTION

    do i = 1, numat
      do j = 1, 3
        coord(j, i) = coord(j, i) + Cos (i*j*.1d0) * 3.0d-9
      end do
    end do

    !
    !       Determing how many boxes are there
    !

    xmax = coord(1, 1)
    ymax = coord(2, 1)
    zmax = coord(3, 1)
    xmin = coord(1, 1)
    ymin = coord(2, 1)
    zmin = coord(3, 1)
    r_max = srad(1)

    do i = 2, numat
      if (coord(1, i) > xmax) then
        xmax = coord(1, i)
      elseif (coord(1, i) < xmin) then
        xmin = coord(1, i)
      end if

      if (coord(2, i) > ymax) then
        ymax = coord(2, i)
      elseif (coord(2, i) < ymin) then
        ymin = coord(2, i)
      end if

      if (coord(3, i) > zmax) then
        zmax = coord(3, i)
      elseif (coord(3, i) < zmin) then
        zmin = coord(3, i)
      end if

      if (srad(i) > r_max) then
        r_max = srad(i)
      end if
    end do

    r_max = r_max + rsolv
    r_max = r_max + r_max   ! cell size

    o_2r = 1 / r_max

    nx = Int ((xmax - xmin) * o_2r)
    if (nx < 1) nx = 1
    if (xmin + nx * r_max < xmax) nx = nx + 1

    ny = Int ((ymax - ymin) * o_2r)
    if (ny < 1) ny = 1
    if (ymin + ny * r_max < ymax) ny = ny + 1

    nz = Int ((zmax - zmin) * o_2r)
    if (nz < 1) nz = 1
    if (zmin + nz * r_max < zmax) nz = nz + 1

    maxrs = 100 * numat
    allocate (atoms_in_box(0:nx*ny*nz), s_box(0:numat-1), &
         & box_number(numat), nn(3, numat), isort(maxrs), &
         & ipsrs(maxrs), nipa(numat), din(Max (numat, 1082)), &
         & stat=alloc_stat)

    if (alloc_stat /= 0) then
      write (iw, *) " Can't allocate memory in COSCANZ"
      go to 100
    end if

    atoms_in_box = 0

    !
    !       Assign atoms to boxes. Boxes are numbered from 0 to n-1
    !

    do i = 1, numat
      ix = Int ((coord(1, i) - xmin) * o_2r)
      if (ix < 0) then
        ix = 0
      elseif (ix >= nx) then
        ix = nx - 1
      end if

      iy = Int ((coord(2, i) - ymin) * o_2r)
      if (iy < 0) then
        iy = 0
      elseif (iy >= ny) then
        iy = ny - 1
      end if

      iz = Int ((coord(3, i) - zmin) * o_2r)
      if (iz < 0) then
        iz = 0
      elseif (iz >= nz) then
        iz = nz - 1
      end if

      box_number(i) = ix + (iy + iz *ny) *nx
      atoms_in_box(box_number(i)) = atoms_in_box(box_number(i)) + 1
    end do

    !
    !       Now we need to collect atoms in the same box together.
    !       The easiest way to do it is by sorting ...
    !

    do i = 0, nx*ny*nz - 1
      atoms_in_box(i+1) = atoms_in_box(i+1) + atoms_in_box(i)
    end do

    do i = 1, numat
      atoms_in_box(box_number(i)) = atoms_in_box(box_number(i)) - 1
      s_box(atoms_in_box(box_number(i))) = i
    end do

    !
    !       It is time to intersect ...
    !
    ilipa = 0
    k = -ny
    do iz = 0, nz - 1
      k = k + ny
      do iy = 0, ny - 1
        l = (iy + k) *nx
        do ix = 0, nx - 1
          i1 = atoms_in_box(l)
          l = l + 1
          i2 = atoms_in_box(l)
          do ii = i1, i2 - 1
            i = s_box(ii)
            ri = srad(i)
            r = ri + rsolv
            rr = r + rsolv
            xa(1:3) = coord(1:3, i)
            izfrom = Max (iz-1, 0)
            izto   = Min (iz+2, nz) - 1
            k1 = (izfrom-1) * ny
            do iz1 = izfrom, izto
              iyfrom = Max (iy-1, 0)
              iyto   = Min (iy+2, ny) - 1
              k1 = k1 + ny
              do iy1 = iyfrom, iyto
                l1 = (iy1 + k1)*nx
                ixfrom = Max (ix-1, 0)
                ixto   = Min (ix+2, nx) - 1
                l1 = l1 + ixfrom - 1
                do ix1 = ixfrom, ixto
                  l1 = l1 + 1
                  do jj = atoms_in_box(l1), atoms_in_box(l1+1) - 1
                    j = s_box(jj)
                    if (j /= i) then
                      dist = (xa(1) - coord(1, j))**2 + &
                           & (xa(2) - coord(2, j))**2 + &
                           & (xa(3) - coord(3, j))**2
                      if (dist < (rr+srad(j))**2) then
                        ilipa = ilipa + 1
                      end if
                    end if
                  end do
                end do   ! ix1
              end do   ! iy1
            end do   ! iz1
          end do
        end do   ! ix
      end do   ! iy
    end do   ! iz
    allocate (temp(ilipa), lipa(ilipa), stat=alloc_stat)
    if (alloc_stat /= 0) then
      write (iw, *) " Can't allocate memory in COSCANZ"
      go to 100
    end if
   !
    ilipa = 0
    k = -ny
    do iz = 0, nz - 1
      k = k + ny
      do iy = 0, ny - 1
        l = (iy + k) *nx
        do ix = 0, nx - 1
          i1 = atoms_in_box(l)
          l = l + 1
          i2 = atoms_in_box(l)

          do ii = i1, i2 - 1
            i = s_box(ii)

            nipa(i) = 0
            ri = srad(i)
            r = ri + rsolv
            rr = r + rsolv
            xa(1:3) = coord(1:3, i)

            izfrom = Max (iz-1, 0)
            izto   = Min (iz+2, nz) - 1
            k1 = (izfrom-1) * ny
            do iz1 = izfrom, izto
              iyfrom = Max (iy-1, 0)
              iyto   = Min (iy+2, ny) - 1
              k1 = k1 + ny
              do iy1 = iyfrom, iyto
                l1 = (iy1 + k1)*nx
                ixfrom = Max (ix-1, 0)
                ixto   = Min (ix+2, nx) - 1
                l1 = l1 + ixfrom - 1
                do ix1 = ixfrom, ixto
                  l1 = l1 + 1

                  do jj = atoms_in_box(l1), atoms_in_box(l1+1) - 1
                    j = s_box(jj)

                    if (j /= i) then
                      dist = (xa(1) - coord(1, j))**2 + &
                           & (xa(2) - coord(2, j))**2 + &
                           & (xa(3) - coord(3, j))**2
                      if (dist < (rr+srad(j))**2) then
                        ! ---
                        ! For some stupid reasons, the code
                        ! below selecting nn atoms depends
                        ! on the order of atoms (and the
                        ! solvation energy of a molecule
                        ! depends on the order of atoms!)
                        ! so to make the code identical to
                        ! original Mopac Cosmo, instead of
                        ! the simple
                        !   nipa(i) = nipa(i) + 1
                        !   ilipa = ilipa + 1
                        !   lipa(ilipa) = j
                        ! we use the following:

                        i0 = ilipa + 1
                        do ik = ilipa, ilipa-nipa(i)+1, -1
                          if (lipa(ik) > j) then
                            lipa(i0) = lipa(ik)
                            i0 = i0 - 1
                          else
                            exit
                          end if
                        end do

                        nipa(i) = nipa(i) + 1
                        ilipa = ilipa + 1
                        lipa(i0) = j
                      end if
                    end if
                  end do
!                  if (ilipa >= 70*numat) then
!                    exit
!                  end if
                end do   ! ix1
!                if (ilipa >= 70*numat) then
!                  exit
!                end if
              end do   ! iy1
!              if (ilipa >= 70*numat) then
!                exit
!              end if
            end do   ! iz1
!            if (ilipa >= 70*numat) then
!              exit
!            end if
          end do
!          if (ilipa >= 70*numat) then
!            exit
!          end if
        end do   ! ix
!        if (ilipa >= 70*numat) then
!          exit
!        end if
      end do   ! iy
!      if (ilipa >= 70*numat) then
!        exit
!      end if
    end do   ! iz
    !
    !       At this stage, array lipa need to be reodered, so
    !       it correspond atoms in order 1...n, instead of
    !       s_box(0) ... s_box(n-1)
    !

    box_number(1) = 0
    do i = 2, numat
      box_number(i) = box_number(i-1) + nipa(i-1)     ! list starts
    end do

    indi = 0
    do ii = 0, numat - 1
      i = s_box(ii)
      do jj = 1, nipa(i)
        temp(box_number(i) + jj) = lipa(indi + jj)
      end do

      indi = indi + nipa(i)
    end do

    lipa(1:indi) = temp(1:indi)

    !
    !       Main loop
    !
    inset = 1
    ilipa = 0
    nps = 0
    area = 0.d0
    cosvol = 0.d0

    do i = 1, numat
      nps0 = nps + 1
      ri = srad(i)
      r = ri + rsolv
      rr = r + rsolv
      ri2 = ri * ri
      xa(1:3) = coord(1:3, i)

      !
      !       SEARCH FOR 3 NEAREST NEIGHBOR ATOMS
      !

      dist1 = 1.d20
      dist2 = 1.d20
      dist3 = 1.d20
      nn(1, i) = 0
      nn(2, i) = 0
      nn(3, i) = 0

      do jj = 1, nipa(i)
        j = lipa(ilipa + jj)
        dist = (xa(1) - coord(1, j))**2 + &
             & (xa(2) - coord(2, j))**2 + &
             & (xa(3) - coord(3, j))**2

        if (dist+0.05d0 < dist3) then
          dist3 = dist
          nn(3, i) = j
        end if

        if (dist3+0.05d0 < dist2) then
          dist = dist2
          dist2 = dist3
          dist3 = dist
          nn(3, i) = nn(2, i)
          nn(2, i) = j
        end if

        if (dist2+0.05d0 < dist1) then
          dist = dist1
          dist1 = dist2
          dist2 = dist
          nn(2, i) = nn(1, i)
          nn(1, i) = j
        end if
      end do

      ilipa = ilipa + nipa(i)

      !
      !       BUILD NEW TRANSFORMATION MATRIX
      !

      if (nn(1, i) == 0) then
        tm(1, 1, i) = 1.d0
        tm(1, 2, i) = 0.d0
        tm(1, 3, i) = 0.d0
      else
        dist1 = (xa(1) - coord(1, nn(1, i)))**2 + &
             & (xa(2) - coord(2, nn(1, i)))**2 + &
             & (xa(3) - coord(3, nn(1, i)))**2
        dist = 1.d0 / Sqrt (dist1)

        tm(1, 1, i) = (coord(1, nn(1, i))-xa(1)) * dist
        tm(1, 2, i) = (coord(2, nn(1, i))-xa(2)) * dist
        tm(1, 3, i) = (coord(3, nn(1, i))-xa(3)) * dist
      end if

      do
        if (nn(2, i) == 0) then
          tm(2, 1, i) = -tm(1, 2, i)
          tm(2, 2, i) = tm(1, 1, i)
          tm(2, 3, i) = 0.d0
          exit
        else
          dist2 = (xa(1) - coord(1, nn(2, i)))**2 + &
               & (xa(2) - coord(2, nn(2, i)))**2 + &
               & (xa(3) - coord(3, nn(2, i)))**2
          dist = 1 / Sqrt (dist2)

          xx(1) = (coord(1, nn(2, i))-xa(1)) * dist
          xx(2) = (coord(2, nn(2, i))-xa(2)) * dist
          xx(3) = (coord(3, nn(2, i))-xa(3)) * dist
          sp = xx(1)*tm(1, 1, i) + xx(2)*tm(1, 2, i) + xx(3)*tm(1, 3, i)

          if (sp*sp > 0.99d0) then
            nn(2, i) = nn(3, i)
            nn(3, i) = 0
            dist2 = dist3
          else
            sininv = 1 / Sqrt (1.d0-sp*sp)
            tm(2, 1, i) = (xx(1)-sp*tm(1, 1, i)) * sininv
            tm(2, 2, i) = (xx(2)-sp*tm(1, 2, i)) * sininv
            tm(2, 3, i) = (xx(3)-sp*tm(1, 3, i)) * sininv
            exit
          end if
        end if
      end do

      tm(3, 1, i) = tm(1, 2, i)*tm(2, 3, i) - tm(2, 2, i)*tm(1, 3, i)
      tm(3, 2, i) = tm(1, 3, i)*tm(2, 1, i) - tm(2, 3, i)*tm(1, 1, i)
      tm(3, 3, i) = tm(1, 1, i)*tm(2, 2, i) - tm(2, 1, i)*tm(1, 2, i)

      !
      !       TRANSFORM DIRVEC ACCORDING TO TM
      !
      do j = 1, 1082
        xx(1) = dirvec(1, j)
        xx(2) = dirvec(2, j)
        xx(3) = dirvec(3, j)
        dirtm(1, j) = xx(1)*tm(1, 1, i) + xx(2)*tm(2, 1, i) + &
             & xx(3)*tm(3, 1, i)
        dirtm(2, j) = xx(1)*tm(1, 2, i) + xx(2)*tm(2, 2, i) + &
             & xx(3)*tm(3, 2, i)
        dirtm(3, j) = xx(1)*tm(1, 3, i) + xx(2)*tm(2, 3, i) + &
             & xx(3)*tm(3, 3, i)
      end do

      !
      !       FIND THE POINTS OF THE BASIC GRID ON THE SAS
      !

      narea = 0
      loop: do j = 1, 1082
        din(j) = .false.
        xx(1) = xa(1) + dirtm(1, j) *r
        xx(2) = xa(2) + dirtm(2, j) *r
        xx(3) = xa(3) + dirtm(3, j) *r
        !
        !             --- WE NEED ONLY TRY THOSE ATOMS INTERSECTING ATOM I
        !
        do ik = ilipa - nipa(i) + 1, ilipa
          k = lipa(ik)
          dist = (xx(1) - coord(1, k))**2 + &
               & (xx(2) - coord(2, k))**2 + &
               & (xx(3) - coord(3, k))**2
          dist = Sqrt (dist) - rsolv - srad(k)
!
! jjps 8-22-2022:  An intermittent error caused the solvent-accessible surface to
! change size significantly, e.g., by 6 square Angstroms in a 1300 atom protein,
! on making a small change, 0.1 Angstroms total, in geometry.
! This error was traced to a special case when "dist" was almost zero.  The error
! was corrected by replacing the original test "if (dist < 0) cycle loop" with
! the text on the next line.  Surprisingly, this appears to have fixed the bug!
!
          if (dist < -0.001d0) cycle loop
        end do

        narea = narea + 1
        cosvol = cosvol + ri2 * dirvec(4, j) * &
             & (dirtm(1, j)*xa(1)+dirtm(2, j)*xa(2) + dirtm(3, j)*xa(3)+ri)
        area = area + ri2 * dirvec(4, j)
        din(j) = .true.
      end do loop

      if (narea /= 0) then
        !
        !  IF HYDROGEN, THEN USE THE SMALLER SET OF POINTS (NORMALLY 12 POINTS)
        !  OTHERWISE, USE THE LARGER SET (NORMALLY 42 POINTS)
        !
        i0 = 1
        if (nat(i) == 1) then
          i0 = 2
        end if

        jmax = n0(i0)
        i0 = (i0-1) * n0(1)
        do j = 1, jmax
          nps = nps + 1
          if (nps > lenabc) then
            call mopend ("NPS IS GREATER THAN &
                 &LENABC-USE SMALLER NSPA")
            return
          else
            iatsp(nps) = i
            xx(1) = dirsm(1, i0+j)
            xx(2) = dirsm(2, i0+j)
            xx(3) = dirsm(3, i0+j)
            cosurf(1, nps) = xx(1)*tm(1, 1, i) + &
                 & xx(2)*tm(2, 1, i) + &
                 & xx(3)*tm(3, 1, i)
            cosurf(2, nps) = xx(1)*tm(1, 2, i) + &
                 & xx(2)*tm(2, 2, i) + &
                 & xx(3)*tm(3, 2, i)
            cosurf(3, nps) = xx(1)*tm(1, 3, i) + &
                 & xx(2)*tm(2, 3, i) + &
                 & xx(3)*tm(3, 3, i)
          end if
        end do

        niter = 0
        do
          niter = niter + 1
          do ips = nps0, nps
            nar_loc(ips) = 0
            phinet(ips, 1) = 0.d0
            phinet(ips, 2) = 0.d0
            phinet(ips, 3) = 0.d0
          end do

          do j = 1, 1082
            if (din(j)) then
              ipm = nps0
              spm = -1.d0
              x1 = dirtm(1, j)
              x2 = dirtm(2, j)
              x3 = dirtm(3, j)
              do ips = nps0, nps
                sp = x1*cosurf(1, ips)+ x2*cosurf(2, ips)+ x3*cosurf(3, ips)
                if (sp >= spm) then
                  spm = sp*(1.d0+1.d-14)
                  ipm = ips
                end if
              end do

              iseg(j) = ipm
              nar_loc(ipm) = nar_loc(ipm) + 1
              phinet(ipm, 1)=phinet(ipm, 1) + dirtm(1, j) * dirvec(4, j)
              phinet(ipm, 2)=phinet(ipm, 2) + dirtm(2, j) * dirvec(4, j)
              phinet(ipm, 3)=phinet(ipm, 3) + dirtm(3, j) * dirvec(4, j)
            end if
          end do

          ips = nps0 - 1
          loop1: do
            ips = ips + 1
            do while (nar_loc(ips) == 0)
              niter = 1
              nps = nps - 1
              if (ips > nps) exit loop1
              do jps = ips, nps
                nar_loc(jps) = nar_loc(jps+1)
                phinet(jps, 1) = phinet(jps+1, 1)
                phinet(jps, 2) = phinet(jps+1, 2)
                phinet(jps, 3) = phinet(jps+1, 3)
              end do
            end do

            dists = phinet(ips, 1) ** 2 + &
                 & phinet(ips, 2) ** 2 + &
                 & phinet(ips, 3) ** 2

            dists = Max (dists, 1.d-20)
            dist = 1.d0 / Sqrt (dists)
            cosurf(1, ips) = phinet(ips, 1) * dist
            cosurf(2, ips) = phinet(ips, 2) * dist
            cosurf(3, ips) = phinet(ips, 3) * dist
            if (ips >= nps) exit
          end do loop1

          if (niter >= 2) exit
        end do
        !
        !        NOW ALL SEGMENTS ARE FINALLY DEFINED AND THE ASSOCIATED
        !        BASIC GRID POINTS ARE CLOSE-PACKED
        !
        do ips = nps0, nps
          nsetf(ips) = inset
          inset = inset + nar_loc(ips)
          nar_loc(ips) = 0
          cosurf(4, ips) = 0.d0
          cosurf(1, ips) = cosurf(1, ips) * ri + xa(1)
          cOsurf(2, ips) = cosurf(2, ips) * ri + xa(2)
          cosurf(3, ips) = cosurf(3, ips) * ri + xa(3)
        end do

        do j = 1, 1082
          if (din(j)) then
            ipm = iseg(j)
            nara = nar_loc(ipm)
            nset(nsetf(ipm)+nara) = j
            nar_loc(ipm) = nara + 1
            cosurf(4, ipm) = cosurf(4, ipm) + dirvec(4, j) * ri2
          end if
        end do
      end if
    end do

    !
    !       HERE THE CONSTRUCTION FOR A SINGLE ATOM ENDS
    !

    do i = 1, numat
      din(i) = .true.
    end do

    do j = 1, nps
      din(iatsp(j)) = .false.
    End do
    s_box(0:numat-1) = 0
    do i=1, nps
      j = iatsp(i) -1
      s_box(j) = s_box(j) + Int(cosurf(4, i))
    end do

    do i=1, numat
      s_box(i-1) = Int(s_box(i-1) / (srad(i))**2 * (srad(i) + rsolv)**2)
    end do

 !   write(iw, *) ' Cosmo atomic radii and accessible surface area '

 !   do i=1, numat
  !    write(iw, '(I10, F12.4, I10)') i, srad(i), s_box(i-1)
  !  end do
    !
    !     NOW THE SEGMENT FORMATION ON ALL ATOMS IS FINISHED
    !     NOW THE CLOSURE OF THE CONCAVE REGIONS OF THE SURFACE WILL BE DONE
    !

    if (ioldcv == 0) then
      i = max(numat, 1082)
      call surclo (coord, nipa, lipa, din, i, rsc, isort, &
           &      ipsrs, nipsrs, nat, srad, cosurf, iatsp, &
           &      nar_loc, nsetf, isude, sude, maxrs)
    end if

    cosvol = cosvol / 3

    npoints(1:numat) = 0
    do i=1, nps
      j = iatsp(i)
      npoints(j) = npoints(j) + 1
    end do
    !
    !       At this point, the npoints array has number of elements
    !       for each atom. Transform it to the list ...
    !
    j = npoints(1)
    npoints(1) = 1
    do i = 2, numat
      k = npoints(i)
      npoints(i) = npoints(i-1) + j
      j = k
    end do
    npoints(numat+1) = npoints(numat) + j

    arat(1:numat) = 0.d0
    do i = 1, nps
      j = iatsp(i)
      arat(j) = arat(j) + cosurf(4, i)
    end do


    deallocate (atoms_in_box, s_box, box_number, din, &
         & nn, isort, ipsrs, nipa, lipa, temp, stat=alloc_stat)

    if (alloc_stat /= 0) then
      write (iw, *) " Deallocation error in CoscanZ"
    end if

    if (allocated (r_vec)) then
      deallocate (r_vec, p_vec, q_vec, z_vec, &
           & a_diag, stat=alloc_stat)
      if  (alloc_stat /= 0) then
        write (iw, *) " Deallocation error in CoscanZ (2)"
        goto 100
      end if
    end if

    allocate (r_vec(nps), p_vec(nps), q_vec(Max (numat, nps)), &
         & z_vec(nps), a_diag(nps), stat=alloc_stat)

    if (alloc_stat /= 0) then
      write (iw, *) " Allocation error in CoscanZ"
    end if

100 continue

  end subroutine coscanz

  subroutine aq_vec (coord, srad, numat, cosurf, nps, nar_csm, nsetf, &
       & nset, rsc, nipsrs, iatsp, tm, ioldcv, maxrs, lenabc, a_diag, q, v)
    !
    ! Computes A*q product using N**2 algorithm
    !
    use cosmo_C, only : dirvec, disex2
    !
    !.. Implicit Declarations ..
    implicit none
    !
    !.. Formal Arguments ..
    integer, intent(in) :: numat, nps, ioldcv, maxrs, lenabc

    double precision, dimension(3, numat), intent(in) :: coord
    double precision, dimension(numat), intent(in) :: srad
    double precision, dimension(4, nps), intent(in) :: cosurf
    integer, dimension(lenabc+1), intent(in) :: nar_csm, nsetf, iatsp
    integer, dimension(*), intent(in) :: nset
    double precision, dimension(4, maxrs), intent(in) :: rsc
    integer, dimension(lenabc), intent(in) :: nipsrs
    double precision, dimension(3, 3, numat), intent(in) :: tm
    double precision, dimension(nps), intent(in) :: a_diag
    double precision, dimension(nps), intent(in) :: q

    double precision, dimension(nps), intent(out) :: v

    !
    ! .. Parameters
    !
    !.. Local Scalars ..
    integer :: i, j, ii, jj, k, l, nfl1, nfl2, ijpos
    double precision :: aa, d2, ri, rj, x1, x2, x3, x4
    !
    !.. Local Arrays ..
    double precision, dimension(4, 500, 2) :: finel
    double precision xi(3), xj(3), xa(3), xb(3)
    v = 0
    ijpos = 0

    do ii = 1, nps
      i = iatsp(ii)
      xa(1:3) = cosurf(1:3, ii)

      if (.not. compute_a_part) then
        ri = srad(i)
        xi(1:3) = coord(1:3, i)
        call mfinel (ii, 1, finel, nar_csm, nsetf, nset, rsc, &
           & nipsrs, dirvec, tm(1, 1, i), xi, ri, &
           & nfl1, ioldcv, maxrs, lenabc, numat)
      end if

      do jj = 1, ii-1
        j = iatsp(jj)
        xb(1:3) = cosurf(1:3, jj)

        d2 = (xb(1)-xa(1))**2 + (xb(2)-xa(2))**2 + &
             & (xb(3)-xa(3))**2

        if (d2 > disex2) then
          aa = 1.d0 / Sqrt (d2)
        else
          if (compute_a_part) then
            ijpos = ijpos +1
            aa = a_part(ijpos)
          else
            xj(1:3) = coord(1:3, j)
            rj = srad(j)
            call mfinel (jj, 2, finel, nar_csm, nsetf, &
               & nset, rsc, nipsrs, dirvec, &
               & tm(1, 1, j), xj, rj, nfl2, &
               & ioldcv, maxrs, lenabc, numat)

            aa = 0.d0
            do k = 1, nfl1
              x1 = finel(1, k, 1)
              x2 = finel(2, k, 1)
              x3 = finel(3, k, 1)
              x4 = finel(4, k, 1)

              do l = 1, nfl2
                aa = aa + x4 * finel(4, l, 2) / &
                   & Sqrt ((x1-finel(1, l, 2))**2 + &
                   &       (x2-finel(2, l, 2))**2 + &
                   &       (x3-finel(3, l, 2))**2)
              end do
            end do

            aa = aa / (cosurf(4, ii) * cosurf(4, jj))
          end if
        end if
        v(ii) = v(ii) + aa * q(jj)
        v(jj) = v(jj) + aa * q(ii)
      end do
    end do
    v = v + a_diag * q
  end subroutine aq_vec

  subroutine simulate_aq_vec (coord, srad, numat, cosurf, nps, nar_csm, nsetf, &
       & nset, rsc, nipsrs, iatsp, tm, ioldcv, maxrs, lenabc, compute, npos)
    !
    ! Computes number of short-range elements in A matrix and
    ! computes the short-range A matrix elements if requested
    ! using N**2 algorithm (off-diagonal elements only)
    !
    use cosmo_C, only : dirvec, disex2

    !
    !.. Implicit Declarations ..
    implicit none
    !
    !.. Formal Arguments ..
    integer, intent(in) :: numat, nps, ioldcv, maxrs, lenabc

    double precision, dimension(3, numat), intent(in) :: coord
    double precision, dimension(numat), intent(in) :: srad
    double precision, dimension(4, nps), intent(in) :: cosurf
    integer, dimension(lenabc+1), intent(in) :: nar_csm, nsetf, iatsp
    integer, dimension(*), intent(in) :: nset
    double precision, dimension(4, maxrs), intent(in) :: rsc
    integer, dimension(lenabc), intent(in) :: nipsrs
    double precision, dimension(3, 3, numat), intent(in) :: tm
    logical, intent(in) :: compute
    integer, intent(out) :: npos
    integer :: i, j, ii, jj, k, l, nfl1, nfl2
    double precision :: aa, d2, ri, rj, x1, x2, x3, x4
    double precision, dimension(4, 500, 2) :: finel
    double precision xi(3), xj(3), xa(3), xb(3)
    npos = 0

    do ii = 1, nps
      i = iatsp(ii)
      ri = srad(i)
      xi(1:3) = coord(1:3, i)
      xa(1:3) = cosurf(1:3, ii)

      if(compute) then
        call mfinel (ii, 1, finel, nar_csm, nsetf, nset, rsc, &
           & nipsrs, dirvec, tm(1, 1, i), xi, ri, &
           & nfl1, ioldcv, maxrs, lenabc, numat)
      end if

      do jj = 1, ii -1
        j = iatsp(jj)
        xj(1:3) = coord(1:3, j)
        xb(1:3) = cosurf(1:3, jj)

        d2 = (xb(1)-xa(1))**2 + (xb(2)-xa(2))**2 + &
             & (xb(3)-xa(3))**2

        if (d2 <= disex2) then
          npos = npos +1
          if(compute) then
            j = iatsp(jj)
            rj = srad(j)
            call mfinel (jj, 2, finel, nar_csm, nsetf, &
               & nset, rsc, nipsrs, dirvec, &
               & tm(1, 1, j), xj, rj, nfl2, &
               & ioldcv, maxrs, lenabc, numat)

            aa = 0.d0
            do k = 1, nfl1
              x1 = finel(1, k, 1)
              x2 = finel(2, k, 1)
              x3 = finel(3, k, 1)
              x4 = finel(4, k, 1)

              do l = 1, nfl2
                aa = aa + x4 * finel(4, l, 2) / &
                   & Sqrt ((x1-finel(1, l, 2))**2 + &
                   &       (x2-finel(2, l, 2))**2 + &
                   &       (x3-finel(3, l, 2))**2)
              end do
            end do

            aa = aa / (cosurf(4, ii) * cosurf(4, jj))

            a_part(npos) = aa
          end if
        end if
      end do
    end do
  end subroutine simulate_aq_vec

  subroutine aq_dir_int (ind1, n1, ind2, n2, c, nc, q, r, same)
    !
    ! Computes part of the A*q product for near surface elements
    !

    use common_arrays_C, only: coord
    use molkst_C, only : numat
    use cosmo_C, only: dirvec, lenabc, ioldcv, disex2, cosurf, srad, &
      nar_csm, iatsp, nsetf
    implicit none
    integer, intent(in) :: n1, n2, nc
    integer, intent(in), dimension(*) :: ind1, ind2
    double precision, intent(in), dimension(nc, *) :: c
    double precision, intent(in), dimension(*) :: q
    double precision, intent(inout), dimension(*) :: r
    logical, intent(in) :: same
    integer :: i, j, ii, jj, i3, j3, k, l, nfl1, nfl2, maxrs
    integer, save :: ijpos
    double precision :: aa, d2, ri, rj, x1, x2, x3, x4
    double precision, dimension(4, 500, 2) :: finel
    double precision xi(3), xj(3), xa(3), xb(3)
!
! Dummy statement to "use" c
!
    ri = c(1,1)
    if (new_iteration) then
      ijpos = 0
      new_iteration = .false.
    end if

    maxrs = 60*numat

    if (same) then
      do i3=1, n1
        ii = ind1(i3)
        i = iatsp(ii)
        ri = srad(i)
        xi(1:3) = coord(1:3, i)
        xa(1:3) = cosurf(1:3, ii)

        if (.not. compute_a_part) then
          call mfinel (ii, 1, finel, nar_csm, nsetf, nset, rsc, &
             & nipsrs, dirvec, tm(1, 1, i), xi, ri, &
             & nfl1, ioldcv, maxrs, lenabc, numat)
        end if

        do j3 = 1, i3 - 1
          jj = ind1(j3)
          j = iatsp(jj)
          xj(1:3) = coord(1:3, j)
          xb(1:3) = cosurf(1:3, jj)

          d2= (xb(1)-xa(1))**2 + (xb(2)-xa(2))**2 + (xb(3)-xa(3))**2

          if (d2 > disex2) then
            aa = 1.d0 / Sqrt (d2)
          else
            if (compute_a_part) then
              ijpos = ijpos +1
              aa = a_part(ijpos)
            else
              rj = srad(j)
              call mfinel (jj, 2, finel, nar_csm, nsetf, &
                 & nset, rsc, nipsrs, dirvec, &
                 & tm(1, 1, j), xj, rj, nfl2, &
                 & ioldcv, maxrs, lenabc, numat)

              aa = 0.d0
              do k = 1, nfl1
                x1 = finel(1, k, 1)
                x2 = finel(2, k, 1)
                x3 = finel(3, k, 1)
                x4 = finel(4, k, 1)

                do l = 1, nfl2
                  aa = aa + x4 * finel(4, l, 2) / &
                     & Sqrt ((x1-finel(1, l, 2))**2+ &
                     & (x2-finel(2, l, 2))**2+ &
                     & (x3-finel(3, l, 2))**2)
                end do
              end do

              aa = aa / (cosurf(4, ii) * cosurf(4, jj))
            end if
          end if

          r(ii) = r(ii) + aa * q(jj)
          r(jj) = r(jj) + aa * q(ii)
        end do
      end do
    else
      do i3 = 1, n1
        ii = ind1(i3)
        i = iatsp(ii)
        ri = srad(i)
        xi(1:3) = coord(1:3, i)
        xa(1:3) = cosurf(1:3, ii)

        if (.not. compute_a_part) then
          call mfinel (ii, 1, finel, nar_csm, nsetf, nset, rsc, &
             & nipsrs, dirvec, tm(1, 1, i), xi, ri, &
             & nfl1, ioldcv, maxrs, lenabc, numat)
        end if

        do j3 = 1, n2
          jj = ind2(j3)
          j = iatsp(jj)
          xj(1:3) = coord(1:3, j)
          xb(1:3) = cosurf(1:3, jj)

          d2 = (xb(1)-xa(1))**2 + (xb(2)-xa(2))**2 + (xb(3)-xa(3))**2

          if (d2 > disex2) then
            aa = 1.d0 / Sqrt (d2)
          else
            if (compute_a_part) then
              ijpos = ijpos +1
              aa = a_part(ijpos)
            else
              rj = srad(j)
              call mfinel (jj, 2, finel, nar_csm, nsetf, &
                 & nset, rsc, nipsrs, dirvec, &
                 & tm(1, 1, j), xj, rj, nfl2, &
                 & ioldcv, maxrs, lenabc, numat)

              aa = 0.d0
              do k = 1, nfl1
                x1 = finel(1, k, 1)
                x2 = finel(2, k, 1)
                x3 = finel(3, k, 1)
                x4 = finel(4, k, 1)

                do l = 1, nfl2
                  aa = aa + x4 * finel(4, l, 2) / &
                     & Sqrt ((x1-finel(1, l, 2))**2+ &
                     & (x2-finel(2, l, 2))**2+ &
                     & (x3-finel(3, l, 2))**2)
                end do
              end do

              aa = aa / (cosurf(4, ii) * cosurf(4, jj))
            end if
          end if

          r(ii) = r(ii) + aa * q(jj)
          r(jj) = r(jj) + aa * q(ii)
        end do
      end do
    end if
  end subroutine aq_dir_int

  subroutine simulate_aq_dir_int (ind1, n1, ind2, n2, c, nc, same, compute, &
                & npos)

    !
    ! Computes number of short-range elements in A matrix and
    ! computes the short-range A matrix elements if requested
    !

    use common_arrays_C, only: coord
    use molkst_C, only : numat
    use cosmo_C, only: dirvec, lenabc, ioldcv, disex2, cosurf, srad, &
      iatsp, nar_csm, nsetf
    implicit none
    integer, intent(in) :: n1, n2, nc
    integer, intent(in), dimension(*) :: ind1, ind2
    double precision, intent(in), dimension(nc, *) :: c
    logical, intent(in) :: same, compute
    integer, intent (inout) :: npos
    integer :: i, j, ii, jj, i3, j3, k, l, nfl1, nfl2, maxrs
    double precision :: aa, d2, ri, rj, x1, x2, x3, x4
    double precision, dimension(4, 500, 2) :: finel
    double precision xi(3), xj(3), xa(3), xb(3)
!
! Dummy statement to "use" c
!
    ri = c(1,1)
    maxrs = 60*numat

    if (same) then
      do i3=1, n1
        ii = ind1(i3)
        i = iatsp(ii)
        ri = srad(i)
        xi(1:3) = coord(1:3, i)
        xa(1:3) = cosurf(1:3, ii)

        if (compute) then
          call mfinel (ii, 1, finel, nar_csm, nsetf, nset, rsc, &
             & nipsrs, dirvec, tm(1, 1, i), xi, ri, &
             & nfl1, ioldcv, maxrs, lenabc, numat)
        end if

        do j3 = 1, i3 -1
          jj = ind1(j3)
          j = iatsp(jj)
          xj(1:3) = coord(1:3, j)
          xb(1:3) = cosurf(1:3, jj)

          d2= (xb(1)-xa(1))**2 + (xb(2)-xa(2))**2 + (xb(3)-xa(3))**2

          if (d2 <= disex2) then
            npos = npos +1
            if (compute) then
              rj = srad(j)
              call mfinel (jj, 2, finel, nar_csm, nsetf, &
                 & nset, rsc, nipsrs, dirvec, &
                 & tm(1, 1, j), xj, rj, nfl2, &
                 & ioldcv, maxrs, lenabc, numat)

              aa = 0.d0
              do k = 1, nfl1
                x1 = finel(1, k, 1)
                x2 = finel(2, k, 1)
                x3 = finel(3, k, 1)
                x4 = finel(4, k, 1)

                do l = 1, nfl2
                  aa = aa + x4 * finel(4, l, 2) / &
                       & Sqrt ((x1-finel(1, l, 2))**2+ &
                       & (x2-finel(2, l, 2))**2+ &
                       & (x3-finel(3, l, 2))**2)
                end do
              end do

              aa = aa / (cosurf(4, ii) * cosurf(4, jj))

              a_part(npos) = aa
            end if
          end if
        end do
      end do
    else
      do i3 = 1, n1
        ii = ind1(i3)
        i = iatsp(ii)
        ri = srad(i)
        xi(1:3) = coord(1:3, i)
        xa(1:3) = cosurf(1:3, ii)
        if (compute) then
          call mfinel (ii, 1, finel, nar_csm, nsetf, nset, rsc, &
             & nipsrs, dirvec, tm(1, 1, i), xi, ri, &
             & nfl1, ioldcv, maxrs, lenabc, numat)
        end if
        do j3 = 1, n2
          jj = ind2(j3)
          j = iatsp(jj)
          xj(1:3) = coord(1:3, j)
          xb(1:3) = cosurf(1:3, jj)

          d2 = (xb(1)-xa(1))**2 + (xb(2)-xa(2))**2 + (xb(3)-xa(3))**2

          if (d2 <= disex2) then
            npos = npos +1

            if (compute) then
              rj = srad(j)
              call mfinel (jj, 2, finel, nar_csm, nsetf, &
                 & nset, rsc, nipsrs, dirvec, &
                 & tm(1, 1, j), xj, rj, nfl2, &
                 & ioldcv, maxrs, lenabc, numat)

              aa = 0.d0
              do k = 1, nfl1
                x1 = finel(1, k, 1)
                x2 = finel(2, k, 1)
                x3 = finel(3, k, 1)
                x4 = finel(4, k, 1)

                do l = 1, nfl2
                  aa = aa + x4 * finel(4, l, 2) / &
                     & Sqrt ((x1-finel(1, l, 2))**2+ &
                     & (x2-finel(2, l, 2))**2+ &
                     & (x3-finel(3, l, 2))**2)
                end do
              end do

              aa = aa / (cosurf(4, ii) * cosurf(4, jj))

              a_part(npos) = aa
            end if
          end if
        end do
      end do
    end if
  end subroutine simulate_aq_dir_int

  subroutine precond_aq_dir_int (ind1, n1, ind2, n2, c, nc, same, compute, &
                & npos)
    !
    ! Computes block preconditioner
    !
    use cosmo_C, only : iatsp, cosurf, disex2
    implicit none
    integer, intent(in) :: n1, n2, nc
    integer, intent(in), dimension(*) :: ind1, ind2
    double precision, intent(in), dimension(nc, *) :: c
    logical, intent(in) :: same, compute
    integer, intent (inout) :: npos
    integer :: i, j, ii, jj, i3, j3, ips, jps
    double precision :: aa, d2
    double precision xa(3), xb(3)
!
! dummy statement to "use" compute and c
!
    if (compute) aa = c(1,1)
    m_vec = 0.d0
    if (same) then
      do i3=1, n1
        ii = ind1(i3)
        i = iatsp(ii)
        xa(1:3) = cosurf(1:3, ii)

        do j3 = 1, i3 -1
          jj = ind1(j3)
          j = iatsp(jj)
          xb(1:3) = cosurf(1:3, jj)

          d2 = (xb(1)-xa(1))**2 + (xb(2)-xa(2))**2 + &
              & (xb(3)-xa(3))**2

          if (d2 > disex2) then
            aa = 1.d0 / Sqrt (d2)
          else
            npos = npos +1
            aa = a_part(npos)
          end if
          if (i == j) then ! same atom
            ips = ii - npoints(i) +1
            jps = jj - npoints(i) +1

            if (ips >= jps) then
              m_vec(iblock_pos(i) -1 + (ips*(ips-1))/2 + jps) = aa
            else
              m_vec(iblock_pos(i) -1 + (jps*(jps-1))/2 + ips) = aa
            end if

          end if

        end do
      end do

    else

      do i3 = 1, n1
        ii = ind1(i3)
        i = iatsp(ii)
        xa(1:3) = cosurf(1:3, ii)

        do j3 = 1, n2
          jj = ind2(j3)
          j = iatsp(jj)
          xb(1:3) = cosurf(1:3, jj)

          d2 = (xb(1)-xa(1))**2 + (xb(2)-xa(2))**2 + &
              & (xb(3)-xa(3))**2

          if (d2 > disex2) then
            aa = 1.d0 / Sqrt (d2)
          else
            npos = npos +1
            aa = a_part(npos)
          end if

          if (i == j) then ! same atom

            ips = ii - npoints(i) +1
            jps = jj - npoints(i) +1
            if (ips >= jps) then
              m_vec(iblock_pos(i) -1 + (ips*(ips-1))/2 + jps) = aa
            else
              m_vec(iblock_pos(i) -1 + (jps*(jps-1))/2 + ips) = aa
            end if

          end if

        end do
      end do
    end if

  end subroutine precond_aq_dir_int

  subroutine aq_mult_int (ind, n, psi, p, x0, y0, z0, y_norm, pmn, c, n3, &
       & q, r)
    !
    ! Computes part of the A*q product for far surface elements,
    ! using local expansions
    !
    implicit none
    integer, intent(in) :: n, p, n3
    integer, intent(in), dimension(*) :: ind
    complex (kind=dcmplx_kind), dimension(-p:p, 0:p), intent(in) :: psi
    double precision, dimension(-p:p, 0:p), intent(in) :: y_norm
    double precision, dimension(-p:p, 0:p), intent(out) :: pmn
    double precision, intent(in), dimension(n3, *) :: c
    double precision, intent(in), dimension(*) :: q
    double precision, intent(inout), dimension(*) :: r
    double precision, intent(in) :: x0, y0, z0
    integer :: i, j, k
    double precision :: dx, dy, dz, d, phi, cos_t, t
    complex (kind=dcmplx_kind) :: tc
    intrinsic Real, Cmplx, Atan2, Sqrt, Exp
!
! Dummy statement to "use" q
!
    dx = q(1)
    do i = 1, n
      dx = c(1, ind(i)) - x0
      dy = c(2, ind(i)) - y0
      dz = c(3, ind(i)) - z0

      d = Sqrt (dx**2 + dy**2 + dz**2)
      cos_t = dz / d
      phi = Atan2 (dy, dx)

      call get_legendre (p, cos_t, pmn)

      t = 0.d0
      do j = 0, p
        tc = (0.d0, 0.d0)

        do k = 1, j
          tc = tc + y_norm(k, j) * pmn(k, j) * d**j &
               & * psi(k, j) * Exp (Cmplx (0.d0, k*phi, dcmplx_kind))
        end do

        t = t + 2.d0 * Real (tc) + pmn(0, j) * d**j * Real (psi(0, j))
      end do

      r(ind(i)) = r(ind(i)) + t

    end do

  end subroutine aq_mult_int

  subroutine aq_far_int (ind, n, psi, p, x0, y0, z0, y_norm, pmn, c, nc, &
       & q, r)
    !
    ! Computes part of the A*q product for far surface elements,
    ! using multipole expansions
    !
    implicit none
    integer, intent(in) :: n, p, nc
    integer, intent(in), dimension(*) :: ind
    complex (kind=dcmplx_kind), dimension(-p:p, 0:p), intent(in) :: psi
    double precision, dimension(-p:p, 0:p), intent(in) :: y_norm
    double precision, dimension(-p:p, 0:p), intent(out) :: pmn
    double precision, intent(in), dimension(nc, *) :: c
    double precision, intent(in), dimension(*) :: q
    double precision, intent(inout), dimension(*) :: r
    double precision, intent(in) :: x0, y0, z0
    integer :: i, j, k
    double precision :: dx, dy, dz, d, phi, cos_t, t
    complex (kind=dcmplx_kind) :: tc
    intrinsic Real, Cmplx, Atan2, Exp, Sqrt
!
! Dummy statement to "use" q
!
    dx = q(1)
    do i = 1, n
      dx =  c(1, ind(i)) - x0
      dy =  c(2, ind(i)) - y0
      dz =  c(3, ind(i)) - z0

      d = Sqrt (dx**2 + dy**2 + dz**2)
      cos_t = dz / d
      phi = Atan2 (dy, dx)

      call get_legendre (p, cos_t, pmn)

      t = 0
      do j = 0, p
        tc = (0.d0, 0.d0)
        do k =  1, j
          tc = tc + y_norm(k, j) * pmn(k, j) / d**(j+1) &
               *psi(k, j) * Exp (Cmplx (0, k*phi, dcmplx_kind))

        end do
        t = t + 2.d0 * Real (tc) + pmn(0, j) / d**(j+1) * Real (psi(0, j))
      end do

      r(ind(i)) = r(ind(i)) + t
    end do

  end subroutine aq_far_int

  subroutine bz_far_int (ind, n, psi, p, x0, y0, z0, y_norm, pmn, c, nc, q, r)
    !
    ! Computes part of the B*q product for far surface elements,
    ! using multipole expansions
    !


    use cosmo_C, only: cosurf
    implicit none
    integer, intent(in) :: n, p, nc
    integer, intent(in), dimension(*) :: ind
    complex (kind=dcmplx_kind), dimension(-p:p, 0:p), intent(in) :: psi
    double precision, dimension(-p:p, 0:p), intent(in) :: y_norm
    double precision, dimension(-p:p, 0:p), intent(out) :: pmn
    double precision, intent(in), dimension(nc, *) :: c
    double precision, intent(in), dimension(*) :: q
    double precision, intent(inout), dimension(*) :: r
    double precision, intent(in) :: x0, y0, z0
    integer :: i, ii, j, jj, k
    double precision :: dx, dy, dz, d, phi, cos_t, t
    complex (kind=dcmplx_kind) :: tc
    intrinsic Real, Cmplx, Atan2, Sqrt, Exp
!
! Dummy statement to "use" q and c
!
    dx = q(1) + c(1,1)
    do ii = 1, n
      i = ind(ii)  ! atom number
      do jj = npoints(i), npoints(i+1) - 1

        dx =  cosurf(1, jj) - x0
        dy =  cosurf(2, jj) - y0
        dz =  cosurf(3, jj) - z0

        d = Sqrt (dx**2 + dy**2 + dz**2)
        cos_t = dz / d
        phi = Atan2 (dy, dx)

        call get_legendre (p, cos_t, pmn)

        t = 0.d0
        do j = 0, p
          tc = (0.d0, 0.d0)
          do k = 1, j
            tc = tc + y_norm(k, j) * pmn(k, j) / d**(j+1) &
                 & * psi(k, j) * Exp (Cmplx (0.d0, k*phi, dcmplx_kind))

          end do
          t = t + 2.d0 * Real (tc) + pmn(0, j) / d**(j+1) * Real(psi(0, j))
        end do

        r(jj) = r(jj) + t
      end do
    end do

  end subroutine bz_far_int

  subroutine bz_mult_int (ind, n, psi, p, x0, y0, z0, y_norm, pmn, c, nc, &
       & q, r)
    !
    ! Computes part of the B*q product for far surface elements,
    ! using local expansions
    !

    use cosmo_C, only: cosurf

    !
    !.. Implicit Declarations ..
    implicit none
    !
    !.. Formal Arguments ..
    integer, intent(in) :: n, p, nc
    integer, intent(in), dimension(*) :: ind
    complex (kind=dcmplx_kind), dimension(-p:p, 0:p), intent(in) :: psi
    double precision, dimension(-p:p, 0:p), intent(in) :: y_norm
    double precision, dimension(-p:p, 0:p), intent(out) :: pmn
    double precision, intent(in), dimension(nc, *) :: c
    double precision, intent(in), dimension(*) :: q
    double precision, intent(inout), dimension(*) :: r
    double precision, intent(in) :: x0, y0, z0
    !
    !.. Local Scalars ..
    integer :: i, ii, j, jj, k
    double precision :: dx, dy, dz, d, phi, cos_t, t
    complex (kind=dcmplx_kind) :: tc
    !
    !.. Intrinsic Functions ..
    intrinsic Real, Atan2, Exp, Sqrt, Cmplx
!
! Dummy statement to "use" q and c
!
    dx = q(1) + c(1,1)
    do ii = 1, n
      i = ind(ii)  ! atom number
      do jj = npoints(i), npoints(i+1) - 1

        dx =  cosurf(1, jj) - x0
        dy =  cosurf(2, jj) - y0
        dz =  cosurf(3, jj) - z0

        d = Sqrt (dx**2 + dy**2 + dz**2)
        cos_t = dz / d
        phi = Atan2 (dy, dx)

        call get_legendre (p, cos_t, pmn)

        t = 0.d0
        do j = 0, p
          tc = (0.d0, 0.d0)
          do k = 1, j
            tc = tc + y_norm(k, j) * pmn(k, j) * d**j &
                 & * psi(k, j) * Exp (Cmplx (0.d0, k*phi, dcmplx_kind))

          end do
          t = t + 2.d0 * Real (tc) + pmn(0, j) * d**j * Real (psi(0, j))
        end do

        r(jj) = r(jj) + t
      end do
    end do

  end subroutine bz_mult_int

  subroutine bp_dir_int (ind1, n1, ind2, n2, c, nc, q, r, same)
    !
    ! Computes part of the B*q product for near surface elements
    !

    use common_arrays_C, only: nfirst, nlast, nat, p
    use cosmo_C, only : cosurf
    use parameters_C, only: tore

    !
    !.. Implicit Declarations ..
    implicit none
    !
    !.. Formal Arguments ..
    integer, intent(in) :: n1, n2, nc
    integer, intent(in), dimension(*) :: ind1, ind2
    double precision, intent(in), dimension(nc, *) :: c
    double precision, intent(in), dimension(*) :: q
    double precision, intent(inout), dimension(*) :: r
    logical, intent(in) :: same
    !
    !.. Local Scalars ..
    integer :: i3, j3, ii, i, j, jj, k, l, m, nao
    double precision :: t
    !
    !.. Local Arrays ..
    double precision :: w(45)
    !
    ! ... Executable Statements ...
!
! Dummy statement to "use" q
!
    t = q(1)
    if(same) then         ! --- points vs. atoms in the same box
      do i3 = 1, n1
        ii = ind1(i3)
        do i = npoints(ii), npoints(ii+1) - 1
          t = 0
          do j3 = 1, n1
            j = ind1(j3)
            nao = nlast(j)-nfirst(j)
            call get_bvec (cosurf(1, i), c(1, j), nao, nat(j), w)

            jj = ijbo_diag(j)

            m = 0
            do k = 1, nao+1
              do l = 1, k-1
                m = m + 1
                t = t - 2*p(jj+m) * w(m)
              end do
              m = m + 1
              t = t - p(jj+m) * w(m)
            end do

            t = t + tore(nat(j)) * w(1)
          end do

          r(i) = r(i) + t
        end do
      end do
    else

      ! --- (1) Points in the 1 box with atoms in 2 box

      do i3 = 1, n1
        ii = ind1(i3)
        do i = npoints(ii), npoints(ii+1) - 1
          t = 0
          do j3 = 1, n2
            j = ind2(j3)
            nao = nlast(j)-nfirst(j)
            call get_bvec (cosurf(1, i), c(1, j), nao, nat(j), w)

            jj = ijbo_diag(j)

            m = 0
            do k = 1, nao+1
              do l = 1, k-1
                m = m + 1
                t = t - 2*p(jj+m) * w(m)
              end do
              m = m + 1
              t = t - p(jj+m) * w(m)
            end do
            t = t + tore(nat(j)) * w(1)
          end do

          r(i) = r(i) + t
        end do
      end do

      ! --- (2) Points in the 2 box with atoms in 1 box

      do i3 = 1, n2
        ii = ind2(i3)
        do i = npoints(ii), npoints(ii+1) - 1
          t = 0
          do j3 = 1, n1
            j = ind1(j3)
            nao = nlast(j)-nfirst(j)
            call get_bvec (cosurf(1, i), c(1, j), nao, nat(j), w)

            jj = ijbo_diag(j)

            m = 0
            do k = 1, nao+1
              do l = 1, k-1
                m = m + 1
                t = t - 2*p(jj+m) * w(m)
              end do
              m = m + 1
              t = t - p(jj+m) * w(m)
            end do
            t = t + tore(nat(j)) * w(1)
          end do

          r(i) = r(i) + t
        end do
      end do
    end if

  end subroutine bp_dir_int

  subroutine cgm_solve (x, n, b, r, p, q, z, m, new_x, new_points)
    !
    !       Conjugated gradient with preconditioner solver
    !       for the system of linear equations
    !       From: "Templates for the solution of linear systems:
    !               Building blocks for iterative methods",
    !               SIAM, Philadelphia, 1994
    !

    use molkst_C, only : numat
    use common_arrays_C, only: coord
    use cosmo_C, only : iatsp, nar_csm, cosurf, srad, nsetf, lenabc, ioldcv

    implicit none
    integer, intent(in) :: n
    logical, intent(in) :: new_x, new_points
    double precision, dimension(n), intent(in) :: b
    double precision, dimension(n), intent(inout) :: x
    double precision, dimension(n), intent(out) :: r, p, q, z
    double precision, dimension(*), intent(out) :: m
    integer, parameter :: max_iters = 100
    double precision, parameter :: stop_tol = 1.d-6
    double precision, parameter :: start_tol = 1.d-2
    integer :: i, maxrs, ierr
    double precision :: t, tt, rho, rho_old, alpha, beta
    double precision, save :: current_tol
    double precision, external :: ddot
    rho_old = 0.d0
    new_iteration = .true.

    maxrs = 60*numat

    if (new_points) then
      current_tol = start_tol
!
      call amat_diag (coord, srad, numat, cosurf, n, nar_csm, nsetf, &
         & nset, rsc, nipsrs, iatsp, tm, ioldcv, maxrs, lenabc, a_diag)
    end if

    if (n > na1max .or. n > na2max) then
      call set_tesselation (surface_handle, ierr)
      if (ierr /= 0) then
        call mopend("Internal error in cgm_solve")
        return
      end if

      if (compute_a_part .and. new_points) then
        i = count_short_ints  (cosurf, 4, simulate_aq_dir_int, .true.)
      end if

    else

      if (compute_a_part .and. new_points) then
        call simulate_aq_vec (coord, srad, numat, cosurf, n, nar_csm, &
       & nsetf, nset, rsc, nipsrs, iatsp, tm, ioldcv, maxrs, lenabc, &
       & .true., i)
      end if


    end if
    if (new_points) then
      call precondition  (cosurf, n, iatsp, numat, a_diag, m)
    end if

    if(new_x) then        ! initial guess
      call precondition_solve (m, x, b, n)
    end if

    if (n > na2max) then

      call afmm (cosurf, 4, n, x, q, n, aq_dir_int, aq_mult_int)
      q = q + a_diag * x

    else if (n > na1max) then

      call simple_mm (cosurf, 4, n, x, q, aq_dir_int, aq_far_int)
      q = q + a_diag * x

    else

      call aq_vec (coord, srad, numat, cosurf, n, &
           & nar_csm, nsetf, nset, rsc, nipsrs, iatsp, &
           & tm, ioldcv, maxrs, lenabc, a_diag, x, q)
    end if

    r = b - q

    if (c_proc < 0.2d0) then
      t = start_tol
    else if (c_proc < 0.4d0) then
      t = start_tol * 0.1d0
    else if (c_proc < 0.6d0) then
      t = start_tol * 0.01d0
    else if (c_proc < 0.8d0) then
      t = start_tol * 0.001d0
    else
      t = stop_tol
    end if

    t = Min (t, current_tol)
    current_tol = t
!
!   write(26,*) ' Tol = ', t
!
    do i = 1, max_iters

      new_iteration = .true.

      call precondition_solve (m, z, r, n)

      rho = ddot (n, r, 1, z, 1)

      if (i == 1) then
        p = z
      else
        beta = rho / rho_old
        p = z + beta * p
      end if

      if (n > na2max) then
        call afmm (cosurf, 4, n, p, q, n, aq_dir_int, aq_mult_int)
        q = q + a_diag * p
      elseif (n > na1max) then
        call simple_mm (cosurf, 4, n, p, q, aq_dir_int, aq_far_int)
        q = q + a_diag * p
      else
        call aq_vec (coord, srad, numat, cosurf, n, &
             & nar_csm, nsetf, nset, rsc, nipsrs, &
             & iatsp, tm, ioldcv, maxrs, lenabc, a_diag, p, q)
      end if

      alpha = rho / ddot (n, p, 1, q, 1)
      x = x + alpha * p
      r = r - alpha * q

      tt = some_norm(r, n)
!
!     write(26,*) ' Micro It. ', i, tt
!
      if (tt < t) then
   !     write (iw, *) i, " microiterations, err = ", tt
        exit
      end if

      rho_old = rho
    end do

  ! if (tt > t) then
  !    write (iw, *) " Warning! cgm_solve finished with error = ", tt
  !  end if

  end subroutine cgm_solve

  double precision function some_norm (v, n)

    !
    !.. Implicit Declarations ..
    implicit none
    !
    !.. Formal Arguments ..
    integer, intent(in) :: n
    double precision, dimension(n), intent(in) :: v
    !
    !.. Local Scalars ..
    integer :: i
    double precision :: t
    !
    !.. Intrinsic Functions ..
    intrinsic Abs
    !
    ! ... Executable Statements ...

    t = Abs (v(1))
    do i = 2, n
      if (t < Abs (v(i))) then
        t = Abs (v(i))
      end if
    end do

    some_norm = t

  end function some_norm

  subroutine precondition (cosurf, nps, iatsp, numat, a_diag, m)

    use cosmo_C, only : disex2
    use chanel_C, only : iw

    !
    !.. Implicit Declarations ..
    implicit none
    !
    !.. Formal Arguments ..
    integer, intent(in) :: nps, numat
    double precision, dimension(4, nps), intent(in) :: cosurf
    integer, dimension(*), intent(in) :: iatsp
    double precision, dimension(nps), intent(in) :: a_diag

    double precision, dimension(*), intent(out) :: m
    !
    !.. Local Scalars ..
    integer :: i, j, ii, jj, ijpos, ips, jps, ierr
    double precision :: aa, d2
    !
    !.. Local Arrays ..
    double precision xa(3), xb(3)
  !  integer, external :: count_short_ints
    !
    ! ... Executable Statements ...


    if (.not. use_a_blocks) then
      m(1:nps) = 1.d0 / a_diag(1:nps)
      return
    end if

    !
    ! Build block-preconditioner
    !

    if (nps > na2max .or. nps > na1max) then
      i = count_short_ints  (cosurf, 4, precond_aq_dir_int, .false.)
      do i=1, nps
        j = iatsp(i)
        ii = iblock_pos(j) -1
        jj = i - npoints(j) +1
        m(ii + (jj*(jj+1))/2) = a_diag(i)
      end do

    else

      ijpos = 0

      do ii = 1, nps
        i = iatsp(ii)
        xa(1:3) = cosurf(1:3, ii)

        do jj = 1, ii
          j = iatsp(jj)
          xb(1:3) = cosurf(1:3, jj)

          d2 = (xb(1)-xa(1))**2 + (xb(2)-xa(2))**2 + &
              & (xb(3)-xa(3))**2

          if (d2 <= disex2 .and. ii /= jj) then
            ijpos = ijpos +1
          end if

          if (i == j) then ! same atom

            if (ii == jj) then ! same point
               aa = a_diag(jj)
            else if (d2 > disex2) then
               aa = 1.d0 / Sqrt (d2)
            else
               aa = a_part(ijpos)
            end if

            ips = ii - npoints(i) +1
            jps = jj - npoints(i) +1


            m(iblock_pos(i) -1 + (ips*(ips-1))/2 + jps) = aa

          end if

        end do
      end do

    end if

    !
    ! Invert blocks
    !

    do i=1, numat
       ii = npoints(i+1) - npoints(i)
       if (ii > 0) then
         ijpos = iblock_pos(i)
         do j=1,ii
           do jj=1,j
             a_block(jj, j) = m(ijpos)
             ijpos = ijpos +1
           end do
         end do

         call dpotrf ('U', ii, a_block, max_block_size, ierr)
      !   if (ierr /= 0) then
      !     write(iw,*) ' Error in decomposition of atom ', i, ' = ', ierr
      !     call mopend(' Internal error')
      !     return
      !  end if

         call dpotri ('U', ii, a_block, max_block_size, ierr)
         if (ierr /= 0) then
           write(iw,*) ' Error in Matrix invert ', i, ' = ', ierr
           call mopend(' Internal error')
           return
         end if

       ! !
       ! ! === Test ===
       ! !
       !
       ! do j=1,ii
       !   do jj=1,j
       !     a_block( j, jj ) = a_block( jj, j )
       !   end do
       ! end do
       !
       ! do j=1,ii
       !   call mult_triangle_vec( m( iblock_pos(i)), a_block(1,j), ii, &
       !        q_vec )
       !
       !   do jj=1,ii
       !     if( jj == j ) then
       !        if( abs( q_vec( jj ) -1.d0 ) > 1.d-6 ) then
       !           write(26,*) ' Error diag', i, j, jj
       !        end if
       !     else
       !        if( abs( q_vec( jj )) > 1.d-6 ) then
       !           write(26,*) ' Error off-diag', i, j, jj
       !        end if
       !     end if
       !   end do
       ! end do
       !
       ! !
       ! ! === end Test ===
       ! !

         ijpos = iblock_pos(i)
         do j=1,ii
           do jj=1,j
             m(ijpos) = a_block(jj, j)
             ijpos = ijpos +1
           end do
         end do

       end if
     end do

  end subroutine precondition

  subroutine amat_diag (coord, srad, numat, cosurf, nps, nar_csm, nsetf, &
       & nset, rsc, nipsrs, iatsp, tm, ioldcv, maxrs, lenabc, m)
    !
    ! Computes diagonal elements of the A matrix
    !

    use cosmo_C, only : dirvec
    use funcon_C, only: pi
    implicit none
    integer, intent(in) :: numat, nps, ioldcv, maxrs, lenabc

    double precision, dimension(3, numat), intent(in) :: coord
    double precision, dimension(numat), intent(in) :: srad
    double precision, dimension(4, nps), intent(in) :: cosurf
    integer, dimension(lenabc+1), intent(in) :: nar_csm, nsetf, iatsp
    integer, dimension(*), intent(in) :: nset
    double precision, dimension(4, maxrs), intent(in) :: rsc
    integer, dimension(lenabc), intent(in) :: nipsrs
    double precision, dimension(3, 3, numat), intent(in) :: tm

    double precision, dimension(nps), intent(out) :: m
    integer :: i, ii, k, l, nfl1
    double precision :: aa, fdiag, ri, x1, x2, x3, x4
    double precision, dimension(4, 500, 2) :: finel
    double precision xi(3)
    fdiag = 2.1d0 * Sqrt (Pi)

    do ii = 1, nps
      i = iatsp(ii)
      ri = srad(i)
      xi(1:3) = coord(1:3, i)

      call mfinel (ii, 1, finel, nar_csm, nsetf, nset, rsc, &
           & nipsrs, dirvec, tm(1, 1, i), xi, ri, &
           & nfl1, ioldcv, maxrs, lenabc, numat)
      !
      !       Diagonal element
      !
      aa = 0.d0
      do k = 1, nfl1
        aa = aa + fdiag * Sqrt (finel(4, k, 1)**3)
        x1 = finel(1, k, 1)
        x2 = finel(2, k, 1)
        x3 = finel(3, k, 1)
        x4 = finel(4, k, 1)

        do l = 1, k - 1
          aa = aa + 2 * x4 * finel(4, l, 1) / &
               &  Sqrt ((x1-finel(1, l, 1))**2 + &
               & (x2-finel(2, l, 1))**2 + &
               & (x3-finel(3, l, 1))**2)
        end do
      end do

      m(ii) = aa / (cosurf(4, ii) **2)
    end do

  end subroutine amat_diag

  subroutine precondition_solve (m, z, r, n)

    use molkst_C, only : numat

    !
    !.. Implicit Declarations ..
    implicit none
    integer, intent (in) :: n
    double precision, dimension(n), intent(in) :: r
    double precision, dimension(*), intent(in) :: m
    double precision, dimension(n), intent(out) :: z
    integer :: i, ndim, ipos
    if (use_a_blocks) then
      ipos = 1
      do i=1, numat
        ndim = npoints(i+1) - npoints(i)
        if (ndim > 0) then
          call mult_triangle_vec (m(iblock_pos(i)) , r(ipos) , ndim, &
            & z(ipos))
          ipos = ipos + ndim
        end if
      end do
    else
      z(1:n) = r(1:n) * m(1:n)
    end if

  end subroutine precondition_solve

  subroutine mult_triangle_vec (t, x, n, y)
    !
    !.. Implicit Declarations ..
    implicit none

    !
    !.. Formal Arguments ..
    integer, intent (in) :: n
    double precision, dimension(*), intent(in) :: t
    double precision, dimension(n), intent(in) :: x
    double precision, dimension(n), intent(out) :: y
    !
    !.. Local Scalars ..
    integer :: i, j, ii, ij
    double precision :: s
    !
    ! ... Executable Statements ...

        ii = 0   ! previous diagonal element
        do i=1, n
          s = 0.d0
          ij = ii
          do j=1, i
            ij = ij +1
            s = s + t(ij) * x(j)
          end do
          ii = ij
          do j=i+1, n
            ij = ij + j -1
            s = s + t(ij) * x(j)
          end do
          y(i) = s
        end do

  end subroutine mult_triangle_vec


  subroutine addnucz (phinet, qscnet, qdenet)
    use molkst_C, only: lm61, numat
    use cosmo_C, only : lenabc, nps, idenat
    use parameters_C, only: tore
    use common_arrays_C, only : nat
    !
    !.. Implicit Declarations ..
    implicit none
    !
    !.. Formal Arguments ..
    double precision, dimension(lenabc+1, 3), intent (inout) :: phinet, &
         & qscnet
    double precision, dimension(lm61, 3), intent (inout) :: qdenet
    !
    !.. Local Scalars ..
    integer :: i
    !
    ! ... Executable Statements ...

    phinet(1:nps, 1) = 0.d0
    qscnet(1:nps, 1) = 0.d0
    qdenet(1:lm61, 1) = 0.d0

    do i = 1, numat
       qdenet(idenat(i), 1) = tore(nat(i))
    end do

  end subroutine addnucz

  subroutine addfckz
    !
    ! Adds correction to the Fock matrix due to charges on surface
    ! Computes corrections to total energy
    !
    use molkst_C, only: lm61, numat
    use cosmo_C, only : nps, solv_energy, fepsi, cosurf, iatsp, phinet, &
      qscnet, qdenet, ipiden, gden, qscat, ediel
    use parameters_C, only: tore
    use funcon_C, only: a0, ev

    use common_arrays_C, only: coord, nfirst, nlast, nat, f, p
    !
    !.. Implicit Declarations ..
    implicit none
    !
    !.. Local Scalars ..
    integer :: i, iat, im, j, ii, nao, nnn, ierr
    double precision :: fcon, qsc3, s1, t
    !.. Local Arrays ..
    integer :: idiag(9) = (/ 1, 3, 6, 10, 15, 21, 28, 36, 45 /)
    double precision :: w(45), v(45)
    !
    !.. Intrinsic Functions ..
    intrinsic Sum
    !
    ! ... Executable Statements ...

    fcon = a0 * ev
    qscat(1:numat) = 0.d0
    phinet(1:nps, 2) = 0.d0
    if (numat > nb1max .or. numat > nb2max) then
      do i = 1, numat
        ii = ijbo_diag(i)
        j = nlast(i) - nfirst(i) + 1
        s1 = Sum(p(ii + idiag(1:j)))
        q_vec(i) = -s1 + tore(nat(i))
      end do

      call set_tesselation (atom_handle, ierr)
      if (ierr == -1000) return ! dummy use of ierr
    end if

    if (numat > nb2max) then
      call afmm (coord, 3, numat, q_vec, phinet(1, 2), nps, &
           & bp_dir_int, bz_mult_int)

    elseif (numat > nb1max) then
      call simple_mm (coord, 3, numat, q_vec, phinet(1, 2), &
           & bp_dir_int, bz_far_int)

    else
      call bz_vec (phinet(1, 1))
      call bpnew_vec (phinet(1, 2))
      phinet(1:nps, 2) = phinet(1:nps, 1) + phinet(1:nps, 2)
      phinet(1:nps, 1) = 0
      !!stop
    end if

    phinet(1:nps, 3) = phinet(1:nps, 2)

    if (new_surface) then
      call cgm_solve (qscnet(1, 2), nps, phinet(1, 2), &
           & r_vec, p_vec, q_vec, z_vec, m_vec, .true., .true.)
      new_surface = .false.
    else
      qscnet(1:nps, 2) = qscnet(1:nps, 2) / (-fepsi)
      call cgm_solve (qscnet(1, 2), nps, phinet(1, 2), &
           & r_vec, p_vec, q_vec, z_vec, m_vec, .false., .false.)
    end if

    ediel = 0.d0
    s1 = 0.0d0
    do i = 1, nps
      iat = iatsp(i)
      qscnet(i, 2) = -fepsi * qscnet(i, 2)
      qsc3 = qscnet(i, 1) + qscnet(i, 2)
      qscnet(i, 3) = qsc3
      ediel = ediel + qsc3 * phinet(i, 3)
      s1 = s1 + qscnet(i, 1)
      qscat(iat) = qscat(iat) + qsc3
    end do
    ediel = ediel * fcon / 2

    if (numat > nb1max .or. numat > nb2max) then
      call set_tesselation (atom_handle, ierr)
      qdenet(1:lm61, 2) = 0
    end if

    if (numat > nb2max) then
      call afmm (coord, 3, numat, qscnet(1, 2), qdenet(1, 2), lm61, &
           & fock_dir_int, fock_mult_int, sphere_f_multipoles)

      f(ipiden(1:lm61)) = f(ipiden(1:lm61)) - qdenet(1:lm61, 2) * fcon

      s1 = 0.d0
      do i = 1, lm61
        s1 = s1 + qdenet(i, 2) * p(ipiden(i)) * gden(i)
      end do
      s1 = -s1 *fcon

    elseif (numat > nb1max) then
      call simple_mm (coord, 3, numat, qscnet(1, 2), qdenet(1, 2), &
           & fock_dir_int, fock_far_int, sphere_f_multipoles)

      f(ipiden(1:lm61)) = f(ipiden(1:lm61)) - &
           & qdenet(1:lm61, 2) * fcon

      s1 = 0.d0
      do i = 1, lm61
        s1 = s1 + qdenet(i, 2) * p(ipiden(i)) * gden(i)
      end do
      s1 = -s1 * fcon

    else

      s1 = 0.d0

      im = 0
      do i = 1, numat
        nao = nlast(i)-nfirst(i)
        ii = ijbo_diag(i)

        if (nao == 0) then ! S - element

          v(1) = 0.d0
          do j=1, nps
            t = 1.d0 / Sqrt (( cosurf(1,j) - coord(1,i) )**2 &
                         &  +( cosurf(2,j) - coord(2,i) )**2 &
                         &  +( cosurf(3,j) - coord(3,i) )**2)
            v(1) = v(1) + qscnet(j, 2) * t
          end do

          v(1) = -v(1) * fcon

          f(ii+1) = f(ii+1) + v(1)

          im = im + 1
          s1 = s1 + v(1) * p(ii+1) * gden(im)

        else

          nnn = ((nao+2)*(nao+1)) / 2
          v(1:nnn) = 0

          do j = 1, nps
            call get_bvec (cosurf(1, j), coord(1, i), nao, nat(i), w)
            v(1:nnn) = v(1:nnn) + w(1:nnn) * qscnet(j, 2)
          end do

          v(1:nnn) = -v(1:nnn) * fcon

          do j = 1, nnn
            f(ii+j) = f(ii+j) + v(j)

            im = im + 1
            s1 = s1 + v(j) * p(ii+j) * gden(im)
          end do

        end if
      end do
    end if

    solv_energy = s1 * 0.5d0 + ediel

    ! CALCULATE QDENET FROM DENSITY MATRIX

    do i = 1, lm61
      qdenet(i, 2) = gden(i) * p(ipiden(i))
      qdenet(i, 3) = qdenet(i, 2) + qdenet(i, 1)
    end do

  end subroutine addfckz

  subroutine fock_dir_int (ind1, n1, ind2, n2, c, nc, q, r, same)
    !
    ! Adds correction to the Fock matrix due to near charges on surface
    !

    use common_arrays_C, only: nfirst, nlast, nat
    use cosmo_C, only : qscnet, cosurf

    !
    !.. Implicit Declarations ..
    implicit none
    !
    !.. Formal Arguments ..
    integer, intent(in) :: n1, n2, nc
    integer, intent(in), dimension(*) :: ind1, ind2
    double precision, intent(in), dimension(nc, *) :: c
    double precision, intent(in), dimension(*) :: q
    double precision, intent(inout), dimension(*) :: r
    logical, intent(in) :: same
    !
    !.. Local Scalars ..
    integer :: i3, j3, ii, i, j, jj, nao, nnn
    !
    !.. Local Arrays ..
    double precision :: v(45), w(45)

!
! Dummy statement to "use" q
!
    v(1) = q(1)
    if(same) then         ! --- Atoms vs. points in the same box
      do i3 = 1, n1
        i = ind1(i3)
        nao = nlast(i)-nfirst(i)
        nnn = ((nao+2)*(nao+1)) / 2
        v(1:nnn) = 0
        ii = iatom_pos(i)

        do j3 = 1, n1
          jj = ind1(j3)

          do j = npoints(jj), npoints(jj+1) - 1
            call get_bvec (cosurf(1, j), c(1, i), nao, nat(i), w)
            v(1:nnn) = v(1:nnn) + w(1:nnn) * qscnet(j, 2)
          end do
        end do

        do j = 1, nnn
          r(ii+j) = r(ii+j) + v(j)
        end do
      end do
    else
      ! --- (1) Atoms n the 1 box with points in 2 box

      do i3 = 1, n1
        i = ind1(i3)
        nao = nlast(i)-nfirst(i)
        nnn = ((nao+2)*(nao+1)) / 2
        v(1:nnn) = 0
        ii = iatom_pos(i)

        do j3 = 1, n2
          jj = ind2(j3)

          do j = npoints(jj), npoints(jj+1) - 1
            call get_bvec (cosurf(1, j), c(1, i), nao, nat(i), w)

            v(1:nnn) = v(1:nnn) + w(1:nnn) * qscnet(j, 2)
          end do
        end do

        do j = 1, nnn
          r(ii+j) = r(ii+j) + v(j)
        end do
      end do

      ! --- (2) Atoms in the 2 box with points in 1 box

      do i3 = 1, n2
        i = ind2(i3)
        nao = nlast(i)-nfirst(i)
        nnn = ((nao+2)*(nao+1)) / 2
        v(1:nnn) = 0
        ii = iatom_pos(i)

        do j3 = 1, n1
          jj = ind1(j3)

          do j = npoints(jj), npoints(jj+1) - 1
            call get_bvec (cosurf(1, j), c(1, i), nao, nat(i), w)

            v(1:nnn) = v(1:nnn) + w(1:nnn) * qscnet(j, 2)
          end do
        end do

        do j = 1, nnn
          r(ii+j) = r(ii+j) + v(j)
        end do
      end do
    end if

  end subroutine fock_dir_int

  subroutine fock_mult_int (ind, n, psi, p, x0, y0, z0, y_norm, pmn, c, &
       & nc, q, r)
    !
    ! Adds correction to the Fock matrix due to far charges on surface
    ! using local expansions
    !

    use common_arrays_C, only: nfirst, nlast
    !
    !.. Implicit Declarations ..
    implicit none
    !
    !.. Formal Arguments ..
    integer, intent(in) :: n, p, nc
    integer, intent(in), dimension(*) :: ind
    complex (kind=dcmplx_kind), dimension(-p:p, 0:p), intent(in) :: psi
    double precision, dimension(-p:p, 0:p), intent(in) :: y_norm
    double precision, dimension(-p:p, 0:p), intent(out) :: pmn
    double precision, intent(in), dimension(nc, *) :: c
    double precision, intent(in), dimension(*) :: q
    double precision, intent(inout), dimension(*) :: r
    double precision, intent(in) :: x0, y0, z0
    !
    !.. Local Scalars ..
    integer :: i, ii, i3, j, k, nnn
    double precision :: dx, dy, dz, d, phi, cos_t, t
    complex (kind=dcmplx_kind) :: tc
    !
    !.. Local Arrays ..
    integer :: idiag(9) = (/ 1, 3, 6, 10, 15, 21, 28, 36, 45 /)
    !
    !.. Intrinsic Functions ..
    intrinsic Real, Cmplx, Atan2, Exp, Sqrt
    !
    ! ... Executable Statements ...

!
! Dummy statement to "use" q
!
    dx = q(1)
    do i3 = 1, n
      i = ind(i3)  ! atom number
      ii = iatom_pos(i)
      nnn = nlast(i) - nfirst(i) + 1

      dx =  c(1, i) - x0
      dy =  c(2, i) - y0
      dz =  c(3, i) - z0

      d = Sqrt (dx**2 + dy**2 + dz**2)
      cos_t = dz / d
      phi = Atan2 (dy, dx)

      call get_legendre (p, cos_t, pmn)

      t = 0.d0
      do j = 0, p
        tc = Cmplx (0.d0, 0.d0, dcmplx_kind)
        do k = 1, j
          tc = tc + y_norm(k, j) * pmn(k, j) * d**j &
               & *psi(k, j) * Exp (Cmplx (0.d0, k*phi, dcmplx_kind))

        end do

        t = t + 2.d0 * Real (tc) + pmn(0, j) * d**j * Real (psi(0, j))
      end do

      r(ii + idiag(1:nnn)) = r(ii + idiag(1:nnn)) + t

    end do

  end subroutine fock_mult_int

  subroutine fock_far_int (ind, n, psi, p, x0, y0, z0, y_norm, pmn, c, &
       & nc, q, r)
    !
    ! Adds correction to the Fock matrix due to far charges on surface
    ! using multipole expansions
    !

    use common_arrays_C, only: nfirst, nlast

    !
    !.. Implicit Declarations ..
    implicit none
    !
    !.. Formal Arguments ..
    integer, intent(in) :: n, p, nc
    integer, intent(in), dimension(*) :: ind
    complex (kind=dcmplx_kind), dimension(-p:p, 0:p), intent(in) :: psi
    double precision, dimension(-p:p, 0:p), intent(in) :: y_norm
    double precision, dimension(-p:p, 0:p), intent(out) :: pmn
    double precision, intent(in), dimension(nc, *) :: c
    double precision, intent(in), dimension(*) :: q
    double precision, intent(inout), dimension(*) :: r
    double precision, intent(in) :: x0, y0, z0
    !
    !.. Local Scalars ..
    integer :: i, ii, i3, j, k, nnn
    double precision :: dx, dy, dz, d, phi, cos_t, t
    complex (kind=dcmplx_kind) :: tc
    !
    !.. Local Arrays ..
    integer :: idiag(9) = (/ 1, 3, 6, 10, 15, 21, 28, 36, 45 /)
    !
    !.. Intrinsic Functions ..
    intrinsic Real, Cmplx, Atan2, Sqrt, Exp
    !
    ! ... Executable Statements ...

!
! Dummy statement to "use" q
!
    t = q(1)
    do i3 = 1, n
      i = ind(i3)  ! atom number
      ii = iatom_pos(i)
      nnn = nlast(i) - nfirst(i) + 1

      dx =  c(1, i) - x0
      dy =  c(2, i) - y0
      dz =  c(3, i) - z0

      d = Sqrt (dx**2 + dy**2 + dz**2)
      cos_t = dz / d
      phi = Atan2 (dy, dx)

      call get_legendre (p, cos_t, pmn)

      t = 0.d0
      do j = 0, p
        tc = (0.d0, 0.d0)
        do k = 1, j
          tc = tc + y_norm(k, j) * pmn(k, j) / d**(j+1) &
               & * psi(k, j) * Exp (Cmplx (0.d0, k*phi, dcmplx_kind))

        end do

        t = t + 2.d0 * Real (tc) + pmn(0, j) / d**(j+1) * Real (psi(0, j))
      end do

      r(ii + idiag(1:nnn)) = r(ii + idiag(1:nnn)) + t

    end do

  end subroutine fock_far_int

  subroutine sphere_f_multipoles (x0, y0, z0, ind, n, pmn, y_norm, &
       & c, nc, q, psi, p)
    !
    ! Computes multipole expansions of surface charges
    !

    use cosmo_C, only: qscnet, cosurf

    !
    !.. Implicit Declarations ..
    implicit none
    !
    !.. Formal Arguments ..
    double precision, intent(in) :: x0, y0, z0
    integer, intent(in) :: n, nc, p
    integer, dimension(*), intent(in) :: ind
    double precision, dimension(nc, *), intent(in) :: c
    double precision, dimension(*), intent(in) :: q
    double precision, dimension(-p:p, 0:p), intent(in) :: y_norm
    double precision, dimension(-p:p, 0:p), intent(inout) :: pmn

    complex (kind=dcmplx_kind), dimension(-p:p, 0:p), intent(inout) :: psi
    !
    !.. Local Scalars ..
    double precision :: dx, dy, dz, t, tt, r, cos_t, phi, qq
    complex (kind=dcmplx_kind) :: tc, tc1
    integer :: j, j1, k, n1, m1
    !
    !.. Intrinsic Functions ..
    intrinsic Sqrt, Atan2, Cmplx, Exp, Conjg
    !
    ! ... Executable Statements ...

!
! Dummy statement to "use" q and c
!
    t = q(1) + c(1,1)
    do j = 1, n
      j1 = ind(j)   ! atom number
      do k = npoints(j1), npoints(j1+1) - 1
        dx = cosurf(1, k) - x0
        dy = cosurf(2, k) - y0
        dz = cosurf(3, k) - z0
        qq = qscnet(k, 2)
        r = Sqrt (dx**2 + dy**2 + dz**2)
        cos_t = dz / r
        phi = Atan2 (dy, dx)

        call get_legendre (p, cos_t, pmn)

        psi(0, 0) = psi(0, 0) + Cmplx(qq, 0.d0, dcmplx_kind)

        t = 1.d0
        do n1 = 1, p
          t = t * r
          do m1 = 0, n1
            tc = Exp (Cmplx (0.d0, -m1*phi, dcmplx_kind))
            tt = y_norm(-m1, n1) * pmn(-m1, n1) * t * qq
            tc1 = tc * tt
            psi(m1, n1) = psi(m1, n1) + tc1
          end do
        end do

        do n1 = 1, p
          do m1 = 1, n1
            psi(-m1, n1) = Conjg (psi(m1, n1))
          end do
        end do

      end do
    end do

  end subroutine sphere_f_multipoles



  subroutine am1dft_solve (qsct, phit, nps)
    ! Wrapper for cgm_solve to calculate AM1DFT correction in writmn.

    !.. Implicit Declarations ..
    implicit none

    !
    !.. Formal Arguments ..
    integer, intent(in) :: nps
    double precision, dimension(nps), intent(inout) :: qsct
    double precision, dimension(nps), intent(in) :: phit
    !
    !
    ! ... Executable Statements ...

    call cgm_solve (qsct, nps, phit, &
         & r_vec, p_vec, q_vec, z_vec, m_vec, .false., .false.)
  end subroutine am1dft_solve

end module linear_cosmo
