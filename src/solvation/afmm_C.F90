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

module afmm_C
  !
  !
  !       Adaptive multipole algorithm
  !       Based on the paper
  !       H. Cheng, L. Greengard, and V. Rokhlin
  !       "A fast adaptive multipole algorithm in three dimensions",
  !       J. Comp. Phys., 155, 468-498 (1999)
  !
  !       But use spherical harmonic addition theorem for
  !       Multipole -> Local transformation, no rotations
  !
  !       L1 list - modified
  !
  use chanel_C, only: iw
  implicit none

  public :: afmm_ini, prepare_tesselations, set_tesselation, &
       & divide_box, afmm, simple_mm, get_legendre, &
       & count_short_ints

  private

  logical, parameter :: DEBUG = .false.
  logical, parameter :: PRINT = .false.

  double precision, parameter ::  small = 1.d-3   ! Tolerance for
                                                  ! for box ident.

  integer, parameter ::   p = 3         ! Multipole expansion level >= 2
  integer, parameter :: dcmplx_kind = Kind ((1.d0, 1.d0))

  type Box
     integer :: parent, level
     integer :: n_p, p_start
     integer :: n_childs, n_l1, n_l2
     integer :: n_up, n_down, n_north, n_south, n_east, n_west
     integer, dimension(8) :: childs
     integer, dimension(28) :: l1
     integer, dimension(189) :: l2

     double precision :: x0, y0, z0

     complex(kind=dcmplx_kind), dimension(-p:p, 0:p) :: phi, psi
  end type Box

  type(Box), dimension(:), pointer :: b
  double precision, dimension(:), pointer :: d
  integer, dimension(:), pointer :: level, ind, ind1

  integer :: nbox, nlev
  integer, dimension(:), allocatable :: points_index

  double precision, dimension(-p:p, 0:p) :: pmn, y_norm, amn
  double precision, dimension(0:2*p) :: fact

  type Tesselation
     type(Box), dimension(:), pointer :: b
     double precision, dimension(:), pointer :: d
     integer, dimension(:), pointer :: level, ind, ind1
     integer :: nbox, nlev, n
  end type Tesselation

  integer :: max_tess, current_tess

  Type(Tesselation), dimension(:), allocatable :: tess

contains

  subroutine prepare_tesselations (n1, n2, ind, n, ierr)

    !
    !.. Implicit Declarations ..
    implicit none
    !
    !.. Formal Arguments ..
    integer, intent(in) :: n1, n2, n
    integer, intent(out) :: ierr
    integer, dimension(n), intent(out) :: ind
    !
    !.. Local Scalars ..
    integer :: i, j
    !
    !.. Intrinsic Functions ..
    intrinsic Allocated, Associated, Max
    !
    ! ... Executable Statements ...

    if (Allocated (tess)) then
      do j = 1, max_tess
        if (Associated (tess(j)%b) ) then
          deallocate (tess(j)%b, tess(j)%d, tess(j)%level, tess(j)%ind, &
               & tess(j)%ind1, stat = i)

          if (i /= 0) then
            ierr = -3
            return
          end if
        end if
      end do

      deallocate (tess, points_index, stat = i)
      if (i /= 0) then
        ierr = -2
        return
      end if
    end if

    allocate (tess(n), points_index(Max (n1, n2)), stat = i)

    if (i /= 0) then
      ierr = -1
      return
    end if

    do j = 1, n
      nullify(tess(j)%b, tess(j)%d, tess(j)%level, tess(j)%ind, &
           & tess(j)%ind1)
    end do

    do i = 1, n
      ind(i) = i
    end do

    max_tess = n
    current_tess = 0
    ierr = 0

  end subroutine prepare_tesselations

  subroutine set_tesselation (handle, ierr)
    !
    !.. Implicit Declarations ..
    implicit none
    !
    !.. Formal Arguments ..
    integer, intent(in) :: handle
    integer, intent(out) :: ierr
    !
    ! ... Executable Statements ...

    if (handle > 0 .and. handle <= max_tess) then
      ierr = 0

      if (current_tess /= handle) then
        current_tess = handle
        b => tess(handle)%b
        d => tess(handle)%d
        level => tess(handle)%level
        ind => tess(handle)%ind
        ind1 => tess(handle)%ind1
        nbox = tess(handle)%nbox
        nlev = tess(handle)%nlev
      end if

      points_index(1:tess(handle)%n) = ind(1:tess(handle)%n)

    else
      ierr = -1
    end if


  end subroutine set_tesselation

  subroutine store_tesselation (handle, n)
    !
    !.. Implicit Declarations ..
    implicit none
    !
    !.. Formal Arguments ..
    integer, intent(in) :: handle, n
    !
    ! ... Executable Statements ...

    current_tess = handle
    tess(handle)%b => b
    tess(handle)%d => d
    tess(handle)%level => level
    tess(handle)%ind => ind
    tess(handle)%ind1 => ind1
    tess(handle)%nbox = nbox
    tess(handle)%nlev = nlev
    tess(handle)%n = n

  end subroutine store_tesselation

  subroutine divide_box (coord, n1, n, min_d, min_p, max_tes, handle, ierr)
    implicit none
    integer, intent (in) :: n1, n, handle
    double precision, intent (in), dimension(n1, n) :: coord
    double precision, intent (in) :: min_d   ! Minimum box size
    integer, intent (in) :: min_p   ! Minimum points in cell
    integer, intent (in) :: max_tes ! Max tesselation size
                                    ! (if < 3, then compute)

    integer, intent(out) :: ierr  ! error, 0 if everything OK
    double precision :: xmax, xmin, ymax, ymin, zmax, zmin, &
         &  dx, dy, dz, x0, y0, z0, t, r
    integer :: max_boxs, max_lev
    integer :: i, i0, i1, j, j0, k, m, status
    logical :: division
    intrinsic Sqrt, Abs, Dble, Int, Log, Max, Min
    !
    !       Step -1. Find initial box
    !       ======
    !
    if (n1 < 3) then
      ierr = -10
      return
    end if

    xmax = coord(1, 1)
    xmin = coord(1, 1)
    ymax = coord(2, 1)
    ymin = coord(2, 1)
    zmax = coord(3, 1)
    zmin = coord(3, 1)

    do i = 2, n
      if (xmax < coord(1, i)) then
        xmax = coord(1, i)
      elseif (xmin > coord(1, i)) then
        xmin = coord(1, i)
      end if

      if (ymax < coord(2, i)) then
        ymax = coord(2, i)
      elseif (ymin > coord(2, i)) then
        ymin = coord(2, i)
      end if

      if (zmax < coord(3, i)) then
        zmax = coord(3, i)
      elseif (zmin > coord(3, i)) then
        zmin = coord(3, i)
      end if
    end do

    dx = xmax - xmin
    dy = ymax - ymin
    dz = zmax - zmin

    t = dx
    if (t < dy) then
      t = dy
    end if
    if (t < dz) then
      t = dz
    end if

    ! --- adjust the initial box

    x0 = (xmax + xmin) /2
    y0 = (ymax + ymin) /2
    z0 = (zmax + zmin) /2
    t = t /2

    xmin = x0 - t
    xmax = x0 + t
    ymin = y0 - t
    ymax = y0 + t
    zmin = z0 - t
    zmax = z0 + t

    if (DEBUG) then
      write (iw, *) " Adjusted, t = ", t
      write (iw, *) "xmax, xmin", xmax, xmin
      write (iw, *) "ymax, ymin", ymax, ymin
      write (iw, *) "zmax, zmin", zmax, zmin
    end if

    if (max_tes < 3) then
      max_lev = Int (Min (Log(t / min_d) / Log(2.d0) + 1.d0, &
           & Log(n / Dble (min_p)) / Log(8.d0) + 2.d0))
    else
      max_lev = max_tes
    end if

    if (max_lev < 3) max_lev = 3

    max_boxs = 0
    j = 1
    do i = 1, max_lev
      max_boxs = max_boxs + j*8
      j = j * 8
    end do

    if (PRINT) then
      write (iw, *) " Estimated number of splits =", max_lev
      write (iw, *) " Estimated number of boxes =", max_boxs
    end if

    allocate (b(max_boxs), d(max_lev), level(max_lev+1), &
         & ind(n+1), ind1(Max(n+1, 189*6)), stat = status)

    if (status /= 0) then
      ierr = -2
      return
    end if

    !
    !       Step 0. Split the initial box
    !       ======
    !
    do i = 1, n
      ind(i) = i
    end do

    call split_box (0, 1, x0, y0, z0, t*0.5, coord, n1, n, 0, 0, m)

    d(1) = t
    level(1) = 1
    level(2) = m + 1
    nbox = m

    do i = 1, max_lev - 1
      division = .false.
      t = t / 2
      do j = level(i), level(i+1) - 1
        if (b(j)%n_p > min_p) then
          call split_box (j, i+1, b(j)%x0, b(j)%y0, b(j)%z0, &
               &    t*0.5, coord, n1, b(j)%n_p, nbox, b(j)%p_start - 1, m)
          division = .true.
          b(j)%n_childs = m
          do k = 1, m
            b(j)%childs(k) = nbox + k
          end do
          nbox = nbox + m
        end if
      end do

      if (division) then
        nlev = i+1
        d(i+1) = t
        level(i+2) = nbox + 1

        do j = level(i+1), level(i+2) - 1
          b(j)%p_start = b(j)%p_start + b(b(j)%parent)%p_start - 1
        end do
      else
        nlev = i
        exit
      end if
    end do

    if (PRINT) then
      write (iw, *) " Actual number of splits =", nlev
      write (iw, *) " Actual number of boxes =", nbox
      do i=1, nlev
        write (iw, *) " Level ", i, "Size = ", d(i)
      end do
    end if

    !
    !       Create lists for each box
    !
    !       l1(b) - set consisting of all boxes adjustent to b.
    !       Let c be in l1(b). Then
    !               if c has children it will be in l1(b)
    !               if c is childless, then only c with numbers < b
    !                       will be in the l1 list of b.
    !       This is done to avoid double counting on step 7.
    !       List l1 is computed for parent boxes, and contains all the boxes
    !       (Note 1. So the list 1 has boxes on the same split level,
    !               as box b, similar to list 2).
    !       (Note 2. in this version the definition of list 1 includes
    !               lists 1 and 3 of the original paper)
    !

    do i = 1, nlev
      t = d(i) + small
      do j = level(i), level(i+1) - 1
        do k = level(i), level(i+1) - 1

          if (j == k) cycle

          dx = Abs (b(j)%x0 - b(k)%x0)
          dy = Abs (b(j)%y0 - b(k)%y0)
          dz = Abs (b(j)%z0 - b(k)%z0)

          if (dx < t .and. dy < t .and. dz < t) then
            if (b(j)%n_childs /= 0) then
              b(j)%n_l1 = b(j)%n_l1 + 1
              b(j)%l1(b(j)%n_l1) = k
            elseif (j > k .or. b(k)%n_childs /= 0) then
              b(j)%n_l1 = b(j)%n_l1 + 1
              b(j)%l1(b(j)%n_l1) = k
            end if
          end if
        end do
      end do
    end do
    !
    !       l2(b) - set consisting of all children of boxes
    !       adjustent to parent of b, but separated from b.
    !       (Note. The list 2 is empty for boxes on level 1)
    !
    do i = level(2), nbox
      t = d(b(i)%level) + small

      i0 = b(i)%parent
      do j = 1, b(i0)%n_l1
        j0 = b(i0)%l1(j)
        do k = 1, b(j0)%n_childs
          i1 = b(j0)%childs(k)

          dx = Abs (b(i)%x0 - b(i1)%x0)
          dy = Abs (b(i)%y0 - b(i1)%y0)
          dz = Abs (b(i)%z0 - b(i1)%z0)

          if (.not. (dx < t .and. dy < t .and. dz < t)) then
            b(i)%n_l2 = b(i)%n_l2 + 1
            b(i)%l2(b(i)%n_l2) = i1
          end if
        end do
      end do
    end do

    if (DEBUG) then
      write (iw, *) " Boxes:"
      ind1 = 0
      do i = 1, nbox
        write (iw, *) "Box, level, x0, y0, z0", i, &
             & b(i)%level, b(i)%x0, b(i)%y0, b(i)%z0
        write (iw, *) "n_childs, n_l1, n_l2, n_p", &
             & b(i)%n_childs, b(i)%n_l1, b(i)%n_l2, b(i)%n_p

        write (iw, *) "parent: ", b(i)%parent

        write (iw, *) " Childs list:"
        write (iw, *) (b(i)%childs(k), k = 1, b(i)%n_childs)

        write (iw, *) " L1 list:"
        write (iw, *) (b(i)%l1(k), k = 1, b(i)%n_l1)

        write (iw, *) " L2 list:"
        write (iw, *) (b(i)%l2(k), k = 1, b(i)%n_l2)

        if (b(i)%n_childs == 0) then
          write (iw, *) " Points:"
          write (iw, *) (ind(b(i)%p_start + k), k = 0, b(i)%n_p - 1)

          do k = 0, b(i)%n_p - 1
            j = b(i)%p_start + k
            if (ind1(j) /= 0) then
              write (iw, *) "--- ERROR in points ---"
            else
              ind1(j) = ind(j)
            end if
          end do
        end if
      end do

      do i = 1, n
        if (ind1(i) == 0) then
          write (iw, *) " +++ Error in points "
        end if
      end do

      do i = 1, nlev
        write (iw, *) " Boxes on level ", i, " : ", &
             & level(i), level(i+1) - 1
      end do
    end if

    if (DEBUG) then

    !
    ! --- Test for points assignment
    !

      do i = 1, nbox
        dx = d(b(i)%level) * 0.5
        r = Sqrt(3.d0) * dx           ! max distance

        do j = 0, b(i)%n_p - 1
          k = ind(b(i)%p_start + j)
          dx = coord(1, k) - b(i)%x0
          dy = coord(2, k) - b(i)%y0
          dz = coord(3, k) - b(i)%z0
          t = Sqrt (dx**2 + dy**2 + dz**2)

          if (t >= r) then
            write (iw, *) " *2* Error", i, j
          end if
        end do
      end do
    end if

    ierr = 0

    call store_tesselation (handle, n)

  end subroutine divide_box

  function count_short_ints ( coord, nc, simulate_dir_int, compute ) result(n)
    !
    !.. Implicit Declarations ..
    implicit none

    !
    !.. Formal Arguments ..
    integer, intent(in) :: nc
    double precision, intent (in), dimension(nc, *) :: coord
    logical, intent( in ) :: compute

    interface
       subroutine simulate_dir_int (ind1, n1, ind2, n2, c, nc, same, &
                & compute, n )
         integer, intent(in) :: n1, n2, nc
         integer, intent(in), dimension(*) :: ind1, ind2
         double precision, intent(in), dimension(nc, *) :: c
         logical, intent(in) :: same, compute
         integer, intent(inout) :: n
       end subroutine simulate_dir_int
    end interface
    !
    !.. Local Scalars ..
    integer :: n, i, j, j1
    !
    ! ... Executable Statements ...

    n = 0

    do i = 1, nbox
      if (i == 34) then
        j1 = 0
      end if
      if (b(i)%n_childs == 0) then
        call simulate_dir_int (points_index(b(i)%p_start), b(i)%n_p, &
             & points_index(b(i)%p_start), b(i)%n_p, &
             & coord, nc, .true., compute, n)

        do j1 = 1, b(i)%n_l1
          j = b(i)%l1(j1)
          call simulate_dir_int (&
               & points_index(b(i)%p_start), b(i)%n_p, &
               & points_index(b(j)%p_start), b(j)%n_p, &
               & coord, nc, .false., compute, n)

        end do
      end if
    end do

  end function count_short_ints

  subroutine afmm (coord, nc, n, q, res, res_size, dir_int, mult_int, &
       & tricky_multipoles )

    !
    !.. Implicit Declarations ..
    implicit none
    !
    !.. Formal Arguments ..
    integer, intent(in) :: nc, n, res_size
    double precision, intent (in), dimension(n) :: q
    double precision, intent (in), dimension(nc, n) :: coord
    double precision, intent (out), dimension(res_size) :: res

    optional :: tricky_multipoles

    interface
       subroutine dir_int (ind1, n1, ind2, n2, c, nc, q, r, same)
         integer, intent(in) :: n1, n2, nc
         integer, intent(in), dimension(*) :: ind1, ind2
         double precision, intent(in), dimension(nc, *) :: c
         double precision, intent(in), dimension(*) :: q
         double precision, intent(inout), dimension(*) :: r
         logical, intent(in) :: same
       end subroutine dir_int

       subroutine mult_int (ind, n, psi, p, x0, y0, z0, y_norm, pmn, &
            & c, n3, q, r)

         integer, intent(in) :: n, p, n3
         integer, intent(in), dimension(*) :: ind
         complex (kind=Kind ((1.d0, 1.d0))), dimension(-p:p, 0:p), &
            & intent(in) :: psi
         double precision, dimension(-p:p, 0:p), intent(in) :: y_norm
         double precision, dimension(-p:p, 0:p), intent(out) :: pmn
         double precision, intent(in), dimension(n3, *) :: c
         double precision, intent(in), dimension(*) :: q
         double precision, intent(inout), dimension(*) :: r
         double precision, intent(in) :: x0, y0, z0
       end subroutine mult_int

       subroutine tricky_multipoles (&
            & x0, y0, z0, ind, n, pmn, y_norm, c, nc, q, psi, p)

         implicit none

         double precision, intent(in) :: x0, y0, z0
         integer, intent(in) :: n, nc, p
         integer, dimension(*), intent(in) :: ind
         double precision, dimension(nc, *), intent(in) :: c
         double precision, dimension(*), intent(in) :: q
         double precision, dimension(-p:p, 0:p), intent(in) :: y_norm
         double precision, dimension(-p:p, 0:p), intent(inout) :: pmn

         complex (kind=Kind ((1.d0, 1.d0))), dimension(-p:p, 0:p), &
            & intent(inout) :: psi

       end subroutine tricky_multipoles

    end interface
    double precision :: t, tt
    integer :: i, j, j1
    intrinsic Present
    if (Present (tricky_multipoles)) then
      call afmm_step1 (coord, nc, n, q, tricky_multipoles)
    else
      call afmm_step1 (coord, nc, n, q)
    end if

    call afmm_step2
    call afmm_step35
    res = 0

    !
    !       Step 7. Calculate the potential at each charge in a box from
    !       ======  direct interactions
    !
    t = 0
    tt = n
    tt = tt*(tt-1.d0) * 0.5d0
    do i = 1, nbox
      if (b(i)%n_childs == 0) then

        call dir_int (points_index(b(i)%p_start), b(i)%n_p, &
             & points_index(b(i)%p_start), b(i)%n_p, &
             & coord, nc, q, res, .true.)

        t = t + (b(i)%n_p * (b(i)%n_p - 1)) /2

        do j1 = 1, b(i)%n_l1
          j = b(i)%l1(j1)
          call dir_int (&
               & points_index(b(i)%p_start), b(i)%n_p, &
               & points_index(b(j)%p_start), b(j)%n_p, &
               & coord, nc, q, res, .false.)

          t = t + b(i)%n_p * b(j)%n_p
        end do
      end if

    end do
    !
    !       Step 6. Calculate the potential at each charge in a box from
    !       ======  it's local expansion
    !
    do i = 1, nbox
      if (b(i)%n_childs == 0) then
        call mult_int (points_index(b(i)%p_start), b(i)%n_p, &
             & b(i)%psi, p, b(i)%x0, b(i)%y0, b(i)%z0, &
             & y_norm, pmn, coord, nc, q, res)
      end if
    end do
  end subroutine afmm

  subroutine afmm_step1 (coord, nc, n, q, tricky_multipoles)
    !
    !.. Implicit Declarations ..
    implicit none
    !
    !.. Formal Arguments ..
    integer, intent(in) :: nc, n
    double precision, intent (in), dimension(n) :: q
    double precision, intent (in), dimension(nc, n) :: coord

    optional :: tricky_multipoles

    interface

       subroutine tricky_multipoles (x0, y0, z0, ind, n, pmn, y_norm, &
            & c, nc, q, psi, p)

         implicit none

         double precision, intent(in) :: x0, y0, z0
         integer, intent(in) :: n, nc, p
         integer, dimension(*), intent(in) :: ind
         double precision, dimension(nc, *), intent(in) :: c
         double precision, dimension(*), intent(in) :: q
         double precision, dimension(-p:p, 0:p), intent(in) :: y_norm
         double precision, dimension(-p:p, 0:p), intent(inout) :: pmn

         complex (kind=Kind ((1.d0, 1.d0))), dimension(-p:p, 0:p), &
              & intent(inout) :: psi

       end subroutine tricky_multipoles

    end interface
    integer :: i
    intrinsic Present
    !
    ! ... Executable Statements ...

    do i = 1, nbox
      b(i)%phi = (0.d0, 0.d0)
      b(i)%psi = (0.d0, 0.d0)
    end do
    !
    !  Step 1. For each childless box form a multipole expansion about it's
    !  ======  centre from all it's charges
    !
    do i = 1, nbox
      if (b(i)%n_childs == 0) then
        if (Present (tricky_multipoles)) then
          call tricky_multipoles (&
               & b(i)%x0, b(i)%y0, b(i)%z0, points_index(b(i)%p_start), &
               & b(i)%n_p, pmn, y_norm, coord, nc, q, b(i)%phi, p)
        else
          call normal_multipoles (&
               & b(i)%x0, b(i)%y0, b(i)%z0, points_index(b(i)%p_start), &
               & b(i)%n_p, pmn, y_norm, coord, nc, q, b(i)%phi, p)
        end if
      end if
    end do
  end subroutine afmm_step1

  subroutine normal_multipoles (x0, y0, z0, ind, n, pmn, y_norm, &
       & c, nc, q, psi, p)
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
    double precision :: dx, dy, dz, t, tt, r, cos_t, phi
    complex (kind=dcmplx_kind) :: tc, tc1
    integer :: j, k, n1, m1
    !
    !.. Intrinsic Functions ..
    intrinsic Sqrt, Atan2, Exp, Cmplx, Conjg
    !
    ! ... Executable Statements ...

    do j = 1, n
      k = ind(j)
      dx = c(1, k) - x0
      dy = c(2, k) - y0
      dz = c(3, k) - z0
      r = Sqrt (dx**2 + dy**2 + dz**2)
      cos_t = dz / r
      phi = Atan2 (dy, dx)

      call get_legendre (p, cos_t, pmn)

      psi(0, 0) = psi(0, 0) + Cmplx (q(k), 0.d0, dcmplx_kind)

      t = 1.d0
      do n1 = 1, p
        t = t * r
        do m1 = 0, n1
          tc = Exp (Cmplx (0.d0, -m1*phi, dcmplx_kind))
          tt = y_norm(-m1, n1) * pmn(-m1, n1) * t * q(k)
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

  end subroutine normal_multipoles

  subroutine afmm_step2
    implicit none
    double precision :: dx, dy, dz, tt, r, cos_t, phi
    complex (kind=dcmplx_kind) :: tc, tc1
    integer :: i, j, j1, k, m, n1, m1, n2, m2
    intrinsic Sqrt, Exp, Cmplx, Atan2, Abs, Conjg

    !
    !       Step 2. For each parent box form a multipole expansion
    !       ======  around it's centre by merging multipole expansion
    !               from it's children
    !
    do i = nlev-1, 2, -1
      do j = level(i), level(i+1) - 1
        do k = 1, b(j)%n_childs
          j1 = b(j)%childs(k)

          dx = b(j1)%x0 - b(j)%x0
          dy = b(j1)%y0 - b(j)%y0
          dz = b(j1)%z0 - b(j)%z0

          r = Sqrt (dx**2 + dy**2 + dz**2)
          cos_t = dz / r
          phi = Atan2 (dy, dx)

          call get_legendre (p, cos_t, pmn)

          do n1 = 0, p
           do m1 = 0, n1
              tc = Cmplx (0.d0, 0.d0, dcmplx_kind)
              do n2 = 0, n1
                do m2 = -n2, n2

                  if (Abs (m1-m2) > n1-n2) then
                    cycle
                  end if

                  m = Abs (m1) - Abs (m2) - Abs (m1-m2)

                  tt = amn(m2, n2) * amn(m1-m2, n1-n2) &
                       & / amn(m1, n1) * r**n2 &
                       & * y_norm(-m2, n2) * pmn(-m2, n2)

                  tc1 = Exp (Cmplx (0, -m2*phi, dcmplx_kind)) &
                       & * b(j1)%phi(m1-m2, n1-n2) &
                       & * (Cmplx (0, 1, dcmplx_kind)) **m

                  tc = tc + tt * tc1

                end do
              end do

              b(j)%phi(m1, n1) = b(j)%phi(m1, n1) + tc
            end do
          end do

        end do

        do n1 = 1, p
          do m1 = 1, n1
            b(j)%phi(-m1, n1) = Conjg (b(j)%phi(m1, n1))
          end do
        end do
      end do
    end do
  end subroutine afmm_step2

  subroutine afmm_step35
    implicit none
    double precision :: dx, dy, dz, tt, r, cos_t, phi
    complex (kind=dcmplx_kind) :: tc, tc1
    integer :: i, j, j1, m, n1, m1, n2, m2
    intrinsic Sqrt, Cmplx, Conjg, Exp, Atan2, Abs
    !
    ! ... Executable Statements ...

    !
    !       Steps 3 and 8. Interaction of not very well separated points
    !       =============
    !
    !       These step is now combined into step 7 due to different
    !       definition of the list 1.
    !

    !
    !       Step 4. Move multipole expansion to boxes in l2 list.
    !       ======
    !

    do i = level(2), nbox
      do j1 = 1, b(i)%n_l2
        j = b(i)%l2(j1)

        dx = b(i)%x0 - b(j)%x0
        dy = b(i)%y0 - b(j)%y0
        dz = b(i)%z0 - b(j)%z0

        r = Sqrt (dx**2 + dy**2 + dz**2)
        cos_t = dz / r
        if( Abs( dx ) < 1.d-6 .and. Abs( dy ) < 1.d-6 ) then
           phi = 0.d0
        else
           phi = Atan2 (dy, dx)
        end if


        call get_legendre (p, cos_t, pmn)

        do n1 = 0, p
          do m1 = 0, n1

            tc = (0.d0, 0.d0)

            do n2 = 0, p

              if (n1+n2 > p) then
                cycle
              end if

              do m2 = -n2, n2

                if (Abs (m2-m1) > n1+n2) then
                  cycle
                end if

                m = Abs (m1-m2) - Abs (m1) - Abs (m2)

                tt = amn(m2, n2) * amn(m1, n1) &
                     & / (amn(m2-m1, n1+n2) * r**(n1+n2+1)) &
                     & * y_norm(m2-m1, n1+n2) * pmn(m2-m1, n1+n2) / (-1)**n2

                tc1 = Exp (Cmplx (0, (m2-m1)*phi, dcmplx_kind)) &
                     & * b(i)%phi(m2, n2) &
                     & * (Cmplx (0, 1, dcmplx_kind)) **m

                tc = tc + tt * tc1
              end do
            end do

            b(j)%psi(m1, n1) = b(j)%psi(m1, n1) + tc
          end do
        end do

        do n1 = 1, p
          do m1 = 1, n1
            b(j)%psi(-m1, n1) = Conjg (b(j)%psi(m1, n1))
          end do
        end do

      end do
    end do
    !
    !       Step 5. For each parent shift the centre of it's local expansion
    !       ======  to it's children.
    !

    do i = level(2), nbox
      do j1 = 1, b(i)%n_childs
        j = b(i)%childs(j1)

        dx = b(i)%x0 - b(j)%x0
        dy = b(i)%y0 - b(j)%y0
        dz = b(i)%z0 - b(j)%z0

        r = Sqrt (dx**2 + dy**2 + dz**2)
        cos_t = dz / r
        phi = Atan2 (dy, dx)

        call get_legendre (p, cos_t, pmn)

        do n1 = 0, p
          do m1 = -n1, n1
            tc = (0.d0, 0.d0)
            do n2 = n1, p
              do m2 = -n2, n2

                if (Abs (m2-m1) > n2-n1) then
                  cycle
                end if

                m = Abs (m2) - Abs (m2-m1) - Abs (m1)

                tt = amn(m2-m1, n2-n1) * amn(m1, n1) &
                     & / amn(m2, n2) * r**(n2-n1) &
                     & * y_norm(m2-m1, n2-n1) * pmn(m2-m1, n2-n1) &
                     & / (-1)**(n2+n1)

                tc1 = Exp (Cmplx (0, (m2-m1)*phi, dcmplx_kind)) &
                     & * b(i)%psi(m2, n2) * (Cmplx (0, 1, dcmplx_kind))**m

                tc = tc + tt * tc1
              end do
            end do

            b(j)%psi(m1, n1) = b(j)%psi(m1, n1) + tc
          end do
        end do
      end do
    end do
  end subroutine afmm_step35

  subroutine simple_mm (coord, nc, n, q, res, dir_int, far_int, &
       & tricky_multipoles)

    !
    !.. Implicit Declarations ..
    implicit none
    !
    !.. Formal Arguments ..
    integer, intent(in) :: nc, n
    double precision, intent (in), dimension(n) :: q
    double precision, intent (in), dimension(nc, n) :: coord
    double precision, intent (out), dimension(n) :: res

    optional :: tricky_multipoles

    interface
       subroutine dir_int (ind1, n1, ind2, n2, c, nc, q, r, same)
         integer, intent(in) :: n1, n2, nc
         integer, intent(in), dimension(*) :: ind1, ind2
         double precision, intent(in), dimension(nc, *) :: c
         double precision, intent(in), dimension(*) :: q
         double precision, intent(inout), dimension(*) :: r
         logical, intent(in) :: same
       end subroutine dir_int

       subroutine far_int (ind, n, phi, p, x0, y0, z0, &
            & y_norm, pmn, c, nc, q, r)

         integer, intent(in) :: n, p, nc
         integer, intent(in), dimension(*) :: ind
         complex (kind=Kind ((1.d0, 1.d0))), dimension(-p:p, 0:p), &
              &intent(in) :: phi
         double precision, dimension(-p:p, 0:p), intent(in) :: y_norm
         double precision, dimension(-p:p, 0:p), intent(out) :: pmn
         double precision, intent(in), dimension(nc, *) :: c
         double precision, intent(in), dimension(*) :: q
         double precision, intent(inout), dimension(*) :: r
         double precision, intent(in) :: x0, y0, z0
       end subroutine far_int

       subroutine tricky_multipoles (&
            & x0, y0, z0, ind, n, pmn, y_norm, c, nc, q, psi, p)

         implicit none

         double precision, intent(in) :: x0, y0, z0
         integer, intent(in) :: n, nc, p
         integer, dimension(*), intent(in) :: ind
         double precision, dimension(nc, *), intent(in) :: c
         double precision, dimension(*), intent(in) :: q
         double precision, dimension(-p:p, 0:p), intent(in) :: y_norm
         double precision, dimension(-p:p, 0:p), intent(inout) :: pmn

         complex (kind=Kind ((1.d0, 1.d0))), dimension(-p:p, 0:p), &
              & intent(inout) :: psi

       end subroutine tricky_multipoles

    end interface
    double precision :: t, tt
    integer :: i, j, j1
    intrinsic Present
    !
    ! ... Executable Statements ...

    res = 0

    if (Present (tricky_multipoles)) then
      call afmm_step1 (coord, nc, n, q, tricky_multipoles)
    else
      call afmm_step1 (coord, nc, n, q)
    end if

    call afmm_step2

    !
    !       Step 3. Compute far interactions
    !
    do i = level(2), nbox
      do j1 = 1, b(i)%n_l2
        j = b(i)%l2(j1)
        call far_int (points_index(b(i)%p_start), b(i)%n_p, &
             & b(j)%phi, p, b(j)%x0, b(j)%y0, b(j)%z0, &
             & y_norm, pmn, coord, nc, q, res)
      end do
    end do


    !
    !       Step 4. Calculate the potential at each charge in a box from
    !       ======  direct interactions
    !


    t = 0
    tt = n
    tt = tt*(tt-1.d0) * 0.5d0
    do i = 1, nbox
      if (b(i)%n_childs == 0) then
        call dir_int (points_index(b(i)%p_start), b(i)%n_p, &
             & points_index(b(i)%p_start), b(i)%n_p, &
             & coord, nc, q, res, .true.)

        t = t + (b(i)%n_p * (b(i)%n_p - 1)) /2

        do j1 = 1, b(i)%n_l1
          j = b(i)%l1(j1)
          call dir_int (&
               & points_index(b(i)%p_start), b(i)%n_p, &
               & points_index(b(j)%p_start), b(j)%n_p, &
               & coord, nc, q, res, .false.)

          t = t + b(i)%n_p * b(j)%n_p
        end do
      end if
    end do
  end subroutine simple_mm

  subroutine split_box (parent, level, x0, y0, z0, d, coord, nc, n, &
       & b_pos, ind_pos, m)
    !
    !       Splits box to 8 smaller boxes & ignore empty boxes
    !

    !
    !.. Implicit Declarations ..
    implicit none
    !
    !.. Formal Arguments ..
    integer, intent(in) :: n, nc, parent, level, b_pos, ind_pos
    double precision, intent(in) :: x0, y0, z0, d
    double precision, intent(in), dimension(nc, *) :: coord
    integer, intent(out) :: m
    !
    !.. Local Scalars ..
    integer :: i, n1, n2, n3, n4, n5, n6
    !
    ! ... Executable Statements ...

    m = 0

    call split_index (1, coord, nc, x0, ind_pos, n, n1, n2)

    if (n1 > 0) then
      call split_index (2, coord, nc, y0, ind_pos, n1, n3, n4)

      if (n3 > 0) then
        call split_index (3, coord, nc, z0, ind_pos, n3, n5, n6)

        if (n5 > 0) then
          m = m + 1
          b(b_pos + m)%x0 = x0 - d
          b(b_pos + m)%y0 = y0 - d
          b(b_pos + m)%z0 = z0 - d
          b(b_pos + m)%p_start = 1
          b(b_pos + m)%n_p = n5
        end if

        if (n6 > 0) then
          m = m + 1
          b(b_pos + m)%x0 = x0 - d
          b(b_pos + m)%y0 = y0 - d
          b(b_pos + m)%z0 = z0 + d
          b(b_pos + m)%p_start = 1 + n5
          b(b_pos + m)%n_p = n6
        end if
      end if

      if (n4 > 0) then
        call split_index (3, coord, nc, z0, ind_pos+n3, n4, n5, n6)

        if (n5 > 0) then
          m = m + 1
          b(b_pos + m)%x0 = x0 - d
          b(b_pos + m)%y0 = y0 + d
          b(b_pos + m)%z0 = z0 - d
          b(b_pos + m)%p_start = n3 + 1
          b(b_pos + m)%n_p = n5
        end if

        if (n6 > 0) then
          m = m + 1
          b(b_pos + m)%x0 = x0 - d
          b(b_pos + m)%y0 = y0 + d
          b(b_pos + m)%z0 = z0 + d
          b(b_pos + m)%p_start = n3 + 1 + n5
          b(b_pos + m)%n_p = n6
        end if
      end if
    end if

    if (n2 > 0) then
      call split_index (2, coord, nc, y0, ind_pos+n1, n2, n3, n4)

      if (n3 > 0) then
        call split_index (3, coord, nc, z0, ind_pos+n1, n3, n5, n6)

        if (n5 > 0) then
          m = m + 1
          b(b_pos + m)%x0 = x0 + d
          b(b_pos + m)%y0 = y0 - d
          b(b_pos + m)%z0 = z0 - d
          b(b_pos + m)%p_start = n1 + 1
          b(b_pos + m)%n_p = n5
        end if

        if (n6 > 0) then
          m = m + 1
          b(b_pos + m)%x0 = x0 + d
          b(b_pos + m)%y0 = y0 - d
          b(b_pos + m)%z0 = z0 + d
          b(b_pos + m)%p_start = n1 + 1 + n5
          b(b_pos + m)%n_p = n6
        end if
      end if

      if (n4 > 0) then
        call split_index (3, coord, nc, z0, ind_pos+n1+n3, n4, n5, n6)

        if (n5 > 0) then
          m = m + 1
          b(b_pos + m)%x0 = x0 + d
          b(b_pos + m)%y0 = y0 + d
          b(b_pos + m)%z0 = z0 - d
          b(b_pos + m)%p_start = n1 + n3 + 1
          b(b_pos + m)%n_p = n5
        end if

        if (n6 > 0) then
          m = m + 1
          b(b_pos + m)%x0 = x0 + d
          b(b_pos + m)%y0 = y0 + d
          b(b_pos + m)%z0 = z0 + d
          b(b_pos + m)%p_start = n1 + n3 + 1 + n5
          b(b_pos + m)%n_p = n6
        end if
      end if
    end if

    do i = b_pos + 1, b_pos + m
      b(i)%parent = parent
      b(i)%level = level
      b(i)%n_childs = 0
      b(i)%n_l1 = 0
      b(i)%n_l2 = 0
      b(i)%phi = (0.d0, 0.d0)
      b(i)%psi = (0.d0, 0.d0)
    end do

  end subroutine split_box

  subroutine split_index (k, coord, nc, g, ind_pos, n, n1, n2)
    !
    !.. Implicit Declarations ..
    implicit none
    !
    !.. Formal Arguments ..
    integer, intent(in) :: k, n, nc, ind_pos
    double precision, intent(in), dimension(nc, *) :: coord
    double precision, intent(in) :: g
    integer, intent (out) :: n1, n2

    !
    !.. Local Scalars ..
    integer :: i
    !
    ! ... Executable Statements ...

    n1 = 0
    n2 = n + 1

    do i = 1, n
      if (coord(k, ind(ind_pos+i)) < g) then
        n1 = n1 + 1
        ind1(n1) = ind(ind_pos+i)
      else
        n2 = n2 - 1
        ind1(n2) = ind(ind_pos+i)
      end if
    end do

    n2 = n - n1

    do i = 1, n
      ind(ind_pos+i) = ind1(i)
    end do

  end subroutine split_index

  subroutine get_legendre (n_max, x, pmn)
    !
    !       Computes array of Legendre functions
    !

    !
    !.. Implicit Declarations ..
    implicit none
    !
    !.. Formal Arguments ..
    integer, intent (in) :: n_max
    double precision, intent (inout) :: x
    double precision, intent (out), dimension(-n_max:n_max, 0:n_max) :: pmn

    !
    !.. Local Scalars ..
    double precision :: x2, t, s, s1
    integer :: n, m
    !
    !.. Intrinsic Functions ..
    intrinsic Sqrt, Dble
    !
    ! ... Executable Statements ...

    x2 = x**2
    if (x2 >= 1.d0) then
      if (x > 0) then
        x = 1.d0
      else
        x = -1.d0
      end if
      x2 = 1.d0
      t = 0.d0
    else
      t = Sqrt (1.d0 - x2)
    end if

    pmn(0, 0) = 1.d0
    pmn(0, 1) = x
    pmn(1, 1) = -t
    pmn(0, 2) = 0.5d0 * (3.d0*x2 - 1.d0)
    pmn(1, 2) = -3.d0*x*t
    pmn(2, 2) = 3.d0*t*t

    s = pmn(2, 2)
    do n = 3, n_max
      s = s * (2*n - 1) * t
      pmn(n, n) = s * (-1)**n
    end do

    do n = 3, n_max
      do m = n-1, 0, -1
        s = 1.d0 / (Dble(n-m))
        s1 = x*(2*n-1)*pmn(m, n-1)
        if (m <= n-2) then
          s1 = s1 - (n+m-1)*pmn(m, n-2)
        end if
        pmn(m, n) = s1 * s
      end do
    end do

    do n = 1, n_max
      do m = 1, n
        pmn(-m, n) = pmn(m, n)
      end do
    end do

  end subroutine get_legendre

  subroutine afmm_ini
    !
    !.. Implicit Declarations ..
    implicit none

    !
    !.. Local Scalars ..
    integer :: i, j
    double precision :: t
    !
    !.. Intrinsic Functions ..
    intrinsic Sqrt
    !
    ! ... Executable Statements ...

    ! --- factorials

    t = 1.d0
    fact(0) = t
    fact(1) = t

    do i = 2, 2*p
      t = t * i
      fact(i) = t
    end do

    ! --- Spherical harmonics norm

    y_norm(0, 0) = 1.d0
    do i = 1, p
      do j = 0, i
        t = Sqrt (fact(i-j) / fact(i+j))
        y_norm(j, i) = t
        y_norm(-j, i) = t
      end do
    end do

    ! --- amn

    t = 1.d0
    amn(0, 0) = 1.d0
    do i = 1, p
      t = -t
      do j = -i, i
        amn(j, i) = t / Sqrt (fact(i-j) * fact(i+j))
      end do
    end do

  end subroutine afmm_ini

end module afmm_C
