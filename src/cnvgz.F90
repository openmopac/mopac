subroutine cnvgz (pnew, p, p1, p2, p3, niter, idiag)
    use molkst_C, only: norbs, mpack
    use MOZYME_C, only : use_three_point_extrap, pmax
    implicit none
    integer, intent (in) :: niter
    integer, dimension (norbs), intent (in) :: idiag ! Pointers to diagonal elements
    double precision, dimension (norbs), intent (inout) :: p1, p2, p3
    double precision, dimension (mpack), intent (inout) :: p, pnew
    integer :: i, j
    double precision :: damp, fac, faca, facb, sa
    intrinsic Abs, Max, Min, Mod, Sign, Sqrt
   !
   ! Save the diagonal of the current and previous density for later use
   !
    do i = 1, norbs
      j = idiag(i)
      p3(i) = pnew(j)
      p2(i) = p(j)
    end do
!
!   Calculate the maximum and RMS change in the density matrix
!
    pmax = 0.0d0
    do i = 1, mpack
      sa = Abs (pnew(i) - p(i))
      pmax = Max (pmax, sa)
    end do
    if (use_three_point_extrap) then
!
!   Three-point extrapolation in use. Extrapolation is used on
!   every third SCF cycle
!
      if (Mod (niter, 3) == 0) then
        faca = 0.0d0
        facb = 0.0d0
        do i = 1, norbs
          sa = Abs (p3(i) - p2(i))
          faca = faca + sa ** 2
          facb = facb + (p3(i)-2.0d0*p2(i)+p1(i)) ** 2
        end do
!
        if (facb > 0.d0 .and. faca < (100.d0*facb)) then
          fac = Sqrt (faca/facb)
          pnew = pnew + fac * (pnew-p)
        end if
!
      end if
!
!   From iteration 4 on, the change in the diagonal elements is 
!   limited to 'damp'
!
      if (niter > 3) then
        damp = 0.05d0
        if (pmax > damp) then
          do i = 1, norbs
            j = idiag(i)
            if (Abs (p3(i)-p2(i)) > damp) then
              pnew(j) = p2(i) + Sign (damp, p3(i)-p2(i))
              pnew(j) = Min (2.0d0, Max (pnew(j), 0.0d0))
            end if
          end do
        end if
      end if
    end if ! three-point extrapolation
!
!  Save the density for use in the next iteration
!
    p1 = p2
    p = pnew
end subroutine cnvgz
