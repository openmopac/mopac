subroutine chrge_for_MOZYME (p, q)
!
! Evaluate the atomic partial charge.
!
    use molkst_C, only: mpack, numat
    use MOZYME_C, only : iorbs
    implicit none
    double precision, dimension (mpack), intent (in) :: p
    double precision, dimension (numat), intent (out) :: q
!
    integer :: i, ii, j
    double precision :: sum
    integer, external :: ijbo
!
    do i = 1, numat
      ii = ijbo (i, i)
      sum = 0.d0
      do j = 1, iorbs(i)
        ii = ii + j
        sum = sum + p(ii)
      end do
      q(i) = sum
    end do
end subroutine chrge_for_MOZYME
