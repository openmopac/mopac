subroutine eimp ()
!
!   Store in the s-s location of array p the value of the Fock terms relating
!   atom A and atom B, for all pairs of atoms.  This quantity will be used by 
!   the diagonalizer
!
  use molkst_C, only: numat
  use common_arrays_C, only : p, f
  use MOZYME_C, only : iorbs
  implicit none
  integer :: i, j, k, l, m
  double precision :: sum
  integer, external :: ijbo
  do i = 1, numat
    do j = 1, i - 1
      k = ijbo (i, j)
      if (k >= 0) then
        l = iorbs(i) * iorbs(j)
        if (l /= 0) then
          l = k + l
          sum = 0.d0
          do m = k + 1, l
            sum = sum + f(m) ** 2
          end do
          p(k+1) = sum
        end if
      end if
    end do
  end do
end subroutine eimp
