subroutine rapid0 (loop)
!
!
    use param_global_C, only : valvar, numvar, xparamp, locvar
!
!
    implicit none
    integer, intent (in) :: loop
!----------------------------------------------------------------
    integer :: i
    double precision :: funct1
!----------------------------------------------------------------
  !
  !   DELTAS will hold the CHANGE IN VALUE of the parameters
  !
    do i = 1, numvar
      xparamp(i) = 0.d0
    end do
  !
  !  Optimize the parameters
  !
    call rapid1 (loop, xparamp, numvar, funct1)
  !
  !  Update the values of the parameters
  !
    do i = 1, numvar
      valvar(i) = valvar(i) - xparamp(i)
      if(locvar(1,i) > 3 .and. locvar(1,i) < 7) valvar(i) = max(0.05d0, valvar(i))
      call update (locvar(1, i), locvar(2, i), valvar(i), 0.d0)
    end do
end subroutine rapid0
