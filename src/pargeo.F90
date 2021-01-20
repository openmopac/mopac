double precision function pargeo (grad, wtgeo, refgeo, loc, diffs, ndif)
  !***********************************************************************
  !
  !   PARGEO  1:  Computes the gradient error.
  !
  !***********************************************************************
    use molkst_C, only : nvar
    use param_global_C, only : power, wtz
    implicit none
    double precision, dimension (300), intent (in) :: grad
    double precision, intent (in) :: wtgeo
    character (len=12), dimension (300), intent (in) :: refgeo
    integer, dimension (2, 300), intent (in) :: loc
    double precision, dimension (300), intent (inout) :: diffs
    integer, intent (out) :: ndif
    integer :: i, lim
    double precision :: relative_weight, sum
    intrinsic Min
!------------------------------------------------------------------------
    lim = Min (100, nvar)
    ndif = 0
    sum = 0.d0
    do i = 1, lim
      if (refgeo(i) /= " ") then
        ndif = ndif + 1
        if (loc(2, i) == 1) then
                                      !
          relative_weight = 1.d0      !  Relative weight of bond-length
                                      !
        else if (loc(2, i) == 2) then !
                                      !
          relative_weight = 1.d0      !  Relative weight of bond-angle
                                      !
        else if (loc(2,i) == 3) then  !
                                      !
          relative_weight = 1.d0      !  Relative weight of dihedral
                                      !                                                                                                                                                                                                                                                                                                                                                                                                                                                             
        end if
        diffs(ndif) = grad(i) * wtgeo * relative_weight * wtz
        sum = sum + Abs(diffs(ndif)) ** power
      end if
    end do
    pargeo = sum
end function pargeo
