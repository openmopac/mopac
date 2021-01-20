subroutine depfn (fns, diffns, loop, loch, ndep, ihref, nfns, numvar)
    use param_global_C, only : weight, hofcal
    implicit none
    integer, intent (in) :: loop, ndep, nfns, numvar
    integer, dimension (loop), intent (in) :: loch
    integer, dimension (ndep), intent (in) :: ihref
    double precision, dimension (nfns), intent (inout) :: fns
    double precision, dimension (numvar, nfns), intent (inout) :: diffns
    integer :: i, iref, j, k, l
  !
  ! ... Executable Statements ...
  !
  !
  !   Modify FNS and DIFFNS to allow for dependent molecules
  !
  !
  !    LOOP = Current molecule number (The H.o.F. of this molecule
  !           depends on the H.o.F. os the molecules in IHREF)
  !
  !    IREF = Address of current molecule in FNS and in DIFFNS
  !
  !    J    = Index of the dependent molecule.
  !
  !    K    = Address of dependent molecule in FNS and in DIFFNS.
    iref = loch(loop)
    do i = 1, ndep
      j = ihref(i)
      k = loch(j)
    !
    !  Calculate the new error in the H.o.F in one step.
    !
      fns(iref) = fns(iref) - hofcal(j) * weight(1, loop)
    !
    !  Calculate the new derivatives:
    !
      do l = 1, numvar
        diffns(l, iref) = diffns(l, iref) - diffns(l, k) / weight(1, j) * &
       & weight(1, loop)
      end do
    end do
end subroutine depfn
