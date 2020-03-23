double precision function pardip (coord, nat)
  !**********************************************************************
  !
  !   PARDIP  1:  Computes the dipole moment.
  !
  !**********************************************************************
    use molkst_C, only : numat
    use parameters_c, only : tore
    use common_arrays_C, only : p, q
  !
  !.. Implicit Declarations ..
    implicit none
  !
  !.. Formal Arguments ..
    double precision, dimension (3, numat), intent (inout) :: coord
    integer, dimension (numat), intent (in) :: nat
  !
  !.. Local Scalars ..
    integer :: i, l
  !
  !.. Local Arrays ..
    double precision, dimension (3) :: dumy
  !
  !.. External Calls ..
    external chrge
  !
  !.. External Functions ..
    double precision, external :: dipole
  !
  ! ... Executable Statements ...
  !
  !
  !  Calculate the total number of electrons on each atom
  !
    call chrge (p, q)
  !
  !  Calculate the charge on each atom
  !
    do i = 1, numat
      l = nat(i)
      q(i) = tore(l) - q(i)
    end do
  !
  !  Calculate the dipole
  !
    pardip = dipole (p, coord, dumy, 1)
end function pardip
