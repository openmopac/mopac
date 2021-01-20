subroutine ionout (ions, m, num)
    use molkst_C, only: natoms, numat, maxtxt
    use chanel_C, only: iw
    use elemts_C, only: elemnt
    use common_arrays_C, only : na, nb, nc, txtatm, labels
    implicit none
    integer, intent (in) :: m, num
    integer, dimension (4, num), intent (in) :: ions
!
    integer, dimension (:), allocatable :: nlabels
    character (len=40) :: blank
    logical :: Int
    integer :: i, j, klim, nabc11, il, iu
    intrinsic Min
    allocate(nlabels(natoms))
    j = 0
    do i = 1, natoms
      if (labels(i) /= 99) j = j + 1
      nlabels(i) = j
    end do
    j = Min (3, numat)
    do i = 1, j
      if (na(i) /= 0) then
        Int = .true.
        klim = 3
        go to 1000
      end if
    end do
    Int = .false.
    klim = 0
1000 blank = " "
    i = maxtxt + 5
    if (i /= 5) then
      i = i + 5
    end if
    write (iw,*)
    if (Int) then
      write (iw, "('     ION   ATOM   WITH   TYPE',A,' CONNECTIVITY')") blank (:i)
      write (iw, "('            No.  DUMMYS')")
    else
      write (iw, "('     ION   ATOM   WITH   TYPE')") 
      write (iw, "('            No.  DUMMYS')")
    end if
   
    nabc11 = na(1)
    na(1) = 0
    if (m /= 3) then
      if (maxtxt /= 0) then
        do i = 1, num
          j = ions(m, i)
          if (txtatm(j)(1:1) == "(") then
          il = 2
          iu = maxtxt - 1
        else
          il = 1
          iu = maxtxt
        end if
          write (iw, "(3I7,3X,A2,'(',A,')',3I7)") i, nlabels(j), j, &
         & elemnt (labels(j)), txtatm(j) (il:iu), na(j), nb(j), nc(j)
        end do
      else
        do i = 1, num
          j = ions(m, i)
          write (iw, "(3I7,3X,A2,3I7)") i,  nlabels(j), j, &
         & elemnt (labels(j)), na(j), nb(j), nc(j)
        end do
      end if
    else if (maxtxt /= 0) then
      do i = 1, num
        j = ions(m, i)
        if (txtatm(j)(1:1) == "(") then
          il = 2
          iu = maxtxt - 1
        else
          il = 1
          iu = maxtxt
        end if
        if (klim == 0) then
          write (iw, "(3I7,3X,A2,'(',A,')',A,SP,I2)") i, nlabels(j), j, elemnt &
         & (labels(j)), txtatm (j) (il:iu), "  CHARGE: ", ions (4, i)
        else
          write (iw, "(3I7,3X,A2,'(',A,')',A,SP,I2,3I7)") i, nlabels(j), j, elemnt &
         & (labels(j)), txtatm (j) (il:iu), "  CHARGE: ", ions (4, i), &
         &  na(j), nb(j), nc(j)
        end if
      end do
    else
      do i = 1, num
        j = ions(m, i)
        if (klim == 0) then
          write (iw, "(3I7,3X,A2,A,SP,I2)") i, nlabels(j), j, elemnt (labels(j)), &
         & "  CHARGE: ", ions (4, i)
        else
          write (iw, "(3I7,3X,A2,3I7,A,SP,I2)") i, nlabels(j), j, elemnt (labels(j)), &
         &  na(j), nb(j), nc(j), "  CHARGE: ", ions (4, i)
        end if
      end do
    end if
    na(1) = nabc11
    deallocate(nlabels)
end subroutine ionout
