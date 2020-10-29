subroutine newflg ()
    use molkst_C, only: natoms, numat, nvar
    use common_arrays_C, only : coord, geo, txtatm, na, nb, nc, loc, xparam
    use chanel_C, only: iw
    implicit none
    integer :: i, j, k
!
    if (numat /= natoms) then
      write (iw,*) " NEWGEO CAN ONLY BE USED IF THERE ARE NO DUMMY ATOMS,", &
     & " OR IF 'RESEQ' IS USED"
      call mopend ("NEWGEO cannot be used here")
    end if
   !
   !    Check that there are internal coordinates
   !
    do i = 1, numat
      if (na(i) /= 0) go to 1000
    end do
    write (iw,*)
    write (iw,*) "  There is a bug in MOZYME.  For unknown reasons,"
    write (iw,*) &
    &"  NEWGEO cannot be used here. (All the coordinates are Cartesian already!"
    call mopend ("There is a bug in MOZYME")
1000 k = 0
    do i = 1, numat
      if (txtatm(i)(:6) /= "ATOM  ") cycle
      if (txtatm(i)(13:15) /= " N " .and. txtatm(i)(13:15) /= " C ") cycle
!
!  CONVERT ATOM COORDINATES FROM INTERNAL TO CARTESIAN
!
      na(i) = 0
      nb(i) = 0
      nc(i) = 0
      do j = 1, 3
        geo(j, i) = coord(j, i)
      end do
      k = 1
    end do
    do i = 1, nvar
      xparam(i) = geo(loc(2, i), loc(1, i))
    end do
    if (k == 0) then
      write (iw,*) " WARNING!  NO BACKBONE ATOMS IDENTIFIED"
    end if
end subroutine newflg
