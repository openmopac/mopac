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
!    Force all coordinates to be internal
!
    call xyzint (coord, numat, na, nb, nc, 1.D0, geo)
!
!  Convert coordinates to Cartesian
!
    call gmetry(geo, coord)
    k = 0
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
    if (k == 0) call mopend("KEYWORD ""NEWGEO"" USED, BUT NO BACKBONE ATOMS IDENTIFIED")
    return
end subroutine newflg
