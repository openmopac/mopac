subroutine savgeo (loop, geo, na, nb, nc, xparam, loc)
!
    use param_global_C, only : contrl, nas, nbs, ncs, geos, &
    names, refher !, ifiles_8
!
    use molkst_C, only : natoms, nvar, keywrd, title, refkey, escf, &
      line, ncomments
!
    use common_arrays_C, only : all_comments

    implicit none
    integer, intent (in) :: loop
    integer, dimension (natoms), intent (in) :: na, nb, nc
    integer, dimension (2, nvar), intent (in) :: loc
    double precision, dimension (nvar), intent (in) :: xparam
    double precision, dimension (3, natoms), intent (inout) :: geo
!--------------------------------------------------------------------
    character (len=100) :: dirnew
    character (len=100) :: name
    logical :: opend, lsav !, exists
    integer :: i, iatm, igeo, j, k, iw
    character :: num*1
    intrinsic Index
    double precision, external :: reada
    character, external :: get_a_name*300
    save :: iatm, igeo
!--------------------------------------------------------------------
  !
  !   Check - should the geometry be saved.  Criterion: NEW_REF exists.
  !   "lsav" is TRUE if the geometry is to be saved.
      lsav = .false.
      iw = 14
      k = Index (contrl, " NEW_REF=")
      if (k /= 0) then
        dirnew = get_a_name(contrl(k + 9:), len_trim(contrl(k + 9:)))
        k = len_trim(dirnew)
        if (dirnew(k:k) /= "/")then
!
! The directory name for the new reference data needs a "/"
!
          k = k + 1
          dirnew(k:k) = "/"
        end if
        do i = 80, 2, -1
          if (names(loop)(i:i) /= " ") exit
        end do
        name = trim(dirnew) // names (loop) (:i) // ".mop"
        inquire (unit=iw, opened=opend)
        if (opend) then
          close (unit=iw, status="KEEP")
        end if
        call add_path(dirnew)
! the directory feature of inquire is specific to the Intel compiler and not standard Fortran, ignoring this for now
!        inquire (directory=trim(dirnew) , exist = exists)
!        if (.not. exists) then
!          write(ifiles_8,'(//10x,a)')"Folder """//trim(dirnew)//""" specified by NEW_REF does not exist"
!          write(ifiles_8,'(/10x,a)')"(Either edit or delete ""NEW_REF"" or create the folder.)"          
!          stop
!        end if
        call add_path(name)
        open (iw, status="UNKNOWN", file=name, blank="ZERO", err=1000)
        write (iw, "(A)", err=1000) " Test"
        rewind (iw)
        lsav = .true.
      end if
!
!  keywrd is 248 characters long, refkey is 360 characters long
!
 1000 refkey(1) = trim(keywrd)
     do i = 1, nvar
      geo(loc(2, i), loc(1, i)) = xparam(i)
     end do
     line = title
     if (Index (contrl, " UPDATE ") /= 0) then
!
! Replace the old reference heat of formation with the calculated HoF.
!
       i = index(title, "H=")
       if (i /= 0) then
         do j = i + 2, i + 20
           if ((ichar(title(j:j)) < ichar("0") .or. ichar(title(j:j)) > ichar("9")) .and. &
              ichar(title(j:j)) /= ichar("-") .and. ichar(title(j:j)) /= ichar(".")) exit 
         end do
         line(100:) = title(j:) ! j = end of number
         if (line(i + 3:i + 3) == "-") then
           line(i + 3:99) = " "
         else
            line(i + 2:99) = " "
         end if
         if (title(j:j) == "+") then
!
! Re-set escf to difference between calculated error (refher) and that supplied 
!
           escf = reada(title, i + 2) + refher
         end if
         if (escf > 0.d0) then
           i = max(int(log10(escf)),0)
           num = char(ichar("5") + i)
         else
           i = max(int(log10(-escf)),0)
           num = char(ichar("6") + i)
         end if
         write(title,"(a,f"//num//".3,a)")trim(line(:99)), escf, trim(line(100:199))         
       end if
     end if
  !
  !  Only write out the geometry if the geometry has been optimized.
  !
    if (lsav) then
      if (Index (contrl, " EXPORT ") /= 0) title = " "
      if (index(contrl," LET") /= 0) then
          i = index(refkey(1)," GEO-OK")
          if (i /= 0) refkey(1)(i:i+7) = " "
        end if
      ncomments = 3
      all_comments(1) = "*"
      write(all_comments(2),"('*', f10.3)") escf
      all_comments(3) = "*"
      call geout (iw)
      ncomments = 0
      close (iw)
    end if
    if (loop == 1) then
      igeo = 0
      iatm = 0
    end if
    do i = 1, natoms
      igeo = igeo + 1
      do j = 1, 3
        geos(j, igeo) = geo(j, i)
      end do
    end do
    do i = 1, natoms
      iatm = iatm + 1
      nas(iatm) = na(i)
      nbs(iatm) = nb(i)
      ncs(iatm) = nc(i)
    end do
end subroutine savgeo
