! Molecular Orbital PACkage (MOPAC)
! Copyright (C) 2021, Virginia Polytechnic Institute and State University
!
! MOPAC is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! MOPAC is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

    subroutine getsym (locpar, idepfn, locdep, depmul)
    use molkst_C, only: natoms, ndep, moperr, id
    use chanel_C, only : iw, ir
    use common_arrays_C, only : na
    use molkst_C, only : line
   !
   !.. Implicit Declarations ..
    implicit none
   !
   !.. Formal Arguments ..
    integer, dimension (3*natoms), intent (inout) :: idepfn, locpar
    integer, dimension (3*natoms), intent (inout) :: locdep
    double precision, dimension (natoms), intent (out) :: depmul
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------


!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(100) :: ivalue
      integer :: n, nvalue, i, ll, j, l, nerror, n_used, i_sym(4)
      double precision :: sum
      double precision, dimension(100) :: value
      character , dimension(19) :: texti*60, textx*60
      character, dimension(19,2) :: text*60
      character, dimension(38) :: used*60
      logical :: ok, sym_pdb

      save text, i_sym
!-----------------------------------------------

      equivalence (text(1,1), texti), (text(1,2), textx)
      data i_sym/3, 6, 8, 12/
      data texti/ &
        ' BOND LENGTH    IS SET EQUAL TO THE REFERENCE BOND LENGTH   ', &
        ' BOND ANGLE     IS SET EQUAL TO THE REFERENCE BOND ANGLE    ', &
        ' DIHEDRAL ANGLE IS SET EQUAL TO THE REFERENCE DIHEDRAL ANGLE', &
        ' DIHEDRAL ANGLE VARIES AS  90 DEGREES - REFERENCE DIHEDRAL  ', &
        ' DIHEDRAL ANGLE VARIES AS  90 DEGREES + REFERENCE DIHEDRAL  ', &
        ' DIHEDRAL ANGLE VARIES AS 120 DEGREES - REFERENCE DIHEDRAL  ', &
        ' DIHEDRAL ANGLE VARIES AS 120 DEGREES + REFERENCE DIHEDRAL  ', &
        ' DIHEDRAL ANGLE VARIES AS 180 DEGREES - REFERENCE DIHEDRAL  ', &
        ' DIHEDRAL ANGLE VARIES AS 180 DEGREES + REFERENCE DIHEDRAL  ', &
        ' DIHEDRAL ANGLE VARIES AS 240 DEGREES - REFERENCE DIHEDRAL  ', &
        ' DIHEDRAL ANGLE VARIES AS 240 DEGREES + REFERENCE DIHEDRAL  ', &
        ' DIHEDRAL ANGLE VARIES AS 270 DEGREES - REFERENCE DIHEDRAL  ', &
        ' DIHEDRAL ANGLE VARIES AS 270 DEGREES + REFERENCE DIHEDRAL  ', &
        ' DIHEDRAL ANGLE VARIES AS - REFERENCE DIHEDRAL              ', &
        ' BOND LENGTH VARIES AS HALF THE REFERENCE BOND LENGTH       ', &
        ' BOND ANGLE VARIES AS HALF THE REFERENCE BOND ANGLE         ', &
        ' BOND ANGLE VARIES AS 180 DEGREES - REFERENCE BOND ANGLE    ', &
        ' DO NOT USE - USE SYMMETY FUNCTION 19 INSTEAD               ', &
        ' BOND LENGTH IS A MULTIPLE OF THE REFERENCE BOND LENGTH     '/
      data textx/ &
        ' X COORDINATE IS SET EQUAL TO   THE REFERENCE X COORDINATE  ', &
        ' Y COORDINATE IS SET EQUAL TO   THE REFERENCE Y COORDINATE  ', &
        ' Z COORDINATE IS SET EQUAL TO   THE REFERENCE Z COORDINATE  ', &
        ' X COORDINATE IS SET EQUAL TO - THE REFERENCE X COORDINATE  ', &
        ' Y COORDINATE IS SET EQUAL TO - THE REFERENCE Y COORDINATE  ', &
        ' Z COORDINATE IS SET EQUAL TO - THE REFERENCE Z COORDINATE  ', &
        ' X COORDINATE IS SET EQUAL TO   THE REFERENCE Y COORDINATE  ', &
        ' Y COORDINATE IS SET EQUAL TO   THE REFERENCE Z COORDINATE  ', &
        ' Z COORDINATE IS SET EQUAL TO   THE REFERENCE X COORDINATE  ', &
        ' X COORDINATE IS SET EQUAL TO - THE REFERENCE Y COORDINATE  ', &
        ' Y COORDINATE IS SET EQUAL TO - THE REFERENCE Z COORDINATE  ', &
        ' Z COORDINATE IS SET EQUAL TO - THE REFERENCE X COORDINATE  ', &
        ' X COORDINATE IS SET EQUAL TO   THE REFERENCE Z COORDINATE  ', &
        ' Y COORDINATE IS SET EQUAL TO   THE REFERENCE X COORDINATE  ', &
        ' Z COORDINATE IS SET EQUAL TO   THE REFERENCE Y COORDINATE  ', &
        ' X COORDINATE IS SET EQUAL TO - THE REFERENCE Z COORDINATE  ', &
        ' Y COORDINATE IS SET EQUAL TO - THE REFERENCE X COORDINATE  ', &
        ' Z COORDINATE IS SET EQUAL TO - THE REFERENCE Y COORDINATE  ', &
        ' NOT USED                                                   '/
      n_used = 0
      nerror = 0
      sym_pdb = .false.
!
! TITLE OUTPUT
!
      write (iw, 10)
   10 format(/,/,/,20x,'PARAMETER DEPENDENCE DATA'/,/,&
        '        REFERENCE ATOM      FUNCTION NO.    DEPENDENT ATOM(S)')
!
! INPUT SYMMETRY : FUNCTION, REFERANCE PARAMETER, AND DEPENDENT ATOMS
!
      n = 0
      if (ndep > 0) then
        j = 2
        do i = 1, ndep
          j = j + 1
          ivalue(j) = locdep(i)
          ok = (i == ndep)
          if (.not. ok) ok = (locpar(i) /= locpar(i + 1) .or. idepfn(i) /= idepfn(i + 1))
          if (ok) then
            write (iw, "(i13, i19, i16, 6 i5, 10(/, 43 x, 7 i5))") &
            locpar(i), idepfn(i), (ivalue(l), l=3, j)
            j = 2
          end if
        end do
        goto 90
      end if
      depmul(1) = 0.D0
   20 continue
      read (ir, '(A)', iostat=i, end=90) line
      if (i/=0) then
        return
      endif
      call upcase(line, len_trim(line))
      do
        i =index(line, '"')
        if (i == 0) exit
        call txt_to_atom_no(line, i, .false.)
        call  l_control("CONTROL_SYM_in_PDB", len_trim("CONTROL_SYM_in_PDB"), 1)
        sym_pdb = .true.
        if (moperr) return
      end do
      call nuchar (line, len_trim(line), value, nvalue)
!   INTEGER VALUES
      do i = 1, nvalue
        ivalue(i) = nint(value(i))
      end do
!   FILL THE LOCDEP ARRAY
      if (nvalue==0 .or. Abs(value(3)) < 1.d-20) go to 90
      if (ivalue(2) == 19) then
        if (na(ivalue(1)) == 0) then
            !
            !  Not allowed: a Cartesian coordinate cannot use function 19
            !
          write (iw,*) "Atom ", ivalue (1), " is Cartesian.  " // &
                     & "Function 19 cannot be used here."
          call mopend ("Error in Symmetry Data")
          return
        end if
        do i = 4, nvalue
          if (ivalue(i) == 0) exit
          ndep = ndep + 1
          locdep(ndep) = ivalue(i)
          locpar(ndep) = ivalue(1)
          idepfn(ndep) = 19
          n = n + 1
!
!  Check:  Is multiplier a square root of a rational ratio?
!
          sum = value(3) ** 2
          do ll = 1, 4
            l = i_sym(ll)
            j = Nint (sum*l)
            if (Abs (j-sum*l) < l*1.d-4) then
              value(3) = Sqrt ((1.d0*j)/l)
              exit
            end if
          end do
          depmul(n) = value(3)
        end do
      else
        if (na(ivalue(1)) /= 0 .and. ivalue(2) == 18) then
            !
            !  Not allowed: an internal coordinate cannot use function 18
            !
          write (iw,*) "Atom ", ivalue (1), " is internal.  " // &
                     & "Function 18 cannot be used here."
          call mopend ("Error in Symmetry Data")
          return
        end if
        do i = 3, nvalue
          if (ivalue(i) == 0) exit
          ndep = ndep + 1
          locdep(ndep) = ivalue(i)
          locpar(ndep) = ivalue(1)
          idepfn(ndep) = ivalue(2)
          if (ivalue(i) > natoms + id) then
            nerror = 1
          end if
        end do
      end if
      ll = i - 1
      if (ivalue(2) == 19) then
          write (iw,'(i13, i13, f13.8, i9, 6 i5, 10(/, 43 x, 7 i5))') &
          ivalue(1), ivalue(2), value(3), (ivalue(j), j=4, ll)
      else
        if (sym_pdb) then
          line = " "
          call atom_no_to_txt(ivalue(1), line(9:))
          write(line(len_trim(line) + 1:),'(i13)')ivalue(2)
          call atom_no_to_txt(ivalue(3), line(len_trim(line) + 11:))
          do j = 4, ll
            call atom_no_to_txt(ivalue(j), line(len_trim(line) + 2:))
          end do
          write(iw, '(a)')trim(line)
        else
          write (iw, "(i13, i19, i16, 6 i5, 10(/, 43 x, 7 i5))") ivalue(1), ivalue(2), (ivalue(j), j=3, ll)
        end if
      end if
      if (na(ivalue(1)) == 0) then
        i = 2
      else
        i = 1
      end if
      line = text(ivalue(2), i)
      do i = 1, n_used
        if (used(i) == line) exit
      end do
      if (i > n_used) then
        n_used = n_used + 1
        used(n_used) = line(:len_trim(line))
      end if
      if (nerror == 1) then
        call mopend("A SYMMETRY FUNCTION IS USED TO DEFINE A NON-EXISTENT ATOM")
        return
      end if
    goto 20
!
! CLEAN UP
   90 continue
      write (iw, 100)
  100 format(/,10x,'   DESCRIPTIONS OF THE FUNCTIONS USED')
      do j = 1, 18
        do i = 1, n_used
          if (used(i) == texti(j)) go to 120
        end do
        cycle
  120   continue
         write (iw, 130) j, used(i)
  130   format(i4,5x,a)
      end do
      do j = 1, 18
        do i = 1, n_used
          if (used(i) == textx(j)) go to 121
        end do
        cycle
  121   continue
         write (iw, 130) j, used(i)
      end do
      write(iw,*)
      return
      end subroutine getsym
