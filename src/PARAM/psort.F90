! Molecular Orbital PACkage (MOPAC)
! Copyright 2021 Virginia Polytechnic Institute and State University
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!    http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

subroutine psort (refers, mode, nmols, reftxt, allref, ifile, l, lp)
    use param_global_C, only : source, maxrfs, contrl
    implicit none
    logical, intent (in) :: lp
    integer, intent (in) :: ifile, mode, nmols
    integer, intent (out) :: l
    character (len=5), dimension (nmols), intent (out) :: reftxt
    character (len=8), dimension (nmols), intent (in) :: refers
    character (len=13), dimension (maxrfs), intent (inout) :: allref
    integer, parameter :: nrefs = 2
    character :: line*120, iw*300, location_of_references(nrefs)*22
    integer :: i, ii, j, m
    logical :: first = .true., exists
    intrinsic :: Index
    data location_of_references / &
    "./REFERENCES.txt     ", &
    "M:/utility/REFERENCES" /
!----------------------------------------------------------------------
    if (mode /= 1) then
      l = 0
      do i = 1, nmols
        if (refers(i) == " ") then
          reftxt (i) = "  "
        else
          do j = 1, l
            if (allref(j) (6:13) == refers(i)) go to 1000
          end do
          l = l + 1
          j = l
          allref(j) (6:13) = refers(i)
1000      reftxt(i) = " "
          write(reftxt(i),'(i5)')j
          allref(j) (1:5) = reftxt(i)
        end if
      end do
      if ( .not. lp) return
    end if
    i = index(contrl, 'REFERENCES="') + 12
    if (i /= 12) then
      j = index(contrl(i:), '"') + i - 2
      iw = trim(contrl(i:j))
      call add_path(iw)
      inquire (file=trim(iw), exist = exists)
    else
      do i = 1, nrefs
        inquire (file=trim(location_of_references(i)), exist = exists)
        if (exists) then
          iw = trim(location_of_references(i))
          exit
        end if
      end do
    end if
    if (exists) then
      open (14, file=trim(iw), blank="ZERO")
      rewind (14)
    end if
    m = l
    do i = 1, l
      source (3, i) = "      *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+* +&
     &* "
      source (2, i) = "      *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+* +&
     &* "
      source (1, i) = " " // allref (i) (5:12) // " NO REFERENCE *+*+*+*+" // &
     & "*+*+*+*+*+*+*+*+*+*"
    end do
    do ii = 1, 10000
      read (14, "(A)", end=1200, err=1200) line
!
! Having read in one reference, is that reference used by this run?
! Check mnemonics only
!
      do j = 1, l
        i = Index(line, " "//allref(j) (6:13))
        if (i /= 0) go to 1100
      end do
!
! Have references become messed up?  References mnemonics should have less than 10 characters
!
      if (line(11:) /= " ") then
        if(first) then
          first = .false.
          write(ifile,*)" Some mnemonics apear to be out-of sequence.  These include:"
        end if
        write(ifile,"(a)")line
      end if
!
! Reference is NOT used, therefore skip over next three lines
!
      read (14, "(A)", end=1200, err=1200) line, line, line
      cycle
1100  read (14, "(A)", end=1200, err=1200) source (1, j), source (2, j), &
     & source(3, j)
      m = m - 1
    end do
1200 close (14)
    do i = 1, l
      write (ifile,*)
      iw = source(1, i)
      if (source(2, i) /= " ") then
        iw(len_trim(iw) + 2:) = source(2, i)
        if (source(3, i) /= " ") then
          iw(len_trim(iw) + 2:) = source(3, i)
        end if
      end if
      if (index(iw, " NO REFERENCE") == 0) then
        write (ifile, "(A)")  allref (i) (1:5) // ":" // trim(iw)
      else
        write (ifile, "(A)")  allref (i) (1:5) // ": " // trim(iw(3:))
      end if
    end do
end subroutine psort
