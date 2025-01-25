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

subroutine point (nvecs, fcells)
  use common_common, only : id, ir, iw, line, Sym_Oper, jobnam, data_set_name, l_read_dat
  implicit none
  integer, intent(in) :: nvecs
  complex, dimension(nvecs, nvecs) :: v
  double precision, dimension(nvecs, *), intent(in) :: fcells
!
  character(len=80) :: input_line
  logical :: leadsp
  integer :: i, j, itype, nvalue
  integer, dimension(20) :: istart
  double precision, dimension(3) :: xyzk, xyzxyz
  double precision, dimension(2*nvecs) :: eigs
  double precision, external :: reada
  do
    Sym_Oper = " "
    if (l_read_dat) write(iw, '(a)') " Please type point in k-space or a blank line to quit"
    if (l_read_dat) write(iw, "(A,I2,A)") " Format is", id, " numbers, free format"
    read(ir, "(A)", end = 1000, err = 1000) line
    if (line == " ") goto 1000
    i = index(line, '{')
    if (i /= 0) then
      j = index(line, "}")
      Sym_Oper = line(i + 1:j - 1)
    end if
    line(30:) = ": Please type point in k-space or a blank line to quit"
    call write_keystrokes(line, len_trim(line))
    read(line, '(a)') input_line
    leadsp = .true.
    nvalue = 0
    do i = 1, 80
      if (leadsp .and. input_line(i:i)/=" ") then
        nvalue = nvalue + 1
        istart(nvalue) = i
      end if
      leadsp = (input_line(i:i)==" ")
    end do
    if( nvalue < id .and. nvalue > 0) then
      i = ichar(input_line(istart(1):istart(1)))
      if (i >= ichar("0") .and. i <= ichar("9")) &
        write(iw,'(a,i2,a)')" Exactly",id," points are needed" 
      cycle    
    end if
    if (nvalue > 3) then
      line = input_line(istart(4):)
    else
      line = " "
    end if
    if (input_line == " ") stop
    do i = 1, id
      xyzk(i) = reada (input_line, istart(i))
    end do
    do i = 1, 3
      xyzxyz(i) = xyzk(i)
    end do
    if (id /= 1) then
      call rot (xyzk)
    end if
    call upcase (input_line)
    !
    !   Option to decide how to represent operations which involve non-
    !   primitive translations
    !
    itype = 1
    if (Index (input_line, " RAW") /= 0) then
      itype = 1
    end if
    if (Index (input_line, " ONE") /= 0) then
      itype = 2
    end if
    if (Index (input_line, " BYK") /= 0) then
      itype = 3
    end if
    if (input_line == " ") then
      stop
    end if
    call kpoint (fcells, nvecs, xyzk(1), xyzk(2), xyzk(3), v, eigs)
    line = Sym_Oper
    call solir (eigs, v, nvecs, xyzk, itype, xyzxyz)
  end do
1000 continue
  call graphics(0.0, 0.0, 100)
  if (l_read_dat) then
    write(iw,'(/5x,a)')"A copy of this output can be found in the file: '"//trim(jobnam)//".txt'"
    write(iw,'(/5x,a)')"Press ""RETURN"" to quit"
    read(iw,'(a)', iostat = i) input_line
  else
    write(iw,'(/5x,a)')"A copy of this output can be found in the file: '"//trim(data_set_name)//".txt'"
    write(iw,'(/5x,a)')"Press ""RETURN"" to quit"
    read(iw,'(a)', iostat = i) input_line
  end if
  stop
end subroutine point
