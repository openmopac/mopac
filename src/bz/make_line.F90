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

subroutine make_line (nvecs, fcells)
!
! 
!
  use common_common, only : id, ir, iw, data_set_name, mr1, mr2, mr3, phonon, iw_new, jobnam, line, &
    l_read_dat
  implicit none
  integer, intent(in) :: nvecs
  double precision, dimension(nvecs, nvecs*(mr1*mr2*mr3 + 1)), intent(in) :: fcells
!
  complex, allocatable, dimension(:) :: v
  double precision, dimension(:,:), allocatable :: alleig, all_points
  integer :: i, j, k, l, m, mm, npts, point
  real :: x, y, scale=0.5, old_dist, ratio, bottom
  double precision :: dist, stepy, xkold = 0.d0, del_point(4)
  double precision :: emax, emin, erange, step, stepx, stepz, xk, xk1, &
       & xk2, xofset, yk, yk1, yk2, zk, zk1, maximum, minimum, reset_zero = -0.2d0, last_x
  double precision :: ykold = 0.d0, zk2, zkold = 0.d0, distance = 0.d0
  double precision :: eigs(2*nvecs), xyzk(3,2)
  integer, parameter :: npoints = 4
  character :: num*1, char_points(3)*3, symbols*100, line_1*100, &
    points(npoints)*6 =  (/"GAMMA", "DELTA", "LAMBDA", "SIGMA" /), &
    pts(npoints)*3 = (/"G G", "D G", "L G", "S G"/)
  logical :: opend
  save
! 
! ... Executable Statements ...
! 
  allocate(alleig(nvecs, 0:10000))
  allocate(v(nvecs**2)) 
  allocate(all_points(4, 0:10000))
  call graphics (0.0, 0.0, 6)
  xofset = 0.d0
  last_x = 0.d0
  emin = 100.d0
  emax = -100.d0
  point = 0
  data_set_name = trim(data_set_name)//".dat"
  do
    if (l_read_dat) then
      write (iw, '(a)') " Please type starting and ending points in k-space"
      write (iw, "(A,I2,A)") " Format is", 2*id, " numbers, free format, all on one line"
      write (iw, *) " To quit, type 'quit'."
    end if
    read (ir, '(a80)', end = 1000, err = 1000) line
    i = len_trim(line)
!
!  Work out the starting, middle, and terminal point symbols for the current line.
!
    line_1 = trim(line)
    call upcase(line)
    symbols = trim(line)
    char_points = " "
    l = 0
    m = 2
    do k = 1, 3
      symbols = " "
      do j = m - 1, len_trim(line)
        if (line(j:j) >= "A" .and. line(j:j) <= "Z") then
          do m = j, len_trim(line) + 1
            if (line(m:m) == " ") then
              symbols = " "//line(j:m)
              line(j:m) = " "
              exit
            end if
          end do
        end if        
        if (symbols /= " ") exit
      end do   
      if (len_trim(symbols) == 3 .and. (symbols(3:3) >= "0" .and. symbols(3:3) <= "9" .or. &
        symbols(3:3) == """" .or. symbols(3:3) == "'") .or. len_trim(symbols) == 2) then
        l = l + 1
        char_points(l) = symbols(2:3)//"R"
      else       
        do j = 1, npoints
          mm = index(symbols, " "//trim(points(j)))
          if (mm /= 0) exit
        end do    
        if (j <= npoints) then
          l = l + 1
          char_points(l) = pts(j)
          mm = mm + 1 + len_trim(points(j))
          if (symbols(mm:mm) == "'") char_points(l)(2:2) = char(162) ! Greek symbol for "'"
        end if
      end if
    end do
    if (l == 2) then
      char_points(3) = char_points(2)
      char_points(2) = " "
    end if
!
! Manipulate the starting and ending points of the line in k-space.
!
    line = trim(line_1)
    if (i < 30) write(line(30:),'(a, i2, a)')": Starting and ending points in k-space, ", 2*id, " numbers, ""quit"" to quit"
    call write_keystrokes(line, len_trim(line))    
    read(line,*, end = 1000, err = 1000) ((xyzk(i, j), i = 1, id), j = 1, 2)
    zk2 = xyzk(3,2)
    yk2 = xyzk(2,2)
    xk2 = xyzk(1,2)
    zk1 = xyzk(3,1)
    yk1 = xyzk(2,1)
    xk1 = xyzk(1,1)
    old_dist = real(Sqrt ((xk1 - xk2)**2 + (yk1 - yk2)**2 + (zk1 - zk2)**2))
    if (id /= 1) then
      do i = 1, 2
        call rot (xyzk(1, i))
      end do
    end if
    zk2 = xyzk(3,2)
    yk2 = xyzk(2,2)
    xk2 = xyzk(1,2)
    zk1 = xyzk(3,1)
    yk1 = xyzk(2,1)
    xk1 = xyzk(1,1)
    dist = Sqrt ((xk1 - xk2)**2 + (yk1 - yk2)**2 + (zk1 - zk2)**2)
    if (dist < 1.d-5) exit
    ratio = real(old_dist/dist)
    if (xofset /= 0.d0) then
      if (Abs (xk1 - xkold) > 1.d-2 .or. Abs (yk1 - ykold) > 1.d-2 .or. Abs (zk1 - zkold) > 1.d-2) then
!
!  Discontinuity in Band Structure - put in a small gap
!
        xofset = xofset + 0.1d0
        distance = distance + 0.1d0
      end if
    end if
    npts = Nint (300*dist)
    stepx = (xk2 - xk1) / npts
    stepy = (yk2 - yk1) / npts
    stepz = (zk2 - zk1) / npts
    xk = xk1 - stepx
    yk = yk1 - stepy
    zk = zk1 - stepz
    do i = 0, npts
      xk = xk + stepx
      yk = yk + stepy
      zk = zk + stepz
      call kpoint (fcells, nvecs, xk, yk, zk, v, eigs)
      point = point + 1
      do j = 1, nvecs
        alleig(j, point) = eigs(j)
      end do
      if (i > 0)  distance = distance + old_dist/npts
      all_points(1,point) = distance
      all_points(2,point) = xk
      all_points(3,point) = yk
      all_points(4,point) = zk
      if (xofset == 0.d0) then
        do j = 1, nvecs
          emin = Min (emin, eigs(j))
          emax = Max (emax, eigs(j))
        end do
      end if
    end do
    erange = 0.505d0*(emax - emin)
    step = dist / npts
    if (id == 1) scale = scale*4.0
!
!  Write scale on left-side of the walk
!
    if (x < 0.01) call write_scale(emin, emax, erange, reset_zero, bottom)
    x = Real (xofset)
!
!  Left-hand side vertical bar.
!
    call graphics (x*scale, -0.05, 2)
    call graphics (x*scale, bottom, 3)
!
!  Write the Brillouin zone symbol
!
    write(line,'(a)')char_points(1)
    call graphics (x*scale - 0.01, bottom, 96) ! Symbol of Brillouin-zone point
    last_x = x
!
!   Draw the plot.
!
    do i = 1, nvecs
      y = 1.95 - Real ((alleig(i, point - npts) - emin + reset_zero)/erange)
      x = Real (0.d0 + xofset)
      call graphics (x*scale,  y, 2)
      do j = 1,npts
      k = point - npts + j 
      y = 1.95 - Real ((alleig(i, k) - emin + reset_zero)/erange)
      x = Real (j*step*ratio + xofset)
      call graphics (x*scale,  y, 3)
      end do
    end do
!
!  Normal right-hand side vertical bar.
!
    call graphics (x*scale, -0.05, 2)
    call graphics (x*scale, bottom, 3)
    write(line,'(a)')char_points(2)
    call graphics (real(x + last_x)*0.5*scale - 0.01, bottom, 96) ! Symbol of Brillouin-zone point
    write(line,'(a)')char_points(3)
    call graphics (x*scale - 0.01, bottom, 96) ! Symbol of Brillouin-zone point
!
!   Draw the bottom line
!
    x = Real (xofset)
    call graphics (x*scale, bottom, 2)
    x = x + Real (dist)*ratio
    call graphics (x*scale, bottom, 3)
    last_x = xofset
    xofset = xofset + dist*ratio
    xkold = xk2
    ykold = yk2
    zkold = zk2
  end do
1000 x = Real (xofset)
  maximum = -1.d9
  minimum = 1.d9
  do j = 1,point
    maximum = max(maximum, alleig(nvecs,j))
    minimum = min(minimum, alleig(1,j))
  end do
  i = int(maximum - minimum)
  if (i > 100) then !  Work out step-size
    j = 20
  else if (i > 50) then
    j = 10
  else if (i > 20) then
    j = 5
  else 
    j = 2
  end if
  maximum = int(maximum/j)*j + j 
  minimum = minimum + 10000
  minimum = int(minimum/j)*j 
  minimum = minimum - 10000
  if (phonon) then
    num = "1"
  else
    num = "3"
  end if
  del_point(:) = all_points(:,2) - all_points(:,1)
!
! Write everything to a file, then stop
!
  inquire(unit=iw_new, opened=opend) 
  if (opend) close (iw_new)
98 open (unit = iw_new, file = trim(jobnam)//".txt", form = "FORMATTED", status = "UNKNOWN", iostat = i)
  if (i /= 0) then
    write(iw,*)" File: '"//trim(jobnam)//".txt' cannot be opened."
    call sleep(1)
    goto 98
  end if
  write(iw_new,'(f8.3,f8.'//num//')',iostat = i)0.0, maximum
  if (i /= 0) then
    write(iw,*)" Cannot write to file: '"//trim(jobnam)//".txt'."
    call sleep(1)
    goto 98
  end if
  do j = 1, point
    if (j > 1) then
      step = abs(del_point(2) - all_points(2,j) + all_points(2,j - 1)) + &
             abs(del_point(3) - all_points(3,j) + all_points(3,j - 1)) + &
             abs(del_point(4) - all_points(4,j) + all_points(4,j - 1)) 
      if (step > 0.0001d0) then
        if (step > 0.01d0) then
          write(iw_new,'(f8.3,100f8.'//num//')')all_points(1,j - 1),minimum
          write(iw_new,'(f8.3,100f8.'//num//')')all_points(1,j - 1),maximum
          write(iw_new,'(f8.3,100f8.'//num//')')all_points(1,j - 1),minimum
          write(iw_new,'(f8.3,100f8.'//num//')')all_points(1,j - 1),minimum
          write(iw_new,'(f8.3,100f8.'//num//')')all_points(1,j),minimum
        end if
        write(iw_new,'(f8.3,100f8.'//num//')')all_points(1,j),minimum
        write(iw_new,'(f8.3,100f8.'//num//')')all_points(1,j),maximum
        write(iw_new,'(f8.3,100f8.'//num//')')all_points(1,j),minimum
        del_point(:) = all_points(:,j + 1) - all_points(:,j)
      end if
    end if
    write(iw_new,'(f8.3,100f8.'//num//')')all_points(1,j),minimum, (alleig(i,j), i = 1, nvecs)
  end do
  write(iw_new,'(f8.3,100f8.'//num//')')all_points(1,j - 1),maximum
  write(iw,'(/5x,a)')"A copy of this output can be found in the file: '"//trim(jobnam)//".txt'"
  write(iw,'(/5x,a)')"Press ""RETURN"" to quit"
  read(iw,'(a)', iostat = i) num
  call graphics(0.0, 0.0, 100)
  stop
  end subroutine make_line
  
  
  subroutine write_scale(emin, emax, erange, reset_zero, bottom)
  use common_common, only : line
  implicit none
  double precision, intent (in) :: emin, emax, reset_zero, erange
  real, intent (out) :: bottom
  integer, parameter :: nsteps = 12
  real :: lower_bound, upper_bound, zero, y, marker, step, range, &
  steps(nsteps) = (/5000.0, 2000.0, 1000.0, 500.0, 200.0, 100.0, 50.0, 20.0, 10.0, 5.0, 2.0, 1.0/)     
  character :: num*1
  integer :: i, npts, num_steps, mini
    lower_bound = real(emin)
    upper_bound = real(emax)
    zero = real(reset_zero)
    range = upper_bound - lower_bound
    num_steps = 10
    do i = 1, nsteps
      if (range*2.5 > steps(i)) then
        step = steps(i)/num_steps
        marker = nint(lower_bound/step)*step - step
        exit
      end if
    end do
    range = real(erange)
    marker = nint(lower_bound/step)*step 
    bottom = 1.95 - (marker - lower_bound + zero)/range
    marker = marker - step
    do i = 0, 20
      marker = marker + step
      if (marker > upper_bound) exit
      num = char(ichar("2") + int(log10(max(abs(marker), 1.0))))
      write(line,'(i'//num//')') nint(marker)
      mini = int(max(log10(abs(marker) + 0.1), 0.0))
      y = 1.95 - (marker - lower_bound + zero)/range
      call graphics (-0.06 - mini*0.03, y - 0.04, 96) ! Write scale
      call graphics (-0.01, y, 2) ! Write tic
      call graphics (-0.00, y, 3) ! Write tic
    end do
  return
  end subroutine write_scale
