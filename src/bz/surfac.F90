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

subroutine surfac (nvecs, fcells)
!
  use common_common, only : id, phonon, ir, iw, jobnam, data_set_name, iw_new, line, &
    l_read_dat, top_l, top_r, bottom_l, bottom_r, grid, isurf
  implicit none
!
!  npts should be odd, so that the system has a center.
!
  integer, parameter :: npts = 81
  integer, intent(in) :: nvecs
  double precision, dimension(nvecs, *), intent(in) :: fcells
!
  complex, allocatable, dimension(:) :: v
  integer :: i, i99, iok, irow, j, k, mpts, npts2, nmax, nmin, l
  double precision :: ca, cb, cntval, cx, cy, cz, px, py, pz, r1, range, &
       & sa, sb, size, step, x, x2, xmax, xmin, xy, y2, eigs(2*nvecs), &
       xyzc(3), xyzp(3), top(3), bottom(3), stepx, stepy, xyz(3), &
       temp_ima(6),  maxima(6,1000), minima(6,1000), store_xyzc(3)
  logical :: first = .true., opend, exists
  character :: q*1, graph_name*200
  double precision, allocatable :: onebnd(:)
  external sleep
! 
!.. Equivalences .. 
  equivalence (xyzp(3), pz), (xyzp(2), py), (xyzp(1), px)
  equivalence (xyzc(3), cz), (xyzc(2), cy), (xyzc(1), cx)
  allocate(grid(npts,npts,nvecs), v(nvecs**2))
  allocate(onebnd((npts*2)**2))
! 
! ... Executable Statements ...
! 
  step = 1.0d0 / (npts-1)
  iok = 1
  inquire (file=jobnam(:len_trim(jobnam))//".bnd", exist = exists)
  if (exists) then
    write (iw, *) "Start a new surface [1], or use an old one [2]?"
    read(6, '(a)',  iostat=i) line
    if (.not. l_read_dat .and. index(line, ": Enter k-space") /= 0) then
      rewind (ir)
      read(ir, '(a)',  iostat=i) line
    else       
      line(20:) = ": Start a new surface [1], or use an old one [2]?"
      call write_keystrokes(line, len_trim(line)) 
      read(line, *) iok
      if (iok == 1 .and. .not. l_read_dat) then
        rewind (ir)
        read(ir,*,  iostat=i) line
      end if 
    end if        
  end if
  call graphics (0.0, 0.0, 6)
  if (iok == 2) then
    i = len_trim(jobnam) 
    open (unit = 17, file = jobnam(:i)//".bnd", status = "OLD", &
         & form = "FORMATTED", iostat = i99)
    if (i99 /= 0) then
      write (iw, *) " Failed to read data. "
      write(iw,*)" Press return to exit"
      read(5,'(a)') q
      call graphics(0.0, 0.0, 100)
      stop
    end if
    do i = 1, 5
      read(17,*)q
    end do
    read(17,"(3f12.5)")px, py, pz
    read(17,"(3f12.5)")store_xyzc
    read(17,"(f12.5)")size
    write(iw,"(3f12.5,a)")px, py, pz, " = center of plot"
    write(iw,"(3f12.5,a)")store_xyzc, " = vector perpendicular to plot"
    write(iw,"(/,f12.5,a)")size, " = length from center to the middle of one edge, in reciprocal space units"
    read(17,*)q
    read(17,*)top_l
    read(17,*)top_r
    read(17,*)bottom_l
    read(17,*)bottom_r
    read(17,*)i,j
    read (17,"(8f10.4)")  (((grid(i, j, k), i = 1, npts), j = 1, npts), k = 1,nvecs)
  else
    do k = 1, 10
      write (iw, '(a,i2,a)') "Enter k-space coordinates for the center of the plot,", id, " numbers"
      read (6, "(Q,A)", iostat = j)i, line 
      read(line,*, iostat = j) px, py, pz
      if (j == 0) exit
      i = len_trim(line)
      line(20:) = ": Enter k-space coordinates for the center of the plot"
      call write_keystrokes(line, len_trim(line)) 
      read(line,*, iostat = j)(xyzp(i), i = 1, id)
      if (j == 0) exit
    end do
    do k = 1, 10
      if (id == 3) then
        write (iw, '(a,i2,a)') "Enter vector perpendicular to plot,", id, " numbers"
        read(6, '(a80)', iostat = j) line
        i = len_trim(line)
        line(20:) = ": Enter vector perpendicular to plot"
        call write_keystrokes(line, len_trim(line)) 
        read(line, *, iostat = j) (xyzc(i), i = 1, id)
        if (j == 0) exit
      else
        xyzc(1) = 0.d0
        xyzc(2) = 0.d0
        xyzc(3) = 1.d0
        exit
      end if
    end do
    store_xyzc = xyzc
    if (abs(xyzc(1)) + abs(xyzc(2)) + abs(xyzc(3)) < 1.d-4) &
      write(iw,'(a,/)')" By default, the vector (0.0,0.0,1.0) will be used"
    do k = 1, 10
      write(iw, '(a)') "Enter length from center to the middle of one edge in "
      write(iw,'(a)') "reciprocal space units (typically 1.0)"
      read(6, '(a80)', iostat = j) line
       i = len_trim(line)
       line(20:) = ": Enter length from center to the middle of one edge in reciprocal space units"
      call write_keystrokes(line, len_trim(line)) 
      read(line, *, iostat = j) size
      if (size < 1.d-4) cycle
      if (j == 0) exit
    end do
!
! CALCULATE DIRECTION COSINES AND SINES.
!
    if (Abs (cx) + Abs (cy) < 1.d-5) cz = 1.d0
    xy = cx*cx + cy*cy
    if (xy > 1.d-15) then
      r1 = Sqrt (xy + cz*cz)
      xy = Sqrt (xy)
      ca = cx/xy
      cb = cz/r1
      sa = cy/xy
      sb = xy/r1
    else
      ca = 1.d0
      cb = 1.d0
      sa = 0.d0
      sb = 0.d0
    end if
    step = 2.d0*size/(npts - 1.d0)
    do i = 1, npts
      x2 = -size + (i - 1)*step
      if (mod(i,10) == 0) write (iw, *) "Raster line", i, " done,", npts-i, " left"
      do j = 1, npts
        y2 = -size + (j-1)*step
        cx = ca*cb*x2 - sa*y2 + px
        cy = sa*cb*x2 + ca*y2 + py
        cz = -sb*x2 + pz
        call rot (xyzc)
        call kpoint (fcells, nvecs, cx, cy, cz, v, eigs)
        do k = 1, nvecs
          grid(npts - i + 1,  j, k) = eigs(k)
        end do
      end do
    end do
    i = len_trim(jobnam)
    if (abs(px) > 1.d4 .or. abs(px) > 1.d4 .or. abs(px) > 1.d4) then
      write(iw,'(//10x,a)')"Data for writing the file "//jobnam(1:i)//".bnd are faulty"
      write(iw,'(10x,a)')"The current job will be halted as no useful results can be made"
      write(iw,'(/5x,a)')"Press ""RETURN"" to quit"
      read(6,'(a)', iostat = i99) q
      call graphics(0.0, 0.0, 100)
      stop
    end if
    bottom_r(1) = ca*cb*( size) - sa*(-size) + px
    bottom_r(2) = sa*cb*( size) + ca*(-size) + py
    bottom_r(3) =   -sb*( size)              + pz
    bottom_l(1) = ca*cb*(-size) - sa*(-size) + px
    bottom_l(2) = sa*cb*(-size) + ca*(-size) + py
    bottom_l(3) =   -sb*(-size)              + pz
    top_r(1) = ca*cb*( size) - sa*( size) + px
    top_r(2) = sa*cb*( size) + ca*( size) + py
    top_r(3) =   -sb*( size)              + pz
    top_l(1) = ca*cb*(-size) - sa*( size) + px
    top_l(2) = sa*cb*(-size) + ca*( size) + py
    top_l(3) =   -sb*(-size)              + pz
    open (unit = 17, file = jobnam(:i)//".bnd", status = "UNKNOWN", &
         & form = "FORMATTED", iostat = i99)
    write(17,"(2a)", iostat = i99)" Brillouin Zone Cross Sections for system: ",jobnam(1:i)
    if (i99 == 0) then
      if( phonon ) then
        write(17,'(A)')" Values of Energy Level, in cm**(-1) for the ", &
        & " Different Vibrational Bands (generated by BZ)"
      else
        write(17,'(a)')" Values of Energy Level, in eV for the ", &
        & " Different Electronic Bands (generated by BZ)"
      end if
      write(17,"(a,i4,a)")" There are",npts," points on each edge of the square"
      write(17,"(a,i4,a)")" There are",nvecs," squares (surfaces)"
      write(17,"(3f12.5,a)")px, py, pz, " = center of plot"
      write(17,"(3f12.5,a)")store_xyzc, " = vector perpendicular to plot"
      write(17,"(f12.5,a)")2*size, " = length of one edge in reciprocal space units"
      write(17,*)" (Gamma to nearest other Gamma is one unit)"
      write(17,"(3f12.5,a)")top_l, " = top left"
      write(17,"(3f12.5,a)")top_r, " = top right"
      write(17,"(3f12.5,a)")bottom_l, " = bottom left"
      write(17,"(3f12.5,a)")bottom_r, " = bottom right"
      write(17,*)npts, nvecs
      write (17,"(8f10.4)")  (((grid(i, j, k), i = 1, npts), j = 1, npts), k = 1, nvecs)
      write(17,"(3f12.5,a)")top_l, " = top left"
      write(17,"(3f12.5,a)")top_r, " = top right"
      write(17,"(3f12.5,a)")bottom_l, " = bottom left"
      write(17,"(3f12.5,a)")bottom_r, " = bottom right"
      write(17,*)npts, nvecs
      write (17,"(8f10.4)")  (((grid(i, j, k), i = 1, npts), j = 1, npts), k = 1, nvecs)
    end if
  end if
  q = char(ichar("1") + int(log10(nvecs*1.0)))
  outer_loop: do 
    do
      if(l_read_dat) write (iw, "(A,I"//q//",a)") "Pick an energy surface in the range 1 - ", nvecs, &
        ", 0 to stop"
      read(ir, '(a80)', iostat = j) line
      if (j /= 0) goto 50
      if (.not. l_read_dat .and. index(line, ": Pick an") == 0) cycle
      i = len_trim(line)
      write(line(20:),'(a,i'//q//',a)')": Pick an energy surface in the range 1 - ", nvecs, ", 0 to stop"
      call write_keystrokes(line, len_trim(line)) 
      read(line, *, end = 50, err = 50) isurf
      if (isurf > nvecs .or. isurf < 1) then
        if (isurf /= 0) then
          write(iw,'(a)') "To stop, enter ""0"", otherwise enter the surface to be drawn"
          cycle
        else
          exit outer_loop
        end if
      end if
      if(l_read_dat) write (iw, "(//,A,I"//q//",/)") " Drawing energy surface No.",isurf
!
!   Extract band from grid array
!
      npts2 = 2*npts - 1
      onebnd = 0.d0
      do i = npts, 1, -1
        do j = npts, 1, -1
          onebnd((2*i - 2)*npts2 + 2*j - 1) = grid(i, j, isurf)
        end do
      end do
!
!  Find turning points
!
      nmax = 0
      nmin = 0
      do i = 2, npts - 1
        do j = 2, npts - 1
          x = grid(i,j,isurf)
!
!  Work hard to avoid spurious maxima and minima
!
          if ( x > grid(i-1,j,isurf)   .and. x > grid(i+1,j,isurf) .and. &
          &    x > grid(i,j-1,isurf)   .and. x > grid(i,j+1,isurf)) then
          if ( x > grid(i-1,j+1,isurf) .and. x > grid(i+1,j-1,isurf) .and. &
             & x > grid(i-1,j-1,isurf) .and. x > grid(i+1,j+1,isurf)) then
!
!  Point x is near to a maximum
!
              x2 = 0.5d0*(grid(i-1,j,isurf) - grid(i+1,j,isurf)) / &
              & (grid(i-1,j,isurf) + grid(i+1,j,isurf) - 2.d0*grid(i,j,isurf))
              stepx = 1.d0 - (i-1+x2)/(npts-1)
              y2 = 0.5d0*(grid(i,j-1,isurf) - grid(i,j+1,isurf)) / &
              & (grid(i,j-1,isurf) + grid(i,j+1,isurf) - 2.d0*grid(i,j,isurf))
              stepy = 1.d0 - (j-1+x2)/(npts-1)
              nmax = Min(nmax + 1, 1000)
              maxima(1:3,nmax) = top_l(:) - stepx*(top_l(:) - top_r(:)) - stepy*(top_l(:) - bottom_l(:))
              maxima(4,nmax) = x + 0.5d0*(grid(i,j-1,isurf) - grid(i,j+1,isurf))*x2 &
              & + 0.5d0*(grid(i,j-1,isurf) + grid(i,j+1,isurf) - 2.d0*grid(i,j,isurf))*x2**2
            end if
          end if
   
          if ( x < grid(i-1,j,isurf)   .and. x < grid(i+1,j,isurf) .and. &
          &    x < grid(i,j-1,isurf)   .and. x < grid(i,j+1,isurf)) then
            if ( x < grid(i-1,j+1,isurf) .and. x < grid(i+1,j-1,isurf) .and. &
             &   x < grid(i-1,j-1,isurf) .and. x < grid(i+1,j+1,isurf)) then
!
!  Point x is near to a minimum
!
              x2 = 0.5d0*(grid(i-1,j,isurf) - grid(i+1,j,isurf)) / &
              & (grid(i-1,j,isurf) + grid(i+1,j,isurf) - 2.d0*grid(i,j,isurf))
              stepx = 1.d0 - (i-1+x2)/(npts-1)
              y2 = 0.5d0*(grid(i,j-1,isurf) - grid(i,j+1,isurf)) / &
              & (grid(i,j-1,isurf) + grid(i,j+1,isurf) - 2.d0*grid(i,j,isurf))
              stepy = 1.d0 - (j-1+x2)/(npts-1)
              nmin = Min(nmin + 1, 1000)
              minima(1:3,nmin) = top_l(:) - stepx*(top_l(:) - top_r(:)) - stepy*(top_l(:) - bottom_l(:)) 
              minima(4,nmin) = x + 0.5d0*(grid(i,j-1,isurf) - grid(i,j+1,isurf))*x2 &
              & + 0.5d0*(grid(i,j-1,isurf) + grid(i,j+1,isurf) - 2.d0*grid(i,j,isurf))*x2**2
              minima(5,nmin) = stepx
              minima(6,nmin) = stepy
              continue
            end if
          endif
        end do
      end do
!
!  Eliminate maxima that are near to each other.
!
      do
        l = 0
        do i = 2, nmax
          if(maxima(4,i) < -999.d0) cycle
          do j = 1, i - 1
            if(maxima(4,j) < -999.d0) cycle
            step = (maxima(1,i) -  maxima(1,j))**2 + &
                 & (maxima(2,i) -  maxima(2,j))**2 + &
                 & (maxima(3,i) -  maxima(3,j))**2
            if (step < 0.025d0) then
              l = k
              if (maxima(4,i) > maxima(4,j)) then
                maxima(4,j) = -1000.d0
              else
                maxima(4,i) = -1000.d0
              end if                   
            end if
          end do
        end do
        if (l == 0) exit
      end do
      j = 0
      do i = 1, nmax
        if (maxima(4,i) > -999.d0) then
          j = j + 1
          maxima(1:4,j) = maxima(1:4,i)
        end if
      end do
      nmax = j
!
!  Eliminate minima that are near to each other.
!
      do
        l = 0
        do i = 2, nmin
          if(minima(4,i) < -999.d0) cycle
          do j = 1, i - 1
            if(minima(4,j) < -999.d0) cycle
            step = (minima(1,i) -  minima(1,j))**2 + &
                 & (minima(2,i) -  minima(2,j))**2 + &
                 & (minima(3,i) -  minima(3,j))**2
            if (step < 0.025d0) then
              l = k
              if (minima(4,i) < minima(4,j)) then
                minima(4,j) = -1000.d0
              else
                minima(4,i) = -1000.d0
              end if                   
            end if
          end do
        end do
        if (l == 0) exit
      end do
      j = 0
      do i = 1, nmin
        if (minima(4,i) > -999.d0) then
          j = j + 1
          minima(:,j) = minima(:,i)
        end if
      end do
!
!  Sort into decreasing order
!
      nmin = j
      do i = 1, nmax
        do j = i + 1, nmax
          if(maxima(4,j) > maxima(4,i)) then
            temp_ima(:) = maxima(:,j)
            maxima(:,j) = maxima(:,i)
            maxima(:,i) = temp_ima(:)
          end if
        end do
      end do
!
! Sort into increasing order
!
      do i = 1, nmin
        do j = i + 1, nmin
          if(minima(4,j) < minima(4,i)) then
            temp_ima(:) = minima(:,j)
            minima(:,j) = minima(:,i)
            minima(:,i) = temp_ima(:)
          end if
        end do
      end do
      write(line,'(a,a,i2.2,a)')trim(data_set_name)," surface number ",isurf," notes.txt"
      inquire(unit=iw_new, opened=opend) 
      if (opend) close (iw_new)
      open(unit = iw_new, file = trim(line))
      if( nmax > 0 )then
        if(nmax == 1) then
          write(iw,'(/19x,a)')"Maximum in plot"
          write(iw_new,'(/19x,a)', iostat = i99)"Maximum in plot"
        else
          write(iw,'(/19x,a)')"Maxima in plot"
          write(iw_new,'(/19x,a)', iostat = i99)"Maxima in plot"
        end if
        write(iw,"(a)")"       kx          ky          kz           Value"
        write(iw_new,"(a)", iostat = i99)"       kx          ky          kz           Value"
        do i = 1,nmax
          write(iw,    "(3f12.4,f14.4)")maxima(1:4,i)
          write(iw_new,"(3f12.4,f14.4)", iostat = i99)maxima(1:4,i)
        end do
      end if
      if( nmin > 0 ) then
        if(nmin == 1) then
          write(iw,'(/19x,a)')"Minimum in plot"
          write(iw_new,'(/19x,a)', iostat = i99)"Minimum in plot"
        else
          write(iw,'(/19x,a)')"Minima in plot"
          write(iw_new,'(/19x,a)', iostat = i99)"Minima in plot"
        end if
        write(iw,"(a)")"       kx          ky          kz          Value"
        write(iw_new,"(a)", iostat = i99)"       kx          ky          kz          Value"
        do i = 1,nmin
          write(iw,"(3f12.4,f14.4,f18.4,f14.4)")minima(:4,i)
          write(iw_new,"(3f12.4,f14.4)", iostat = i99)minima(1:4,i)
        end do
      end if
!
!   Double up mesh.  Array was NPTS, increase it to 2*NPTS-1
!
      mpts = npts*2 - 1
      do i = 1, mpts, 2
        irow = (i-1) * mpts
        onebnd(irow+2) = 0.5d0 * (onebnd(irow+1)+onebnd(irow+3))
        onebnd(i*mpts-1) = 0.5d0 * (onebnd(i*mpts)+onebnd(i*mpts-2))
        do j = 3, mpts-4, 2
          onebnd(irow+j+1) = 0.0625d0 * &
               & (-onebnd(irow+j-2)-onebnd(irow+j+4) + &
               & 9.d0 * (onebnd(irow+j)+onebnd(irow+j+2)))
        end do
      end do
      npts2 = mpts * mpts
      do i = 1, mpts
        onebnd(i+mpts) = 0.5d0 * (onebnd(i)+onebnd(i+mpts*2))
        onebnd(i+npts2-2*mpts) = 0.5d0 * &
             & (onebnd(i+npts2-mpts)+onebnd(i+npts2-3*mpts))
        do j = 3, mpts-4, 2
          onebnd(i+j*mpts) = 0.0625d0 * &
               & (-onebnd(i+(j-1)*mpts)-onebnd(i+(j+3)*mpts) + &
               & 9.d0*(onebnd(i+ (j-1)*mpts)+onebnd(i+ (j+1)*mpts)))
        end do
      end do
!
!   Determine maximum and minimum contours, and step size.
!
      xmin = 1.d10
      xmax = -1.d10
      do i = 1, npts2
        xmin = Min (xmin, onebnd(i))
        xmax = Max (xmax, onebnd(i))
      end do
      range = xmax - xmin
!
!   Determine a suitable step-size
!
      i = Int(Log10(range) + 10.d0) - 10
      step = range*10.d0**(-i)
      i = i-1
      if (step > 5) then
        step = 2.d0*10.d0**i
      else if (step > 2.d0) then
        step = 1.d0*10.d0**i
      else
        step = 0.5d0*10.d0**i
      end if
      cntval = step * Int (xmin/step) - step
      if (first) then
! 
!  One-time setup
!
       first = .false.
      end if
!
!  INITIALIZE & SET-UP FOR NEXT PICTURE
!

      call graphics (0.0, 0.0, 6)
      call graphics (0.0, 0.0, 99)
      write(iw,"(/,3f12.5,a)")top_l, " = top left"
      write(iw,"(3f12.5,a)")top_r, " = top right"
      write(iw,"(3f12.5,a)")bottom_l, " = bottom left"
      write(iw,"(3f12.5,a)")bottom_r, " = bottom right"
      write(iw,"(3f12.5,a)")px, py, pz, " = center"
      write(iw_new,"(/,3f12.5,a)", iostat = i99)top_l, " = top left"
      write(iw_new,"(3f12.5,a)", iostat = i99)top_r, " = top right"
      write(iw_new,"(3f12.5,a)", iostat = i99)bottom_l, " = bottom left"
      write(iw_new,"(3f12.5,a)", iostat = i99)bottom_r, " = bottom right"
      write(iw_new,"(3f12.5,a)", iostat = i99)px, py, pz, " = center"

      if(phonon) then
        write(iw,"(/5x,a,f9.5,a)")"Contour interval = ",step," cm**(-1)"
        write(iw_new,"(/5x,a,f9.5,a)", iostat = i99)"Contour interval = ",step," cm**(-1)"
      else
        write(iw,"(/5x,a,f8.5,a)")"Contour interval = ",step," eV"
        write(iw_new,"(/5x,a,f8.5,a)", iostat = i99)"Contour interval = ",step," eV"
      end if
      write(iw,'(/5x,a,/5x,a,/)')"To exit, click on a point in the window ""Program BZ""", &
      "that is outside the contour map"
      call graphics (0.0, 0.0, 98)  ! Select black
      call graphics (0.0, 0.0, 2)
      call graphics (0.0, 2.0, 3)
      call graphics (2.0, 2.0, 3)
      call graphics (2.0, 0.0, 3)
      call graphics (0.0, 0.0, 3)
      write(line,'(a,i3)')"Band:",isurf
      call graphics (0.90, -0.05, 96) ! Write band number
      write(line,'(a,3f8.3,a)')"(",top_l,")"
      call graphics (0.0, -0.05, 96)  ! Write coordinate of top left
      write(line,'(a,3f8.3,a)')"(",top_r,")"
      call graphics (1.50, -0.05, 96) ! Write coordinate of top right
      write(line,'(a,3f8.3,a)')"(",bottom_l,")"
      call graphics (0.0, 2.02, 96)   ! Write coordinate of bottom left
      write(line,'(a,3f8.3,a)')"(",bottom_r,")"
      call graphics (1.50, 2.02, 96)  ! Write coordinate of bottom right
!
      do i = 1, 1000
        cntval = cntval + step
         if ( cntval > xmax) exit
!
!  Select new color
!
        call graphics(i*1.0, 0.0, 99)
        call cntour (onebnd, mpts, mpts, mpts, 0.d0, cntval, 1, 2, 3, x)
      end do
!
! Draw a box in black
!
      
      call graphics (0.0, 0.0, 98)  ! Select black
      call graphics (0.0, 0.0, 2)
      call graphics (0.0, 2.0, 3)
      call graphics (2.0, 2.0, 3)
      call graphics (2.0, 0.0, 3)
      call graphics (0.0, 0.0, 3)
      do 
        if (line(1:4) == "EXIT") exit
        call sleep(1)
      end do
!
!  Write raw data in a form suitable for use in Microsoft Excel
!
      write(graph_name,'(a,a,i2.2,a)')trim(data_set_name)," surface number ",isurf,".txt"
      write(iw,'(/5x,a)')"Writing suface to file """//trim(graph_name)//""""      
      inquire(unit=18, opened=opend) 
      if (opend) close (18)
      open(unit = 18, file = graph_name)
      write(18,'(8x,200f8.3)')(i*1.d0, i = 1,mpts)
      do i = 1, mpts
        write(18,'(201f8.3)')i*1.d0,(onebnd((i-1)*mpts + j), j = 1, mpts)
      end do
      call sleep(2)
    end do    
  end do outer_loop
50 continue
  if (ir == 99) then
    write(iw,'(/5x,a)')"Press ""RETURN"" to quit"
    read(iw,'(a)', iostat = i99) q
  end if
  call graphics(0.0, 0.0, 100)
  stop
end subroutine surfac
