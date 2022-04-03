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

  SUBROUTINE graphics(xx, yy, Action ) 
!
! Subroutine "graphics" looks after all graphical operations
!
! On input xx = x-coordinate of the point from the non-graphics side of the program 
!          yy = y-coordinate of the point from the non-graphics side of the program
!
! Action: An integer with value 1, 2, 3, 6, 99, 100, idicating which graphics task is to be done.
!
  USE IFQWIN   
  use common_common, only : line, xscale, yscale, xoffset, yoffset, top_left_x, top_left_y
  implicit none
  interface
  subroutine Details_of_cursor_point(unit, mouseevent, keystate, MouseXpos,MouseYpos)
      INTEGER unit
      INTEGER mouseevent
      INTEGER keystate
      INTEGER MouseXpos
      INTEGER MouseYpos
    end subroutine 
  end interface
  integer, intent (in) :: Action
  real, intent (in) :: xx, yy
  !
  ! Local
  !
 ! real :: xoffset,          &  ! Define the relative position on the monitor
 !         yoffset,          &  ! Define the relative position on the monitor
 !         xscale  =  450.0, &  ! Convert from non-graphics coordinates to graphics coordinates
 !         yscale  =  450.0, &  ! Convert from non-graphics coordinates to graphics coordinates
  real ::        x, y
  integer(4) ::   i, j, k, l, col(3), all_colors(50)
  INTEGER(2) :: color, pick_color
  TYPE (windowconfig) :: myscreen
  TYPE (qwinfo)      :: FrameSize
  TYPE (xycoord) :: xy
  character :: type_font*30, num1*1, num2*1
  logical :: status, l_font = .false.
!
! All colors used for lines, expressed as three fractions per color, blue,green,red
!
  real, dimension (3,50) :: colors = (/& 
& 0.0,0.0,1.0,  0.0,0.1,0.9,  0.0,0.2,0.8, 0.0,0.3,0.7, 0.0,0.4,0.6, &
& 0.0,0.5,0.5,  0.0,0.6,0.4,  0.0,0.7,0.3, 0.0,0.8,0.2, 0.0,0.9,0.1, &
& 0.0,1.0,0.0,  0.0,1.0,0.1,  0.0,1.0,0.2, 0.0,1.0,0.3, 0.0,1.0,0.4, &
& 0.0,1.0,0.5,  0.0,1.0,0.6,  0.0,1.0,0.7, 0.0,1.0,0.8, 0.0,1.0,0.9, &
& 1.0,0.0,1.0,  1.0,0.1,0.9,  1.0,0.2,0.8, 1.0,0.3,0.7, 1.0,0.4,0.6, &
& 1.0,0.5,0.5,  1.0,0.6,0.4,  1.0,0.7,0.3, 1.0,0.8,0.2, 1.0,0.9,0.1, &
& 1.0,1.0,0.0,  1.0,1.0,0.1,  1.0,1.0,0.2, 1.0,1.0,0.3, 1.0,1.0,0.4, &
& 1.0,1.0,0.5,  1.0,1.0,0.6,  1.0,1.0,0.7, 1.0,1.0,0.8, 1.0,1.0,0.9, &
& 1.0,0.0,1.0,  1.0,0.1,1.0,  1.0,0.2,1.0, 1.0,0.3,1.0, 1.0,0.4,1.0, &
& 1.0,0.5,1.0,  1.0,0.6,1.0,  1.0,0.7,1.0, 1.0,0.8,1.0, 1.0,0.9,1.0  &
 /)
save
  select case (Action)
!
!
  case (1) !  Run one time only, to set up graphics initialization
!
!  Create graphics window frame to fill the whole screen
!
    Status=GETWSIZEQQ(QWIN$FRAMEWINDOW,QWIN$SIZEMAX,FrameSize)
    Status=SETWSIZEQQ(QWIN$FRAMEWINDOW,FrameSize)
!
!  Create  child window to hold text
!

    OPEN (UNIT = 6, FILE = 'USER', IOFOCUS = .TRUE.)
    myscreen.numxpixels = int2(0.40*FrameSize.W)  ! Width of the text window
    myscreen.numypixels = int2(0.75*FrameSize.H)  ! Height of the text window
    myscreen.numtextcols=-1
    myscreen.numtextrows=-1
    myscreen.numcolors=-1
    myscreen.fontsize=-1
    myscreen.title = "Text "C ! Title of the graphics screen
    status = SETWINDOWCONFIG(myscreen) ! Set configuration of the graphics window
    Framesize.X = 0
    Status=SETWSIZEQQ(6,FrameSize)
!
!  Background colors are set in six bytes: BBGGRR. White is FFFFFF and black is 000000.
!  Colow child window white
!
    i=setbkcolorrgb(#FFFFFF)   !! white background
    i=settextcolorrgb(#000000) !! black text
    call clearscreen($GCLEARSCREEN)
!
!  Create child window to hold graphics
!
    OPEN (UNIT = 10, FILE = 'USER', IOFOCUS = .TRUE.)

    myscreen.numxpixels = int2(0.60*FrameSize.W)   ! Width of the graphics window
    myscreen.numypixels = int2(0.90*FrameSize.H)  ! Height of the graphics window
    myscreen.numtextcols=-1
    myscreen.numtextrows=-1
    myscreen.numcolors=-1
    myscreen.fontsize=-1
    myscreen.title = "Program BZ "C ! Title of the graphics screen
    status = SETWINDOWCONFIG(myscreen) ! Set configuration of the graphics window
!
!  Position child window for graphics
!
    Framesize.TYPE = QWIN$SET
    Framesize.X = 0.053*Framesize.W     ! Position graphics window to the right of the screen
    top_left_x  = Framesize.X
    top_left_y  = Framesize.Y
    Status=SETWSIZEQQ(10,FrameSize)
    i = SETBKCOLORRGB(Z'FFFFFF')
    call clearscreen($GCLEARSCREEN)
    Status = REGISTERMOUSEEVENT (10, MOUSE$LBUTTONDOWN, Details_of_cursor_point)
    xoffset = 70.0
    yoffset = 40.0
    xscale  =  450.0*framesize.W/1920.0 ! Convert from non-graphics coordinates to graphics coordinates
    yscale  =  450.0*framesize.H/1133.0 ! Convert from non-graphics coordinates to graphics coordinates

    pick_color = 0
!
!  Set colors used in lines as four-byte integers. 
!
!  Byte 1: Red
!  Byte 2: Green
!  Byte 3: Blue
!  Byte 4: not used
!
    do i = 1,50
      do l = 1,3
      j = nint(colors(l,i)*255)
      col(l) = j
      end do
      all_colors(i) = col(1) + 256*(col(2) + 256*col(3))
    end do      
!
! Construct part of format for text
!
    i = nint((18*framesize.H)/1133.0)
    num1 = char(ichar("1") + i/10)
    j = nint((10*framesize.W)/1920.0)
    num2 = char(ichar("1") + j/10)
!
! Set up fonts
!
    i = INITIALIZEFONTS()
    return
!
!
  case (2) ! Move to new point, don't draw a line
    i = SETACTIVEQQ(10)
    i = FOCUSQQ(10)
    x = xx*xscale + xoffset
    y = yy*yscale + yoffset
    CALL MOVETO(int2(x), int2(y), xy ) 
    return
!
!
  case (3) ! Move to new point drawing a line from the old point
    x = xx*xscale + xoffset
    y = yy*yscale + yoffset
    i = LINETO(int2(x), int2(y)) 
    return
!
!
  case (6)  ! Run each time a new picture is drawn, initialize a picture.  
    i = SETACTIVEQQ(10)
    i = FOCUSQQ(10)
    CALL SETVIEWORG( INT2(0),  INT2( 0 ), xy )
    call clearscreen($GCLEARSCREEN)
    i = SETCOLORRGB(Z'000000')
    pick_color = 0
    return
!
!
  case (96)  ! Write text to graph in Arial and Symbol
    i = SETCOLORRGB(Z'000000')
    write(type_font,'(a,i'//num1//',a,i'//num2//',a)')"t'Arial'h" ,i, "w", j, "vb'"
    if (len_trim(line) == 3) then
      if (line(3:3) == "G") write(type_font,'(a,i'//num1//',a,i'//num2//',a)')"t'Symbol'h" ,i, "w", j, "vb'"
      if (line(3:3) == "G" .or. line(3:3) == "R")line(3:3) = " "
    end if
    i = SETFONT(trim(type_font))
    x = xx*xscale + xoffset
    y = yy*yscale + yoffset
    CALL MOVETO(int2(x), int2(y), xy )
    call OUTGTEXT(trim(line))
    return
!
!
  case (97)  ! Set focus to channel 6

    i = FOCUSQQ(6)
    return
!
!
  case (98)  ! Select black

    i = SETACTIVEQQ(10)
    i = FOCUSQQ(10)
    i = SETCOLORRGB(Z'000000')
    return
!
!
    case (99)  ! Select color
    pick_color = pick_color + 1
    if (pick_color == 50) pick_color = 1 
    i = SETCOLORRGB(all_colors(pick_color))
    return
!
!
  case (100)  ! Run one time only, to gracefully terminate the job
    i = SETEXITQQ(QWIN$EXITNOPERSIST)
    return
  end select
END SUBROUTINE graphics