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

subroutine paths()
!-----------------------------------------------
!   Follow a reaction path in which the values of the coordinate are supplied by the user
!-----------------------------------------------
    use molkst_C, only : iflepo, numat, keywrd, nvar, tleft, &
          & time0, escf, norbs, moperr, line, nl_atoms
    use maps_C, only : lparam, react, latom, rxn_coord
    use common_arrays_C, only : geo, xparam, na, hesinv, l_atom, nat, coord
    use ef_C, only : alparm, x0, x1, x2, iloop
    use chanel_C, only : iw0, iw, ires, restart_fn, ixyz, xyz_fn
    use elemts_C, only : elemnt
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
    implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    integer :: maxcyc, i, j, ii, lpr, npts, iw00, percent = 0, ipdb = 14, imodel = 0
    double precision :: x3, c3, cc1, cc2, cb1, cb2, delf0, delf1, aconst, bconst, cconst, c1
    logical :: lef, debug
    double precision, external :: reada, seconds
!-----------------------------------------------

!***********************************************************************
!
!   PATH FOLLOWS A REACTION COORDINATE.   THE REACTION COORDINATE IS ON
!        ATOM LATOM, AND IS A DISTANCE IF LPARAM=1,
!                           AN ANGLE   IF LPARAM=2,
!                           AN DIHEDRALIF LPARAM=3.
!
!***********************************************************************
    iloop = 1! Counter for reaction coordinate steps
    iw00 = iw0
    debug = (index(keywrd, " DEBUG") /= 0)
    if (index(keywrd, " PDBOUT") /= 0) &
      open(unit = ipdb, file = xyz_fn(:len_trim(xyz_fn) - 3)//"pdb")
    open(unit=ixyz, file=xyz_fn)
    if (.not. debug) iw0 = -1
    if (allocated(hesinv)) deallocate(hesinv)
    allocate (hesinv(nvar*nvar))
    hesinv = 0.0D00

! Set a maximum of reaction coordinate steps (BIGCYCLES)
    maxcyc = 100000
    if (index(keywrd,' BIGCYCLES') /= 0) &
          maxcyc = nint(reada(keywrd,index(keywrd,' BIGCYCLES')))

! Use LBFGS optimizer as default, EF on demand
    lef = (index(keywrd,' LBFGS') == 0 .and. nvar <101)

    if (allocated(alparm)) deallocate(alparm)
    allocate(alparm(3, nvar))
    alparm = 0.d0

! Initialize if RESTART
    CheckRestart: if (lef) then
        write (iw, '(''  ABOUT TO ENTER EF FROM PATH'')')
        if (index(keywrd,'RESTAR') /= 0) then
            open(unit=ires, file=restart_fn, status='UNKNOWN', &
                  & form='UNFORMATTED', position='asis')
            rewind ires
            read (ires, end=120, err=120) i, j
            if (norbs /= j .or. numat /= i) then
                call mopend("Restart file read in does not match current data set")
                goto 99
            end if
            read (ires, end=120, err=120) ((alparm(j,i),j=1,3),i=1,nvar)
            read (ires, end=120, err=120) iloop, x0, x1, x2
            close (ires)
        end if

    else
        write (iw, '(''  ABOUT TO ENTER LBFGS FROM PATH'')')
        if (index(keywrd,'RESTAR') /= 0) then
          call mopend(" RESTART is not supported in a path calculation unless EF is used")
          write(iw,'(10x, "(Edit the previous PATH output to make a new job that"'// &
          '"starts the path where the previous path stopped)")')
          return
        end if

    end if CheckRestart

! Conversion factor for angles
    if (lparam /= 1 .and. na(latom) /= 0) then
        c1 = 57.29577951308232D0
    else
        c1 = 1.D0
    end if

    FirstReactionStep: if (iloop <= 1) then
        time0 = seconds(1)
        if (maxcyc == 0) tleft = -100.D0
        if (lef) then
            call ef (xparam, escf)
        else
            call lbfgs (xparam, escf)
        end if
        i = index(keywrd,'RESTAR')
        if (i /= 0) keywrd(i:i+6) = " "
        if (iw00 > -1) then
            write (line, '('' :'',F16.5,F16.6)') geo(lparam,latom)*c1, escf
            call to_screen(line)
        end if
        if (moperr) goto 99
        if (iflepo == (-1)) goto 99
        write (iw, '(''  OPTIMIZED VALUES OF PARAMETERS, INITIAL POINT'')')
        if (index(keywrd, " PDBOUT") /= 0) then
          imodel = imodel + 1
          write(ipdb,'(a,i7)')"MODEL",imodel
          call pdbout(ipdb)
          write(ipdb,'(a)')"ENDMDL"
        end if
        call writmo
        time0 = seconds(1)
    end if FirstReactionStep

    NextReactionStep: if (iloop <= 2) then
        geo(lparam,latom) = react(2)

        SecondReactionStep: if (iloop == 1) then
            x0 = react(1)
            x1 = x0
            x2 = react(2)
            if (x2 < (-100.D0)) call mopend ('Error in PATHS')
            if (x2 < (-100.D0)) goto 99
            alparm(2,:nvar) = xparam(:nvar)
            alparm(1,:nvar) = xparam(:nvar)
            iloop = 2
        end if SecondReactionStep

        if (maxcyc == 1) tleft = -100.D0
        if (lef) then
            call ef (xparam, escf)
        else
            call lbfgs (xparam, escf)
        end if
        if (iw00 > -1) then
            write (line, '('' :'',F16.5,F16.6)') geo(lparam,latom)*c1, escf
            call to_screen(line)
        end if
        if (iflepo == (-1)) goto 99
        rxn_coord = react(2)
        if (lparam > 1 .and. na(latom) > 0) then
            rxn_coord = rxn_coord*57.29577951308232D0
            write (iw, &
                  '(1X,16(''*****''),2/17X,''REACTION COORDINATE = ''       ,F12.4,2X,A10,1&
                  &9X,2/1X,16(''*****''))') rxn_coord, 'DEGREES   '
        else
            write (iw, &
                  '(1X,16(''*****''),2/17X,''REACTION COORDINATE = ''       ,F12.4,2X,A10,1&
                  &9X,2/1X,16(''*****''))') rxn_coord, 'ANGSTROMS   '
        end if

        call writmo
!
!  Write out "xyz" file
!
        write(ixyz,"(i6,a)") nl_atoms," "
        write(ixyz,*)"PATH "
        do i = 1, numat
          if (l_atom(i)) write(ixyz,"(3x,a2,3f15.5)")elemnt(nat(i)), (coord(j,i),j=1,3)
        end do
        if (index(keywrd, " PDBOUT") /= 0) then
          imodel = imodel + 1
          write(ipdb,'(a,i7)')"MODEL",imodel
          call pdbout(ipdb)
          write(ipdb,'(a)')"ENDMDL"
        end if
        time0 = seconds(1)
        alparm(3,:nvar) = xparam(:nvar)
        if (iloop == 2) iloop = 3
    end if NextReactionStep

! Find number of reaction path steps
    lpr = iloop
    do npts = 1,10000
        if (react(npts) < -100.D0) exit
    end do

    MainLoop: do ii = lpr, npts - 1
        iloop = ii
        if (iloop - lpr > maxcyc - 3) then
            tleft = -100.D0
        end if
        rxn_coord = react(iloop)
        if (lparam > 1 .and. na(latom) > 0) then
            rxn_coord = rxn_coord*57.29577951308232D0
            write (iw, &
                  '(1X,16(''*****''),2/17X,''REACTION COORDINATE = '',F12.4,2X,A10,19X,2/1X,16(''*****''))') &
                  rxn_coord, 'DEGREES   '
        else
            write (iw, &
                  '(1X,16(''*****''),2/17X,''REACTION COORDINATE = '',F12.4,2X,A10,19X,2/1X,16(''*****''))') &
                  rxn_coord, 'ANGSTROMS   '
        end if

! Determine type of interpolation
        x3 = react(iloop)
        c3 = (x0*x0 - x1*x1)*(x1 - x2) - (x1*x1 - x2*x2)*(x0 - x1)
        if (abs(c3) < 1.D-8) then
 ! Linear interpolation
            cc1 = 0.D0
            cc2 = 0.D0
        else
 ! Quadratic interpolation
            cc1 = (x1 - x2)/c3
            cc2 = (x0 - x1)/c3
        end if
        cb1 = 1.D0/(x1 - x2)
        cb2 = (x1*x1 - x2*x2)*cb1

! Calculate the interpolated coordinates
        do i = 1, nvar
            delf0 = alparm(1,i) - alparm(2,i)
            delf1 = alparm(2,i) - alparm(3,i)
            aconst = cc1*delf0 - cc2*delf1
            bconst = cb1*delf1 - aconst*cb2
            cconst = alparm(3,i) - bconst*x2 - aconst*x2**2
            xparam(i) = cconst + bconst*x3 + aconst*x3**2
            alparm(1,i) = alparm(2,i)
            alparm(2,i) = alparm(3,i)
        end do

! Check that the guessed geometry is not too absurd
        CheckGeo: do i = 1, nvar
            if (abs(xparam(i)-alparm(3,i)) > 0.2D0) then
     ! Too large deviation, resort to old geometry
                write (iw, &
                      '('' GEOMETRY TOO UNSTABLE FOR EXTRAPOLATION TO BE USED'',/,  &
                      & " - THE LAST GEOMETRY IS BEING USED TO START THE NEXT",'' CALCULATION'')')
                xparam(:nvar) = alparm(3,:nvar)
                exit CheckGeo
            end if
        end do CheckGeo

        x0 = x1
        x1 = x2
        x2 = x3
        geo(lparam,latom) = react(iloop)
       if (lef) then
            call ef (xparam, escf)
        else
            call lbfgs (xparam, escf)
        end if
        if (iw00 > -1) then
            i = nint((100.0*iloop)/npts)
            if (i /= percent) then
                percent = i
                write(line,"(i4,a)")percent, "% of Reaction Coordinate done"
                call to_screen(line)
            end if
            write (line, '('' :'',F16.5,F16.6)') geo(lparam,latom)*c1, escf
            call to_screen(line)
        end if
        if (iflepo == (-1)) goto 99
        if (index(keywrd, " PDBOUT") /= 0) then
          imodel = imodel + 1
          write(ipdb,'(a,i7)')"MODEL",imodel
          call pdbout(ipdb)
          write(ipdb,'(a)')"ENDMDL"
        end if
        call writmo
        time0 = seconds(1)
        alparm(3,:nvar) = xparam(:nvar)
    end do MainLoop

    iw0 = iw00
    goto 99
120 continue
    write (iw, '(2/10X,''Restart file is corrupt!'')')
    call mopend ('Restart file is corrupt!')
99  continue
    if (allocated(alparm)) deallocate(alparm)
    return
end subroutine paths
