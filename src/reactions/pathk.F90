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

      subroutine pathk()
!-----------------------------------------------
!   Follow a reaction path in which the step-size is a constant,
!   Keywords STEP and POINTS are used
!-----------------------------------------------
      use molkst_C, only : iflepo, numat, keywrd, nvar, tleft, escf, &
      line, norbs, numcal, id, natoms, nl_atoms, time0, ncomments
      use maps_C, only : rxn_coord, lparam, react, latom, kloop
      use common_arrays_C, only : geo, xparam, profil, na, l_atom, coord, nat, &
        loc, grad, all_comments, p, q, labels
      use chanel_C, only : iw0, iw, ires, iarc, restart_fn, archive_fn, &
        ixyz, xyz_fn
      use elemts_C, only : elemnt
      USE parameters_C, only : tore
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(20) :: mdfp
      integer :: npts, maxcyc, i, lloop, iloop, l, m, k, iw00, percent, ipdb, &
        imodel, ixyz1
      double precision :: step, degree, c1, cputot, cpu1, cpu2, cpu3, stepc1, factor, dip, &
        dipvec(3), xdfp(20),  gd(3*numat), xlast(3*numat)
      logical :: use_lbfgs, opend, scale, debug, l_dipole
      character :: num1*1, num2*1
      double precision, external :: dipole, reada, seconds
!-----------------------------------------------
      imodel = 0
      ipdb = 14
      ixyz1 = 0
      percent = 0
      use_lbfgs = (index(keywrd,' LBFGS') /=0 .or. nvar > 100 .and. index(keywrd,' EF') == 0)
      step = reada(keywrd,index(keywrd,'STEP') + 5)
      debug = (index(keywrd, " DEBUG") /= 0)
      l_dipole = (index(keywrd, " DIPOLE") /= 0)
      npts = nint(reada(keywrd,index(keywrd,'POINT') + 6))
!
!  THE SMALLEST VALUE IN THE PATH IS
!      REACT(1) DEGREE OR GEO(LPARAM,LATOM) RADIANS
!
      degree = 57.29577951308232D0
      if (lparam /= 1 .and. na(latom) /= 0) then
        step = step/degree
        c1 = degree
      else
        c1 = 1.D0
      end if
      if (index(keywrd, " PDBOUT") /= 0) &
      open(unit = ipdb, file = xyz_fn(:len_trim(xyz_fn) - 3)//"pdb")
      if (index(keywrd, " MINI") /= 0) then
        loc = 0
        nvar = 0
        do i = 1, natoms
          do l = 1, 3
            if (lparam /= l .or. latom /= i) then
              nvar = nvar + 1
              loc(1,nvar) = i
              loc(2,nvar) = l
              xparam(nvar) = geo(l,i)
            end if
          end do
        end do
        if (allocated(grad)) deallocate(grad)
        allocate(grad(nvar))
      end if
      open(unit=ixyz, file=xyz_fn)
      if (l_dipole) then
      ixyz1 = ixyz + 1
      open(unit=ixyz1, file=xyz_fn(:len_trim(xyz_fn) - 4)//" for dipole.xyz")
    end if
!
      kloop = 1
      maxcyc = 100000
      if (index(keywrd,' BIGCYCLES') /= 0) maxcyc = nint(reada(keywrd,index(&
        keywrd,' BIGCYCLES')))
      cputot = 0.0D0
      rxn_coord = geo(lparam,latom)
      if (allocated(profil)) deallocate(profil)
      if (allocated(react)) deallocate(react)
      allocate(profil(npts + 1), react(npts + 1))
      profil = 0.D0
      react(1) = geo(lparam,latom)
      if (use_lbfgs) then
        write (iw, '(/10x,''ABOUT TO ENTER L-BFGS FROM PATHK'')')
        if (index(keywrd,'RESTAR') /= 0) then
          mdfp(9) = 0 !   This section is almost certainly faulty!
          gd = 0.d0
          xlast = 0.d0
          xdfp = 0.d0
          call dfpsav (cputot, xparam, gd, xlast, escf, mdfp, xdfp)
          write (iw, '(2/10X,'' RESTARTING AT POINT'',I3)') kloop
        end if
      else
        write (iw, '(''  ABOUT TO ENTER EF FROM PATHK'')')
        if (index(keywrd,'RESTAR') /= 0) then
          open(unit=ires, file=restart_fn, status='UNKNOWN', form=&
            'UNFORMATTED', position='asis')
          rewind ires
          read (ires, end=60, err=60) i, l
          if (norbs /= l .or. numat /= i) then
              call mopend("Restart file read in does not match current data set")
              goto 99
          end if
          read (ires, err=60) kloop
          read (ires, err=60) rxn_coord
          read (ires, err=60) (profil(i),i=1,kloop)
          close (ires)
          write (iw, '(2/10X,'' RESTARTING AT POINT'',I3)') kloop
        end if
      end if
!
      geo(lparam,latom) = rxn_coord
      lloop = kloop
      scale = .false.
      if (id == 1) then
        do i = 1, natoms
          if (na(i) /= 0) exit
        end do
        scale = (i > natoms .and.latom == natoms)
      end if
      iw00 = iw0
      if (.not. debug) iw0 = -1
      if (index(keywrd, " HTML") /= 0) then
        if (index(keywrd,' DIPOLE') /= 0) call write_path_html(2)
        call write_path_html(1)
      end if
      do iloop = kloop, npts
        if (iloop - lloop >= maxcyc) tleft = -100.D0
        time0 = seconds(1)
        cpu1 = seconds(2)
        if (iloop > 1 .and. scale) then
          factor = geo(lparam,latom)/rxn_coord
          do i = 1, latom - 1
            xparam((i - 1)*3 +lparam) = xparam((i - 1)*3 +lparam)*factor
          end do
        end if
        rxn_coord = geo(lparam,latom)
        numcal = numcal + 1
        if (use_lbfgs) then
          write(iw,'(/10x,a)')"Geometry optimization using L-BFGS"
          if (iloop == 1 .and.iw00 > -1) then
            if (.not. debug) iw0 = 0 ! Temporarily allow writing to screen
            call to_screen(" Geometry optimization using L-BFGS")
            if (.not. debug) iw0 = -1
          end if
          call lbfgs (xparam, escf)
        else
          write(iw,'(/10x,a)')"Geometry optimization using EF"
          if (iloop == 1 .and.iw00 > -1) then
            if (.not. debug) iw0 = 0 ! Temporarily allow writing to screen
            call to_screen(" Geometry optimization using EF")
            if (.not. debug) iw0 = -1
          end if
          call ef (xparam, escf)
        end if
        i = index(keywrd,'RESTAR')
        if (i /= 0) keywrd(i:i+6) = " "
        i = index(keywrd,'OLDENS')
        if (i /= 0) keywrd(i:i+5) = " "
        if (iflepo == (-1) .or. tleft < 0.d0 ) goto 99
        kloop = kloop + 1
        cpu2 = seconds(2)
        cpu3 = cpu2 - cpu1
        cputot = cputot + cpu3
        profil(iloop) = escf
        write (iw, '(/''          VARIABLE        FUNCTION'')')
        write (iw, '('' :'',F16.5,F16.6)') geo(lparam,latom)*c1, escf
        if (iw00 > -1) then
          if (.not. debug) iw0 = 0 ! Temporarily allow writing to screen
          i = nint((100.0*iloop)/npts)
          if (i /= percent) then
            percent = i
            write(line,"('' :'',F16.4,F16.4,i9,a)")geo(lparam,latom)*c1, escf, &
              percent, "% of Reaction Coordinate done"
          else
            write(line,"('' :'',F16.4,F16.4,i9,a)")geo(lparam,latom)*c1, escf
          end if
          call to_screen(line)
        end if
        call to_screen("To_file: Reaction path")
        call geout (iw)
        if (index(keywrd, " PRTXYZ") /= 0) then
          write (iw, '(29X,''CARTESIAN COORDINATES '',/)')
          l = 0 
          do i = 1, natoms
            if (labels(i) == 99 .or. labels(i) == 107) cycle
            l = l + 1
            write (iw, '(I4,3X,A2,3x, 3F16.9)') l, elemnt(labels(i)), (coord(k,l),k=1,3)
          end do
        end if
        if (index(keywrd, " PDBOUT") /= 0) then
          imodel = imodel + 1
          write(line,'(F13.5)') escf
          do i = 1, 12
            if (line(i:i) /= " ") exit
          end do
          write(ipdb,'(a, i9, 2x, a)')"MODEL",imodel, trim(line(i:))
          if (ncomments == 0) ncomments = 1
          write(all_comments(1),'(a,f12.3,a)')" HEADER Heat of Formation =", escf, " Kcal/mol"
          call pdbout(ipdb)
          write(ipdb,'(a)')"ENDMDL"
        end if
        geo(lparam,latom) = geo(lparam,latom) + step
!
!  Write out "xyz" file
!
        write(ixyz,"(i6,a)") nl_atoms," "
        num1 = char(ichar("1") + int(log10(iloop*1.01)))
        factor = abs(escf)
        num2 = char(max(ichar("0"), ichar("0") + min(9, int(log10(factor)))))
        write(ixyz,'(a, i'//num1//', a, f1'//num2//'.5, a)')"Profile.", iloop, &
        " HEAT OF FORMATION =", escf, " KCAL "
        do i = 1, numat
          if (l_atom(i)) write(ixyz,"(3x,a2,3f15.5)")elemnt(nat(i)), (coord(l,i),l=1,3)
        end do
        if (l_dipole) then
          call chrge (p, q)
          q(:numat) = tore(labels(:numat)) - q(:numat)
          dip = dipole(p, coord, dipvec,0)
          write(ixyz1,"(i6,a)") nl_atoms," "
          write(ixyz1,'(a, i'//num1//', a, f1'//num2//'.5, a)')"Profile.", iloop, &
          " DIPOLE =", dip, " DEBYE"
          do i = 1, numat
            if (l_atom(i)) write(ixyz1,"(3x,a2,3f15.5)")elemnt(nat(i)), (coord(l,i),l=1,3)
          end do
        end if
        if (.not. debug) iw0 = -1
      end do
      if (cputot > 1.d7) cputot = cputot - 1.d7
      react(1) = react(1)*c1
      stepc1 = step*c1
      do i = 2, npts
        react(i) = react(i-1) + stepc1
      end do
      write (iw, &
      '(/16X,''POINTS ON REACTION PATH '',/16X,''AND CORRESPONDING HEATS'',2/)')
      inquire(unit=iarc, opened=opend)
      if (opend) close(unit=iarc, status='KEEP')
      open(unit=iarc, file=archive_fn, status='UNKNOWN', position=&
        'asis')
      write (iarc, 30)
      call wrttxt (iarc)
   30 format(' ARCHIVE FILE FOR PATH CALCULATION'/,&
        'A PROFILE OF COORDINATES - HEATS'/)
      write (iarc, '(/'' TOTAL JOB TIME : '',F10.3/)') cputot
!
      l = npts/8
      m = npts - l*8
      if (l >= 1) then
        do k = 0, l - 1
          write (iw, '(9F17.8)') (react(i),i=k*8 + 1,k*8 + 8)
          write (iw, '(9F17.8,/)') (profil(i),i=k*8 + 1,k*8 + 8)
        end do
      end if
      if (m > 0) then
        write (iw, '(9F17.8)') (react(i),i=l*8 + 1,l*8 + m)
        write (iw, '(9F17.8,/)') (profil(i),i=l*8 + 1,l*8 + m)
      end if
      do i = 1, npts
        write(iarc,'(2f17.8)')react(i), profil(i)
      end do
      if (.not. debug) iw0 = iw00
      goto 99
   60 continue
      write (iw, '(A,I3,A)') ' ERROR DETECTED DURING READ FROM CHANNEL', ires, &
        ' IN SUBROUTINE PATHK'
      call mopend ('ERROR DETECTED DURING READ IN SUBROUTINE PATHK')
  99  if (allocated(profil)) deallocate(profil)
      if (allocated(react)) deallocate(react)
      return
  end subroutine pathk
  subroutine write_path_html(mode)
    use chanel_C, only: input_fn
    use molkst_C, only : line, koment, escf, title, backslash, keywrd, l_normal_html
    implicit none
    integer, intent (in) :: mode
    logical :: exists, l_pdb
    integer :: iprt=27, i, j
    double precision :: store_escf
    character :: suffix*3
    l_normal_html = .false.
    if (mode == 1) then
      line = input_fn(:len_trim(input_fn) - 4)//"html"
    else
      line = input_fn(:len_trim(input_fn) - 5)//" for dipole.html"
    end if
    open(unit=iprt, file=trim(line))
    write(iprt,"(a)")"<!DOCTYPE html> "
    write(iprt,"(a)")"<html>"
    write(iprt,"(a)")"<head>"
    write(iprt,"(a)")"<!--"
    write(iprt,"(a)")""
    write(iprt,"(a)")"Template for JSmol reaction path: JSmol window plus HoF energy plot."
    write(iprt,"(a)")""
    write(iprt,"(a)")"-->"
    write(iprt,"(a)")""
    write(iprt,"(a)")"<style>"
    write(iprt,"(a)")"@media print { .noprint {display:none} .printonly {display:block}}"
    write(iprt,"(a)")"@media screen {.noprint {display:block} .printonly {display:none}}"
    write(iprt,"(a)")"</style>"
    write(iprt,"(a)")"<script type=""text/javascript"" src=""../jsmol/JSmol.min.js""></script>"
    write(iprt,"(a)")"<script type=""text/javascript"" src=""../jsmol/js/Jmol2.js""></script>"
    write(iprt,"(a)")"<script type=""text/javascript"" src=""../jsmol/flot/jquery.flot2.js""></script>"
    write(iprt,"(a)")""
    write(iprt,"(a)")"<script type=""text/javascript"">"
    write(iprt,"(a)")"function roundoff(x,ndec){"
    write(iprt,"(a)")"//round x to ndec decimal places"
    write(iprt,"(a)")"if(x==0)return 0"
    write(iprt,"(a)")"if(ndec==0)return Math.round(x)"
    write(iprt,"(a)")"var neg=(x<0?""-"":"""")"
    write(iprt,"(a)")"var xs=Math.abs(x)+"""""
    write(iprt,"(a)")"var i=(xs.indexOf(""E"") & xs.indexOf(""e""))"
    write(iprt,"(a)")"i=xs.indexOf(""."")"
    write(iprt,"(a)")"if (i<0) {xs=xs+""."";i=xs.indexOf(""."");}"
    write(iprt,"(a)")"xs=xs+""000000000"""
    write(iprt,"(a)")"var s="".""+xs.substring(i+1+ndec,xs.length)"
    write(iprt,"(a)")"xs=xs.substring(0,i)+xs.substring(i+1,i+1+ndec)"
    write(iprt,"(a)")"var add1=(xs.charAt(0)==""0"")"
    write(iprt,"(a)")"if(add1)xs=""1""+xs"
    write(iprt,"(a)")"xs=eval(xs)+Math.round(eval(s))+"""""
    write(iprt,"(a)")"if(add1)xs=xs.substring(1,xs.length)"
    write(iprt,"(a)")"xs=xs.substring(0,xs.length-ndec)+"".""+xs.substring(xs.length-ndec,xs.length)"
    write(iprt,"(a)")"if(xs.substring(0,1)==""."")xs=""0""+xs"
    write(iprt,"(a)")"return neg+xs"
    write(iprt,"(a)")"}"
    write(iprt,"(a)")""
    suffix = "xyz"
    l_pdb = (index(keywrd, " PDBOUT") /= 0)
    if (l_pdb) suffix = "pdb"
    if (mode == 1) then
      line = input_fn(:len_trim(input_fn) - 4)//suffix
    else
      line = input_fn(:len_trim(input_fn) - 5)//" for dipole."//suffix
    end if
    do i = len_trim(line), 1, -1
      if (line(i:i) == "/" .or. line(i:i) == backslash) exit
    end do
    line = line(i+1:)
    write(iprt,"(a)")"modelFile = """//backslash//"'"//trim(line)//backslash//"'; center {visible}; zoom 0;"""
    write(iprt,"(a)")"var appletPrintable = (navigator.appName != ""Netscape"")"
    write(iprt,"(a)")""
    write(iprt,"(a)")"// The following code is black magic - don't mess with it!"
    write(iprt,"(a)")""
    write(iprt,"(a)")"function setImage() {if (appletPrintable)return; var image = jmolGetPropertyAsString(""image"")"
    write(iprt,"(a)")"var html = '<img src=""data:image/jpeg;base64,'+image+'"" />'; "
    write(iprt,"(a)")"document.getElementById(""imagediv"").innerHTML = html }"
    write(iprt,"(a)")""
    write(iprt,"(a)")"function doPrintAll() { setImage(); window.print()}"
    write(iprt,"(a)")""
    write(iprt,"(a)")"$(function () { if (!appletPrintable)$(""#appletdiv"").addClass(""noprint"")})"
    write(iprt,"(a)")""
    write(iprt,"(a)")"var theplot = null   "
    write(iprt,"(a)")""
    write(iprt,"(a)")"var data = []"
    write(iprt,"(a)")""
    write(iprt,"(a)")"function plotEnergies(a,b,c,d,e) {"
    write(iprt,"(a)")"if (c == ""zapped"")return"
    write(iprt,"(a)")""
    write(iprt,"(a)")"setImage()"
    write(iprt,"(a)")""
    write(iprt,"(a)")"data = []"
    write(iprt,"(a)")"var A = []"
    write(iprt,"(a)")"var nplots = 1"
    write(iprt,"(a)")"Info = jmolGetPropertyAsArray(""auxiliaryInfo.models"");  "
    write(iprt,"(a)")"var modelCount = Info.length"
    write(iprt,"(a)")" "
    write(iprt,"(a)")"// Gather the data we want for each model."
    write(iprt,"(a)")"//   Note that Flot allows additional element data other than just x and y -- we use this in the callback."
    write(iprt,"(a)")"//   We build an array: [x,y,modelnumber,label]."
    write(iprt,"(a)")""
    if (l_pdb) then
      write(iprt,"(a)")"var x = Jmol.evaluateVar(jmolApplet0, ""show('file').lines.find('MODEL')"");"
      write(iprt,"(a)")"var e = []; for (var i = 0;i < x.length; i++){e.push(+x[i].substring(16))};"
    end if
    write(iprt,"(a)")"for (var i = 0; i < modelCount; i++) {"
    write(iprt,"(a)")"  var modelnumber = 0 + Info[i].modelNumber"
    write(iprt,"(a)")"  var name = Info[i].name"
    write(iprt,"(a)")"  var Properties = Info[i].modelProperties"
    if (l_pdb) then
      write(iprt,"(a)")"  var energy =  e[i]"
    else
      write(iprt,"(a)")"  var energy =  parseFloat(name.substring((name.toLowerCase() + "" kc"")."// &
       "split(""kc"")[0].lastIndexOf(""="") + 1)); //parse the name to pull out the energy"
    end if
    if (mode == 2) then
      write(iprt,"(a)")"  var label =  'Model = ' + modelnumber + ', Dipole = ' + roundoff(energy,3) + ' Debye'"
    else
      write(iprt,"(a)")"  var label =  'Model = ' + modelnumber + ', Energy = ' + roundoff(energy,3) + ' Kcal/mol'"
    end if
    write(iprt,"(a)")"  A.push([i+1,0 + energy,modelnumber,label])"
    write(iprt,"(a)")"  }"
    write(iprt,"(a)")""
    write(iprt,"(a)")"// add that data to the array"
    write(iprt,"(a)")""
    write(iprt,"(a)")"data.push(A)"
    write(iprt,"(a)")""
    write(iprt,"(a)")"// select flot options:"
    write(iprt,"(a)")""
    write(iprt,"(a)")"var options = {"
    write(iprt,"(a)")"  lines: { show: true },"
    write(iprt,"(a)")"  points: { show: true, color:""#FF0000"" },"
    write(iprt,"(a)")"  selection: { mode: (nplots == 1 ? ""x"" : ""xy""), hoverMode: (nplots == 1 ? ""x"" : ""xy"") },"
    write(iprt,"(a)")"  grid: { hoverable: true, clickable: true, hoverDelay: 10, hoverDelayDefault: 10 }"
    write(iprt,"(a)")"  }"
    write(iprt,"(a)")""
    write(iprt,"(a)")"// draw the plot"
    write(iprt,"(a)")""
    write(iprt,"(a)")"theplot = $.plot($(""#plotarea""), data, options)"
    write(iprt,"(a)")""
    write(iprt,"(a)")""
    write(iprt,"(a)")""
    write(iprt,"(a)")"// globals for callback"
    write(iprt,"(a)")""
    write(iprt,"(a)")" previousPoint = null"
    write(iprt,"(a)")" item0 = {datapoint:A[0]}"
    write(iprt,"(a)")""
    write(iprt,"(a)")" $(""#plotarea"").unbind(""plothover plotclick"", null)"
    write(iprt,"(a)")" $(""#plotarea"").bind(""plothover"", plotHoverCallback)"
    write(iprt,"(a)")" $(""#plotarea"").bind(""plotclick"", plotClickCallback)"
    write(iprt,"(a)")""
    write(iprt,"(a)")" // execute an initial plotClickCallback:"
    write(iprt,"(a)")""
    write(iprt,"(a)")""
    write(iprt,"(a)")" setTimeout('plotClickCallback(null,null,item0);"// &
      "iamready=true;jmolScript(""animation mode loop 0 0;animation on"")',100)"
    write(iprt,"(a)")""
    write(iprt,"(a)")""
    write(iprt,"(a)")" "
    write(iprt,"(a)")"}"
    write(iprt,"(a)")""
    write(iprt,"(a)")"var iamready = false"
    write(iprt,"(a)")""
    write(iprt,"(a)")"function doHighlight(app, modelIndex) {"
    write(iprt,"(a)")" if (!iamready)return"
    write(iprt,"(a)")" theplot.unhighlight(0,-1)"
    write(iprt,"(a)")""
    write(iprt,"(a)")"  modelIndex = 0 + modelIndex // required for JavaScript"
    write(iprt,"(a)")""
    write(iprt,"(a)")" theplot.highlight(0, modelIndex)"
    write(iprt,"(a)")" var label = data[0][modelIndex][3]"
    write(iprt,"(a)")"        setTimeout('jmolScript(""set echo top left;echo ' + label+'"")',100)"
    write(iprt,"(a)")"}"
    write(iprt,"(a)")""
    write(iprt,"(a)")"var item0"
    write(iprt,"(a)")"var previousPoint = null"
    write(iprt,"(a)")""
    write(iprt,"(a)")"function plotClickCallback(event, pos, item) {"
    write(iprt,"(a)")" // We're getting back  Flot ""event"" that returns the [x,y,modelIndex,label] point clicked."
    write(iprt,"(a)")""
    write(iprt,"(a)")" if (!item)return"
    write(iprt,"(a)")" var model = item.datapoint[2]"
    write(iprt,"(a)")" var label = item.datapoint[3]"
    write(iprt,"(a)")" var script = ' model 1.'+model+';font echo 16; set echo top left;echo ' + label"
    write(iprt,"(a)")""
    write(iprt,"(a)")" // It's important to use jmolScriptWait here, otherwise Jmol "
    write(iprt,"(a)")" // will impose a 100 ms wait of its own in processing the script queue."
    write(iprt,"(a)")""
    write(iprt,"(a)")" jmolScriptWait(script)"
    write(iprt,"(a)")"}"
    write(iprt,"(a)")""
    write(iprt,"(a)")"function plotHoverCallback(event, pos, item) {"
    write(iprt,"(a)")""
    write(iprt,"(a)")" // for hover, we set the ""tooltip"" div -- "
    write(iprt,"(a)")""
    write(iprt,"(a)")" // $(""#xxxx"") is a jQuery wrapper for a div. $(""#xxxx"")[0] is the div itself"
    write(iprt,"(a)")""
    write(iprt,"(a)")" if (item) {"
    write(iprt,"(a)")"  if (previousPoint != item.datapoint) {"
    write(iprt,"(a)")"   $(""#tooltip"").remove()"
    write(iprt,"(a)")"   previousPoint = item.datapoint "
    write(iprt,"(a)")"   var y = item.datapoint[1]"
    write(iprt,"(a)")"   var model = item.datapoint[2]"
    if (mode == 2) then
      write(iprt,"(a)")"   var label = ""&nbsp;&nbsp;Model ""+ model + "", Dipole = "" + y +"" Debye"""
    else
      write(iprt,"(a)")"   var label = ""&nbsp;&nbsp;Model ""+ model + "", Energy = "" + y +"" Kcal/mol)"""
    end if
    write(iprt,"(a)")"   showTooltip(item.pageX, item.pageY + 10, label, pos)"
    write(iprt,"(a)")"  }"
    write(iprt,"(a)")""
    write(iprt,"(a)")"  if (pos.canvasY > 350)plotClickCallback(event, pos, item)"
    write(iprt,"(a)")" } else {"
    write(iprt,"(a)")"  $(""#tooltip"").remove()"
    write(iprt,"(a)")"  previousPoint = null"
    write(iprt,"(a)")" }"
    write(iprt,"(a)")"}"
    write(iprt,"(a)")""
    write(iprt,"(a)")"function showTooltip(x, y, contents, pos) {"
    write(iprt,"(a)")""
    write(iprt,"(a)")" // from the plot hover callback -- create the tooltip div and set its content and style"
    write(iprt,"(a)")""
    write(iprt,"(a)")" if (pos.canvasY > 340) y += (340 - pos.canvasY)"
    write(iprt,"(a)")" $('<div id=""tooltip"">' + contents + '</div>').css( {"
    write(iprt,"(a)")"  position: 'absolute',"
    write(iprt,"(a)")"  display: 'none',"
    write(iprt,"(a)")"  top: y + 5,"
    write(iprt,"(a)")"  left: x + 5,"
    write(iprt,"(a)")"  border: '1px solid #fdd',"
    write(iprt,"(a)")"  padding: '2px',"
    write(iprt,"(a)")"  'background-color': '#fee',"
    write(iprt,"(a)")"  opacity: 0.80"
    write(iprt,"(a)")"  }).appendTo(""body"").fadeIn(200);"
    write(iprt,"(a)")"}"
    write(iprt,"(a)")""
    write(iprt,"(a)")""
    write(iprt,"(a)")"</script>"
    write(iprt,"(a)")"</head>"
    write(iprt,"(a)")"<body>"
    write(iprt,"(a)")""
    write(iprt,"(a)")""
    if (len_trim(koment) == 0 .or. (index(koment(:8), " NULL") /= 0)) then
      write(iprt,"(a)") "<HTML><HEAD><TITLE>"//input_fn(:len_trim(input_fn) - 5)//"</TITLE></HEAD>"
      write(iprt,"(a)") "<h1 align=""center"">"//input_fn(:len_trim(input_fn) - 5)//"</h1>"
    else
      write(iprt,"(a)") "<HTML><HEAD><TITLE>"//trim(koment)//"</TITLE></HEAD>"
      write(iprt,"(a)") "<h1 align=""center"">"//trim(koment)//"</h1>"
    end if
     if (len_trim(title) /= 0 .and. (index(title(:8), " NULL") == 0) ) &
      write(iprt,"(a)") '<h2 align="center">'//trim(title)//'</h2>'
    store_escf = escf
    escf = 0.d0
    call write_data_to_html(iprt)
    escf = store_escf
    write(iprt,"(a)")"<table><tr><td bgcolor=lightblue align=center>"
    write(iprt,"(a)")"<div id=""imagediv"" class=""printonly""></div>"
    write(iprt,"(a)")"<div id=""appletdiv"" style=""width:450;height:470;text-align:center"">"
    write(iprt,"(a)")"<script type=""text/javascript"">"
    write(iprt,"(a)")""
    write(iprt,"(a)")"  Jmol.Info.j2sPath = ""../jsmol/j2s"""
    write(iprt,"(a)")"  "
    write(iprt,"(a)")" jmolInitialize(""java"",""JmolAppletSigned0.jar"")"
    write(iprt,"(a)")" jmolSetAppletColor(""lightblue"");"
    write(iprt,"(a)")" jmolApplet(600, ""set antialiasDisplay;set loadStructCallback 'plotEnergies';"// &
      "set animFrameCallback 'doHighlight'; load "" + modelFile);"
    write(iprt,"(a)")"jmolScript(""set measurementUnits ANGSTROMS; set "// &
      "bondRadiusMilliAngstroms (50); spacefill 15%; "")</script>"
    write(iprt,"(a)")"</script>"
    write(iprt,"(a)")"</div>"
    write(iprt,"(a)")"animation"
    write(iprt,"(a)")" <a href='javascript:jmolScript(""model first"")'>first</a> "
    write(iprt,"(a)")" <a href='javascript:jmolScript(""model last"")'>last</a> "
    write(iprt,"(a)")" <a href='javascript:jmolScript(""model prev"")'>previous</a> "
    write(iprt,"(a)")" <a href='javascript:jmolScript(""model next"")'>next</a> "
    write(iprt,"(a)")" <a href='javascript:jmolScript(""animation mode loop 0 0;animation play"")'>loop</a> "
    write(iprt,"(a)")" <a href='javascript:jmolScript(""animation mode palindrome 0 0;animation play"")'>palindrome</a> "
    write(iprt,"(a)")" <a href='javascript:jmolScript(""animation off"")'>off</a> &nbsp;&nbsp; "
    line = input_fn(:len_trim(input_fn) - 4)//"txt"
    do i = len_trim(line), 1, -1
      if (line(i:i) == "/" .or. line(i:i) == backslash) exit
    end do
    line = line(i+1:)
    write(iprt,"(a)")" <a href='javascript:jmolScript(""script "//backslash//"""common.txt" &
      //backslash//""";"")'>Common</a>&nbsp;&nbsp;"
    write(iprt,"(a)")" <a href='javascript:jmolScript(""script "//backslash//""""//trim(line) &
      //backslash//""";"")'> Script</a>"
    write(iprt,"(a)")"<br>"
    write(iprt,"(a)")"</td>"
    write(iprt,"(a)")"<td bgcolor=lightblue align=center>"
    write(iprt,"(a)")"<table><tr><td width=40></td><td>"
    write(iprt,"(a)")"<div id=""graphdiv"" style=""background-color:lightblue"">"
    write(iprt,"(a)")" <div id=""plottitle""></div>"
    write(iprt,"(a)")" <br>"
    write(iprt,"(a)")" <div id=""plotarea"" style=""width:600px;height:600px;background-color:lightblue""></div>"
    write(iprt,"(a)")"</div>"
    write(iprt,"(a)")"</td></tr></table>"
    write(iprt,"(a)")"</td>"
    write(iprt,"(a)")"</tr>"
    write(iprt,"(a)")"<tr>"
    write(iprt,"(a)")"<td valign=top>"
    write(iprt,"(a)")"<span class=""noprint"">"
    write(iprt,"(a)")"Click on a point, or scan along the baseline to move rapidly among the models."
    write(iprt,"(a)")"<br><br>"
    write(iprt,"(a)")""
    write(iprt,"(a)")"</span>"
    write(iprt,"(a)")"</td>"
    write(iprt,"(a)")"</tr>"
    write(iprt,"(a)")"<tr><td colspan=2>"
    write(iprt,"(a)")""
    write(iprt,"(a)")""
    write(iprt,"(a)")"</td></tr></table>"
    write(iprt,"(a)")""
    write(iprt,"(a)")"</body>"
    write(iprt,"(a)")"</html>"
!
!  Write out a simple script file
!
    close (iprt)
    line = line(:len_trim(line) - 3)
    call add_path(line)
    inquire (file=trim(line)//"txt", exist = exists)
    if (.not. exists) then
      open(unit=iprt, file=trim(line)//"txt")
      i = 0
      do j = 1, len_trim(line)
        if (line(j:j) == "/" .or. line(j:j) == backslash) i = j
      end do
      if (i /= 0) line = line(i + 1:)
      write(iprt,"(a)")"#","# Script for use with the HTML file """//trim(line)//"html""","#"
    end if
    close (iprt)
  end subroutine write_path_html
