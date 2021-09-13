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

subroutine partab
!
    use param_global_C, only : power, ifiles_8, refers, names, is_a_ref, &
    contrl, refgeo, refher, refder, refcer, reftot, refger, maxmol, &
    maxrfs, nmols, refhof, refips, molnum, wtgeo, wthof, refdip, &
    wtdip, wtips
!
    use molkst_C, only :  natoms, nvar, tleft, jobnam, escf, &
    & gnorm, last, koment, moperr, norbs, cp298 => temp_1, s298 => temp_2
!
    use chanel_C, only : ir, iw, ilog
    use common_arrays_C, only : labels, loc, xparam, nat, &
    na, nb, nc, c, cb, geometry => geo
!
!------------------------------------------------------------------------
!
    implicit none
    character (len=5) :: wrk2
    character (len=8) :: els
    character (len=40) :: formla
    character (len=96) :: wrk
    logical :: ldips = .false., lgeos = .false., lheat = .false., lips = &
   & .false., large, opend, cp, one_scf
    integer :: i, itime, itxt, j, k, loop,  nerr, ndmols, nimols, nbmols, &
    & nhmols, namols
    double precision :: sum, aveher, aveder, aveier, aveber, aveaer, &
    & rmsher, rmsder, rmsier, rmsger, sumher, sumder, sumier, sumger, &
    & sumerr
    character (len=20), dimension (300) :: ttype
    character (len=5), dimension (maxmol, 4) :: reftxt
    character (len=13), dimension (maxrfs, 4) :: allref
    integer, dimension (4) :: ipfile
    integer, dimension (4) :: nzs
    integer, dimension (4) :: sortl
    double precision :: errors(300), yparam(300), dummy(2)
    double precision, external :: seconds
    intrinsic Abs, Index, Int, Min
    data ipfile / 16, 17, 18, 19 /
!------------------------------------------------------------------------
   large = (Index (contrl, "LARGE") /= 0)
   cp = (Index (contrl, " CP ") /= 0)
   i = len_trim(jobnam)
    open (unit=ipfile(1), form="FORMATTED", status="UNKNOWN", &
   & file=jobnam(:i)//".heats")
    open (unit=ipfile(2), form="FORMATTED", status="UNKNOWN", &
   & file=jobnam(:i)//".ips")
    open (unit=ipfile(3), form="FORMATTED", status="UNKNOWN", &
   & file=jobnam(:i)//".dips")
    open (unit=ipfile(4), form="FORMATTED", status="UNKNOWN", &
   & file=jobnam(:i)//".geos")
   allref(1,1) = " "
   reftxt(1,1) = " "
    do i = 1, 4
      call psort (refers(1, i), 0, nmols, reftxt(1, i), allref(1, i), &
     & ipfile(i), sortl(i), .false. )
    end do
    write (ipfile(1), "(' Formula        Molecule              ',26x,'  &
   &Heat of Formation   Diff. Ref.')")
    write (ipfile(1), "(67X,'Exp.     Calc.')")
    write (ipfile(2), "(' Formula        Molecule              ',20x,'Io&
   &nization Potential  Diff. Ref.')")
    write (ipfile(2), "(67X,'Exp.     Calc.')")
    write (ipfile(4), "(' Formula        Molecule              ',6x,' Geome&
   &tric        Exp.    Calc.   Ref.')")
    write (ipfile(4), "(45X,'Parameter')")
    write (ipfile(3), "(' Formula        Molecule              ','            &
   &     Dipole         Diff. Ref.')")
    write (ipfile(3), "(67X,'Exp.     Calc.')")
    write (ifiles_8, "(a,26x,2A)") "   Molecule     ","         Time  Heat      ", "Dipole   I.P.     Geometry     Total"
    itime = Int (seconds (1))
     if (large) then
        inquire (unit=iw, opened=opend)
        if (.not. opend) then
        i = Index (jobnam, " ") - 1
        open (unit=iw, file=jobnam(:i)//".arc", status="UNKNOWN")
        rewind (iw)
        end if
        else
        open (unit=iw, status="SCRATCH", form="FORMATTED")
      !
      !  Standard MOPAC output is NOT wanted, therefore send it to scratch.
      !
        inquire (unit=iw, opened=opend)
        if (opend) then
          close (unit=iw, status="DELETE")
        end if
        open (unit=iw, status="SCRATCH", form="FORMATTED")
      !
      !  Standard MOPAC LOG files are NOT wanted, therefore send
      !  them to scratch.
      !
        inquire (unit=ilog, opened=opend)
        if (opend) then
          close (unit=ilog, status="DELETE")
        end if
        open (unit=ilog, status="SCRATCH", form="FORMATTED")
      end if
    nhmols = 0
    ndmols = 0
    nimols = 0
    nbmols = 0
    namols = 0
    aveher = 0.d0
    aveder = 0.d0
    aveier = 0.d0
    aveaer = 0.d0
    aveber = 0.d0
    rmsher = 0.d0
    rmsder = 0.d0
    rmsier = 0.d0
    rmsger = 0.d0
    sumher = 0.d0
    sumder = 0.d0
    sumier = 0.d0
    sumger = 0.d0
    one_scf = (index(contrl, "PKA") + index(contrl, "1SCF")/= 0)
    do loop = 1, nmols  
      molnum = molnum + 1
    !
    !  Restore all information for molecule number LOOP
    !
      call calpar
      call getmol (loop)
      deallocate(c, cb)
      allocate(c(norbs, norbs), cb(norbs, norbs))
      write(iw,*)trim(koment)
      gnorm = 0
      last = 0
      call empire (formla, nat, els)
      moperr = .false.
      tleft = 1.d7
    !
    !  Compute the value of all errors.  Gradient limit is set
    !  to -0.10 to ensure that geometries are OK.
    !

      if (one_scf) then
        call compfg (xparam, .true., escf, .true., dummy, .false.)      
      else
        call optgeo (xparam, yparam, nvar, refgeo, -0.1d0)
      end if
      call savgeo (loop, geometry, na, nb, nc, xparam, loc)
      call parfg (errors, ttype, nerr, loop, .false.)
      if (cp) then
        cp298 = 0.d0
        s298 = 0.d0                                                                                                                                                                                                                                                                                                         
        call force()
        refher = cp298 - refhof
        refder = s298 - refdip
      end if
      endfile (ifiles_8)
      backspace (ifiles_8)
    !
    !  Write out information on this system
    !
      i = Int (seconds (1))
      if (is_a_ref(loop)) then
      write (ifiles_8, "(A48,I5)") koment, i - itime
        cycle
      end if
      write (ifiles_8, "(A48,I5,F9.3,F10.3,F7.2,F12.3,F12.3)") koment, i - itime, refher, &
        refder, refcer, refger, refher + refder + refcer + refger
      itime = i
    !
    !  Write out the results in a suitable format.
    !
    !                       Heat of Formation
    !
      if ((Abs (refher) > 1.d-20 .or. natoms == 1 .and. refhof > 0.d0) .and. (index(koment, " TS") /= 0 .or. .not. one_SCF)) then
        write (ipfile(1), "(A15,A47,F9.2,F10.2,F9.2,1X,A5,3X,A,5X,A)") formla, &
       & koment, refhof, refher + refhof, refher, reftxt(loop, 1), els, trim(names(loop))
   !     if (method_pm6) write(51,"(a)")koment(:len_trim(koment))
        lheat = .true.
        nhmols = nhmols + 1
        aveher = aveher + Abs(refher)
          rmsher = rmsher + Abs(refher) ** power
      end if
    !
    !                       Ionization Potential
    !
      if (Abs (refcer) > 1.d-20) then
        write (ipfile(2), "(A15,A47,F9.2,F10.2,F9.2,1X,A5,3X,A,5X,A)") formla, &
       & koment, refips, refcer + refips, refcer, reftxt(loop, 2), els, trim(names(loop))
  !      if (method_pm6) write(52,"(a)")koment(:len_trim(koment))
        lips = .true.
        nimols = nimols + 1
        aveier = aveier + Abs(refcer)
        rmsier = rmsier + Abs(refcer) ** power
      end if
    !
    !                       Dipole Moment
    !
      if (Abs (refder) > 1.d-20) then
        write (ipfile(3), "(A15,A47,F9.2,F10.2,F9.2,1X,A5,3X,A,5X,A)") formla, &
       & koment, refdip, refdip + refder, refder, reftxt(loop, 3), els, trim(names(loop))
   !     if (method_pm6) write(53,"(a)")koment(:len_trim(koment))
        ldips = .true.
        ndmols = ndmols + 1
        aveder = aveder + Abs(refder)
        rmsder = rmsder + Abs(refder) ** power
      end if
    !
    !                       Geometric Variable
    !
      if (refger > 1.d-20) then
        j = Min (300, nvar)
        itxt = 0
        do i = 1, j
         if (refgeo(i) /= " ") then
!
!  Put into nzs the atomic numbers of the atoms used in the geometric variable
!
!   k = atom number of atom of geometric variable
!   j = 1 - bond length, 2 - angle, 3 - dihedral
!
          k = loc(1, i)
          j = loc(2, i) + 1
          select case (j)
          case (2)  !  Bond length
            if (labels(na(k)) > labels(k)) then
              nzs(2) = labels(na(k))
              nzs(1) = labels(k)
            else
              nzs(2) = labels(k)
              nzs(1) = labels(na(k))
            end if
          case (3)
            if (labels(nb(k)) > labels(k)) then
              nzs(3) = labels(nb(k))
              nzs(2) = labels(na(k))
              nzs(1) = labels(k)
            else
              nzs(3) = labels(k)
              nzs(2) = labels(na(k))
              nzs(1) = labels(nb(k))
            end if
          case (4)
          if (labels(nc(k)) > labels(k)) then
              nzs(4) = labels(nc(k))
              nzs(3) = labels(nb(k))
              nzs(2) = labels(na(k))
              nzs(1) = labels(k)
            else              
              nzs(4) = labels(k)
              nzs(3) = labels(na(k))
              nzs(2) = labels(nb(k))
              nzs(1) = labels(nc(k))
            end if
          end select
          sum = 1.d0
          if (loc(2, i) > 1) then
            sum = 180 / 3.141592653598d0
          end if
          if (itxt == 0 .and. refgeo(i) /= " ") then
            wrk = formla(:15) // koment(1:81)
            wrk2 = reftxt(loop, 4)
            itxt = 1
          else
            wrk = " "
            wrk2 = " "
          end if
!
!  xparam:  experimental reference
!  yparam:  calculated geometric parameter
!
              if (sum > 5.d0) then
              write (ipfile(4), "(A47,A12,F8.3,F8.3,1X,A5,7X,F8.3,4I3)") wrk, &
             & refgeo(i), xparam(i) * sum, yparam(i) * sum, wrk2, &
             & (yparam(i)-xparam(i)) * sum, (nzs(k), k=1, j)
  !           if (method_pm6)    write(50,"(a)")wrk(16:len_trim(wrk))
             namols = namols + 1
             aveaer = aveaer + Abs((yparam(i)-xparam(i)) * sum)
            else
              write (ipfile(4), "(A47,A12,2F8.3,1X,A5,F7.3,8X,4I3)") wrk, &
             & refgeo(i), xparam(i) * sum, yparam(i) * sum, wrk2, &
             & (yparam(i)-xparam(i)) * sum, (nzs(k), k=1, j)
  !           if (method_pm6) write(50,"(a)")wrk(16:len_trim(wrk))
             nbmols = nbmols + 1
             aveber = aveber + Abs((yparam(i)-xparam(i)) * sum)
            end if
            rmsger = rmsger + Abs(errors(i)/wtgeo) ** power
            lgeos = .true.
          end if
        end do
      end if
        sumher = sumher + Abs(refher*wthof) ** power
        sumder = sumder + Abs(refder*wtdip) ** power
        sumier = sumier + Abs(refcer*wtips) ** power
        sumger = sumger + refger
        reftot = Abs(refher*wthof) ** power + Abs(refder*wtdip) ** power + &
        & Abs(refcer*wtips) ** power + refger
        write(iw,*)trim(koment)//" end"
    end do
!
!  Use iostat, but ignore the result - causes "graceful" exit
!
    if ( .not. lheat) then
      close (unit=ipfile(1), status="DELETE", iostat=i)
    end if
    if ( .not. lips) then
      close (unit=ipfile(2), status="DELETE", iostat=i)
    end if
    if ( .not. ldips) then
      close (unit=ipfile(3), status="DELETE", iostat=i)
    end if
    if ( .not. lgeos) then
      close (unit=ipfile(4), status="DELETE", iostat=i)
    end if
    close (unit=ir, status="DELETE", iostat=i)
    if( .not. large) close (unit=iw, status="DELETE") 
  !
  !  Summarize the errors
  !
    sumerr = sumher + sumder + sumier + sumger
    write (ifiles_8, "(25x,A23,F14.3,F9.2,F8.2,F12.3,F12.3)") " TOTALS", sumher, &
   & sumder, sumier, sumger, sumerr
    aveher = aveher / (nhmols+0.000001d0)
    aveder = aveder / (ndmols+0.000001d0)
    aveier = aveier / (nimols+1.d-4)
    
    write (ifiles_8, "(25x,'  UNSIGNED AVE. ERROR   ',i4,f8.2,i4,f6.2,i3,f5.2,i6,f6.3)")&
    & nhmols,  aveher, ndmols, aveder, nimols, aveier, namols + nbmols, &
    & (aveaer*3.141592653598d0/180 + aveber)/(namols + nbmols + 1.d-4)
    aveber = aveber / (nbmols + 1.d-4)
    aveaer = aveaer / (namols + 1.d-4)
    write(ifiles_8,'(/,a,/30x,a)')"     Average Errors for Various Quantities", &
    &"No. in set"
    if (nhmols > 0) write(ifiles_8,"(a,f7.2,a,i6)")' Delta Hf:   ',&
    & aveher,"  kcal/mol",nhmols
    if (ndmols > 0) write(ifiles_8,"(a,f7.2,a,i6)")' Dipole:     ',&
    & aveder,"  Debye   ",ndmols
    if (nimols > 0) write(ifiles_8,"(a,f7.2,a,i6)")' I.P.:       ',&
    & aveier,"  eV      ",nimols
    if (nbmols > 0) write(ifiles_8,"(a,f8.3,a,i5)")' Bond length:',&
    & aveber," Angstroms",nbmols
    if (namols > 0) write(ifiles_8,"(a,f7.2,a,i6)")' Angle:      ',&
    & aveaer,"  Degrees ",namols
    write (ifiles_8, "(/,2A)") " All systems calculated; " // "now generating re&
   &ferences"
    do i = 1, 4
      inquire(unit=ipfile(i), opened=opend) 
      if ( .not. opend) cycle
      write(ipfile(i),"(20a)")("****",j=1,20)
      call psort (refers(1, i), 1, nmols, reftxt(1, i), allref(1, i), &
     & ipfile(i), sortl(i), .false.)
    end do
    write (ifiles_8, "(2A)") " References generated"
end subroutine partab
