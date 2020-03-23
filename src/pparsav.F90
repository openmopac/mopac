  subroutine pparsav(save_parameters)
    use param_global_C, only : valvar, valold, toplim, &
    botlim, numvar, ifiles_8, locvar, contrl, fnsnew, penalty
!
    use chanel_C, only : iext
    use molkst_C, only : jobnam, line
    use parameters_C, only : partyp
    use elemts_C, only : elemnt
    implicit none
    logical :: save_parameters
!-----------------------------------------------------------------------
    integer :: i, j, k
    double precision :: penalty_fn, sum
    character :: elemnt2*2, num*1
    intrinsic Index
!-----------------------------------------------------------------------
    if(save_parameters) then
      if (numvar == 0) return
      k = len_trim(jobnam)
      i = Index(contrl, "NEW_RP=")
      if (i /= 0) then
      j = Index(contrl(i:)," ") + i - 2
      if(contrl(j:j) /= "/") then
        j = j + 1
        contrl(j:j) = "/"
      end if
      line = contrl(i+13:j)//jobnam (:k) // "." // "rp"
      else
      line = jobnam (:k) // "." // "rp"
      end if
      open (unit=iext, form="FORMATTED", status="UNKNOWN", file=line)
      rewind (iext)
      write(iext,"('*',/,2a,/,'*')")'* Parameter    Element     New value   ',&
      & '    Change   Limits: Low    High      Penalty     Gradient'
      write(ifiles_8,"(2a)")'  Parameter    Element     New value   ',&
      & '    Change   Limits: Low    High      Penalty     Gradient'
      do i = 1, numvar
        j = locvar(2, i)/200
        if(j /= 0) then
          elemnt2 = elemnt(j)
           if(elemnt2(1:1) == " ")elemnt2 = elemnt2(2:2)
          j = locvar(2,i) - j*200
        else
          elemnt2 = " "
          j = locvar(2,i)
        endif
        k = i
        call lockit(valvar(i), k)
        penalty_fn = penalty*((Max(0.d0, valvar(i)-toplim(i))+Min(0.d0, valvar(i)-botlim(i))))**2
        if (locvar(1, i) == 41) then
          num = "1"
          if (j > 9) num = "2"
          write(elemnt2,'(i'//num//')')j
          write (line, "(4X,'PAR',a2,11X,2F15.8,F12.2,F8.2,f12.2,f14.2)") elemnt2, &
          &             valvar(i), valvar(i) - valold(i), botlim(i), toplim(i), &
          & penalty_fn,fnsnew(i)
        else
          write (line, "(4X,A7,5X,A2,2X,2F15.8,F12.2,F8.2,f12.2,f14.2)") partyp(locvar(1, i))//elemnt2, &
          & elemnt(j), valvar(i), valvar(i) - valold(i), botlim(i), toplim(i), &
          & penalty_fn,fnsnew(i)
        endif
        if (valvar(i)-toplim(i) < 0.d0 .and. valvar(i)-botlim(i) > 0.d0) line(52:82) = " "
        write (iext, "(a)") " "//trim(line)
      end do
    end if
!
!  For parameters that should always be a maximum or a minimum, modify
!  the limits to impose a mild force on the parameter
!
    do i = 1, numvar
      if (botlim(i) < -1.d3) then
!
! This parameter should be as negative as possible.  Set the upper bound
! slightly lower than the actual parameter. (difference is such that the force = botlim*0.0001)
!
        toplim(i) = valvar(i) - 0.000158114*sqrt(-botlim(i))
      end if
      if (toplim(i) > 1.d3) then
!
! This parameter should be as positive as possible.  Set the lower bound
! slightly higher than the actual parameter.
!
        botlim(i) = valvar(i) + 0.000158114*sqrt(toplim(i))
      end if
    end do
    sum = 0.d0
    do i = 1, numvar
      j = locvar(2, i)/200
      if (j /= 0) then
        elemnt2 = elemnt(j)
        if(elemnt2(1:1) == " ")elemnt2 = elemnt2(2:2)
        j = locvar(2,i) - j*200
      else
         elemnt2 = " "
         j = locvar(2,i)
      endif
      k = i
      call lockit(valvar(i), k)
      penalty_fn = penalty*((Max(0.d0, valvar(i)-toplim(i))+Min(0.d0, valvar(i)-botlim(i))))**2
      if (locvar(1, i) == 41) then
        num = "1"
        if (j > 9) num = "2"
        write(elemnt2,'(i'//num//')')j
        write (line, "(4X,'PAR',a2,11X,2F15.8,F12.2,F8.2,f12.2,f14.2)") elemnt2, &
        &             valvar(i), valvar(i) - valold(i), botlim(i), toplim(i), &
        & penalty_fn,fnsnew(i)
      else
        write (line, "(4X,A7,5X,A2,2X,2F15.8,F12.2,F8.2,f12.2,f14.2)") partyp(locvar(1, i))//elemnt2, &
        & elemnt(j), valvar(i), valvar(i) - valold(i), botlim(i), toplim(i), &
        & penalty_fn,fnsnew(i)
      endif
      if (valvar(i)-toplim(i) < 0.d0 .and. valvar(i)-botlim(i) > 0.d0) line(52:82) = " "
      write (ifiles_8, "(a)") " "//trim(line)
    end do
    do i = 1, numvar
      valold(i) = valvar(i)
    end do
    close (iext)
    if (sum > 1.d-6) write(ifiles_8,"(a,f12.2)")' Total penalty function:', sum
    if(save_parameters) write (ifiles_8,*) " PARAMETERS DUMPED O.K."
  end subroutine pparsav
