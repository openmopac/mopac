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

subroutine rapid1 (loop, xparam, numvar_loc, funct1)
!
!
    use param_global_C, only : ifiles_8, parab, contrl, nfns, error, &
    factor, maxpms
    integer, intent (inout) :: numvar_loc
    integer, intent (in) :: loop
    double precision, dimension (numvar_loc), intent (inout) :: xparam
    double precision, intent (out) :: funct1
    double precision, external :: reada
 !-------------------------------------------------------------------------
    logical :: first = .true., fix_parab
    logical :: okc, okf, reset, first_print
    integer :: i, ihdim, ii, ik, iloop, j, jnrst, k, lnstop
    double precision :: alpha, bsmvf, cncadd, Cos, cosine_loc, del, dell, &
   & dott, funct, funct2, ggd, ggggg, gnorm_loc, pmstep, pnlast, pnorm, &
   & renorm, rootv, rst, s, smval, startf, startg, sy, y, yhy
    double precision, dimension (:), allocatable :: hesinv
    double precision, dimension (maxpms) :: gg
    double precision, dimension (maxpms) :: glast
    double precision, dimension (maxpms) :: gd, grad, pvect, xd
    double precision, dimension (maxpms) :: xlast
    double precision, dimension (maxpms) :: xold
    external finish
    double precision, external :: dot
    intrinsic Abs, Max, Min, Sign, Sqrt
    save :: first, fix_parab, del, dell, rootv, rst
    data xlast / maxpms * 0.d0 /
    data glast / maxpms * 0.d0 /
!--------------------------------------------------------------------------
    ggd = 0.d0
    if (first) then
      do i = 1, numvar_loc
        xold(i) = 1.d0
      end do
      first = .false.
      parab = 10000.d0
      i = Index(contrl," PARAB=")
      if(i /= 0)then
      parab = Reada(contrl,i)
      fix_parab = .true.
      else
      fix_parab = .false.
      end if
    !
    !   THE FOLLOWING CONSTANTS SHOULD BE SET BY THE USER.
    !
      rst = 0.05d0
      dell = 0.0000001d0
      del = dell
    !
    !    THESE CONSTANTS SHOULD BE SET BY THE PROGRAM.
    !
      rootv = Sqrt (numvar_loc+1.d-5)
    end if
  !
  !     AND FINALLY, THE FOLLOWING CONSTANTS ARE CALCULATED.
  !
    ihdim = (numvar_loc*(numvar_loc+1)) / 2
    allocate (hesinv(ihdim))
    cncadd = Min (1.0d00/rootv, 0.15d0)
  !
  !     FIRST, WE INITIALISE THE VARIABLES.
  !
    lnstop = 1
    alpha = 0.01d00
    pnorm = .10d00
    jnrst = 0
    Cos = 0.0d00
    iloop = 0
    first_print = .true.
  !
  ! CALCULATE THE VALUE OF THE FUNCTION -> FUNCT1, AND GRADIENTS -> GRAD.
  ! NORMAL SET-UP OF FUNCT1 AND GRAD, DONE ONCE ONLY.
  !
    call rapid2 (xparam, funct1, grad, .true.)
    startg = Sqrt (dot(grad, grad, numvar_loc))
    if (startg >= 1.0d-4) then
      startf = funct1
      do i = 1, numvar_loc
        gd(i) = grad(i)
      end do
    !     *
    !     START OF EACH ITERATION CYCLE ...
    !     *
    !
      reset = .false.
      do
        gnorm_loc = Sqrt (dot(grad, grad, numvar_loc))
        jnrst = jnrst + 1
        if (lnstop /= 1 .and. Cos > rst .and. jnrst < 30) then
        !
        !     *
        !     UPDATE VARIABLE-METRIC MATRIX
        !     *
        !
          sy = 0.0d00
          yhy = 0.0d00
          do i = 1, numvar_loc
            s = 0.0d00
            do k = 1, numvar_loc
              ik = i + numvar_loc * (k-1) - ((k*(k-1))/2)
              if (k > i) then
                ik = k + numvar_loc * (i-1) - ((i*(i-1))/2)
              end if
              s = s + hesinv(ik) * (grad(k)-glast(k))
            end do
            gg(i) = s
            y = grad(i) - glast(i)
            yhy = yhy + gg(i) * y
            sy = sy + (xparam(i)-xlast(i)) * y
          end do
          if (abs(sy) > 1.d-10 .and. abs(yhy) > 1.d-10) then
          do i = 1, numvar_loc
            y = xparam(i) - xlast(i)
            do k = i, numvar_loc
              ik = k + numvar_loc * (i-1) - ((i*(i-1))/2)
              hesinv(ik) = hesinv(ik) + y * (xparam(k)-xlast(k)) / sy - gg &
             & (i) * gg(k) / yhy
            end do
          end do
          end if
        else
          reset = .true.
          do i = 1, numvar_loc
            xd(i) = xparam(i) - Sign (del, grad(i))
          end do
        !
        ! THIS CALL OF COMPFG IS USED TO CALCULATE THE SECOND-ORDER MATRIX IN H
        ! IF THE NEW POINT HAPPENS TO IMPROVE THE RESULT, THEN IT IS KEPT.
        ! OTHERWISE IT IS SCRAPPED, BUT STILL THE SECOND-ORDER MATRIX IS O.K.
        !
          call rapid2 (xd, funct2, gd, .true.)
          do i = 1, ihdim
            hesinv(i) = 0.0d00
          end do
          do i = 1, numvar_loc
            ii = i + numvar_loc * (i-1) - ((i*(i-1))/2)
            ggggg = grad(i) - gd(i)
            if (Abs (ggggg) >= 1.d-12) then
              ggd = Abs (grad(i))
              if (funct2 < funct1) then
                ggd = Abs (gd(i))
              end if
              hesinv(ii) = Sign (del, grad(i)) / ggggg
              if (hesinv(ii) >= 0.0d00 .or. ggd >= 1.d-12) then
                if (hesinv(ii) < 0.0d00) then
                  hesinv(ii) = 6.d0 * del / ggd
                end if
                go to 1000
              end if
            end if
            hesinv(ii) = 0.01d00
1000        if (ggd < 1.d-12) then
              ggd = 1.d-12
            end if
            pmstep = Abs (0.1d0/ggd)
            if (hesinv(ii) > pmstep) then
              hesinv(ii) = pmstep
            end if
          end do
          jnrst = 0
          if (funct2 < funct1) then
            funct1 = funct2
            gnorm_loc = 0.0d00
            do i = 1, numvar_loc
              xparam(i) = xd(i)
              grad(i) = gd(i)
              gnorm_loc = gnorm_loc + grad(i) ** 2
            end do
            gnorm_loc = Sqrt (gnorm_loc)
          end if
        end if
        do
        !
        !     *
        !     ESTABLISH NEW SEARCH DIRECTION
        !     *
          pnlast = pnorm
          pnorm = 0.0d00
          dott = 0.0d00
          do k = 1, numvar_loc
            s = 0.0d00
            do i = 1, numvar_loc
              ik = i + numvar_loc * (k-1) - ((k*(k-1))/2)
              if (k > i) then
                ik = k + numvar_loc * (i-1) - ((i*(i-1))/2)
              end if
              s = s - hesinv(ik) * grad(i)
            end do
            pvect(k) = s
            pnorm = pnorm + pvect(k) ** 2
            dott = dott + pvect(k) * grad(k)
          end do
          pnorm = Sqrt (pnorm)
          Cos = -dott / (pnorm*gnorm_loc)
          if (jnrst == 0) exit
          if (Cos > cncadd) then
            if (Cos > rst) then
              i = 0
              exit
            end if
          end if
          pnorm = pnlast
          reset = .true.
          do i = 1, numvar_loc
            xd(i) = xparam(i) - Sign (del, grad(i))
          end do
        !
        ! THIS CALL OF COMPFG IS USED TO CALCULATE THE SECOND-ORDER MATRIX IN H
        ! IF THE NEW POINT HAPPENS TO IMPROVE THE RESULT, THEN IT IS KEPT.
        ! OTHERWISE IT IS SCRAPPED, BUT STILL THE SECOND-ORDER MATRIX IS O.K.
        !
          call rapid2 (xd, funct2, gd, .true.)
          do i = 1, ihdim
            hesinv(i) = 0.0d00
          end do
          do i = 1, numvar_loc
            ii = i + numvar_loc * (i-1) - ((i*(i-1))/2)
            ggggg = grad(i) - gd(i)
            if (Abs (ggggg) >= 1.d-12) then
              ggd = Abs (grad(i))
              if (funct2 < funct1) then
                ggd = Abs (gd(i))
              end if
              hesinv(ii) = Sign (del, grad(i)) / ggggg
              if (hesinv(ii) >= 0.0d00 .or. ggd >= 1.d-12) then
                if (hesinv(ii) < 0.0d00) then
                  hesinv(ii) = 6.d0 * del / ggd
                end if
                go to 1100
              end if
            end if
            hesinv(ii) = 0.01d00
1100        if (ggd < 1.d-12) then
              ggd = 1.d-12
            end if
            pmstep = Abs (0.1d0/ggd)
            if (hesinv(ii) > pmstep) then
              hesinv(ii) = pmstep
            end if
          end do
          jnrst = 0
          if (funct2 < funct1) then
            funct1 = funct2
            gnorm_loc = 0.0d00
            do i = 1, numvar_loc
              xparam(i) = xd(i)
              grad(i) = gd(i)
              gnorm_loc = gnorm_loc + grad(i) ** 2
            end do
            gnorm_loc = Sqrt (gnorm_loc)
          end if
        end do
        lnstop = 0
        alpha = alpha * pnlast / pnorm
        do i = 1, numvar_loc
          glast(i) = grad(i)
          xlast(i) = xparam(i)
        end do
        if (jnrst == 0) then
          alpha = 0.1d00
        end if
        iloop = iloop + 1
     !   if (Abs (smval-funct1) < 1.d-5) exit
        if (Mod(iloop,200) == 0) then
          if (first_print) then
            write(ifiles_8,*)" Mini-cycle  Error Function  Gradient"
            write(ifiles_8,'(i7,f16.2, f14.2)')0,startf,startg
             first_print = .false.
          end if
          write(ifiles_8,'(i7,f16.2, f14.2)')iloop,funct1,gnorm_loc
          endfile (ifiles_8)
          backspace (ifiles_8)
        end if
        if (jnrst /= 0 .and. gnorm_loc < 0.01d0*startg .or. iloop > 2000) then
          i = 0
          exit
        end if
        smval = funct1
        call rapid3 (xparam, alpha, pvect, numvar_loc, funct1, okf, okc)
      !   WE WANT ACCURATE DERIVATIVES AT THIS POINT
      !
      !   RAPID3 DOES NOT GENERATE ANY DERIVATIVES, THEREFORE RAPID2 MUST BE
      !   CALLED TO END THE SEARCH
      !
        if (reset) then
          do j = 1, numvar_loc
            grad(j) = 0.d0
          end do
          reset = .false.
        end if
        call rapid2 (xparam, funct1, grad, .true.)
        if (numvar_loc == 1) then
          gnorm_loc = Abs(grad(1))
          i = 0
          exit
        end if
        bsmvf = Abs (smval-funct1)
        if (bsmvf > 10.d00) then
          Cos = 0.0d00
        end if
        del = 0.0000002d00
        if (bsmvf > 1.0d00) then
          del = dell / 2.0d00
        end if
      !
      ! END OF ITERATION LOOP, EVERYTHING IS STILL O.K. SO GO TO
      ! NEXT CYCLE OF OPTIMIZATION.
      !
        if (bsmvf > 5.0d00) then
          del = dell
        end if
        if (Cos < rst) then
          do i = 1, numvar_loc
            gd(i) = 0.5d0
          end do
        end if
      end do
    !
    !   RAPID1 IS ENDING PROPERLY. THIS IS IMMEDIATELY BEFORE THE RETURN.
    !
      call rapid2 (xparam, funct, grad, .false.)
      do i = 1, numvar_loc
        funct = funct - xparam(i) ** 2 * parab
      end do
    !
    !  Alter damping factor using the cosine of the step with the old step
    !
      renorm = Sqrt (dot(xparam, xparam, numvar_loc)*dot(xold, xold, numvar_loc))
      cosine_loc = dot (xparam, xold, numvar_loc) / renorm
      do i = 1, numvar_loc
        xold(i) = xparam(i)
      end do
!
!  Set factor to 1.d0 in order to calculate "raw" ssq and grad.
!
      error(1:nfns) = factor(1:nfns)
      factor(1:nfns) = 1.d0
      pvect = 0.d0
      call rapid2(pvect, startf, grad,.true.)
  !    startg = Sqrt(dot(grad,grad,numvar_loc))
      call rapid2(xparam, funct, grad, .false.)
  !    write(ifiles_8,'(5f12.6)')xparam(1:5),grad(1:5)
      factor(1:nfns) = error(1:nfns)
   !   gnorm_loc = Sqrt(dot(grad,grad,numvar_loc))
      if(.not. fix_parab) then
      if (cosine_loc > 0.6d0) then
        parab = Max (parab/1.3d0, -10.d0)
      else
        parab = Min (Max(10.d0,parab*1.3d0), 5000000.d0)
      end if
      end if
      write (ifiles_8, "(i4,A,F12.5,A,F17.5,a,i4)") loop, " COSINE ", cosine_loc, " PARAB", &
     & parab, " CYCLES:", iloop
      write (ifiles_8, "(I4,A,2F14.2,A,2F15.1)") loop, " START SSQ:", &
     & startf, funct, "   GRAD:", startg, gnorm_loc
     if(startg < 5.d-1) then
       write(ifiles_8,*) "Optimization Successful"
       call finish
     end if
      if (startf-funct >=1.d-1) return
    else
    write(ifiles_8,*) "Optimization Successful"
    call finish
    end if
end subroutine rapid1
