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

subroutine rapid3 (xparam, step, pvect, numvar, funct, okf, okc)
    use param_global_C, only : maxpms
    implicit none
    logical, intent (out) :: okc, okf
    integer, intent (inout) :: numvar
    double precision, intent (inout) :: funct, step
    double precision, dimension (numvar), intent (inout) :: pvect, xparam
!---------------------------------------------------------------------------
    logical :: first = .true.
    integer :: center, i, ictr, iquit, left, maxlin, right
    double precision :: aabs, alfs, alpha, beta, drop, energy, eps, estor, &
   & fin, fmax, fmin, funold, gamma, pabs, s, sqstor, ssqlst, stlast, tee, &
   & tiny, xcrit, xmaxm, xminm, xxm, ymaxst
    double precision, dimension (3) :: phi, vt, grad
    double precision, dimension (maxpms) :: xstor
    external exchng
    intrinsic Abs, Max, Min, Sign
    save :: first, drop, xmaxm, i, maxlin, xcrit, ymaxst, eps, tee, energy
!---------------------------------------------------------------------------
    if (first) then
      first = .false.
      drop = 0.002d0
      xmaxm = 0.4d0
      i = 2
      step = 1.d0
      maxlin = 15
      if (numvar == 1) maxlin = 30
      xcrit = 0.0001d0
      ymaxst = 0.4d0
      eps = 10 ** (-i)
      tee = eps
      energy = 0.d0
    end if
    xmaxm = 0.d0
    step = Max (0.d0, Min (100.d0, step))
    do i = 1, numvar
      pabs = Abs (pvect(i))
      xmaxm = Max (xmaxm, pabs)
    end do
    xminm = xmaxm
    xmaxm = ymaxst / xmaxm
    fin = funct
    ssqlst = funct
    iquit = 0
    phi(1) = funct
    vt(1) = 0.0d00
    vt(2) = step / 4.0d00
    if (vt(2) > xmaxm) then
      vt(2) = xmaxm
    end if
    fmax = funct
    fmin = funct
    step = vt(2)
    do i = 1, numvar
      xparam(i) = xparam(i) + step * pvect(i)
    end do
    call rapid2 (xparam, funct, grad, .false.)
    phi(2) = funct
    if (phi(2) > fmax) then
      fmax = phi(2)
    end if
    if (phi(2) < fmin) then
      fmin = phi(2)
    end if
    call exchng (phi(2), sqstor, energy, estor, xparam, xstor, step, alfs, &
   & numvar)
    if (phi(1) <= phi(2)) then
      vt(3) = -vt(2)
      left = 3
      center = 1
      right = 2
    else
      vt(3) = 2.0d00 * vt(2)
      left = 1
      center = 2
      right = 3
    end if
    stlast = vt(3)
    step = stlast - step
    do i = 1, numvar
      xparam(i) = xparam(i) + step * pvect(i)
    end do
    call rapid2 (xparam, funct, grad, .false.)
    if (funct > fmax) then
      fmax = funct
    end if
    if (funct < fmin) then
      fmin = funct
    end if
    if (funct < sqstor) then
      call exchng (funct, sqstor, energy, estor, xparam, xstor, step, alfs, &
     & numvar)
    end if
    if (funct < fin) then
      iquit = 1
    end if
    phi(3) = funct
    okc = .true.
    do ictr = 3, maxlin
      alpha = vt(2) - vt(3)
      beta = vt(3) - vt(1)
      gamma = vt(1) - vt(2)
      if (Abs (alpha*beta*gamma) <= 1.d-20) go to 1000
      alpha = - (phi(1)*alpha+phi(2)*beta+phi(3)*gamma) / (alpha*beta*gamma)
      beta = ((phi(1)-phi(2))/gamma) - alpha * (vt(1)+vt(2))
      if (alpha > 0.0d0) then
        step = -beta / (2.0d00*alpha)
        s = step - stlast
        xxm = 2.0d00 * xmaxm
        if (Abs (s) > xxm) then
          s = Sign (xxm, s) * (1.d0+0.01d0*(xxm/s))
        end if
        step = s + stlast
      else
        if (phi(right) > phi(left)) then
          step = 3.0d00 * vt(left) - 2.0d00 * vt(center)
        else
          step = 3.0d00 * vt(right) - 2.0d00 * vt(center)
        end if
        s = step - stlast
        if (Abs (s) > xmaxm) then
          s = Sign (xmaxm, s) * (1.d0+0.01d0*(xmaxm/s))
        end if
        step = s + stlast
      end if
      if (ictr > 3) then
        aabs = Abs (s*xminm)
        if (aabs < xcrit .and. ictr > 4) go to 1000
      end if
      do i = 1, numvar
        xparam(i) = xparam(i) + s * pvect(i)
      end do
      funold = funct
      call rapid2 (xparam, funct, grad, .false.)
      if (funct > fmax) then
        fmax = funct
      end if
      if (funct < fmin) then
        fmin = funct
      end if
      if (funct < sqstor) then
        call exchng (funct, sqstor, energy, estor, xparam, xstor, step, alfs, &
       & numvar)
      end if
      if (funct < fin) then
        iquit = 1
      end if
    !
    ! TEST TO EXIT FROM RAPID3 IF NOT DROPPING IN VALUE OF FUNCTION FAST.
    !
      tiny = Max ((ssqlst-fmin)*0.2d0, drop)
      tiny = Min (tiny, 0.5d0)
      if (numvar == 1) tiny = 0.01d0*tiny
      if (Abs (funold-funct) < tiny .and. iquit == 1 .and. ictr > 4) go to 1000
      if ((Abs(step-stlast) <= eps*Abs(step+stlast)+tee) .and. (iquit == 1) .and. ictr > 4) &
     & go to 1000
      stlast = step
      if ((step > vt(right)) .or. (step > vt(center) .and. funct < &
     & phi(center)) .or. (step > vt(left) .and. step < vt(center) .and. funct &
     & > phi(center))) then
        vt(left) = step
        phi(left) = funct
      else
        vt(right) = step
        phi(right) = funct
      end if
      if (vt(center) >= vt(right)) then
        i = center
        center = right
        right = i
      end if
      if (vt(left) >= vt(center)) then
        i = left
        left = center
        center = i
      end if
      if (vt(center) >= vt(right)) then
        i = center
        center = right
        right = i
      end if
    end do
    okc = .false.
1000 call exchng (sqstor, funct, estor, energy, xstor, xparam, alfs, step, &
   & numvar)
    okf = (funct < ssqlst)
    if (funct >= ssqlst) return
    if (step < 0.0d0) then
      step = -step
      do i = 1, numvar
        pvect(i) = -pvect(i)
      end do
    end if
end subroutine rapid3
