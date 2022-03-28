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

      subroutine linmin(xparam, alpha, pvect, nvar, funct, okf, ic, dott)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use molkst_C, only : moperr, numcal, keywrd, cosine, line
      use common_arrays_C, only : loc
      use chanel_C, only : iw
      use funcon_C, only : pi
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: nvar
      integer , intent(inout) :: ic
      double precision  :: alpha
      double precision, intent (out)  :: funct
      double precision , intent(in) :: dott
      logical , intent(inout) :: okf
      double precision  :: xparam(nvar)
      double precision , intent(inout) :: pvect(nvar)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: left, right, center, icalcn, i, maxlin, nsame, ictr
      double precision, dimension(3) :: phi
      double precision, dimension(4) :: vt, grad ! grad is not used in linmin
      double precision, dimension(nvar) :: xstor, xparef
      double precision :: alpold, xmaxm, delta1, delta2, xnear, ymaxst, pabs, &
        fin, ssqlst, fmax, fmin, sqstor, energy = 0.d0, estor, alfs, beta, gamma, s&
        , xxm, funold, hlast, drop, sum
      logical :: print, diis

      save print, alpold, icalcn, xmaxm, i, delta1, delta2, xnear, maxlin, &
        ymaxst
!-----------------------------------------------
!*********************************************************************
!
!  LINMIN DOES A LINE MINIMISATION.
!
!  ON INPUT:  XPARAM = STARTING COORDINATE OF SEARCH.
!             ALPHA  = STEP SIZE FOR INITIATING SEARCH.
!             PVECT  = DIRECTION OF SEARCH.
!             NVAR   = NUMBER OF VARIABLES IN XPARAM.
!             FUNCT  = INITIAL VALUE OF THE FUNCTION TO BE MINIMIZED.
!             ISOK   = NOT IMPORTANT.
!             COSINE = COSINE OF ANGLE OF CURRENT AND PREVIOUS GRADIENT.
!
!  ON OUTPUT: XPARAM = COORDINATE OF MINIMUM OF FUNCTI0N.
!             ALPHA  = NEW STEP SIZE, USED IN NEXT CALL OF LINMIN.
!             FUNCT  = FINAL, MINIMUM VALUE OF THE FUNCTION.
!             OKF    = TRUE IF LINMIN IMPROVED FUNCT, FALSE OTHERWISE.
!
!**********************************************************************
!
!  THE FOLLOWING COMMON IS USED TO FIND OUT IF A NON-VARIATIONALLY
!  OPTIMIZED WAVE-FUNCTION IS BEING USED.
      data icalcn/ 0/
      data alpold/ 0.D0/
      if (icalcn /= numcal) then
        xmaxm = 0.4D0
        delta2 = 0.001D0
        if (index(keywrd,'NOTH') == 0) then
          delta1 = 0.5D0
        else
          delta1 = 0.1D0
        end if
        alpha = 1.D0
        maxlin = 15
        xnear = 1.D-4
        if (nvar == 1) then
          pvect(1) = 0.01D0
          alpha = 1.D0
          xnear = 1.D-5
          delta1 = 0.0005D0
          delta2 = 0.0001D0
          if (index(keywrd,'PREC') /= 0) then
            delta1 = 0.D0
            delta2 = delta2*0.01D0
          end if
          maxlin = 50
        end if
        cosine = 99.99D0
!
        ymaxst = 0.4D0
        print = index(keywrd,'LINMIN') /= 0
        icalcn = numcal
      end if
      nsame = 0
      xparef(:nvar) = xparam(:nvar)
      xmaxm = 0.D0
      do i = 1, nvar
        pabs = abs(pvect(i))
        xmaxm = max(xmaxm,pabs)
      end do
      xmaxm = ymaxst/xmaxm
      if (nvar == 1) call compfg (xparam, .TRUE., funct, .TRUE., grad, .FALSE.)
      if (moperr) return
      fin = funct
      ssqlst = funct
      diis = ic==1 .and. nvar>1
      phi(1) = funct
      alpha = 1.D0
      vt(1) = 0.0D00
      vt(2) = alpha
      vt(2) = min(xmaxm,vt(2))
      fmax = funct
      fmin = funct
      alpha = vt(2)
      i = 1
      if (nvar > 0) then
        xparam(:nvar) = xparef(:nvar) + alpha*pvect(:nvar)
        i = nvar + 1
      end if
      call compfg (xparam, .TRUE., funct, .TRUE., grad, .FALSE.)
      phi(2) = funct
      if (moperr) return
      fmax = dmax1(phi(2),fmax)
      fmin = min(phi(2),fmin)
      call exchng (phi(2), sqstor, energy, estor, xparam, xstor, alpha, alfs, &
        nvar)
      if (.not.diis) then
        if (nvar > 1) then
!
!   CALCULATE A NEW ALPHA BASED ON THIEL'S FORMULA
!
          alpha = -alpha**2*dott/(2.D0*(phi(2)-ssqlst-alpha*dott))
          alpha = min(2.D0,alpha)
        else
          if (phi(2) < phi(1)) then
            alpha = 2*alpha
          else
            alpha = -alpha
          end if
        end if
!#      IF(PRINT)WRITE(IW,'(3(A,F12.6))')' ESTIMATED DROP:',DOTT*0.5D0,
!#     1'  ACTUAL: ',PHI(2)-SSQLST, '  PREDICTED ALPHA',ALPHA
        okf = okf .or. phi(2)<ssqlst
        if (delta1 > 0.3D0) then
!
!  THIEL'S TESTS # 18 AND 19
!
          if (okf .and. alpha<2.D0) go to 190
        end if
        vt(3) = alpha
        if (vt(3) <= 1.D0) then
          left = 3
          center = 1
          right = 2
        else
          left = 1
          center = 2
          right = 3
        end if
        i = 1
        if (nvar > 0) then
          xparam(:nvar) = xparef(:nvar) + alpha*pvect(:nvar)
          i = nvar + 1
        end if
        call compfg (xparam, .TRUE., funct, .TRUE., grad, .FALSE.)
        if (moperr) return
        fmax = dmax1(funct,fmax)
        fmin = min(funct,fmin)
        if (funct < sqstor) call exchng (funct, sqstor, energy, estor, xparam, &
          xstor, alpha, alfs, nvar)
        okf = okf .or. funct<fin
        phi(3) = funct
        if (print) write (iw, 50) vt(1), phi(1), phi(1) - fin, vt(2), phi(2), &
          phi(2) - fin, vt(3), phi(3), phi(3) - fin
          if (print) then
          call to_screen(' ---QLINMN ')
          write(line,"(5x,'LEFT    ...',f17.8,2f17.7)")vt(1), phi(1), phi(1) - fin
          call to_screen(line)
          write(line,"(5x,'CENTER  ...',f17.8,2f17.7)") phi(2), phi(2) - fin, vt(3)
          call to_screen(line)
          write(line,"(5x,'RIGHT   ...',f17.8,2f17.7)")vt(3), phi(3), phi(3) - fin
          call to_screen(line)
          end if
   50   format(' ---QLINMN ',/,5x,'LEFT    ...',f17.8,2f17.7,/,5x,'CENTER  ...',&
          f17.8,2f17.7,/,5x,'RIGHT   ...',f17.8,2f17.7,/)
        alpold = 0.D0
        l180: do ictr = 3, maxlin
          alpha = vt(2) - vt(3)
          beta = vt(3) - vt(1)
          gamma = vt(1) - vt(2)
          if (abs(alpha*beta*gamma) > xnear) then
            alpha = -(phi(1)*alpha+phi(2)*beta+phi(3)*gamma)/(alpha*beta*gamma)
          else
!
!   FINISH BECAUSE TWO POINTS CALCULATED ARE VERY CLOSE TOGETHER
!
            exit  l180
          end if
          beta = (phi(1)-phi(2))/gamma - alpha*(vt(1)+vt(2))
          if (alpha <= 0.D0) then
            if (phi(right) <= phi(left)) then
              alpha = 3.0D00*vt(right) - 2.0D00*vt(center)
            else
              alpha = 3.0D00*vt(left) - 2.0D00*vt(center)
            end if
            s = alpha - alpold
            if (abs(s) > xmaxm) s = sign(xmaxm,s)*(1 + 0.01D0*(xmaxm/s))
            alpha = s + alpold
          else
            alpha = -beta/(2.0D00*alpha)
            s = alpha - alpold
            xxm = 2.0D00*xmaxm
            if (abs(s) > xxm) s = sign(xxm,s)*(1 + 0.01D0*(xxm/s))
            alpha = s + alpold
          end if
!
!   FINISH IF CALCULATED POINT IS NEAR TO POINT ALREADY CALCULATED
!
          do i = 1, 3
            if (.not.(abs(alpha - vt(i))<delta1*(1.D0+vt(i)) .and. okf)) cycle
            exit  l180
          end do
          i = 1
          if (nvar > 0) then
            xparam(:nvar) = xparef(:nvar) + alpha*pvect(:nvar)
            i = nvar + 1
          end if
          funold = funct
          call compfg (xparam, .TRUE., funct, .TRUE., grad, .FALSE.)
          if (moperr) return
          fmax = dmax1(funct,fmax)
          fmin = min(funct,fmin)
          if (funct < sqstor) call exchng (funct, sqstor, energy, estor, xparam&
            , xstor, alpha, alfs, nvar)
          okf = okf .or. funct<fin
          if (print) then
            write(line,"(5x,'LEFT   ...',f17.8,2f17.7)")vt(left), phi(left), phi(left) - fin
            call to_screen(line)
            write(line,"(5x,'CENTER ...',f17.8,2f17.7)") vt(center), phi(center), phi(center) - fin
            call to_screen(line)
            write(line,"(5x,'RIGHT  ...',f17.8,2f17.7)")vt(right), phi(right), phi(right) - fin
            call to_screen(line)
            if (nvar == 1) then
              if (loc(2,1) == 1) then
                sum = 1.d0
              else
                sum = 180/pi
              end if
              write(line,"(5x,'NEW    ...',f17.8,4f17.7)")alpha, funct, funct - fin, xparam(1)*sum
              write (iw, 130) vt(left), phi(left), phi(left) - fin, &
                              vt(center), phi(center), phi(center) - fin, &
                              vt(right), phi(right), phi(right) - fin, &
                              alpha, funct, funct - fin, xparam(1)*sum
            else
              write(line,"(5x,'NEW    ...',f17.8,2f17.7)")alpha, funct, funct - fin
              write (iw, 130) vt(left), phi(left), phi(left) - fin, &
                              vt(center), phi(center), phi(center) - fin, &
                              vt(right), phi(right), phi(right) - fin, &
                              alpha, funct, funct - fin
            end if
            call to_screen(line)
          end if
  130     format(5x,'LEFT    ...',f17.8,2f17.7,/,5x,'CENTER  ...',f17.8,2f17.7, &
            & /,5x,'RIGHT   ...',f17.8,2f17.7,/,5x,'NEW     ...',f17.8,3f17.7,/)
!
! TEST TO EXIT FROM LINMIN IF NOT DROPPING IN VALUE OF FUNCTION FAST.
!
          if (abs(funold - funct)<delta2 .and. okf) exit  l180
          alpold = alpha
          if (alpha>vt(right) .or. alpha>vt(center) .and. funct<phi(center)&
             .or. alpha>vt(left) .and. alpha<vt(center) .and. funct>phi(center)&
            ) go to 140
!
!   IF SYSTEM IS CREEPING UP ON A DISTANT HIGH VALUE, DISCARD
!   DISTANT HIGH VALUE.
!
          if (nsame > 4) go to 140
          nsame = nsame + 1
          vt(right) = alpha
          phi(right) = funct
          go to 150
  140     continue
          vt(left) = alpha
          nsame = 0
          phi(left) = funct
  150     continue
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
          if (vt(center) < vt(right)) cycle  l180
          i = center
          center = right
          right = i
        end do l180
      end if
  190 continue
      ic = 2
      if (abs(estor - energy) < 1.D-12) ic = 1
      hlast = funct - fin
      drop = sqstor - fin
      call exchng (sqstor, funct, estor, energy, xstor, xparam, alfs, alpha, &
        nvar)
!
!  Check that the density matrix is O.K., or almost O.K.
!
      if (hlast /= 0.D0 .and. hlast > 0.5D0*drop) then
        call compfg (xparam, .TRUE., funct, .TRUE., grad, .FALSE.)
        if (moperr) return
      end if
      okf = funct < ssqlst .or. diis
      if (funct >= ssqlst) return
      if (alpha < 0.D0) then
        alpha = -alpha
        i = 1
        if (nvar > 0) then
          pvect(:nvar) = -pvect(:nvar)
          i = nvar + 1
        end if
      end if
      return
!
!
      end subroutine linmin
