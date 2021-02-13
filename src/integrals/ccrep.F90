       subroutine ccrep(ni, nj, r, enuclr, gab) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      use parameters_C, only : alp, tore, guess1, guess2, guess3, alpb, xfac, &
        par1, par2, par3, par4
      use funcon_C, only : a0
      use molkst_C, only : method_mndod, method_pm6, method_am1, &
        method_pm7, method_PM8
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: ni 
      integer , intent(in) :: nj 
      double precision , intent(inout) :: r 
      double precision , intent(in) :: gab
      double precision , intent(out) :: enuclr 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nt, ig, i, j
    !  logical :: opend
      double precision :: alpni, alpnj, enuc, abond, fff, scale, eni, enj, ax
!-----------------------------------------------
!     CONVERT TO ANGSTROM AND INITIALIZE VARIABLES.
      r = r*a0 
      alpni = alp(ni) 
      alpnj = alp(nj) 
!
!     CALCULATE REPULSIVE TERM. 
!
      enuc = tore(ni)*tore(nj)*gab 
      
!
! Get bond parameters if defined.
!
      if (ni < 101 .and. nj < 101) then
        fff = xfac(ni,nj)        
      else
        fff = 0.d0
      end if
      if (method_PM7) then
        if (abs(fff) < 1.d-5) then
          if (ni < 99 .and. nj < 99) then
!
!  This is a poor guess for the value of alpb(ni,nj) and xfac(ni,nj)
!
            xfac(ni,nj) = 0.5d0*(xfac(ni,ni) + xfac(nj,nj))
            alpb(ni,nj) = 0.5d0*(alpb(ni,ni) + alpb(nj,nj))
            fff = xfac(ni,nj)
          else
            fff = 0.d0  !  unreal atoms - set to an arbitary value
          end if
        end if
      end if
      if (abs(fff) > 1.d-5) then
     !
     ! Bond parameters defined
     !
        abond = alpb(ni,nj)
        if (abond < 1.d-6) abond = 1.2d0
        if (.not. method_mndod) then
          if (method_pm6 .or. method_pm7 .or. method_pm8) then
            scale = 1.0d0 + 2.d0 * fff * Exp (-abond*(r + 0.0003*r**6)) 
  !
  !   Put all special handling codes here.  To date, these are:
  !
  !    If an O-H, C-H, or N-H interaction, use a different core-core term
  !   
            i = max(ni, nj)
            j = min(ni, nj)
            select case (j) 
            case (1)             
              select case (i)
              case (1)
              case (6:7)
                scale = 1.0d0 + 2.d0 * fff * Exp (-abond*r**2) 
              case (8)
!
! Slow O - H term
!
                scale = 1.0d0 + 2.d0 * fff * Exp (-abond*r**2) - par3*exp(-par4*r*2) ! par3*exp(-par4*r*2) used by PM7
              end select 
            case (6)
              select case (i)
              case (6)
                scale = scale + par1*exp(-par2*r) ! To correct C-C triple bond HoF 
              end select 
            case (7) 
              select case (i)
              case (7)
              end select 
            case (8)
             select case (i)
             case (14)
!
!  The maximum gradient is, to a first approximation, given by: exp(-ax*(r - R0 + 1/(2.d0*ax))**2) where
!  R0 is the defined distance.
!  For R0 = 3.6 Angstroms, energy in kcal/mol per Si-O interaction is (11952/r = 3320) times premultiplier.
!
                ax = 1.0d0
                scale = scale  -0.7d-3*exp(-(r - 2.9d0)**2) ! To correct Si-O weak long-range interaction 
              end select
            end select
          else  ! Not PM6
            if (method_am1 .and. (ni == 42 .and. nj == 1 .or. ni == 1 .and. nj == 42)) then
!
!  AM1-d Mo-H interaction
!
                scale = 1.0d0 + r * 2.d0 * fff * Exp (-abond*r)
            else
!
!  All other interactions
!
              scale = 1.0d0 + 2.d0 * fff * Exp (-abond*r)
            end if
          end if
        else
         !
         ! MNDO/d interactions
         !
         ! Alpb(ip) is used as alpha parameter for one element, and
         ! the standard alpha parameter is used for the other element.
          if (ni == nj) then
            scale = 1.0d0 + 2.d0 * Exp (-abond*r)
          else
            select case (nj)
            case(11, 12, 13)
              scale = 1.0d0 + Exp (-abond*r) + Exp (-alp(ni)*r)
            case default
              scale = 1.0d0 + Exp (-abond*r) + Exp (-alp(nj)*r)
            end select
          end if
        end if
       !
       ! Multiply monopole term by scaling factor
       !
        enuclr = enuc * scale
      else 
        abond = 0.d0
        if (method_pm6 .or. method_pm7 .or. method_pm8) then
          if (ni .gt. 56 .and. ni  .lt. 72 &
         .or. nj .gt. 56 .and. nj  .lt. 72 )then
            scale = 10.d0*exp((-3.d0*r))
          else
            scale = 10.d0*exp((-2.18d0*r)) ! This is a generic core-core term.
          end if
          eni = 0.d0
          enj = 0.d0
        else
          eni = exp((-alpni*r)) 
          enj = exp((-alpnj*r)) 
          scale = eni + enj 
        end if
!
!  This is almost certainly dead code.
        nt = ni + nj 
        if (nt == 8 .or. nt == 9) then 
          if (ni == 7 .or. ni == 8) scale = scale + (r - 1.D0)*eni 
          if (nj == 7 .or. nj == 8) scale = scale + (r - 1.D0)*enj 
        end if 
!
!  End of probable dead code
!
        enuclr = abs(scale*enuc) + enuc
      end if
      scale = 0.d0 
      if (method_pm6 .or. method_pm7 .or. method_pm8) then
  !
  !  VdW term
  !
        ax = guess2(ni,1)*(r - guess3(ni,1))**2 
        if (ax < 25.D0) scale = scale + tore(ni)*tore(nj)/r*guess1(ni,1)*exp((-ax)) 
        ax = guess2(nj,1)*(r - guess3(nj,1))**2 
        if (ax < 25.D0) scale = scale + tore(ni)*tore(nj)/r*guess1(nj,1)*exp((-ax))
        if (abond > 1.d-4) then
          i = 0
        else
          i = 4
        end if
      else
          i = 4
          if (fff > 1.d-4) i = 0
      end if
      do ig = 1, i 
        if (abs(guess1(ni,ig)) > 0.D0) then 
          ax = guess2(ni,ig)*(r - guess3(ni,ig))**2 
          if (ax <= 25.D0) scale = scale + tore(ni)*tore(nj)/r*guess1(ni,ig)*&
            exp((-ax)) 
        end if 
        if (abs(guess1(nj,ig)) <= 0.D0) cycle  
        ax = guess2(nj,ig)*(r - guess3(nj,ig))**2 
        if (ax > 25.D0) cycle  
        scale = scale + tore(ni)*tore(nj)/r*guess1(nj,ig)*exp((-ax)) 
      end do  
      enuclr = enuclr + scale 
      if (method_pm6 .or. method_pm7 .or. method_pm8) then
!
!  The next term is the unpolarizable core - unpolarizable core interaction
!  It should have no effect on heats of formation or other properties.
!  Its purpose is to allow systems where the atoms are forced together
!  to be realistic.  The form of the expression is the Lennard-Jones "12" part
!  in the "12 - 6" potential. 
!
!  The multiplier "1.d-8" was chosen to make the function negligible at
!  normal bonding distances. "ax" is proportional to the interatomic distance 
!  divided by the sum of the two covalent radii.
!
        ax = r/(ni**0.3333d0 + nj**0.3333d0)
        if (ax < 3.d0) then
          scale = 1.d-8/ax**12
          enuclr = enuclr + min(scale, 1.d5)
        end if
      end if  
      return  
      end subroutine ccrep 
