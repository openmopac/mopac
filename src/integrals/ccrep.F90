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

      subroutine ccrep(ni, nj, r, enuclr, gab)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use parameters_C, only : alp, tore, guess1, guess2, guess3, alpb, xfac, &
        par1, par2, par3, par4
      use funcon_C, only : a0
      use molkst_C, only : method_mndod, method_pm6, method_am1, &
        method_pm7, method_pm8, method_pm6_org
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
        if (method_pm6_org) then
          call ccrep_PM6_ORG(ni, nj, r, fff, abond, scale)
        else if (.not. method_mndod) then
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
        if (method_pm6 .or. method_pm7 .or. method_pm8 .or. method_pm6_org) then
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
      if (method_pm6 .or. method_pm7 .or. method_pm8 .or. method_pm6_org) then
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
      ! AM1 B-H, B-C, & B-halogen corrections
      else if (method_am1 .and. (ni == 5 .or. nj == 5) .and. &
          (ni == 1 .or. nj == 1 .or. ni == 6 .or. nj == 6 .or. ni == 9 .or. nj == 9 .or. &
           ni == 17 .or. nj == 17 .or. ni == 35 .or. nj == 35 .or. ni == 53 .or. nj == 53)) then
          if (ni == 1 .or. nj == 1) then
            ax = 10.0d0*(r - 0.832586d0)**2
            if (ax <= 25.D0) scale = scale + tore(ni)*tore(nj)/r*0.412253d0*exp(-ax)
            ax = 6.0d0*(r - 1.186220d0)**2
            if (ax <= 25.D0) scale = scale + tore(ni)*tore(nj)/r*(-0.149917d0)*exp(-ax)
          else if (ni == 6 .or. nj == 6) then
            ax = 8.0d0*(r - 1.063995d0)**2
            if (ax <= 25.D0) scale = scale + tore(ni)*tore(nj)/r*0.261751d0*exp(-ax)
            ax = 5.0d0*(r - 1.936492d0)**2
            if (ax <= 25.D0) scale = scale + tore(ni)*tore(nj)/r*0.050275d0*exp(-ax)
          else
            ax = 9.0d0*(r - 0.819351d0)**2
            if (ax <= 25.D0) scale = scale + tore(ni)*tore(nj)/r*0.359244d0*exp(-ax)
            ax = 9.0d0*(r - 1.574414d0)**2
            if (ax <= 25.D0) scale = scale + tore(ni)*tore(nj)/r*0.074729d0*exp(-ax)
          end if
          if (ni == 5) then
            do ig = 1, 4
              ax = guess2(nj,ig)*(r - guess3(nj,ig))**2
              if (ax <= 25.D0) scale = scale + tore(ni)*tore(nj)/r*guess1(nj,ig)*exp(-ax)
            end do
          else
            do ig = 1, 4
              ax = guess2(ni,ig)*(r - guess3(ni,ig))**2
              if (ax <= 25.D0) scale = scale + tore(ni)*tore(nj)/r*guess1(ni,ig)*exp(-ax)
            end do
          end if
          i = 0
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
      if (method_pm6 .or. method_pm7 .or. method_pm8 .or. method_pm6_org) then
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
!
  subroutine ccrep_PM6_ORG(ni, nj, r, fff, abond, scale)
!
!  Core-core scaling factor for method PM6-ORG
!
  use parameters_C, only : par1, par2, par3, par4, par5, par10, par11, par12, par13, par14, par15, &
    par16, par17, par18, par19, par20, par21, par22, par23, par24, par25, par26, par27, par28, &
    par29, par30, par31, par32, par33, par34, par35, par36, par37, par38, par39, par40, par41, par42, par43, &
    par44, par45, par46
  implicit none
  integer, intent (in) :: ni, nj
  double precision, intent (in) :: r, fff, abond
  double precision, intent (out) :: scale
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: i, j
  double precision :: sum = 0.01d0
  save :: sum
  scale = 1.0d0 + 2.d0 * fff * Exp (-abond*(r + 0.0003*r**6)) 
!
!   Put all special handling codes here.
!   
  i = max(ni, nj)
  j = min(ni, nj)
  select case (j) 
!
!  Hydrogen with everything
!
  case (1)             
    select case (i)
    case (1)
!
! Steric repulsion term for H - H interaction. Uses: PAR16, PAR17, and PAR18
!
      if (r - par18 > 0.d0) then
        scale = scale + sum*par16*exp(-par17*(r-par18)**2) 
      else
        scale = scale + sum*par16
      end if 
    case (6)
!
! Steric repulsion term for C - H interaction. Uses: PAR19, PAR11, and PAR12
!
      if (r - par12 > 0.d0) then
        scale = 1.0d0 + 2.d0 * fff * Exp (-abond*r**2) + sum*par19*exp(-par11*(r-par12)**2) 
      else
        scale = 1.0d0 + 2.d0 * fff * Exp (-abond*r**2) + sum*par19
      end if 
    case (7)   
!
! Steric repulsion term for N - H interaction. Uses: PAR38, PAR39, and PAR40
!
      if (r - par40 > 0.d0) then
        scale = 1.0d0 + 2.d0 * fff * Exp (-abond*r**2) + sum*par38*exp(-par39*(r-par40)**2) 
      else
        scale = 1.0d0 + 2.d0 * fff * Exp (-abond*r**2) + sum*par38
      end if 
    case (8)
!
! Steric repulsion term for O - H interaction. Uses: PAR3, PAR4, and PAR5
!
      if (r - par5 > 0.d0) then
        scale = 1.0d0 + 2.d0 * fff * Exp (-abond*r**2) + sum*par3*exp(-par4*(r-par5)**2) 
      else
        scale = 1.0d0 + 2.d0 * fff * Exp (-abond*r**2) + sum*par3
      end if      
    case (16)
!
! Steric term for S - H interaction. Uses: PAR44, PAR45, and PAR46
!
     if (r - par46 > 0.d0) then
        scale = scale + sum*par44*exp(-par45*(r-par46)**2) 
      else
        scale = scale + sum*par44
      end if  
    end select
!
!  Carbon 
!
  case (6)
    select case (i)
    case (6)
      scale = scale + par1*exp(-par2*r) ! To correct C-C triple bond HoF 
!
! Steric repulsion term for C - C interaction. Uses: PAR13, PAR14, and PAR15
!
      if (r - par15 > 0.d0) then
        scale = scale + sum*par13*exp(-par14*(r-par15)**2) 
      else
        scale = scale + sum*par13
      end if 
    case (7)
!
! Steric repulsion term for N - C interaction. Uses: PAR35, PAR36, and PAR37
!      
      if (r - par37 > 0.d0) then
        scale = scale + sum*par35*exp(-par36*(r-par37)**2) 
      else
        scale = scale + sum*par35
      end if      
    case (8)
!
! Steric repulsion term for C - O interaction. Uses: PAR20, PAR21, and PAR22
!
      if (r - par22 > 0.d0) then
        scale = scale + sum*par20*exp(-par21*(r-par22)**2) 
      else
        scale = scale + sum*par20
      end if    
    case (16)
!
! Steric repulsion term for S - C interaction. Uses: PAR29, PAR30, and PAR31
!
      if (r - par31 > 0.d0) then
        scale = scale + sum*par29*exp(-par30*(r-par31)**2) 
      else
        scale = scale + sum*par29
      end if      
    end select 
!
!  Nitrogen 
!
  case (7) 
    select case (i)
    case (8)
!
! Steric repulsion term for N - O interaction. Uses: PAR26, PAR27, and PAR28
!      
      if (r - par28 > 0.d0) then
        scale = scale + sum*par26*exp(-par27*(r-par28)**2) 
      else
        scale = scale + sum*par26
      end if      
   case (16)
!
! Steric repulsion term for S - N interaction. Uses: PAR41, PAR42, and PAR43
!
      if (r - par43 > 0.d0) then
        scale = scale + sum*par41*exp(-par42*(r-par43)**2) 
      else
        scale = scale + sum*par41
      end if
    end select 
!
!  Oxygen 
!
  case (8)
    select case (i)
    case (8)
!
! Steric repulsion term for O - O interaction. Uses: PAR32, PAR33, and PAR34
!
      if (r - par34 > 0.d0) then
        scale = scale + sum*par32*exp(-par33*(r-par34)**2) 
      else
        scale = scale + sum*par32
      end if     
    case (14)
!
!  For R0 = 3.6 Angstroms, energy in kcal/mol per Si-O interaction is (11952/r = 3320) times premultiplier.
!
      scale = scale  -0.7d-3*exp(-(r - 2.9d0)**2) ! To correct Si-O weak long-range interaction 
    case (16)
!
! Steric repulsion term for S - O interaction. Uses: PAR23, PAR24, and PAR25
!
      if (r - par25 > 0.d0) then
        scale = scale + sum*par23*exp(-par24*(r-par25)**2) 
      else
        scale = scale + sum*par23
      end if      
    end select 
  end select
  return
  end subroutine ccrep_PM6_ORG
