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

      subroutine mpcsyb(chr, kchrge, eionis, dip)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use molkst_C, only : numat, norbs, nclose, nalpha, nbeta, escf, &
      keywrd
      use common_arrays_C, only : coord, eigs
      use chanel_C, only : isyb, syb_fn
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: kchrge
      double precision , intent(in) :: eionis, chr(numat)
      double precision , intent(inout) :: dip
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, i1, i2, nfilled
!-----------------------------------------------
   open(unit = isyb, file = syb_fn)
!  Write out the charge flag and number of atoms
      write (16, '(2I4)', err=30) 1, numat
!  Write out the coordinates and charges
      do i = 1, numat
        write (16, '(4F12.6)', err=30) (coord(j,i),j=1,3), chr(i)
      end do
      nfilled = max(nclose, nalpha, nbeta)
      i1 = max(1,nfilled - 1)
      i2 = min(norbs,nfilled + 2)
!
!  Write out the 2 highest and 2 lowest orbital energies
!
      write (16, 20, err=30) (eigs(j),j=i1,i2), nfilled
   20 format(4f12.6,2x,i4,2x,'HOMOs,LUMOs,# of occupied MOs')
!
!  Write out the Heat of Formation and Ionisation Potential
!
      write (16, '(2F12.6,4X,''HF and IP'')', err=30) escf, eionis
!
!  Write out the Dipole Moment
!
      if (kchrge /= 0) dip = 0.0
      write (16, '(I4,F10.3,''  Charge,Dipole Moment'')', err=30) kchrge, dip
      if (index(keywrd," MULL") /= 0) then
        call mpcpop(1)
      else
        call mpcpop(0)
      end if
      close(unit = isyb, status = "keep")
      return
   30 continue
      write (6, '(A)') 'Error writing SYBYL MOPAC output'
      return
      end subroutine mpcsyb
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      subroutine mpcpop(icok)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE molkst_C, only : numat, keywrd
      use common_arrays_C, only : nfirst, nlast, nat, pb, chrg
      use parameters_C, only : tore
      use chanel_C, only : iw, isyb
      use elemts_C, only : elemnt
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: icok
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, if, il, j, k
      double precision, dimension(numat) :: pop
      double precision :: sum
      double precision, allocatable :: Pot(:), Dipole_correction(:)
      double precision :: Ene
!-----------------------------------------------
!
! This subroutine calculates the total Mulliken populations on the
!   atoms by summing the diagonal elements from the  Mulliken
!   population analysis.
!
      if (icok > -1) write (isyb, '(I4,5X,'' MULLIKEN POPULATION AND CHARGE'')', err=40) icok
!
! ICOK = 1  ==> SYBYL present - do analysis, and write MOPAC and SYBYL output
! ICOK = 0  ==> KEYWORD SYBYL present, but MULLIK absent - don't do analysis
! ICOK = -1 ==> MULLIK present but SYBYL absent - do analysis, but don't write SYBYL output
!
      if (allocated(chrg)) deallocate(chrg)
      allocate(chrg(numat))
      if (icok /= 0) then
        do i = 1, numat
          if = nfirst(i)
          il = nlast(i)
          sum = 0.0
          pop(i) = 0.0
          chrg(i) = 0.0
          do j = if, il
!
!    Diagonal element of mulliken matrix
!
            sum = sum + pb((j*(j+1))/2)
          end do
          k = nat(i)
!
!    Mulliken population for i'th atom
!
          pop(i) = sum
          chrg(i) = tore(k) - pop(i)
        end do
        write (iw, '(3/8X,''MULLIKEN POPULATIONS AND CHARGES'',/)')
        write (iw, '(6X,''NO.  ATOM   POPULATION      CHARGE'')')
        write (iw, '(5X,I4,3X,A2,F13.6,F14.6)') &
         (j, elemnt(nat(j)), pop(j), chrg(j), j = 1, numat)

         if (index(keywrd, " EF") /= 0) then
!
!  Keyword EF is used here to activate the CPE correction.  When CPE is fully integrated, the activation can be removed
!
          allocate(Pot(numat), Dipole_correction(3*numat))
 !         do i = 1, numat !  Use 4 digits after decimal point, to mimic stand-alone code
 !           chrg(i) = nint(chrg(i)*1.d4)/1.d4
 !         end do
          call CPE_energy(Ene, Pot, Dipole_correction)
          write(iw, '(//10x,a,/)') "Chemical-Potential Equalization"
          write(iw, '(10x,a, F20.16, a)')  'Energy                    = ', Ene ,                     "   kcal/mol"
          write(iw, '(10x,a, F20.16, a)')  "Pot(1)                    = ", Pot(1),                   "   eV"
          write(iw, '(10x,a, F20.16, a)')  "Dipole correction(1)      = ", Dipole_correction(1),     "   au"
        end if

        if (icok > -1) write (isyb, "(2f12.6)", err=40) pop(:numat), chrg(:numat)
        call to_screen("To_file: Mulliken")
      end if
      return
   40 continue
      write (iw, '(A)') 'Error writing SYBYL Mulliken population output'
      return
      end subroutine mpcpop
!
