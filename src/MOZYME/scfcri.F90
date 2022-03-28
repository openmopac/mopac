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

subroutine scfcri (selcon)
    use molkst_C, only: numcal, keywrd, efield
    use chanel_C, only: iw
 !   use common_convrg, only: scfref
   !
   !.. Implicit Declarations ..
    implicit none
   !
   !.. Formal Arguments ..
    double precision, intent (inout) :: selcon
   !
   !.. Local Scalars ..
    logical :: precis
    integer :: i
    integer :: icalcn = 0
    double precision, save  :: scfcrt, scfref
    double precision, external :: reada
    if (icalcn /= numcal) then
      icalcn = numcal
      !
      !   DETERMINE THE SELF-CONSISTENCY CRITERION
      !
      ! SCFCRT IS MACHINE-PRECISION DEPENDENT
      !
      !   IF SCFCRT IS CHANGED, THEN ALSO CHANGE DEFAULT IN WRTKEY
      !
      scfcrt = 1.d-2
      !
      !  THE USER CAN STATE THE SCF CRITERION, IF DESIRED.
      !
      i = Index (keywrd, " TS") + Index (keywrd, " FORCETS") + Index (keywrd, " IRC=")
      if (i /= 0) scfcrt = 2.0d-3
      precis = (Index (keywrd, " PRECIS") /= 0)
      i = Index (keywrd, " RELSCF")
      if (i /= 0) then
        scfcrt = reada (keywrd, i) * scfcrt
        write (iw, "('  SCF CRITERION =',G14.4)") scfcrt
        if (scfcrt < 1.d-5) then
10000     format (//2 x, " THERE IS A RISK OF INFINITE LOOPING WITH", &
         & " THE SCFCRT LESS THAN 1.D-5")
          write (iw, 10000)
        end if
      end if
      i = Index (keywrd, " SCFCRT")
      if (i /= 0) then
        scfcrt = reada (keywrd, i)
        write (iw, "('  SCF CRITERION =',G14.4)") scfcrt
        if (scfcrt < 1.d-5) then
          write (iw, "(//2x,' THERE IS A RISK OF INFINITE LOOPING WITH', &
         & ' THE SCFCRT LESS THAN 1.D-5')")
        end if
      end if
      !
      !   SELF-CONSISTENCY CRITERIA: SELCON IS IN KCAL/MOL
      !   IF GNORM IS LARGE, MAKE SELCON BIGGER
      !
      if (precis) then
        scfcrt = scfcrt * 0.01d0
      end if
      !
      ! For polarizability calculations, the default convergence should be
      ! tightened
      !
      if (Index (keywrd, " POLAR") /= 0 .and. scfref == 0.0d0) then
        scfcrt = 1.d-4
      end if
      selcon = scfcrt
      scfref = scfcrt
    else
      !
      !  IF POLARIZATION IS BEING CALCULATED, TIGHTEN SCF CRITERION
      !
      if (Abs (efield(1))+Abs(efield(2))+Abs(efield(3)) > 1.d-6) then
        selcon = 1.d-4
      end if
      if (scfref /= 0.d0) then
        selcon = scfref
      end if
    end if
    return
end subroutine scfcri
