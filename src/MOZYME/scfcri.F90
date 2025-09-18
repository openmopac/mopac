! Molecular Orbital PACkage (MOPAC)
! Copyright 2021 Virginia Polytechnic Institute and State University
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!    http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

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
      precis = (Index (keywrd, " PRECISE") /= 0)
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
