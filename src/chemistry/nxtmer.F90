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

  subroutine nxtmer (iatom, nbackb)
    use common_arrays_C, only : nat, ibonds, nbonds
    implicit none
    integer, intent (in) :: iatom
    integer, dimension (4), intent (inout) :: nbackb
!
    integer :: i, ii, iii, j, jatom, jj, jjofco = 0, jofco = 0, k, kk, l, ihcr, ico, ico1
    ico = 0
    ico1 = 0
    ihcr = 0
!
    jatom = nbackb(4)
    do i = 1, nbonds(iatom)
       if (nat(ibonds(i, iatom)) == 6) ihcr = ibonds(i, iatom)
    end do
   !
   !   START WITH NITROGEN
   !
    outer_loop: do i = 1, nbonds(iatom)
      !
      !    IS THERE A CARBON ATTACHED TO THE NITROGEN?
      !
      if (nat(ibonds(i, iatom)) == 6) then
        l = ibonds(i, iatom)
        jj = 0
        jjofco = 0
        do k = 1, nbonds(l)
            !
            !    IS THERE A CARBON ATTACHED TO THE CARBON ATTACHED TO THE NITROGEN?
            !
          if (nat(ibonds(k, l)) == 6) then
            j = ibonds(k, l)
            do ii = 1, nbonds(j)
!
!   IS THERE AN OXYGEN ATOM THAT IS BONDED ONLY TO THE CARBON ATTACHED TO THE
!   CARBON ATTACHED TO THE NITROGEN?
!
              if (nat(ibonds(ii, j)) == 8 .and. nbonds(ibonds(ii, j)) == 1) go to 1000
            end do
            cycle
1000        jofco = ibonds(ii, j)
            ico1 = ico
            ico = j
            iii = ii
            ihcr = l
            do ii = 1, nbonds(j)
                  !
                  !   IS THERE A NITROGEN ATTACHED TO THE
                  !   CARBON ATTACHED TO THE CARBON ATTACHED TO THE NITROGEN?
                  !
              if (nat(ibonds(ii, j)) == 7) go to 1010
            end do
            kk = 0
            do ii = 1, nbonds(j)
                  !
                  !   ARE THERE TWO OXYGEN ATOMS ATTACHED TO THE
                  !   CARBON ATTACHED TO THE CARBON ATTACHED TO THE NITROGEN?
                  !
              if (nat(ibonds(ii, j)) == 8) kk = kk + 1
            end do
            if (kk == 2) then
!
! A -COO group detected
!
              jatom = 0
              exit outer_loop
            end if
               !
               !   NO BOND TO NITROGEN, THEREFORE NOT BACKBONE
               !
            kk = 0
            do ii = 1, nbonds(j)
              if (nat(ibonds(ii, j)) == 6) then
                kk = kk + 1
              end if
            end do
               !
               !  IF THERE IS MORE THAN ONE CARBON ATOM ATTACHED TO THE
               !  SECOND CARBON ATOM THEN IT CANNOT BE A BACKBONE ATOM
               !
            if (kk == 1) then
              jj = j
              jjofco = ibonds(iii, j)
            end if
            cycle
1010        jatom = ibonds(ii, j)
            do ii = 1, nbonds(jatom)
              if (nat(ibonds(ii, jatom)) == 6) exit outer_loop
            end do
          end if
        end do
         !
         !  IN ORDER TO GET HERE, THE BACKBONE HAS NOT BEEN RECOGNIZED,
         !  THEREFORE SET CARBON OF C=O TO JJ
         !
        if (jj /= 0) then
          j = jj
        end if
        if (jjofco /= 0) then
          jofco = jjofco
        end if
      end if
    end do outer_loop
   !
   !    IATOM = STARTING NITROGEN
   !    L     = CARBON OF -(HCR)-
   !    J     = CARBON OF -(C=O)-
   !    JATOM = NITROGEN OF -(C=O)-NH-
   !    JOFCO = OXYGEN OF -(C=O)-
    if (ico1 /= 0) then
!
! Distinguish a terminal threonine with "CO" instead of "COOH"
!
      do i = 1, nbonds(ico1)
        if (ibonds(i,ico1) == jofco) ico = ico1
      end do
    end if
    nbackb(1) = ihcr
    nbackb(2) = ico
    nbackb(3) = jofco
    nbackb(4) = jatom
    return
  end subroutine nxtmer
