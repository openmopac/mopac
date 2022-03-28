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

subroutine chkion (ox_calc, n_lone_pairs, atom_charge)
   !***********************************************************************
   !
   !   CHKION works out which atoms are ionized.  By the time CHKION is
   !          called, all the obvious elements of the Lewis structure have
   !   been determined (bonds, groups such as -NH3(+), etc.).
   !
   !   An atom is ionized if it has unused atomic orbitals. Using tables of
   !   normal oxidation states (ox_ref_sp) and hypervalent oxidation states
   !   (ox_ref_spd1), the degree of ionization is worked out.  To do this,
   !   the starting oxidation state (ox_calc) is compared to that expected.
   !
   !   This should work for all elements in their common oxidation states.
   !
   !   A limitation is that the sum of atomic number and oxidation state
   !   must be even.  This means that Cr must be Cr(II) and not Cr(III).
   !
   !***********************************************************************
    use chanel_C, only: iw
    use parameters_C, only: tore
    use common_arrays_C, only : nat, nbonds, ibonds, txtatm
    use molkst_C, only : numat, nelecs, line, keywrd
    use elemts_C, only : elemnt
    use atomradii_C, only: is_metal
    use elemts_C, only : cap_elemnt
    use MOZYME_C, only: ions, Lewis_tot, Lewis_elem, &
      iz, ib
    implicit none
    character, intent (in) :: atom_charge(numat)
    integer, intent (inout) :: n_lone_pairs, ox_calc(numat)
    integer :: i, j, k, jj, loop, m, mm, outer_loop, ni, dummy, nocc_local
    integer :: ox_ref(107), ox_ref_sp(107), ox_ref_spd1(107), ox_ref_spd2(107)
    character :: num(10)*1
    double precision, external :: reada
!
!           H           Default formal oxidation states                     He
!           Li Be                                            B  C  N  O  F  Ne
!           Na Mg                                            Al Si P  S  Cl Ar
!           K  Ca Sc            Ti V  Cr Mn Fe Co Ni Cu Zn   Ga Ge As Se Br Kr
!           Rb Sr Y             Zr Nb Mo Tc Ru Rh Pd Ag Cd   In Sn Sb Te I  Xe
!           Cs Ba La Ce-Lu      Hf Ta W  Re Os Ir Pt Au Hg   Tl Pb Bi Po At Rn
!           Fr Ra Ac Th Pa U    Np Pu Am Cm Bk Cf            Cb ++ +  -- -  Tv
!                                      "s" shell
!                              Isolated atom (Most of these will never exist, e.g. C)
    data ox_ref &
        &/ 1,                                                                0, &!    2
        &  1, 2,                                              3, 4, 3,-2,-1, 0, &!   10
        &  1, 2,                                              3, 4, 3,-2,-1, 0, &!   18
        &  1, 2, 3,              4, 3, 2, 3, 2, 3, 2, 1, 2,   3, 4, 3,-2,-1, 0, &!   36
        &  1, 2, 3,              4, 3, 2, 3, 2, 3, 2, 1, 2,   3, 4, 3,-2,-1, 0, &!   54
        &  1, 2, 3, 14*3,        4, 3, 2, 3, 2, 3, 2, 1, 2,   3, 2, 3,-2,-1, 0, &!   86
        &  1, 2, 3, 4, 2, 2,     4, 3, 2, 3, 2, 3, 2, 0, 0,   0, 0, 0, 0, 0, 0 /
!
!                                Atom normal oxidation state.
!  Be very careful in changing these!  They are designed to suit the largest number
!  of common oxidation states.  Also, they affect the sign of the hypervalents.
!  Note in particular, the signs of Groups V, VI, and VII.
!
    data ox_ref_sp &
        &/ 1,                                                                0, &!    2
        &  1, 2,                                              3, 4, 3,-2,-1, 0, &!   10
        &  1, 2,                                              3, 4, 3,-2,-1, 0, &!   18
        &  1, 2, 3,              4, 3, 2, 3, 2, 3, 2, 1, 2,   3, 4, 3,-2, 1, 0, &!   36
        &  1, 2, 3,              4, 3, 2, 3, 2, 3, 2, 1, 2,   3, 4, 3,-2, 1, 0, &!   54
        &  1, 2, 3, 5*0,3*2,6*2, 4, 3, 2, 3, 2, 3, 2, 1, 2,   3, 2, 3,-2, 0, 0, &!   86
        &  1, 2, 3, 4, 2, 2,     4, 3, 2, 3, 2, 3, 2, 0, 0,   0, 0, 0, 0, 0, 0 /
!
!                                Atom lower hypervalent oxidation state
!
    data ox_ref_spd1 &
        &/ 1,                                                                0, &!    2
        &  1, 2,                                              3, 4, 3,-2,-1, 0, &!   10
        &  1, 2,                                              3, 4, 5, 6, 3, 0, &!   18
        &  1, 2, 3,              4, 5, 6, 7, 2, 3, 2, 1, 2,   3, 4, 5, 6, 3, 0, &!   36
        &  1, 2, 3,              4, 5, 6, 7, 2, 3, 2, 1, 2,   3, 4, 5, 6, 3, 0, &!   54
        &  1, 2, 3, 5*0,3*2,6*2, 4, 5, 6, 7, 2, 3, 2, 1, 2,   3, 2, 3, 4, 0, 0, &!   86
        &  1, 2, 3, 4, 2, 2,     4, 5, 6, 7, 2, 3, 2, 0, 0,   0, 0, 0, 0, 0, 0 /
!
!                                Atom higher hypervalent oxidation state
!
    data ox_ref_spd2 &
        &/ 1,                                                                0, &!    2
        &  1, 2,                                              3, 4,-3, 2,-1, 0, &!   10
        &  1, 2,                                              3, 4,-5, 6, 7, 0, &!   18
        &  1, 2, 2,              4, 3, 2, 3, 2, 3, 2, 1, 2,   3, 4,-5, 6, 5, 0, &!   36
        &  1, 2, 2,              4, 3, 2, 3, 2, 3, 2, 1, 2,   3, 4,-5, 6, 7, 0, &!   54
        &  1, 2, 2, 5*0,3*2,6*2, 4, 3, 2, 3, 2, 3, 2, 1, 2,   3, 4,-3, 2, 0, 0, &!   86
        &  1, 2, 2, 4, 2, 2,     4, 3, 2, 3, 2, 3, 2, 0, 0,   0, 0, 0, 0, 0, 0 /
!
!  If instructed to, change oxidation states
!
    i = index(keywrd," METAL")
    if (i /= 0) then
      j = index(keywrd(i:), ") ") + i
      if (i /= j) then
        line = keywrd(i + 7:j)
        do k = 1, 83
          if (is_metal(k)) then
            if (cap_elemnt(k)(2:2) == " ") then
              jj = 1
            else
              jj = 2
            end if
            m = index(line,cap_elemnt(k)(:jj))
            if (m /= 0) then
              jj = jj + m
              if (line(jj:jj) == "(") then
              m = nint(reada(line, m))
              if (mod(m - ox_ref(k), 2) == 1) then
                write(iw,'(10x, a, i3)')"Formal oxidation state of "//elemnt(k)//" is:   ", ox_ref(k)
                write(iw,'(10x,a, i3)')"Requested oxidation state of "//elemnt(k)//" is:", m
                if (mod(ox_ref(k), 2) == 1) then
                write(iw,'(10x,a, i3)')"Alternate oxidation states must be odd, e.g, -3, -1, +1, +3, +5, etc."
                else
                  write(iw,'(10x,a, i3)')"Alternate oxidation states must be even, e.g, -4, -2, +0, +2, +4, etc."
                end if
                call mopend("OXIDATION STATE REQUESTED IS NOT ALLOWED")
                return
              end if
              ox_ref(k) = m
              ox_ref_spd1(k) = ox_ref(k)
              end if
            end if
          end if
        end do
      end if
    end if
!
!  ALL ORBITALS SHOULD HAVE BEEN USED BY THIS TIME.  SOME ORBITALS
!  HAVE NOT BEEN USED: CHECK TO SEE IF THEY SHOULD BE FILLED OR
!  BE EMPTY.
!
!  First, work out the sign of the oxidation state of the atoms
!
    do i = 1, numat
      k = ions(i)
      if (nbonds(i) == 0) then
        jj = sign(1, ox_ref(nat(i)))    ! Isolated atom or ion
      else
        jj = sign(1, ox_ref_sp(nat(i))) ! Lowest valence state
      end if
      do j = 1, Lewis_tot
        if (Lewis_elem(1,j) /= 0 .and. Lewis_elem(2,j) /= 0) then
          if (Lewis_elem(1,j) == i) k = k + jj
          if (Lewis_elem(2,j) == i) k = k + jj
        end if
      end do
      ox_calc(i) = k
      if (ib(i) == 0) then
        ions(i) = iz(i)
        iz(i) = 0
      end if
    end do
!
!   Write in MOZYME(Cr(IV))
!

    do outer_loop = 1, 2
!
!  The first cycle of outer_loop assignes all obvious charges,
!  the second cycle assignes the more subtle charges
!  such as those on C2H5
!
      atom_loop: do i = 1, numat
        loop = 0
        do
          if (ib(i) <= 0) cycle atom_loop
          loop = loop + 1
          if (loop > 10) goto 1040
          if (nat(i) == 6) then
  !
  !  Carbon is particularly difficult: decide its charge based on its environment
  !
            if (outer_loop == 1) then
              if (iz(i) > 1) then
                call add_Lewis_element(i,0,0, n_lone_pairs)  !  Lone pair
                if (iz(i) == 0 .and. ib(i) == 1) then
                  call add_Lewis_element(0,i,0, dummy)  ! Virtual lone pair
                end if
              else
                do m = 1, nbonds(i)
                  if (ions(ibonds(m,i)) > 0) then ! Adjacent atom is a cation
                    call add_Lewis_element(i,0,-1, n_lone_pairs)    ! Lone pair
                    cycle atom_loop
                  else if (ions(ibonds(m,i)) < 0) then ! Adjacent atom is an anion
                    call add_Lewis_element(0,i,1, dummy)    ! Cation
                    cycle atom_loop
                  end if
                end do
  !
  !     If there is a nearby atom of Group V, make the atom a cation
  !
                do m = 1, nbonds(i)
                  mm = ibonds(m, i)
                     !
                     !   IS ANY ATOM ATTACHED TO THE ATOM OF GROUP 5?
                     !
                  if (Nint(tore(nat(mm))) == 5) go to 1010
                  if (Nint(tore(nat(mm))) == 6) go to 1012
                  do j = 1, nbonds(mm)
                    jj = ibonds(j, mm)
                        !
                        !   IS ANY ATOM ATTACHED TO ANY ATOM ATTACHED
                        !   TO THE ATOM OF GROUP 5?
                        !
                    if (Nint(tore(nat(jj))) == 5) go to 1010
                  end do
                end do
                goto 99
  1010        call add_Lewis_element(0,i,1, dummy)    ! Cation
                cycle atom_loop
  1012        call add_Lewis_element(i,0,-1, n_lone_pairs)    ! Lone pair
                cycle atom_loop
      !
      !     If there is a nearby atom of Group VI, make the atom an anion
      !
    99          do m = 1, nbonds(i)
                  mm = ibonds(m, i)
                     !
                     !   IS ANY ATOM ATTACHED TO THE ATOM OF GROUP 5?
                     !
                  if (Nint(tore(nat(mm))) == 6) go to 1011
                  do j = 1, nbonds(mm)
                    jj = ibonds(j, mm)
                        !
                        !   IS ANY ATOM ATTACHED TO ANY ATOM ATTACHED
                        !   TO THE ATOM OF GROUP 6?
                        !
                    if (Nint(tore(nat(jj))) == 6) go to 1011
                  end do
                end do
                cycle atom_loop
  1011          call add_Lewis_element(i,0, -1, n_lone_pairs)   ! Lone pair
                cycle atom_loop
              end if
            else
  !
  !  To get here, there are no obvious indications of the charge on
  !  carbon, so postpone assigning the charge until after all other
  !  charges are known.
  !
              call add_Lewis_element(-i,0, 0, dummy) ! Anion or cation - decide later
              cycle atom_loop
            end if
          else  !  Atom is not a carbon - use generic rules
            if (loop == 1) then
!
!  Adjust oxidation state to match that expected.
!
              ni = nat(i)
              do
                if (atom_charge(i) == " ") then
                  if (Abs(ox_calc(i)) > Abs(ox_ref_spd1(ni))) then
                    jj = ox_ref_spd2(ni) - ox_calc(i)
                  else if (Abs(ox_calc(i)) > Abs(ox_ref_sp(ni))) then
                    jj = ox_ref_spd1(ni) - ox_calc(i)
                  else if (nbonds(i) > 0) then
                    jj = ox_ref_sp(ni) - ox_calc(i)
                  else
                    jj = ox_ref(ni) - ox_calc(i)
                  end if
                else if (atom_charge(i) == "+") then
                  jj = 1
                else if (atom_charge(i) == "-") then
                  jj = -1
                else if (atom_charge(i) == "0") then
                  jj = 0
                end if
!
!  jj is the presumed charge on the atom.
!  "presumed" = difference between calculated oxidation state and default oxidation state,
!               or charge specified by the user.
!
                if (jj == 0) exit
                if (jj >= 2 .and. nat(i) == 7) then
                  jj = -2 ! Azide
                end if
                if (ib(i) == 0) exit
                  if (jj > 1) then
                  call add_Lewis_element(0,i,2, dummy)  ! di-cation
                  ox_calc(i) = ox_calc(i) + 2
                 else if (jj == 1) then
                  call add_Lewis_element(0,i,1, dummy)  ! Cation
                   ox_calc(i) = ox_calc(i) + 1
                 else if (jj == -1) then
                  call add_Lewis_element(i,0,-1, n_lone_pairs) ! Anion
                  ox_calc(i) = ox_calc(i) - 1
                 else if (jj < -1) then
                  call add_Lewis_element(i,0,-2, n_lone_pairs) ! di-Anion
                  ox_calc(i) = ox_calc(i) - 2
                end if
              end do
            else  if (iz(i) > 1) then
              call add_Lewis_element(i,0,0, n_lone_pairs)  ! Lone pair
            else
              call add_Lewis_element(0,i,0, dummy)  ! Virtual lone pair
            end if
          end if
        end do
      end do atom_loop
    end do
    nocc_local = 0
    do i = 1, Lewis_tot
      if (Lewis_elem(1,i) /= 0) nocc_local = nocc_local + 1
    end do
    mm = nocc_local - nelecs/2
!
!  mm is the number of undetermined Lewis elements that should be
!  virtual.  These elements are always charged carbon atoms.
!
    do i = 1, mm
      do j = 1, Lewis_tot
        jj = Lewis_elem(1,j)
        if (jj < 0) then
          jj = -jj              ! Convert undetermined charged carbon
          Lewis_elem(2,j) = jj   ! into a cation
          Lewis_elem(1,j) = 0
          iz(jj) = iz(jj) - 1
          if (ib(jj) == 1) then
            write (iw,'(/10x,a,i6,a)') " WARNING: Carbon atom ", jj, " has an unexpected valency!   PDB Label: """// &
            txtatm(jj)//""""
            if (nbonds(jj) > 0) then
              num = "1"
              do k = 1, nbonds(jj)
                num(k) = char(ichar("1") +int(log10(ibonds(k,jj) + 0.05)))
              end do
              line = '(10x,a,a,i'//num(1)//',3x,a,i'//num(2)//',3x,a,i'//num(3)
              line = trim(line)//',3x,a,i'//num(4)//',3x,a,i'//num(5)//',3x,a,i'//num(6)
              line = trim(line)//',3x,a,i'//num(7)//',3x,a,i'//num(8)//',3x,a,i'//num(9)//',3x,a,i'//num(10)//',a)'
              write (iw,trim(line))" It is connected to:  ",(elemnt(nat(ibonds(k,jj))),ibonds(k,jj), k = 1, nbonds(jj)), "only"
            else
              write (iw,'(a)')" It is not connected to any other atom."
            end if
            write(iw,"(/)")
          end if
!
! Convert all remaining unused atomic orbitals into virtual lone pairs
!
          do
            if (ib(jj) <= 0) exit
            call add_Lewis_element(0,jj,0, dummy)  ! Virtual lone pair
          end do
          exit
        end if
      end do
    end do
!
!   Convert all remaining undetermined charged carbon atoms into anions
!
    do j = 1, Lewis_tot
      jj = Lewis_elem(1,j)
      if (jj < 0) then
        jj = -jj
        Lewis_elem(1,j) = jj
        n_lone_pairs = n_lone_pairs + 1
        if (iz(jj) == 1) iz(jj) = 0
        if (ib(jj) == 1) then
!
! create a virtual lone pair
!
          write (iw,'(/10x,a,i6,a)') " WARNING: Carbon atom ", jj, " has an unexpected valency!   PDB Label: """// &
            txtatm(jj)//""""
          if (nbonds(jj) > 0) then
            num = "1"
            do k = 1, nbonds(jj)
              num(k) = char(ichar("1") +int(log10(ibonds(k,jj) + 0.05)))
            end do
            line = '(10x,a,a,i'//num(1)//',3x,a,i'//num(2)//',3x,a,i'//num(3)//',3x,a,i'//num(4)//',3x,a,i'
            line = trim(line)//num(5)//',3x,a,i'//num(6)//',3x,a,i'//num(7)//',3x,a,i'//num(8)
            line = trim(line)//',3x,a,i'//num(9)//',3x,a,i'//num(10)//',a)'
            write (iw,trim(line))" It is connected to:  ",(elemnt(nat(ibonds(k,jj))),ibonds(k,jj), k = 1, nbonds(jj)), "only"
          else
            write (iw,'(a)')" It is not connected to any other atom."
          end if
          write(iw,"(/)")
          call add_Lewis_element(0,jj,0, dummy)  ! Virtual lone pair
          ib(jj) = 0
        end if
      end if
    end do

    return
1040 write (iw,*) " CHARGE ON ATOM", i, " UNREASONABLE"
    call mopend ("Charge on an atom is unreasonable")
end subroutine chkion
