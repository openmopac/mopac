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

  subroutine wrt_diffs
!
! Given two geometries, work out the differences in bond-lengths
! and print the resulting values, largest first
!
! The two geometries are in geo and geoa
!
   use chanel_C, only : iw, job_fn
   use molkst_C, only : numat, refkey, line, geo_dat_name, maxtxt
!
   use elemts_C, only : elemnt
!
   use common_arrays_C, only : geo, geoa, nat, nbonds, ibonds, txtatm, coord
!
   implicit none
   integer :: i, j, j1, k, l, m, n, npairs, ijbonds(15)
   integer, allocatable :: pairs(:,:)
   logical :: first = .TRUE.
   double precision :: sum, sum_a, sum_b, sum_aa, sum_bb
   double precision, allocatable :: dist(:,:)
   integer, allocatable :: jbonds(:,:), njbonds(:)
   integer, external :: quoted
     allocate(pairs(2,numat*10), dist(3,numat*10), jbonds(15,numat), njbonds(numat))
!
! Work out connectivity twice - the sequence of atoms on one set might be
! different to that in the other set.
!
     coord(:,:numat) = geoa(:,:numat)
     call set_up_dentate()
     jbonds(:,:numat) = ibonds(:,:numat)
     njbonds(:numat) = nbonds(:numat)
     coord(:,:numat) = geo(:,:numat)
     call set_up_dentate()
!
!  Decide on bonds to an atom by whichever system has the more bonds to that atom
!
     do i = 1, numat
       if (njbonds(i) > nbonds(i)) then
         nbonds(i) = njbonds(i)
         ibonds(:nbonds(i),i) = jbonds(:nbonds(i),i)
       end if
       if (nbonds(i) > njbonds(i)) then
         njbonds(i) = nbonds(i)
        jbonds(:nbonds(i),i) = ibonds(:nbonds(i),i)
       end if
     end do
     npairs = 0
     do i = 1, numat
       do j = 1, nbonds(i)
         k = ibonds(j,i)
!
! Decide which atom(s) in the first system matche(s) atom "k" in the second system
! If the choice is unambiguous, use that atom, otherwise use all atoms that are bonded to atom "i"
! and then pick the best match.
!
         m = 0
         n = 0
         do j1 = 1, nbonds(i)
           l = jbonds(j1,i)
           if (txtatm(l)(12:) == txtatm(k)(12:)) then
             m = l
             n = n + 1
           end if
         end do
         if (n == 1) then
!
!  The label on one and only one atom matches atom "i"
!
           ijbonds(1) = m
           m = 1
         else
           m = nbonds(i)
           ijbonds(:m) = jbonds(:m,i)
         end if
         sum = 10.d0
         sum_aa = 0.d0
         sum_bb = 0.d0
         do j1 = 1, m
           l = ijbonds(j1)
           if (l > i .or. k > i .or. k == 0 .or. l == 0) cycle
           if (nat(k) /= nat(l)) cycle
           if (maxtxt == 26) then
             if (txtatm(k)(13:16) /= txtatm(l)(13:16)) cycle
           end if
           sum_a = sqrt((geo(1,i) - geo(1,k))**2 + (geo(2,i) - geo(2,k))**2 + (geo(3,i) - geo(3,k))**2)
           sum_b = sqrt((geoa(1,i) - geoa(1,l))**2 + (geoa(2,i) - geoa(2,l))**2 + (geoa(3,i) - geoa(3,l))**2)
           if (sum > abs(sum_a - sum_b)) then
             sum = abs(sum_a - sum_b)
             sum_aa = sum_a
             sum_bb = sum_b
           end if
         end do
         if (sum < 9.d0) then
           npairs = npairs + 1
           pairs(1, npairs) = i
           pairs(2, npairs) = k
           dist(1,npairs) = sum
           dist(2,npairs) = sum_aa
           dist(3,npairs) = sum_bb
         end if
       end do
     end do
    !
    ! Write out differences,largest first
    !
    do i = 1, min(25, npairs)
      sum = 0.d0
      k = 0
      do j = 1, npairs
        if (sum < dist(1,j)) then
          sum = dist(1,j)
          k = j
        end if
      end do
      if (sum < 0.01d0) exit
      j = pairs(1,k)
      l = pairs(2,k)
      if (nat(l) == nat(j) .and. j > l .or. nat(l) > nat(j)) then
        j = pairs(2,k)
        l = pairs(1,k)
      end if
      if (first) then
        if (trim(job_fn) == trim(geo_dat_name)) then
          line = "dataset"
        else
          line = "GEO_DAT"
        end if
        write(iw,'(/28x,a)')" Differences between bond-lengths for the two geometries"
        write(iw,'(/,a)')"           Diff.                 Atom 1           "//&
          &"            Atom 2                  Bond length      Bond length"
        write(iw,'(85x,a)')" in "//trim(line)//"      in GEO_REF"
        first = .false.
      end if
      write(iw,'(i4, SP, f12.3, S, 4x, a, 3x, a, f11.3, 5x, f12.3)') &
        i, dist(2,k) - dist(3,k), elemnt(nat(j))//"("//trim(txtatm(j))//")", &
        elemnt(nat(l))//"("//trim(txtatm(l))//")", dist(2,k), dist(3,k)
      dist(1,k) = -10.d0
    end do
    if ( .not. first) then
      i = quoted(" GEO_DAT")
      if (i > 0) then
        do k = 1, 6
          line = " "//trim(refkey(k))
          call upcase(line, len_trim(line))
          i = index(line," GEO_DAT")
          if (i /= 0) exit
        end do
        j = index(refkey(k)(i + 10:),'"') + index(refkey(k)(i + 10:),"'") + i + 8
        line = refkey(k)(i + 9:j)
        write(iw,'(/3x,a)')" GEO_DAT:  '"//trim(line)//"'"
      else
        write(iw,'(/3x,a)')" Data set: '"//trim(job_fn)//"'"
      end if
      do k = 1, 6
        line = " "//trim(refkey(k))
        call upcase(line, len_trim(line))
        i = index(line," GEO_REF")
        if (i /= 0) exit
      end do
      j = index(refkey(k)(i + 10:),'"') + index(refkey(k)(i + 10:),"'") + i + 8
      line = refkey(k)(i + 9:j)
      write(iw,'(3x,a)')" GEO_REF:  '"//trim(line)//"'"
    end if
  end subroutine wrt_diffs
  subroutine analyze_h_bonds
!
!  Using two geometries, compare the hydrogen bonds to find:
!  (A) How many hydrogen bonds are the same.
!  (B) Which hydrogen bonds are in the first geometry and not in the second.
!  (C) Which hydrogen bonds are in the second geometry and not in the first.
    use molkst_C, only : numat,  P_Hbonds, numcal, geo_ref_name, geo_dat_name, &
    E_disp, E_hb, E_hh, line
!
    use common_arrays_C, only : coord, geoa, H_energy, H_txt
!
    use chanel_C, only : iw, job_fn
    implicit none
    double precision :: sum, E_1_disp, E_2_disp, E_1_hb, E_2_hb, E_1_hh, E_2_hh, sum_a
    integer :: set_1_P_Hbonds, set_2_P_Hbonds, i, j, k, l, n_common
    logical :: l_prt
    double precision, allocatable :: set_1_H_energy(:), set_2_H_energy(:)
    logical, allocatable :: l_used_1(:), l_used_2(:)
    character, allocatable :: set_1_H_txt(:)*62, set_2_H_txt(:)*62, set_1_H_dist(:)*5, set_2_H_dist(:)*5
      allocate(set_1_H_txt(numat), set_2_H_txt(numat), set_1_H_energy(numat), set_2_H_energy(numat), &
        l_used_1(numat), l_used_2(numat), set_1_H_dist(numat), set_2_H_dist(numat))
      call l_control("PRT DISP(1.0) SILENT", len("PRT DISP(1.0) SILENT"), 1)
      call get_H_bonds
      E_1_disp = E_disp
      E_1_hb = E_hb
      E_1_hh = E_hh
      set_1_P_Hbonds = P_Hbonds
      do i = 1, P_Hbonds
        set_1_H_txt(i) = H_txt(i)(:31)//trim(H_txt(i)(37:68))
        set_1_H_dist(i) = H_txt(i)(32:36)
      end do
      if (P_Hbonds > 0) set_1_H_energy(:P_Hbonds) = H_energy(:P_Hbonds)
      coord(:,:numat) = geoa(:,:numat)
      numcal = numcal + 1
      call get_H_bonds
      E_2_disp = E_disp
      E_2_hb = E_hb
      E_2_hh = E_hh
      set_2_P_Hbonds = P_Hbonds
      do i = 1, P_Hbonds
        set_2_H_txt(i) = H_txt(i)(:31)//trim(H_txt(i)(37:68))
        set_2_H_dist(i) = H_txt(i)(32:36)
      end do
      if (P_Hbonds > 0) set_2_H_energy(:P_Hbonds) = H_energy(:P_Hbonds)
      i = max(len_trim(geo_ref_name), len_trim(geo_dat_name))
      write(iw,'(//20x,a)')"Analysis of Non-Covalent Interactions"
!
!  How many hydrogen bonds are the same?
!
      l_used_1 = .false.
      l_used_2 = .false.
      n_common = 0
      do i = 1, set_1_P_Hbonds
        do j = 1, set_2_P_Hbonds
          if (set_1_H_txt(i) == set_2_H_txt(j)) then
            l_used_1(i) = .true.
            l_used_2(j) = .true.
            n_common = n_common + 1
          end if
        end do
      end do
      if (trim(job_fn) == trim(geo_dat_name)) then
        line = "dataset"
      else
        line = "GEO_DAT"
      end if
      write(iw,'(/10x,a, f9.2, a)')"Total non-covalent energy of "//trim(line)//" system: ", &
        E_1_disp + E_1_hb + E_1_hh, " Kcal/mol"
      write(iw,'(10x,a, f9.2, a)')"Total non-covalent energy of GEO_REF system: ", &
        E_2_disp + E_2_hb + E_2_hh, " Kcal/mol"
      sum = E_1_disp + E_1_hb + E_1_hh - E_2_disp - E_2_hb - E_2_hh
      write(iw,'(43x, a, f9.2, a)')"Difference: ", sum, " Kcal/mol"
      if (n_common > 0) write(iw,'(/10x,a, i5)')"Number of hydrogen bonds common to both systems:", n_common
!
!  Which hydrogen bonds are in the data set but not in the reference?
!
      j = 0
      l_prt = .true.
      do i = 1, set_1_P_Hbonds
        if (.not. l_used_1(i)) then
          j = j + 1
          set_1_H_txt(j) = set_1_H_txt(i)
          set_1_H_energy(j) = set_1_H_energy(i)
          set_1_H_dist(j) = set_1_H_dist(i)
        end if
      end do
      sum_a = 0.d0
      do l = 1, j
        k = 0
        sum = 10.d0
        do i = 1, j
          if (set_1_H_energy(i) < sum) then
            k = i
            sum = set_1_H_energy(i)
          end if
        end do
        if (sum > -0.8d0) exit
        if (set_1_H_dist(k) > "3.000") cycle
        if (l_prt) then
          if (trim(job_fn) == trim(geo_dat_name)) then
            line = "dataset"
          else
            line = "GEO_DAT"
          end if
          write(iw,'(//28x,a)')"Hydrogen bonds in "//trim(line)//" but not in GEO_REF"
          write(iw,'(/15x,a,21x,a,12x,a,5x,a,/)')"Donor atom","Hydrogen atom","H-bond length(A)","Energy    Sum"
          l_prt = .false.
        end if
        sum_a = sum_a + set_1_H_energy(k)
        write(iw,'(i4,3x,a, 7x, a, f17.3, f9.3)')l, set_1_H_txt(k), set_1_H_dist(k), set_1_H_energy(k), sum_a
        set_1_H_energy(k) = 5.d0
      end do
      j = 0
      l_prt = .true.
      do i = 1, set_2_P_Hbonds
        if (.not. l_used_2(i)) then
          j = j + 1
          set_2_H_txt(j) = set_2_H_txt(i)
          set_2_H_energy(j) = set_2_H_energy(i)
          set_2_H_dist(j) = set_2_H_dist(i)
        end if
      end do
      sum_a = 0.d0
      do l = 1, j
        k = 0
        sum = 10.d0
        do i = 1, j
          if (set_2_H_energy(i) < sum) then
            k = i
            sum = set_2_H_energy(i)
          end if
        end do
        if (sum > -0.8d0) exit
        if (set_2_H_dist(k) > "3.000") cycle
        if (l_prt) then
          if (trim(job_fn) == trim(geo_dat_name)) then
            line = "dataset"
          else
            line = "GEO_DAT"
          end if
          write(iw,'(//28x,a)')"Hydrogen bonds in GEO_REF but not in "//trim(line)
          write(iw,'(/15x,a,21x,a,12x,a,5x,a,/)')"Donor atom","Hydrogen atom","H-bond length(A)","Energy    Sum"
          l_prt = .false.
        end if
        sum_a = sum_a + set_2_H_energy(k)
        write(iw,'(i4, 3x, a, 7x, a, f17.3, f9.3)')l, set_2_H_txt(k), set_2_H_dist(k), set_2_H_energy(k), sum_a
        set_2_H_energy(k) = 5.d0
      end do
      return
  end subroutine analyze_h_bonds
  subroutine get_H_bonds
    logical :: L_grad = .false.
    double precision :: correction
    call post_scf_corrections(correction, l_grad)
    return
  end subroutine get_H_bonds
