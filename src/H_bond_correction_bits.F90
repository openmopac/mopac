logical function connected(atom_i, atom_j, criterion)
!
!   "connected" is true if atom_i and atom_j are within sqrt(criterion) Angstroms
!   of each other, false otherwise.
!
!   if "connected" then:
!
!   If the system is infinite (polymer, layer or solid),
!   then cell_ijk will hold the unit cell translation indices that move atom_j near to
!   atom_i, zero therwise.
!
!   Vab is the postition of atom_i relative to atom_j
!
!   Rab is distance from atom_i to atom_j
!
!   if not "connected" then Rab, Vab, and cell_ijk are meaningless.
!
  use molkst_C, only : id, l11, l21, l31, Rab
  use common_arrays_C, only: tvec, coord, cell_ijk, Vab
  implicit none
  double precision ::  criterion
  integer :: atom_i, atom_j, ii, jj, kk
  double precision :: V_ab(3), R_ab
!
  if (id == 0) then
    Vab = coord(:,atom_i) - coord(:,atom_j)
    Rab = Vab(1)**2 + Vab(2)**2 + Vab(3)**2 
  else
    Rab = 1.d8
    do ii = -l11, l11 
      do jj = -l21, l21 
        do kk = -l31, l31 
          V_ab = coord(:,atom_i) - coord(:,atom_j) + tvec(:,1)*ii + tvec(:,2)*jj + tvec(:,3)*kk 
          R_ab = V_ab(1)**2 + V_ab(2)**2 + V_ab(3)**2 
          if (R_ab < Rab) then
            Rab = R_ab
            Vab = V_ab
            cell_ijk(1) = ii
            cell_ijk(2) = jj
            cell_ijk(3) = kk
          end if                  
        end do 
      end do 
    end do 
  end if 
  connected = (Rab < criterion) 
  if (connected) Rab = sqrt(Rab)
  return       
end function connected
subroutine find_XH_bonds(acc, nacc, h_b, nhb)
!
!  Locate all hydrogen atoms involved in O-H and N-H bonds
!  On exit:
!  nhb:   Number of hydrogen atoms involved in O-H or N-H bonds
!  h_b:   Atom numbers of the hydrogen atoms
!
  use molkst_C, only : numat, method_PM7, method_pm6_dh_plus
  use common_arrays_C, only: nat
  implicit none
  integer :: nacc, nhb

  integer :: acc(numat), h_b(numat)
  integer :: i, j, is
  double precision :: RAH
  logical :: used(numat)
  logical, external :: connected
    if (method_pm6_dh_plus) then
      RAH = 1.4d0
      is = 8
    else if (method_PM7) then
      RAH = 1.4d0
      is = 8  ! Change to 16 ASAP
    else
      RAH = 1.15d0
      is = 16
    end if
    used = .false.
    nacc = 0
    nhb  = 0
    do i = 1, numat
      if (nat(i) ==  7 .or. nat(i) ==  8 .or. nat(i) ==  is) then
        nacc = nacc + 1
        acc(nacc) = i
        do j = 1, numat
          if (nat(j) ==  1 .and. .not. used(j)) then
            if (connected(i, j, RAH**2)) then
              nhb = nhb + 1
              h_b(nhb) = j
              used(j) = .true.
            end if
          end if
        end do
      end if
    end do
end subroutine find_XH_bonds
subroutine find_H__Y_bonds(acc_a, nacc_a, acc_b, nacc_b, bonding_a_h, nb_a_h, hblist1, hblist2, hblist3, max_h_bonds, nrpairs)
!
!  Find all sets of three atoms that form a hydrogen bond
!
!   acc_a:       All O and N atoms 
!   bonding_a_h: Hydrogen atoms bonded to acc_a, not in any particular order
!
  use chanel_C, only : iw
  use funcon_C, only : pi
  use molkst_C, only : keywrd, method_pm7
  implicit none
  integer :: nacc_a, nacc_b, nb_a_h, nrpairs, max_h_bonds
  integer :: acc_a(nacc_a), acc_b(nacc_b), bonding_a_h(nb_a_h), hblist1(max_h_bonds), hblist2(max_h_bonds), hblist3(max_h_bonds)
  integer:: ii, i, jj, j, kk, k, i1
  double precision :: RAH, cutoff
  logical, external :: connected
  double precision, external :: angle
  if (index(keywrd, "PM6-DH+") /= 0 ) then
      RAH = 1.4d0
      cutoff = 10.d0
    else if (method_PM7) then
      RAH = 1.4d0 
      cutoff = 7.d0
    else
      RAH = 1.15d0
      cutoff = 7.d0
    end if
    do ii = 1, nacc_a
      i = acc_a(ii)   !  i = Acceptor atom bonded to H  
      do jj = 1, nb_a_h
        j = bonding_a_h(jj)   !  j = Hydrogen bonded to acceptor atom
        if (connected(i, j, RAH**2)) then
          do kk = 1, nacc_b
            k = acc_b(kk)  
            if (k /= i) then
              if (connected(k, j, cutoff**2)) then
                if (angle(k,j,i) > pi*0.5d0) then
!
!  Eliminate bonds of type O(n) - H - O(m) if O(m) - H - O(n) exists
!
                  do i1 = 1, nrpairs
                    if (hblist2(i1) /= j) cycle
                    if (hblist1(i1) /= k) cycle
                    if (hblist3(i1) /= i) cycle
                    exit
                  end do
                  if (i1 /= nrpairs + 1) cycle
                  do i1 = 1, nrpairs
                    if (hblist2(i1) /= j) cycle
                    if (hblist1(i1) /= i) cycle
                    if (hblist3(i1) /= k) cycle
                    exit
                  end do
                  if (i1 /= nrpairs + 1) cycle                  
                  nrpairs = nrpairs + 1 !  k = Distant acceptor atom
                  if (nrpairs > max_h_bonds) then
                    write(iw,'(a)')" The default array size for hydrogen bonds is too small"
                    write(iw,'(a,i6,a,i6)')" Array size:", max_h_bonds,", estimated size needed:", &
                      nint(float(max_h_bonds)*float(nacc_a)/float(ii))
                    if (index(keywrd, "PM6-DH+") /= 0 ) then
                      write(iw,'(a)')" If possible, use PM6-DH+=(text)"
                      call web_message(iw,"PM6_DH_plus.html")
                    else
                      write(iw,'(a)')" If possible, use PM6-DH2=(text)"
                      call web_message(iw,"PM6_DH2.html")
                    end if
                    call mopend("The default array size for hydrogen bonds is too small")
                    nrpairs = nrpairs - 1
                    return
                  end if
                  hblist3(nrpairs) = k ! O or N hydrogen-bonded to hydrogen j
                  hblist2(nrpairs) = j ! H singly bonded to i
                  hblist1(nrpairs) = i ! O or N of single bond to H
                end if
              end if
            end if
          end do                
        end if          
      end do
    end do
  end subroutine find_H__Y_bonds
 
  function truncation(R, limit, spread)
!
!  Truncation has the values:
!
!  limit, when R < limit - spread
!  greater than R when R is in the range: limit - spread < R < limit + spread
!  R, when R > limit + spread
!
!  At limit + spread, the slope of truncation is 1, i.e., equal to the slope for R > limit + spread
!  At limit - spread, the slope is zero, i.e., equal to the slope of the constant "limit"
!
  implicit none
  double precision, intent(in) :: R, limit, spread
  double precision :: truncation, a, b
    a = limit - spread
    b = limit + spread  
    if (R < b) then
      if (R < a) then
        truncation = limit
      else      
        truncation = limit + (limit - a)/(a - b)**2*(R - a)**2
      end if
    else
      truncation = R
    end if
    return
  end function truncation



  subroutine all_h_bonds(hblist1, hblist2, hblist3, max_h_bonds, nrpairs)
!
!  all_h_bonds detects all potential hydrogen bonds.  A hydrogen bond is a set of three atoms,
!  an oxygen or nitrogen, a hydrogen, and an oxygen or nitrogen atom.  The hydrogen atom must
! be within RAH Angstroms of either an oxygen or a nitrogen atom.
!
!  On exit:
!
!    nrpairs: number of pairs of atoms
!    hblist1: Atom number of oxygen or nitrogen attached to hydrogen
!    hblist2: Atom number of hydrogen atom involved in hydrogen bonding
!    hblist3: Atom number ofoxygen or nitrogen hydrogen bonded to the hydrogen atom
!    
!
  use common_arrays_C, only: bonding_a_h, bonding_b_h, acceptor_a, acceptor_b
  use molkst_C, only : numat, line, moperr
  use chanel_C, only : iw
  use to_screen_I
  implicit none
    integer, intent (in) :: max_h_bonds
    integer, intent (out) :: hblist1(max_h_bonds), hblist2(max_h_bonds), hblist3(max_h_bonds), nrpairs
!
!  Local variables
!
    integer :: i, nacceptor_a, nbonding_a_h
!
!  Work out list of of potential hydrogen bonds
!
!   First, find all N-H and O-H groups 
!
    if (allocated(acceptor_a))   deallocate(acceptor_a)
    if (allocated(acceptor_b))   deallocate(acceptor_b)
    if (allocated(bonding_a_h))  deallocate(bonding_a_h)
    if (allocated(bonding_b_h))  deallocate(bonding_b_h)
    allocate (acceptor_a(numat*2), acceptor_b(numat*2), bonding_a_h(numat*2), &
      bonding_b_h(numat*2), stat=i)
    if (i /= 0) then
      line = " Cannot allocate arrays for hydrogen bonds"
      write(iw,*)" Cannot allocate arrays for hydrogen bonds"
      call to_screen(trim(line))
      call mopend(trim(line))
      return
    end if
    call find_XH_bonds(acceptor_a, nacceptor_a, bonding_a_h, nbonding_a_h)
!
!  acceptor_a now holds the list of acceptor O and N atoms (for the first fragment)
!  acceptor_b now holds the list of acceptor O and N atoms (for the second fragment, if it exists)
!  bonding_a_h now holds the list of H atoms attached to an O or N (for the first fragment)
!  bonding_b_h now holds the list of H atoms attached to an O or N (for the second fragment, if it exists)
!  
!
!  Find if there is a hydrogen bond
!     
    nrpairs = 0
    call find_H__Y_bonds(acceptor_a, nacceptor_a, acceptor_a, nacceptor_a, &
      bonding_a_h, nbonding_a_h, hblist1, hblist2, hblist3, max_h_bonds, nrpairs)
    if (moperr) return
  end subroutine all_h_bonds

  double precision function distance(a, b)
    use common_arrays_C, only: coord, tvec
    use molkst_C, only : id, l1u, l2u, l3u
    implicit none
      integer :: a,  b
      integer :: ik, jk, kl
      double precision :: coord1(3)
      if (id == 0) then
        distance = sqrt((coord(1,a) - coord(1,b))**2 + (coord(2,a) - coord(2,b))**2 + (coord(3,a) - coord(3,b))**2)
      else
        distance = 1.d6
        do ik = -l1u, l1u
          do jk = -l2u, l2u
            do kl = -l3u, l3u
              coord1(:) = coord(:, a) + tvec(:, 1)*ik + tvec(:, 2)*jk + tvec(:, 3)*kl
              distance = Min (distance, (coord1(1)-coord(1, b))**2 &
                                      + (coord1(2)-coord(2, b))**2 &
                                      + (coord1(3)-coord(3, b))**2)
            end do
          end do
        end do
        distance = sqrt(distance)
      end if    
  end function distance
  double precision function angle(a, b, c) 
    use common_arrays_C, only: coord
      integer :: a,  b,  c
      call bangle(coord, a, b, c, angle)
  end
  double precision function torsion(i, j, k, l)
      use common_arrays_C, only: coord
      integer :: i,  j,  k,  l
      call dihed (coord, i, j, k, l, torsion)
  end

  double precision function bonding(x, y, covrad)
    use common_arrays_C, only: nat
    implicit none
    integer,  INTENT(IN) :: x, y 
    double precision :: covrad(94)
      bonding = covrad(nat(x)) + covrad(nat(y)) 
    return
  end function bonding
  
  subroutine prt_hbonds(D, H, A, energy)
  use common_arrays_C, only: nat, txtatm, H_txt, H_energy
  use molkst_C, only : numat, keywrd, numcal, P_hbonds, maxtxt
  use elemts_C, only: elemnt
  implicit none
  integer, intent (in) :: D, H, A
  double precision, intent (in) :: energy
!
!  Local
!
  integer :: icalcn = -1
  double precision :: cutoff, sum1
  double precision, external :: distance, reada
  logical :: prt_first
  save :: prt_first, cutoff, icalcn
    if (icalcn /= numcal) then
      icalcn = numcal
      if (allocated(H_txt)) deallocate (H_txt, H_energy)
      prt_first = .true.
      allocate(H_txt(numat), H_energy(numat))
      P_hbonds = 0
    end if
    if (prt_first) then
      prt_first = .false.
      cutoff = reada(keywrd, index(keywrd," DISP("))
      cutoff = -abs(cutoff)
      P_Hbonds = 0
    end if
    if (energy > -0.5d0) return
    sum1 = distance(D,H)
    P_Hbonds = min(numat, P_Hbonds + 1)
    if (txtatm(D) == " ") write(txtatm(D),'(a,i5,3x,a)')"Atom No.:", D, elemnt(nat(D))
    if (txtatm(H) == " ") write(txtatm(H),'(a,i5,3x,a)')"Atom No.:", H, elemnt(nat(H))
    if (txtatm(A) == " ") write(txtatm(A),'(a,i5,3x,a)')"Atom No.:", A, elemnt(nat(A))
    if (maxtxt == 0) then
      write(H_txt(P_Hbonds),'(2x,a,f15.3,6x,a,12x,a,10x, f7.2,a)')trim(txtatm(D)), sum1,  &
      trim(txtatm(H)), trim(txtatm(A)), energy, " Kcal/mol"
    else      
      write(H_txt(P_Hbonds),'(a,2x,f6.3,3x,a,3x,a,3x, f7.2,a)')""""//txtatm(D)(:maxtxt)//"""", sum1,  &
      """"//txtatm(H)(:maxtxt)//"""", """"//txtatm(A)(:maxtxt)//"""", energy, " Kcal/mol"
    end if
    H_energy(P_Hbonds) = energy
    return
  end subroutine prt_hbonds

