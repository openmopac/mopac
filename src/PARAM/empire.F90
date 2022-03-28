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

subroutine empire (formla, nat, els)
    use molkst_C, only : numat, keywrd, natoms
    use common_arrays_C, only : labels
    implicit none
    character (len=40), intent (out) :: formla
    integer, dimension (numat), intent (in) :: nat
    character (len=8), intent (out) :: els
    logical :: first = .true.
    integer :: i, j
    integer :: maxele
    integer :: l, nc, ne, mers(3), k, nele(109), m, nelements
    character (len=2), dimension (107) :: elemnt, lowel
    character (len=3), dimension (14) :: presnt
    character :: num
    integer, dimension (10) :: ntype
    integer, dimension (107) :: nmap, nos, nosort, numz
    double precision, external :: reada
    save :: first, nmap,nos, nosort, numz, elemnt, lowel, presnt, maxele
!
!.. Data Declarations ..
!
!  Sequence taken from the CRC Handbook page B-30
!
    data lowel /       &
 & "Cs","Rb","K ","Na","Li", & !  Group I
 & "Ba","Sr","Ca","Mg","Be", & !  Group IIA
 & "La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Tm","Yb","Lu", &
 &      "Y ","Sc",           & !  Group IIIB
 & "Hf","Zr","Ti",           & !  Group IVB
 & "Ta","Nb","V ",           & !  Group VB
 & "W ","Mo","Cr",           & !  Group VIB
 & "Re","Tc","Mn",           & !  Group VIIB
 & "Pt","Ir","Os",           & !  This is the sequence in CRC B-30
 & "Pd","Rh","Ru","Ni","Co","Fe","Au",      &
 & "Ag","Cu","Hg","Cd","Zn","Tl","In","Ga","Al","Pb","Sn","Ge","Bi", &
 & "B ","Si","C ","Sb","As","P ","N ","H ","Te","Se","S ","I ","Br", &
 & "Cl","O ","F ","He","Ne","Ar","Kr","Xe",25*"XX"/

    data elemnt / "H ", "He", "Li", "Be", "B ", "C ", "N ", "O ", "F ", "Ne", &
   & "Na", "Mg", "Al", "Si", "P ", "S ", "Cl", "Ar", "K ", "Ca", "Sc", "Ti", "&
   &V ", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br&
   &", "Kr", "Rb", "Sr", "Y ", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", &
   & "Cd", "In", "Sn", "Sb", "Te", "I ", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "&
   &Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf&
   &", "Ta", "W ", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", &
   & "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U ", "Np", "Pu", "Am", "Cm", "&
   &Bk", "Cf", "XX", "Fm", "Md", "Cb", "++", " +", "--", " -", "Tv" /
  !
  ! ... Executable Statements ...
  !
    if (first) then
      first = .false.
    !
    !  CONSTRUCT THE MAP TO DETERMINE THE ORDER IN WHICH THE ELEMENTS WILL
    !  BE PRINTED IN THE EMPIRIC FORMULA
    !
      do i = 1, 107
        if (lowel(i) =="XX") exit
        do j = 1, 107
          if (elemnt(j) == lowel(i)) go to 1000
        end do
        go to 1100
1000    nmap(i) = j
      end do
      maxele = i - 1
      go to 1200
1100  write (6,*) "   FAULT:", lowel (i), i
      stop
    end if
1200 els = " "
    do i = 1, maxele
      if (lowel(i) (2:2) == " ") then
        numz(i) = 1
      else
        numz(i) = 2
      end if
    end do
  !
  !  INITIALIZE ELEMENT ARRAYS
  !
    do i = 1, 107
      nos(i) = 0
    end do
  !
  !  READ FORMULA AND WORK OUT ELEMENTAL COMPOSITION
  !
    do i = 1, numat
      j = nat(i)
      nos(j) = nos(j) + 1
    end do
  !
  !   SEQUENCE ELEMENTS TO LOOK NICE
  !
    do i = 1, maxele
      nosort(i) = nos(nmap(i))
    end do
  !
  !  REDUCE IT TO AN EMPIRICAL FORMULA
  !
    nelements = 0
    do i = 1, maxele
      if (nosort(i) /= 0) then
        nelements = nelements + 1
        ntype(nelements) = i
        presnt(nelements) = lowel(i)
      end if
    end do
    j = Index (keywrd, " MERS")
    if (j /= 0) then
      mers = 0
      k = 0
      i = Index (keywrd(j + 1:), " ") + j
      do l = 1, 3
        j = j + k
        if (l > 1 .and. k == 0) exit
        mers(l) = Nint (reada (keywrd(j:), 1))
        k = Index (keywrd(j:i), ",")
      end do
      i = Index(keywrd, " Z=")
      if (i /= 0) then
        k = Nint(reada(keywrd,i))
        k = mers(1)*mers(2)*mers(3)*k
        if (Index (keywrd, " BCC") /= 0) k = k/2
      else
        nele = 0
        do i = 1, natoms - 3
          nele(labels(i)) = nele(labels(i)) + 1
        end do
        j = 0
        do i = 1, 98
          if (nele(i) > 0) then
            j = j + 1
            nele(j) = nele(i)
          end if
        end do
        k = 1000
        do i = 1, j
          if(nele(i) < k) k = nele(i)
        end do
  !
  !  k is the smallest number of atoms of any element in the formula
  !
        do i = 1, 100
          m = 0
          do l = 1, j
            if (Abs((i*nele(l))/k - (i*1.d0*nele(l))/k) > 1.d-5) m = 1
          end do
          if (m == 0) exit
        end do
        k = k/i
      end if
      nosort = nosort/max(k,1)
    end if
!
!   CONVERT IT TO A CHARACTER FORMULA
!
    formla = " "
    ne = 2
    nc = 2
    do i = 1, nelements
      formla(nc:) = presnt(i) (:numz(ntype(i)))
      nc = nc + numz(ntype(i))
      els(ne:) = presnt(i) (:numz(ntype(i)))
      ne = ne + numz(ntype(i))
      if (nosort(ntype(i)) /= 1) then
        if (nosort(ntype(i)) == 0) then
   !       write(*,'(a)')" Z definitely wrong in "//trim(title)
    !      write(31,'(a)')" Z definitely wrong in "//trim(title)
        else
          j = int(log10(nosort(ntype(i))*1.0)) + 1
          num = char(ichar("0") + j)
          write(formla(nc:),"(i"//num//")")nosort(ntype(i))
          nc = len_trim(formla) + 1
        end if
      end if
    end do
end subroutine empire
