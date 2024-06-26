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

module common_elemts
    implicit none
    character (len=2), dimension (107) :: elemnt
        data elemnt / " H", "He", "Li", "Be", " B", " C", " N", " O", " F", "Ne", &
         &"Na", "Mg", "Al", "Si", " P", " S", "Cl", "Ar", " K", "Ca", "Sc", &
         &"Ti", " V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", &
         &"As", "Se", "Br", "Kr", "Rb", "Sr", " Y", "Zr", "Nb", "Mo", "Tc", &
         &"Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", " I", "Xe", &
         &"Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", &
         &"Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", " W", "Re", "Os", &
         &"Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", &
         &"Ra", "Ac", "Th", "Pa", " U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", &
         &"XX", "Fm", "Md", "Cb", "++", " +", "--", " -", "Tv" /
end module common_elemts
