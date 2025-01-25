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
