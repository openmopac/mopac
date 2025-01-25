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

      module elemts_C
      character, dimension(107) :: elemnt*2, cap_elemnt*2
      character (len=12), dimension(107) :: atom_names
      data elemnt/ ' H', 'He', 'Li', 'Be', ' B', ' C', ' N', ' O', ' F', 'Ne', &
        'Na', 'Mg', 'Al', 'Si', ' P', ' S', 'Cl', 'Ar', ' K', 'Ca', 'Sc', 'Ti'&
        , ' V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', &
        'Se', 'Br', 'Kr', 'Rb', 'Sr', ' Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh'&
        , 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', ' I', 'Xe', 'Cs', 'Ba', &
        'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er'&
        , 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', ' W', 'Re', 'Os', 'Ir', 'Pt', 'Au', &
        'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa'&
        , ' U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Mi', 'XX', '+3', '-3', 'Cb', &
        '++', ' +', '--', ' -', 'Tv'/
      data cap_elemnt / "H ", "HE", "LI", "BE", "B ", "C ", "N ", "O ", "F ", "NE", &
        & "NA", "MG", "AL", "SI", "P ", "S ", "CL", "AR", "K ", "CA", "SC", "TI", &
        & "V ", "CR", "MN", "FE", "CO", "NI", "CU", "ZN", "GA", "GE", "AS", "SE", "BR", &
        & "KR", "RB", "SR", "Y ", "ZR", "NB", "MO", "TC", "RU", "RH", "PD", "AG", &
        & "CD", "IN", "SN", "SB", "TE", "I ", "XE", "CS", "BA", "LA", "CE", "PR", &
        & "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU", "HF", &
        & "TA", " W", "RE", "OS", "IR", "PT", "AU", "HG", "TL", "PB", "BI", "PO", &
        & "AT", "RN", "FR", "RA", "AC", "TH", "PA", " U", "NP", "PU", "AM", "CM", &
        & "BK", "CF", "XX", "+3", "-3", "CB", '++', ' +', '--', ' -', 'TV'/
               data atom_names / &
   & "    Hydrogen", "      Helium", "     Lithium", "   Beryllium", &
   & "       Boron", "      Carbon", "    Nitrogen", "      Oxygen", &
   & "    Fluorine", "        Neon", "      Sodium", "   Magnesium", &
   & "    Aluminum", "     Silicon", "  Phosphorus", "      Sulfur", &
   & "    Chlorine", "       Argon", "   Potassium", "     Calcium", &
   & "    Scandium", "    Titanium", "    Vanadium", "    Chromium", &
   & "   Manganese", "        Iron", "      Cobalt", "      Nickel", &
   & "      Copper", "        Zinc", "     Gallium", "   Germanium", &
   & "     Arsenic", "    Selenium", "     Bromine", "     Krypton", &
   & "    Rubidium", "   Strontium", "     Yttrium", "   Zirconium", &
   & "     Niobium", "  Molybdenum", "  Technetium", "   Ruthenium", &
   & "     Rhodium", "   Palladium", "      Silver", "     Cadmium", &
   & "      Indium", "         Tin", "    Antimony", "   Tellurium", &
   & "      Iodine", "       Xenon", "      Cesium", "      Barium", &
   & "   Lanthanum", "      Cerium", "Praseodymium", "   Neodymium", &
   & "  Promethium", "    Samarium", "    Europium", "  Gadolinium", &
   & "     Terbium", "  Dysprosium", "     Holmium", "      Erbium", &
   & "     Thulium", "   Ytterbium", "    Lutetium", &
   & "     Hafnium", "    Tantalum", "    Tungsten", &
   & "     Rhenium", "      Osmium", "     Iridium", "    Platinum", &
   & "        Gold", "     Mercury", "    Thallium", "        Lead", &
   & "     Bismuth", "    Polonium", "    Astatine", "       Radon", &
   & "    Francium", "      Radium", "    Actinium", "     Thorium", &
   & "Protactinium", "     Uranium", "   Neptunium", "   Plutonium", &
   & "   Americium", "      Curium", "   Berkelium", "     Mithril", &
   & "  Dummy atom", &
   & "  3+ Sparkle", "  3- Sparkle", " Capped bond", "  ++ Sparkle", &
   & "   + Sparkle", "  -- Sparkle", "   - Sparkle", "     Tv     "/
      end module elemts_C
