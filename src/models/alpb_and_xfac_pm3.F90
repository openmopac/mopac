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

  subroutine alpb_and_xfac_pm3
    use parameters_C, only : xfac, alpb
    alpb = 0.d0
    xfac = 0.d0
    !
    alpb(11, 1) =   1.800472d0 !      Sodium -     Hydrogen
    xfac(11, 1) =   3.171946d0 !      Sodium -     Hydrogen
    alpb(11, 6) =   1.321600d0 !      Sodium -       Carbon
    xfac(11, 6) =   0.876072d0 !      Sodium -       Carbon
    alpb(11, 7) =   0.999895d0 !      Sodium -     Nitrogen
    xfac(11, 7) =   0.295812d0 !      Sodium -     Nitrogen
    alpb(11, 8) =   2.116028d0 !      Sodium -       Oxygen
    xfac(11, 8) =   3.392634d0 !      Sodium -       Oxygen
    alpb(11, 9) =   1.860719d0 !      Sodium -     Fluorine
    xfac(11, 9) =   1.605474d0 !      Sodium -     Fluorine
    alpb(11,11) =   1.024150d0 !      Sodium -       Sodium
    xfac(11,11) =   1.126080d0 !      Sodium -       Sodium
    !
    alpb(16,11) =   0.999990d0 !      Sulfur -       Sodium
    xfac(16,11) =   0.241710d0 !      Sulfur -       Sodium
    !
    alpb(17,11) =   1.420756d0 !    Chlorine -       Sodium
    xfac(17,11) =   1.168589d0 !    Chlorine -       Sodium
    !
    alpb(19, 1) =   1.027701d0 !   Potassium -     Hydrogen
    xfac(19, 1) =   0.893642d0 !   Potassium -     Hydrogen
    alpb(19, 6) =   1.453975d0 !   Potassium -       Carbon
    xfac(19, 6) =   2.009100d0 !   Potassium -       Carbon
    alpb(19, 7) =   0.999984d0 !   Potassium -     Nitrogen
    xfac(19, 7) =   0.425244d0 !   Potassium -     Nitrogen
    alpb(19, 8) =   1.470651d0 !   Potassium -       Oxygen
    xfac(19, 8) =   0.999945d0 !   Potassium -       Oxygen
    alpb(19, 9) =   0.999898d0 !   Potassium -     Fluorine
    xfac(19, 9) =   0.342063d0 !   Potassium -     Fluorine
    alpb(19,16) =   1.000576d0 !   Potassium -       Sulfur
    xfac(19,16) =   0.584513d0 !   Potassium -       Sulfur
    alpb(19,17) =   0.999989d0 !   Potassium -     Chlorine
    xfac(19,17) =   0.584823d0 !   Potassium -     Chlorine
    alpb(19,19) =   1.532875d0 !   Potassium -    Potassium
    xfac(19,19) =   8.250024d0 !   Potassium -    Potassium
    !
    alpb(20, 1) =   1.427846d0 !     Calcium -     Hydrogen
    xfac(20, 1) =   1.805489d0 !     Calcium -     Hydrogen
    alpb(20, 6) =   0.999993d0 !     Calcium -       Carbon
    xfac(20, 6) =   0.320183d0 !     Calcium -       Carbon
    alpb(20, 7) =   0.999993d0 !     Calcium -     Nitrogen
    xfac(20, 7) =   0.275261d0 !     Calcium -     Nitrogen
    alpb(20, 8) =   1.934512d0 !     Calcium -       Oxygen
    xfac(20, 8) =   0.480203d0 !     Calcium -       Oxygen
    alpb(20, 9) =   1.964914d0 !     Calcium -     Fluorine
    xfac(20, 9) =   0.594612d0 !     Calcium -     Fluorine
    alpb(20,16) =   0.999960d0 !     Calcium -       Sulfur
    xfac(20,16) =   0.212682d0 !     Calcium -       Sulfur
    alpb(20,17) =   1.669363d0 !     Calcium -     Chlorine
    xfac(20,17) =   0.928683d0 !     Calcium -     Chlorine
    alpb(20,20) =   0.999174d0 !     Calcium -      Calcium
    xfac(20,20) =   3.160373d0 !     Calcium -      Calcium
    !
    !
    alpb(35,11) =   1.517965d0 !     Bromine -       Sodium
    xfac(35,11) =   1.741658d0 !     Bromine -       Sodium
    alpb(35,19) =   1.401604d0 !     Bromine -    Potassium
    xfac(35,19) =   2.440790d0 !     Bromine -    Potassium
    alpb(35,20) =   1.509690d0 !     Bromine -      Calcium
    xfac(35,20) =   0.874420d0 !     Bromine -      Calcium
    !
    alpb(37, 1) =   2.066936d0 !    Rubidium -     Hydrogen
    xfac(37, 1) =   9.999919d0 !    Rubidium -     Hydrogen
    alpb(37, 5) =   1.996668d0 !    Rubidium -        Boron
    xfac(37, 5) =  10.000119d0 !    Rubidium -        Boron
    alpb(37, 8) =   1.999948d0 !    Rubidium -       Oxygen
    xfac(37, 8) =   1.817233d0 !    Rubidium -       Oxygen
    alpb(37, 9) =   3.083855d0 !    Rubidium -     Fluorine
    xfac(37, 9) =  10.000006d0 !    Rubidium -     Fluorine
    alpb(37,17) =   2.371423d0 !    Rubidium -     Chlorine
    xfac(37,17) =   9.970607d0 !    Rubidium -     Chlorine
    alpb(37,35) =   2.071332d0 !    Rubidium -      Bromine
    xfac(37,35) =   7.407687d0 !    Rubidium -      Bromine
    alpb(37,37) =   0.539344d0 !    Rubidium -     Rubidium
    xfac(37,37) =   2.654922d0 !    Rubidium -     Rubidium
    !
    alpb(38, 1) =   1.385332d0 !   Strontium -     Hydrogen
    xfac(38, 1) =   7.195639d0 !   Strontium -     Hydrogen
    alpb(38, 6) =   1.392807d0 !   Strontium -       Carbon
    xfac(38, 6) =   5.455084d0 !   Strontium -       Carbon
    alpb(38, 7) =   1.392514d0 !   Strontium -     Nitrogen
    xfac(38, 7) =   5.293460d0 !   Strontium -     Nitrogen
    alpb(38, 8) =   2.471371d0 !   Strontium -       Oxygen
    xfac(38, 8) =   2.147711d0 !   Strontium -       Oxygen
    alpb(38, 9) =   2.525416d0 !   Strontium -     Fluorine
    xfac(38, 9) =   5.359983d0 !   Strontium -     Fluorine
    alpb(38,16) =   1.045944d0 !   Strontium -       Sulfur
    xfac(38,16) =   0.594774d0 !   Strontium -       Sulfur
    alpb(38,17) =   1.473507d0 !   Strontium -     Chlorine
    xfac(38,17) =   1.176805d0 !   Strontium -     Chlorine
    alpb(38,35) =   1.338861d0 !   Strontium -      Bromine
    xfac(38,35) =   1.082828d0 !   Strontium -      Bromine
    alpb(38,38) =   1.393785d0 !   Strontium -    Strontium
    xfac(38,38) =   5.138251d0 !   Strontium -    Strontium
    !
    alpb(53,11) =   1.148784d0 !      Iodine -       Sodium
    xfac(53,11) =   0.747959d0 !      Iodine -       Sodium
    alpb(53,19) =   1.402113d0 !      Iodine -    Potassium
    xfac(53,19) =   2.568324d0 !      Iodine -    Potassium
    alpb(53,20) =   1.513600d0 !      Iodine -      Calcium
    xfac(53,20) =   0.856865d0 !      Iodine -      Calcium
    alpb(53,38) =   1.730393d0 !      Iodine -    Strontium
    xfac(53,38) =   3.855491d0 !      Iodine -    Strontium
    !
    alpb(55, 1) =   1.448617d0 !      Cesium -     Hydrogen
    xfac(55, 1) =   6.721688d0 !      Cesium -     Hydrogen
    alpb(55, 6) =   0.999871d0 !      Cesium -       Carbon
    xfac(55, 6) =   5.579858d0 !      Cesium -       Carbon
    alpb(55, 7) =   0.999559d0 !      Cesium -     Nitrogen
    xfac(55, 7) =   0.242498d0 !      Cesium -     Nitrogen
    alpb(55, 8) =   3.222538d0 !      Cesium -       Oxygen
    xfac(55, 8) =   7.547802d0 !      Cesium -       Oxygen
    alpb(55, 9) =   2.748652d0 !      Cesium -     Fluorine
    xfac(55, 9) =   3.360596d0 !      Cesium -     Fluorine
    alpb(55,15) =   0.999863d0 !      Cesium -   Phosphorus
    xfac(55,15) =   0.099894d0 !      Cesium -   Phosphorus
    alpb(55,16) =   0.999607d0 !      Cesium -       Sulfur
    xfac(55,16) =   0.869429d0 !      Cesium -       Sulfur
    alpb(55,17) =   0.999578d0 !      Cesium -     Chlorine
    xfac(55,17) =   0.288955d0 !      Cesium -     Chlorine
    alpb(55,35) =   0.999651d0 !      Cesium -      Bromine
    xfac(55,35) =   0.350614d0 !      Cesium -      Bromine
    alpb(55,53) =   0.998419d0 !      Cesium -       Iodine
    xfac(55,53) =   0.377054d0 !      Cesium -       Iodine
    alpb(55,55) =   0.997434d0 !      Cesium -       Cesium
    xfac(55,55) =   8.549732d0 !      Cesium -       Cesium
    !
    alpb(56, 1) =   0.999997d0 !      Barium -     Hydrogen
    xfac(56, 1) =   0.099997d0 !      Barium -     Hydrogen
    alpb(56, 6) =   0.999997d0 !      Barium -       Carbon
    xfac(56, 6) =   0.099997d0 !      Barium -       Carbon
    alpb(56, 7) =   0.999997d0 !      Barium -     Nitrogen
    xfac(56, 7) =   0.099997d0 !      Barium -     Nitrogen
    alpb(56, 8) =   1.249857d0 !      Barium -       Oxygen
    xfac(56, 8) =   0.352542d0 !      Barium -       Oxygen
    alpb(56, 9) =   1.886689d0 !      Barium -     Fluorine
    xfac(56, 9) =   1.528347d0 !      Barium -     Fluorine
    alpb(56,16) =   0.999862d0 !      Barium -       Sulfur
    xfac(56,16) =   0.904098d0 !      Barium -       Sulfur
    alpb(56,17) =   1.490059d0 !      Barium -     Chlorine
    xfac(56,17) =   1.846584d0 !      Barium -     Chlorine
    alpb(56,35) =   1.628325d0 !      Barium -      Bromine
    xfac(56,35) =   3.652542d0 !      Barium -      Bromine
    alpb(56,53) =   1.461169d0 !      Barium -       Iodine
    xfac(56,53) =   2.514512d0 !      Barium -       Iodine
    alpb(56,56) =   0.999997d0 !      Barium -       Barium
    xfac(56,56) =   0.099997d0 !      Barium -       Barium
    !
  end subroutine alpb_and_xfac_pm3
