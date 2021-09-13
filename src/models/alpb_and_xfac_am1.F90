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

  subroutine alpb_and_xfac_am1
    use parameters_C, only : xfac, alpb
    xfac = 0.d0
    alpb = 0.d0
    alpb( 3, 1) =   2.975116d0 !     Lithium -     Hydrogen
    xfac( 3, 1) =  10.000006d0 !     Lithium -     Hydrogen
    alpb( 3, 3) =   2.074231d0 !     Lithium -      Lithium
    xfac( 3, 3) =  10.000005d0 !     Lithium -      Lithium
    !
    alpb( 4, 1) =   1.521999d0 !   Beryllium -     Hydrogen
    xfac( 4, 1) =   0.535972d0 !   Beryllium -     Hydrogen
    alpb( 4, 4) =   0.999154d0 !   Beryllium -    Beryllium
    xfac( 4, 4) =   0.441202d0 !   Beryllium -    Beryllium
    !
    alpb( 6, 3) =   1.709195d0 !      Carbon -      Lithium
    xfac( 6, 3) =   1.679411d0 !      Carbon -      Lithium
    alpb( 6, 4) =   0.999473d0 !      Carbon -    Beryllium
    xfac( 6, 4) =   0.132557d0 !      Carbon -    Beryllium
    !
    alpb( 7, 3) =   1.924790d0 !    Nitrogen -      Lithium
    xfac( 7, 3) =   1.762718d0 !    Nitrogen -      Lithium
    alpb( 7, 4) =   0.999523d0 !    Nitrogen -    Beryllium
    xfac( 7, 4) =   0.099561d0 !    Nitrogen -    Beryllium
    !
    alpb( 8, 3) =   1.873415d0 !      Oxygen -      Lithium
    xfac( 8, 3) =   1.294666d0 !      Oxygen -      Lithium
    alpb( 8, 4) =   2.816061d0 !      Oxygen -    Beryllium
    xfac( 8, 4) =   1.104010d0 !      Oxygen -    Beryllium
    !
    alpb( 9, 3) =   1.890271d0 !    Fluorine -      Lithium
    xfac( 9, 3) =   1.542043d0 !    Fluorine -      Lithium
    alpb( 9, 4) =   1.866431d0 !    Fluorine -    Beryllium
    xfac( 9, 4) =   0.411014d0 !    Fluorine -    Beryllium
    !
    alpb(11, 1) =   2.430601d0 !      Sodium -     Hydrogen
    xfac(11, 1) =   9.994724d0 !      Sodium -     Hydrogen
    alpb(11, 6) =   1.711276d0 !      Sodium -       Carbon
    xfac(11, 6) =   2.226935d0 !      Sodium -       Carbon
    alpb(11, 7) =   2.323882d0 !      Sodium -     Nitrogen
    xfac(11, 7) =  10.000002d0 !      Sodium -     Nitrogen
    alpb(11, 8) =   2.161260d0 !      Sodium -       Oxygen
    xfac(11, 8) =   4.057658d0 !      Sodium -       Oxygen
    alpb(11, 9) =   2.056819d0 !      Sodium -     Fluorine
    xfac(11, 9) =   3.274203d0 !      Sodium -     Fluorine
    alpb(11,11) =   1.612164d0 !      Sodium -       Sodium
    xfac(11,11) =   9.944446d0 !      Sodium -       Sodium
    !
    alpb(12, 1) =   2.174856d0 !   Magnesium -     Hydrogen
    xfac(12, 1) =   2.832813d0 !   Magnesium -     Hydrogen
    alpb(12, 6) =   2.216434d0 !   Magnesium -       Carbon
    xfac(12, 6) =   2.837201d0 !   Magnesium -       Carbon
    alpb(12, 7) =   2.091799d0 !   Magnesium -     Nitrogen
    xfac(12, 7) =   1.504138d0 !   Magnesium -     Nitrogen
    alpb(12, 8) =   1.666897d0 !   Magnesium -       Oxygen
    xfac(12, 8) =   0.472916d0 !   Magnesium -       Oxygen
    alpb(12, 9) =   2.142571d0 !   Magnesium -     Fluorine
    xfac(12, 9) =   1.227425d0 !   Magnesium -     Fluorine
    alpb(12,12) =   1.794083d0 !   Magnesium -    Magnesium
    xfac(12,12) =  10.000008d0 !   Magnesium -    Magnesium
    !
    alpb(14,12) =   1.900000d0 !     Silicon -    Magnesium
    xfac(14,12) =   2.000000d0 !     Silicon -    Magnesium
    !
    alpb(15,11) =   1.850000d0 !  Phosphorus -       Sodium
    xfac(15,11) =   2.000000d0 !  Phosphorus -       Sodium
    !
    alpb(16, 3) =   2.693882d0 !      Sulfur -      Lithium
    xfac(16, 3) =  10.000007d0 !      Sulfur -      Lithium
    alpb(16, 4) =   0.999448d0 !      Sulfur -    Beryllium
    xfac(16, 4) =   0.200086d0 !      Sulfur -    Beryllium
    alpb(16,11) =   1.903681d0 !      Sulfur -       Sodium
    xfac(16,11) =   1.053743d0 !      Sulfur -       Sodium
    alpb(16,12) =   1.572524d0 !      Sulfur -    Magnesium
    xfac(16,12) =   0.571504d0 !      Sulfur -    Magnesium
    !
    alpb(17, 3) =   1.779447d0 !    Chlorine -      Lithium
    xfac(17, 3) =   1.436235d0 !    Chlorine -      Lithium
    alpb(17, 4) =   1.382661d0 !    Chlorine -    Beryllium
    xfac(17, 4) =   0.211418d0 !    Chlorine -    Beryllium
    alpb(17,11) =   1.968902d0 !    Chlorine -       Sodium
    xfac(17,11) =   4.053650d0 !    Chlorine -       Sodium
    alpb(17,12) =   1.875713d0 !    Chlorine -    Magnesium
    xfac(17,12) =   1.195314d0 !    Chlorine -    Magnesium
    !
    alpb(19, 1) =   2.224557d0 !   Potassium -     Hydrogen
    xfac(19, 1) =   9.800241d0 !   Potassium -     Hydrogen
    alpb(19, 6) =   1.687938d0 !   Potassium -       Carbon
    xfac(19, 6) =   3.166630d0 !   Potassium -       Carbon
    alpb(19, 7) =   1.938615d0 !   Potassium -     Nitrogen
    xfac(19, 7) =   7.651496d0 !   Potassium -     Nitrogen
    alpb(19, 8) =   1.461751d0 !   Potassium -       Oxygen
    xfac(19, 8) =   1.046264d0 !   Potassium -       Oxygen
    alpb(19, 9) =   1.402823d0 !   Potassium -     Fluorine
    xfac(19, 9) =   0.811417d0 !   Potassium -     Fluorine
    alpb(19,16) =   3.439238d0 !   Potassium -       Sulfur
    xfac(19,16) =   9.817827d0 !   Potassium -       Sulfur
    alpb(19,17) =   2.072076d0 !   Potassium -     Chlorine
    xfac(19,17) =   8.149694d0 !   Potassium -     Chlorine
    alpb(19,19) =   1.267822d0 !   Potassium -    Potassium
    xfac(19,19) =   9.852489d0 !   Potassium -    Potassium
    !
    alpb(20, 1) =   1.593033d0 !     Calcium -     Hydrogen
    xfac(20, 1) =   3.654870d0 !     Calcium -     Hydrogen
    alpb(20, 6) =   1.005258d0 !     Calcium -       Carbon
    xfac(20, 6) =   0.567133d0 !     Calcium -       Carbon
    alpb(20, 7) =   1.003892d0 !     Calcium -     Nitrogen
    xfac(20, 7) =   0.453895d0 !     Calcium -     Nitrogen
    alpb(20, 8) =   2.574888d0 !     Calcium -       Oxygen
    xfac(20, 8) =   2.892181d0 !     Calcium -       Oxygen
    alpb(20, 9) =   2.268048d0 !     Calcium -     Fluorine
    xfac(20, 9) =   1.986328d0 !     Calcium -     Fluorine
    alpb(20,16) =   0.999897d0 !     Calcium -       Sulfur
    xfac(20,16) =   0.335089d0 !     Calcium -       Sulfur
    alpb(20,17) =   1.694685d0 !     Calcium -     Chlorine
    xfac(20,17) =   1.371206d0 !     Calcium -     Chlorine
    alpb(20,20) =   0.999786d0 !     Calcium -      Calcium
    xfac(20,20) =   4.524140d0 !     Calcium -      Calcium
    !
    alpb(35, 3) =   1.886037d0 !     Bromine -      Lithium
    xfac(35, 3) =   2.214447d0 !     Bromine -      Lithium
    alpb(35, 4) =   1.166450d0 !     Bromine -    Beryllium
    xfac(35, 4) =   0.193603d0 !     Bromine -    Beryllium
    alpb(35,11) =   2.189596d0 !     Bromine -       Sodium
    xfac(35,11) =   9.992652d0 !     Bromine -       Sodium
    alpb(35,12) =   1.939738d0 !     Bromine -    Magnesium
    xfac(35,12) =   1.898893d0 !     Bromine -    Magnesium
    alpb(35,19) =   1.997576d0 !     Bromine -    Potassium
    xfac(35,19) =   9.709128d0 !     Bromine -    Potassium
    alpb(35,20) =   1.635628d0 !     Bromine -      Calcium
    xfac(35,20) =   1.769827d0 !     Bromine -      Calcium
    !
    alpb(37, 1) =   1.999930d0 !    Rubidium -     Hydrogen
    xfac(37, 1) =   7.589965d0 !    Rubidium -     Hydrogen
    alpb(37, 5) =   1.999970d0 !    Rubidium -        Boron
    xfac(37, 5) =   6.012783d0 !    Rubidium -        Boron
    alpb(37, 8) =   1.999662d0 !    Rubidium -       Oxygen
    xfac(37, 8) =   2.914922d0 !    Rubidium -       Oxygen
    alpb(37, 9) =   2.914271d0 !    Rubidium -     Fluorine
    xfac(37, 9) =   8.743147d0 !    Rubidium -     Fluorine
    alpb(37,17) =   1.999948d0 !    Rubidium -     Chlorine
    xfac(37,17) =   3.708034d0 !    Rubidium -     Chlorine
    alpb(37,35) =   1.999755d0 !    Rubidium -      Bromine
    xfac(37,35) =   6.443546d0 !    Rubidium -      Bromine
    alpb(37,37) =   1.999854d0 !    Rubidium -     Rubidium
    xfac(37,37) =  10.000003d0 !    Rubidium -     Rubidium
    !
    alpb(38, 1) =   1.491059d0 !   Strontium -     Hydrogen
    xfac(38, 1) =   8.262735d0 !   Strontium -     Hydrogen
    alpb(38, 6) =   3.009340d0 !   Strontium -       Carbon
    xfac(38, 6) =   5.935637d0 !   Strontium -       Carbon
    alpb(38, 7) =   3.004535d0 !   Strontium -     Nitrogen
    xfac(38, 7) =   6.408169d0 !   Strontium -     Nitrogen
    alpb(38, 8) =   2.987339d0 !   Strontium -       Oxygen
    xfac(38, 8) =   5.727702d0 !   Strontium -       Oxygen
    alpb(38, 9) =   1.849351d0 !   Strontium -     Fluorine
    xfac(38, 9) =   1.086554d0 !   Strontium -     Fluorine
    alpb(38,16) =   2.051857d0 !   Strontium -       Sulfur
    xfac(38,16) =   4.276806d0 !   Strontium -       Sulfur
    alpb(38,17) =   2.072171d0 !   Strontium -     Chlorine
    xfac(38,17) =   4.093194d0 !   Strontium -     Chlorine
    alpb(38,35) =   2.168786d0 !   Strontium -      Bromine
    xfac(38,35) =   8.440324d0 !   Strontium -      Bromine
    alpb(38,38) =   2.982057d0 !   Strontium -    Strontium
    xfac(38,38) =   6.335178d0 !   Strontium -    Strontium
!
    alpb(42, 1) =   2.240000d0 !  Molybdenum -     Hydrogen
    xfac(42, 1) =   1.000000d0 !  Molybdenum -     Hydrogen
    alpb(42, 6) =   2.465000d0 !  Molybdenum -       Carbon
    xfac(42, 6) =   2.500000d0 !  Molybdenum -       Carbon
    alpb(42, 7) =   2.324000d0 !  Molybdenum -     Nitrogen
    xfac(42, 7) =   1.500000d0 !  Molybdenum -     Nitrogen
    alpb(42, 8) =   2.492000d0 !  Molybdenum -       Oxygen
    xfac(42, 8) =   1.500000d0 !  Molybdenum -       Oxygen
    alpb(42, 9) =   2.485000d0 !  Molybdenum -     Fluorine
    xfac(42, 9) =   1.500000d0 !  Molybdenum -     Fluorine
    alpb(42,15) =   3.000000d0 !  Molybdenum -   Phosphorus
    xfac(42,15) =   6.000000d0 !  Molybdenum -   Phosphorus
    alpb(42,16) =   2.290000d0 !  Molybdenum -       Sulfur
    xfac(42,16) =   1.500000d0 !  Molybdenum -       Sulfur
    alpb(42,17) =   2.500000d0 !  Molybdenum -     Chlorine
    xfac(42,17) =   2.500000d0 !  Molybdenum -     Chlorine
    alpb(42,24) =   1.955826d0 !  Molybdenum -     Chromium
    xfac(42,24) =   1.677542d0 !  Molybdenum -     Chromium
    alpb(42,35) =   2.130000d0 !  Molybdenum -      Bromine
    xfac(42,35) =   1.500000d0 !  Molybdenum -      Bromine
    alpb(42,42) =   2.550000d0 !  Molybdenum -   Molybdenum
    xfac(42,42) =   6.000000d0 !  Molybdenum -   Molybdenum

    !
    alpb(53, 3) =   2.648765d0 !      Iodine -      Lithium
    xfac(53, 3) =   9.992209d0 !      Iodine -      Lithium
    alpb(53, 4) =   0.999470d0 !      Iodine -    Beryllium
    xfac(53, 4) =   0.101571d0 !      Iodine -    Beryllium
    alpb(53,11) =   2.130872d0 !      Iodine -       Sodium
    xfac(53,11) =  10.000013d0 !      Iodine -       Sodium
    alpb(53,12) =   2.685633d0 !      Iodine -    Magnesium
    xfac(53,12) =   9.956335d0 !      Iodine -    Magnesium
    alpb(53,19) =   2.062884d0 !      Iodine -    Potassium
    xfac(53,19) =   9.900715d0 !      Iodine -    Potassium
    alpb(53,20) =   2.164489d0 !      Iodine -      Calcium
    xfac(53,20) =   7.940431d0 !      Iodine -      Calcium
    alpb(53,38) =   2.011009d0 !      Iodine -    Strontium
    xfac(53,38) =   6.634610d0 !      Iodine -    Strontium
    alpb(53,42) =   2.150000d0 !      Iodine -   Molybdenum
    xfac(53,42) =   1.500000d0 !      Iodine -   Molybdenum
    !
    !
    alpb(55, 1) =   0.997604d0 !      Cesium -     Hydrogen
    xfac(55, 1) =   1.985832d0 !      Cesium -     Hydrogen
    alpb(55, 6) =   0.999375d0 !      Cesium -       Carbon
    xfac(55, 6) =   0.099375d0 !      Cesium -       Carbon
    alpb(55, 7) =   0.998759d0 !      Cesium -     Nitrogen
    xfac(55, 7) =   0.205926d0 !      Cesium -     Nitrogen
    alpb(55, 8) =   2.249702d0 !      Cesium -       Oxygen
    xfac(55, 8) =   1.307467d0 !      Cesium -       Oxygen
    alpb(55, 9) =   0.999854d0 !      Cesium -     Fluorine
    xfac(55, 9) =   0.119889d0 !      Cesium -     Fluorine
    alpb(55,15) =   0.999228d0 !      Cesium -   Phosphorus
    xfac(55,15) =   0.099806d0 !      Cesium -   Phosphorus
    alpb(55,16) =   0.997626d0 !      Cesium -       Sulfur
    xfac(55,16) =   0.766755d0 !      Cesium -       Sulfur
    alpb(55,17) =   1.061905d0 !      Cesium -     Chlorine
    xfac(55,17) =   0.324568d0 !      Cesium -     Chlorine
    alpb(55,35) =   0.999274d0 !      Cesium -      Bromine
    xfac(55,35) =   0.416656d0 !      Cesium -      Bromine
    alpb(55,53) =   1.135145d0 !      Cesium -       Iodine
    xfac(55,53) =   0.642690d0 !      Cesium -       Iodine
    alpb(55,55) =   0.994764d0 !      Cesium -       Cesium
    xfac(55,55) =   3.263516d0 !      Cesium -       Cesium
    !
    alpb(56, 1) =   1.081595d0 !      Barium -     Hydrogen
    xfac(56, 1) =   0.108779d0 !      Barium -     Hydrogen
    alpb(56, 6) =   1.005494d0 !      Barium -       Carbon
    xfac(56, 6) =   0.109472d0 !      Barium -       Carbon
    alpb(56, 7) =   1.004182d0 !      Barium -     Nitrogen
    xfac(56, 7) =   0.109951d0 !      Barium -     Nitrogen
    alpb(56, 8) =   1.630584d0 !      Barium -       Oxygen
    xfac(56, 8) =   0.465945d0 !      Barium -       Oxygen
    alpb(56, 9) =   2.230982d0 !      Barium -     Fluorine
    xfac(56, 9) =   2.476341d0 !      Barium -     Fluorine
    alpb(56,16) =   1.000116d0 !      Barium -       Sulfur
    xfac(56,16) =   0.714303d0 !      Barium -       Sulfur
    alpb(56,17) =   1.341436d0 !      Barium -     Chlorine
    xfac(56,17) =   1.041642d0 !      Barium -     Chlorine
    alpb(56,35) =   1.427436d0 !      Barium -      Bromine
    xfac(56,35) =   2.384411d0 !      Barium -      Bromine
    alpb(56,53) =   1.303416d0 !      Barium -       Iodine
    xfac(56,53) =   2.097065d0 !      Barium -       Iodine
    alpb(56,56) =   1.110416d0 !      Barium -       Barium
    xfac(56,56) =   0.101890d0 !      Barium -       Barium
    !
  end subroutine alpb_and_xfac_am1
