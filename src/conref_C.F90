  module conref_C 
    double precision, dimension(2,10) :: fpcref 
!
!   Fundamental physical constants
!
!   fpcref(1,*) = 1998 Codata  fpcref(2,8) = old constants
!
    data fpcref(1, 1) / 1.602176462d0 /         ! Elementary charge, in Coulombs*10**19
    data fpcref(2, 1) / 1.60217733d0 /
    data fpcref(1, 2) / 14.399643d0 /           ! a0*(one au in eV)
    data fpcref(2, 2) / 14.399d0 /
    data fpcref(1, 3) / 0.5291772083d0 /        ! a0
    data fpcref(2, 3) / 0.529167d0 /
    data fpcref(1, 4) / 27.2113834d0 /          ! one au in eV
    data fpcref(2, 4) / 27.21d0 /
    data fpcref(1, 5) / 1.9872065d0 /           ! "R", the gas constant, in calories
    data fpcref(2, 5) / 1.98726d0 /
    data fpcref(1, 6) / 6.6260755d-27 /         ! Planck's constant in erg-seconds
    data fpcref(2, 6) / 6.626d-27 /
    data fpcref(1, 7) / 1.3806503d-16 /         ! Boltzmann's constant
    data fpcref(2, 7) / 1.3807d-16 /
    data fpcref(1, 8) / 2.99792458d+10 /        ! Speed of light, in cm/sec
    data fpcref(2, 8) / 2.99776d+10 /
    data fpcref(1, 9) / 23.060529d0 /           ! 1eV in kcal/mol = (ref(1)*ref(10)*10**(-3))/4.184
    data fpcref(2, 9) / 23.061d0 /
    data fpcref(1, 10) / 6.0221367d23 /         ! Avogadro's number (1/mole)
    data fpcref(2, 10) / 6.02205d23 /
  end module conref_C

