  subroutine alpb_and_xfac_mndod
    use parameters_C, only : xfac, alpb
    xfac = 0.d0
    alpb = 0.d0
    alpb(11, 1) =  1.05225212D0  !     Sodium - Hydrogen
    xfac(11, 1) =  1.00000000d0  !     Sodium - Hydrogen   
    alpb(11, 6) =  1.05225212D0  !     Sodium - Carbon
    xfac(11, 6) =  1.00000000d0  !     Sodium - Carbon
    alpb(12, 1) =  1.35052992D0  !     Magnesium - Hydrogen
    xfac(12, 1) =  1.00000000d0  !     Magnesium - Hydrogen
    alpb(12, 6) =  1.48172071D0  !     Magnesium - Hydrogen
    xfac(12, 6) =  1.00000000d0  !     Magnesium - Hydrogen 
    alpb(16,12) =  1.48172071D0  !     Sulfur - Magnesium
    xfac(16,12) =  1.00000000d0  !     Sulfur - Magnesium 
    alpb(13, 1) =  1.38788000D0  !     Aluminum - Hydrogen
    xfac(13, 1) =  1.00000000d0  !     Aluminum - Hydrogen
    alpb(13, 6) =  1.38788000D0  !     Aluminum - Carbon
    xfac(13, 6) =  1.00000000d0  !     Aluminum - Carbon  
    alpb(13,13) =  1.38788000D0  !     Aluminum - Aluminum
    xfac(13,13) =  1.00000000d0  !     Aluminum - Aluminum 
  end subroutine alpb_and_xfac_mndod

