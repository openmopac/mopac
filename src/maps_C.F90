      module maps_C 
 
      double precision, dimension(:), allocatable  :: &
   &  surf,    &
   &  react

      integer :: &
  &  ijlp    , & !  
  &  ilp,      & !  
  &  jlp,      & !  
  &  jlp1,     & !  
  &  ione,     & !
  &  latom,    & !
  &  lparam,   & !
  &  kloop,    &
  &  latom1,   &
  &  lpara1,   &
  &  latom2,   &
  &  lpara2
      double precision :: &
  &  rxn_coord1,  &
  &  rxn_coord2,  &
  &  rxn_coord,   &
  &  rc_escf,     &  !  Reaction coordinate Heat of Formation
  &  rc_dipo,     &  !  Reaction coordinate dipole
  &  ekin,        &  !  Kinetic energy
  &  dummy
      end module maps_C 
