! Example of MOPAC Fortran API usage
program mopac_api_example
  use mopac_api_f
  implicit none
  type(mopac_system_f) :: water
  type(mopac_state_f) :: state
  type(mopac_properties_f) :: output
  integer :: i

  ! define a water molecule (derived type has default values)  
  water%natom = 3
  allocate(water%atom(water%natom))
  allocate(water%coord(water%natom*3))
  water%atom(1) = 1
  water%coord(1) = -0.02110d0
  water%coord(2) = -0.00200d0
  water%coord(3) = 0.00000d0
  water%atom(2) = 8
  water%coord(4) = 0.83450d0
  water%coord(5) = 0.45190d0
  water%coord(6) = 0.00000d0
  water%atom(3) = 1
  water%coord(7) = 1.47690d0
  water%coord(8) = -0.27300d0
  water%coord(9) = 0.00000d0

  ! self-consistent field (SCF) calculation
  call mopac_scf_f(water, state, output)
  write(*,*) "PM7 heat of water =", output%heat, "kcal/mol"
  if(output%nerror > 0)then
    do i=1,output%nerror
       write(*,*) "SCF ERROR #", i, ":", output%error_msg(i)
    end do
    stop 1
  end if

  ! SCF calculation in COSMO solvent
  water%epsilon = 78.4d0
  call mopac_scf_f(water, state, output)
  if(output%nerror > 0)then
    do i=1,output%nerror
       write(*,*) "COSMO ERROR #", i, ":", output%error_msg(i)
    end do
    stop 1 
  end if
  write(*,*) "PM7 heat of water in COSMO solvent =", output%heat, "kcal/mol"

  ! geometry relaxation
  water%epsilon = 1.0d0 ! reset dielectric constant back to vacuum
  water%natom_move = 3
  call mopac_relax_f(water, state, output)
  if(output%nerror > 0)then
    do i=1,output%nerror
       write(*,*) "RELAX ERROR #", i, ":", output%error_msg(i)
    end do
    stop 1
  end if
  write(*,*) "geometry relaxation:"
  do i=1,water%natom*3
    write(*,*) i, ":", water%coord(i), "->", output%coord_update(i)
    water%coord(i) = output%coord_update(i)
  end do

  ! vibrational calculation of relaxed geometry
  call mopac_vibe_f(water, state, output)
  if(output%nerror > 0)then
    do i=1,output%nerror
       write(*,*) "VIBE ERROR #", i, ":", output%error_msg(i)
    end do
    stop 1
  end if
  write(*,*) "vibrational modes:"
  do i=1,water%natom_move*3
    write(*,*) i, ":", output%freq(i), "1/cm"
  end do
end program mopac_api_example
