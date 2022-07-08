program API_example
  use mopac_API
  implicit none
  integer, parameter :: Natom1=3, Natom2=5
  double precision :: XYZ1(3, Natom1), dXYZ1(3, Natom1)
  double precision :: charges1(Natom1),BO1(Natom1,Natom1),energy
  double precision :: XYZ2(3, Natom2), dXYZ2(3, Natom2)
  character(3) :: labels1(Natom1), labels2(Natom2)
  integer :: i,j


  ! perform calculations for two different system sizes; 3 atom H2O and 5 atom NH4+


  ! H2O geometry:
  !  O  0.000000  0.000000  0.000000
  !  H  0.758602  0.000000  0.504284
  !  H  0.758602  0.000000  -0.504284 

  ! define labels
  labels1(1) = 'O'
  labels1(2) = 'H'
  labels1(3) = 'H'

  ! define geometry - dimension of XYZ1 is (Natom, 3)
  XYZ1(:,1) = (/0.d0, 0.d0, 0.d0/)
  XYZ1(:,2) = (/0.758602d0, 0.000000d0,  0.504284d0/)
  XYZ1(:,3) = (/0.758602d0, 0.000000d0,  -0.504284d0/)

  ! define keywords
  keywords='RM1 1SCF XYZ GEO-OK NOSYM'

  ! initialize mopac
  call mopac_initialize(Natom1, XYZ1, labels1)

  ! can now obtain energy and gradient
  write(*,*) 'Water test system'
  write(*,*)
  call mopac_return_energy_and_grad(Natom1, XYZ1, energy, dXYZ1)
  write(*,*) 'First geometry'
  write(*,*) 'The energy is:',energy
  write(*,*) "The gradient is:"
  do i=1,Natom1
     write(*,"(f18.8,$)") (dXYZ1(j,i), j=1,3)
     write(*,*)
  end do
  write(*,*)
  call mopac_return_charges(Natom1, charges1)
  call mopac_return_bonds(Natom1, BO1)
  write(*,*) 'The charges are: ', charges1
  write(*,*) 'Bond orders:'
  do i = 1, Natom1
     write(*,"(f18.8,$)") BO1(i,:)
     write(*,*)
  end do
  write(*,*)
  write(*,*)

  ! change geometry and obtain a new energy+grad
  XYZ1(1,3) = XYZ1(1,3)+.1d0
  call mopac_return_energy_and_grad(Natom1, XYZ1, energy, dXYZ1)
  write(*,*) 'Second geometry, with old density as start'
  write(*,*) 'The energy is:',energy
  write(*,*) "The gradient is:"
  do i=1,Natom1
     write(*,"(f18.8,$)") (dXYZ1(j,i), j=1,3)
     write(*,*)
  end do
  write(*,*)

  ! this result will be slightly different from a direct call to MOPAC
  ! because it is using the density from the previous calculation
  ! the API allows for the density to be reset
  call mopac_return_energy_and_grad(Natom1, XYZ1, energy, dXYZ1,reset_density=.true.)
  write(*,*) 'Second geometry with reset density'
  write(*,*) 'The energy is:',energy
  write(*,*) "The gradient is:"
  do i=1,Natom1
     write(*,"(f18.8,$)") (dXYZ1(j,i), j=1,3)
     write(*,*)
  end do
  write(*,*)

  ! now move to a different system
  ! NH4+ geometry
  !N 	0.0000 	0.0000 	0.0000
  !H 	0.6276 	0.6276 	0.6276
  !H 	0.6276 	-0.6276	-0.6276
  !H 	-0.6276 0.6276 	-0.6276
  !H 	-0.6276 -0.6276 0.6276 
  write(*,*)
  write(*,*) 'NH4+ test system'

  ! define new labels
  labels2(1) = 'N'
  labels2(2) = 'H'
  labels2(3) = 'H'
  labels2(4) = 'H'
  labels2(5) = 'H'

  ! and geometry
  XYZ2(:,1) = 0.d0
  XYZ2(:,2) = (/0.6276, 0.6276, 0.6276/)
  XYZ2(:,3) = (/0.6276, -0.6276, -0.6276/)
  XYZ2(:,4) = (/-0.6276, 0.6276, -0.6276/)
  XYZ2(:,5) = (/-0.6276, -0.6276, 0.6276/)

  ! and keywords
  keywords='RM1 CHARGE=1 1SCF XYZ GEO-OK NOSYM'


  ! need to re-initialize mopac when the system size/charge changes
  call mopac_initialize(Natom2, XYZ2, labels2)

  ! obtain the energy and grad
  call mopac_return_energy_and_grad(Natom2, XYZ2, energy, dXYZ2)
  write(*,*) 'The energy is:',energy
  write(*,*) "The gradient is:"
  do i=1,Natom2
     write(*,"(f18.8,$)") (dXYZ2(j,i), j=1,3)
     write(*,*)
  end do
  
end program API_example
