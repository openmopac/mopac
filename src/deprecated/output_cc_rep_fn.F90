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

subroutine output_cc_rep_fn
!
!  Output the core-core repulsion function in a form suitable for EXCEL plotting.
!
  use Param_global_C, only : ifiles_8
  use funcon_C, only: a0, ev
  use parameters_C, only: tore, po
  implicit none
  integer :: i, j, ni, nj, pair
  double precision :: r, gab
  double precision, dimension (100,10) :: vab
  double precision, dimension (100) :: rab
  integer, parameter  :: np = 7
  integer, dimension (2,np) :: nij 
  data nij /1,8,6,8,8,8,7,7,1,1,6,1,6,6/
  do pair = 1, np
    do i = 5,60
      r = 0.1d0 * i / a0
      ni = nij(1,pair)
      nj = nij(2,pair)
        gab = ev / Sqrt (r*r/(a0*a0)+ (po(9, ni)+po(9, nj))**2)
      call ccrep (ni, nj, r, vab(i,pair), gab)
    
   !
   ! Multiply by core charges
      vab(i,pair) = 23.06*(vab(i,pair) - tore(ni) * tore(nj) * gab)
       rab(i) = r
    end do
  end do
  write(ifiles_8,"(a)")" Distance    O-H         C-O         O-O         N-N         H-H         C-H         C-C"
  do i = 5, 60
    write(ifiles_8,"(f6.2,9f12.3)")rab(i),(vab(i,j), j=1,np)
  end do
end subroutine output_cc_rep_fn
