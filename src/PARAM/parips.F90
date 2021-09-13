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

double precision function parips (eigs, loop)
    use molkst_C, only : norbs, nopen, nclose, nalpha
    use param_global_C, only : lions
    implicit none
    integer, intent(in) :: loop
    double precision, dimension (norbs), intent (in) :: eigs
    if(nclose /= 0) then
!
!  Only allow for higher IPs for closed-shell systems
!
      parips = -eigs(lions(loop))
      if (nopen /= nclose) then
           parips = Min (parips,-eigs(nopen))
      end if
    elseif (nalpha > 0) then
      parips = -eigs(nalpha)
    else
      parips = 0.d0
    end if
end function parips
