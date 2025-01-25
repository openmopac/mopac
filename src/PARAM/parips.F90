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
