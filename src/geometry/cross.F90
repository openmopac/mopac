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

! vector cross product
subroutine cross(in1, in2, out)
   implicit none
   double precision, dimension(3), intent(in) :: in1, in2
   double precision, dimension(3), intent(out) :: out

   out(1) = in1(2)*in2(3) - in1(3)*in2(2)
   out(2) = in1(3)*in2(1) - in1(1)*in2(3)
   out(3) = in1(1)*in2(2) - in1(2)*in2(1)
end subroutine cross
