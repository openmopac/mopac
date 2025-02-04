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
 
      SUBROUTINE FINISH
!
!   MOPEND SHUTS ALL FILES WHICH MAY HAVE BEEN OPENED
!        AND THEN STARTS A RAPID RETURN TO THE MAIN SEGMENT
!            
      use chanel_C, only: iend, end_fn
      use param_global_C, only : ifiles_8, large
      implicit none
      logical :: exists
      integer :: i
      inquire (file = end_fn, exist = exists)
      if (exists) then
        open(unit=iend, file=end_fn, iostat=i)
        if (i == -100) return
        close(iend, status = 'delete', iostat=i)
        if (i == -100) return
      end if
      end_fn = end_fn(:len_trim(end_fn) - 3)//"res"
      inquire (file = end_fn, exist = exists)
      if (exists) then
        open(unit=iend, file=end_fn, iostat=i)
        if (i == -100) return
        close(iend, status = 'delete', iostat=i)
        if (i == -100) return
      end if
      if (large) write (ifiles_8, '(/,'' == PARAM DONE =='')')
      stop
      END
