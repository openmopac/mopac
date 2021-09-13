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

 
  subroutine pulay_for_gpu(f, p, n, fppf, fock, emat, &
                        & lfock, nfock, msize, start, pl) 
      USE vast_kind_param, ONLY:  double 
      integer  :: n,iopc 
      integer , intent(inout) :: lfock 
      integer , intent(inout) :: nfock 
      integer , intent(in) :: msize 
      real(double) :: pl 
      end
  
subroutine gemma_cublas(txt, txt_2, n, nocc, n_1, one, fmo, n_2, vector , mdim, zero, fck, norbs)
  end
