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

subroutine adjvec (cvecb, ncvb, icvecb, nib, nncb, ncb_loc, nnb, ncvecb, lmob, &
& iorbs, cveca, ncva, icveca, nia, nnca, nca_loc, nna_loc, ncveca, lmoa, beta, &
& iused, sumtot)
    use molkst_C, only: numat, norbs
    use MOZYME_C, only : thresh
    implicit none
    integer, intent (in) :: lmoa, lmob, ncva, ncvb, nia, nib, nna_loc, nnb
    double precision, intent (in) :: beta
    double precision, intent (inout) :: sumtot
    integer, dimension (nia), intent (in) :: icveca
    integer, dimension (nib), intent (inout) :: icvecb
    integer, dimension (nna_loc), intent (in) :: nca_loc, ncveca, nnca
    integer, dimension (nnb), intent (in) :: ncvecb, nncb
    integer, dimension (nnb), intent (inout) :: ncb_loc
    integer, dimension (numat), intent (inout) :: iused
    integer, dimension (norbs), intent (in) :: iorbs
    double precision, dimension (ncva), intent (in) :: cveca
    double precision, dimension (ncvb), intent (inout) :: cvecb
!
    integer :: ii, la, lb, llim, mia, mla, mlb, mlim, mlla, mllb
    double precision :: cutoff, sum
    mlla = 0
    cutoff = thresh * 1.d1
    if (Abs (beta) < cutoff) return
    sumtot = sumtot + Abs (beta)
   !
   !   IDENTIFY THE ATOMS IN THE SECOND SET OF LMOs.  ATOMS NOT USED ARE
   !   FLAGGED BY '-1' IN IUSED.
   !
    do la = nnca(lmoa) + 1, nnca(lmoa) + nca_loc(lmoa)
      iused(icveca(la)) = -1
    end do
    mlb = ncvecb(lmob)
    if (lmob == nnb) then
      llim = nib
    else
      llim = nncb(lmob+1)
    end if
    if (lmob == nnb) then
      mlim = ncvb - 4
    else
      mlim = ncvecb(lmob+1) - 4
    end if
    do lb = nncb(lmob) + 1, nncb(lmob) + ncb_loc(lmob)
      ii = icvecb(lb)
      iused(ii) = mlb
      mlb = mlb + iorbs(ii)
    end do
    mla = ncveca(lmoa)
   !
   !      ROTATION OF SECOND VECTOR TO MAKE IT ORTHOGONAL TO THE FIRST.
   !
    do la = nnca(lmoa) + 1, nnca(lmoa) + nca_loc(lmoa)
      mia = icveca(la)
      if (iused(mia) >= 0) then
        mllb = iused(mia)
         !
         !   ROTATION INVOLVING ATOMS WHICH ARE COMMON TO BOTH LMOs
         !
        do mlla = mla + 1, mla + iorbs(mia)
          mllb = mllb + 1
          cvecb(mllb) = cvecb(mllb) - beta * cveca(mlla)
        end do
      else
         !
         !   ATOM 'MIA' DOES NOT EXIST IN LMOB.  CHECK IF IT SHOULD EXIST.
         !   IF IT SHOULD EXIST, MAKE IT EXIST.
         !
        sum = 0.d0
        do mlla = mla + 1, mla + iorbs(mia)
          sum = sum + cveca(mlla) ** 2
        end do
        if (beta**2*sum > cutoff) then
          if (ncb_loc(lmob) < llim .and. mlb < mlim) then
            ncb_loc(lmob) = ncb_loc(lmob) + 1
            icvecb(nncb(lmob)+ncb_loc(lmob)) = mia
            iused(mia) = mlb
            do mlla = mla + 1, mla + iorbs(mia)
              mlb = mlb + 1
              cvecb(mlb) = -beta * cveca(mlla)
            end do
          end if
        end if
      end if
      mla = mla + iorbs(mia)
    end do
    if (mlla /= -1) return
!
!   Use in debugging - sum is the value of the overlap after orthogonalization
!
    sum = 0
    mla = ncveca(lmoa)
    do la = nnca(lmoa) + 1, nnca(lmoa) + nca_loc(lmoa)
      mia = icveca(la)
      if (iused(mia) >= 0) then
        mllb = iused(mia)
        do mlla = mla + 1, mla + iorbs(mia)
          mllb = mllb + 1
          sum = sum + cvecb(mllb)*cveca(mlla)
        end do
      end if
      mla = mla + iorbs(mia)
    end do
end subroutine adjvec
