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

   subroutine dihed(xyz, i, j, k, l, angle)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
   use funcon_C, only : pi
   use common_arrays_C, only: tvec
   use molkst_C, only : id, l11, l21, l31
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
   implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
   integer , intent(in) :: i
   integer , intent(in) :: j
   integer , intent(in) :: k
   integer , intent(in) :: l
   double precision, intent (out)  :: angle
   double precision , intent(in) :: xyz(3,*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
   double precision :: xi1, xj1, xl1, yi1, yj1, yl1, zi1, zj1, zl1, dist, cosa, &
     ddd, yxdist, xi2, xl2, yi2, yl2, costh, sinth, cosph, sinph, yj2, yi3, &
     yl3, rmin_ik, rmin_jk, rmin_lk, Vab_ik(3), Vab_jk(3), Vab_lk(3), Vab(3), &
     Rab
   integer :: ii, jj, kk
!-----------------------------------------------
!********************************************************************
!
!      DIHED CALCULATES THE DIHEDRAL ANGLE BETWEEN ATOMS I, J, K,
!            AND L.  THE CARTESIAN COORDINATES OF THESE ATOMS
!            ARE IN ARRAY XYZ.
!
!     DIHED IS A MODIFIED VERSION OF A SUBROUTINE OF THE SAME NAME
!           WHICH WAS WRITTEN BY DR. W. THIEL IN 1973.
!
!********************************************************************
  if (id == 0) then
     xi1 = xyz(1,i) - xyz(1,k)
     xj1 = xyz(1,j) - xyz(1,k)
     xl1 = xyz(1,l) - xyz(1,k)
     yi1 = xyz(2,i) - xyz(2,k)
     yj1 = xyz(2,j) - xyz(2,k)
     yl1 = xyz(2,l) - xyz(2,k)
     zi1 = xyz(3,i) - xyz(3,k)
     zj1 = xyz(3,j) - xyz(3,k)
     zl1 = xyz(3,l) - xyz(3,k)
   else
     rmin_ik = 1.d8
     rmin_jk = 1.d8
     rmin_lk = 1.d8
     Vab_ik = 0.d0
     Vab_jk = 0.d0
     Vab_lk = 0.d0
     do ii = -l11, l11
       do jj = -l21, l21
         do kk = -l31, l31
           Vab = xyz(:,i) - xyz(:,k) + tvec(:,1)*ii + tvec(:,2)*jj + tvec(:,3)*kk
           Rab = Vab(1)**2 + Vab(2)**2 + Vab(3)**2
           if (Rab < rmin_ik) then
             rmin_ik = Rab
             Vab_ik = Vab
           end if
           Vab = xyz(:,j) - xyz(:,k) + tvec(:,1)*ii + tvec(:,2)*jj + tvec(:,3)*kk
           Rab = Vab(1)**2 + Vab(2)**2 + Vab(3)**2
           if (Rab < rmin_jk) then
             rmin_jk = Rab
             Vab_jk = Vab
           end if
           Vab = xyz(:,l) - xyz(:,k) + tvec(:,1)*ii + tvec(:,2)*jj + tvec(:,3)*kk
           Rab = Vab(1)**2 + Vab(2)**2 + Vab(3)**2
           if (Rab < rmin_lk) then
             rmin_lk = Rab
             Vab_lk = Vab
           end if
         end do
       end do
     end do
     xi1 = Vab_ik(1)
     xj1 = Vab_jk(1)
     xl1 = Vab_lk(1)
     yi1 = Vab_ik(2)
     yj1 = Vab_jk(2)
     yl1 = Vab_lk(2)
     zi1 = Vab_ik(3)
     zj1 = Vab_jk(3)
     zl1 = Vab_lk(3)
   end if
!      ROTATE AROUND Z AXIS TO PUT KJ ALONG Y AXIS
   dist = sqrt(xj1*xj1 + yj1*yj1 + zj1*zj1)
   cosa = zj1/dist
   cosa = min(1.0D0,cosa)
   cosa = dmax1(-1.0D0,cosa)
   ddd = 1.0D0 - cosa**2
   if (ddd <= 0.0D0) go to 10
   yxdist = dist*sqrt(ddd)
   if (yxdist > 1.0D-6) go to 20
    10 continue
   xi2 = xi1
   xl2 = xl1
   yi2 = yi1
   yl2 = yl1
   costh = cosa
   sinth = 0.D0
   go to 30
    20 continue
   cosph = yj1/yxdist
   sinph = xj1/yxdist
   xi2 = xi1*cosph - yi1*sinph
   xl2 = xl1*cosph - yl1*sinph
   yi2 = xi1*sinph + yi1*cosph
   yj2 = xj1*sinph + yj1*cosph
   yl2 = xl1*sinph + yl1*cosph
!      ROTATE KJ AROUND THE X AXIS SO KJ LIES ALONG THE Z AXIS
   costh = cosa
   sinth = yj2/dist
    30 continue
   yi3 = yi2*costh - zi1*sinth
   yl3 = yl2*costh - zl1*sinth
   call dang (xl2, yl3, xi2, yi3, angle)
!     6.2831853  IS 2 * 3.1415926535 = 180 DEGREE
   if (angle < 0.) angle = pi*2.d0 + angle
   if (angle >= 6.28318530717959D0) angle = 0.D0
   return
   end subroutine dihed
