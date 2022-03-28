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

      module reimers_C
! Full module added by RMG
      integer :: n, na, nb2, nham, ix, iy, &
                 nlocl, nstrt, ndump, ion, nstop, &
                 nshell, nconf, multci, n2phot, nese, &
                 nciout, nciouv, nwch, ntrmet, nirreg, nte, &
                 nrep, ifroti = 0, nptg, nr, ndtype, &
                 ncore, nol, noh, nvl, nvh, nov, neltot, nunp, &
                 mfg(24), isok(80), nprin(80), nbfa(80), ndconf, &
                 ixprd(8,8), ig1(74,8), ig2(74,8), ig3(74,8), &
                 ig4(74,8), matind(50000), ifgfac(11), &
                 kflag(3), jflag(3), mb, ic, id, lenf, nfrag, &
                 ncoupb1, ncoupv1, ncoupb2, ncoupv2, &
                 ncisym, icisym(8), nmp2, nmp2t, &
                 ifst, afst, nrange(8), ifbinc(0:8), &
                 ndiff, io(2), iv(2), irrtyp, ns(8), nslwr(8), &
                 nsupr(8), nfact(0:10), nci(8), &
                 norb(4), norbl(4), norbh(4), nex, mspn, moass, &
                 kass(8,4), nold, nvhd, nolc, nvhc

      character*2   ::  iatsym(80)
      character*3   ::  kind(9)
      character*80  ::  line
      double precision  ::  alpha(9,9)
      character*3   ::  nmptg(8), nmrep(8,8), cisym(8)
      character*80  ::  ident
      character*8   ::  ldcomp(3)
      character*40  ::  filenm
      character*8, allocatable :: lab(:)
      logical*1     ::  fastci, ixeof, sing, dbls, mrci

      double precision :: emaxci, vnn, dipsym, value, &
                 zeta(80), zetad(2,80), zetawt(2,80), &
                 weight(80), zcorea(80), betaa(3,80), fg(24,80), &
                 f0(2,2) = 0.d0, tomk,g(74,8), &
                 fssig, fpsig, fppi, fdsig, fdpi, fddel, fall, pp, &
                 bincoe(465), fact(30), a(35), b(35), d(8), e(8), &
                 fspdf(4,4), ss(4), fintfa(6), &
                 avec1(0:4), bvec1(0:4), avec2(0:4), bvec2(0:4), &
                 dia(0:4), tot, edef = 0.d0, ecore, rr, &
                 dq(3), dda(3), dtot(3), pol(6), nel(4), ef(3) = 0.d0


      integer, allocatable :: natm(:), nbf(:), ibf(:), iat(:), &
                 nbt(:), nprn(:), itrmet(:), natt(:), nstr(:), &
                 istr(:,:), nsym(:), ivv(:), npsn(:), mofrag(:), &
                 nfocc(:), nfvir(:), icifrag(:), imp2d(:), imp2s(:), &
                 isc(:), krefd(:), vv(:), iwk(:,:), &
                 istate(:), nseig(:), icif(:), jcif(:), ncif(:), &
                 nbtmo(:), iconf(:), nspn(:)

      double precision, allocatable :: vca(:,:), vcb(:,:), occfr(:), &
                 x(:), y(:), z(:), zcore(:), gamma(:,:), beta(:), &
                 betao(:), s(:), r(:,:), ppg(:,:), pg(:,:), dd(:,:), &
                 ff(:,:), aa(:), cc0(:,:), dm(:,:), xz(:,:), ci(:,:), &
                 evalmo(:), evalci(:), spintr(:,:), cimatr(:), &
                 dmci(:,:), stwt(:), e0(:), wk1(:), wk2(:), tr1(:), &
                 wk0(:,:), wk3(:), occ(:), qgs(:), dipgs(:,:), &
                 qbfcore(:,:), qbf(:,:), qcore(:), q0(:), dip(:,:), &
                 dtmp(:,:)

      double precision, allocatable :: eec(:), ee2(:)

      logical*1, allocatable :: aocc(:,:,:), bocc(:,:,:), aor1(:), &
                 bor1(:), aor(:), bor(:), aos(:), bos(:), aod(:), &
                 bod(:), ao1(:,:), bo1(:,:)

      integer :: ind(0:8), nintg, nbeta

      integer       ::  mis2, mis1, mip2, mip1, mid2, mid1, mid0, &
     &                  mic2, mic1, mic0, &
     &                  mf0ss, mf0sd, mf0dd, mg1sp, mf2pp, &
     &                  mg2sd, mg1pd, mf2pd, mg3pd, mf2dd, mf4dd, &
     &                  mr1sppd, mr2sddd, mr2sdpp

      double precision :: au2ev,au2ang,au2cm,debye

      data au2ev,au2ang,au2cm,debye / &
       27.2114D0, 0.529177D0, 219475.3D0, 4.80294D0 /
      data nmptg /' C1',' CS',' C2','C2V','D2H',' CI','C2H',' D2'/
      data nmrep / 8*'   ',              ' A''',' A"',6*'   ',&
     &             '  A','  B',6*' ',    ' A1',' A2',' B1',' B2',4*' ',&
     &             ' Ag',' Au','B1g','B1u','B2g','B2u','B3g','B3u',&
     &             '  g','  u',6*' ',    ' Ag',' Au',' Bg',' Bu',4*' ',&
     &             '  A',' B1',' B2',    ' B3',4*' ' /
      data ixprd /1,2,3,4,5,6,7,8,2,1,4,3,6,5,8,7,3,4,1,2,7,8,5,6,4,3,2,&
     &1,8,7,6,5,5,6,7,8,1,2,3,4,6,5,8,7,2,1,4,3,7,8,5,6,3,4,1,2,8,7,6,5,&
     &4,3,2,1/
      data (iatsym(i),i=1,80) &
     &            /" H", "He", "Li", "Be", " B", " C", " N", " O", " F", "Ne", &
     & "Na", "Mg", "Al", "Si", "P ", "S ", "Cl", "Ar", "K ", "Ca", "Sc", "Ti", &
     & " V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", &
     & "Kr", "Rb", "Sr", "Y ", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", &
     & "Cd", "In", "Sn", "Sb", "Te", " I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", &
     & "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", &
     & "Hf", "Ta", "W ", "Re", "Os", "Ir", "Pt", "Au", "Hg" /
      data kind    /' S ', ' Px',' Py',' Pz', 'Dz2','xy2','Dxy',&
     &              'Dxz','Dyz'/

      data mfg / 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,&
     &           21,22,23,24/
      data ind  /1,1,1,1,2,2,2,2,2/
      data tomk / 1.2D0 /
      data ifgfac / 3,25,5,15,35,245,49,441,1,1,1 /
      data fintfa / 1.D0,1.267D0,0.585D0,3*1.D0/
      data              ic,id/1,2/
      data ldcomp       /'X Dipole','Y Dipole','Z Dipole'/

      equivalence (mfg(1),mis2), (mfg(2),mis1), (mfg(3),mip2), &
     & (mfg(4),mip1),     (mfg(5),mid2),     (mfg(6),mid1), &
     & (mfg(7),mid0),     (mfg(8),mic2),     (mfg(9),mic1), &
     & (mfg(10),mic0),    (mfg(11),mf0ss),   (mfg(12),mf0sd), &
     & (mfg(13),mf0dd),   (mfg(14),mg1sp),   (mfg(15),mf2pp), &
     & (mfg(16),mg2sd),   (mfg(17),mg1pd),   (mfg(18),mf2pd), &
     & (mfg(19),mg3pd),   (mfg(20),mf2dd),   (mfg(21),mf4dd), &
     & (mfg(22),mr1sppd), (mfg(23),mr2sddd), (mfg(24),mr2sdpp)

      equivalence (fintfa(1),fssig),(fintfa(2),fpsig), &
     &            (fintfa(3),fppi), (fintfa(4),fdsig), &
     &            (fintfa(5),fdpi), (fintfa(6),fddel)

      end module reimers_C
