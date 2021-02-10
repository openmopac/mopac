      subroutine symtrz(vects, eigs, itype, geteig) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE molkst_C, only : numat, norbs, keywrd, moperr
      use meci_C, only : lab, maxci
      use common_arrays_C, only : nfirst, nlast, coord, nat, atmass
      USE symmetry_C, only :  jndex, nclass, elem, namo, jelem
      USE chanel_C, only : iw 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use molsym_I 
      use makopr_I 
      use mult33_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent (in)  :: itype 
      logical, intent(in) :: geteig 
      double precision  :: vects(*) 
      double precision  :: eigs(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k, ierror, i, j, nvecs, nat_store(numat)
      double precision, dimension(3,numat) :: cotim 
      double precision, dimension(3,3) :: r 
!-----------------------------------------------
!**************************************************************
!                                                             *
!     DETERMINE POINT GROUP & SYMMETRIZE ORBITALS             *
!                                                             *
!**************************************************************
!
!    On input:   COORD  = Cartesian coordinates of atoms in system
!                VECTS  = Eigenvectors
!                ITYPE  = 1 for Molecular Orbitals
!                         2 for vibrational normal coordinates
!                         3 for State Functions.
!                GETEIG = .TRUE. if eigenvector analysis wanted.
!    On output   NAME   = Point Group Symbol as 4-character string
!                NAMO   = Irreducible Representations of Eigenvectors
!                         as 4-character strings
!                JNDEX  = Principal Quantum Number of each I.R.
!
!
!**************************************************************
!
!   Store cartesian geometry
!
     
      cotim(:,1:numat) = coord(:,1:numat) 
!
!    Identify Molecular Point Group
!
      nat_store(1:numat) = nat(1:numat)
!
!  Change the atomic number to a function of the isotope.
!  The multiplier "20" was chosen so that an isotope that was different by
!  0.05 amu would be considered a different type of atom.
!  Do not make the multiplier much larger because a large value would upset
!  the detection of symmetry operations.
!
      if (size(jelem) < numat*20) return
      do i = 1, numat
        nat(i) = nat(i)*20 + Nint(atmass(i)*20.d0) 
      end do
      call molsym (cotim, ierror, r) 
      nat(:numat) = nat_store
      if (moperr) return  
      if (allocated(namo)) deallocate (namo)
      if (allocated(jndex)) deallocate (jndex)
      i = max(norbs, maxci + 20,3*numat, lab)
      allocate(namo(i), jndex(i))
      j = size(jndex)
      namo = ' '
      do i = 1, j 
        jndex(i) = i 
      end do 
!
!   A Point-Group has been positively identified.  Now to
!   make the symmetry-operations (as vectors in ELEM).
!
      if ((geteig .or. itype == 2) .and. ierror == 0) call makopr (numat, cotim, ierror, r) 
!
!    Return to the user, if wanted, the symmetry operations for
!    the system.  The operations are relative to the cartesian
!    coordinate system supplied in COORD.  
!
      if (itype == 2 .and. .not.geteig) then 
        if (index(keywrd,'SYMTRZ') /= 0) then 
          write (iw, *) ' Symmetry Operations in SYMTRZ' 
          do i = 1, nclass 
            write (iw, *) ' Operation:', i 
            write (iw, '(3F12.6)') ((elem(j,k,i),j=1,3),k=1,3) 
          end do 
          write (iw, *) ' Orientation Matrix' 
          write (iw, '(3F12.6)') r 
        end if 
        do i = 2, nclass 
          call mult33 (r, i) 
        end do 
        if (index(keywrd,'SYMTRZ') /= 0) then 
          write (iw, *) ' Symmetry Operations for FORCE calculation' 
          write (iw, *) ' Operation:', i 
          do i = 1, nclass 
            write (iw, '(3F12.6)') ((elem(j,k,i),j = 1,3),k = 1,3) 
          end do 
        end if 
      end if 
!
!    Characterize Eigenvectors
!
      if (ierror == 0 .and. geteig) then 
        if (itype == 2) then 
!
!   In vibrational analyses, each atom has three modes, x, y, and z.
!
      jndex(:numat) = 3 
          nvecs = 3*numat 
        else 
!
!    In Molecular Orbital analyses, atoms have bases sets of
!    atomic orbitals.
!
          jndex(:numat) = nlast(:numat) - nfirst(:numat) + 1 
          nvecs = norbs 
        end if 
        if (nvecs >= 1) then 
          i = nvecs 
          if (itype == 3) i = lab 
          call symoir (itype, vects, eigs, nvecs, r, i) 
        end if 
      end if 
      return  
      end subroutine symtrz 
